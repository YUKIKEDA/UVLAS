#include "pchip_boundary.hpp"

#include "pchip_curve_plan.hpp"
#include "pchip_slot_private.hpp"

#include <algorithm>
#include <bit>
#include <type_traits>

namespace uvlas::aero::detail::pchip {
namespace boundary_detail {
using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::supplement;

struct Visits {
    std::uint64_t limit, count = 0;
    const std::stop_token& stop;
    core::FailureDetails details;

    core::Result<void> next() noexcept {
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        if (count == limit)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        ++count;
        return {};
    }
};

core::Result<BitBound> shape(RationalView value, RationalLimits limits,
                             core::FailureDetails details) noexcept {
    const auto words = limits.input_bits / 32 + (limits.input_bits % 32 != 0 ? 1U : 0U);
    if (value.numerator.size() > words || value.denominator.size() > words)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    if (value.denominator.empty() || value.denominator.back() == 0 ||
        (!value.numerator.empty() && value.numerator.back() == 0) ||
        (value.numerator.empty() ? value.sign != 0 : (value.sign != -1 && value.sign != 1)))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    CheckedCount checked;
    const auto bits = [&](std::span<const std::uint32_t> magnitude) {
        return magnitude.empty()
                   ? std::uint64_t{0}
                   : checked.add(checked.mul(magnitude.size() - 1, 32),
                                 static_cast<std::uint64_t>(std::bit_width(magnitude.back())));
    };
    const BitBound bound{bits(value.numerator), bits(value.denominator)};
    if (!checked.ok || std::max(bound.numerator, bound.denominator) > limits.input_bits)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    return bound;
}

// The counter uses the same grid already admitted by the coefficient forecast.
struct BoundedCoefficients {
    SampleGrid grid;
    BitBound coefficient_width;
    [[nodiscard]] bool active() const noexcept {
        return true;
    }
    [[nodiscard]] std::uint64_t nodes() const noexcept {
        return grid.times.size();
    }
    [[nodiscard]] std::uint64_t components() const noexcept {
        return grid.components;
    }
    [[nodiscard]] RationalView time(std::uint64_t node) const noexcept {
        return grid.times[static_cast<std::size_t>(node)];
    }
    [[nodiscard]] RationalView value(std::uint64_t node, std::uint64_t component) const noexcept {
        return grid.values[static_cast<std::size_t>(node * grid.components + component)];
    }
    [[nodiscard]] BitBound coefficient(std::uint64_t, std::uint64_t, std::uint64_t) const noexcept {
        return coefficient_width;
    }
};

template <class Arithmetic, class Source>
class Program {
    static constexpr bool counting = std::is_same_v<Arithmetic, Counter>;
    using Value = std::conditional_t<counting, BitBound, Slot>;

public:
    Program(Arithmetic& arithmetic, Visits& visits, RationalLimits limits,
            const Source& source) noexcept
        : arithmetic_(arithmetic),
          visits_(visits),
          limits_(limits),
          source_(source),
          last_(source.nodes() - 1) {}

    core::Result<void> check(const BoundaryRequest& request) noexcept {
        auto begin = read(request.begin);
        if (!begin)
            return std::unexpected(begin.error());
        auto end = read(request.end);
        if (!end)
            return std::unexpected(end.error());
        const auto* periodic = std::get_if<PeriodicBoundary>(&request.policy);
        if (periodic != nullptr) {
            auto period = read(periodic->period);
            if (!period)
                return std::unexpected(period.error());
            return domain(*begin, *end, request.policy, &*period);
        }
        return domain(*begin, *end, request.policy, nullptr);
    }

private:
    Arithmetic& arithmetic_;
    Visits& visits_;
    RationalLimits limits_;
    const Source& source_;
    std::uint64_t last_;

    core::Result<Value> read(RationalView value) noexcept {
        const auto visited = visits_.next();
        if (!visited)
            return std::unexpected(visited.error());
        if constexpr (counting) {
            const auto bound = shape(value, limits_, visits_.details);
            if (!bound)
                return std::unexpected(bound.error());
            return arithmetic_.input(*bound);
        } else {
            return arithmetic_.input(value);
        }
    }

    core::Result<Value> read(core::Result<RationalView> value) noexcept {
        if (!value)
            return std::unexpected(supplement(value.error(), visits_.details));
        return read(*value);
    }

    core::Result<Value> read(BitBound bound) noexcept
        requires counting
    {
        return visits_.next().and_then([&] { return arithmetic_.input(bound); });
    }

    core::Result<Value> integer(std::int64_t value) noexcept {
        return visits_.next().and_then([&] { return arithmetic_.integer(value); });
    }

    core::Result<Value> binary(RationalOp op, const Value& a, const Value& b) noexcept {
        return visits_.next().and_then([&] { return arithmetic_.binary(op, a, b); });
    }

    core::Result<std::strong_ordering> compare(const Value& a, const Value& b) noexcept {
        const auto visited = visits_.next();
        if (!visited)
            return std::unexpected(visited.error());
        if constexpr (counting) {
            return arithmetic_.compare(a, b).transform([] { return std::strong_ordering::equal; });
        } else {
            return arithmetic_.compare(a, b);
        }
    }

    core::Result<void> domain(const Value& begin, const Value& end,
                              const decltype(BoundaryRequest::policy)& policy,
                              const Value* period) noexcept {
        const auto order = compare(begin, end);
        if (!order)
            return std::unexpected(order.error());
        if (*order > 0)
            return invalid();

        const auto first = read(source_.time(0));
        if (!first)
            return std::unexpected(first.error());
        const auto last = read(source_.time(last_));
        if (!last)
            return std::unexpected(last.error());

        if (period != nullptr)
            return seam(*period, *first, *last);

        const auto left = compare(begin, *first);
        if (!left)
            return std::unexpected(left.error());
        const auto right = compare(end, *last);
        if (!right)
            return std::unexpected(right.error());
        if (std::holds_alternative<RejectBoundary>(policy))
            return *left < 0 || *right > 0 ? invalid() : core::Result<void>{};

        const auto through_first = compare(end, *first);
        if (!through_first)
            return std::unexpected(through_first.error());
        const auto through_last = compare(begin, *last);
        if (!through_last)
            return std::unexpected(through_last.error());
        // Counting treats comparisons as equal, covering both endpoints for Hold.
        const auto first_contact = *left <= 0 && *through_first >= 0;
        const auto last_contact = *through_last <= 0 && *right >= 0;
        if (!first_contact && !last_contact)
            return {};
        return hold_slopes(first_contact, last_contact);
    }

    core::Result<void> hold_slopes(bool first_contact, bool last_contact) noexcept {
        const auto zero = integer(0);
        if (!zero)
            return std::unexpected(zero.error());
        for (std::uint64_t component = 0; component < source_.components(); ++component) {
            if (first_contact) {
                const auto slope = read(source_.coefficient(0, component, 1));
                if (!slope)
                    return std::unexpected(slope.error());
                const auto valid = equal(*slope, *zero);
                if (!valid)
                    return std::unexpected(valid.error());
            }
            if (last_contact) {
                const auto slope = last_u_derivative(component);
                if (!slope)
                    return std::unexpected(slope.error());
                const auto valid = equal(*slope, *zero);
                if (!valid)
                    return std::unexpected(valid.error());
            }
        }
        return {};
    }

    core::Result<void> seam(const Value& period, const Value& first, const Value& last) noexcept {
        // The owned grid has increasing knots, so matching its width also proves positivity.
        const auto width = binary(RationalOp::Subtract, last, first);
        if (!width)
            return std::unexpected(width.error());
        const auto period_matches = equal(period, *width);
        if (!period_matches)
            return std::unexpected(period_matches.error());
        const auto second = read(source_.time(1));
        if (!second)
            return std::unexpected(second.error());
        const auto previous = read(source_.time(last_ - 1));
        if (!previous)
            return std::unexpected(previous.error());
        const auto first_width = binary(RationalOp::Subtract, *second, first);
        if (!first_width)
            return std::unexpected(first_width.error());
        const auto last_width = binary(RationalOp::Subtract, last, *previous);
        if (!last_width)
            return std::unexpected(last_width.error());

        for (std::uint64_t component = 0; component < source_.components(); ++component) {
            const auto first_value = read(source_.value(0, component));
            if (!first_value)
                return std::unexpected(first_value.error());
            const auto last_value = read(source_.value(last_, component));
            if (!last_value)
                return std::unexpected(last_value.error());
            const auto value_matches = equal(*first_value, *last_value);
            if (!value_matches)
                return std::unexpected(value_matches.error());
            // Compare physical derivatives by their positive interval widths.
            const auto first_u = read(source_.coefficient(0, component, 1));
            if (!first_u)
                return std::unexpected(first_u.error());
            const auto last_u = last_u_derivative(component);
            if (!last_u)
                return std::unexpected(last_u.error());
            const auto first_scaled = binary(RationalOp::Multiply, *first_u, *last_width);
            if (!first_scaled)
                return std::unexpected(first_scaled.error());
            const auto last_scaled = binary(RationalOp::Multiply, *last_u, *first_width);
            if (!last_scaled)
                return std::unexpected(last_scaled.error());
            const auto slope_matches = equal(*first_scaled, *last_scaled);
            if (!slope_matches)
                return std::unexpected(slope_matches.error());
        }
        return {};
    }

    core::Result<Value> last_u_derivative(std::uint64_t component) noexcept {
        // Horner derivative at u=1: c1 + 2*c2 + 3*c3.
        const auto c1 = read(source_.coefficient(last_ - 1, component, 1));
        if (!c1)
            return std::unexpected(c1.error());
        const auto twice_c2 = scaled_last_coefficient(component, 2);
        if (!twice_c2)
            return std::unexpected(twice_c2.error());
        const auto three_c3 = scaled_last_coefficient(component, 3);
        if (!three_c3)
            return std::unexpected(three_c3.error());
        return binary(RationalOp::Add, *c1, *twice_c2).and_then([&](const Value& sum) {
            return binary(RationalOp::Add, sum, *three_c3);
        });
    }

    core::Result<Value> scaled_last_coefficient(std::uint64_t component,
                                                std::uint64_t power) noexcept {
        const auto factor = integer(static_cast<std::int64_t>(power));
        if (!factor)
            return std::unexpected(factor.error());
        return read(source_.coefficient(last_ - 1, component, power))
            .and_then([&](const Value& coefficient) {
                return binary(RationalOp::Multiply, coefficient, *factor);
            });
    }

    core::Result<void> equal(const Value& a, const Value& b) noexcept {
        return compare(a, b).and_then(
            [&](const auto ordering) { return ordering == 0 ? core::Result<void>{} : invalid(); });
    }

    [[nodiscard]] core::Result<void> invalid() const noexcept {
        return std::unexpected(failure(core::FailureCode::InvalidInput, visits_.details));
    }
};

struct Builder {
    struct Lease {
        BoundaryWorkspace& workspace;
        ~Lease() {
            workspace.busy_.clear(std::memory_order_release);
        }
    };

    template <class Source>
    static core::Result<BoundaryPlan> forecast(const Source& source, const BoundaryRequest& request,
                                               RationalLimits limits, const std::stop_token& stop,
                                               core::FailureDetails details,
                                               std::uint64_t derived_bits = 0) noexcept {
        if (!source.active() || source.nodes() < 2 || source.components() == 0)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        // Check the caller's shapes before a requested stop. Canonical GCD checks are
        // still part of admitted execution.
        const auto* periodic = std::get_if<PeriodicBoundary>(&request.policy);
        const auto arguments = std::array{request.begin, request.end,
                                          periodic == nullptr ? RationalView{} : periodic->period};
        const std::uint64_t argument_count = periodic == nullptr ? 2U : 3U;
        for (std::size_t i = 0; i < argument_count; ++i) {
            const auto bound = shape(arguments[i], limits, details);
            if (!bound)
                return std::unexpected(bound.error());
        }
        Visits scan{limits.steps, 0, stop, details};
        auto arithmetic_limits = limits;
        arithmetic_limits.input_bits = std::max(limits.input_bits, derived_bits);
        Counter counter{arithmetic_limits, details};
        Program program{counter, scan, limits, source};
        const auto checked = program.check(request);
        if (!checked)
            return std::unexpected(checked.error());
        const auto arithmetic = counter.finish();
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        CheckedCount count;
        const auto outer = count.add(count.mul(scan.count, 2), argument_count + 1);
        const auto work = count.add(outer, arithmetic->reserved_steps());
        if (!count.ok || work > limits.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        BoundaryPlan result;
        result.arithmetic_ = *arithmetic;
        result.scan_ = count.add(scan.count, argument_count);
        result.work_ = work;
        return result;
    }

    static core::Result<CurveNumericsPlan> forecast_curve(const PreparationRequest& source,
                                                          RationalLimits limits,
                                                          const std::stop_token& stop,
                                                          core::FailureDetails details,
                                                          bool retain_original) noexcept {
        const auto coefficients =
            forecast_coefficients(source.samples, limits, stop, details, retain_original);
        if (!coefficients)
            return std::unexpected(coefficients.error());
        const auto bits = coefficients->arithmetic().integer_bits();
        const BoundedCoefficients bounded{source.samples, {bits, bits}};
        auto boundary_limits = limits;
        boundary_limits.steps -= coefficients->reserved_steps();
        const auto boundary =
            forecast(bounded, source.boundary, boundary_limits, stop, details, bits);
        if (!boundary)
            return std::unexpected(boundary.error());
        CheckedCount count;
        const auto bytes = count.add(coefficients->payload_bytes(), boundary->payload_bytes());
        const auto forecasts =
            count.add(coefficients->preflight_steps(), boundary->preflight_steps());
        const auto calls = count.add(coefficients->reserved_steps(), boundary->reserved_steps());
        const auto work = count.add(retain_original ? calls : count.add(forecasts, calls), 1);
        if (!count.ok || work > limits.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));

        CurveNumericsPlan result;
        result.coefficients_ = *coefficients;
        result.boundary_ = *boundary;
        result.boundary_limits_ = boundary_limits;
        result.boundary_limits_.input_bits = std::max(limits.input_bits, bits);
        result.bytes_ = bytes;
        result.work_ = work;
        result.preflight_ = forecasts + 1;
        return result;
    }

    static core::Result<BoundaryUsage>
    check(const CoefficientBlock& source, const BoundaryRequest& request,
          BoundaryWorkspace& workspace, RationalLimits limits, core::ResourcePlan& memory,
          const std::stop_token& stop, core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto plan = forecast(source, request, limits, stop, details);
        if (!plan)
            return std::unexpected(plan.error());
        auto arithmetic = Workspace::create(plan->arithmetic_, memory, details, stop);
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        return run(source, request, *plan, *arithmetic, limits, stop, details, plan->scan_);
    }

    static core::Result<BoundaryUsage>
    run(const CoefficientBlock& source, const BoundaryRequest& request, const BoundaryPlan& plan,
        Workspace& arithmetic, RationalLimits limits, const std::stop_token& stop,
        core::FailureDetails details, std::uint64_t preflight = 0) noexcept {
        auto session = arithmetic.begin(stop, details);
        if (!session)
            return std::unexpected(session.error());
        const auto outer = plan.work_ - plan.arithmetic_.reserved_steps();
        Visits visits{outer, preflight, stop, details};
        Program program{*session, visits, limits, source};
        const auto checked = program.check(request);
        if (!checked)
            return std::unexpected(checked.error());
        const auto completion = visits.next();
        if (!completion)
            return std::unexpected(completion.error());
        const auto usage = session->usage();
        if (!usage)
            return std::unexpected(usage.error());
        return BoundaryUsage{usage->steps + visits.count, usage->slots};
    }
};
} // namespace boundary_detail

core::Result<BoundaryPlan> forecast_boundary(const CoefficientBlock& source,
                                             BoundaryRequest request, RationalLimits limits,
                                             const std::stop_token& stop,
                                             core::FailureDetails details) noexcept {
    return boundary_detail::Builder::forecast(source, request, limits, stop, details);
}

core::Result<CurveNumericsPlan> forecast_curve_numerics(const PchipPreparationInput& input,
                                                        RationalLimits limits,
                                                        const std::stop_token& stop,
                                                        core::FailureDetails details) noexcept {
    const auto source = input.view();
    if (!source)
        return std::unexpected(slot_detail::supplement(source.error(), details));
    return boundary_detail::Builder::forecast_curve(*source, limits, stop, details, false);
}

core::Result<BoundaryUsage>
check_prepared_boundary(const CoefficientBlock& source, BoundaryRequest request,
                        const BoundaryPlan& plan, Workspace& arithmetic, RationalLimits limits,
                        const std::stop_token& stop, core::FailureDetails details) noexcept {
    return boundary_detail::Builder::run(source, request, plan, arithmetic, limits, stop, details);
}

core::Result<BoundaryUsage> check_boundary(const CoefficientBlock& source, BoundaryRequest request,
                                           BoundaryWorkspace& workspace, RationalLimits limits,
                                           core::ResourcePlan& memory, const std::stop_token& stop,
                                           core::FailureDetails details) noexcept {
    return boundary_detail::Builder::check(source, request, workspace, limits, memory, stop,
                                           details);
}
core::Result<CurveNumericsPlan> forecast_curve_numerics(const PreparationRequest& source,
                                                        RationalLimits limits,
                                                        const std::stop_token& stop,
                                                        core::FailureDetails details) noexcept {
    return boundary_detail::Builder::forecast_curve(source, limits, stop, details, true);
}
} // namespace uvlas::aero::detail::pchip
