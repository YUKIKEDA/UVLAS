#include "backend_memory_detail.hpp"
#include "compensated.hpp"
#include "gmres_residual_detail.hpp"
#include "linear_solve_detail.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <uvlas/numerics/gmres_solve.hpp>

namespace uvlas::numerics::linear {
namespace {
using core::FailureCode;
using core::Index;
using core::Result;
constexpr auto operation = core::Operation::SolveLinear;
core::FailureDetails known(const Context& context) noexcept {
    core::FailureDetails result{};
    result.model = context.cache.model;
    result.trial = context.trial;
    result.numerical_generation = context.cache.matrix;
    return result;
}
auto fail(FailureCode code, const Context& context) noexcept {
    return std::unexpected(core::Failure{code, operation, known(context)});
}
core::Failure internal_failure(core::Failure failure, const Context& context) noexcept {
    auto& d = failure.details;
    if (!d.model)
        d.model = context.cache.model;
    if (!d.trial)
        d.trial = context.trial;
    if (!d.numerical_generation)
        d.numerical_generation = context.cache.matrix;
    return failure;
}
const DenseReal& high(const PrecisionBuffer<double>& x) noexcept {
    if (const auto* binary = std::get_if<Binary64<double>>(&x.storage()))
        return binary->values;
    return std::get_if<TwoWord<double>>(&x.storage())->hi;
}
const DenseReal* low(const PrecisionBuffer<double>& x) noexcept {
    const auto* words = std::get_if<TwoWord<double>>(&x.storage());
    return words ? &words->lo : nullptr;
}
Limits shared_limits(const GmresLimits& limits) noexcept {
    return {limits.maximum_order, limits.maximum_rhs_columns, 0, 0, limits.backend_memory};
}
Result<void> check_limits(const GmresLimits& limits, const Context& context) noexcept {
    if (limits.maximum_order < 0 || limits.maximum_rhs_columns < 0 || limits.restart <= 0)
        return fail(FailureCode::InvalidInput, context);
    return {};
}
Result<void> state_matches(const GmresOperator& op, const FixedPreconditioner& p,
                           const Context& context, const FixedPreconditionerIdentity& fixed,
                           Index n) noexcept {
    if (op.context() != context || p.identity() != fixed || op.order() != n || p.order() != n)
        return fail(FailureCode::StaleGeneration, context);
    return {};
}
Result<double> norm(std::span<const double> x, const Context& context) noexcept {
    double value = 0;
    for (const double item : x) {
        if (!std::isfinite(item))
            return fail(FailureCode::Nonfinite, context);
        value = std::hypot(value, item);
    }
    if (!std::isfinite(value))
        return fail(FailureCode::PrecisionUnresolved, context);
    return value;
}
Result<double> dot(std::span<const double> a, std::span<const double> b,
                   const Context& context) noexcept {
    uvlas::numerics::detail::CompensatedAccumulator sum;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double product = a[i] * b[i];
        if (!std::isfinite(product))
            return fail(FailureCode::PrecisionUnresolved, context);
        auto r = sum.add(product);
        if (!r)
            return std::unexpected(internal_failure(r.error(), context));
    }
    auto value = sum.result();
    if (!value)
        return std::unexpected(internal_failure(value.error(), context));
    // Search-space arithmetic is approximate; only the independent original-system gate certifies.
    return value->value;
}
Result<Index> product(Index a, Index b) noexcept {
    return core::checked_extent_product(std::array<Index, 2>{a, b}, operation);
}
} // namespace
namespace detail {
class GmresFactory {
public:
    static GmresFailure failure(core::Failure value, GmresDiagnostics diagnostics) noexcept {
        return {value, diagnostics};
    }
    static GmresSolution success(Solution<double> solution, GmresDiagnostics diagnostics) noexcept {
        return {std::move(solution), diagnostics};
    }
    template <class Value, class Function>
    static GmresResult<Value> guard(const Context& context, GmresDiagnostics& diagnostics,
                                    Function&& action) noexcept {
        auto result = core::guard_entry(operation, FailureCode::LinearUnconverged, known(context),
                                        std::forward<Function>(action));
        if (!result)
            return std::unexpected(failure(internal_failure(result.error(), context), diagnostics));
        return std::move(*result);
    }
};
struct GmresData {
    GmresData(Context c, FixedPreconditionerIdentity p, GmresLimits limits,
              BackendAdmission memory) noexcept
        : context(c),
          preconditioner(p),
          capacity(limits),
          admission(memory),
          lease(core::PlannedBuffer<std::byte>{}) {}
    Context context;
    FixedPreconditionerIdentity preconditioner;
    GmresLimits capacity;
    BackendAdmission admission;
    MetadataBuffers scales;
    core::Workspace lease;
    Index width = 0;
    Index order = 0;
    core::PlannedBuffer<double> basis, hessenberg, cosines, sines, g, y;
    core::PlannedBuffer<double> residual, residual_error, denominator, rhs_hi, rhs_lo, word_hi,
        word_lo, tail_hi, tail_lo, error_hi, error_lo, weights, work, temp, delta, proposal;
};
} // namespace detail
namespace {
bool same_scales(const Scales& scales, const detail::MetadataBuffers& saved) noexcept {
    if (scales.unknown.size() != saved.unknown_scales.size())
        return false;
    for (std::size_t i = 0; i < scales.unknown.size(); ++i)
        if (std::bit_cast<std::uint64_t>(scales.unknown[i]) !=
                std::bit_cast<std::uint64_t>(saved.unknown_scales.values()[i]) ||
            std::bit_cast<std::uint64_t>(scales.residual[i]) !=
                std::bit_cast<std::uint64_t>(saved.residual_scales.values()[i]) ||
            scales.unknown_units[i] != saved.unknown_units.values()[i] ||
            scales.residual_units[i] != saved.residual_units.values()[i])
            return false;
    return true;
}
class Runner {
public:
    Runner(detail::GmresData& data, const GmresOperator& a, const FixedPreconditioner& p,
           const GmresLimits& limits, GmresDiagnostics& diagnostics) noexcept
        : d(data), op(a), preconditioner(p), limits_(limits), stats(diagnostics) {}
    [[nodiscard]] Result<void> stable() const noexcept {
        return state_matches(op, preconditioner, d.context, d.preconditioner, d.order);
    }
    Result<void> admit(std::uint64_t requested) noexcept {
        if (requested > limits_.maximum_operator_actions - stats.operator_actions) {
            stats.stop = GmresStopReason::OperatorLimit;
            auto details = known(d.context);
            if (requested <= std::numeric_limits<std::uint64_t>::max() - stats.operator_actions)
                details.resources =
                    core::ResourceUsage{core::ResourceKind::OperatorApplications,
                                        stats.operator_actions + requested, stats.operator_actions};
            return std::unexpected(
                core::Failure{FailureCode::LinearUnconverged, operation, details});
        }
        return {};
    }
    Result<void> charge(bool absolute, bool p) noexcept {
        auto fixed = stable();
        if (!fixed)
            return fixed;
        const auto maximum = std::numeric_limits<std::uint64_t>::max();
        if (stats.total_action_callbacks == maximum ||
            (p ? stats.preconditioner_actions == maximum : stats.operator_actions == maximum) ||
            (absolute ? stats.absolute_actions == maximum : !p && stats.forward_actions == maximum))
            return fail(FailureCode::ResourceExceeded, d.context);
        if (!p) {
            auto room = admit(1);
            if (!room)
                return room;
        }
        ++stats.total_action_callbacks;
        if (p)
            ++stats.preconditioner_actions;
        else {
            ++stats.operator_actions;
            if (absolute)
                ++stats.absolute_actions;
            else
                ++stats.forward_actions;
        }
        return {};
    }
    Result<void> action(std::span<const double> x, std::span<double> estimate,
                        std::span<double> tail, std::span<double> error) noexcept {
        for (const double value : x)
            if (!std::isfinite(value))
                return fail(FailureCode::PrecisionUnresolved, d.context);
        auto charged = charge(false, false);
        if (!charged)
            return charged;
        stats.stop = GmresStopReason::ActionFailure;
        std::fill(estimate.begin(), estimate.end(), std::numeric_limits<double>::quiet_NaN());
        std::fill(error.begin(), error.end(), std::numeric_limits<double>::quiet_NaN());
        std::fill(tail.begin(), tail.end(), std::numeric_limits<double>::quiet_NaN());
        auto result = core::guard_entry(operation, FailureCode::LinearUnconverged, known(d.context),
                                        [&] { return op.apply_word(x, estimate, tail, error); });
        auto after = stable();
        if (!result)
            return result;
        if (!after)
            return after;
        for (std::size_t i = 0; i < estimate.size(); ++i) {
            if (!std::isfinite(estimate[i]) || !std::isfinite(tail[i]) || !std::isfinite(error[i]))
                return fail(FailureCode::Nonfinite, d.context);
            if (error[i] < 0)
                return fail(FailureCode::InvalidInput, d.context);
        }
        stats.stop = GmresStopReason::ArithmeticFailure;
        return {};
    }
    Result<void> absolute(std::span<const double> weights, std::span<double> result) noexcept {
        auto charged = charge(true, false);
        if (!charged)
            return charged;
        stats.stop = GmresStopReason::ActionFailure;
        std::fill(result.begin(), result.end(), std::numeric_limits<double>::quiet_NaN());
        auto outcome =
            core::guard_entry(operation, FailureCode::LinearUnconverged, known(d.context),
                              [&] { return op.absolute_action_lower(weights, result); });
        auto after = stable();
        if (!outcome)
            return outcome;
        if (!after)
            return after;
        for (const double value : result) {
            if (!std::isfinite(value))
                return fail(FailureCode::Nonfinite, d.context);
            if (value < 0)
                return fail(FailureCode::InvalidInput, d.context);
        }
        stats.stop = GmresStopReason::ArithmeticFailure;
        return {};
    }
    Result<void> precondition(std::span<const double> residual,
                              std::span<double> correction) noexcept {
        for (const double value : residual)
            if (!std::isfinite(value))
                return fail(FailureCode::PrecisionUnresolved, d.context);
        auto charged = charge(false, true);
        if (!charged)
            return charged;
        stats.stop = GmresStopReason::ActionFailure;
        std::fill(correction.begin(), correction.end(), std::numeric_limits<double>::quiet_NaN());
        auto result = core::guard_entry(operation, FailureCode::LinearUnconverged, known(d.context),
                                        [&] { return preconditioner.apply(residual, correction); });
        auto after = stable();
        if (!result)
            return result;
        if (!after)
            return after;
        for (const double value : correction)
            if (!std::isfinite(value))
                return fail(FailureCode::Nonfinite, d.context);
        stats.stop = GmresStopReason::ArithmeticFailure;
        return {};
    }
    void qualification() noexcept {
        stats.stop = GmresStopReason::QualificationFailure;
    }
    detail::GmresData& d;

private:
    const GmresOperator& op;
    const FixedPreconditioner& preconditioner;
    const GmresLimits& limits_;
    GmresDiagnostics& stats;
};
} // namespace
GmresPrepared::GmresPrepared(std::unique_ptr<detail::GmresData> data) noexcept
    : data_(std::move(data)) {}
GmresPrepared::GmresPrepared(GmresPrepared&&) noexcept = default;
GmresPrepared& GmresPrepared::operator=(GmresPrepared&&) noexcept = default;
GmresPrepared::~GmresPrepared() = default;
GmresResult<GmresPrepared> GmresPrepared::prepare(const GmresOperator& op,
                                                  const FixedPreconditioner& p,
                                                  const Scales& scales, const GmresLimits& limits,
                                                  core::ResourcePlan& plan) noexcept {
    GmresDiagnostics stats{};
    const auto context = op.context();
    return detail::GmresFactory::guard<GmresPrepared>(
        context, stats, [&]() -> Result<GmresPrepared> {
            const auto valid_limits = check_limits(limits, context);
            if (!valid_limits)
                return std::unexpected(valid_limits.error());
            const auto n = op.order();
            if (n < 0 || p.order() != n)
                return fail(FailureCode::InvalidInput, context);
            const auto fixed = p.identity();
            if (fixed.target != context)
                return fail(FailureCode::StaleGeneration, context);
            auto valid = detail::check_scale_inputs(n, scales, context);
            if (!valid)
                return std::unexpected(internal_failure(valid.error(), context));
            valid = detail::check_capacity_inputs(n, limits.maximum_rhs_columns, 0,
                                                  shared_limits(limits));
            if (!valid)
                return std::unexpected(internal_failure(valid.error(), context));
            if (n > 0) {
                auto a = op.validate_frozen_state();
                if (!a)
                    return std::unexpected(a.error());
                auto b = p.validate_frozen_state();
                if (!b)
                    return std::unexpected(b.error());
                auto after = state_matches(op, p, context, fixed, n);
                if (!after)
                    return std::unexpected(after.error());
            }
            auto memory = admit_backend_memory(n == 0 ? detail::BackendPlanFactory::make(0, 0)
                                                      : detail::BackendPlanFactory::make(0),
                                               limits.backend_memory);
            if (!memory)
                return std::unexpected(internal_failure(memory.error(), context));
            auto data = std::make_unique<detail::GmresData>(context, fixed, limits, *memory);
            data->order = n;
            data->width = std::min(n, limits.restart);
            auto metadata = detail::SolutionFactory::prepare_metadata(n, 0, scales, {}, plan);
            if (!metadata)
                return std::unexpected(internal_failure(metadata.error(), context));
            data->scales = std::move(*metadata);
            if (n > 0) {
                if (data->width == std::numeric_limits<Index>::max())
                    return fail(FailureCode::ResourceExceeded, context);
                auto basis_size = product(n, data->width + 1),
                     h_size = product(data->width + 1, data->width);
                if (!basis_size)
                    return std::unexpected(internal_failure(basis_size.error(), context));
                if (!h_size)
                    return std::unexpected(internal_failure(h_size.error(), context));
                for (auto* buffer :
                     {&data->residual, &data->residual_error, &data->denominator, &data->rhs_hi,
                      &data->rhs_lo, &data->word_hi, &data->word_lo, &data->tail_hi, &data->tail_lo,
                      &data->error_hi, &data->error_lo, &data->weights, &data->work, &data->temp,
                      &data->delta, &data->proposal}) {
                    auto allocation = plan.allocate<double>(n, operation);
                    if (!allocation)
                        return std::unexpected(internal_failure(allocation.error(), context));
                    *buffer = std::move(*allocation);
                }
                auto basis = plan.allocate<double>(*basis_size, operation),
                     h = plan.allocate<double>(*h_size, operation),
                     g = plan.allocate<double>(data->width + 1, operation);
                if (!basis)
                    return std::unexpected(internal_failure(basis.error(), context));
                if (!h)
                    return std::unexpected(internal_failure(h.error(), context));
                if (!g)
                    return std::unexpected(internal_failure(g.error(), context));
                data->basis = std::move(*basis);
                data->hessenberg = std::move(*h);
                data->g = std::move(*g);
                for (auto* buffer : {&data->cosines, &data->sines, &data->y}) {
                    auto allocation = plan.allocate<double>(data->width, operation);
                    if (!allocation)
                        return std::unexpected(internal_failure(allocation.error(), context));
                    *buffer = std::move(*allocation);
                }
            }
            return GmresPrepared{std::move(data)};
        });
}
} // namespace uvlas::numerics::linear
namespace uvlas::numerics::linear {
namespace {
Result<ColumnMetrics> evaluate(Runner& run, std::span<const double> xhi,
                               std::span<const double> xlo, const Scales& scales,
                               const Accuracy& accuracy) noexcept {
    auto& d = run.d;
    const auto n = static_cast<std::size_t>(d.order);
    bool use_hi = false, use_lo = false, use_abs = false;
    for (std::size_t i = 0; i < n; ++i) {
        use_hi |= xhi[i] != 0;
        use_lo |= xlo[i] != 0;
        auto magnitude = detail::actual_magnitude_bounds(xhi[i], xlo[i]);
        if (!magnitude)
            return std::unexpected(magnitude.error());
        d.weights.values()[i] = magnitude->lower;
        use_abs |= magnitude->lower != 0;
    }
    auto admitted =
        run.admit(static_cast<std::uint64_t>(use_hi) + static_cast<std::uint64_t>(use_lo) +
                  static_cast<std::uint64_t>(use_abs));
    if (!admitted)
        return std::unexpected(admitted.error());
    for (auto* b :
         {&d.word_hi, &d.word_lo, &d.tail_hi, &d.tail_lo, &d.error_hi, &d.error_lo, &d.denominator})
        std::fill(b->values().begin(), b->values().end(), 0.0);
    if (use_hi) {
        auto r = run.action(xhi, d.word_hi.values(), d.tail_hi.values(), d.error_hi.values());
        if (!r)
            return std::unexpected(r.error());
    }
    if (use_lo) {
        auto r = run.action(xlo, d.word_lo.values(), d.tail_lo.values(), d.error_lo.values());
        if (!r)
            return std::unexpected(r.error());
    }
    if (use_abs) {
        auto r = run.absolute(d.weights.values(), d.denominator.values());
        if (!r)
            return std::unexpected(r.error());
    }
    run.qualification();
    auto combined = detail::combine_gmres_residual(
        d.order, d.rhs_hi.values(), d.rhs_lo.values(), d.word_hi.values(), d.word_lo.values(),
        d.tail_hi.values(), d.tail_lo.values(), d.error_hi.values(), d.error_lo.values(),
        d.residual.values(), d.residual_error.values());
    if (!combined)
        return std::unexpected(combined.error());
    for (std::size_t i = 0; i < n; ++i) {
        auto b = detail::actual_magnitude_bounds(d.rhs_hi.values()[i], d.rhs_lo.values()[i]);
        if (!b)
            return std::unexpected(b.error());
        d.denominator.values()[i] =
            detail::add_magnitudes({d.denominator.values()[i], d.denominator.values()[i]}, *b)
                .lower;
    }
    return detail::evaluate_column_metrics<double>(
        d.order, d.residual.values(), d.residual_error.values(), d.denominator.values(),
        d.rhs_hi.values(), d.rhs_lo.values(), scales.residual, accuracy);
}

// Left preconditioning acts in physical coordinates; the basis uses the declared
// unknown scales. Hessenberg residuals generate candidates, never certificates.
Result<void> cycle(Runner& run, std::span<double> hi, std::span<double> lo, const Scales& scales,
                   bool& two, RetentionPath& path, GmresDiagnostics& diagnostics) noexcept {
    auto& d = run.d;
    const auto n = static_cast<std::size_t>(d.order);
    const auto m = static_cast<std::size_t>(d.width), ld = m + 1;
    auto available = run.admit(1);
    if (!available)
        return available;
    auto p = run.precondition(d.residual.values(), d.temp.values());
    if (!p)
        return p;
    for (std::size_t i = 0; i < n; ++i)
        d.work.values()[i] = d.temp.values()[i] / scales.unknown[i];
    auto beta = norm(d.work.values(), d.context);
    if (!beta)
        return std::unexpected(beta.error());
    if (*beta == 0) {
        diagnostics.stop = GmresStopReason::KrylovBreakdown;
        return fail(FailureCode::LinearUnconverged, d.context);
    }
    std::fill(d.hessenberg.values().begin(), d.hessenberg.values().end(), 0.0);
    std::fill(d.g.values().begin(), d.g.values().end(), 0.0);
    d.g.values()[0] = *beta;
    for (std::size_t i = 0; i < n; ++i)
        d.basis.values()[i] = d.work.values()[i] / *beta;
    std::size_t steps = 0;
    bool breakdown = false;
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i < n; ++i)
            d.work.values()[i] = d.basis.values()[j * n + i] * scales.unknown[i];
        auto a = run.action(d.work.values(), d.word_hi.values(), d.tail_hi.values(),
                            d.error_hi.values());
        if (!a)
            return a;
        for (std::size_t i = 0; i < n; ++i)
            d.word_hi.values()[i] += d.tail_hi.values()[i];
        auto pre = run.precondition(d.word_hi.values(), d.temp.values());
        if (!pre)
            return pre;
        for (std::size_t i = 0; i < n; ++i)
            d.work.values()[i] = d.temp.values()[i] / scales.unknown[i];
        // Two deterministic MGS passes reduce loss of basis orthogonality.
        for (unsigned pass = 0; pass < 2; ++pass)
            for (std::size_t k = 0; k <= j; ++k) {
                auto q = std::span<const double>{d.basis.values().subspan(k * n, n)};
                auto h = dot(q, d.work.values(), d.context);
                if (!h)
                    return std::unexpected(h.error());
                d.hessenberg.values()[j * ld + k] += *h;
                for (std::size_t i = 0; i < n; ++i)
                    d.work.values()[i] -= *h * q[i];
            }
        auto hnext = norm(d.work.values(), d.context);
        if (!hnext)
            return std::unexpected(hnext.error());
        d.hessenberg.values()[j * ld + j + 1] = *hnext;
        breakdown = *hnext == 0;
        if (!breakdown)
            for (std::size_t i = 0; i < n; ++i)
                d.basis.values()[(j + 1) * n + i] = d.work.values()[i] / *hnext;
        for (std::size_t k = 0; k < j; ++k) {
            auto& upper = d.hessenberg.values()[j * ld + k];
            auto& lower = d.hessenberg.values()[j * ld + k + 1];
            const double u = d.cosines.values()[k] * upper + d.sines.values()[k] * lower;
            lower = -d.sines.values()[k] * upper + d.cosines.values()[k] * lower;
            upper = u;
        }
        auto& diagonal = d.hessenberg.values()[j * ld + j];
        auto& sub = d.hessenberg.values()[j * ld + j + 1];
        const double radius = std::hypot(diagonal, sub);
        if (!std::isfinite(radius))
            return fail(FailureCode::PrecisionUnresolved, d.context);
        if (radius == 0) {
            diagnostics.stop = GmresStopReason::KrylovBreakdown;
            return fail(FailureCode::LinearUnconverged, d.context);
        }
        d.cosines.values()[j] = diagonal / radius;
        d.sines.values()[j] = sub / radius;
        diagonal = radius;
        sub = 0;
        d.g.values()[j + 1] = -d.sines.values()[j] * d.g.values()[j];
        d.g.values()[j] *= d.cosines.values()[j];
        ++diagnostics.arnoldi_steps;
        steps = j + 1;
        if (breakdown)
            break;
    }
    for (std::size_t ii = steps; ii-- > 0;) {
        double value = d.g.values()[ii];
        for (std::size_t j = ii + 1; j < steps; ++j)
            value -= d.hessenberg.values()[j * ld + ii] * d.y.values()[j];
        const double pivot = d.hessenberg.values()[ii * ld + ii];
        if (pivot == 0)
            return fail(FailureCode::LinearUnconverged, d.context);
        d.y.values()[ii] = value / pivot;
        if (!std::isfinite(d.y.values()[ii]))
            return fail(FailureCode::PrecisionUnresolved, d.context);
    }
    bool unchanged = true;
    for (std::size_t i = 0; i < n; ++i) {
        uvlas::numerics::detail::CompensatedAccumulator sum;
        for (std::size_t j = 0; j < steps; ++j) {
            auto r = sum.product(d.basis.values()[j * n + i], d.y.values()[j], 1);
            if (!r)
                return r;
        }
        auto value = sum.result();
        if (!value)
            return std::unexpected(value.error());
        d.delta.values()[i] = value->value * scales.unknown[i];
        d.proposal.values()[i] = hi[i] + d.delta.values()[i];
        if (!std::isfinite(d.proposal.values()[i]))
            return fail(FailureCode::PrecisionUnresolved, d.context);
        unchanged &= std::bit_cast<std::uint64_t>(d.proposal.values()[i]) ==
                     std::bit_cast<std::uint64_t>(hi[i]);
    }
    if (!two && unchanged) {
        two = true;
        path = RetentionPath::PromotedTwoWord;
    }
    bool words_unchanged = true;
    for (std::size_t i = 0; i < n; ++i) {
        const auto oldhi = std::bit_cast<std::uint64_t>(hi[i]),
                   oldlo = std::bit_cast<std::uint64_t>(lo[i]);
        if (two) {
            auto r = uvlas::numerics::detail::apply_two_word_correction(hi[i], lo[i],
                                                                        d.delta.values()[i]);
            if (!r)
                return r;
        } else
            hi[i] = d.proposal.values()[i];
        words_unchanged &= oldhi == std::bit_cast<std::uint64_t>(hi[i]) &&
                           oldlo == std::bit_cast<std::uint64_t>(lo[i]);
    }
    // Exact unchanged words imply that repeating this complete fixed cycle cannot advance.
    if (words_unchanged) {
        diagnostics.stop = GmresStopReason::Stagnated;
        return fail(FailureCode::LinearUnconverged, d.context);
    }

    return {};
}
} // namespace
} // namespace uvlas::numerics::linear
namespace uvlas::numerics::linear {
GmresResult<GmresSolution> solve_linear(const GmresOperator& op, const FixedPreconditioner& p,
                                        const Rhs<double>& rhs, const Scales& scales,
                                        const Precision<double>& precision,
                                        const Accuracy& accuracy, const GmresLimits& limits,
                                        GmresPrepared& prepared,
                                        core::ResourcePlan& plan) noexcept {
    GmresDiagnostics diagnostics;
    return detail::GmresFactory::guard<GmresSolution>(
        op.context(), diagnostics, [&]() -> Result<GmresSolution> {
            const auto context = op.context();
            if (!prepared.data_)
                return fail(FailureCode::InvalidInput, context);
            auto& d = *prepared.data_;
            const Index n = d.order, columns = rhs.values().shape().columns;
            auto checked = check_limits(limits, context);
            if (!checked)
                return std::unexpected(checked.error());
            auto state = state_matches(op, p, d.context, d.preconditioner, n);
            if (!state)
                return std::unexpected(state.error());
            auto valid = detail::check_solve_inputs(n, 0, d.context, rhs, scales, precision,
                                                    accuracy, shared_limits(limits));
            if (!valid)
                return std::unexpected(valid.error());
            auto capacity = detail::check_capacity_inputs(n, columns, 0, shared_limits(d.capacity));
            if (!capacity)
                return std::unexpected(capacity.error());
            if (limits.restart > d.capacity.restart || !same_scales(scales, d.scales))
                return fail(FailureCode::StaleGeneration, context);
            auto admitted = admit_backend_memory(d.admission.plan(), limits.backend_memory);
            if (!admitted)
                return std::unexpected(admitted.error());
            auto lease = d.lease.acquire();
            if (!lease)
                return std::unexpected(lease.error());
            if (n > 0 && columns > 0) {
                auto a = op.validate_frozen_state();
                if (!a)
                    return std::unexpected(a.error());
                auto pc = p.validate_frozen_state();
                if (!pc)
                    return std::unexpected(pc.error());
                state = state_matches(op, p, d.context, d.preconditioner, n);
                if (!state)
                    return std::unexpected(state.error());
            }
            auto metadata = detail::SolutionFactory::prepare_metadata(
                n, columns, scales, precision.required_checks, plan);
            if (!metadata)
                return std::unexpected(metadata.error());
            auto hi =
                DenseReal::copy(rhs.values().shape(), high(rhs.values()).values(), plan, operation);
            if (!hi)
                return std::unexpected(hi.error());
            auto lo =
                DenseReal::copy(rhs.values().shape(), high(rhs.values()).values(), plan, operation);
            if (!lo)
                return std::unexpected(lo.error());
            std::fill(lo->values().begin(), lo->values().end(), 0.0);
            const auto* projection = std::get_if<Binary64Projection>(&precision.returned);
            const auto* retained_request = std::get_if<RetainedPrecision>(&precision.returned);
            if (!projection && !retained_request)
                return fail(FailureCode::PrecisionUnsupported, context);
            std::optional<DenseReal> projected;
            if (projection) {
                auto value = DenseReal::copy(rhs.values().shape(), high(rhs.values()).values(),
                                             plan, operation);
                if (!value)
                    return std::unexpected(value.error());
                projected.emplace(std::move(*value));
            }
            bool two = rhs.values().precision() == core::Precision::TwoWord;
            RetentionPath path = two ? RetentionPath::InitiallyTwoWord : RetentionPath::Binary64;
            const auto* rhs_low = low(rhs.values());
            Runner runner{d, op, p, limits, diagnostics};
            const auto size = static_cast<std::size_t>(n);
            // A call may request a smaller restart width without reallocating preparation.
            const Index prepared_width = d.width;
            d.width = std::min(n, limits.restart);
            struct Restore {
                Index& target;
                Index value;
                ~Restore() {
                    target = value;
                }
            };
            const Restore restore{d.width, prepared_width};
            for (Index c = 0; c < columns; ++c) {
                const auto offset = static_cast<std::size_t>(c * hi->shape().leading_dimension);
                auto xhi = hi->values().subspan(offset, size),
                     xlo = lo->values().subspan(offset, size);
                std::fill(xhi.begin(), xhi.end(), 0.0);
                for (std::size_t i = 0; i < size; ++i) {
                    const auto input =
                        static_cast<std::size_t>(c * rhs.values().shape().leading_dimension) + i;
                    d.rhs_hi.values()[i] = high(rhs.values()).values()[input];
                    d.rhs_lo.values()[i] = rhs_low ? rhs_low->values()[input] : 0;
                }
                auto metrics = evaluate(runner, xhi, xlo, scales, accuracy);
                if (!metrics)
                    return std::unexpected(metrics.error());
                bool first_cycle = true;
                while (!detail::meets_accuracy(*metrics, accuracy)) {
                    if (!first_cycle) {
                        auto available = runner.admit(1);
                        if (!available)
                            return std::unexpected(available.error());
                        if (diagnostics.restarts == std::numeric_limits<std::uint64_t>::max())
                            return fail(FailureCode::ResourceExceeded, context);
                        ++diagnostics.restarts;
                    }
                    first_cycle = false;
                    auto step = cycle(runner, xhi, xlo, scales, two, path, diagnostics);
                    if (!step)
                        return std::unexpected(step.error());
                    metrics = evaluate(runner, xhi, xlo, scales, accuracy);
                    if (!metrics)
                        return std::unexpected(metrics.error());
                }
                metadata->retained_columns.values()[static_cast<std::size_t>(c)] = *metrics;
                if (projected) {
                    auto target = projected->values().subspan(offset, size);
                    for (std::size_t i = 0; i < size; ++i) {
                        target[i] = xhi[i] + xlo[i];
                        if (!std::isfinite(target[i])) {
                            diagnostics.stop = GmresStopReason::QualificationFailure;
                            return fail(FailureCode::PrecisionUnresolved, context);
                        }
                    }
                    std::fill(d.proposal.values().begin(), d.proposal.values().end(), 0.0);
                    auto projected_metrics =
                        evaluate(runner, target, d.proposal.values(), scales, accuracy);
                    if (!projected_metrics) {
                        if (projected_metrics.error().code == FailureCode::LinearUnconverged &&
                            diagnostics.stop == GmresStopReason::OperatorLimit) {
                            auto failure = projected_metrics.error();
                            failure.code = FailureCode::PrecisionUnresolved;
                            return std::unexpected(failure);
                        }
                        return std::unexpected(projected_metrics.error());
                    }
                    if (!detail::meets_accuracy(*projected_metrics, accuracy)) {
                        diagnostics.stop = GmresStopReason::QualificationFailure;
                        return fail(FailureCode::PrecisionUnresolved, context);
                    }
                    metadata->returned_columns.values()[static_cast<std::size_t>(c)] =
                        *projected_metrics;
                }
            }
            const auto generation =
                projection ? projection->retained_generation : retained_request->generation;
            auto retained =
                two ? PrecisionBuffer<double>::two_word(std::move(*hi), std::move(*lo), generation)
                    : PrecisionBuffer<double>::binary64(std::move(*hi), generation);
            if (!retained)
                return std::unexpected(retained.error());
            std::optional<PrecisionBuffer<double>> projected_buffer;
            if (projected) {
                auto value = PrecisionBuffer<double>::binary64(std::move(*projected),
                                                               projection->target_generation);
                if (!value)
                    return std::unexpected(value.error());
                projected_buffer.emplace(std::move(*value));
            }
            diagnostics.stop = GmresStopReason::QualificationFailure;
            auto result = detail::SolutionFactory::finalize(
                std::move(*retained), std::move(projected_buffer), precision, std::move(*metadata),
                detail::Finalization{d.context, rhs.cache_request().rhs, rhs.values().generation(),
                                     accuracy, d.admission, *admitted, 0, path, NotEvaluated{}});
            if (!result)
                return std::unexpected(result.error());
            auto final_state = runner.stable();
            if (!final_state)
                return std::unexpected(final_state.error());
            diagnostics.stop = GmresStopReason::Qualified;
            return detail::GmresFactory::success(std::move(*result), diagnostics);
        });
}
} // namespace uvlas::numerics::linear
