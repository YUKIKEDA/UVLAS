#include "pchip_curve.hpp"

#include "pchip_curve_plan.hpp"
#include "pchip_domain_route.hpp"
#include "pchip_knot_locator.hpp"
#include "pchip_slot_private.hpp"

#include <algorithm>
#include <array>

namespace uvlas::aero::detail::pchip {
namespace curve_detail {
using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::operation;
using slot_detail::supplement;

struct CurveOwner {
    CurveOwner(const SharedPreparationInput& original, CoefficientBlock values,
               std::uint64_t steps) noexcept
        : input(original), coefficients(std::move(values)), work(steps) {}

    SharedPreparationInput input;
    CoefficientBlock coefficients;
    std::uint64_t work;
};

using PendingCurve = std::optional<CurveOwner>;

struct CurveHolder {
    explicit CurveHolder(core::PlannedBuffer<PendingCurve> buffer) noexcept
        : storage(std::move(buffer)) {}
    core::PlannedBuffer<PendingCurve> storage;
};

struct JetOwner {
    std::shared_ptr<const CurveOwner> parent;
    core::PlannedBuffer<RationalView> records;
    core::PlannedBuffer<std::uint32_t> limbs;
    PchipJetLocation location{std::nullopt, RequestedSide::Automatic, JetCandidate::Only};
    PchipInputStamp stamp{PchipIssuerId{0}, PchipInputGeneration{0}};
    JetOrder order = JetOrder::Value;
    std::uint64_t components = 0, work = 0;
};

constexpr std::array<std::uint32_t, 1> one{1};
constexpr RationalView zero{0, {}, one};
constexpr std::uint64_t wrapper_work = 8;
constexpr std::uint64_t allocator_estimate = 4 * sizeof(void*);
constexpr std::uint64_t curve_metadata = sizeof(CurveHolder) + 2 * allocator_estimate;
constexpr std::uint64_t preparation_arrays =
    7; // Input header, records, limbs, dimensions, arena three.
constexpr std::uint64_t numerical_arrays =
    10; // Coefficient arena three plus four arrays; C1 arena three.

struct Meter {
    std::uint64_t cap = 0;
    const std::stop_token& stop;
    core::FailureDetails details;
    std::uint64_t used = 0;

    [[nodiscard]] std::uint64_t remaining() const noexcept {
        return cap - used;
    }

    [[nodiscard]] core::Result<void> visit(std::uint64_t steps = 1) noexcept {
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        if (steps > remaining())
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        used += steps;
        return {};
    }
};

[[nodiscard]] core::Result<void> admit_metadata(numerics::linear::BackendMemoryLimit limit,
                                                std::uint64_t bytes,
                                                core::FailureDetails details) noexcept {
    const auto selection = limit.selection();
    const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&selection);
    if (estimate == nullptr)
        return std::unexpected(failure(core::FailureCode::Unsupported, details));
    if (estimate->bytes < bytes)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    return {};
}

[[nodiscard]] core::Result<void> admit_payload(const CurveAdmission& forecast,
                                               const core::ResourcePlan& memory,
                                               core::FailureDetails details) noexcept {
    if (forecast.payload_bytes > memory.remaining_bytes()) {
        const auto already = memory.byte_limit() - memory.remaining_bytes();
        CheckedCount checked;
        const auto needed = checked.add(already, forecast.payload_bytes);
        if (!details.resources && checked.ok)
            details.resources = core::ResourceUsage{core::ResourceKind::Bytes, needed, already};
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    }
    return {};
}

template <class T>
[[nodiscard]] core::Result<core::PlannedBuffer<T>>
allocate(core::ResourcePlan& memory, std::uint64_t count, core::FailureDetails details) noexcept {
    if (count > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    auto result = memory.allocate<T>(static_cast<core::Index>(count), operation);
    if (!result)
        return std::unexpected(supplement(result.error(), details));
    return result;
}

[[nodiscard]] core::Result<BitBound> shape(RationalView value, std::uint64_t limit,
                                           core::FailureDetails details) noexcept {
    const auto words = limit / 32 + (limit % 32 == 0 ? 0U : 1U);
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
    const BitBound result{bits(value.numerator), bits(value.denominator)};
    if (!checked.ok || std::max(result.numerator, result.denominator) > limit)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    return result;
}

[[nodiscard]] BitBound cover(BitBound left, BitBound right) noexcept {
    return {std::max(left.numerator, right.numerator),
            std::max(left.denominator, right.denominator)};
}

[[nodiscard]] core::Result<PreparationRequest> input_view(const SharedPreparationInput& owner,
                                                          core::FailureDetails details) noexcept {
    const auto source = owner.input();
    if (!source)
        return std::unexpected(supplement(source.error(), details));
    const auto view = source->get().view();
    if (!view)
        return std::unexpected(supplement(view.error(), details));
    return view;
}

[[nodiscard]] std::uint64_t coefficient_count(const CoefficientBlock& coefficients,
                                              CheckedCount& checked) noexcept {
    return checked.mul(checked.mul(coefficients.nodes() - 1, coefficients.components()), 4);
}

[[nodiscard]] core::Result<RationalView> flat_coefficient(const CoefficientBlock& coefficients,
                                                          std::uint64_t flat) noexcept {
    return coefficients.coefficient(flat / 4 / coefficients.components(),
                                    (flat / 4) % coefficients.components(), flat % 4);
}

struct PreparePlan {
    CurveAdmission all;
    CurveNumericsPlan numerics;
    std::optional<SharedPreparationPlan> shared;
    std::optional<PreparationPlan> input;
    std::uint64_t preflight = 0;
};

struct QualificationStorage {
    std::shared_ptr<CurveHolder> holder;
    CoefficientStorage coefficients;
    Workspace boundary;
};

using DomainForecast = std::variant<DomainPlan, PeriodicTimePlan>;
using DomainStorage = std::variant<Workspace, PeriodicStorage>;

struct EvaluationPlan {
    CurveAdmission all;
    RationalLimits derived;
    DomainForecast domain;
    LocatorPlan locator;
    IntervalPlan interval;
    std::uint64_t records = 0, limbs = 0, agreement = 0, preflight = 0;
    std::size_t candidates = 0;
};

struct EvaluationStorage {
    DomainStorage domain;
    Workspace locator;
    std::array<std::optional<IntervalStorage>, 2> interval;
};

struct RestorePlan {
    PreparePlan prepare;
    std::uint64_t preflight = 0, comparison = 0, limbs = 0;
};

struct Route {
    RationalView normalized;
    std::optional<PeriodicTime> periodic;
    std::optional<IntervalRequest> left, right;
    std::optional<std::uint64_t> constant;
    SideLocation location = SideLocation::Single;
    std::uint64_t steps = 0;
};

struct Selected {
    std::optional<IntervalValues> values;
    std::optional<std::uint64_t> constant;
    PchipJetLocation location{std::nullopt, RequestedSide::Automatic, JetCandidate::Only};
    std::uint64_t steps = 0;
};

struct Builder {
    struct Lease {
        CurveWorkspace& workspace;
        ~Lease() {
            workspace.busy_.clear(std::memory_order_release);
        }
    };

    static bool acquire(CurveWorkspace& workspace) noexcept {
        return !workspace.busy_.test_and_set(std::memory_order_acquire);
    }

    static core::Result<PreparePlan> plan_prepare(const PreparationRequest& input,
                                                  CurveLimits limits, const std::stop_token& stop,
                                                  core::FailureDetails details,
                                                  bool copy_input = true) noexcept {
        std::uint64_t input_bytes = 0, input_steps = 0, input_metadata = 0;
        std::uint64_t input_preflight = 0, input_bits = 0;
        std::optional<SharedPreparationPlan> shared_plan;
        std::optional<PreparationPlan> input_plan;
        auto numerical_limits = limits.arithmetic;
        if (copy_input) {
            const auto source = forecast_preparation_input(input, limits.arithmetic, stop, details);
            if (!source)
                return std::unexpected(source.error());
            input_plan.emplace(*source);
            const auto shared = forecast_shared_preparation(
                {limits.metadata, numerical_limits.steps - source->reserved_steps()}, stop,
                details);
            if (!shared)
                return std::unexpected(shared.error());
            shared_plan.emplace(*shared);
            CheckedCount checked;
            input_bytes = checked.add(source->payload_bytes(), shared->payload_bytes());
            input_preflight = source->preflight_steps();
            input_bits = source->arithmetic().integer_bits();
            input_steps = checked.add(source->reserved_steps(), shared->reserved_steps());
            input_metadata =
                checked.add(preparation_arrays * allocator_estimate,
                            checked.add(shared->metadata_estimate_bytes(), allocator_estimate));
            if (!checked.ok)
                return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
            numerical_limits.steps -= input_steps;
        }
        const auto numerics = forecast_curve_numerics(input, numerical_limits, stop, details);
        if (!numerics)
            return std::unexpected(numerics.error());
        CheckedCount checked;
        const CurveAdmission result{
            checked.add(checked.add(input_bytes, numerics->payload_bytes()), sizeof(PendingCurve)),
            checked.add(input_metadata, numerical_arrays * allocator_estimate + curve_metadata),
            checked.add(checked.add(input_steps, numerics->reserved_steps()), wrapper_work),
            std::max({input_bits, numerics->coefficients().arithmetic().integer_bits(),
                      numerics->boundary().arithmetic().integer_bits()})};
        const auto preflight = checked.add(input_preflight, numerics->preflight_steps());
        if (!checked.ok || result.reserved_steps > limits.arithmetic.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        const auto metadata =
            admit_metadata(limits.metadata, result.metadata_estimate_bytes, details);
        if (!metadata)
            return std::unexpected(metadata.error());
        return PreparePlan{result, *numerics, shared_plan, input_plan, preflight};
    }

    static core::Result<QualificationStorage>
    provision_qualification(const PreparePlan& plan, core::ResourcePlan& memory,
                            const std::stop_token& stop, core::FailureDetails details) noexcept {
        auto storage = allocate<PendingCurve>(memory, 1, details);
        if (!storage)
            return std::unexpected(storage.error());
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        auto holder =
            core::guard_entry(operation, core::FailureCode::InvalidInput, details,
                              [&]() -> core::Result<std::shared_ptr<CurveHolder>> {
                                  return std::make_shared<CurveHolder>(std::move(*storage));
                              });
        if (!holder)
            return std::unexpected(holder.error());
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        auto coefficients =
            prepare_coefficient_storage(plan.numerics.coefficients(), memory, stop, details);
        if (!coefficients)
            return std::unexpected(coefficients.error());
        auto boundary =
            Workspace::create(plan.numerics.boundary().arithmetic(), memory, details, stop);
        if (!boundary)
            return std::unexpected(boundary.error());
        return QualificationStorage{std::move(*holder), std::move(*coefficients),
                                    std::move(*boundary)};
    }

    static core::Result<PchipCurveReady>
    qualify(const SharedPreparationInput& source, const PreparePlan& plan,
            QualificationStorage storage, CurveLimits limits, const std::stop_token& stop,
            core::FailureDetails details, std::uint64_t earlier_work,
            std::span<const RationalView> saved = {}, std::uint64_t comparison = 0) noexcept {
        const auto input = input_view(source, details);
        if (!input)
            return std::unexpected(input.error());
        auto coefficients = prepare_staged_coefficients(source, plan.numerics.coefficients(),
                                                        std::move(storage.coefficients),
                                                        limits.arithmetic, stop, details);
        if (!coefficients)
            return std::unexpected(coefficients.error());
        const auto boundary = check_prepared_boundary(
            *coefficients, input->boundary, plan.numerics.boundary(), storage.boundary,
            plan.numerics.boundary_limits(), stop, details);
        if (!boundary)
            return std::unexpected(boundary.error());
        Meter meter{comparison, stop, details};
        if (comparison != 0) {
            for (std::size_t i = 0; i < saved.size(); ++i) {
                const auto expected = flat_coefficient(*coefficients, i);
                if (!expected)
                    return std::unexpected(supplement(expected.error(), details));
                const auto compared =
                    equal_record(*expected, saved[i], meter, core::FailureCode::InvalidInput);
                if (!compared)
                    return std::unexpected(compared.error());
            }
        }
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        const auto used = earlier_work + plan.preflight + coefficients->used_steps() +
                          boundary->steps + meter.used + wrapper_work;
        const auto& result =
            storage.holder->storage.values()[0].emplace(source, std::move(*coefficients), used);
        return PchipCurveReady{
            std::shared_ptr<const CurveOwner>{std::move(storage.holder), &result}};
    }

    static core::Result<PchipCurveReady> prepare(const PreparationRequest& input,
                                                 CurveWorkspace& workspace, CurveLimits limits,
                                                 core::ResourcePlan& memory,
                                                 const std::stop_token& stop,
                                                 core::FailureDetails details) noexcept {
        if (!acquire(workspace))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto plan = plan_prepare(input, limits, stop, details);
        if (!plan)
            return std::unexpected(plan.error());
        const auto admitted = admit_payload(plan->all, memory, details);
        if (!admitted)
            return std::unexpected(admitted.error());
        if (!plan->shared || !plan->input)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        auto shared = prepare_shared_preparation_storage(*plan->shared, memory, stop, details);
        if (!shared)
            return std::unexpected(shared.error());
        auto storage = provision_qualification(*plan, memory, stop, details);
        if (!storage)
            return std::unexpected(storage.error());
        // Input preparation performs the remaining allocations, then owns and checks its snapshot.
        auto original = make_planned_preparation_input(input, *plan->input, limits.arithmetic,
                                                       memory, stop, details);
        if (!original)
            return std::unexpected(original.error());
        const auto input_work = original->used_steps();
        auto source =
            publish_shared_preparation(std::move(*original), std::move(*shared), stop, details);
        if (!source)
            return std::unexpected(source.error());
        return qualify(*source, *plan, std::move(*storage), limits, stop, details,
                       input_work + plan->shared->reserved_steps());
    }

    static core::Result<void> equal_record(RationalView expected, RationalView actual, Meter& meter,
                                           core::FailureCode mismatch) noexcept {
        const auto ready = meter.visit();
        if (!ready)
            return std::unexpected(ready.error());
        if (expected.sign != actual.sign || expected.numerator.size() != actual.numerator.size() ||
            expected.denominator.size() != actual.denominator.size())
            return std::unexpected(failure(mismatch, meter.details));
        const auto equal_words = [&](auto a, auto b) -> core::Result<void> {
            for (std::size_t i = 0; i < a.size(); ++i) {
                const auto next = meter.visit();
                if (!next)
                    return std::unexpected(next.error());
                if (a[i] != b[i])
                    return std::unexpected(failure(mismatch, meter.details));
            }
            return {};
        };
        const auto numerator = equal_words(expected.numerator, actual.numerator);
        if (!numerator)
            return std::unexpected(numerator.error());
        return equal_words(expected.denominator, actual.denominator);
    }

    static core::Result<std::uint64_t> match_expected(const CurveOwner& owner,
                                                      const SharedPreparationInput& expected,
                                                      CurveLimits limits,
                                                      const std::stop_token& stop,
                                                      core::FailureDetails details) noexcept {
        if (!expected.active())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        if (owner.input.shares_input_with(expected)) {
            Meter meter{limits.arithmetic.steps, stop, details};
            const auto found = meter.visit();
            if (!found)
                return std::unexpected(found.error());
            return meter.used;
        }
        const auto actual_input = owner.input.input();
        const auto expected_input = expected.input();
        if (!actual_input || !expected_input)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        return match_preparation_input(actual_input->get(), expected_input->get(),
                                       limits.arithmetic.steps, stop, details);
    }

    static core::Result<EvaluationPlan> plan_evaluation(const PchipCurveReady& curve,
                                                        const SharedPreparationInput& expected,
                                                        PchipQuery query, CurveLimits limits,
                                                        const std::stop_token& stop,
                                                        core::FailureDetails details) noexcept {
        if (!curve.active())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto side = decide_side(SideLocation::Single, query.side, query.order, details);
        if (!side)
            return std::unexpected(side.error());
        const auto time = shape(query.time, limits.arithmetic.input_bits, details);
        if (!time)
            return std::unexpected(time.error());
        const auto matched = match_expected(*curve.owner_, expected, limits, stop, details);
        if (!matched)
            return std::unexpected(matched.error());
        const auto input = curve.input();
        if (!input)
            return std::unexpected(supplement(input.error(), details));
        const auto& coefficients = curve.owner_->coefficients;
        auto remaining = limits.arithmetic;
        remaining.steps -= *matched;
        auto normalized = *time;
        auto interval_time = *time;
        DomainForecast domain_plan;
        std::uint64_t domain_bytes = 0, domain_steps = 0, domain_preflight = 0, domain_bits = 0;
        if (const auto* periodic = std::get_if<PeriodicBoundary>(&input->boundary.policy)) {
            const auto domain = forecast_periodic_time(coefficients, query.time, periodic->period,
                                                       remaining, stop, details);
            if (!domain)
                return std::unexpected(domain.error());
            normalized = domain->normalized_bound();
            const auto endpoint = coefficients.time(coefficients.nodes() - 1);
            if (!endpoint)
                return std::unexpected(supplement(endpoint.error(), details));
            const auto end = shape(*endpoint, limits.arithmetic.input_bits, details);
            if (!end)
                return std::unexpected(end.error());
            interval_time = cover(normalized, *end); // The left seam uses its own endpoint time.
            domain_bytes = domain->payload_bytes();
            domain_steps = domain->reserved_steps();
            domain_preflight = domain->preflight_steps();
            domain_bits = domain->arithmetic().integer_bits();
            domain_plan = *domain;
        } else {
            const auto domain = forecast_domain(
                coefficients, query.time, domain_policy(input->boundary), remaining, stop, details);
            if (!domain)
                return std::unexpected(domain.error());
            domain_bytes = domain->payload_bytes();
            domain_steps = domain->reserved_steps();
            domain_preflight = domain->preflight_steps();
            domain_bits = domain->arithmetic().integer_bits();
            domain_plan = *domain;
        }
        auto derived = limits.arithmetic;
        derived.input_bits = std::max(derived.input_bits, derived.work_bits);
        remaining.input_bits = derived.input_bits;
        remaining.steps -= domain_steps;
        const auto locator =
            forecast_location_bound(coefficients, normalized, remaining, stop, details);
        if (!locator)
            return std::unexpected(locator.error());
        remaining.steps -= locator->reserved_steps();
        const auto interval = forecast_interval_bound(coefficients, interval_time, query.order,
                                                      remaining, stop, details);
        if (!interval)
            return std::unexpected(interval.error());

        CheckedCount count;
        const std::size_t candidates =
            query.order == JetOrder::Second && query.side == RequestedSide::Automatic ? 2 : 1;
        // The interval upper bound covers both words of every derivative output. Automatic
        // may compare two jets at a join or compare the inside derivatives to held zeros.
        const auto words = interval->arithmetic().words();
        const auto components = coefficients.components();
        const auto records =
            count.add(count.mul(components, 1 + static_cast<std::uint64_t>(query.order)), 2);
        const auto per_record = count.mul(words, 2);
        const auto interval_work = count.add(
            interval->preflight_steps(),
            count.mul(candidates, interval->reserved_steps() - interval->local_preflight_steps()));
        if (!count.ok || interval_work > remaining.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        std::uint64_t packed_words = words;
        Meter hold_scan{remaining.steps - interval_work, stop, details};
        if (std::holds_alternative<HoldBoundary>(input->boundary.policy)) {
            for (const auto node : {std::uint64_t{0}, coefficients.nodes() - 1}) {
                for (std::uint64_t c = 0; c < components; ++c) {
                    const auto entered = hold_scan.visit();
                    if (!entered)
                        return std::unexpected(entered.error());
                    const auto value = coefficients.value(node, c);
                    if (!value)
                        return std::unexpected(supplement(value.error(), details));
                    packed_words = std::max(
                        {packed_words, value->numerator.size(), value->denominator.size()});
                }
            }
        }
        const auto limbs = count.mul(records, count.mul(packed_words, 2));
        const auto agreement = count.add(count.mul(components, count.add(per_record, 2)), 1);
        const auto pack_steps = count.add(records, limbs);
        const auto structural = count.add(count.add(pack_steps, agreement), wrapper_work);
        CurveAdmission total;
        total.payload_bytes =
            count.add(count.add(domain_bytes, locator->payload_bytes()),
                      count.add(count.mul(candidates, interval->payload_bytes()),
                                count.add(sizeof(JetOwner),
                                          count.add(count.mul(records, sizeof(RationalView)),
                                                    count.mul(limbs, sizeof(std::uint32_t))))));
        // Locator arena three, each interval arena three plus four arrays, and jet three.
        const std::uint64_t route_allocations =
            std::holds_alternative<PeriodicTimePlan>(domain_plan) ? 4 : 3;
        total.metadata_estimate_bytes =
            (route_allocations + 3 + candidates * 7 + 3) * allocator_estimate;
        const auto planning =
            count.add(count.add(*matched, hold_scan.used),
                      count.add(domain_preflight, count.add(locator->preflight_steps(),
                                                            interval->preflight_steps())));
        // Staged execution reuses the plan, so only the initial planning pass is counted.
        const auto execution = count.add(
            count.add(domain_steps - domain_preflight,
                      locator->reserved_steps() - locator->preflight_steps()),
            count.mul(candidates, interval->reserved_steps() - interval->local_preflight_steps()));
        total.maximum_integer_bits = std::max({domain_bits, locator->arithmetic().integer_bits(),
                                               interval->arithmetic().integer_bits()});
        total.reserved_steps = count.add(count.add(execution, planning), structural);
        constexpr auto max_index =
            static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
        if (!count.ok || total.reserved_steps > limits.arithmetic.steps || records > max_index ||
            limbs > max_index)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        const auto metadata =
            admit_metadata(limits.metadata, total.metadata_estimate_bytes, details);
        if (!metadata)
            return std::unexpected(metadata.error());
        return EvaluationPlan{total,   derived, domain_plan, *locator, *interval,
                              records, limbs,   agreement,   planning, candidates};
    }

    static DomainPolicy domain_policy(const BoundaryRequest& boundary) noexcept {
        return std::holds_alternative<HoldBoundary>(boundary.policy) ? DomainPolicy::Hold
                                                                     : DomainPolicy::Reject;
    }

    static core::Result<DomainStorage> provision_domain(const DomainForecast& plan,
                                                        core::ResourcePlan& memory,
                                                        const std::stop_token& stop,
                                                        core::FailureDetails details) noexcept {
        if (const auto* periodic = std::get_if<PeriodicTimePlan>(&plan)) {
            auto storage = prepare_periodic_storage(*periodic, memory, stop, details);
            if (!storage)
                return std::unexpected(storage.error());
            return DomainStorage{std::in_place_type<PeriodicStorage>, std::move(*storage)};
        }
        const auto* domain = std::get_if<DomainPlan>(&plan);
        if (domain == nullptr)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        auto storage = prepare_domain_storage(*domain, memory, stop, details);
        if (!storage)
            return std::unexpected(storage.error());
        return DomainStorage{std::in_place_type<Workspace>, std::move(*storage)};
    }

    static core::Result<EvaluationStorage>
    provision_evaluation(const EvaluationPlan& plan, core::ResourcePlan& memory,
                         const std::stop_token& stop, core::FailureDetails details) noexcept {
        auto domain = provision_domain(plan.domain, memory, stop, details);
        if (!domain)
            return std::unexpected(domain.error());
        auto locator = prepare_locator_storage(plan.locator, memory, stop, details);
        if (!locator)
            return std::unexpected(locator.error());
        EvaluationStorage storage{std::move(*domain), std::move(*locator), {}};
        for (std::size_t i = 0; i < plan.candidates; ++i) {
            auto interval = prepare_interval_storage(plan.interval, memory, stop, details);
            if (!interval)
                return std::unexpected(interval.error());
            storage.interval[i].emplace(std::move(*interval));
        }
        return storage;
    }

    static core::Result<Route> route(const CoefficientBlock& coefficients,
                                     const BoundaryRequest& boundary, PchipQuery query,
                                     const EvaluationPlan& plan, EvaluationStorage& storage,
                                     const std::stop_token& stop,
                                     core::FailureDetails details) noexcept {
        Route result;
        result.normalized = query.time;
        if (const auto* period = std::get_if<PeriodicBoundary>(&boundary.policy)) {
            const auto* planned = std::get_if<PeriodicTimePlan>(&plan.domain);
            auto* prepared = std::get_if<PeriodicStorage>(&storage.domain);
            if (planned == nullptr || prepared == nullptr)
                return std::unexpected(failure(core::FailureCode::InvalidInput, details));
            auto normal =
                normalize_prepared_periodic_time(coefficients, query.time, period->period, *planned,
                                                 std::move(*prepared), stop, details);
            if (!normal)
                return std::unexpected(normal.error());
            result.periodic.emplace(std::move(*normal));
            const auto t = result.periodic->time();
            if (!t)
                return std::unexpected(supplement(t.error(), details));
            result.normalized = *t;
            result.steps = result.periodic->used_steps();
            if (result.periodic->is_seam()) {
                const auto last = coefficients.time(coefficients.nodes() - 1);
                if (!last)
                    return std::unexpected(supplement(last.error(), details));
                result.location = SideLocation::C1Join;
                result.left = IntervalRequest{coefficients.nodes() - 2, *last, query.order};
                result.right = IntervalRequest{0, result.normalized, query.order};
                return result;
            }
        } else {
            const auto* planned = std::get_if<DomainPlan>(&plan.domain);
            auto* prepared = std::get_if<Workspace>(&storage.domain);
            if (planned == nullptr || prepared == nullptr)
                return std::unexpected(failure(core::FailureCode::InvalidInput, details));
            const auto domain =
                route_prepared_domain(coefficients, query.time, domain_policy(boundary), *planned,
                                      std::move(*prepared), stop, details);
            if (!domain)
                return std::unexpected(domain.error());
            result.steps = domain->used_steps();
            if (domain->has_constant()) {
                const auto endpoint = domain->constant_node();
                if (!endpoint)
                    return std::unexpected(supplement(endpoint.error(), details));
                result.constant = *endpoint;
                if (domain->relation() == DomainRelation::AtFirst)
                    result.location = SideLocation::HoldAtFirst;
                if (domain->relation() == DomainRelation::AtLast)
                    result.location = SideLocation::HoldAtLast;
            }
            if (!domain->has_inside())
                return result; // Exterior Hold, the only candidate is its constant.
        }
        const auto location = locate_prepared_time(coefficients, result.normalized, plan.locator,
                                                   std::move(storage.locator), stop, details);
        if (!location)
            return std::unexpected(location.error());
        result.steps += location->used_steps();
        if (location->has_left()) {
            const auto interval = location->left_interval();
            if (!interval)
                return std::unexpected(supplement(interval.error(), details));
            result.left = IntervalRequest{*interval, result.normalized, query.order};
        }
        if (location->has_right()) {
            const auto interval = location->right_interval();
            if (!interval)
                return std::unexpected(supplement(interval.error(), details));
            result.right = IntervalRequest{*interval, result.normalized, query.order};
        }
        if (location->is_knot() && result.left && result.right)
            result.location = SideLocation::C1Join;
        return result;
    }

    static core::Result<Selected> select(const CoefficientBlock& coefficients, const Route& route,
                                         PchipQuery query, const EvaluationPlan& plan,
                                         EvaluationStorage& storage, const std::stop_token& stop,
                                         core::FailureDetails details) noexcept {
        const auto decision = decide_side(route.location, query.side, query.order, details);
        if (!decision)
            return std::unexpected(decision.error());
        std::size_t next = 0;
        const auto evaluate =
            [&](const std::optional<IntervalRequest>& request) -> core::Result<IntervalValues> {
            if (!request || next >= plan.candidates)
                return std::unexpected(failure(core::FailureCode::InvalidInput, details));
            auto& prepared = storage.interval[next++];
            if (!prepared)
                return std::unexpected(failure(core::FailureCode::InvalidInput, details));
            return evaluate_staged_interval(coefficients, *request, plan.interval,
                                            std::move(*prepared), plan.derived, stop, details);
        };
        Selected result;
        result.location.requested = query.side;
        if (const auto* candidate = std::get_if<JetCandidate>(&*decision)) {
            result.location.selected = *candidate;
            const auto& chosen = *candidate == JetCandidate::Left    ? route.left
                                 : *candidate == JetCandidate::Right ? route.right
                                 : route.left                        ? route.left
                                                                     : route.right;
            if (!chosen) {
                if (!route.constant)
                    return std::unexpected(failure(core::FailureCode::InvalidInput, details));
                result.constant = route.constant;
                return result;
            }
            auto values = evaluate(chosen);
            if (!values)
                return std::unexpected(values.error());
            result.location.interval = chosen->interval;
            result.steps = values->used_steps();
            result.values.emplace(std::move(*values));
            return result;
        }
        if (const auto* both = std::get_if<CompareSecondDerivatives>(&*decision)) {
            auto left = evaluate(route.left);
            if (!left)
                return std::unexpected(left.error());
            auto right = evaluate(route.right);
            if (!right)
                return std::unexpected(right.error());
            const auto agreement = agree_side(*both, *left, *right, plan.agreement, stop, details);
            if (!agreement)
                return std::unexpected(agreement.error());
            const auto interval = right->interval();
            if (!interval)
                return std::unexpected(supplement(interval.error(), details));
            result.location.interval = *interval;
            result.location.selected = agreement->selected;
            result.steps = left->used_steps() + right->used_steps() + agreement->used_steps;
            result.values.emplace(std::move(*right));
            return result;
        }
        const auto* held = std::get_if<CompareHeldDerivatives>(&*decision);
        if (held == nullptr)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        auto inside = evaluate(held->inside == JetCandidate::Left ? route.left : route.right);
        if (!inside)
            return std::unexpected(inside.error());
        const auto agreement = agree_side(*held, *inside, plan.agreement, stop, details);
        if (!agreement)
            return std::unexpected(agreement.error());
        const auto interval = inside->interval();
        if (!interval)
            return std::unexpected(supplement(interval.error(), details));
        result.location.interval = *interval;
        result.location.selected = agreement->selected;
        result.steps = inside->used_steps() + agreement->used_steps;
        result.values.emplace(std::move(*inside));
        return result;
    }

    static core::Result<void> append_record(RationalView& record, std::span<std::uint32_t> limbs,
                                            RationalView source, std::size_t& cursor,
                                            Meter& meter) noexcept {
        const auto entered = meter.visit();
        if (!entered)
            return std::unexpected(entered.error());
        const auto copy = [&](std::span<const std::uint32_t> words)
            -> core::Result<std::span<const std::uint32_t>> {
            const auto remaining = limbs.size() - cursor;
            if (words.size() > remaining)
                return std::unexpected(failure(core::FailureCode::ResourceExceeded, meter.details));
            auto destination = limbs.subspan(cursor, words.size());
            for (std::size_t i = 0; i < words.size(); ++i) {
                const auto next = meter.visit();
                if (!next)
                    return std::unexpected(next.error());
                destination[i] = words[i];
            }
            cursor += words.size();
            return destination;
        };
        const auto numerator = copy(source.numerator);
        if (!numerator)
            return std::unexpected(numerator.error());
        const auto denominator = copy(source.denominator);
        if (!denominator)
            return std::unexpected(denominator.error());
        record = {source.sign, *numerator, *denominator};
        return {};
    }

    static core::Result<RationalView> read_selected(const Selected& selected,
                                                    const CoefficientBlock& coefficients,
                                                    std::uint64_t component, JetOrder derivative,
                                                    core::FailureDetails details) noexcept {
        if (selected.values)
            return selected.values->derivative(component, derivative);
        if (!selected.constant)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        return derivative == JetOrder::Value ? coefficients.value(*selected.constant, component)
                                             : core::Result<RationalView>{zero};
    }

    static core::Result<ExactPchipJet>
    evaluate(const PchipCurveReady& curve, const SharedPreparationInput& expected, PchipQuery query,
             CurveWorkspace& workspace, CurveLimits limits, core::ResourcePlan& memory,
             const std::stop_token& stop, core::FailureDetails details) noexcept {
        if (!acquire(workspace))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto plan = plan_evaluation(curve, expected, query, limits, stop, details);
        if (!plan)
            return std::unexpected(plan.error());
        const auto admitted = admit_payload(plan->all, memory, details);
        if (!admitted)
            return std::unexpected(admitted.error());
        auto header = allocate<JetOwner>(memory, 1, details);
        if (!header)
            return std::unexpected(header.error());
        auto records = allocate<RationalView>(memory, plan->records, details);
        if (!records)
            return std::unexpected(records.error());
        auto limbs = allocate<std::uint32_t>(memory, plan->limbs, details);
        if (!limbs)
            return std::unexpected(limbs.error());
        auto storage = provision_evaluation(*plan, memory, stop, details);
        if (!storage)
            return std::unexpected(storage.error());
        auto& result = header->values()[0];
        result.records = std::move(*records);
        result.limbs = std::move(*limbs);
        Meter meter{plan->all.reserved_steps, stop, details};
        std::size_t cursor = 0;
        const auto original = append_record(result.records.values()[0], result.limbs.values(),
                                            query.time, cursor, meter);
        if (!original)
            return std::unexpected(original.error());
        query.time = result.records.values()[0]; // Route only the owned snapshot.
        const auto input = curve.input();
        if (!input)
            return std::unexpected(supplement(input.error(), details));
        const auto& coefficients = curve.owner_->coefficients;
        const auto routed =
            route(coefficients, input->boundary, query, *plan, *storage, stop, details);
        if (!routed)
            return std::unexpected(routed.error());
        auto selected = select(coefficients, *routed, query, *plan, *storage, stop, details);
        if (!selected)
            return std::unexpected(selected.error());
        const auto normalized = append_record(result.records.values()[1], result.limbs.values(),
                                              routed->normalized, cursor, meter);
        if (!normalized)
            return std::unexpected(normalized.error());
        std::size_t record_number = 2;
        for (std::uint64_t c = 0; c < coefficients.components(); ++c) {
            for (std::uint64_t d = 0; d <= static_cast<std::uint64_t>(query.order); ++d) {
                const auto value =
                    read_selected(*selected, coefficients, c, static_cast<JetOrder>(d), details);
                if (!value)
                    return std::unexpected(supplement(value.error(), details));
                const auto copied = append_record(result.records.values()[record_number++],
                                                  result.limbs.values(), *value, cursor, meter);
                if (!copied)
                    return std::unexpected(copied.error());
            }
        }
        const auto completed = meter.visit(wrapper_work);
        if (!completed)
            return std::unexpected(completed.error());
        result.parent = curve.owner_;
        result.stamp = input->stamp;
        result.order = query.order;
        result.components = coefficients.components();
        result.location = selected->location;
        result.work = plan->preflight + meter.used + routed->steps + selected->steps;
        return ExactPchipJet{std::move(*header)};
    }

    static core::Result<CurveAdmission> plan_snapshot(const PchipCurveReady& curve,
                                                      CurveLimits limits,
                                                      const std::stop_token& stop,
                                                      core::FailureDetails details) noexcept {
        if (!curve.active())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        CheckedCount checked;
        const auto records = coefficient_count(curve.owner_->coefficients, checked);
        const CurveAdmission result{checked.mul(records, sizeof(RationalView)), allocator_estimate,
                                    checked.add(records, wrapper_work)};
        if (!checked.ok || result.reserved_steps > limits.arithmetic.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        const auto metadata =
            admit_metadata(limits.metadata, result.metadata_estimate_bytes, details);
        if (!metadata)
            return std::unexpected(metadata.error());
        return result;
    }

    static core::Result<PchipSnapshot> snapshot(const PchipCurveReady& curve,
                                                CurveWorkspace& workspace, CurveLimits limits,
                                                core::ResourcePlan& memory,
                                                const std::stop_token& stop,
                                                core::FailureDetails details) noexcept {
        if (!acquire(workspace))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto plan = plan_snapshot(curve, limits, stop, details);
        if (!plan)
            return std::unexpected(plan.error());
        const auto admitted = admit_payload(*plan, memory, details);
        if (!admitted)
            return std::unexpected(admitted.error());
        const auto records = plan->payload_bytes / sizeof(RationalView);
        auto saved = allocate<RationalView>(memory, records, details);
        if (!saved)
            return std::unexpected(saved.error());
        Meter meter{plan->reserved_steps, stop, details};
        for (std::size_t i = 0; i < saved->size(); ++i) {
            const auto visited = meter.visit();
            if (!visited)
                return std::unexpected(visited.error());
            const auto value = flat_coefficient(curve.owner_->coefficients, i);
            if (!value)
                return std::unexpected(supplement(value.error(), details));
            saved->values()[i] = *value;
        }
        const auto finished = meter.visit(wrapper_work);
        if (!finished)
            return std::unexpected(finished.error());
        return PchipSnapshot{curve.owner_, std::move(*saved)};
    }

    static core::Result<void> compare_inputs(const PreparationRequest& actual,
                                             const PreparationRequest& expected,
                                             Meter& meter) noexcept {
        const auto visited = meter.visit();
        if (!visited)
            return std::unexpected(visited.error());
        constexpr auto mismatch = core::FailureCode::StaleGeneration;
        if (actual.stamp != expected.stamp)
            return std::unexpected(failure(mismatch, meter.details));
        CheckedCount checked;
        const auto value_count =
            checked.mul(actual.samples.times.size(), actual.samples.components);
        if (!checked.ok)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, meter.details));
        if (actual.samples.times.size() < 2 || actual.samples.components == 0 ||
            value_count != actual.samples.values.size() ||
            actual.samples.components != actual.dimensions.size())
            return std::unexpected(failure(core::FailureCode::InvalidInput, meter.details));

        if (actual.samples.components != expected.samples.components ||
            actual.dimensions.size() != expected.dimensions.size() ||
            actual.boundary.policy.index() != expected.boundary.policy.index())
            return std::unexpected(failure(mismatch, meter.details));

        for (std::size_t component = 0; component < expected.dimensions.size(); ++component) {
            const auto next = meter.visit();
            if (!next)
                return std::unexpected(next.error());
            if (actual.dimensions[component].powers != expected.dimensions[component].powers)
                return std::unexpected(failure(mismatch, meter.details));
        }
        using Records = std::span<const RationalView>;
        const auto compare_arrays = [&](Records a, Records b) -> core::Result<void> {
            if (a.size() != b.size())
                return std::unexpected(failure(mismatch, meter.details));
            for (std::size_t i = 0; i < a.size(); ++i) {
                const auto same = equal_record(a[i], b[i], meter, mismatch);
                if (!same)
                    return std::unexpected(same.error());
            }
            return {};
        };
        for (const auto& [record, saved] :
             std::array{std::pair{expected.samples.times, actual.samples.times},
                        std::pair{expected.samples.values, actual.samples.values},
                        std::pair{expected.events, actual.events}}) {
            const auto matched = compare_arrays(record, saved);
            if (!matched)
                return std::unexpected(matched.error());
        }
        const auto begin =
            equal_record(expected.boundary.begin, actual.boundary.begin, meter, mismatch);
        if (!begin)
            return std::unexpected(begin.error());
        const auto end = equal_record(expected.boundary.end, actual.boundary.end, meter, mismatch);
        if (!end)
            return std::unexpected(end.error());
        if (const auto* period = std::get_if<PeriodicBoundary>(&expected.boundary.policy)) {
            const auto* saved = std::get_if<PeriodicBoundary>(&actual.boundary.policy);
            if (saved == nullptr)
                return std::unexpected(failure(mismatch, meter.details));
            return equal_record(period->period, saved->period, meter, mismatch);
        }
        return {};
    }

    static core::Result<RestorePlan> plan_restore(const PchipSnapshotView& saved,
                                                  const SharedPreparationInput& expected,
                                                  CurveLimits limits, const std::stop_token& stop,
                                                  core::FailureDetails details) noexcept {
        const auto current = input_view(expected, details);
        if (!current)
            return std::unexpected(current.error());
        Meter meter{limits.arithmetic.steps, stop, details};
        const auto matched = compare_inputs(saved.input, *current, meter);
        if (!matched)
            return std::unexpected(matched.error());
        auto remaining = limits;
        remaining.arithmetic.steps = meter.remaining();
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        if (remaining.arithmetic.steps == 0)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        if (saved.formula_version != CoefficientBlock::formula_version())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        CheckedCount checked;
        const auto count = checked.mul(
            checked.mul(current->samples.times.size() - 1, current->samples.components), 4);
        if (!checked.ok)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        if (count != saved.coefficients.size())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto prepare = plan_prepare(*current, remaining, stop, details, false);
        if (!prepare)
            return std::unexpected(prepare.error());
        const auto width = prepare->numerics.coefficients().arithmetic().integer_bits();
        Meter saved_scan{remaining.arithmetic.steps - prepare->all.reserved_steps, stop, details};
        std::uint64_t limbs = 0;
        for (const auto coefficient : saved.coefficients) {
            const auto visited = saved_scan.visit();
            if (!visited)
                return std::unexpected(visited.error());
            const auto bound =
                shape(coefficient, std::max(width, limits.arithmetic.input_bits), details);
            if (!bound)
                return std::unexpected(bound.error());
            limbs = checked.add(
                limbs, checked.add(coefficient.numerator.size(), coefficient.denominator.size()));
        }
        const auto comparison = checked.add(count, limbs);
        const auto structural_work = checked.add(meter.used, saved_scan.used);
        RestorePlan result{*prepare, structural_work, comparison, limbs};
        auto& admission = result.prepare.all;
        admission.payload_bytes = checked.add(
            admission.payload_bytes, checked.add(checked.mul(count, sizeof(RationalView)),
                                                 checked.mul(limbs, sizeof(std::uint32_t))));
        admission.metadata_estimate_bytes =
            checked.add(admission.metadata_estimate_bytes, 2 * allocator_estimate);
        admission.reserved_steps = checked.add(
            admission.reserved_steps, checked.add(structural_work, checked.mul(2, comparison)));
        constexpr auto max_index =
            static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
        if (!checked.ok || admission.reserved_steps > limits.arithmetic.steps ||
            count > max_index || limbs > max_index)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        const auto metadata =
            admit_metadata(limits.metadata, admission.metadata_estimate_bytes, details);
        if (!metadata)
            return std::unexpected(metadata.error());
        return result;
    }

    static core::Result<PchipCurveReady>
    restore(const PchipSnapshotView& saved, const SharedPreparationInput& expected,
            CurveWorkspace& workspace, CurveLimits limits, core::ResourcePlan& memory,
            const std::stop_token& stop, core::FailureDetails details) noexcept {
        if (!acquire(workspace))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto plan = plan_restore(saved, expected, limits, stop, details);
        if (!plan)
            return std::unexpected(plan.error());
        const auto admitted = admit_payload(plan->prepare.all, memory, details);
        if (!admitted)
            return std::unexpected(admitted.error());
        auto coefficients = allocate<RationalView>(memory, saved.coefficients.size(), details);
        if (!coefficients)
            return std::unexpected(coefficients.error());
        auto limbs = allocate<std::uint32_t>(memory, plan->limbs, details);
        if (!limbs)
            return std::unexpected(limbs.error());
        auto storage = provision_qualification(plan->prepare, memory, stop, details);
        if (!storage)
            return std::unexpected(storage.error());
        Meter meter{plan->comparison, stop, details};
        std::size_t cursor = 0;
        for (std::size_t i = 0; i < saved.coefficients.size(); ++i) {
            const auto copied = append_record(coefficients->values()[i], limbs->values(),
                                              saved.coefficients[i], cursor, meter);
            if (!copied)
                return std::unexpected(copied.error());
        }
        return qualify(expected, plan->prepare, std::move(*storage), limits, stop, details,
                       plan->preflight + meter.used, coefficients->values(), plan->comparison);
    }
};
} // namespace curve_detail

PchipSnapshot::PchipSnapshot(std::shared_ptr<const curve_detail::CurveOwner> owner,
                             core::PlannedBuffer<RationalView> coefficients) noexcept
    : owner_(std::move(owner)), coefficients_(std::move(coefficients)) {}

core::Result<PchipSnapshotView> PchipSnapshot::view() const& noexcept {
    if (!owner_)
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    const auto input = curve_detail::input_view(owner_->input, {});
    if (!input)
        return std::unexpected(input.error());
    return PchipSnapshotView{*input, CoefficientBlock::formula_version(), coefficients_.values()};
}

bool PchipCurveReady::active() const noexcept {
    return owner_ != nullptr;
}
core::Result<PreparationRequest> PchipCurveReady::input() const& noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return curve_detail::input_view(owner_->input, {});
}
core::Result<SharedPreparationInput> PchipCurveReady::retain_input() const noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return owner_->input;
}
core::Result<std::uint64_t> PchipCurveReady::formula_version() const noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return CoefficientBlock::formula_version();
}
core::Result<RationalView> PchipCurveReady::coefficient(std::uint64_t interval,
                                                        std::uint64_t component,
                                                        std::uint64_t power) const& noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return owner_->coefficients.coefficient(interval, component, power);
}
std::uint64_t PchipCurveReady::used_steps() const noexcept {
    return active() ? owner_->work : 0;
}

ExactPchipJet::ExactPchipJet(core::PlannedBuffer<curve_detail::JetOwner> owner) noexcept
    : owner_(std::move(owner)) {}
ExactPchipJet::ExactPchipJet(ExactPchipJet&&) noexcept = default;
ExactPchipJet& ExactPchipJet::operator=(ExactPchipJet&&) noexcept = default;
ExactPchipJet::~ExactPchipJet() = default;
bool ExactPchipJet::active() const noexcept {
    return !owner_.values().empty();
}
std::uint64_t ExactPchipJet::components() const noexcept {
    return active() ? owner_.values()[0].components : 0;
}
std::uint64_t ExactPchipJet::used_steps() const noexcept {
    return active() ? owner_.values()[0].work : 0;
}
core::Result<PchipInputStamp> ExactPchipJet::stamp() const noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return owner_.values()[0].stamp;
}
core::Result<std::uint64_t> ExactPchipJet::formula_version() const noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return CoefficientBlock::formula_version();
}
core::Result<PchipJetLocation> ExactPchipJet::location() const noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return owner_.values()[0].location;
}
core::Result<JetOrder> ExactPchipJet::order() const noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return owner_.values()[0].order;
}
core::Result<RationalView> ExactPchipJet::original_time() const& noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return owner_.values()[0].records.values()[0];
}
core::Result<RationalView> ExactPchipJet::normalized_time() const& noexcept {
    if (!active())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    return owner_.values()[0].records.values()[1];
}
core::Result<RationalView> ExactPchipJet::derivative(std::uint64_t component,
                                                     JetOrder derivative) const& noexcept {
    if (!active() || component >= components())
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    const auto& owner = owner_.values()[0];
    if (derivative > owner.order)
        return std::unexpected(curve_detail::failure(core::FailureCode::InvalidInput));
    const auto fields = 1 + static_cast<std::uint64_t>(owner.order);
    return owner.records.values()[static_cast<std::size_t>(2 + component * fields +
                                                           static_cast<std::uint64_t>(derivative))];
}

core::Result<CurveAdmission> forecast_pchip_prepare(const PreparationRequest& input,
                                                    CurveLimits limits, const std::stop_token& stop,
                                                    core::FailureDetails details) noexcept {
    const auto plan = curve_detail::Builder::plan_prepare(input, limits, stop, details);
    if (!plan)
        return std::unexpected(plan.error());
    return plan->all;
}
core::Result<PchipCurveReady> prepare_pchip_curve(const PreparationRequest& input,
                                                  CurveWorkspace& workspace, CurveLimits limits,
                                                  core::ResourcePlan& memory,
                                                  const std::stop_token& stop,
                                                  core::FailureDetails details) noexcept {
    return curve_detail::Builder::prepare(input, workspace, limits, memory, stop, details);
}
core::Result<CurveAdmission> forecast_pchip_evaluation(const PchipCurveReady& curve,
                                                       const SharedPreparationInput& expected,
                                                       PchipQuery query, CurveLimits limits,
                                                       const std::stop_token& stop,
                                                       core::FailureDetails details) noexcept {
    const auto plan =
        curve_detail::Builder::plan_evaluation(curve, expected, query, limits, stop, details);
    if (!plan)
        return std::unexpected(plan.error());
    return plan->all;
}
core::Result<ExactPchipJet> evaluate_pchip_curve(const PchipCurveReady& curve,
                                                 const SharedPreparationInput& expected,
                                                 PchipQuery query, CurveWorkspace& workspace,
                                                 CurveLimits limits, core::ResourcePlan& memory,
                                                 const std::stop_token& stop,
                                                 core::FailureDetails details) noexcept {
    return curve_detail::Builder::evaluate(curve, expected, query, workspace, limits, memory, stop,
                                           details);
}
core::Result<CurveAdmission> forecast_pchip_snapshot(const PchipCurveReady& curve,
                                                     CurveLimits limits,
                                                     const std::stop_token& stop,
                                                     core::FailureDetails details) noexcept {
    return curve_detail::Builder::plan_snapshot(curve, limits, stop, details);
}
core::Result<PchipSnapshot> snapshot_pchip_curve(const PchipCurveReady& curve,
                                                 CurveWorkspace& workspace, CurveLimits limits,
                                                 core::ResourcePlan& memory,
                                                 const std::stop_token& stop,
                                                 core::FailureDetails details) noexcept {
    return curve_detail::Builder::snapshot(curve, workspace, limits, memory, stop, details);
}
core::Result<CurveAdmission> forecast_pchip_restore(const PchipSnapshotView& saved,
                                                    const SharedPreparationInput& expected,
                                                    CurveLimits limits, const std::stop_token& stop,
                                                    core::FailureDetails details) noexcept {
    const auto plan = curve_detail::Builder::plan_restore(saved, expected, limits, stop, details);
    if (!plan)
        return std::unexpected(plan.error());
    return plan->prepare.all;
}
core::Result<PchipCurveReady>
restore_pchip_curve(const PchipSnapshotView& saved, const SharedPreparationInput& expected,
                    CurveWorkspace& workspace, CurveLimits limits, core::ResourcePlan& memory,
                    const std::stop_token& stop, core::FailureDetails details) noexcept {
    return curve_detail::Builder::restore(saved, expected, workspace, limits, memory, stop,
                                          details);
}
} // namespace uvlas::aero::detail::pchip
