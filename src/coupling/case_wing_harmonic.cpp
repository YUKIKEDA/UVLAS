#include "case_wing_harmonic.hpp"

#include "io/si_rational_projection.hpp"

#include <algorithm>
#include <cmath>
#include <memory_resource>
#include <vector>

namespace uvlas::coupling::detail {
struct CaseWingHarmonicPreparation {
    static CaseHarmonicWake bind(const CaseWingGeometry& owner,
                                 std::array<io::SiEnclosure, 3> velocity) {
        return CaseHarmonicWake{owner, std::move(velocity)};
    }
};

namespace {
namespace a = aero::detail;
namespace p = a::pchip;
namespace pg = structure::detail;
namespace s = io::schema;
namespace si = io::detail;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_mul;
using pg::partition_take;
using enum core::FailureCode;
using Math = pg::ReferenceArithmetic;
using Exact = structure::ExactRational;
using XYZ = std::array<Exact, 3>;
constexpr auto operation = core::Operation::EvaluateAero;

template <class T>
core::PlannedBuffer<T> allocate(std::uint64_t count, Math& m) {
    return partition_take(m.plan.allocate<T>(partition_index(count), operation));
}

template <class T>
class BorrowedInputs {
public:
    BorrowedInputs(std::size_t count, Math& m)
        : storage_(allocate<std::byte>(
              partition_add(partition_mul(count, sizeof(T)), alignof(T) - 1), m)),
          arena_(storage_.values().data(), storage_.size(), std::pmr::null_memory_resource()),
          values_(&arena_) {
        values_.reserve(count);
    }

    void add(const T& value) {
        values_.push_back(value);
    }
    [[nodiscard]] std::span<const T> values() const noexcept {
        return values_;
    }

private:
    core::PlannedBuffer<std::byte> storage_;
    std::pmr::monotonic_buffer_resource arena_;
    std::pmr::vector<T> values_;
};

template <class Record, class Variant>
Record select(const Variant& value) {
    if (const auto* row = std::get_if<Record>(&value))
        return *row;
    partition_fail(Unsupported);
}

template <class Field>
Exact number(const io::SiCase& source, Field field, Math& m) {
    m.charge();
    const auto recipe = partition_take(source.quantity(field)).recipe;
    if (recipe.pi_power)
        partition_fail(InvalidInput);
    return partition_take(io::materialize_rational(recipe, m.leaf_limits, m.plan));
}

template <class Fields>
XYZ vector(const io::SiCase& source, Fields fields, Math& m) {
    if (fields.size() != 3)
        partition_fail(InvalidInput);
    return {number(source, partition_take(fields.at(0)), m),
            number(source, partition_take(fields.at(1)), m),
            number(source, partition_take(fields.at(2)), m)};
}

double unscale(double word, std::int64_t exponent) {
    if (word == 0)
        return word;
    if (exponent < -2200 || exponent > 2200)
        partition_fail(PrecisionUnresolved);
    const auto shift = static_cast<int>(exponent);
    const auto physical = std::scalbn(word, shift);
    if (!std::isfinite(physical) || std::scalbn(physical, -shift) != word)
        partition_fail(PrecisionUnresolved);
    return physical;
}

struct Projection {
    const CaseWingHarmonicProjection& policy;
    std::uint64_t remaining;

    io::SiEnclosure apply(const Exact& value, const CaseFlowComponentBudget& budget, Math& m) {
        const auto encode = [](p::RationalView word) -> si::SiRationalWordsView {
            return {word.sign, word.numerator, word.denominator};
        };
        const auto words =
            partition_take(si::copy_rational_words(value, m.plan, m.limits.cancellation));
        const auto measured = partition_take(si::project_rational_words(
            words.view(), {encode(budget.scale), encode(budget.absolute), encode(budget.relative)},
            policy.representation, {policy.limits.input_integer_bits, remaining},
            m.limits.cancellation, {}));
        remaining -= measured.used_steps;
        return measured.value;
    }

    a::ScalarWords scalar(const Exact& value, const CaseFlowComponentBudget& budget, Math& m) {
        const auto measured = apply(value, budget, m);
        if (const auto* single = std::get_if<io::ScaledBinary64>(&measured))
            return a::SingleWord{unscale(single->center, single->exponent2)};
        const auto pair = std::get<io::ScaledTwoWord>(measured);
        return a::DoubleWord{unscale(pair.high, pair.exponent2), unscale(pair.low, pair.exponent2)};
    }

    NonpenetrationVector vector(const XYZ& value,
                                const std::array<CaseFlowComponentBudget, 3>& budgets, Math& m) {
        return {scalar(value[0], budgets[0], m), scalar(value[1], budgets[1], m),
                scalar(value[2], budgets[2], m)};
    }
};

s::ObservationAero observation(const s::Case& declared, std::uint64_t id, std::string_view kind,
                               Math& m) {
    const auto values = declared.observations;
    for (std::uint64_t i = 0; i < values.size(); ++i) {
        m.visit(1);
        const auto value = partition_take(values.at(partition_index(i))).record();
        const auto* aero = std::get_if<s::ObservationAero>(&value);
        if (aero != nullptr && aero->id.unsigned_value() == id) {
            const auto basis = select<s::BasisWorldMaterial>(aero->basis.record());
            if (basis.kind.text() != "world" || aero->quantity.text() != kind)
                partition_fail(InvalidInput);
            return *aero;
        }
    }
    partition_fail(InvalidInput);
}

bool contains(const s::ObservationAero& observation, std::uint64_t surface, Math& m) {
    for (std::uint64_t i = 0; i < observation.surfaces.size(); ++i) {
        m.visit(1);
        if (partition_take(observation.surfaces.at(partition_index(i))).unsigned_value() == surface)
            return true;
    }
    return false;
}

s::ObservationAero group(const io::SiCase& source, const s::Case& declared,
                         CaseWingHarmonicGroup ids, const XYZ& origin, Math& m) {
    const auto force = observation(declared, ids.force_observation, "force", m);
    const auto moment = observation(declared, ids.moment_observation, "moment", m);
    if (force.surfaces.size() != moment.surfaces.size() ||
        vector(source, force.origin, m) != origin || vector(source, moment.origin, m) != origin)
        partition_fail(InvalidInput);
    for (std::uint64_t i = 0; i < force.surfaces.size(); ++i) {
        m.visit(1);
        const auto surface = partition_take(force.surfaces.at(partition_index(i))).unsigned_value();
        if (!contains(moment, surface, m))
            partition_fail(InvalidInput);
    }
    return force;
}

CaseWingHarmonicResult
run(const io::SiCase& source, const CaseWingGeometry& wing, const a::GeometryContext& context,
    const CaseWingHarmonicRequest& request, const CaseWingHarmonicProjection& projection,
    const structure::PgReferenceLimits& exact, const DiscreteNonpenetrationLimits& limits,
    const DiscreteCirculationRequest& solving, core::ResourcePlan& memory) {
    Math m{exact, memory};
    m.cancel();
    partition_take(wing.validate_for(source, context, {}));
    const auto wake = wing.straight_wake();
    if (!wake || wing.wake_bands().size() != wing.trailing_bands().size() ||
        wing.wake_bands().empty())
        partition_fail(StaleGeneration);
    m.admit();
    const auto model = partition_take(partition_take(source.source_case()).get().root());
    const auto flow = select<s::FlowKnown>(model.flow.record());
    const auto& velocity = wake->velocity_metres_per_second;
    const auto origin = vector(source, flow.reference.record().moment_origin, m);
    const std::array groups{group(source, model, request.observations[0], origin, m),
                            group(source, model, request.observations[1], origin, m)};
    const auto reference_surfaces = model.model.record().aero_surfaces;
    for (std::uint64_t i = 0; i < reference_surfaces.size(); ++i) {
        m.visit(1);
        const auto surface = partition_take(reference_surfaces.at(partition_index(i))).record();
        if (!std::holds_alternative<s::AeroAttachmentFixed>(surface.attachment.record()))
            partition_fail(Unsupported);
    }
    BorrowedInputs<DiscreteCollocation> rows{wing.cells().size(), m};
    auto sources = allocate<HarmonicLoadSource>(wing.cells().size(), m);
    for (std::size_t i = 0; i < wing.cells().size(); ++i) {
        m.visit(1);
        const auto& cell = wing.cells()[i];
        const bool first = contains(groups[0], cell.surface, m);
        const bool second = contains(groups[1], cell.surface, m);
        if (first == second)
            partition_fail(InvalidInput);
        sources.values()[i] = {cell.face, cell.surface, first ? 0U : 1U};
        rows.add(partition_take(wing.collocation(i)));
    }
    BorrowedInputs<DiscreteWakePanel> strips{wing.wake_bands().size(), m};
    for (std::size_t i = 0; i < wing.wake_bands().size(); ++i) {
        m.visit(1);
        const auto& band = wing.wake_bands()[i];
        const auto& transfer = wing.trailing_bands()[i].transfer;
        if (band.trailing_edge != transfer.edge || band.bound_face != transfer.bound_face)
            partition_fail(StaleGeneration);
        strips.add({band.facet, AttachedWakeFront{}, SteadyWakeClosure{}});
    }
    const auto& contact = partition_take(wing.global()).get();
    const auto model_generation = context.model.value();
    const DiscreteNonpenetrationInput input{
        contact,
        rows.values(),
        strips.values(),
        {context.diagnostics.parent.value_or(core::AcceptedStateId{0}), {}},
        {model_generation, context.current.value()},
        {p::PchipIssuerId{model_generation}, p::PchipInputGeneration{context.current.value()}},
        p::RequestedSide::Automatic};
    Projection put{projection, projection.limits.steps};
    const StraightWakeHarmonicInput phase{
        CaseWingHarmonicPreparation::bind(wing,
                                          {put.apply(velocity[0], projection.velocity[0], m),
                                           put.apply(velocity[1], projection.velocity[1], m),
                                           put.apply(velocity[2], projection.velocity[2], m)}),
        request.age_panel_starts_seconds, request.omega_radians_per_second,
        request.maximum_age_samples};
    const HarmonicLoadInput load{put.scalar(number(source, flow.density, m), projection.density, m),
                                 put.vector(origin, projection.origin, m)};
    auto solved = partition_take(
        solve_discrete_harmonic_batch(input, phase, {request.drives, load, sources.values()},
                                      limits, solving, memory, exact.cancellation));
    return {partition_take(wing.source()), request.observations, std::move(solved)};
}
} // namespace

core::Result<CaseWingHarmonicResult> solve_case_wing_harmonic(
    const io::SiCase& source, const CaseWingGeometry& wing, const a::GeometryContext& context,
    const CaseWingHarmonicRequest& request, const CaseWingHarmonicProjection& projection,
    const structure::PgReferenceLimits& exact, const DiscreteNonpenetrationLimits& limits,
    const DiscreteCirculationRequest& solving, core::ResourcePlan& memory) noexcept {
    return core::guard_entry(operation, PrecisionUnresolved, context.diagnostics,
                             [&]() -> core::Result<CaseWingHarmonicResult> {
                                 try {
                                     return run(source, wing, context, request, projection, exact,
                                                limits, solving, memory);
                                 } catch (const pg::PartitionAbort& error) {
                                     return std::unexpected(core::Failure{
                                         error.failure.code, operation, error.failure.details});
                                 }
                             });
}
} // namespace uvlas::coupling::detail
