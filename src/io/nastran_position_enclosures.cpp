#include "nastran_position_enclosures.hpp"

#include "../aero/certified_trigonometry_private.hpp"
#include "../structure/pg_rank_storage.hpp"
#include "../structure/pg_reference_arithmetic.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

namespace m = aero::detail::p1;
namespace g = structure::detail;
using g::partition_fail;
using g::partition_index;
using g::partition_take;
using Rational = structure::ExactRational;
using Code = core::FailureCode;
using Reason = NastranPositionFailure;
using Position = std::array<m::Interval, 3>;
using FixedInteger = aero::detail::ExactInteger;
constexpr auto operation = core::Operation::PrepareModel;

[[nodiscard]] g::RankFraction<FixedInteger> fixed_fraction(const Rational& source, m::Work& work) {
    work.step();
    constexpr auto width = static_cast<std::uint64_t>(std::numeric_limits<g::RankWide>::digits);
    if (std::max(source.numerator_bits(), source.denominator_bits()) > width) {
        partition_fail(Code::ResourceExceeded);
    }
    return g::rank_fraction<FixedInteger>(source);
}

[[nodiscard]] m::Interval ratio(const FixedInteger& numerator, const FixedInteger& denominator,
                                m::Work& work) {
    return m::divide(m::exact(numerator, 0), m::exact(denominator, 0), work);
}

[[nodiscard]] Rational retain(const m::Bound& input, const structure::ExactLimits& exact,
                              core::ResourcePlan& plan, m::Work& work) {
    work.step();
    if (input.mantissa == 0) {
        return {};
    }
    constexpr auto width = std::numeric_limits<g::RankWide>::digits;
    if (input.power < -width + 1 ||
        static_cast<std::int64_t>(m::bits(input.mantissa)) + input.power > width) {
        partition_fail(Code::ResourceExceeded);
    }
    const g::RankWide units{input.mantissa};
    if (input.power < 0) {
        return g::rank_materialize_dyadic(units, static_cast<std::uint32_t>(-input.power), exact,
                                          plan);
    }
    auto shifted = g::rank_abs(units) << static_cast<unsigned>(input.power);
    if (units < 0) {
        shifted = -shifted;
    }
    return g::rank_exact(g::RankFraction<g::RankWide>{std::move(shifted), g::RankWide{1}}, exact,
                         plan);
}

[[nodiscard]] Position subtract(const Position& first, const Position& second, m::Work& work) {
    return {m::subtract(first[0], second[0], work), m::subtract(first[1], second[1], work),
            m::subtract(first[2], second[2], work)};
}

[[nodiscard]] Position cross(const Position& first, const Position& second, m::Work& work) {
    Position output;
    for (std::size_t i = 0; i < output.size(); ++i) {
        const auto j = (i + 1) % 3, k = (i + 2) % 3;
        output[i] = m::subtract(m::multiply(first[j], second[k], work),
                                m::multiply(first[k], second[j], work), work);
    }
    return output;
}

[[nodiscard]] Position normalized(const Position& direction, m::Work& work) {
    m::Interval squared;
    for (const auto& value : direction) {
        squared = m::add(squared, m::multiply(value, value, work), work);
    }
    if (squared.upper.mantissa == 0) {
        partition_fail(Code::GeometryInvalid);
    }
    if (squared.lower.mantissa <= 0) {
        partition_fail(Code::PrecisionUnresolved);
    }
    const auto length = m::square_root(squared, work);
    return {m::divide(direction[0], length, work), m::divide(direction[1], length, work),
            m::divide(direction[2], length, work)};
}

[[nodiscard]] std::optional<NastranBulkPosition>
source_position(const std::optional<NastranBulkLexeme>& field) noexcept {
    return field ? std::optional{field->position} : std::nullopt;
}

class Transform {
public:
    Transform(const NastranCoordinateView& deck, const structure::ExactLimits& exact,
              core::ResourcePlan& plan, m::Work& work, NastranPositionDiagnostic& active) noexcept
        : deck_(deck), exact_(exact), plan_(plan), work_(work), active_(active) {}

    [[nodiscard]] NastranGridEnclosure grid(NastranGridId id, const m::Interval& unit) {
        active_ = {Reason::Selection, id, {}, {}};
        work_.step();
        const auto source = std::ranges::find_if(deck_.grids, [&](const NastranGridPoint& grid) {
            work_.step();
            return grid.id == id;
        });
        if (source == deck_.grids.end()) {
            partition_fail(Code::InvalidInput);
        }

        const auto* frame = coordinate(source->position_system);
        active_.reason = Reason::Tuple;
        active_.coordinate = source->position_system.id;
        auto position = cartesian(source->position, kind(frame));
        while (frame) {
            work_.step();
            const auto* parent = coordinate(frame->reference);
            active_.coordinate = frame->id;
            active_.reason = Reason::FrameOrigin;
            const auto origin = cartesian(frame->a, kind(parent));

            active_.reason = Reason::FrameAxis;
            const auto point_b = cartesian(frame->b, kind(parent));
            active_.position = source_position(frame->b[0].field).value_or(frame->name.position);
            const auto z = normalized(subtract(point_b, origin, work_), work_);

            active_.reason = Reason::FramePlane;
            const auto point_c = cartesian(frame->c, kind(parent));
            active_.position = source_position(frame->c[0].field).value_or(frame->name.position);
            const auto y = normalized(cross(z, subtract(point_c, origin, work_), work_), work_);
            const auto x = cross(y, z, work_);
            Position parent_coordinates;
            for (std::size_t i = 0; i < position.size(); ++i) {
                parent_coordinates[i] =
                    m::add(origin[i],
                           m::add(m::add(m::multiply(x[i], position[0], work_),
                                         m::multiply(y[i], position[1], work_), work_),
                                  m::multiply(z[i], position[2], work_), work_),
                           work_);
            }
            position = std::move(parent_coordinates);
            frame = parent;
        }

        active_.reason = Reason::Conversion;
        active_.coordinate = source->position_system.id;
        active_.position = source->id_field.position;
        NastranGridEnclosure result;
        result.grid = id;
        result.cp = source->position_system.id;
        result.source = source->id_field.position;
        result.cp_field = source_position(source->position_system.field);
        for (std::size_t i = 0; i < position.size(); ++i) {
            result.tuple_fields[i] = source_position(source->position[i].field);
            const auto metres = m::multiply(position[i], unit, work_);
            result.metres[i] = {retain(metres.lower, exact_, plan_, work_),
                                retain(metres.upper, exact_, plan_, work_)};
        }
        return result;
    }

    [[nodiscard]] m::Interval number(const Rational& value) const {
        const auto value_fraction = fixed_fraction(value, work_);
        return ratio(value_fraction.numerator, value_fraction.denominator, work_);
    }

private:
    [[nodiscard]] const NastranCoordinateDefinition*
    coordinate(const NastranCoordinateReference& reference) {
        if (reference.id.value() == 0) {
            return nullptr;
        }
        active_.reason = Reason::CoordinateDefinition;
        active_.coordinate = reference.id;
        active_.position = source_position(reference.field);
        const auto found =
            std::ranges::find_if(deck_.coordinates, [&](const NastranCoordinateDefinition& row) {
                work_.step();
                return row.id == reference.id;
            });
        if (found == deck_.coordinates.end()) {
            partition_fail(Code::InvalidInput);
        }
        return &*found;
    }

    [[nodiscard]] static NastranCoordinateKind
    kind(const NastranCoordinateDefinition* definition) noexcept {
        return definition ? definition->kind : NastranCoordinateKind::Rectangular;
    }

    [[nodiscard]] Rational read(const NastranCoordinateValue& field) {
        active_.position = source_position(field.field);
        work_.step();
        return partition_take(Rational::from_decimal(field.decimal(), exact_, plan_));
    }

    [[nodiscard]] m::SinCos angle(const NastranCoordinateValue& field) {
        const auto exact = read(field);
        const auto value = fixed_fraction(exact, work_);
        work_.step();
        const auto wrapped = value.numerator % (value.denominator * 360);
        const auto phase = ratio(wrapped, value.denominator * 180, work_);
        if (!pi_) {
            pi_ = m::pi_constant(work_);
        }
        return m::sine_cosine(phase, m::PhaseUnits::PiMultiple, *pi_, work_);
    }

    [[nodiscard]] Position cartesian(const std::array<NastranCoordinateValue, 3>& raw,
                                     NastranCoordinateKind frame_kind) {
        const auto first = number(read(raw[0]));
        if (frame_kind == NastranCoordinateKind::Rectangular) {
            return {first, number(read(raw[1])), number(read(raw[2]))};
        }
        const auto theta = angle(raw[1]);
        if (frame_kind == NastranCoordinateKind::Cylindrical) {
            return {m::multiply(first, theta.cosine, work_), m::multiply(first, theta.sine, work_),
                    number(read(raw[2]))};
        }
        const auto phi = angle(raw[2]);
        const auto radial_xy = m::multiply(first, theta.sine, work_);
        return {m::multiply(radial_xy, phi.cosine, work_), m::multiply(radial_xy, phi.sine, work_),
                m::multiply(first, theta.cosine, work_)};
    }

    const NastranCoordinateView& deck_;
    const structure::ExactLimits& exact_;
    core::ResourcePlan& plan_;
    m::Work& work_;
    NastranPositionDiagnostic& active_;
    std::optional<m::Interval> pi_;
};

} // namespace

NastranPositionEnclosures::NastranPositionEnclosures(
    RootSnapshotLocation source, SourceDigest digest, SiLengthUnit length,
    core::PlannedBuffer<NastranGridEnclosure>&& grids, std::uint64_t work) noexcept
    : source_(source),
      digest_(digest),
      length_unit_(length),
      grids_(std::move(grids)),
      work_(work) {}

NastranPositionEnclosures::NastranPositionEnclosures(NastranPositionEnclosures&& other) noexcept
    : source_(std::exchange(other.source_, std::nullopt)),
      digest_(other.digest_),
      length_unit_(other.length_unit_),
      grids_(std::move(other.grids_)),
      work_(other.work_) {}

std::optional<NastranPositionView> NastranPositionEnclosures::view() const& noexcept {
    if (!source_) {
        return std::nullopt;
    }
    return NastranPositionView{*source_, digest_, length_unit_, grids_.values(), work_};
}

core::Result<NastranPositionEnclosures> NastranPositionEnclosures::prepare_view(
    const NastranCoordinateView* input, SiLengthUnit length,
    std::span<const NastranGridId> requested, const NastranPositionLimits& limits,
    core::ResourcePlan& plan, std::optional<NastranPositionDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    NastranPositionDiagnostic active;
    auto result = g::partition_entry([&]() -> core::Result<NastranPositionEnclosures> {
        try {
            if (limits.stop.stop_requested()) {
                partition_fail(Code::Cancelled);
            }
            if (!input) {
                partition_fail(Code::InvalidInput);
            }
            if (requested.size() > limits.maximum_positions) {
                partition_fail(Code::ResourceExceeded);
            }
            if (requested.empty()) {
                return NastranPositionEnclosures{
                    input->bulk.location, input->bulk.digest, length, {}, 0};
            }
            if (g::partition_mul(requested.size(), sizeof(NastranGridEnclosure)) >
                plan.remaining_bytes()) {
                partition_fail(Code::ResourceExceeded);
            }
            m::Work work{0, limits.maximum_work, limits.stop};
            const structure::PgReferenceLimits admission{limits.exact, limits.maximum_work,
                                                         limits.maximum_work, limits.stop};
            g::ReferenceArithmetic budget{admission, plan};
            budget.admit(g::partition_add(16, g::partition_mul(6, requested.size())));
            Transform transform{*input, budget.leaf_limits, plan, work, active};
            const SiUnits units{length, SiMassUnit::Kilogram, SiTimeUnit::Second,
                                SiAngleUnit::Radian};
            const SiRecipeView scale{SiDecimalView{1, "1", 0}, {}, units, {{1, 0, 0, 0}}, 0};
            active.reason = Reason::Conversion;
            work.step();
            const auto unit = transform.number(
                partition_take(materialize_rational(scale, budget.leaf_limits, plan)));
            auto rows = partition_take(
                plan.allocate<NastranGridEnclosure>(partition_index(requested.size()), operation));
            for (std::size_t i = 0; i < requested.size(); ++i) {
                rows.values()[i] = transform.grid(requested[i], unit);
            }
            if (limits.stop.stop_requested()) {
                partition_fail(Code::Cancelled);
            }
            return NastranPositionEnclosures{input->bulk.location, input->bulk.digest, length,
                                             std::move(rows), work.used};
        } catch (const m::ArithmeticFailure& failure) {
            partition_fail(failure.code);
        } catch (const std::overflow_error&) {
            partition_fail(Code::PrecisionUnresolved);
        } catch (const std::range_error&) {
            partition_fail(Code::PrecisionUnresolved);
        }
    });
    if (!result && diagnostic && result.error().code != Code::Cancelled) {
        *diagnostic = active;
    }
    return result;
}

core::Result<NastranPositionEnclosures>
NastranPositionEnclosures::prepare(const NastranCoordinateCards& input, SiLengthUnit length,
                                   std::span<const NastranGridId> grids,
                                   const NastranPositionLimits& limits, core::ResourcePlan& plan,
                                   std::optional<NastranPositionDiagnostic>* diagnostic) noexcept {
    const auto source = input.view();
    return prepare_view(source ? &*source : nullptr, length, grids, limits, plan, diagnostic);
}

core::Result<NastranPositionEnclosures>
NastranPositionEnclosures::prepare(const NastranShellCards& input, SiLengthUnit length,
                                   std::span<const NastranGridId> grids,
                                   const NastranPositionLimits& limits, core::ResourcePlan& plan,
                                   std::optional<NastranPositionDiagnostic>* diagnostic) noexcept {
    const auto source = input.view();
    return prepare_view(source ? &source->input : nullptr, length, grids, limits, plan, diagnostic);
}

} // namespace uvlas::io::detail
