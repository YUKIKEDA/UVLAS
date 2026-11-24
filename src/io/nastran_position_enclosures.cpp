#include "nastran_position_enclosures.hpp"

#include "../aero/certified_trigonometry_private.hpp"
#include "../structure/pg_reference_arithmetic.hpp"
#include "nastran_isotropic_shell_cards.hpp"
#include "nastran_position_arithmetic.hpp"

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
using nastran_directed::cross;
using nastran_directed::normalized;
using nastran_directed::Position;
using nastran_directed::subtract;
constexpr auto operation = core::Operation::PrepareModel;

[[nodiscard]] std::optional<NastranBulkPosition>
source_position(const std::optional<NastranBulkLexeme>& field) noexcept {
    return field ? std::optional{field->position} : std::nullopt;
}

[[nodiscard]] bool selects(const NastranSupportTerm& term, NastranGridId grid) {
    if (const auto* point = std::get_if<NastranSupportPoint>(&term.target)) {
        return point->grid == grid;
    }
    const auto& range = std::get<NastranSupportRange>(term.target);
    return grid >= range.first && grid <= range.last;
}

class Transform {
public:
    Transform(const NastranCoordinateView& deck, const structure::ExactLimits& exact,
              core::ResourcePlan& plan, m::Work& work, NastranPositionDiagnostic& active) noexcept
        : deck_(deck), exact_(exact), plan_(plan), work_(work), active_(active) {}

    [[nodiscard]] const NastranGridPoint& select(NastranGridId id) {
        active_ = {Reason::Selection, id, {}, {}};
        work_.step();
        const auto source = std::ranges::find_if(deck_.grids, [&](const NastranGridPoint& grid) {
            work_.step();
            return grid.id == id;
        });
        if (source == deck_.grids.end()) {
            partition_fail(Code::InvalidInput);
        }
        return *source;
    }

    [[nodiscard]] NastranGridEnclosure grid(NastranGridId id, const m::Interval& unit) {
        const auto& source = select(id);
        auto position = basic_position(source);
        active_.reason = Reason::Conversion;
        active_.coordinate = source.position_system.id;
        active_.position = source.id_field.position;
        NastranGridEnclosure result;
        result.grid = id;
        result.cp = source.position_system.id;
        result.source = source.id_field.position;
        result.cp_field = source_position(source.position_system.field);
        for (std::size_t i = 0; i < position.size(); ++i) {
            result.tuple_fields[i] = source_position(source.position[i].field);
            position[i] = m::multiply(position[i], unit, work_);
        }
        result.metres = retain_axis(position);
        return result;
    }

    [[nodiscard]] std::array<NastranComponentAxis, 3> component_axes(const NastranGridPoint& point,
                                                                     std::uint8_t components) {
        const auto* definition = coordinate(point.displacement_system);
        const auto cd_kind = kind(definition);
        auto physical = reference_frame(definition);
        if (cd_kind == NastranCoordinateKind::Rectangular) {
            return {retain_axis(physical.axes[0]), retain_axis(physical.axes[1]),
                    retain_axis(physical.axes[2])};
        }

        constexpr unsigned transverse_components =
            (1U << 0U) | (1U << 1U) | (1U << 3U) | (1U << 4U);
        if (cd_kind == NastranCoordinateKind::Cylindrical &&
            (components & transverse_components) == 0) {
            // At a pole the axial components are the same in the cylindrical
            // and defining rectangular frames. Do not infer either transverse axis.
            return {NastranComponentAxis{}, NastranComponentAxis{}, retain_axis(physical.axes[2])};
        }

        Position position;
        if (point.position_system.id == point.displacement_system.id) {
            active_.reason = Reason::Tuple;
            active_.coordinate = point.position_system.id;
            position = cartesian(point.position, cd_kind);
        } else {
            const auto relative = subtract(basic_position(point), physical.origin, work_);
            for (std::size_t i = 0; i < position.size(); ++i) {
                position[i] = dot(relative, physical.axes[i]);
            }
        }

        active_.reason = Reason::DisplacementAxis;
        active_.coordinate = point.displacement_system.id;
        active_.position = source_position(point.displacement_system.field);
        const Position transverse{position[0], position[1], {}};
        const auto radius2 = nastran_directed::norm_squared(transverse, work_);
        if (radius2.upper.mantissa == 0) {
            partition_fail(Code::Unsupported);
        }
        if (radius2.lower.mantissa <= 0) {
            partition_fail(Code::PrecisionUnresolved);
        }
        const auto equatorial = normalized(transverse, work_);
        const Position azimuth{m::subtract({}, equatorial[1], work_), equatorial[0], {}};
        if (cd_kind == NastranCoordinateKind::Cylindrical) {
            return {retain_axis(rotate(physical.axes, equatorial)),
                    retain_axis(rotate(physical.axes, azimuth)), retain_axis(physical.axes[2])};
        }

        const auto radial = normalized(position, work_);
        const auto meridional = cross(azimuth, radial, work_);
        return {retain_axis(rotate(physical.axes, radial)),
                retain_axis(rotate(physical.axes, meridional)),
                retain_axis(rotate(physical.axes, azimuth))};
    }

private:
    using Axes = std::array<Position, 3>;

    struct RectangularFrame {
        Position origin;
        Axes axes;
    };

    [[nodiscard]] Position rotate(const Axes& axes, const Position& position) {
        Position result;
        for (std::size_t i = 0; i < result.size(); ++i) {
            result[i] = m::add(m::add(m::multiply(axes[0][i], position[0], work_),
                                      m::multiply(axes[1][i], position[1], work_), work_),
                               m::multiply(axes[2][i], position[2], work_), work_);
        }
        return result;
    }

    [[nodiscard]] m::Interval dot(const Position& a, const Position& b) {
        return m::add(m::add(m::multiply(a[0], b[0], work_), m::multiply(a[1], b[1], work_), work_),
                      m::multiply(a[2], b[2], work_), work_);
    }

    [[nodiscard]] Position translate(const Position& point, const Position& origin) {
        return {m::add(point[0], origin[0], work_), m::add(point[1], origin[1], work_),
                m::add(point[2], origin[2], work_)};
    }

    [[nodiscard]] NastranComponentAxis retain_axis(const Position& axis) {
        NastranComponentAxis result;
        for (std::size_t i = 0; i < axis.size(); ++i) {
            result[i] = {nastran_directed::retain(axis[i].lower, exact_, plan_, work_),
                         nastran_directed::retain(axis[i].upper, exact_, plan_, work_)};
        }
        return result;
    }

    [[nodiscard]] RectangularFrame read_frame(const NastranCoordinateDefinition& definition,
                                              NastranCoordinateKind parent) {
        work_.step();
        active_.coordinate = definition.id;
        active_.reason = Reason::FrameOrigin;
        auto origin = cartesian(definition.a, parent);

        active_.reason = Reason::FrameAxis;
        const auto point_b = cartesian(definition.b, parent);
        active_.position =
            source_position(definition.b[0].field).value_or(definition.name.position);
        auto z = normalized(subtract(point_b, origin, work_), work_);

        active_.reason = Reason::FramePlane;
        const auto point_c = cartesian(definition.c, parent);
        active_.position =
            source_position(definition.c[0].field).value_or(definition.name.position);
        auto y = normalized(cross(z, subtract(point_c, origin, work_), work_), work_);
        auto x = cross(y, z, work_);
        return {std::move(origin), {std::move(x), std::move(y), std::move(z)}};
    }

    [[nodiscard]] Position basic_position(const NastranGridPoint& source) {
        const auto* frame = coordinate(source.position_system);
        active_.reason = Reason::Tuple;
        active_.coordinate = source.position_system.id;
        auto position = cartesian(source.position, kind(frame));
        while (frame) {
            const auto* parent = coordinate(frame->reference);
            const auto in_parent = read_frame(*frame, kind(parent));
            position = translate(rotate(in_parent.axes, position), in_parent.origin);
            frame = parent;
        }
        return position;
    }

    [[nodiscard]] RectangularFrame reference_frame(const NastranCoordinateDefinition* leaf) {
        const auto one = m::integer(1);
        RectangularFrame result;
        for (std::size_t i = 0; i < result.axes.size(); ++i) {
            result.axes[i][i] = one;
        }
        for (const auto* frame = leaf; frame;) {
            const auto* parent = coordinate(frame->reference);
            const auto in_parent = read_frame(*frame, kind(parent));
            for (auto& direction : result.axes) {
                direction = rotate(in_parent.axes, direction);
            }
            result.origin = translate(rotate(in_parent.axes, result.origin), in_parent.origin);
            frame = parent;
        }
        return result;
    }

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
        const auto value = nastran_directed::fixed_fraction(exact, work_);
        work_.step();
        const auto wrapped = value.numerator % (value.denominator * 360);
        const auto phase = nastran_directed::ratio(wrapped, value.denominator * 180, work_);
        if (!pi_) {
            pi_ = m::pi_constant(work_);
        }
        return m::sine_cosine(phase, m::PhaseUnits::PiMultiple, *pi_, work_);
    }

    [[nodiscard]] Position cartesian(const std::array<NastranCoordinateValue, 3>& raw,
                                     NastranCoordinateKind frame_kind) {
        const auto first = nastran_directed::number(read(raw[0]), work_);
        if (frame_kind == NastranCoordinateKind::Rectangular) {
            return {first, nastran_directed::number(read(raw[1]), work_),
                    nastran_directed::number(read(raw[2]), work_)};
        }
        const auto theta = angle(raw[1]);
        if (frame_kind == NastranCoordinateKind::Cylindrical) {
            return {m::multiply(first, theta.cosine, work_), m::multiply(first, theta.sine, work_),
                    nastran_directed::number(read(raw[2]), work_)};
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
            const auto unit = nastran_directed::number(
                partition_take(materialize_rational(scale, budget.leaf_limits, plan)), work);
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

core::Result<NastranPositionEnclosures>
NastranPositionEnclosures::prepare(const NastranIsotropicShellCards& input, SiLengthUnit length,
                                   std::span<const NastranGridId> grids,
                                   const NastranPositionLimits& limits, core::ResourcePlan& plan,
                                   std::optional<NastranPositionDiagnostic>* diagnostic) noexcept {
    const auto source = input.view();
    return prepare_view(source ? &source->input.input : nullptr, length, grids, limits, plan,
                        diagnostic);
}

core::Result<NastranGridDisplacementAxes>
prepare_nastran_displacement_axes(const NastranIsotropicShellCards& input, NastranGridId id,
                                  std::uint8_t required_axes,
                                  const NastranSpcDirectionLimits& limits, core::ResourcePlan& plan,
                                  std::optional<NastranPositionDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    NastranPositionDiagnostic active;
    auto result = g::partition_entry([&]() -> core::Result<NastranGridDisplacementAxes> {
        try {
            const auto source = input.view();
            if (!source || required_axes == 0 || required_axes > 7) {
                partition_fail(Code::InvalidInput);
            }
            m::Work work{0, limits.maximum_work, limits.stop};
            const structure::PgReferenceLimits admission{limits.exact, limits.maximum_work,
                                                         limits.maximum_work, limits.stop};
            g::ReferenceArithmetic budget{admission, plan};
            budget.admit(18 + 16);
            Transform transform{source->input.input, budget.leaf_limits, plan, work, active};
            auto axes = transform.component_axes(transform.select(id), required_axes);
            work.step();
            return NastranGridDisplacementAxes{std::move(axes), work.used};
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

NastranSpcDirections::NastranSpcDirections(const NastranSpcDirectionOrigin& origin,
                                           std::uint8_t components,
                                           std::array<NastranComponentAxis, 3>&& axes,
                                           std::uint64_t work) noexcept
    : origin_(origin), components_(components), axes_(std::move(axes)), work_(work) {}

NastranSpcDirections::NastranSpcDirections(NastranSpcDirections&& other) noexcept
    : origin_(std::exchange(other.origin_, std::nullopt)),
      components_(other.components_),
      axes_(std::move(other.axes_)),
      work_(other.work_) {}

const NastranSpcDirectionOrigin* NastranSpcDirections::origin() const& noexcept {
    return origin_ ? &*origin_ : nullptr;
}

std::optional<NastranSpcComponentDirection>
NastranSpcDirections::direction(unsigned component) const& noexcept {
    if (!origin_ || component == 0 || component > 6) {
        return std::nullopt;
    }
    const auto zero_based = component - 1;
    if ((components_ & (1U << zero_based)) == 0) {
        return std::nullopt;
    }
    const auto kind =
        component <= 3 ? NastranSpcComponentKind::Translation : NastranSpcComponentKind::Rotation;
    return NastranSpcComponentDirection{kind, axes_[zero_based % axes_.size()]};
}

std::uint64_t NastranSpcDirections::work() const noexcept {
    return origin_ ? work_ : 0;
}

core::Result<NastranSpcDirections> NastranSpcDirections::prepare(
    const NastranCoordinateCards& coordinates, const NastranSupportCards& supports,
    std::size_t term, NastranGridId grid, const NastranSpcDirectionLimits& limits,
    core::ResourcePlan& plan, std::optional<NastranPositionDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    NastranPositionDiagnostic active;
    auto result = g::partition_entry([&]() -> core::Result<NastranSpcDirections> {
        try {
            if (limits.stop.stop_requested()) {
                partition_fail(Code::Cancelled);
            }
            const auto deck = coordinates.view();
            const auto parsed = supports.view();
            if (!deck || !parsed) {
                partition_fail(Code::InvalidInput);
            }
            if (!deck->bulk.same_snapshot(parsed->source.location, parsed->source.digest)) {
                partition_fail(Code::StaleGeneration);
            }
            active = {Reason::SpcTerm, grid, {}, {}};
            if (term >= parsed->terms.size()) {
                partition_fail(Code::InvalidInput);
            }
            const auto& row = parsed->terms[term];
            active.position = row.component_field.position;
            if (!selects(row, grid)) {
                partition_fail(Code::InvalidInput);
            }

            m::Work work{0, limits.maximum_work, limits.stop};
            const structure::PgReferenceLimits admission{limits.exact, limits.maximum_work,
                                                         limits.maximum_work, limits.stop};
            g::ReferenceArithmetic budget{admission, plan};
            // Each axis retains six exact endpoints. The remaining slots
            // cover decimal parsing and the fixed-size interval transform.
            budget.admit(18 + 16);
            Transform transform{*deck, budget.leaf_limits, plan, work, active};
            const auto& point = transform.select(grid);
            auto axes = transform.component_axes(point, row.components);
            work.step();
            return NastranSpcDirections{{parsed->source, row.set, grid,
                                         point.displacement_system.id, point.id_field.position,
                                         row.component_field.position, row.card.position,
                                         source_position(point.displacement_system.field),
                                         source_position(row.motion_field), row.prescribed_motion},
                                        row.components,
                                        std::move(axes),
                                        work.used};
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

} // namespace uvlas::io::detail
