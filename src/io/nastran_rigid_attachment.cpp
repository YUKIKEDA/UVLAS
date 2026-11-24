#include "nastran_rigid_attachment.hpp"

#include "../structure/pg_reference_polynomials.hpp"
#include "nastran_position_arithmetic.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace uvlas::io::detail {
namespace {

namespace g = structure::detail;
namespace d = nastran_directed;
namespace m = d::m;
using g::partition_fail;
using g::partition_take;
using g::ReferenceArithmetic;
using Rational = structure::ExactRational;
using Code = core::FailureCode;
using Reason = NastranAttachmentFailure;
constexpr auto operation = core::Operation::PrepareModel;
constexpr unsigned translation_mask = 0b000111U;

struct Calculation {
    ReferenceArithmetic& math;
    m::Work& work;
    NastranAttachmentDiagnostic& active;

    [[nodiscard]] NastranCartesianInterval retain(const m::Interval& value) const {
        return {d::retain(value.lower, math.leaf_limits, math.plan, work),
                d::retain(value.upper, math.leaf_limits, math.plan, work)};
    }

    [[nodiscard]] NastranRigidIntervalVector retain(const d::Position& row) const {
        return {retain(row[0]), retain(row[1]), retain(row[2])};
    }

    void grid_site(Reason reason, const NastranRigidGridReference& point) const {
        active.reason = reason;
        active.grid = point.id;
        active.position = point.field.position;
        active.coordinate.reset();
    }

    [[noreturn]] void forward(const core::Failure& failure,
                              const std::optional<NastranPositionDiagnostic>& detail) const {
        if (detail) {
            if (detail->position) {
                active.position = detail->position;
            }
            active.coordinate = detail->coordinate;
        }
        throw g::PartitionAbort{failure};
    }
};

[[nodiscard]] d::Position add(const d::Position& a, const d::Position& b, m::Work& work) {
    return {m::add(a[0], b[0], work), m::add(a[1], b[1], work), m::add(a[2], b[2], work)};
}

[[nodiscard]] d::Position scaled(const d::Position& a, const m::Interval& x, m::Work& work) {
    return {m::multiply(a[0], x, work), m::multiply(a[1], x, work), m::multiply(a[2], x, work)};
}

[[nodiscard]] m::Interval dot(const d::Position& a, const d::Position& b, m::Work& work) {
    return m::add(m::add(m::multiply(a[0], b[0], work), m::multiply(a[1], b[1], work), work),
                  m::multiply(a[2], b[2], work), work);
}

[[nodiscard]] d::Position read(const NastranRigidExactVector& source, m::Work& work) {
    return {d::number(source[0], work), d::number(source[1], work), d::number(source[2], work)};
}

[[nodiscard]] d::Position read(const NastranRigidIntervalVector& source, m::Work& work) {
    d::Position row;
    for (std::size_t i = 0; i < row.size(); ++i) {
        row[i] = {d::number(source[i].lower, work).lower, d::number(source[i].upper, work).upper};
    }
    return row;
}

[[nodiscard]] d::Position grid_position(const NastranShellCards& shells,
                                        const NastranRigidGridReference& point, SiLengthUnit length,
                                        const Calculation& c) {
    c.grid_site(Reason::ReferencePosition, point);
    const std::array selection{point.id};
    std::optional<NastranPositionDiagnostic> problem;
    const auto owner = NastranPositionEnclosures::prepare(
        shells, length, selection, {c.math.leaf_limits, 1, c.work.limit - c.work.used, c.work.stop},
        c.math.plan, &problem);
    if (!owner) {
        c.forward(owner.error(), problem);
    }
    const auto view = owner->view();
    if (!view || view->grids.size() != 1) {
        partition_fail(Code::InvalidInput);
    }
    c.work.used += view->work;
    return read(view->grids.front().metres, c.work);
}

[[nodiscard]] NastranRigidPointReference dependent_reference(const NastranShellCards& shells,
                                                             const NastranRigidGridReference& point,
                                                             const d::Position& master,
                                                             SiLengthUnit length, std::uint8_t axes,
                                                             const Calculation& c) {
    const auto position = grid_position(shells, point, length, c);
    const auto relative = d::subtract(position, master, c.work);
    c.grid_site(Reason::DisplacementAxes, point);
    std::optional<NastranPositionDiagnostic> problem;
    auto directions = prepare_nastran_displacement_axes(
        shells, point.id, axes, {c.math.leaf_limits, c.work.limit - c.work.used, c.work.stop},
        c.math.plan, &problem);
    if (!directions) {
        c.forward(directions.error(), problem);
    }
    c.work.used += directions->work;
    NastranRigidPointReference row{point, c.retain(position), c.retain(relative), {}};
    for (std::size_t i = 0; i < row.component_axes.size(); ++i) {
        if ((axes & (1U << i)) != 0) {
            row.component_axes[i] = std::move(directions->directions[i]);
        }
    }
    return row;
}

[[nodiscard]] std::span<const NastranRigidGridReference>
selected_dependents(const NastranRigidLinkView& owner, const NastranRigidLink& link) {
    if (link.dependent_count == 0 || link.first_dependent > owner.dependents.size() ||
        link.dependent_count > owner.dependents.size() - link.first_dependent) {
        partition_fail(Code::InvalidInput);
    }
    return owner.dependents.subspan(link.first_dependent, link.dependent_count);
}

[[nodiscard]] bool selected(std::uint8_t mask, std::size_t axis, unsigned offset = 0) noexcept {
    return (mask & (1U << (axis + offset))) != 0;
}

void verify_loads(const NastranRigidReferenceView& owner,
                  std::span<const NastranRigidPointLoad> loads, const Calculation& c) {
    if (loads.size() != owner.dependents.size()) {
        c.active.reason = Reason::LoadSelection;
        c.active.position = owner.card.components_field.position;
        partition_fail(Code::InvalidInput);
    }
    for (std::size_t j = 0; j < loads.size(); ++j) {
        c.math.visit(1);
        const auto& point = owner.dependents[j].source;
        const auto& entry = loads[j];
        if (entry.point != point.id) {
            c.grid_site(Reason::LoadSelection, point);
            partition_fail(Code::InvalidInput);
        }
        for (std::size_t i = 0; i < entry.force_n.size(); ++i) {
            if ((!selected(owner.card.dependent_components, i) && entry.force_n[i].sign() != 0) ||
                (!selected(owner.card.dependent_components, i, 3) &&
                 entry.couple_nm[i].sign() != 0)) {
                c.grid_site(Reason::LoadSelection, point);
                c.active.position = owner.card.components_field.position;
                partition_fail(Code::InvalidInput);
            }
        }
    }
}

struct WorldMotion {
    d::Position origin, velocity, spin;
    std::array<d::Position, 3> rotation_rows;
};

[[nodiscard]] WorldMotion import_motion(const NastranRigidMotion& input, const Calculation& c) {
    g::validate_frame(c.math, {input.rotation, input.origin_metres, c.math.constant("1")});
    WorldMotion motion{read(input.origin_metres, c.work),
                       read(input.world_velocity_metres_per_second, c.work),
                       read(input.world_spin_rad_per_second, c.work),
                       {}};
    for (std::size_t i = 0; i < motion.rotation_rows.size(); ++i) {
        for (std::size_t j = 0; j < motion.rotation_rows[i].size(); ++j) {
            motion.rotation_rows[i][j] = d::number(input.rotation[i * 3 + j], c.work);
        }
    }
    return motion;
}

[[nodiscard]] d::Position rotate(const WorldMotion& motion, const d::Position& vector,
                                 m::Work& work) {
    return {dot(motion.rotation_rows[0], vector, work), dot(motion.rotation_rows[1], vector, work),
            dot(motion.rotation_rows[2], vector, work)};
}

struct Wrench {
    d::Position force, moment;
};

[[nodiscard]] NastranRigidProjectedPoint projection(const NastranRigidPointReference& reference,
                                                    const NastranRigidPointLoad& applied,
                                                    std::uint8_t mask, const WorldMotion& motion,
                                                    Wrench& total, const Calculation& c) {
    const auto source_offset = read(reference.offset_from_master_metres, c.work);
    const auto lever = rotate(motion, source_offset, c.work);

    const auto rigid_target = add(motion.origin, lever, c.work);
    const auto velocity = add(motion.velocity, d::cross(motion.spin, lever, c.work), c.work);
    NastranRigidProjectedPoint row{reference.source.id, c.retain(rigid_target), {}, {}, {}, {}};
    if ((mask & (translation_mask << 3U)) == (translation_mask << 3U)) {
        row.rigid_target_orientation.emplace();
    }
    Wrench local;
    for (std::size_t i = 0; i < reference.component_axes.size(); ++i) {
        if (const auto& stored = reference.component_axes[i]) {
            const auto axis = read(*stored, c.work);
            if (selected(mask, i)) {
                row.translation_target[i] = c.retain(dot(axis, rigid_target, c.work));
                row.translational_velocity[i] = c.retain(dot(axis, velocity, c.work));
                const auto effort = d::number(applied.force_n[i], c.work);
                local.force = add(local.force, scaled(axis, effort, c.work), c.work);
            }
            if (selected(mask, i, 3)) {
                row.spin[i] = c.retain(dot(axis, motion.spin, c.work));
                if (auto& orientation = row.rigid_target_orientation) {
                    (*orientation)[i] = c.retain(rotate(motion, axis, c.work));
                }
                const auto effort = d::number(applied.couple_nm[i], c.work);
                local.moment = add(local.moment, scaled(axis, effort, c.work), c.work);
            }
        }
    }
    total.force = add(total.force, local.force, c.work);
    total.moment =
        add(total.moment, add(d::cross(lever, local.force, c.work), local.moment, c.work), c.work);
    return row;
}

} // namespace

NastranRigidAttachment::NastranRigidAttachment(
    NastranRigidLinkSource source, const NastranRigidLink& card,
    NastranRigidIntervalVector&& master,
    core::PlannedBuffer<NastranRigidPointReference>&& dependents,
    structure::PgReferenceCallInfo work, std::uint64_t directed) noexcept
    : source_(source),
      card_(card),
      master_(std::move(master)),
      dependents_(std::move(dependents)),
      exact_work_(work),
      directed_work_(directed) {}

std::optional<NastranRigidReferenceView> NastranRigidAttachment::view() const& noexcept {
    if (dependents_.size() == 0) {
        return std::nullopt;
    }
    return NastranRigidReferenceView{source_,     card_,         master_, dependents_.values(),
                                     exact_work_, directed_work_};
}

core::Result<NastranRigidAttachment> NastranRigidAttachment::prepare(
    const NastranShellCards& shells, const NastranRigidLinkCards& links, NastranElementId id,
    SiLengthUnit length, const NastranRigidAttachmentLimits& limits, core::ResourcePlan& plan,
    std::optional<NastranAttachmentDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    NastranAttachmentDiagnostic active;
    auto result = g::partition_entry([&]() -> core::Result<NastranRigidAttachment> {
        try {
            ReferenceArithmetic math{limits.arithmetic, plan};
            m::Work work{0, limits.maximum_directed_work, limits.arithmetic.cancellation};
            const Calculation c{math, work, active};
            math.cancel();
            const auto model = shells.view();
            const auto deck = links.view();
            if (!model || !deck) {
                partition_fail(Code::InvalidInput);
            }
            if (!model->input.bulk.same_snapshot(deck->source.location, deck->source.digest)) {
                partition_fail(Code::StaleGeneration);
            }
            active = {Reason::Selection, id, {}, {}, {}};
            const auto found = std::ranges::find_if(deck->links, [&](const auto& row) {
                math.visit(1);
                return row.id == id;
            });
            if (found == deck->links.end()) {
                partition_fail(Code::InvalidInput);
            }
            const auto& card = *found;
            const auto dependents = selected_dependents(*deck, card);
            if (dependents.size() > limits.maximum_dependents) {
                partition_fail(Code::ResourceExceeded);
            }
            math.admit(g::partition_add(64, g::partition_mul(30, dependents.size())));
            if (const auto& alpha = card.alpha) {
                active.reason = Reason::Thermal;
                active.position = alpha->field.position;
                math.charge();
                const auto exact = partition_take(
                    Rational::from_decimal(alpha->value.decimal(), math.leaf_limits, plan));
                if (exact.sign() != 0) {
                    partition_fail(Code::Unsupported);
                }
            }

            const auto axes = static_cast<std::uint8_t>(
                (card.dependent_components | (card.dependent_components >> 3U)) & translation_mask);
            const auto master = grid_position(shells, card.independent, length, c);
            auto points = partition_take(plan.allocate<NastranRigidPointReference>(
                g::partition_index(dependents.size()), operation));
            for (std::size_t i = 0; i < dependents.size(); ++i) {
                math.visit(1);
                points.values()[i] =
                    dependent_reference(shells, dependents[i], master, length, axes, c);
            }
            auto original = c.retain(master);
            math.cancel();
            return NastranRigidAttachment{deck->source,      card,        std::move(original),
                                          std::move(points), math.info(), work.used};
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

core::Result<NastranRigidProjection> evaluate_nastran_rigid_attachment(
    const NastranRigidAttachment& input, const NastranRigidMotion& state,
    std::span<const NastranRigidPointLoad> loads, const NastranRigidAttachmentLimits& limits,
    core::ResourcePlan& plan, std::optional<NastranAttachmentDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    NastranAttachmentDiagnostic active;
    auto result = g::partition_entry([&]() -> core::Result<NastranRigidProjection> {
        try {
            ReferenceArithmetic math{limits.arithmetic, plan};
            m::Work work{0, limits.maximum_directed_work, limits.arithmetic.cancellation};
            const Calculation c{math, work, active};
            math.cancel();
            const auto owner = input.view();
            if (!owner) {
                partition_fail(Code::InvalidInput);
            }
            active.link = owner->card.id;
            if (owner->dependents.size() > limits.maximum_dependents) {
                partition_fail(Code::ResourceExceeded);
            }
            verify_loads(*owner, loads, c);
            math.admit(g::partition_add(64, g::partition_mul(42, loads.size())));
            active.reason = Reason::Pose;
            const auto motion = import_motion(state, c);
            auto rows = partition_take(plan.allocate<NastranRigidProjectedPoint>(
                g::partition_index(loads.size()), operation));
            Wrench total;
            for (std::size_t i = 0; i < loads.size(); ++i) {
                math.visit(1);
                c.grid_site(Reason::DisplacementAxes, owner->dependents[i].source);
                rows.values()[i] = projection(owner->dependents[i], loads[i],
                                              owner->card.dependent_components, motion, total, c);
            }
            auto force = c.retain(total.force), moment = c.retain(total.moment);
            math.cancel();
            return NastranRigidProjection{owner->source,    owner->card.id,    std::move(rows),
                                          std::move(force), std::move(moment), math.info(),
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
