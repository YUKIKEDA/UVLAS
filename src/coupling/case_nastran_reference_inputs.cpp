#include "case_nastran_reference_inputs.hpp"

#include "io/si_case_detail.hpp"
#include "structure/pg_reference_polynomials.hpp"

#include <algorithm>
#include <tuple>
#include <type_traits>
#include <utility>

namespace uvlas::coupling::detail {
namespace {

namespace native = io::detail;
using core::Failure;
using core::FailureCode;
constexpr auto operation = core::Operation::PrepareModel;

template <class T, class D>
void record_failure(const core::Result<T>& result, const std::optional<D>& detail,
                    CaseNastranReferenceDiagnostic* diagnostic) {
    if (!result && detail && diagnostic) {
        diagnostic->card = *detail;
    }
}

template <class Owned, io::SchemaNamespace Space, class Prepare>
[[nodiscard]] core::Result<core::PlannedBuffer<std::optional<Owned>>>
prepare_rows(std::span<const native::NastranSidecarElementBinding<Space>> rows,
             const std::stop_token& stop, core::ResourcePlan& plan,
             CaseNastranReferenceDiagnostic* diagnostic, Prepare prepare) {
    if (!std::in_range<core::Index>(rows.size())) {
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
    }
    auto result =
        plan.allocate<std::optional<Owned>>(static_cast<core::Index>(rows.size()), operation);
    if (!result) {
        return std::unexpected(result.error());
    }
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (stop.stop_requested()) {
            return std::unexpected(Failure{FailureCode::Cancelled, operation});
        }
        const auto& row = rows[i];
        auto physical = prepare(row.element);
        if (!physical) {
            if (diagnostic) {
                diagnostic->selection =
                    CaseNastranCardSelection{Space, row.target.value(), row.element, row.node};
            }
            return std::unexpected(physical.error());
        }
        result->values()[i].emplace(Owned{row, std::move(*physical)});
    }
    return result;
}

[[nodiscard]] auto
prepare_beams(const native::NastranSidecar& rows, const native::NastranShellCards& deck,
              const native::NastranIsotropicShellCards& isotropic,
              const CaseNastranReferenceInputLimits& limits, const std::stop_token& stop,
              core::ResourcePlan& plan, CaseNastranReferenceDiagnostic* diagnostic)
    -> core::Result<core::PlannedBuffer<std::optional<CaseNastranReferenceInputs::Beam>>> {
    if (rows.beams.size() == 0) {
        return core::PlannedBuffer<std::optional<CaseNastranReferenceInputs::Beam>>{};
    }
    std::optional<native::NastranBeamDiagnostic> cards_detail;
    const auto cards =
        native::NastranBeamCards::prepare(deck, limits.beam_cards, plan, stop, &cards_detail);
    if (!cards) {
        record_failure(cards, cards_detail, diagnostic);
        return std::unexpected(cards.error());
    }
    std::optional<native::NastranPbeamDiagnostic> properties_detail;
    const auto properties = native::NastranPbeamCards::prepare(
        isotropic, *cards, limits.beam_properties, plan, stop, &properties_detail);
    if (!properties) {
        record_failure(properties, properties_detail, diagnostic);
        return std::unexpected(properties.error());
    }
    auto limits_si = limits.beam;
    limits_si.arithmetic.cancellation = stop;
    return prepare_rows<CaseNastranReferenceInputs::Beam>(
        rows.beams.values(), stop, plan, diagnostic, [&](native::NastranElementId element) {
            std::optional<native::NastranBeamSiDiagnostic> detail;
            auto result = prepare_nastran_si_beam(isotropic, *cards, *properties, element,
                                                  rows.units, limits_si, plan, &detail);
            record_failure(result, detail, diagnostic);
            return result;
        });
}

[[nodiscard]] auto
prepare_masses(const native::NastranSidecar& rows, const native::NastranShellCards& deck,
               const CaseNastranReferenceInputLimits& limits, const std::stop_token& stop,
               core::ResourcePlan& plan, CaseNastranReferenceDiagnostic* diagnostic)
    -> core::Result<core::PlannedBuffer<std::optional<CaseNastranReferenceInputs::Mass>>> {
    if (rows.masses.size() == 0) {
        return core::PlannedBuffer<std::optional<CaseNastranReferenceInputs::Mass>>{};
    }
    std::optional<native::NastranLumpedMassDiagnostic> cards_detail;
    const auto cards =
        native::NastranLumpedMassCards::prepare(deck, limits.mass_cards, plan, stop, &cards_detail);
    if (!cards) {
        record_failure(cards, cards_detail, diagnostic);
        return std::unexpected(cards.error());
    }
    auto limits_si = limits.mass;
    limits_si.arithmetic.cancellation = stop;
    return prepare_rows<CaseNastranReferenceInputs::Mass>(
        rows.masses.values(), stop, plan, diagnostic, [&](native::NastranElementId element) {
            std::optional<native::NastranMassGeometryDiagnostic> detail;
            auto result = prepare_nastran_basic_mass(deck, *cards, element, rows.units, limits_si,
                                                     plan, &detail);
            record_failure(result, detail, diagnostic);
            return result;
        });
}

[[nodiscard]] auto
prepare_joints(const native::NastranSidecar& rows, const native::NastranShellCards& deck,
               const CaseNastranReferenceInputLimits& limits, const std::stop_token& stop,
               core::ResourcePlan& plan, CaseNastranReferenceDiagnostic* diagnostic)
    -> core::Result<core::PlannedBuffer<std::optional<CaseNastranReferenceInputs::Joint>>> {
    if (rows.joints.size() == 0) {
        return core::PlannedBuffer<std::optional<CaseNastranReferenceInputs::Joint>>{};
    }
    std::optional<native::NastranRigidLinkDiagnostic> cards_detail;
    const auto cards =
        native::NastranRigidLinkCards::prepare(deck, limits.rigid_cards, plan, stop, &cards_detail);
    if (!cards) {
        record_failure(cards, cards_detail, diagnostic);
        return std::unexpected(cards.error());
    }
    auto limits_si = limits.rigid;
    limits_si.arithmetic.cancellation = stop;
    return prepare_rows<CaseNastranReferenceInputs::Joint>(
        rows.joints.values(), stop, plan, diagnostic, [&](native::NastranElementId element) {
            std::optional<native::NastranAttachmentDiagnostic> detail;
            auto result = native::NastranRigidAttachment::prepare(
                deck, *cards, element, rows.units.length, limits_si, plan, &detail);
            record_failure(result, detail, diagnostic);
            return result;
        });
}

} // namespace

core::Result<CaseNastranReferenceInputs> CaseNastranReferenceInputs::prepare(
    const io::SiCase& si, const native::ExternalMeshSidecarBinding& bound,
    const native::NastranLaminateCards& quantities, const CasePgReferenceLimits& reference_limits,
    const native::ExternalSidecarPgLimits& external_limits,
    const CaseNastranReferenceInputLimits& limits, core::Workspace& workspace,
    core::ResourcePlan& plan, CaseNastranReferenceDiagnostic* diagnostic) noexcept {
    if (diagnostic) {
        *diagnostic = {};
    }
    const auto& stop = reference_limits.cancellation;
    if (stop.stop_requested()) {
        return std::unexpected(Failure{FailureCode::Cancelled, operation});
    }
    const auto* rows = std::get_if<native::NastranSidecar>(&bound.sidecar().payload());
    const auto* deck = std::get_if<native::NastranShellCards>(&bound.source());
    if (!rows || !deck) {
        return std::unexpected(Failure{FailureCode::Unsupported, operation});
    }
    const auto schema = si.source_case();
    if (!schema) {
        return std::unexpected(schema.error());
    }
    const auto identity = schema->get().identity();
    if (!identity) {
        return std::unexpected(identity.error());
    }
    const auto source = bound.sidecar().source();
    const auto& parent = source.parent;
    if (identity->source != parent.source ||
        identity->contract_revision != parent.contract_revision ||
        identity->representation_revision != parent.representation_revision) {
        return std::unexpected(Failure{FailureCode::StaleGeneration, operation});
    }
    const auto original = deck->view();
    const auto qualified = quantities.view();
    if (!original || !qualified ||
        !original->input.bulk.same_snapshot(qualified->input.input.input.bulk.location,
                                            qualified->input.input.input.bulk.digest)) {
        return std::unexpected(Failure{FailureCode::StaleGeneration, operation});
    }
    auto remaining = limits.maximum_bound_cards;
    for (const auto count : {rows->beams.size(), rows->masses.size(), rows->joints.size()}) {
        if (count > remaining) {
            return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
        }
        remaining -= count;
    }

    CaseNastranReferenceInputs result{source, bound.sidecar().digest(),
                                      original->input.bulk.location, original->input.bulk.digest};
    if (!bound.faces().empty()) {
        auto shells = CaseNastranPlanarPgSet::prepare(
            si, bound, quantities, reference_limits, external_limits, limits.shells, workspace,
            plan, diagnostic ? &diagnostic->shells : nullptr);
        if (!shells) {
            return std::unexpected(shells.error());
        }
        result.shells_.emplace(std::move(*shells));
    }

    auto beams =
        prepare_beams(*rows, *deck, quantities.isotropic(), limits, stop, plan, diagnostic);
    if (!beams) {
        return std::unexpected(beams.error());
    }
    auto masses = prepare_masses(*rows, *deck, limits, stop, plan, diagnostic);
    if (!masses) {
        return std::unexpected(masses.error());
    }
    auto joints = prepare_joints(*rows, *deck, limits, stop, plan, diagnostic);
    if (!joints) {
        return std::unexpected(joints.error());
    }
    if (stop.stop_requested()) {
        return std::unexpected(Failure{FailureCode::Cancelled, operation});
    }
    result.beams_ = std::move(*beams);
    result.masses_ = std::move(*masses);
    result.joints_ = std::move(*joints);
    return result;
}

namespace {

namespace ref = structure::detail;
using MassInput = CaseNastranReferenceInputs::Mass;

template <class CaseMass>
using MassEvidence = std::conditional_t<std::is_same_v<CaseMass, io::schema::MassRigid>,
                                        CaseNastranRigidMass, CaseNastranParticleMass>;

template <class CaseMass>
core::Result<MassEvidence<CaseMass>>
certify_case_mass(const io::SiCase& si, const CaseNastranReferenceInputs& prepared,
                  CaseNastranMassId id, ref::ReferenceArithmetic& math, bool admitted,
                  CaseNastranParticleDiagnostic* diagnostic, const MassInput* selected = nullptr,
                  std::optional<core::Index> case_ordinal = {}) noexcept {
    using ref::partition_take;
    using Reason = CaseNastranParticleFailure;
    using Result = core::Result<MassEvidence<CaseMass>>;
    constexpr bool rigid = std::is_same_v<CaseMass, io::schema::MassRigid>;
    if (diagnostic)
        *diagnostic = {};
    return ref::partition_entry([&]() -> Result {
        CaseNastranParticleDiagnostic detail{};
        const auto mark = [&](Reason reason,
                              std::optional<native::NastranBulkPosition> field = {}) {
            detail.reason = reason;
            detail.field = field;
        };
        const auto reject = [&](FailureCode code) -> Result {
            if (diagnostic)
                *diagnostic = detail;
            return std::unexpected(Failure{code, operation});
        };
        math.cancel();
        const auto& schema = partition_take(si.source_case()).get();
        const auto identity = partition_take(schema.identity());
        const auto sidecar = prepared.sidecar_source();
        const auto deck = prepared.deck_source();
        if (!native::same_identity(identity, sidecar.parent) ||
            !native::same_identity(identity, deck.parent))
            return reject(FailureCode::StaleGeneration);

        mark(Reason::Selection);
        if (!selected) {
            for (const auto& slot : prepared.masses()) {
                math.visit(1);
                if (!slot || slot->binding.target != id)
                    continue;
                detail.element = slot->binding.element;
                if (selected)
                    return reject(FailureCode::InvalidInput);
                selected = &*slot;
            }
        }
        if (!selected)
            return reject(FailureCode::InvalidInput);
        detail.element = selected->binding.element;
        const auto& original = selected->physical;
        if (!native::same_identity(identity, original.source.location.parent) ||
            original.source.location.source != deck.source ||
            original.source.digest != prepared.deck_digest()) {
            mark(Reason::Source);
            return reject(FailureCode::StaleGeneration);
        }

        mark(Reason::Kind, original.card.id_field.position);
        constexpr auto required_kind =
            rigid ? native::NastranBasicMassKind::Rigid : native::NastranBasicMassKind::Particle;
        if (original.kind != required_kind)
            return reject(FailureCode::InvalidInput);
        if constexpr (!rigid) {
            for (std::size_t axis = 0; axis < original.grid_to_centre_metres.size(); ++axis) {
                math.visit(1);
                const auto& lever = original.grid_to_centre_metres[axis];
                if (lever.lower.sign() == 0 && lever.upper.sign() == 0)
                    continue;
                const auto& field = original.card.position[axis].field;
                mark(Reason::Kind, field ? field->position : original.card.grid_field.position);
                return reject(lever.lower.sign() > 0 || lever.upper.sign() < 0
                                  ? FailureCode::InvalidInput
                                  : FailureCode::PrecisionUnresolved);
            }
        }

        const auto model = partition_take(schema.root()).model.record();
        std::optional<CaseMass> case_mass;
        const auto read_case_mass = [&](core::Index i) {
            math.visit(1);
            const auto value = partition_take(model.masses.at(i)).record();
            if (const auto* candidate = std::get_if<CaseMass>(&value);
                candidate && candidate->id.unsigned_value() == id.value())
                case_mass = *candidate;
        };
        if (case_ordinal) {
            if (*case_ordinal >= 0)
                read_case_mass(*case_ordinal);
        } else {
            for (std::size_t i = 0; i < model.masses.size() && !case_mass; ++i)
                read_case_mass(ref::partition_index(i));
        }
        if (!case_mass)
            return reject(FailureCode::InvalidInput);

        const auto point = case_mass->target.record();
        const auto member = point.target.record();
        mark(Reason::Position, original.card.grid_field.position);
        if (member.kind.selected() != decltype(member.kind.selected())::ItemShell ||
            point.side.selected() != decltype(point.side.selected())::ItemInterior)
            return reject(FailureCode::Unsupported);
        const native::MeshShellId shell{member.id.unsigned_value()};
        const auto& set = prepared.shells();
        if (!set)
            return reject(FailureCode::Unsupported);
        const auto find_shell = [&]() -> const CaseNastranPlanarPgSet::Shell* {
            for (std::size_t i = 0; i < set->size(); ++i) {
                math.visit(1);
                const auto& item = partition_take(set->at(ref::partition_index(i))).get();
                if (item.match().planar.shell == shell)
                    return &item;
            }
            return nullptr;
        };
        const auto* certified = find_shell();
        if (!certified)
            return reject(FailureCode::Unsupported);
        const auto& match = certified->match().planar;
        const auto geometry = match.geometry.view();
        const auto patch = partition_take(certified->reference().identity());
        if (!geometry || geometry->patch != patch ||
            !native::same_identity(identity, geometry->source.parent) ||
            geometry->source.source != deck.source || geometry->digest != prepared.deck_digest() ||
            match.sidecar.source != sidecar.source ||
            match.sidecar_digest != prepared.sidecar_digest()) {
            mark(Reason::Source);
            return reject(FailureCode::StaleGeneration);
        }

        if (!admitted)
            math.admit();
        const auto exact = [&](const auto& field) {
            math.charge();
            return partition_take(io::materialize_rational(
                partition_take(si.quantity(field)).recipe, math.leaf_limits, math.plan));
        };
        mark(Reason::Mass,
             original.card.mass.field.transform([](const auto& f) { return f.position; }));
        auto mass = exact(case_mass->mass);
        if (math.compare(mass, original.mass_kg) != 0)
            return reject(FailureCode::InvalidInput);

        mark(Reason::Position, original.card.grid_field.position);
        if (exact(partition_take(point.coordinates.at(2))).sign() != 0)
            return reject(FailureCode::InvalidInput);
        const auto s = exact(partition_take(point.coordinates.at(0)));
        const auto t = exact(partition_take(point.coordinates.at(1)));
        for (const auto& face : geometry->faces) {
            for (const auto& corner : face.ring()) {
                math.visit(1);
                if (corner.node != original.card.grid.value() ||
                    math.compare(corner.chart.s, s) != 0 || math.compare(corner.chart.t, t) != 0)
                    continue;
                if (!corner.nastran) {
                    mark(Reason::Source);
                    return reject(FailureCode::StaleGeneration);
                }
                if (const auto &unshifted = corner.nastran->unshifted_metres;
                    unshifted && !std::ranges::equal(*unshifted, corner.position_metres,
                                                     [&](const auto& a, const auto& b) {
                                                         return math.compare(a, b) == 0;
                                                     }))
                    return reject(FailureCode::InvalidInput);
                CaseNastranMassAttachment attachment{deck,
                                                     sidecar,
                                                     prepared.deck_digest(),
                                                     prepared.sidecar_digest(),
                                                     selected->binding,
                                                     shell,
                                                     patch,
                                                     original.card.grid,
                                                     face.element,
                                                     std::move(mass)};
                if constexpr (!rigid) {
                    return attachment;
                } else {
                    using R = structure::ExactRational;
                    mark(Reason::Frame, original.card.grid_field.position);
                    const auto frame = case_mass->frame.record();
                    std::array<R, 3> origin, offset;
                    std::array<R, 9> rotation, inertia;
                    for (std::size_t i = 0; i < 3; ++i) {
                        const auto row = ref::partition_index(i);
                        origin[i] = exact(partition_take(frame.origin.at(row)));
                        if (math.compare(origin[i], corner.position_metres[i]) != 0)
                            return reject(FailureCode::InvalidInput);
                        offset[i] = exact(partition_take(case_mass->center.at(row)));
                        const auto rotation_row = partition_take(frame.rotation.at(row));
                        const auto inertia_row = partition_take(case_mass->inertia.at(row));
                        for (std::size_t j = 0; j < 3; ++j) {
                            const auto column = ref::partition_index(j);
                            rotation[3 * i + j] = exact(partition_take(rotation_row.at(column)));
                            inertia[3 * i + j] = exact(partition_take(inertia_row.at(column)));
                        }
                    }
                    const auto frame_valid = ref::partition_entry([&]() -> core::Result<void> {
                        const auto one = math.constant("1");
                        ref::validate_frame(math, {rotation, origin, one});
                        return {};
                    });
                    if (!frame_valid)
                        return reject(frame_valid.error().code);

                    std::array<R, 3> centre;
                    const bool absolute =
                        std::holds_alternative<native::NastranBasicMassCentre>(original.card.frame);
                    for (std::size_t i = 0; i < 3; ++i) {
                        const auto& field = original.card.position[i].field;
                        mark(Reason::Centre,
                             field ? field->position : original.card.grid_field.position);
                        auto displacement = R{};
                        for (std::size_t j = 0; j < 3; ++j)
                            displacement =
                                math.add(displacement, math.mul(rotation[3 * i + j], offset[j]));
                        centre[i] = math.add(origin[i], displacement);
                        const auto& field_value = original.card_position_metres[i];
                        const auto relation =
                            absolute ? math.compare(centre[i], field_value)
                                     : math.compare(centre[i], math.add(origin[i], field_value));
                        if (relation != 0)
                            return reject(FailureCode::InvalidInput);
                    }

                    const auto& original_inertia = original.card.inertia;
                    const std::array inertia_fields{
                        original_inertia.i11.field, original_inertia.i21.field,
                        original_inertia.i22.field, original_inertia.i31.field,
                        original_inertia.i32.field, original_inertia.i33.field};
                    std::array<R, 6> centre_inertia;
                    std::size_t packed = 0;
                    for (std::size_t i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j <= i; ++j, ++packed) {
                            const auto& field = inertia_fields[packed];
                            mark(Reason::Inertia,
                                 field ? field->position : original.card.card.position);
                            if (math.compare(inertia[3 * i + j], inertia[3 * j + i]) != 0)
                                return reject(FailureCode::InvalidInput);
                            auto& value = centre_inertia[packed];
                            for (std::size_t k = 0; k < 3; ++k) {
                                for (std::size_t l = 0; l < 3; ++l) {
                                    value = math.add(value, math.mul(math.mul(rotation[3 * i + k],
                                                                              inertia[3 * k + l]),
                                                                     rotation[3 * j + l]));
                                }
                            }
                            if (math.compare(value, original.centre_inertia_kg_m2[packed]) != 0)
                                return reject(FailureCode::InvalidInput);
                        }
                    }
                    return CaseNastranRigidMass{std::move(attachment), std::move(centre),
                                                std::move(centre_inertia)};
                }
            }
        }
        return reject(FailureCode::InvalidInput);
    });
}

} // namespace

core::Result<CaseNastranParticleMass>
certify_case_nastran_particle_mass(const io::SiCase& si, const CaseNastranReferenceInputs& prepared,
                                   CaseNastranMassId id, const structure::PgReferenceLimits& limits,
                                   core::ResourcePlan& plan,
                                   CaseNastranParticleDiagnostic* diagnostic) noexcept {
    ref::ReferenceArithmetic math{limits, plan};
    return certify_case_mass<io::schema::MassParticle>(si, prepared, id, math, false, diagnostic);
}

core::Result<CaseNastranRigidMass>
certify_case_nastran_rigid_mass(const io::SiCase& si, const CaseNastranReferenceInputs& prepared,
                                CaseNastranMassId id, const structure::PgReferenceLimits& limits,
                                core::ResourcePlan& plan,
                                CaseNastranParticleDiagnostic* diagnostic) noexcept {
    ref::ReferenceArithmetic math{limits, plan};
    return certify_case_mass<io::schema::MassRigid>(si, prepared, id, math, false, diagnostic);
}

core::Result<CaseNastranMassSet>
CaseNastranMassSet::certify(const io::SiCase& si, const CaseNastranReferenceInputs& prepared,
                            const structure::PgReferenceLimits& limits, core::ResourcePlan& plan,
                            std::optional<CaseNastranMassSetDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return ref::partition_entry([&]() -> core::Result<CaseNastranMassSet> {
        ref::ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto& schema = ref::partition_take(si.source_case()).get();
        const auto identity = ref::partition_take(schema.identity());
        if (!native::same_identity(identity, prepared.sidecar_source().parent) ||
            !native::same_identity(identity, prepared.deck_source().parent))
            ref::partition_fail(FailureCode::StaleGeneration);

        const auto rows = prepared.masses();
        if (rows.empty())
            return CaseNastranMassSet{core::PlannedBuffer<std::optional<Mass>>{}};
        auto live_rationals = ref::reference_live_rationals;
        for (const auto& row : rows) {
            math.visit(1);
            if (!row)
                ref::partition_fail(FailureCode::InvalidInput);
            // A particle retains mass. A rigid mass also retains three centre
            // coordinates and six centroidal inertia components.
            const auto retained =
                row->physical.kind == native::NastranBasicMassKind::Rigid ? 10U : 1U;
            live_rationals = ref::partition_add(live_rationals, retained);
        }
        math.admit(live_rationals);
        struct MassIndex {
            std::uint64_t id{};
            core::Index ordinal{};
        };
        const auto compares_id = [&](std::uint64_t left, std::uint64_t right) {
            math.visit(1);
            return left < right;
        };
        // Sort indices only. The owning evidence is built in sidecar declaration order.
        {
            auto bindings = ref::partition_take(
                plan.allocate<MassIndex>(ref::partition_index(rows.size()), operation));
            auto order = bindings.values();
            for (std::size_t i = 0; i < rows.size(); ++i) {
                math.visit(1);
                order[i] = {rows[i]->binding.target.value(), ref::partition_index(i)};
            }
            std::ranges::sort(order, [&](const MassIndex& left, const MassIndex& right) {
                math.visit(1);
                return std::tie(left.id, left.ordinal) < std::tie(right.id, right.ordinal);
            });
            for (std::size_t i = 1; i < order.size(); ++i) {
                math.visit(1);
                if (order[i - 1].id != order[i].id)
                    continue;
                const auto ordinal = order[i].ordinal;
                const auto& duplicate = rows[static_cast<std::size_t>(ordinal)]->binding;
                if (diagnostic) {
                    CaseNastranParticleDiagnostic detail;
                    detail.reason = CaseNastranParticleFailure::Selection;
                    detail.element = duplicate.element;
                    *diagnostic = CaseNastranMassSetDiagnostic{{io::SchemaNamespace::Mass,
                                                                duplicate.target.value(),
                                                                duplicate.element, duplicate.node},
                                                               detail};
                }
                ref::partition_fail(FailureCode::InvalidInput);
            }
        }
        const auto model = ref::partition_take(schema.root()).model.record();
        auto cases = ref::partition_take(
            plan.allocate<MassIndex>(ref::partition_index(model.masses.size()), operation));
        auto case_order = cases.values();
        for (std::size_t i = 0; i < case_order.size(); ++i) {
            math.visit(1);
            const auto ordinal = ref::partition_index(i);
            const auto value = ref::partition_take(model.masses.at(ordinal)).record();
            const auto id =
                std::visit([](const auto& item) { return item.id.unsigned_value(); }, value);
            case_order[i] = {id, ordinal};
        }
        std::ranges::sort(case_order, compares_id, &MassIndex::id);
        auto owned = ref::partition_take(
            plan.allocate<std::optional<Mass>>(ref::partition_index(rows.size()), operation));
        const auto match = [&](const auto& row, auto tag) -> core::Result<Mass> {
            CaseNastranParticleDiagnostic detail;
            const auto& binding = row.binding;
            const auto target = binding.target.value();
            const auto found =
                std::ranges::lower_bound(case_order, target, compares_id, &MassIndex::id);
            const core::Index ordinal =
                found != case_order.end() && found->id == target ? found->ordinal : -1;
            auto qualified = certify_case_mass<typename decltype(tag)::type>(
                si, prepared, binding.target, math, true, &detail, &row, ordinal);
            if (!qualified && diagnostic) {
                *diagnostic =
                    CaseNastranMassSetDiagnostic{{io::SchemaNamespace::Mass, binding.target.value(),
                                                  binding.element, binding.node},
                                                 detail};
            }
            if (!qualified)
                return std::unexpected(qualified.error());
            return Mass{std::move(*qualified)};
        };
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto& row = rows[i].value();
            auto result = row.physical.kind == native::NastranBasicMassKind::Rigid
                              ? match(row, std::type_identity<io::schema::MassRigid>{})
                              : match(row, std::type_identity<io::schema::MassParticle>{});
            owned.values()[i].emplace(ref::partition_take(std::move(result)));
        }
        math.cancel();
        return CaseNastranMassSet{std::move(owned)};
    });
}

} // namespace uvlas::coupling::detail
