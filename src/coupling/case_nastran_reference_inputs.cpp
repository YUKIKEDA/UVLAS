#include "case_nastran_reference_inputs.hpp"

#include "io/si_case_detail.hpp"
#include "structure/pg_reference_polynomials.hpp"

#include <algorithm>
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

template <class CaseMass>
using MassEvidence = std::conditional_t<std::is_same_v<CaseMass, io::schema::MassRigid>,
                                        CaseNastranRigidMass, CaseNastranParticleMass>;

template <class CaseMass>
core::Result<MassEvidence<CaseMass>>
certify_case_mass(const io::SiCase& si, const CaseNastranReferenceInputs& prepared,
                  CaseNastranMassId id, const structure::PgReferenceLimits& limits,
                  core::ResourcePlan& plan, CaseNastranParticleDiagnostic* diagnostic) noexcept {
    namespace ref = structure::detail;
    using ref::partition_take;
    using Reason = CaseNastranParticleFailure;
    using Result = core::Result<MassEvidence<CaseMass>>;
    constexpr bool rigid = std::is_same_v<CaseMass, io::schema::MassRigid>;
    if (diagnostic)
        *diagnostic = {};
    return ref::partition_entry([&]() -> Result {
        ref::ReferenceArithmetic math{limits, plan};
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
        const CaseNastranReferenceInputs::Mass* selected = nullptr;
        for (const auto& slot : prepared.masses()) {
            math.visit(1);
            if (!slot || slot->binding.target != id)
                continue;
            detail.element = slot->binding.element;
            if (selected)
                return reject(FailureCode::InvalidInput);
            selected = &*slot;
        }
        if (!selected)
            return reject(FailureCode::InvalidInput);
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
        for (std::size_t i = 0; i < model.masses.size(); ++i) {
            math.visit(1);
            const auto value = partition_take(model.masses.at(ref::partition_index(i))).record();
            if (const auto* candidate = std::get_if<CaseMass>(&value);
                candidate && candidate->id.unsigned_value() == id.value()) {
                case_mass = *candidate;
                break;
            }
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

        math.admit();
        const auto exact = [&](const auto& field) {
            math.charge();
            return partition_take(io::materialize_rational(
                partition_take(si.quantity(field)).recipe, math.leaf_limits, plan));
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
    return certify_case_mass<io::schema::MassParticle>(si, prepared, id, limits, plan, diagnostic);
}

core::Result<CaseNastranRigidMass>
certify_case_nastran_rigid_mass(const io::SiCase& si, const CaseNastranReferenceInputs& prepared,
                                CaseNastranMassId id, const structure::PgReferenceLimits& limits,
                                core::ResourcePlan& plan,
                                CaseNastranParticleDiagnostic* diagnostic) noexcept {
    return certify_case_mass<io::schema::MassRigid>(si, prepared, id, limits, plan, diagnostic);
}

} // namespace uvlas::coupling::detail
