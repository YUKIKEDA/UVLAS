#include "case_nastran_reference_inputs.hpp"

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

} // namespace uvlas::coupling::detail
