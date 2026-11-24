#include "case_external_pg.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::coupling::detail {
namespace imported = io::detail;

struct CaseExternalPgFactory {
    template <class Match>
    static CaseExternalPg<Match> make(structure::PgReferenceGeometry&& reference,
                                      Match&& match) noexcept {
        return CaseExternalPg<Match>{std::move(reference), std::forward<Match>(match)};
    }

    static CaseExternalPg<CaseNastranPgMatch>
    with_material_axes(CaseExternalPg<imported::QualifiedExternalSidecarPg>&& planar,
                       imported::NastranPgMaterialAxes&& axes) noexcept {
        return make(std::move(planar.reference_),
                    CaseNastranPgMatch{std::move(planar.match_), std::move(axes)});
    }
};

core::Result<CaseExternalPg<imported::QualifiedExternalSidecarPg>> prepare_case_external_planar_pg(
    const io::SiCase& si, const imported::ExternalMeshSidecarBinding& bound,
    imported::MeshShellId shell, const imported::NastranLaminateCards* quantities,
    const CasePgReferenceLimits& reference_limits,
    const imported::ExternalSidecarPgLimits& external_limits, core::Workspace& workspace,
    core::ResourcePlan& plan, imported::ExternalSidecarPgDiagnostic* diagnostic) noexcept {
    if (diagnostic)
        *diagnostic = {};
    const auto source = si.source_case();
    if (!source)
        return std::unexpected(source.error());
    const auto& schema = source->get();
    const auto root = schema.root();
    if (!root)
        return std::unexpected(root.error());
    auto reference =
        prepare_case_pg_reference(si, shell.value(), reference_limits, workspace, plan);
    if (!reference)
        return std::unexpected(reference.error());
    const auto unit_fields = root->units.record();
    const io::SiUnits units{static_cast<io::SiLengthUnit>(unit_fields.length.selected()),
                            static_cast<io::SiMassUnit>(unit_fields.mass.selected()),
                            static_cast<io::SiTimeUnit>(unit_fields.time.selected()),
                            static_cast<io::SiAngleUnit>(unit_fields.angle.selected())};
    auto match = imported::qualify_external_sidecar_planar_pg(
        schema, bound, shell, *reference, units, quantities, external_limits, plan,
        reference_limits.cancellation, diagnostic);
    if (!match)
        return std::unexpected(match.error());
    return CaseExternalPgFactory::make(std::move(*reference), std::move(*match));
}

core::Result<CaseExternalPg<CaseNastranPgMatch>> prepare_case_nastran_planar_pg(
    const io::SiCase& si, const imported::ExternalMeshSidecarBinding& bound,
    imported::MeshShellId shell, const imported::NastranLaminateCards& quantities,
    const CasePgReferenceLimits& reference_limits,
    const imported::ExternalSidecarPgLimits& external_limits, std::uint64_t maximum_axis_matches,
    core::Workspace& workspace, core::ResourcePlan& plan,
    CaseNastranPgDiagnostic* diagnostic) noexcept {
    if (diagnostic)
        *diagnostic = {};
    if (!std::holds_alternative<imported::NastranShellCards>(bound.source())) {
        return std::unexpected(
            core::Failure{core::FailureCode::Unsupported, core::Operation::PrepareModel});
    }
    auto planar = prepare_case_external_planar_pg(si, bound, shell, &quantities, reference_limits,
                                                  external_limits, workspace, plan,
                                                  diagnostic ? &diagnostic->planar : nullptr);
    if (!planar)
        return std::unexpected(planar.error());

    imported::NastranPgAxisLimits axis_limits{external_limits.geometry.arithmetic,
                                              maximum_axis_matches};
    axis_limits.arithmetic.cancellation = reference_limits.cancellation;
    auto axes =
        imported::match_nastran_pg_material_axes(si, bound, planar->match(), axis_limits, plan,
                                                 diagnostic ? &diagnostic->material_axes : nullptr);
    if (!axes)
        return std::unexpected(axes.error());
    return CaseExternalPgFactory::with_material_axes(std::move(*planar), std::move(*axes));
}

core::Result<std::reference_wrapper<const CaseNastranPlanarPgSet::Shell>>
CaseNastranPlanarPgSet::at(core::Index index) const& noexcept {
    if (index < 0 || std::cmp_greater_equal(index, shells_.size())) {
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    }
    if (const auto& member = shells_.values()[static_cast<std::size_t>(index)]) {
        return std::cref(*member);
    }
    return std::unexpected(
        core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
}

core::Result<CaseNastranPlanarPgSet> CaseNastranPlanarPgSet::prepare(
    const io::SiCase& si, const imported::ExternalMeshSidecarBinding& bound,
    const imported::NastranLaminateCards& quantities, const CasePgReferenceLimits& reference_limits,
    const imported::ExternalSidecarPgLimits& external_limits, const CaseNastranPgSetLimits& limits,
    core::Workspace& workspace, core::ResourcePlan& plan,
    CaseNastranPgSetDiagnostic* diagnostic) noexcept {
    using Code = core::FailureCode;
    constexpr auto operation = core::Operation::PrepareModel;
    if (diagnostic) {
        *diagnostic = {};
    }
    const auto& stop = reference_limits.cancellation;
    if (stop.stop_requested()) {
        return std::unexpected(core::Failure{Code::Cancelled, operation});
    }
    const auto* native = std::get_if<imported::NastranShellCards>(&bound.source());
    if (!native) {
        return std::unexpected(core::Failure{Code::Unsupported, operation});
    }
    const auto schema = si.source_case();
    if (!schema) {
        return std::unexpected(schema.error());
    }
    const auto identity = schema->get().identity();
    if (!identity) {
        return std::unexpected(identity.error());
    }
    const auto parent = bound.sidecar().source().parent;
    const auto native_view = native->view();
    const auto material_view = quantities.view();
    if (identity->source != parent.source ||
        identity->contract_revision != parent.contract_revision ||
        identity->representation_revision != parent.representation_revision || !native_view ||
        !material_view) {
        return std::unexpected(core::Failure{Code::StaleGeneration, operation});
    }
    const auto& material_bulk = material_view->input.input.input.bulk;
    if (!native_view->input.bulk.same_snapshot(material_bulk.location, material_bulk.digest)) {
        return std::unexpected(core::Failure{Code::StaleGeneration, operation});
    }

    const auto faces = bound.faces();
    if (faces.empty()) {
        return std::unexpected(core::Failure{Code::InvalidInput, operation});
    }
    const auto capacity = std::min<std::uint64_t>(faces.size(), limits.maximum_shells);
    if (!std::in_range<core::Index>(capacity)) {
        return std::unexpected(core::Failure{Code::ResourceExceeded, operation});
    }
    auto selections = plan.allocate<std::size_t>(static_cast<core::Index>(capacity), operation);
    if (!selections) {
        return std::unexpected(selections.error());
    }

    std::size_t count = 0;
    auto work = limits.maximum_selection_work;
    for (std::size_t i = 0; i < faces.size(); ++i) {
        if (stop.stop_requested()) {
            return std::unexpected(core::Failure{Code::Cancelled, operation});
        }
        if (work == 0) {
            return std::unexpected(core::Failure{Code::ResourceExceeded, operation});
        }
        --work;
        const auto& face = faces[i];
        bool selected = false;
        for (const auto previous : selections->values().first(count)) {
            if (stop.stop_requested()) {
                return std::unexpected(core::Failure{Code::Cancelled, operation});
            }
            if (work == 0) {
                return std::unexpected(core::Failure{Code::ResourceExceeded, operation});
            }
            --work;
            if (faces[previous].shell == face.shell) {
                selected = true;
                break;
            }
        }
        if (!selected) {
            if (count == selections->size()) {
                if (diagnostic) {
                    diagnostic->shell = face.shell;
                    diagnostic->first_declaration = face.sidecar_node;
                }
                return std::unexpected(core::Failure{Code::ResourceExceeded, operation});
            }
            selections->values()[count++] = i;
        }
    }
    auto result = plan.allocate<std::optional<Shell>>(static_cast<core::Index>(count), operation);
    if (!result) {
        return std::unexpected(result.error());
    }

    for (std::size_t i = 0; i < count; ++i) {
        const auto& face = faces[selections->values()[i]];
        auto member = prepare_case_nastran_planar_pg(
            si, bound, face.shell, quantities, reference_limits, external_limits,
            limits.maximum_axis_matches_per_shell, workspace, plan,
            diagnostic ? &diagnostic->member : nullptr);
        if (!member) {
            if (diagnostic) {
                diagnostic->shell = face.shell;
                diagnostic->first_declaration = face.sidecar_node;
            }
            return std::unexpected(member.error());
        }
        result->values()[i].emplace(std::move(*member));
    }
    if (stop.stop_requested()) {
        return std::unexpected(core::Failure{Code::Cancelled, operation});
    }
    return CaseNastranPlanarPgSet{std::move(*result)};
}

core::Result<CaseExternalPg<imported::CertifiedExternalSidecarPgBounds>>
prepare_case_external_pg_bounds(const io::SiCase& si,
                                const imported::ExternalMeshSidecarBinding& bound,
                                imported::MeshShellId shell,
                                const CasePgReferenceLimits& reference_limits,
                                const imported::ExternalSidecarPgLimits& external_limits,
                                std::uint64_t maximum_directed_work, core::Workspace& workspace,
                                core::ResourcePlan& plan,
                                imported::ExternalSidecarPgDiagnostic* diagnostic) noexcept {
    if (diagnostic)
        *diagnostic = {};
    auto reference =
        prepare_case_pg_reference(si, shell.value(), reference_limits, workspace, plan);
    if (!reference)
        return std::unexpected(reference.error());
    auto match = imported::certify_external_sidecar_pg_bounds(
        si, bound, shell, *reference, external_limits, maximum_directed_work, plan,
        reference_limits.cancellation, diagnostic);
    if (!match)
        return std::unexpected(match.error());
    return CaseExternalPgFactory::make(std::move(*reference), std::move(*match));
}
} // namespace uvlas::coupling::detail
