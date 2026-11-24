#include "case_external_pg.hpp"

#include <utility>

namespace uvlas::coupling::detail {
namespace imported = io::detail;

struct CaseExternalPgFactory {
    template <class Match>
    static CaseExternalPg<Match> make(structure::PgReferenceGeometry&& reference,
                                      Match&& match) noexcept {
        return CaseExternalPg<Match>{std::move(reference), std::forward<Match>(match)};
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
