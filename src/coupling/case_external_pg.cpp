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
    auto match = imported::qualify_external_sidecar_planar_pg(
        si, bound, shell, quantities, reference_limits, external_limits, workspace, plan,
        reference_limits.cancellation, diagnostic);
    if (!match)
        return std::unexpected(match.error());
    return CaseExternalPgFactory::make(std::move(match->reference), std::move(match->match));
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
    auto match = imported::certify_external_sidecar_pg_bounds(
        si, bound, shell, reference_limits, external_limits, maximum_directed_work, workspace, plan,
        reference_limits.cancellation, diagnostic);
    if (!match)
        return std::unexpected(match.error());
    return CaseExternalPgFactory::make(std::move(match->reference), std::move(match->match));
}
} // namespace uvlas::coupling::detail
