#pragma once

#include "case_pg_reference.hpp"
#include "io/external_mesh_sidecar_pg.hpp"

namespace uvlas::coupling::detail {

template <class Match>
struct CaseExternalPg {
    structure::PgReferenceGeometry reference;
    Match match;
};

/// Prepare the Case shell's own PG reference and verify exact planar correspondence and
/// applicable quantities. The Nastran property owner must share the authenticated deck.
[[nodiscard]] core::Result<CaseExternalPg<io::detail::QualifiedExternalSidecarPg>>
prepare_case_external_planar_pg(const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&,
                                io::detail::MeshShellId, const io::detail::NastranLaminateCards*,
                                const CasePgReferenceLimits&,
                                const io::detail::ExternalSidecarPgLimits&, core::Workspace&,
                                core::ResourcePlan&,
                                io::detail::ExternalSidecarPgDiagnostic* = nullptr) noexcept;

/// Prepare the same reference for whole-face bounds against an explicit conversion budget.
/// The retained observation IDs identify the remaining independent physical checks.
[[nodiscard]] core::Result<CaseExternalPg<io::detail::CertifiedExternalSidecarPgBounds>>
prepare_case_external_pg_bounds(const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&,
                                io::detail::MeshShellId, const CasePgReferenceLimits&,
                                const io::detail::ExternalSidecarPgLimits&,
                                std::uint64_t maximum_directed_work, core::Workspace&,
                                core::ResourcePlan&,
                                io::detail::ExternalSidecarPgDiagnostic* = nullptr) noexcept;

} // namespace uvlas::coupling::detail
