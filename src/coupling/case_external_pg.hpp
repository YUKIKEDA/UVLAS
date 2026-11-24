#pragma once

#include "case_pg_reference.hpp"
#include "io/nastran_pg_material_axes.hpp"

#include <optional>
#include <utility>

namespace uvlas::coupling::detail {

template <class Match>
class CaseExternalPg {
public:
    [[nodiscard]] const structure::PgReferenceGeometry& reference() const& noexcept {
        return reference_;
    }
    [[nodiscard]] const structure::PgReferenceGeometry& reference() const&& = delete;

    [[nodiscard]] const Match& match() const& noexcept {
        return match_;
    }
    [[nodiscard]] const Match& match() const&& = delete;

private:
    friend struct CaseExternalPgFactory;
    CaseExternalPg(structure::PgReferenceGeometry&& reference, Match&& match) noexcept
        : reference_(std::move(reference)), match_(std::move(match)) {}

    structure::PgReferenceGeometry reference_;
    Match match_;
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

struct CaseNastranPgMatch {
    io::detail::QualifiedExternalSidecarPg planar;
    io::detail::NastranPgMaterialAxes material_axes;
};

struct CaseNastranPgDiagnostic {
    io::detail::ExternalSidecarPgDiagnostic planar{};
    std::optional<io::detail::NastranPgAxisDiagnostic> material_axes{};
};

/// Prepare one Nastran shell with exact planar, thickness, layer and material-axis evidence.
/// The authenticated quantity owner and the Case remain readable after success or failure.
/// The Case preparation cancellation token also covers the material-axis comparison.
[[nodiscard]] core::Result<CaseExternalPg<CaseNastranPgMatch>> prepare_case_nastran_planar_pg(
    const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&, io::detail::MeshShellId,
    const io::detail::NastranLaminateCards&, const CasePgReferenceLimits&,
    const io::detail::ExternalSidecarPgLimits&, std::uint64_t maximum_axis_matches,
    core::Workspace&, core::ResourcePlan&, CaseNastranPgDiagnostic* = nullptr) noexcept;

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
