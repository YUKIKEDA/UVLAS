#pragma once

#include "case_external_pg.hpp"
#include "io/nastran_position_enclosures.hpp"

#include <uvlas/structure/pg_full_trace.hpp>
#include <uvlas/structure/pg_prepared_patch.hpp>

namespace uvlas::coupling::detail {

class CaseStaticPgPointPatch;

struct CaseNastranPointSupportLimits {
    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_terms{};
};

enum class CaseNastranPointSupportFailure : std::uint8_t {
    Binding,
    Support,
    Position,
    Motion,
    Components,
};

struct CaseNastranPointSupportDiagnostic {
    CaseNastranPointSupportFailure reason = CaseNastranPointSupportFailure::Binding;
    io::detail::NastranSupportSetId set{0};
    io::detail::NastranGridId grid{0};
    std::optional<core::Index> declaration{};
    std::optional<io::detail::NastranBulkPosition> field{}, cd_field{};
    std::optional<std::uint64_t> element{}, offset_byte{};
    std::optional<io::detail::NastranPositionDiagnostic> direction{};
};

struct CaseNastranPointSupportTerm {
    io::detail::NastranSupportSetId source_set{0};
    std::uint8_t local_components{};
    io::detail::NastranBulkPosition card{}, components{};
    std::optional<io::detail::NastranBulkPosition> motion{};
};

struct CaseNastranPointSupport {
    io::detail::NastranSupportOrigin deck;
    io::detail::RootSnapshotLocation sidecar;
    io::SourceDigest sidecar_digest;
    io::detail::NastranSidecarSupport mapping;
    io::detail::MeshShellId shell{0};
    structure::PgReferenceIdentity patch;
    std::uint64_t element{};
    structure::PgPartitionPoint chart;                ///< Exact (s,t) at zeta=0.
    structure::PgCartesianExactPosition world_metres; ///< Original GRID in basic SI coordinates.
    io::detail::NastranBulkPosition grid_field;
    io::detail::NastranCoordinateId cd{0};
    std::optional<io::detail::NastranBulkPosition> cd_field{};
    std::uint8_t world_components{};
    core::PlannedBuffer<CaseNastranPointSupportTerm> terms;
    structure::PgReferenceCallInfo call_info;
};

/// Exact dimensionless point evaluation coefficients in the same prepared basis.
struct CaseStaticPgPointRow {
    structure::PgPreparedPatchIdentity basis;
    io::detail::MeshSupportId support;
    io::detail::MeshShellId shell;
    std::uint8_t world_components;
    core::PlannedBuffer<structure::PgTraceSparseEntry> terms;
    structure::PgCartesianExactPosition prescribed_metres;
    structure::PgReferenceCallInfo call_info;
};

[[nodiscard]] core::Result<CaseStaticPgPointRow> assemble_case_nastran_static_point_row(
    const CaseNastranPointSupport&, const structure::PgPreparedPatch&,
    const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

/// Same static centre-surface point row when the Case itself supplies its target.
[[nodiscard]] core::Result<CaseStaticPgPointRow> assemble_case_static_pg_point_row(
    const io::SiCase&, const io::schema::SupportConstraint&, const CaseStaticPgPointPatch&,
    const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

// Compare a homogeneous translational SPC selection to a static Case point support
// on an original node of its qualified planar Shell. The result owns the evidence.
[[nodiscard]] core::Result<CaseNastranPointSupport> certify_case_nastran_static_point_support(
    const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&, const CaseNastranPlanarPgSet&,
    io::detail::NastranSupportSetId, io::detail::NastranGridId,
    const CaseNastranPointSupportLimits&, core::ResourcePlan&,
    std::optional<CaseNastranPointSupportDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::coupling::detail
