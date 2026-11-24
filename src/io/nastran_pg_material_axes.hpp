#pragma once

#include "external_mesh_sidecar_pg.hpp"

namespace uvlas::io::detail {

enum class NastranPgAxisFailure : std::uint8_t {
    Source,
    Thickness,
    CaseAxis,
    SourceAxis,
    AngleRepresentation,
    Difference,
};

struct NastranPgAxisDiagnostic {
    NastranPgAxisFailure reason = NastranPgAxisFailure::Source;
    std::optional<std::uint64_t> element{};
    std::optional<NastranPropertyId> property{};
    std::optional<core::LayerId> layer{};
    std::optional<NastranBulkPosition> element_field{};
    std::optional<NastranBulkPosition> ply_field{};
};

struct NastranPgAxisLimits {
    NastranPgAxisLimits(const structure::PgReferenceLimits& exact, std::uint64_t matches) noexcept
        : arithmetic(exact), maximum_matches(matches) {}
    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_matches;
};

struct NastranPgLayerAxis {
    std::uint64_t element{};
    NastranPropertyId property{0};
    core::LayerId layer{0};
    ExternalPgWinding winding = ExternalPgWinding::Forward;
    std::optional<NastranBulkPosition> element_field{};
    std::optional<NastranBulkPosition> ply_field{};
    // Exact rays in basic Cartesian coordinates. The angle is measured from the
    // tangent about the Case normal; it is an exact number of degrees.
    std::array<structure::ExactRational, 3> tangent;
    std::array<structure::ExactRational, 3> normal;
    structure::ExactRational degrees;
};

struct NastranPgMaterialAxes {
    RootSnapshotLocation deck;
    SourceDigest deck_digest;
    RootSnapshotLocation sidecar;
    SourceDigest sidecar_digest;
    structure::PgReferenceIdentity pg;
    MeshShellId shell{0};
    core::PlannedBuffer<NastranPgLayerAxis> matches; // Face order, then Case layer order.
    structure::PgReferenceCallInfo work;
};

// Compare the constant through-depth material triads of selected Nastran faces.
// The Case owner, already admitted geometry and explicit local thickness proof
// must all match. A projected Case axis may use equal general rotations when the
// residual relative rotation is an exact quarter turn. No tolerance is used.
[[nodiscard]] core::Result<NastranPgMaterialAxes>
match_nastran_pg_material_axes(const SiCase&, const ExternalMeshSidecarBinding&,
                               const QualifiedExternalSidecarPg&, const NastranPgAxisLimits&,
                               core::ResourcePlan&,
                               std::optional<NastranPgAxisDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
