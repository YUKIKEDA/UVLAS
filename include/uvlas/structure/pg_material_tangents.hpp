#pragma once

#include <uvlas/structure/pg_prepared_patch.hpp>

namespace uvlas::structure {

/// Exact partials of the fixed material correspondence, in (s,t,zeta) order.
/// The two rows are derivatives with respect to the caller's surface parameters.
using PgMaterialDirections = std::array<std::array<ExactRational, 3>, 2>;

struct PgMaterialTangentProjection {
    core::Precision representation{};
    double maximum_derivative_error{};
};

/// Own two fixed tangent rows per qualified material point, in the prepared
/// patch's Cartesian column order. Valid at a selected interior or face side.
class PgMaterialTangentReference {
public:
    PgMaterialTangentReference(PgMaterialTangentReference&&) noexcept;
    PgMaterialTangentReference& operator=(PgMaterialTangentReference&&) noexcept;
    PgMaterialTangentReference(const PgMaterialTangentReference&) = delete;
    PgMaterialTangentReference& operator=(const PgMaterialTangentReference&) = delete;

    [[nodiscard]] static core::Result<PgMaterialTangentReference>
    prepare(const PgPreparedPatch&, std::span<const PgReferencePoint>,
            std::span<const PgMaterialDirections>, const PgMaterialTangentProjection&,
            const PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] const PgReferenceCallInfo& call_info() const noexcept {
        return call_info_;
    }

private:
    friend core::Result<void>
    apply_pg_material_tangents(const PgMaterialTangentReference&, PgPreparedPatchIdentity,
                               PgVolumePositionDirection, PgCartesianWords,
                               const PgVolumePositionLimits&, core::Workspace&,
                               PgVolumePositionResult&) noexcept;

    PgMaterialTangentReference(PgPreparedPatchIdentity, std::size_t, PgMaterialTangentProjection,
                               PgReferenceCallInfo, core::PlannedBuffer<PgProjectedReferencePoint>,
                               core::PlannedBuffer<double>, core::PlannedBuffer<double>) noexcept;
    PgPreparedPatchIdentity source_;
    std::size_t column_count_;
    PgMaterialTangentProjection projection_;
    PgReferenceCallInfo call_info_;
    core::PlannedBuffer<PgProjectedReferencePoint> points_;
    core::PlannedBuffer<double> high_, low_;
};

/// Apply the selected material directions or their structural transpose.
/// Inputs and outputs keep both retained words when requested.
[[nodiscard]] core::Result<void>
apply_pg_material_tangents(const PgMaterialTangentReference&, PgPreparedPatchIdentity,
                           PgVolumePositionDirection, PgCartesianWords,
                           const PgVolumePositionLimits&, core::Workspace&,
                           PgVolumePositionResult&) noexcept;

} // namespace uvlas::structure
