#pragma once

#include <uvlas/structure/pg_volume_mass_reference.hpp>

namespace uvlas::structure {

/// Select the fixed Cartesian position map or its work-conjugate transpose.
enum class PgVolumePositionDirection : std::uint8_t { ToPoints, TransposeToColumns };

/// Three Cartesian components per basis column or evaluation row, contiguous in
/// x,y,z order. Representation identifies whether the low span is required.
struct PgCartesianWords {
    core::Precision representation;
    std::span<const double> high;
    std::span<const double> low;
};

struct PgVolumePositionLimits {
    std::uint64_t maximum_steps;
    double maximum_absolute_error;
    std::stop_token cancellation;
};

class PgBoundaryPositionReference;
class PgBoundaryTangentReference;
class PgMaterialTangentReference;
class PgPreparedVolumeReference;
struct PgPreparedPatchIdentity;
class PgProjectedReferenceBatch;

/// Reusable owned destination. Its prior words remain valid after any failed
/// action. Stage writes in a separate exclusive Workspace before publishing.
class PgVolumePositionResult {
public:
    PgVolumePositionResult(PgVolumePositionResult&&) noexcept = default;
    PgVolumePositionResult& operator=(PgVolumePositionResult&&) noexcept = default;
    PgVolumePositionResult(const PgVolumePositionResult&) = delete;
    PgVolumePositionResult& operator=(const PgVolumePositionResult&) = delete;

    /// Capacity is a number of scalar components, including all three axes.
    [[nodiscard]] static core::Result<PgVolumePositionResult>
    allocate(core::Index capacity, core::Precision, core::ResourcePlan&) noexcept;

    /// Scratch to stage up to this destination's full capacity.
    [[nodiscard]] core::Index workspace_bytes() const noexcept;

    /// Empty before the first successful action. Borrow until the next
    /// successful action or destruction of this destination.
    [[nodiscard]] PgCartesianWords words() const& noexcept;
    [[nodiscard]] PgCartesianWords words() const&& = delete;

private:
    friend core::Result<void> apply_pg_volume_position(const PgVolumeMassReference&,
                                                       PgVolumePositionDirection, PgCartesianWords,
                                                       const PgVolumePositionLimits&,
                                                       core::Workspace&,
                                                       PgVolumePositionResult&) noexcept;
    friend core::Result<void>
    apply_pg_boundary_position(const PgBoundaryPositionReference&, PgVolumePositionDirection,
                               PgCartesianWords, const PgVolumePositionLimits&, core::Workspace&,
                               PgVolumePositionResult&) noexcept;

    friend core::Result<void>
    apply_pg_boundary_tangents(const PgBoundaryTangentReference&, PgVolumePositionDirection,
                               PgCartesianWords, const PgVolumePositionLimits&, core::Workspace&,
                               PgVolumePositionResult&) noexcept;

    friend core::Result<void>
    apply_pg_prepared_volume_position(const PgPreparedVolumeReference&, PgPreparedPatchIdentity,
                                      PgVolumePositionDirection, PgCartesianWords,
                                      const PgVolumePositionLimits&, core::Workspace&,
                                      PgVolumePositionResult&) noexcept;

    friend core::Result<void>
    apply_pg_prepared_reference_position(const PgProjectedReferenceBatch&, PgPreparedPatchIdentity,
                                         PgVolumePositionDirection, PgCartesianWords,
                                         const PgVolumePositionLimits&, core::Workspace&,
                                         PgVolumePositionResult&) noexcept;

    friend core::Result<void>
    apply_pg_material_tangents(const PgMaterialTangentReference&, PgPreparedPatchIdentity,
                               PgVolumePositionDirection, PgCartesianWords,
                               const PgVolumePositionLimits&, core::Workspace&,
                               PgVolumePositionResult&) noexcept;

    struct MatrixView {
        const PgConditioningCertificate* certificate;
        std::size_t rank, rows;
        std::span<const double> high, low;
        std::size_t prefix;
        double maximum_coefficient_error;
        std::size_t column_stride = 1;
    };

    [[nodiscard]] core::Result<void> apply(const MatrixView&, PgVolumePositionDirection,
                                           PgCartesianWords, const PgVolumePositionLimits&,
                                           core::Workspace&) noexcept;

    PgVolumePositionResult(core::Precision, core::Index, core::PlannedBuffer<double>,
                           core::PlannedBuffer<double>) noexcept;

    core::Precision representation_;
    core::Index workspace_bytes_;
    core::PlannedBuffer<double> high_, low_;
    std::size_t size_ = 0;
};

/// Apply N or its transpose in fixed Cartesian order, using both input and
/// reference words. Verify the returned values against the unrounded N using
/// the preparation's error bound. No weights, density or measures enter.
[[nodiscard]] core::Result<void>
apply_pg_volume_position(const PgVolumeMassReference&, PgVolumePositionDirection, PgCartesianWords,
                         const PgVolumePositionLimits&, core::Workspace&,
                         PgVolumePositionResult&) noexcept;

} // namespace uvlas::structure
