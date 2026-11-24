#pragma once

#include <uvlas/structure/pg_volume_position.hpp>

namespace uvlas::structure {

/// Explicitly select a horizontal material side or one external lateral face.
using PgBoundaryPositionPoint = std::variant<PgReferenceHorizontalPoint, PgReferenceLateralPoint>;

enum class PgBoundaryPositionSide : std::uint8_t { Lower, Upper, Lateral };

struct PgBoundaryPositionProjection {
    core::Precision representation;
    double maximum_shape_error; ///< Dimensionless N, applied to every coefficient.
};

struct PgBoundaryPositionRecord {
    PgReferenceIdentity identity;
    PgCellIndex cell;
    core::Index partition_record{};
    core::LayerId layer{0};
    PgBoundaryPositionSide side{};
};

/// One row per face point, with N in the certified column order. Low is empty
/// for Binary64 and has the same shape as high for TwoWord.
struct PgBoundaryPositionWords {
    std::span<const double> high;
    std::span<const double> low;
};

/// Owns the certified basis and its projected boundary position values.
class PgBoundaryPositionReference {
public:
    PgBoundaryPositionReference(PgBoundaryPositionReference&&) noexcept = default;
    PgBoundaryPositionReference& operator=(PgBoundaryPositionReference&&) noexcept = default;
    PgBoundaryPositionReference(const PgBoundaryPositionReference&) = delete;
    PgBoundaryPositionReference& operator=(const PgBoundaryPositionReference&) = delete;

    /// Transfer the certificate on success. Preserve the input point order,
    /// including coincident points with distinct side or layer identities.
    [[nodiscard]] static core::Result<PgBoundaryPositionReference>
    prepare(PgConditioningCertificate&&, std::span<const PgBoundaryPositionPoint>,
            const PgBoundaryPositionProjection&, const PgReferenceLimits&,
            core::ResourcePlan&) noexcept;

    [[nodiscard]] std::size_t column_count() const noexcept {
        return column_count_;
    }
    [[nodiscard]] PgBoundaryPositionProjection projection() const noexcept {
        return projection_;
    }
    [[nodiscard]] PgReferenceCallInfo call_info() const noexcept {
        return call_info_;
    }

    [[nodiscard]] std::span<const PgBoundaryPositionRecord> points() const& noexcept {
        return points_.values();
    }
    [[nodiscard]] std::span<const PgBoundaryPositionRecord> points() const&& = delete;
    [[nodiscard]] PgBoundaryPositionWords words() const& noexcept {
        return {high_.values(), low_.values()};
    }
    [[nodiscard]] PgBoundaryPositionWords words() const&& = delete;
    [[nodiscard]] const PgConditioningCertificate& conditioning() const& noexcept {
        return conditioning_;
    }
    [[nodiscard]] const PgConditioningCertificate& conditioning() const&& = delete;

private:
    PgBoundaryPositionReference(PgConditioningCertificate&&, std::size_t,
                                PgBoundaryPositionProjection, PgReferenceCallInfo,
                                core::PlannedBuffer<PgBoundaryPositionRecord>,
                                core::PlannedBuffer<double>, core::PlannedBuffer<double>) noexcept;

    PgConditioningCertificate conditioning_;
    std::size_t column_count_;
    PgBoundaryPositionProjection projection_;
    PgReferenceCallInfo call_info_;
    core::PlannedBuffer<PgBoundaryPositionRecord> points_;
    core::PlannedBuffer<double> high_, low_;
};

/// Apply the ordered boundary N to Cartesian coordinates or N^T to Cartesian
/// point forces. The same reusable destination and Workspace serve either action.
[[nodiscard]] core::Result<void>
apply_pg_boundary_position(const PgBoundaryPositionReference&, PgVolumePositionDirection,
                           PgCartesianWords, const PgVolumePositionLimits&, core::Workspace&,
                           PgVolumePositionResult&) noexcept;

} // namespace uvlas::structure
