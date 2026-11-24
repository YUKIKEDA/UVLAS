#pragma once

#include <uvlas/structure/pg_boundary_position.hpp>

namespace uvlas::structure {

struct PgBoundaryTangentProjection {
    core::Precision representation;
    double maximum_derivative_error; ///< Dimensionless, common to both face directions.
};

struct PgBoundaryTangentWords {
    std::span<const double> high;
    std::span<const double> low;
};

/// Owns the projected derivatives and ordered face records. No source borrow is retained.
class PgBoundaryTangentReference {
public:
    PgBoundaryTangentReference(PgBoundaryTangentReference&&) noexcept;
    PgBoundaryTangentReference& operator=(PgBoundaryTangentReference&&) noexcept;
    PgBoundaryTangentReference(const PgBoundaryTangentReference&) = delete;
    PgBoundaryTangentReference& operator=(const PgBoundaryTangentReference&) = delete;

    /// Use (u,v) on a horizontal face and (along,depth) on a lateral face.
    [[nodiscard]] static core::Result<PgBoundaryTangentReference>
    prepare(const PgConditioningCertificate&, std::span<const PgBoundaryPositionPoint>,
            const PgBoundaryTangentProjection&, const PgReferenceLimits&,
            core::ResourcePlan&) noexcept;

    [[nodiscard]] std::size_t column_count() const noexcept {
        return column_count_;
    }
    [[nodiscard]] PgBoundaryTangentProjection projection() const noexcept {
        return projection_;
    }
    [[nodiscard]] PgReferenceCallInfo call_info() const noexcept {
        return call_info_;
    }

    [[nodiscard]] std::span<const PgBoundaryPositionRecord> points() const& noexcept {
        return points_.values();
    }
    [[nodiscard]] std::span<const PgBoundaryPositionRecord> points() const&& = delete;

    /// The first and second direction rows for each point, in that order.
    /// Each row contains all columns in the certified order.
    [[nodiscard]] PgBoundaryTangentWords words() const& noexcept {
        return {high_.values(), low_.values()};
    }
    [[nodiscard]] PgBoundaryTangentWords words() const&& = delete;

private:
    PgBoundaryTangentReference(std::size_t, PgBoundaryTangentProjection, PgReferenceCallInfo,
                               core::PlannedBuffer<PgBoundaryPositionRecord>,
                               core::PlannedBuffer<double>, core::PlannedBuffer<double>) noexcept;

    std::size_t column_count_;
    PgBoundaryTangentProjection projection_;
    PgReferenceCallInfo call_info_;
    core::PlannedBuffer<PgBoundaryPositionRecord> points_;
    core::PlannedBuffer<double> high_, low_;
};

/// ToPoints produces six components per point: xyz for each direction in chart order.
/// TransposeToColumns consumes the same six-component layout as conjugate covectors.
[[nodiscard]] core::Result<void>
apply_pg_boundary_tangents(const PgBoundaryTangentReference&, PgVolumePositionDirection,
                           PgCartesianWords, const PgVolumePositionLimits&, core::Workspace&,
                           PgVolumePositionResult&) noexcept;

} // namespace uvlas::structure
