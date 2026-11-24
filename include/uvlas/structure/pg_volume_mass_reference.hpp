#pragma once

#include <uvlas/core/error_certificate.hpp>
#include <uvlas/structure/pg_conditioning.hpp>

namespace uvlas::structure {

/// Absolute errors of the actual stored representation, in the native units.
struct PgVolumeMassProjection {
    core::Precision representation;
    double maximum_shape_error;                        ///< Dimensionless N.
    double maximum_volume_jacobian_error_metres_cubed; ///< Per mapped reference point.
};

struct PgVolumeMassPoint {
    PgReferenceIdentity identity;
    PgCellIndex cell;
    core::Index partition_record{};
    core::LayerId layer{0};
};

/// Point-major, with each row containing one physical volume Jacobian followed by
/// every conditioned N in the certificate's column order. The low span is empty
/// for Binary64 and the same length as high for TwoWord.
struct PgVolumeMassWords {
    std::span<const double> high;
    std::span<const double> low;
};

/// Fixed volume data for subsequent mass integration at caller-selected points.
/// The basis and the admitted values move together, preventing reinterpretation
/// with a separately prepared transform.
class PgVolumeMassReference {
public:
    PgVolumeMassReference(PgVolumeMassReference&&) noexcept = default;
    PgVolumeMassReference& operator=(PgVolumeMassReference&&) noexcept = default;
    PgVolumeMassReference(const PgVolumeMassReference&) = delete;
    PgVolumeMassReference& operator=(const PgVolumeMassReference&) = delete;

    /// On failure the certificate remains valid. No quadrature weights or
    /// densities are introduced by this geometric preparation.
    [[nodiscard]] static core::Result<PgVolumeMassReference>
    prepare(PgConditioningCertificate&&, std::span<const PgReferenceVolumePoint>,
            const PgVolumeMassProjection&, const PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] std::size_t column_count() const noexcept {
        return column_count_;
    }
    [[nodiscard]] PgVolumeMassProjection projection() const noexcept {
        return projection_;
    }
    [[nodiscard]] PgReferenceCallInfo call_info() const noexcept {
        return call_info_;
    }

    [[nodiscard]] std::span<const PgVolumeMassPoint> points() const& noexcept {
        return points_.values();
    }
    [[nodiscard]] std::span<const PgVolumeMassPoint> points() const&& = delete;
    [[nodiscard]] PgVolumeMassWords words() const& noexcept {
        return {high_.values(), low_.values()};
    }
    [[nodiscard]] PgVolumeMassWords words() const&& = delete;
    [[nodiscard]] const PgConditioningCertificate& conditioning() const& noexcept {
        return conditioning_;
    }
    [[nodiscard]] const PgConditioningCertificate& conditioning() const&& = delete;

private:
    PgVolumeMassReference(PgConditioningCertificate&&, std::size_t, PgVolumeMassProjection,
                          PgReferenceCallInfo, core::PlannedBuffer<PgVolumeMassPoint>,
                          core::PlannedBuffer<double>, core::PlannedBuffer<double>) noexcept;

    PgConditioningCertificate conditioning_;
    std::size_t column_count_;
    PgVolumeMassProjection projection_;
    PgReferenceCallInfo call_info_;
    core::PlannedBuffer<PgVolumeMassPoint> points_;
    core::PlannedBuffer<double> high_, low_;
};

} // namespace uvlas::structure
