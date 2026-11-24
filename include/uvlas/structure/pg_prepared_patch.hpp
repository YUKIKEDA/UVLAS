#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <uvlas/structure/pg_boundary_position.hpp>
#include <uvlas/structure/pg_discrete_gram.hpp>
#include <uvlas/structure/pg_volume_position.hpp>

namespace uvlas::structure {

/// Reference, exact selected basis, coordinate transform, and retained finite rule.
struct PgPreparedPatchIdentity {
    std::array<std::byte, 32> bytes{};
    bool operator==(const PgPreparedPatchIdentity&) const = default;
};

struct PgPreparedPatchInfo {
    PgPreparedPatchIdentity identity;
    PgReferenceIdentity geometry;
    PgConditioningInfo basis;
    PgReferenceCallInfo preparation;
};

using PgCartesianExactPosition = std::array<ExactRational, 3>;

enum class PgProjectedReferenceDomain : std::uint8_t { Volume, Lower, Upper, Lateral };

struct PgProjectedReferencePoint {
    PgCellIndex cell;
    core::Index partition_record{};
    core::LayerId layer{0};
    PgProjectedReferenceDomain domain = PgProjectedReferenceDomain::Volume;
};

/// Componentwise absolute errors in the recorded physical units.
struct PgPreparedReferenceProjection {
    core::Precision representation;
    double maximum_position_error_metres;
    double maximum_shape_error;
    double maximum_gradient_error_per_metre;
    double maximum_volume_error_metres_cubed;
    double maximum_area_error_metres_squared;
};

/// Point-major words. Each row is X0[3], positive measure, then (N, G[3]) per column.
class PgProjectedReferenceBatch {
public:
    PgProjectedReferenceBatch(PgProjectedReferenceBatch&&) noexcept = default;
    PgProjectedReferenceBatch& operator=(PgProjectedReferenceBatch&&) noexcept = default;
    PgProjectedReferenceBatch(const PgProjectedReferenceBatch&) = delete;
    PgProjectedReferenceBatch& operator=(const PgProjectedReferenceBatch&) = delete;

    [[nodiscard]] PgPreparedPatchIdentity source() const noexcept {
        return source_;
    }
    [[nodiscard]] std::size_t column_count() const noexcept {
        return columns_;
    }
    [[nodiscard]] PgPreparedReferenceProjection projection() const noexcept {
        return projection_;
    }
    [[nodiscard]] PgReferenceCallInfo call_info() const noexcept {
        return call_info_;
    }

    [[nodiscard]] std::span<const PgProjectedReferencePoint> points() const& noexcept {
        return points_.values();
    }
    [[nodiscard]] std::span<const PgProjectedReferencePoint> points() const&& = delete;
    [[nodiscard]] std::span<const double> high_words() const& noexcept {
        return high_.values();
    }
    [[nodiscard]] std::span<const double> high_words() const&& = delete;
    /// Empty in Binary64. TwoWord keeps the same layout as high_words().
    [[nodiscard]] std::span<const double> low_words() const& noexcept {
        return low_.values();
    }
    [[nodiscard]] std::span<const double> low_words() const&& = delete;

private:
    friend class PgPreparedPatch;
    PgProjectedReferenceBatch(PgPreparedPatchIdentity, std::size_t, PgPreparedReferenceProjection,
                              PgReferenceCallInfo, core::PlannedBuffer<PgProjectedReferencePoint>,
                              core::PlannedBuffer<double>, core::PlannedBuffer<double>) noexcept;

    PgPreparedPatchIdentity source_;
    std::size_t columns_;
    PgPreparedReferenceProjection projection_;
    PgReferenceCallInfo call_info_;
    core::PlannedBuffer<PgProjectedReferencePoint> points_;
    core::PlannedBuffer<double> high_, low_;
};

/// Use the N entries of a prepared point batch, in the patch's Cartesian column order.
[[nodiscard]] core::Result<void>
apply_pg_prepared_reference_position(const PgProjectedReferenceBatch&, PgPreparedPatchIdentity,
                                     PgVolumePositionDirection, PgCartesianWords,
                                     const PgVolumePositionLimits&, core::Workspace&,
                                     PgVolumePositionResult&) noexcept;

/// Immutable reference patch. q0 is Cartesian metres, by column and axis.
class PgPreparedPatch {
public:
    PgPreparedPatch(PgPreparedPatch&&) noexcept = default;
    PgPreparedPatch& operator=(PgPreparedPatch&&) noexcept = default;
    PgPreparedPatch(const PgPreparedPatch&) = delete;
    PgPreparedPatch& operator=(const PgPreparedPatch&) = delete;

    /// Preserve the caller's accepted conditioning source on failure.
    [[nodiscard]] static core::Result<PgPreparedPatch> prepare(PgConditioningCertificate&&,
                                                               const PgReferenceLimits&,
                                                               core::Workspace&,
                                                               core::ResourcePlan&) noexcept;
    /// Preserve the exact finite rule alongside its certified basis.
    [[nodiscard]] static core::Result<PgPreparedPatch> prepare(PgFixedDiscreteGramCertificate&&,
                                                               const PgReferenceLimits&,
                                                               core::Workspace&,
                                                               core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgPreparedPatchInfo> info() const noexcept;

    [[nodiscard]] std::span<const PgCartesianExactPosition>
    reference_coordinates() const& noexcept {
        return reference_.values();
    }
    [[nodiscard]] std::span<const PgCartesianExactPosition>
    reference_coordinates() const&& = delete;

    /// Exact scalar coordinates of one, usable for translation along any Cartesian axis.
    [[nodiscard]] std::span<const ExactRational> translation() const& noexcept {
        return translation_.values();
    }
    [[nodiscard]] std::span<const ExactRational> translation() const&& = delete;

    [[nodiscard]] const PgConditioningCertificate& conditioning() const& noexcept;
    [[nodiscard]] const PgConditioningCertificate& conditioning() const&& = delete;

    [[nodiscard]] std::optional<std::span<const PgOwnedDiscreteVolumePoint>>
    fixed_rule() const& noexcept;
    [[nodiscard]] std::optional<std::span<const PgOwnedDiscreteVolumePoint>>
    fixed_rule() const&& = delete;

    [[nodiscard]] core::Result<PgConditionedReferenceBatch>
    evaluate_reference(std::span<const PgReferencePoint>, const PgReferenceLimits&,
                       core::ResourcePlan&) const noexcept;

    /// Project the same qualified point order without changing or consuming this patch.
    [[nodiscard]] core::Result<PgProjectedReferenceBatch>
    project_reference(std::span<const PgReferencePoint>, const PgPreparedReferenceProjection&,
                      const PgReferenceLimits&, core::ResourcePlan&) const noexcept;

    /// Project material-side queries without introducing a volume point.
    [[nodiscard]] core::Result<PgProjectedReferenceBatch>
    project_boundary_reference(std::span<const PgBoundaryPositionPoint>,
                               const PgPreparedReferenceProjection&, const PgReferenceLimits&,
                               core::ResourcePlan&) const noexcept;

private:
    template <class Query>
    [[nodiscard]] core::Result<PgProjectedReferenceBatch>
    project_points(std::span<const Query>, const PgPreparedReferenceProjection&,
                   const PgReferenceLimits&, core::ResourcePlan&) const noexcept;

    struct Coordinates;
    using Conditioning = std::variant<PgConditioningCertificate, PgFixedDiscreteGramCertificate>;
    PgPreparedPatch(Conditioning, Coordinates) noexcept;

    [[nodiscard]] static Coordinates derive(const PgConditioningCertificate&,
                                            std::span<const PgOwnedDiscreteVolumePoint>,
                                            const PgReferenceLimits&, core::Workspace&,
                                            core::ResourcePlan&);

    Conditioning conditioning_;
    core::PlannedBuffer<PgCartesianExactPosition> reference_;
    core::PlannedBuffer<ExactRational> translation_;
    PgPreparedPatchInfo info_;
};

} // namespace uvlas::structure
