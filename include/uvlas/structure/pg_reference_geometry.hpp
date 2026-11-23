#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stop_token>
#include <uvlas/structure/pg_coefficients.hpp>
#include <uvlas/structure/pg_integration_partition.hpp>
#include <variant>

namespace uvlas::structure {

/// Three product candidates, consumed only when reference preparation succeeds.
struct PgReferenceInputs {
    PgGraph graph;
    PgCoefficientDictionary dictionary;
    PgIntegrationPartition partition;
};

/// Read-only view of the prepared product owners. Valid while they remain owned.
struct PgReferenceSource {
    const PgGraph& graph;
    const PgCoefficientDictionary& dictionary;
    const PgIntegrationPartition& partition;
};

/// Exact physical frame. Rotation entries are in row-major Cartesian order.
struct PgReferenceFrame {
    const std::array<ExactRational, 9>& rotation;
    const std::array<ExactRational, 3>& origin_metres;
    const ExactRational& length_metres;
};

/// Bounds on one preparation or point query. Zero does not mean unlimited.
struct PgReferenceLimits {
    ExactLimits algebra;
    std::uint64_t maximum_scalar_operations;
    std::uint64_t maximum_partition_visits; ///< Includes graph comparison and its limb work.
    std::stop_token cancellation;
};

/// Reference triangle coordinates satisfy u,v >= 0 and u+v <= 1.
struct PgReferenceVolumePoint {
    core::Index volume;
    const ExactRational& u;
    const ExactRational& v;
    const ExactRational& depth; ///< Fraction in the selected layer, in [0,1].
};

/// Selects the material side via a horizontal face record, including interfaces.
struct PgReferenceHorizontalPoint {
    core::Index face;
    const ExactRational& u;
    const ExactRational& v;
};

/// An external boundary segment in the domain's counterclockwise orientation.
struct PgReferenceLateralPoint {
    core::Index face;
    const ExactRational& along; ///< Fraction along the directed boundary, in [0,1].
    const ExactRational& depth; ///< Fraction in the selected layer, in [0,1].
};

using PgReferencePoint =
    std::variant<PgReferenceVolumePoint, PgReferenceHorizontalPoint, PgReferenceLateralPoint>;

/// Value a+b sqrt(g), where g belongs to the same owned evaluation result.
struct PgReferenceScalar {
    ExactRational rational;
    ExactRational radical;
};

using PgReferenceVector = std::array<PgReferenceScalar, 3>;

/// At the reference triangle and layer coordinates of a volume query (m³).
struct PgReferenceVolumeMeasure {
    PgReferenceScalar jacobian;
};

/// Outward vector (m²) and squared norm (m⁴). Area is its positive square root.
struct PgReferenceAreaMeasure {
    PgReferenceVector outward;
    PgReferenceScalar squared_norm;
};

/// A horizontal face keeps its material side, even at a shared physical interface.
struct PgReferenceHorizontalMeasure {
    PgPartitionFaceSide side;
    PgReferenceAreaMeasure area;
};

/// A lateral face is an outer boundary, with the domain's counterclockwise tangent.
struct PgReferenceLateralMeasure {
    PgReferenceAreaMeasure area;
};

/// Raw dictionary information before activity rank or conditioning.
struct PgReferenceColumnValue {
    PgCoefficientColumn column;
    /// R/P/C rows, derivative columns 0,s,t,ss,st,tt in chart coordinates.
    std::array<std::array<ExactRational, 6>, 3> numerators;
    PgReferenceScalar value;
    PgReferenceVector chart_derivatives; ///< Partial derivatives in s,t,zeta order.
};

/// Content stamp of the exact physical reference; it cannot confer qualification.
struct PgReferenceIdentity {
    std::array<std::byte, 32> bytes{};
    bool operator==(const PgReferenceIdentity&) const = default;
};

/// Successful call work and cumulative owned storage, separate from limb estimates.
struct PgReferenceCallInfo {
    std::uint64_t scalar_operations = 0;
    std::uint64_t partition_visits = 0;
    std::uint64_t allocated_owned_payload_bytes = 0;
};

/// Fully owned point result with one exact radical and physical Cartesian data.
struct PgReferenceValue {
    PgReferenceIdentity identity;
    PgCellIndex cell;
    core::Index piece;
    core::Index thickness_triangle;
    core::Index partition_record;
    core::LayerId layer{0};
    ExactRational s, t, zeta;
    ExactRational full_thickness_metres;
    std::array<ExactRational, 2> thickness_gradient_metres;
    /// Graph height jets by total order 0..4, decreasing s-order within each order.
    std::array<ExactRational, 15> height_chart_jet;
    ExactRational g;
    PgReferenceVector position_metres;
    PgReferenceVector normal;
    std::array<PgReferenceVector, 2> normal_chart_derivatives;
    std::array<PgReferenceVector, 3> jacobian_columns_metres; ///< s,t,zeta columns.
    PgReferenceScalar determinant_metres_cubed;
    std::variant<PgReferenceVolumeMeasure, PgReferenceHorizontalMeasure, PgReferenceLateralMeasure>
        measure;
    std::optional<PgReferenceColumnValue> column;
    PgReferenceCallInfo call_info;
};

/// Owns the exact physical frame and a qualified bound on the local normal tube.
/// No prepared value claims global injectivity or an independently conditioned basis.
class PgReferenceGeometry {
public:
    PgReferenceGeometry(PgReferenceGeometry&&) noexcept = default;
    PgReferenceGeometry& operator=(PgReferenceGeometry&&) noexcept = default;
    PgReferenceGeometry(const PgReferenceGeometry&) = delete;
    PgReferenceGeometry& operator=(const PgReferenceGeometry&) = delete;

    [[nodiscard]] static core::Result<PgReferenceGeometry> prepare(PgReferenceInputs&&,
                                                                   const PgReferenceFrame&,
                                                                   const PgReferenceLimits&,
                                                                   core::ResourcePlan&) noexcept;

    /// Return owned physical values from the reference's actual partition.
    [[nodiscard]] core::Result<PgReferenceValue> evaluate(const PgReferencePoint&,
                                                          std::optional<PgCoefficientColumn>,
                                                          const PgReferenceLimits&,
                                                          core::ResourcePlan&) const noexcept;

    /// Check the Jacobian sign at one point without asserting whole-domain regularity.
    [[nodiscard]] static core::Result<void>
    diagnose_point(const PgReferenceInputs&, const PgReferenceFrame&, const PgReferencePoint&,
                   const PgReferenceLimits&, core::ResourcePlan&) noexcept;

    /// Canonical cross-call identity. A moved-from object cannot supply it.
    [[nodiscard]] core::Result<PgReferenceIdentity> identity() const noexcept;

    /// Counts for the preparation which produced this object.
    [[nodiscard]] PgReferenceCallInfo preparation() const noexcept;

    /// Borrow the actual validated three-owner source; temporary borrows are forbidden.
    [[nodiscard]] core::Result<PgReferenceSource> source() const& noexcept;
    [[nodiscard]] core::Result<PgReferenceSource> source() const&& = delete;

    /// Positive physical determinant bound before mapping a partition record.
    [[nodiscard]] core::Result<const ExactRational*> determinant_lower_bound() const& noexcept;
    [[nodiscard]] core::Result<const ExactRational*> determinant_lower_bound() const&& = delete;

private:
    friend class PgReferenceBasis;
    PgReferenceGeometry(PgReferenceInputs&& source, std::array<ExactRational, 9> rotation,
                        std::array<ExactRational, 3> origin, ExactRational length,
                        ExactRational determinant_bound, PgReferenceIdentity identity,
                        PgReferenceCallInfo info) noexcept;

    PgReferenceInputs inputs_;
    std::array<ExactRational, 9> rotation_;
    std::array<ExactRational, 3> origin_;
    ExactRational length_;
    ExactRational determinant_bound_;
    PgReferenceIdentity identity_;
    PgReferenceCallInfo info_;
};

} // namespace uvlas::structure
