#pragma once

#include <atomic>
#include <uvlas/numerics/eigen_arnoldi.hpp>

namespace uvlas::numerics::eigen::detail {
using ArnoldiComplex = std::complex<double>;
template <class T>
using ArnoldiResult = ::uvlas::numerics::eigen::Result<T>;

// Capacity is explicit and separate from count. Unused storage never enters a
// view, a metric, a consumer, or a certificate. Both result generations own their
// arrays independently when projection is requested.
struct ArnoldiPacket {
    core::PlannedBuffer<ArnoldiComplex> eigenvalues, right, left, right_t, left_t;
    core::PlannedBuffer<std::optional<PairMetrics>> right_pairs, left_pairs;
    core::PlannedBuffer<core::Index> right_mapping, left_mapping;
    std::optional<BasisMetrics> right_basis, left_basis;
    core::Index count = 0;
};
struct ArnoldiWorkspaceData {
    core::Index order = 0, capacity = 0;
    std::atomic_flag busy = ATOMIC_FLAG_INIT;
    core::PlannedBuffer<double> real_input, real_input_low, real_high, real_low, real_error;
    core::PlannedBuffer<ArnoldiComplex> input, input_low, high, low, cluster_action;
    core::PlannedBuffer<ArnoldiComplex> candidate_values, candidate_vectors;
    core::PlannedBuffer<double> error, cluster_error, triangular_error;
    core::PlannedBuffer<core::Index> assigned;
    ArnoldiPacket left_candidates;
};
class ArnoldiWorkspaceAccess {
public:
    [[nodiscard]] static ArnoldiWorkspaceData* data(ArnoldiWorkspace& workspace) noexcept {
        return workspace.data_.get();
    }
};
struct ArnoldiData {
    ArnoldiData(Context original, ActionStamp action_stamp, StateDiagnosticContext reporting,
                linear::ReturnedPrecision precision, ArnoldiAssociation matching,
                ArnoldiSelectionMetadata transform, ArnoldiBackend versions,
                linear::UnitDimensions lambda_units, double lambda_scale, bool request_left,
                bool request_subspace, ArnoldiSettings configured_settings,
                ArnoldiLimits configured_limits) noexcept
        : context(original),
          original_stamp(action_stamp),
          state(reporting),
          returned(precision),
          association(matching),
          selection(transform),
          backend(versions),
          eigenvalue_units(lambda_units),
          eigenvalue_scale(lambda_scale),
          left_requested(request_left),
          subspace_requested(request_subspace),
          settings(configured_settings),
          limits(configured_limits),
          diagnostics({}, {}, configured_settings.nev, configured_settings.ncv, 0, 0, std::nullopt,
                      0, 0, {}) {}
    Context context;
    ActionStamp original_stamp;
    StateDiagnosticContext state;
    linear::ReturnedPrecision returned;
    ArnoldiAssociation association;
    ArnoldiSelectionMetadata selection;
    ArnoldiBackend backend;
    linear::UnitDimensions eigenvalue_units;
    double eigenvalue_scale;
    bool left_requested;
    bool subspace_requested;
    ArnoldiSettings settings;
    ArnoldiLimits limits;
    core::PlannedBuffer<std::optional<linear::ConsumerCheck>> required_checks;
    bool projected = false;
    core::Index order = 0;
    ArnoldiPacket retained, projection;
    core::PlannedBuffer<double> components;
    core::PlannedBuffer<linear::UnitDimensions> units;
    core::PlannedBuffer<ArnoldiComplex> initial;
    core::PlannedBuffer<ArnoldiSensitivity> sensitivities;
    core::PlannedBuffer<std::optional<ArnoldiInverseAssessment>> inverse_assessments;
    core::PlannedBuffer<std::optional<core::ErrorCertificate>> consumer_proofs;
    std::size_t inverse_count = 0, consumer_count = 0;
    ArnoldiDiagnostics diagnostics;
    ArnoldiMeasurements measurements{ArnoldiUnmeasured{}};
    std::optional<linear::BackendAdmission> admission;
};

[[nodiscard]] inline ArnoldiMatrixView
arnoldi_matrix(const core::PlannedBuffer<ArnoldiComplex>& owner, core::Index rows,
               core::Index columns) noexcept {
    return {
        {rows, columns, rows},
        owner.values().first(static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns))};
}
[[nodiscard]] inline ArnoldiCandidateView
arnoldi_view(const ArnoldiData& data, const ArnoldiPacket& packet,
             core::NumericalGeneration generation) noexcept {
    const auto right = arnoldi_matrix(packet.right, data.order, packet.count);
    const auto left = arnoldi_matrix(packet.left, data.left_requested ? data.order : 0,
                                     data.left_requested ? packet.count : 0);
    const ArnoldiVectorView right_view =
        data.subspace_requested
            ? ArnoldiVectorView{ArnoldiSubspaceView{
                  right, arnoldi_matrix(packet.right_t, packet.count, packet.count),
                  packet.right_mapping.values().first(static_cast<std::size_t>(packet.count))}}
            : ArnoldiVectorView{ArnoldiPairView{right}};
    ArnoldiLeftView left_view{NotRequested{}};
    if (data.left_requested) {
        if (data.subspace_requested) {
            left_view = ArnoldiSubspaceView{
                left, arnoldi_matrix(packet.left_t, packet.count, packet.count),
                packet.left_mapping.values().first(static_cast<std::size_t>(packet.count))};
        } else {
            left_view = ArnoldiPairView{left};
        }
    }
    return {arnoldi_matrix(packet.eigenvalues, packet.count, 1),
            right_view,
            left_view,
            data.context,
            data.state,
            generation,
            core::Precision::Binary64};
}
[[nodiscard]] inline core::NumericalGeneration
retained_generation(linear::ReturnedPrecision policy) noexcept {
    if (const auto* retained = std::get_if<linear::RetainedPrecision>(&policy)) {
        return retained->generation;
    }
    return std::get<linear::Binary64Projection>(policy).retained_generation;
}
} // namespace uvlas::numerics::eigen::detail
