#pragma once
#include "backend_memory_detail.hpp"
#include "eigen_metric_bounds.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <atomic>
#include <uvlas/numerics/eigen_full.hpp>

namespace uvlas::numerics::eigen::detail {
using FullComplex = std::complex<double>;
using FullMatrix = Eigen::Matrix<FullComplex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
struct FullIdentity {};
struct FullData {
    FullData(Context c, StateDiagnosticContext s, const FullOutputRequest& request,
             linear::BackendAdmission admission, std::optional<BasisMetricRequest> accuracy)
        : context(c),
          state(s),
          basis_request(accuracy),
          returned(request.returned),
          limits(request.limits),
          backend(admission) {}
    Context context;
    StateDiagnosticContext state;
    std::optional<BasisMetricRequest> basis_request;
    linear::ReturnedPrecision returned;
    FullLimits limits;
    linear::BackendAdmission backend;
    FullKind kind = FullKind::Spectrum;
    core::Index n = 0, k = 0;
    double eigenvalue_scale = 1;
    linear::UnitDimensions eigenvalue_units{};
    core::NumericalGeneration generation{0}, retained_generation{0};
    core::PlannedBuffer<FullComplex> original, values, right, left, triangular, left_triangular;
    core::PlannedBuffer<double> scales, triangular_error, left_triangular_error;
    core::PlannedBuffer<linear::UnitDimensions> units, coefficient_units;
    core::PlannedBuffer<core::Index> mapping;
    core::PlannedBuffer<std::optional<FullPairProof>> pairs, retained_pairs;
    core::PlannedBuffer<std::optional<core::ErrorCertificate>> consumers, retained_consumers;
    std::optional<BasisMetrics> basis, left_basis, retained_basis, retained_left_basis;
    std::optional<FullAssociation> association;
    std::shared_ptr<const FullIdentity> identity;
    AttemptDiagnostics diagnostics;
    std::optional<FullTimings> timings;
    std::uint64_t clock_reads = 0;
    bool timing_overflow = false;
    std::optional<double> independence_upper;
    bool left_requested = false;
    bool real_input = false;
    mutable std::atomic<bool> busy{false};
};
class FullFactory {
public:
    template <StorageScalar Scalar>
    static eigen::Result<FullSpectrum>
    spectrum(const PrecisionBuffer<Scalar>& original, Context context, StateDiagnosticContext state,
             FullScales scales, const FullRequest& request, core::ResourcePlan& plan) noexcept;
    static eigen::Result<FullSelection> selection(const FullSpectrum&, std::span<const core::Index>,
                                                  core::ResourcePlan&) noexcept;
    static eigen::Result<FullPairs> pairs(const FullSpectrum&, const FullSelection&,
                                          const FullPairRequest&, core::ResourcePlan&) noexcept;
    static eigen::Result<FullCluster> cluster(const FullSpectrum&, const FullSelection&,
                                              const FullClusterRequest&,
                                              core::ResourcePlan&) noexcept;
};
} // namespace uvlas::numerics::eigen::detail
