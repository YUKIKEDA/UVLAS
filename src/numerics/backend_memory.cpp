#include <uvlas/numerics/backend_memory.hpp>

namespace uvlas::numerics::linear {
core::Result<BackendAdmission> admit_backend_memory(BackendMemoryPlan plan,
                                                    BackendMemoryLimit limit) noexcept {
    std::uint64_t selected = plan.peak_payload_estimate_bytes();
    std::uint64_t maximum = 0;
    const auto selection = limit.selection();
    if (const auto* forecast = std::get_if<EstimateAdmission>(&selection)) {
        maximum = forecast->bytes;
    } else {
        const auto bound = plan.proven_peak_payload_upper_bound_bytes();
        if (!bound) {
            return std::unexpected(
                core::Failure{core::FailureCode::Unsupported, core::Operation::SolveLinear});
        }
        selected = *bound;
        maximum = std::get_if<RequireProvenBound>(&selection)->bytes;
    }
    if (selected > maximum) {
        // Forecasts and bounds are not observed allocation counts. Do not fill resources.
        return std::unexpected(
            core::Failure{core::FailureCode::ResourceExceeded, core::Operation::SolveLinear});
    }
    return BackendAdmission{plan, limit, selected};
}
} // namespace uvlas::numerics::linear
