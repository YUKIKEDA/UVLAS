#include <uvlas/numerics/eigen_common.hpp>

namespace uvlas::numerics::eigen {
EigenFailure complete_failure(core::Failure cause, const Context& context,
                              StateDiagnosticContext state,
                              AttemptDiagnostics diagnostics) noexcept {
    // A nested action can identify a different implicated component/matrix. Preserve
    // that information; an outer context supplies only fields the action did not know.
    if (!cause.details.model)
        cause.details.model = context.cache.model;
    if (!cause.details.trial)
        cause.details.trial = context.trial;
    if (!cause.details.numerical_generation)
        cause.details.numerical_generation = context.cache.matrix;
    if (!cause.details.parent)
        cause.details.parent = state.parent;
    if (!cause.details.last_accepted)
        cause.details.last_accepted = state.last_accepted;
    return {cause, diagnostics};
}

EigenFailure complete_failure(EigenFailure failure, const Context& context,
                              StateDiagnosticContext state,
                              const AttemptDiagnostics& fallback) noexcept {
    auto& report = failure.diagnostics;
    if (!report.requested_count)
        report.requested_count = fallback.requested_count;
    if (!report.converged_count)
        report.converged_count = fallback.converged_count;
    if (!report.backend_iterations)
        report.backend_iterations = fallback.backend_iterations;
    if (!report.extraction_attempts)
        report.extraction_attempts = fallback.extraction_attempts;
    if (!report.forward_applications)
        report.forward_applications = fallback.forward_applications;
    if (!report.adjoint_applications)
        report.adjoint_applications = fallback.adjoint_applications;
    if (!report.inverse_applications)
        report.inverse_applications = fallback.inverse_applications;
    return complete_failure(failure.cause, context, state, report);
}
} // namespace uvlas::numerics::eigen
