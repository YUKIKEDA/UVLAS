#include <cmath>
#include <optional>
#include <uvlas/core/error_certificate.hpp>

namespace uvlas::core {
namespace {
[[nodiscard]] std::optional<Failure> validate(CertificateContext context, double budget,
                                              double upper_bound) noexcept {
    FailureDetails details;
    details.numerical_generation = context.numerical_generation;
    details.problem = context.observation;
    if (!std::isfinite(budget) || !std::isfinite(upper_bound)) {
        return Failure{FailureCode::Nonfinite, Operation::CertifyError, details};
    }
    if (budget < 0 || upper_bound < 0) {
        return Failure{FailureCode::InvalidInput, Operation::CertifyError, details};
    }
    if (context.representation != Precision::Binary64 &&
        context.representation != Precision::TwoWord) {
        return Failure{FailureCode::PrecisionUnsupported, Operation::CertifyError, details};
    }
    return std::nullopt;
}
} // namespace

Result<ErrorCertificate> ErrorCertificate::unevaluated(CertificateContext context,
                                                       double budget) noexcept {
    if (const auto failure = validate(context, budget, 0)) {
        return std::unexpected(*failure);
    }
    return ErrorCertificate{context, budget, Unevaluated{}};
}

Result<ErrorCertificate> ErrorCertificate::assess(CertificateContext context, double budget,
                                                  double upper_bound) noexcept {
    if (const auto failure = validate(context, budget, upper_bound)) {
        return std::unexpected(*failure);
    }
    // Bounds and budgets refer to the same observation/scope and scale. This
    // records an evaluator's bound; it does not generate or round a solution.
    if (upper_bound <= budget) {
        return ErrorCertificate{context, budget, Achieved{upper_bound}};
    }
    return ErrorCertificate{context, budget, Exceeded{upper_bound}};
}
} // namespace uvlas::core
