#pragma once

#include <cstdint>
#include <uvlas/core/failure.hpp>
#include <variant>

namespace uvlas::core {

/// Representation of the evaluated field; a tag alone does not qualify its consumer (N4-06).
enum class Precision : std::uint8_t {
    Binary64, ///< Field evaluated in binary64.
    TwoWord   ///< Field evaluated in the two-word representation.
};

/// Identity of the observation, scope and numerical state being certified (N4-02).
/// The evaluator is responsible for their correspondence to the supplied error bound.
struct CertificateContext {
    /// Observable defined by the observation plan.
    ObservationId observation;
    /// Evaluation scope defined by that plan; a narrower scope cannot certify a wider one.
    ScopeId scope;
    /// Numerical generation actually evaluated.
    NumericalGeneration numerical_generation;
    /// Representation of the evaluated field, not of the stored double budget.
    Precision representation;
    /// Store explicit evaluation identities; membership and consistency are not checked here.
    constexpr CertificateContext(ObservationId observation_id, ScopeId scope_id,
                                 NumericalGeneration generation, Precision precision) noexcept
        : observation(observation_id),
          scope(scope_id),
          numerical_generation(generation),
          representation(precision) {}
    /// Compare every identity and the representation tag exactly.
    bool operator==(const CertificateContext&) const = default;
};

/// No bound has been evaluated; this state cannot be interpreted as achievement.
struct Unevaluated {
    /// All unevaluated state markers compare equal; context belongs to the certificate.
    bool operator==(const Unevaluated&) const = default;
};
/// The evaluator's nonnegative finite upper bound is at most the certificate budget.
struct Achieved {
    /// Reported conservative bound, in the same scale as the certificate budget.
    double upper_bound;
};
/// The evaluator's nonnegative finite upper bound exceeds the certificate budget.
struct Exceeded {
    /// Reported conservative bound, in the same scale as the certificate budget.
    double upper_bound;
};
/// Explicit assessment state; inspect the alternative instead of testing a success flag.
using ErrorAssessment = std::variant<Unevaluated, Achieved, Exceeded>;

/// Owned, validated record of a supplied error bound against its budget (N4-02).
/// Does not compute or establish conservatism of a bound, verify scope membership,
/// convert a solution, or qualify TwoWord consumers. Those are evaluator responsibilities.
/// Construction is restricted to factories so unevaluated cannot masquerade as achieved.
class ErrorCertificate {
public:
    /// Create a record with no evaluated bound; no achievement is implied.
    /// @param context Observation, scope, numerical generation and field representation.
    /// @param budget Nonnegative finite budget in the observation's prescribed scale.
    /// Zero is allowed and requires exact agreement when later assessed.
    /// @return Owned Unevaluated record, or Nonfinite, InvalidInput,
    /// then PrecisionUnsupported in that validation order. Failures identify the
    /// observation and numerical generation with Operation::CertifyError.
    [[nodiscard]] static Result<ErrorCertificate> unevaluated(CertificateContext context,
                                                              double budget) noexcept;
    /// Classify a supplied conservative bound without an added numerical tolerance.
    /// @param context Identities and field representation actually evaluated by the caller.
    /// @param budget Nonnegative finite budget for this observation and scope.
    /// @param upper_bound Nonnegative finite conservative bound in the same scale as budget.
    /// @return Owned Achieved record when upper_bound <= budget, otherwise an Exceeded
    /// record. Exceeded is a valid record, not a failed factory or a successful analysis.
    /// Invalid inputs return Nonfinite, InvalidInput, then PrecisionUnsupported in that
    /// order, with observation and numerical generation under Operation::CertifyError.
    [[nodiscard]] static Result<ErrorCertificate> assess(CertificateContext context, double budget,
                                                         double upper_bound) noexcept;

    /// Borrow immutable context; the reference is valid only while this object lives
    /// and is not assigned to or moved from.
    [[nodiscard]] const CertificateContext& context() const noexcept {
        return context_;
    }
    /// Return the validated budget by value, without altering its binary64 bits.
    [[nodiscard]] double budget() const noexcept {
        return budget_;
    }
    /// Borrow the explicit state; no bound exists in the Unevaluated alternative.
    /// The reference is valid only while this object lives and is not assigned to or moved from.
    [[nodiscard]] const ErrorAssessment& assessment() const noexcept {
        return assessment_;
    }

private:
    ErrorCertificate(CertificateContext context, double budget, ErrorAssessment assessment) noexcept
        : context_(context), budget_(budget), assessment_(assessment) {}

    CertificateContext context_;
    double budget_;
    ErrorAssessment assessment_;
};

} // namespace uvlas::core
