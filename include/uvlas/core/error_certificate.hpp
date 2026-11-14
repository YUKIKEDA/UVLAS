#pragma once

#include <cstdint>
#include <uvlas/core/failure.hpp>
#include <variant>

namespace uvlas::core {

/// Arithmetic representation used to evaluate the field associated with an error bound.
/// Selecting a tag does not perform a conversion or change the arithmetic of a calculation.
enum class Precision : std::uint8_t {
    Binary64, ///< Field evaluated in binary64.
    TwoWord   ///< Field evaluated using values represented by two floating-point words.
};

/// Identifies the measured quantity, evaluation domain and state to which a bound applies.
/// Supply identities that refer to the calculation that produced the bound. These
/// values are stored without checking them against a model or a registry.
struct CertificateContext {
    /// Measured quantity, such as displacement or force, whose error was bounded.
    ObservationId observation;
    /// Domain or subset covered by the bound; a subset's bound does not cover a larger domain.
    ScopeId scope;
    /// Version of the numerical state evaluated; change it when that state changes.
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

/// Records whether a supplied error upper bound meets a permitted error budget.
/// Both values must describe the same measured quantity, domain, state and scale.
/// For an absolute error use the quantity's physical units for both; for a relative
/// error use the same normalization for both. This class does not choose the scale.
/// The caller must compute a conservative bound: this class checks its numeric
/// validity and compares it to the budget, without recomputing or verifying the bound.
/// Use unevaluated() when no bound is available, and inspect assessment() explicitly.
/// A successful factory call means the record is valid, not that its budget was met.
/// @code{.cpp}
/// using namespace uvlas::core;
/// const CertificateContext context{ObservationId{1}, ScopeId{2},
///                                  NumericalGeneration{3}, Precision::Binary64};
/// const auto record = ErrorCertificate::assess(context, 1e-6, 5e-7);
/// if (record) {
///     if (std::holds_alternative<Achieved>(record->assessment())) {
///         // The reported bound is no larger than the budget.
///     }
/// }
/// @endcode
class ErrorCertificate {
public:
    /// Create a record with no evaluated bound; no achievement is implied.
    /// @param context Observation, scope, numerical generation and field representation.
    /// @param budget Nonnegative finite permitted error, in the scale used for this quantity.
    /// Zero is allowed; only a zero upper bound can meet a zero budget.
    /// @return Owned Unevaluated record. Invalid inputs return Failure: Nonfinite for
    /// a nonfinite budget, InvalidInput for a negative budget, or PrecisionUnsupported
    /// for an unknown representation tag, checked in that order. The failure identifies
    /// the supplied observation and numerical generation with Operation::CertifyError.
    [[nodiscard]] static Result<ErrorCertificate> unevaluated(CertificateContext context,
                                                              double budget) noexcept;
    /// Classify a supplied conservative bound without an added numerical tolerance.
    /// @param context Identities and field representation actually evaluated by the caller.
    /// @param budget Nonnegative finite budget for this observation and scope.
    /// @param upper_bound Nonnegative finite conservative bound in the same scale as budget.
    /// @return Owned Achieved record when upper_bound <= budget, otherwise an Exceeded
    /// record. Exceeded is a valid record, not a failed factory or a successful analysis.
    /// Invalid inputs return Failure: Nonfinite if either number is nonfinite,
    /// InvalidInput if either is negative, or PrecisionUnsupported if the representation
    /// tag is unknown, checked in that order. The failure identifies the supplied
    /// observation and numerical generation with Operation::CertifyError.
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
