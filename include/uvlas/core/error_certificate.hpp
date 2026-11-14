#pragma once

#include <cstdint>
#include <uvlas/core/failure.hpp>
#include <variant>

namespace uvlas::core {

enum class Precision : std::uint8_t { Binary64, TwoWord };

struct CertificateContext {
    ObservationId observation;
    ScopeId scope;
    NumericalGeneration numerical_generation;
    Precision representation;
    constexpr CertificateContext(ObservationId observation_id, ScopeId scope_id,
                                 NumericalGeneration generation, Precision precision) noexcept
        : observation(observation_id),
          scope(scope_id),
          numerical_generation(generation),
          representation(precision) {}
    bool operator==(const CertificateContext&) const = default;
};

struct Unevaluated {
    bool operator==(const Unevaluated&) const = default;
};
struct Achieved {
    double upper_bound;
};
struct Exceeded {
    double upper_bound;
};
using ErrorAssessment = std::variant<Unevaluated, Achieved, Exceeded>;

class ErrorCertificate {
public:
    [[nodiscard]] static Result<ErrorCertificate> unevaluated(CertificateContext context,
                                                              double budget) noexcept;
    [[nodiscard]] static Result<ErrorCertificate> assess(CertificateContext context, double budget,
                                                         double upper_bound) noexcept;

    [[nodiscard]] const CertificateContext& context() const noexcept {
        return context_;
    }
    [[nodiscard]] double budget() const noexcept {
        return budget_;
    }
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
