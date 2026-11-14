#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <limits>
#include <type_traits>
#include <uvlas/core/error_certificate.hpp>

using namespace uvlas::core;

namespace {
[[nodiscard]] CertificateContext context(Precision representation) {
    return {ObservationId{3}, ScopeId{4}, NumericalGeneration{5}, representation};
}
} // namespace

TEST_CASE("unevaluated certificates cannot be interpreted as achieved", "[unit][N4-02]") {
    for (const auto precision : {Precision::Binary64, Precision::TwoWord}) {
        const auto result = ErrorCertificate::unevaluated(context(precision), 1e-10);
        REQUIRE(result.has_value());
        CHECK(std::holds_alternative<Unevaluated>(result->assessment()));
        CHECK(std::get_if<Achieved>(&result->assessment()) == nullptr);
        CHECK(std::get_if<Exceeded>(&result->assessment()) == nullptr);
        CHECK(result->context() == context(precision));
        CHECK(std::bit_cast<std::uint64_t>(result->budget()) ==
              std::bit_cast<std::uint64_t>(1e-10));
    }
    STATIC_REQUIRE_FALSE(std::is_default_constructible_v<ErrorCertificate>);
    STATIC_REQUIRE_FALSE(std::is_convertible_v<ErrorCertificate, bool>);
    STATIC_REQUIRE_FALSE(
        std::is_constructible_v<CertificateContext, ObservationId, ScopeId, NumericalGeneration>);
    STATIC_REQUIRE_FALSE(
        std::is_constructible_v<ErrorCertificate, CertificateContext, double, ErrorAssessment>);
}

TEST_CASE("budget comparison accepts equality and rejects the next representable bound",
          "[unit][N4-02]") {
    for (const auto precision : {Precision::Binary64, Precision::TwoWord}) {
        // N4-02: compare the supplied bound and budget without a hidden tolerance.
        // Both representations label the evaluated field, not a solution conversion.
        for (const double budget : {0.0, std::numeric_limits<double>::denorm_min(), 1e-10, 1.0}) {
            const auto exact = ErrorCertificate::assess(context(precision), budget, budget);
            REQUIRE(exact.has_value());
            CHECK(std::holds_alternative<Achieved>(exact->assessment()));
            const double above = std::nextafter(budget, std::numeric_limits<double>::infinity());
            const auto missed = ErrorCertificate::assess(context(precision), budget, above);
            REQUIRE(missed.has_value());
            CHECK(std::holds_alternative<Exceeded>(missed->assessment()));
            const auto* assessment = std::get_if<Exceeded>(&missed->assessment());
            REQUIRE(assessment != nullptr);
            CHECK(std::bit_cast<std::uint64_t>(assessment->upper_bound) ==
                  std::bit_cast<std::uint64_t>(above));
        }
    }
}

TEST_CASE("certificates reject negative and nonfinite bounds or budgets",
          "[negative][N4-02][N4-08]") {
    constexpr double inf = std::numeric_limits<double>::infinity();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto precision : {Precision::Binary64, Precision::TwoWord}) {
        for (const double value :
             {-1.0, -std::numeric_limits<double>::denorm_min(), inf, -inf, nan}) {
            const auto expected =
                std::isfinite(value) ? FailureCode::InvalidInput : FailureCode::Nonfinite;
            const auto bad_budget = ErrorCertificate::unevaluated(context(precision), value);
            REQUIRE_FALSE(bad_budget.has_value());
            CHECK(bad_budget.error().code == expected);
            const auto assessed_bad_budget = ErrorCertificate::assess(context(precision), value, 0);
            REQUIRE_FALSE(assessed_bad_budget.has_value());
            CHECK(assessed_bad_budget.error().code == expected);
            const auto bad_bound = ErrorCertificate::assess(context(precision), 1e-10, value);
            REQUIRE_FALSE(bad_bound.has_value());
            CHECK(bad_bound.error().code == expected);
            CHECK(bad_bound.error().operation == Operation::CertifyError);
            CHECK(bad_bound.error().details.numerical_generation == NumericalGeneration{5});
            CHECK(bad_bound.error().details.problem == std::optional<ProblemId>{ObservationId{3}});
        }
    }
}

TEST_CASE("certificates reject unsupported representation tags", "[negative][N4-02][N4-08]") {
    const auto result = ErrorCertificate::assess(context(static_cast<Precision>(255)), 1e-10, 0);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == FailureCode::PrecisionUnsupported);
    const auto unassessed =
        ErrorCertificate::unevaluated(context(static_cast<Precision>(255)), 1e-10);
    REQUIRE_FALSE(unassessed.has_value());
    CHECK(unassessed.error().code == FailureCode::PrecisionUnsupported);
    const auto nonfinite = ErrorCertificate::assess(context(static_cast<Precision>(255)), -1,
                                                    std::numeric_limits<double>::quiet_NaN());
    REQUIRE_FALSE(nonfinite.has_value());
    CHECK(nonfinite.error().code == FailureCode::Nonfinite);
}
