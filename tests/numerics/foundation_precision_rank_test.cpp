#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <utility>
#include <uvlas/numerics/dense_solve.hpp>
#include <uvlas/numerics/precision.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
using Rational = boost::multiprecision::cpp_rational;
constexpr auto operation = core::Operation::SolveLinear;
constexpr auto source_generation = core::NumericalGeneration{17};
constexpr double tiny = 0x1p-180;

template <std::size_t N>
[[nodiscard]] DenseReal matrix(DenseShape shape, const std::array<double, N>& values,
                               core::ResourcePlan& plan) {
    auto result = DenseReal::copy(shape, values, plan, operation);
    REQUIRE(result.has_value());
    return std::move(*result);
}

[[nodiscard]] Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            0.25};
}

class AlgebraicVerifier final : public ProjectionVerifier<double> {
public:
    mutable unsigned calls = 0;
    [[nodiscard]] core::Result<ProjectionAssessment>
    verify(const PrecisionBuffer<double>& source, const DenseReal& candidate,
           const ProjectionContext& identity) const override {
        ++calls;
        CHECK(identity.source == source.generation());
        CHECK(identity.target == core::NumericalGeneration{18});
        CHECK(identity.consumer == core::ScopeId{31});
        CHECK(std::ranges::equal(candidate.values(), std::array{1., 1.}));
        // This fixture requests only the algebraic residual checked by the
        // projection entry; it does not invent a physical constraint consumer.
        return ProjectionAssessment::Verified;
    }
};
} // namespace

TEST_CASE("Frozen rational determinant survives where its binary64 matrix loses rank",
          "[N4-01][N4-05][N4-10][verification][negative]") {
    // Frozen precision_hash.cpp uses t=2^-180 and det[[1,1],[1,1+t]].
    // This explicit finite determinant is an independent reference, not a
    // production arbitrary-rational rank routine or a geometric rank certificate.
    const Rational reference_tiny =
        Rational(1) / Rational(boost::multiprecision::cpp_int(1) << 180);
    const Rational one = 1;
    const Rational determinant = one * (one + reference_tiny) - one * one;
    CHECK(determinant == reference_tiny);
    CHECK(determinant > 0);
    CHECK(Rational(tiny) == reference_tiny);
    const double collapsed_entry = (one + reference_tiny).convert_to<double>();
    REQUIRE(collapsed_entry == 1);
    const Rational collapsed_determinant = Rational(collapsed_entry) - one;
    CHECK(collapsed_determinant == 0);

    core::ResourcePlan inputs{1048576}, work{1048576};
    const std::array entries{1., 1., 1., collapsed_entry};
    const auto collapsed = matrix({2, 2, 2}, entries, inputs);
    const auto original_context = context();
    const DenseOperator original{collapsed, original_context};
    const std::array scales{1., 1.};
    const std::array<UnitDimensions, 2> units{};
    const Scales fixed{scales, scales, units, units, original_context};
    const Limits limits{2, 1, 0, 5, EstimateAdmission{1048576}};
    const auto rejected = DensePrepared::prepare(original, fixed, limits, work);
    REQUIRE_FALSE(rejected.has_value());
    // An exact zero LU pivot is unresolved rank; no small-pivot cutoff or
    // certificate of the original rational matrix is inferred from this failure.
    CHECK(rejected.error().code == core::FailureCode::RankUnresolved);
    CHECK(rejected.error().operation == operation);
    CHECK(std::ranges::equal(collapsed.values(), entries));
    CHECK(original.context == original_context);
}

TEST_CASE("Frozen rank-collapse low word is consumed and its unsafe projection is rejected",
          "[N4-01][N4-06][N4-10][verification][negative]") {
    core::ResourcePlan inputs{1048576}, work{1048576};
    const std::array high_values{1., 1.}, low_values{0., tiny};
    auto source =
        PrecisionBuffer<double>::two_word(matrix({2, 1, 2}, high_values, inputs),
                                          matrix({2, 1, 2}, low_values, inputs), source_generation);
    REQUIRE(source.has_value());
    const auto& retained = std::get<TwoWord<double>>(source->storage());
    const Rational second = Rational(retained.hi.values()[1]) + Rational(retained.lo.values()[1]);
    CHECK(second == Rational(1) + Rational(tiny));
    CHECK(second - Rational(retained.hi.values()[0]) == Rational(tiny));

    // The second row of the fixed rational matrix is x=(1,1+t). Its exact
    // row difference C*x=t, C=[-1,1], is also the explicit determinant above.
    const auto difference = matrix({1, 2, 1}, std::array{-1., 1.}, inputs);
    const auto expected = matrix({1, 1, 1}, std::array{tiny}, inputs);
    const auto residual = original_residual(difference, *source, expected, work);
    REQUIRE(residual.has_value());
    CHECK(residual->values()[0] == 0);

    auto rounded = PrecisionBuffer<double>::binary64(matrix({2, 1, 2}, high_values, inputs),
                                                     core::NumericalGeneration{19});
    REQUIRE(rounded.has_value());
    const auto lost_residual = original_residual(difference, *rounded, expected, work);
    REQUIRE(lost_residual.has_value());
    CHECK(lost_residual->values()[0] == tiny);
    // N4-10: a fixed 1e-12 small-system budget. The physical residual scale is
    // the known input t, never a scale chosen from the computed residual.
    const std::array residual_scale{tiny};
    CHECK(Rational(lost_residual->values()[0]) / Rational(residual_scale[0]) == 1);
    const AlgebraicVerifier verifier;
    const auto projected =
        project_binary64(*source, difference, expected, residual_scale, 1e-12,
                         core::NumericalGeneration{18}, core::ScopeId{31}, verifier, work);
    REQUIRE_FALSE(projected.has_value());
    CHECK(projected.error().code == core::FailureCode::PrecisionUnresolved);
    CHECK(projected.error().operation == operation);
    CHECK(verifier.calls == 0);

    // A zero low plane is a positive control: these same high words satisfy
    // C*x=0 and may be projected. An implementation that rejects everything fails.
    auto exact = PrecisionBuffer<double>::two_word(matrix({2, 1, 2}, high_values, inputs),
                                                   matrix({2, 1, 2}, std::array{0., 0.}, inputs),
                                                   source_generation);
    REQUIRE(exact.has_value());
    const auto zero = matrix({1, 1, 1}, std::array{0.}, inputs);
    const auto accepted =
        project_binary64(*exact, difference, zero, residual_scale, 1e-12,
                         core::NumericalGeneration{18}, core::ScopeId{31}, verifier, work);
    REQUIRE(accepted.has_value());
    CHECK(verifier.calls == 1);
    CHECK(std::ranges::equal(accepted->values().values(), high_values));
    CHECK(source->generation() == source_generation);
    CHECK(std::ranges::equal(retained.hi.values(), high_values));
    CHECK(std::ranges::equal(retained.lo.values(), low_values));
    CHECK(std::ranges::equal(difference.values(), std::array{-1., 1.}));
    CHECK(expected.values()[0] == tiny);
}
