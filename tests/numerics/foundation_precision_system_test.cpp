#include "full_precision_fixture.hpp"
#include "full_precision_oracle.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <uvlas/numerics/dense_solve.hpp>
#include <uvlas/numerics/precision.hpp>
#include <variant>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
namespace reference = uvlas::test_support::full_precision;
using Complex = std::complex<double>;
constexpr auto operation = core::Operation::SolveLinear;
constexpr std::size_t order = 4096;
constexpr double relative_budget = 1e-11;
constexpr double backward_budget = 1e-12;
constexpr double observation_budget = 1e-8;
constexpr auto observation_scope = core::ScopeId{4096};
constexpr std::uint64_t payload_budget = std::uint64_t{1} << 30;

template <class T>
[[nodiscard]] T take(core::Result<T> result) {
    if (!result.has_value()) {
        FAIL("Unexpected product failure " << static_cast<unsigned>(result.error().code));
    }
    return std::move(result).value();
}
template <StorageScalar T>
[[nodiscard]] DenseMatrix<T> matrix(core::Index rows, core::Index columns,
                                    std::span<const T> values, core::ResourcePlan& plan) {
    return take(DenseMatrix<T>::copy({rows, columns, rows}, values, plan, operation));
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
[[nodiscard]] std::vector<Complex> column(std::span<const Complex> row_major, std::size_t index) {
    if (row_major.size() != order * 3 || index >= 3)
        throw std::runtime_error("Unexpected full-system column shape");
    std::vector<Complex> result(order);
    for (std::size_t i = 0; i < order; ++i)
        result[i] = row_major[i * 3 + index];
    return result;
}
struct Words {
    std::span<const Complex> high, low;
};
[[nodiscard]] Words words(const PrecisionBuffer<Complex>& value) {
    if (const auto* single = std::get_if<Binary64<Complex>>(&value.storage()))
        return {single->values.values(), {}};
    const auto& pair = std::get<TwoWord<Complex>>(value.storage());
    return {pair.hi.values(), pair.lo.values()};
}
[[nodiscard]] std::vector<Complex> weights(const reference::Fixture& fixture,
                                           const reference::Observation& observation) {
    if (observation.column >= 3)
        throw std::runtime_error("Unknown observation column");
    std::vector<Complex> result(order);
    const auto frequency = 16 * fixture.k;
    if (observation.kind == "force") {
        for (std::size_t i = 0; i < order; ++i)
            result[i] = {fixture.fc[i], frequency * fixture.fu[i]};
    } else if (observation.kind == "moment") {
        const auto pivot = fixture.modes[observation.column].xp;
        for (std::size_t i = 0; i < order; ++i)
            result[i] = {fixture.mc0[i] + pivot * fixture.fc[i],
                         frequency * (fixture.mu0[i] + pivot * fixture.fu[i])};
    } else if (observation.kind == "pressure") {
        const auto lower = observation.lower, upper = observation.upper;
        if (!(upper > lower))
            throw std::runtime_error("Invalid fixed pressure window");
        const auto left = std::upper_bound(fixture.xv.begin(), fixture.xv.end(), lower);
        const auto right = std::lower_bound(fixture.xv.begin(), fixture.xv.end(), upper);
        if (left == fixture.xv.begin() || right == fixture.xv.begin())
            throw std::runtime_error("Pressure window outside the fixed field");
        const auto left_index = static_cast<std::size_t>(left - fixture.xv.begin() - 1);
        const auto right_index = static_cast<std::size_t>(right - fixture.xv.begin() - 1);
        if (left_index >= order || right_index >= order)
            throw std::runtime_error("Pressure trace outside the fixed field");
        for (std::size_t i = 0; i < order; ++i) {
            const auto overlap =
                std::max(0., std::min(fixture.xv[i + 1], upper) - std::max(fixture.xv[i], lower));
            result[i] = {0, frequency * overlap};
        }
        result[left_index] -= Complex{8, 0};
        result[right_index] += Complex{8, 0};
        for (auto& value : result)
            value *= 1.2 / (upper - lower);
    } else {
        throw std::runtime_error("Unknown fixed observation kind");
    }
    return result;
}

// This is a fixture-specific linear map, not a product UVLM load kernel. Its
// reference is the frozen R&D observation, independent of the current solution.
class Observations final : public ConsumerEvaluator<Complex> {
public:
    Observations(const reference::Fixture& fixture, std::size_t selected)
        : fixture_(fixture), selected_(selected) {
        for (std::size_t i = 0; i < fixture_.observations.size(); ++i) {
            if (fixture_.observations[i].column == selected_) {
                indices_.push_back(i);
                coefficients_.push_back(weights(fixture_, fixture_.observations[i]));
                checks_.push_back({CheckKind::Observation, core::ObservationId{1000 + i},
                                   observation_scope, observation_budget});
            }
        }
        if (indices_.size() != 13)
            throw std::runtime_error("Incomplete force/moment/pressure observation group");
        for (std::size_t i = 0; i < indices_.size(); ++i) {
            const auto& observation = fixture_.observations[indices_[i]];
            const auto expected_kind = i == 0 ? "force" : i == 1 ? "moment" : "pressure";
            if (observation.kind != expected_kind ||
                observation.mode != fixture_.modes[selected_].name)
                throw std::runtime_error("Changed frozen observation identity or order");
        }
    }
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] std::span<const ConsumerCheck> checks() const noexcept {
        return checks_;
    }
    [[nodiscard]] unsigned calls(std::size_t index) const {
        return calls_.at(index);
    }
    [[nodiscard]] double error(std::size_t local, Words value) const {
        const auto& expected = fixture_.observations[indices_.at(local)];
        return reference::observation_error_upper(coefficients_.at(local), value.high, value.low,
                                                  expected.value, expected.scale);
    }
    [[nodiscard]] bool all_qualified(Words candidate) const {
        for (std::size_t i = 0; i < indices_.size(); ++i) {
            const auto upper = error(i, candidate);
            if (!std::isfinite(upper) || upper > observation_budget)
                return false;
        }
        return true;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<Complex>&, const PrecisionBuffer<Complex>& candidate,
           const Context&, const ConsumerCheck& request,
           core::CertificateContext expected) const override {
        if (request.kind != CheckKind::Observation || request.scope != observation_scope)
            return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
        for (std::size_t i = 0; i < checks_.size(); ++i) {
            if (request.observation == checks_[i].observation) {
                ++calls_[i];
                return core::ErrorCertificate::assess(expected, request.budget,
                                                      error(i, words(candidate)));
            }
        }
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    }

private:
    const reference::Fixture& fixture_;
    std::size_t selected_;
    std::vector<std::size_t> indices_;
    std::vector<std::vector<Complex>> coefficients_;
    std::vector<ConsumerCheck> checks_;
    mutable std::array<unsigned, 13> calls_{};
};

void check_audit(const reference::Audit& audit) {
    CHECK(audit.l2_squared_ratio <=
          reference::Rational(relative_budget) * reference::Rational(relative_budget));
    CHECK(audit.backward_upper <= reference::Rational(backward_budget));
}
[[nodiscard]] double relative_l2(const reference::Audit& audit) {
    // Reporting only: acceptance above compares the exact squared quantities.
    return std::sqrt(audit.l2_squared_ratio.convert_to<double>());
}
[[nodiscard]] double maximum_observation_error(const Observations& observations, Words value) {
    double result = 0;
    for (std::size_t i = 0; i < observations.checks().size(); ++i)
        result = std::max(result, observations.error(i, value));
    return result;
}
class FullProjection final : public ProjectionVerifier<Complex> {
public:
    FullProjection(std::span<const double> base, std::span<const Complex> update,
                   std::span<const Complex> rhs, const Observations& observations)
        : base_(base), update_(update), rhs_(rhs), observations_(observations) {}
    [[nodiscard]] core::Result<ProjectionAssessment>
    verify(const PrecisionBuffer<Complex>&, const DenseComplex& candidate,
           const ProjectionContext&) const override {
        ++calls;
        const auto audit =
            reference::audit_full_column(base_, update_, rhs_, candidate.values(), {});
        // The component envelope is a conservative sufficient component test;
        // passing it does not establish the sum over the whole column. This
        // callback separately evaluates the required global L2 condition.
        if (audit.l2_squared_ratio >
                reference::Rational(relative_budget) * reference::Rational(relative_budget) ||
            audit.backward_upper > reference::Rational(backward_budget) ||
            !observations_.all_qualified({candidate.values(), {}}))
            return std::unexpected(
                core::Failure{core::FailureCode::PrecisionUnresolved, operation});
        return ProjectionAssessment::Verified;
    }
    mutable unsigned calls = 0;

private:
    std::span<const double> base_;
    std::span<const Complex> update_, rhs_;
    const Observations& observations_;
};

[[nodiscard]] DenseComplex effective_matrix(std::span<const double> base,
                                            std::span<const Complex> update,
                                            core::ResourcePlan& plan) {
    std::vector<Complex> coefficients(order * order);
    for (std::size_t j = 0; j < order; ++j)
        for (std::size_t i = 0; i < order; ++i)
            coefficients[j * order + i] = j + 1 == order
                                              ? Complex{base[j * order + i], 0} + update[i]
                                              : Complex{base[j * order + i], 0};
    return matrix<Complex>(order, order, coefficients, plan);
}
[[nodiscard]] bool unchanged_effective(const DenseComplex& actual, std::span<const double> base,
                                       std::span<const Complex> update) {
    for (std::size_t j = 0; j < order; ++j)
        for (std::size_t i = 0; i < order; ++i) {
            const auto expected = j + 1 == order ? Complex{base[j * order + i], 0} + update[i]
                                                 : Complex{base[j * order + i], 0};
            if (actual.values()[j * order + i] != expected)
                return false;
        }
    return true;
}
} // namespace

TEST_CASE("Full 4096 updated product qualifies every returned RHS and frozen observation",
          "[N4-05][N4-06][N4-10][verification][negative]") {
    const auto fixture = reference::load_fixture(UVLAS_FULL_PRECISION_FIXTURE);
    REQUIRE(fixture.m == static_cast<core::Index>(order));
    const auto base = reference::reconstruct_base(fixture);
    core::ResourcePlan inputs{payload_budget}, workspace{payload_budget};
    const auto original = matrix<double>(order, order, base, inputs);
    const auto update = matrix<Complex>(order, 1, fixture.update, inputs);
    std::vector<Complex> selector(order);
    selector.back() = {1, 0};
    const auto transpose_factor = matrix<Complex>(order, 1, selector, inputs);
    const ColumnUpdatedOperator op{original, update, transpose_factor, context()};
    const std::vector<double> ones(order, 1.);
    const std::vector<UnitDimensions> units(order);
    const Scales scales{ones, ones, units, units, context()};
    const Limits limits{order, 1, 1, 5, EstimateAdmission{payload_budget}};
    auto prepared = take(UpdatedPrepared::prepare(op, scales, limits, workspace));

    for (std::size_t j = 0; j < 3; ++j) {
        CAPTURE(j, fixture.modes[j].name);
        const auto rhs_values = column(fixture.rhs, j);
        const auto rhs_buffer = take(PrecisionBuffer<Complex>::binary64(
            matrix<Complex>(order, 1, rhs_values, inputs), core::NumericalGeneration{20 + j}));
        const Rhs<Complex> rhs{rhs_buffer, context(), core::RhsGeneration{100 + j}};
        const Observations observations{fixture, j};
        const Precision<Complex> precision{RetainedPrecision{core::NumericalGeneration{40 + j}},
                                           observations.checks(), observations};
        const Accuracy accuracy{Norm::L2,
                                reference::conservative_l2_allowance(rhs_values, relative_budget),
                                0,
                                backward_budget,
                                core::ObservationId{50},
                                core::ObservationId{51}};
        core::ResourcePlan output{payload_budget};
        auto result =
            take(solve_linear(op, rhs, scales, precision, accuracy, limits, prepared, output));
        const auto actual = words(result.values());
        const auto audit =
            reference::audit_full_column(base, fixture.update, rhs_values, actual.high, actual.low);
        check_audit(audit);
        std::cout << std::setprecision(17) << "full4096 product column=" << j
                  << " precision=" << (actual.low.empty() ? "Binary64" : "TwoWord")
                  << " relative_l2=" << relative_l2(audit)
                  << " backward_upper=" << reference::to_upper(audit.backward_upper)
                  << " maximum_observation_error="
                  << maximum_observation_error(observations, actual) << '\n';
        // Preserve each completed column in the log during the long full-system run.
        std::cout << std::flush;
        CHECK(result.certificate().consumer_count() == 13);
        for (std::size_t i = 0; i < 13; ++i)
            CHECK(observations.calls(i) > 0);

        // A failure after a successful solve may not mutate the old result or RHS.
        const std::vector<Complex> old_high(actual.high.begin(), actual.high.end());
        const std::vector<Complex> old_low(actual.low.begin(), actual.low.end());
        // A deliberately changed trailing coefficient changes the force map.
        // The consumer must assess this actual candidate, not return a cached
        // Achieved result from the successful solve above.
        auto altered_high = old_high;
        altered_high.back() += Complex{1, 0};
        auto altered_low = old_low;
        if (altered_low.empty())
            altered_low.resize(order);
        const auto altered = take(PrecisionBuffer<Complex>::two_word(
            matrix<Complex>(order, 1, altered_high, inputs),
            matrix<Complex>(order, 1, altered_low, inputs), core::NumericalGeneration{60 + j}));
        const auto& force_check = observations.checks().front();
        const auto missed = take(observations.assess(
            result.values(), altered, context(), force_check,
            core::CertificateContext{force_check.observation, force_check.scope,
                                     altered.generation(), core::Precision::TwoWord}));
        CHECK(std::holds_alternative<core::Exceeded>(missed.assessment()));
        core::ResourcePlan none{0};
        const auto failed =
            solve_linear(op, rhs, scales, precision, accuracy, limits, prepared, none);
        REQUIRE_FALSE(failed.has_value());
        CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
        CHECK(std::ranges::equal(words(result.values()).high, old_high));
        CHECK(std::ranges::equal(words(result.values()).low, old_low));
        CHECK(std::ranges::equal(words(rhs_buffer).high, rhs_values));
        CHECK(rhs_buffer.generation() == core::NumericalGeneration{20 + j});
    }
    CHECK(reference::base_hash(original.values(), order) == fixture.base_sha256);
    CHECK(std::ranges::equal(update.values(), fixture.update));
    CHECK(std::ranges::equal(transpose_factor.values(), selector));
    CHECK(op.context == context());
}

TEST_CASE("Full frozen TwoWord reference preserves its global residual and projection outcomes",
          "[N4-05][N4-06][N4-10][verification][negative]") {
    const auto fixture = reference::load_fixture(UVLAS_FULL_PRECISION_FIXTURE);
    const auto base = reference::reconstruct_base(fixture);
    core::ResourcePlan inputs{payload_budget}, work{payload_budget};
    const auto effective = effective_matrix(base, fixture.update, inputs);
    for (std::size_t j = 0; j < 3; ++j) {
        CAPTURE(j, fixture.modes[j].name);
        const auto rhs_values = column(fixture.rhs, j);
        const auto high = column(fixture.reference_hi, j), low = column(fixture.reference_lo, j);
        const auto source = take(PrecisionBuffer<Complex>::two_word(
            matrix<Complex>(order, 1, high, inputs), matrix<Complex>(order, 1, low, inputs),
            core::NumericalGeneration{17}));
        const auto retained_audit =
            reference::audit_full_column(base, fixture.update, rhs_values, high, low);
        check_audit(retained_audit);
        const Observations observations{fixture, j};
        CHECK(observations.all_qualified(words(source)));

        std::vector<Complex> projected_words(order);
        for (std::size_t i = 0; i < order; ++i)
            projected_words[i] = high[i] + low[i];
        const auto projected_audit =
            reference::audit_full_column(base, fixture.update, rhs_values, projected_words, {});
        std::cout << std::setprecision(17) << "full4096 reference column=" << j
                  << " retained_relative_l2=" << relative_l2(retained_audit)
                  << " projected_relative_l2=" << relative_l2(projected_audit)
                  << " retained_backward_upper="
                  << reference::to_upper(retained_audit.backward_upper)
                  << " retained_maximum_observation_error="
                  << maximum_observation_error(observations, words(source)) << '\n';
        // Preserve each completed column in the log during the long full-system run.
        std::cout << std::flush;
        const auto exact_budget =
            reference::Rational(relative_budget) * reference::Rational(relative_budget);
        CHECK((projected_audit.l2_squared_ratio <= exact_budget) == (j == 2));

        const auto rhs = matrix<Complex>(order, 1, rhs_values, inputs);
        const std::vector<double> residual_scale(
            order, reference::conservative_l2_allowance(rhs_values, 1.));
        const FullProjection verifier{base, fixture.update, rhs_values, observations};
        const auto projected =
            project_binary64(source, effective, rhs, residual_scale, relative_budget,
                             core::NumericalGeneration{18}, core::ScopeId{4096}, verifier, work);
        if (j == 2) {
            REQUIRE(projected.has_value());
            CHECK(verifier.calls == 1);
            CHECK(std::ranges::equal(projected->values().values(), projected_words));
        } else {
            REQUIRE_FALSE(projected.has_value());
            CHECK(projected.error().code == core::FailureCode::PrecisionUnresolved);
        }
        CHECK(std::ranges::equal(words(source).high, high));
        CHECK(std::ranges::equal(words(source).low, low));
        CHECK(source.generation() == core::NumericalGeneration{17});
        CHECK(std::ranges::equal(rhs.values(), rhs_values));
    }
    CHECK(reference::base_hash(base, order) == fixture.base_sha256);
    CHECK(unchanged_effective(effective, base, fixture.update));
}
