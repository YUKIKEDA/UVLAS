#include "sparse_test_support.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <bit>
#include <cmath>
#include <limits>

using namespace sparse_test;

TEST_CASE("Sparse product solves KKT with original constraints and support reaction",
          "[N4-05][N4-03]") {
    Fixture fixture{3};
    auto matrix = fixture.matrix(3, kkt);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{4.0, -1.0, 0.0, 8.0, -2.0, 0.0};
    auto b = fixture.rhs(3, 2, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    const std::array checks{
        ConsumerCheck{CheckKind::Constraint, core::ObservationId{50}, core::ScopeId{51}, 1e-12},
        ConsumerCheck{CheckKind::Reaction, core::ObservationId{52}, core::ScopeId{53}, 1e-12}};
    auto result = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(checks),
                                    accuracy(), limits(), prepared, fixture.plan));
    const std::array expected{1.0, -1.0, 2.0, 2.0, -2.0, 4.0};
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(std::abs(values(result)[i] - expected[i]) <= 1e-12);
    }
    CHECK(result.certificate().consumer_count() == 2);
    CHECK(fixture.evaluator.calls == 2);
    CHECK(result.certificate().rhs_generation() == core::RhsGeneration{41});
    CHECK(result.certificate().columns().size() == 2);
    CHECK(result.certificate().context() == context());
    auto repeated = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                      accuracy(), limits(), prepared, fixture.plan));
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(std::bit_cast<std::uint64_t>(values(result)[i]) ==
              std::bit_cast<std::uint64_t>(values(repeated)[i]));
    }
    fixture.evaluator.reject = true;
    const auto failed = solve_linear(original, rhs, fixture.scales(), fixture.precision(checks),
                                     accuracy(), limits(), prepared, fixture.plan);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().code == core::FailureCode::PrecisionUnresolved);
    CHECK(matrix.values()[0] == 2);
    CHECK(values(result)[0] == values(repeated)[0]);
}

TEST_CASE("Sparse SPD product verifies the original matrix before LLT", "[N4-05]") {
    Fixture fixture{3};
    const std::array<Triplet, 7> entries{
        {{0, 0, 2}, {1, 0, 4}, {2, 0, -6}, {0, 1, 4}, {1, 1, 11}, {0, 2, -6}, {2, 2, 71}}};
    auto matrix = fixture.matrix(3, entries);
    const SparseOperator original{matrix, context()};
    auto diagnostic = take(inspect_sparse_spd(original, limits(), fixture.plan));
    CHECK(diagnostic.status() == SparseSpdStatus::PositiveDefinite);
    CHECK(diagnostic.witness() == nullptr);
    auto prepared = take(SparsePrepared::prepare(
        original, fixture.scales(), SparseMethod::PositiveDefiniteLlt, limits(), fixture.plan));
    const std::array data{-24.0, -18.0, 207.0};
    auto b = fixture.rhs(3, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    auto result = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                    accuracy(Norm::L2), limits(), prepared, fixture.plan));
    const std::array expected{1.0, -2.0, 3.0};
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(std::abs(values(result)[i] - expected[i]) <= 1e-12);
    }
}

TEST_CASE("Sparse diagnostics distinguish proved negative and zero witnesses", "[N4-05]") {
    Fixture fixture{2};
    const std::array<Triplet, 4> negative{{{0, 0, 1}, {1, 0, 2}, {0, 1, 2}, {1, 1, 1}}};
    auto indefinite = fixture.matrix(2, negative);
    const SparseOperator original{indefinite, context()};
    auto report = take(inspect_sparse_spd(original, limits(), fixture.plan));
    CHECK(report.status() == SparseSpdStatus::NegativeCurvature);
    REQUIRE(report.witness() != nullptr);
    REQUIRE(report.quadratic_upper().has_value());
    if (report.quadratic_upper().has_value()) {
        CHECK(report.quadratic_upper().value() < 0);
    }
    const auto rejected = SparsePrepared::prepare(
        original, fixture.scales(), SparseMethod::PositiveDefiniteLlt, limits(), fixture.plan);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::RankUnresolved);
    const std::array<Triplet, 4> zero{{{0, 0, 1}, {1, 0, -1}, {0, 1, -1}, {1, 1, 1}}};
    auto rigid = fixture.matrix(2, zero);
    auto null = take(inspect_sparse_spd(SparseOperator{rigid, context()}, limits(), fixture.plan));
    CHECK(null.status() == SparseSpdStatus::ZeroMode);
    REQUIRE(null.witness() != nullptr);
    CHECK(null.witness()->values()[0] != 0);
}

TEST_CASE("Sparse cache matches actual contents and every operator context", "[N4-03][N4-05]") {
    Fixture fixture{3};
    auto matrix = fixture.matrix(3, kkt);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{4.0, -1.0, 0.0};
    auto b = fixture.rhs(3, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    auto changed_entries = kkt;
    changed_entries[0].value = 4;
    auto changed = fixture.matrix(3, changed_entries);
    auto mismatch = solve_linear(SparseOperator{changed, context()}, rhs, fixture.scales(),
                                 fixture.precision(), accuracy(), limits(), prepared, fixture.plan);
    REQUIRE_FALSE(mismatch);
    CHECK(mismatch.error().code == core::FailureCode::StaleGeneration);
    auto identity = context();
    identity.cache.matrix = core::NumericalGeneration{999};
    auto epoch = solve_linear(SparseOperator{matrix, identity}, rhs, fixture.scales(),
                              fixture.precision(), accuracy(), limits(), prepared, fixture.plan);
    REQUIRE_FALSE(epoch);
    CHECK(epoch.error().code == core::FailureCode::StaleGeneration);
    const Rhs<double> new_rhs{b, context(), core::RhsGeneration{123}};
    auto solved = take(solve_linear(original, new_rhs, fixture.scales(), fixture.precision(),
                                    accuracy(), limits(), prepared, fixture.plan));
    CHECK(solved.certificate().rhs_generation() == core::RhsGeneration{123});
}

TEST_CASE("Sparse retains stalled corrections and keeps strict backend admission explicit",
          "[N4-05][N4-06]") {
    Fixture fixture{1};
    const std::array<Triplet, 1> entries{{{0, 0, 3}}};
    auto matrix = fixture.matrix(1, entries);
    const SparseOperator original{matrix, context()};
    auto strict = limits();
    strict.backend_memory = RequireProvenBound{1000000};
    auto no_bound = SparsePrepared::prepare(original, fixture.scales(), SparseMethod::GeneralLu,
                                            strict, fixture.plan);
    REQUIRE_FALSE(no_bound);
    CHECK(no_bound.error().code == core::FailureCode::Unsupported);
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array high{1.0};
    auto b = fixture.rhs(1, 1, high);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    auto demand = accuracy();
    demand.absolute_tolerance = 1e-30;
    demand.componentwise_backward_tolerance = 1e-30;
    auto retained = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(), demand,
                                      limits(), prepared, fixture.plan));
    REQUIRE(retained.values().precision() == core::Precision::TwoWord);
    const auto& pair = std::get<TwoWord<double>>(retained.values().storage());
    CHECK(pair.lo.values()[0] != 0);
    using Big = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<200>>;
    const Big error = Big{1} - Big{3} * (Big{pair.hi.values()[0]} + Big{pair.lo.values()[0]});
    CHECK(abs(error) <= Big{"1e-30"});
    auto projection = fixture.precision();
    projection.returned =
        Binary64Projection{core::NumericalGeneration{100}, core::NumericalGeneration{101}};
    auto rejected = solve_linear(original, rhs, fixture.scales(), projection, demand, limits(),
                                 prepared, fixture.plan);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::PrecisionUnresolved);
}

TEST_CASE("Sparse empty and zero RHS values qualify without invalid divisions", "[N4-05]") {
    Fixture fixture{0};
    const std::array<Triplet, 0> none{};
    auto matrix = fixture.matrix(0, none);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    auto b = fixture.rhs(0, 2, {});
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    auto result = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                    accuracy(), limits(), prepared, fixture.plan));
    CHECK(result.values().shape().rows == 0);
    CHECK(result.values().shape().columns == 2);
    CHECK(result.certificate().columns().size() == 2);
}

TEST_CASE("Sparse projection and zero-column RHS retain exact identities", "[N4-05][N4-06]") {
    Fixture fixture{3};
    auto matrix = fixture.matrix(3, kkt);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{4.0, -1.0, 0.0};
    auto b = fixture.rhs(3, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    auto request = fixture.precision();
    request.returned =
        Binary64Projection{core::NumericalGeneration{100}, core::NumericalGeneration{101}};
    auto projected = take(solve_linear(original, rhs, fixture.scales(), request, accuracy(),
                                       limits(), prepared, fixture.plan));
    CHECK(projected.values().generation() == core::NumericalGeneration{101});
    CHECK(std::abs(values(projected)[0] - 1) <= 1e-12);
    auto empty = fixture.rhs(3, 0, {});
    const Rhs<double> zero{empty, context(), core::RhsGeneration{42}};
    auto solved = take(solve_linear(original, zero, fixture.scales(), fixture.precision(),
                                    accuracy(), limits(), prepared, fixture.plan));
    CHECK(solved.values().shape().columns == 0);
    CHECK(solved.certificate().columns().empty());
    const std::array zeros{0.0, 0.0, 0.0};
    auto zero_buffer = fixture.rhs(3, 1, zeros);
    const Rhs<double> zero_rhs{zero_buffer, context(), core::RhsGeneration{43}};
    auto exact = accuracy();
    exact.absolute_tolerance = 0;
    exact.relative_tolerance = 0;
    exact.componentwise_backward_tolerance = 0;
    auto zero_solution =
        take(solve_linear(original, zero_rhs, fixture.scales(), fixture.precision(), exact,
                          limits(), prepared, fixture.plan));
    CHECK(std::ranges::all_of(values(zero_solution), [](double value) { return value == 0; }));
}

TEST_CASE("Sparse failure leaves previous solutions and reusable factors intact",
          "[N4-05][N4-03]") {
    Fixture fixture{3};
    auto matrix = fixture.matrix(3, kkt);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{4.0, -1.0, 0.0};
    auto b = fixture.rhs(3, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    const std::array checks{
        ConsumerCheck{CheckKind::Constraint, core::ObservationId{50}, core::ScopeId{51}, 1e-12}};
    fixture.evaluator.wrong_generation = true;
    auto stale = solve_linear(original, rhs, fixture.scales(), fixture.precision(checks),
                              accuracy(), limits(), prepared, fixture.plan);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == core::FailureCode::StaleGeneration);
    fixture.evaluator.wrong_generation = false;
    core::ResourcePlan exhausted{0};
    auto unavailable = solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                    accuracy(), limits(), prepared, exhausted);
    REQUIRE_FALSE(unavailable);
    CHECK(unavailable.error().code == core::FailureCode::ResourceExceeded);
    auto recovered = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                       accuracy(), limits(), prepared, fixture.plan));
    CHECK(std::abs(values(recovered)[2] - 2) <= 1e-12);
    auto moved = std::move(prepared);
    auto inactive = solve_linear(original, rhs, fixture.scales(), fixture.precision(), accuracy(),
                                 limits(), prepared, fixture.plan);
    REQUIRE_FALSE(inactive);
    CHECK(inactive.error().code == core::FailureCode::StaleGeneration);
    auto reused = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                    accuracy(), limits(), moved, fixture.plan));
    CHECK(values(recovered)[2] == values(reused)[2]);
}

TEST_CASE("Sparse LU rejects singular constraints without a small-pivot cutoff", "[N4-05]") {
    Fixture fixture{2};
    {
        Fixture rank_fixture{3};
        const std::array<Triplet, 5> redundant_kkt{
            {{0, 0, 2}, {1, 0, 1}, {2, 0, 2}, {0, 1, 1}, {0, 2, 2}}};
        auto rank_matrix = rank_fixture.matrix(3, redundant_kkt);
        auto rank_failure =
            SparsePrepared::prepare(SparseOperator{rank_matrix, context()}, rank_fixture.scales(),
                                    SparseMethod::GeneralLu, limits(), rank_fixture.plan);
        REQUIRE_FALSE(rank_failure);
        CHECK(rank_failure.error().code == core::FailureCode::RankUnresolved);
    }
    const std::array<Triplet, 4> redundant{{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}};
    auto matrix = fixture.matrix(2, redundant);
    auto rejected = SparsePrepared::prepare(SparseOperator{matrix, context()}, fixture.scales(),
                                            SparseMethod::GeneralLu, limits(), fixture.plan);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::RankUnresolved);
    const double tiny = std::ldexp(1.0, -500);
    const std::array<Triplet, 2> entries{{{0, 0, 1}, {1, 1, tiny}}};
    auto small = fixture.matrix(2, entries);
    const SparseOperator original{small, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{2.0, tiny};
    auto b = fixture.rhs(2, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    auto solved = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                    accuracy(), limits(), prepared, fixture.plan));
    CHECK(values(solved)[0] == 2);
    CHECK(values(solved)[1] == 1);
}

TEST_CASE("Sparse consumer reentry refuses workspace reuse without damaging outer result",
          "[N4-03][N4-05]") {
    Fixture fixture{3};
    auto matrix = fixture.matrix(3, kkt);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{4.0, -1.0, 0.0};
    auto b = fixture.rhs(3, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    class Reentry final : public ConsumerEvaluator<double> {
    public:
        Reentry(const SparseOperator& a, const Rhs<double>& b, Fixture& f, SparsePrepared& p)
            : original(a), rhs(b), fixture(f), prepared(p) {}
        const SparseOperator& original;
        const Rhs<double>& rhs;
        Fixture& fixture;
        SparsePrepared& prepared;
        mutable bool refused = false;
        SupportedPrecision supported_precision() const noexcept override {
            return SupportedPrecision::Binary64Only;
        }
        core::Result<core::ErrorCertificate>
        assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>&, const Context&,
               const ConsumerCheck& request, core::CertificateContext expected) const override {
            auto nested = solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                       accuracy(), limits(), prepared, fixture.plan);
            refused = !nested && nested.error().code == core::FailureCode::InvalidInput;
            return core::ErrorCertificate::assess(expected, request.budget, 0);
        }
    };
    const Reentry evaluator{original, rhs, fixture, prepared};
    const std::array checks{
        ConsumerCheck{CheckKind::Observation, core::ObservationId{50}, core::ScopeId{51}, 1e-12}};
    const Precision<double> request{RetainedPrecision{core::NumericalGeneration{100}}, checks,
                                    evaluator};
    auto result = take(solve_linear(original, rhs, fixture.scales(), request, accuracy(), limits(),
                                    prepared, fixture.plan));
    CHECK(evaluator.refused);
    CHECK(std::abs(values(result)[0] - 1) <= 1e-12);
    auto next = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(), accuracy(),
                                  limits(), prepared, fixture.plan));
    CHECK(values(result)[0] == values(next)[0]);
}

TEST_CASE("Sparse TwoWord RHS low plane affects retained solution and consumer proof",
          "[N4-05][N4-06]") {
    Fixture fixture{1};
    const std::array<Triplet, 1> entries{{{0, 0, 1}}};
    auto matrix = fixture.matrix(1, entries);
    const SparseOperator original{matrix, context()};
    const std::array high{1.0}, low{std::ldexp(1.0, -80)};
    class LowWordEvaluator final : public ConsumerEvaluator<double> {
    public:
        explicit LowWordEvaluator(double expected) : expected_low(expected) {}
        double expected_low;
        mutable unsigned calls = 0;
        SupportedPrecision supported_precision() const noexcept override {
            return SupportedPrecision::Binary64AndTwoWord;
        }
        core::Result<core::ErrorCertificate>
        assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>& candidate,
               const Context&, const ConsumerCheck& request,
               core::CertificateContext expected) const override {
            ++calls;
            const auto* pair = std::get_if<TwoWord<double>>(&candidate.storage());
            const double upper =
                pair && pair->lo.values()[0] == expected_low ? 0 : 2 * request.budget;
            return core::ErrorCertificate::assess(expected, request.budget, upper);
        }
    };
    const LowWordEvaluator evaluator{low[0]};
    const std::array checks{
        ConsumerCheck{CheckKind::Observation, core::ObservationId{55}, core::ScopeId{56}, 1e-30}};
    const Precision<double> request{RetainedPrecision{core::NumericalGeneration{100}}, checks,
                                    evaluator};

    for (const auto method : {SparseMethod::GeneralLu, SparseMethod::PositiveDefiniteLlt}) {
        auto prepared = take(
            SparsePrepared::prepare(original, fixture.scales(), method, limits(), fixture.plan));
        auto hi = take(DenseReal::copy({1, 1, 1}, high, fixture.plan, operation));
        auto lo = take(DenseReal::copy({1, 1, 1}, low, fixture.plan, operation));
        auto input = take(PrecisionBuffer<double>::two_word(std::move(hi), std::move(lo),
                                                            core::NumericalGeneration{40}));
        const Rhs<double> rhs{input, context(), core::RhsGeneration{41}};
        auto demand = accuracy();
        demand.absolute_tolerance = 1e-30;
        demand.componentwise_backward_tolerance = 1e-30;
        auto retained = take(solve_linear(original, rhs, fixture.scales(), request, demand,
                                          limits(), prepared, fixture.plan));
        REQUIRE(retained.values().precision() == core::Precision::TwoWord);
        const auto& pair = std::get<TwoWord<double>>(retained.values().storage());
        CHECK(pair.hi.values()[0] == 1);
        CHECK(pair.lo.values()[0] == low[0]);
    }
    CHECK(evaluator.calls == 2);
}

TEST_CASE("Sparse last correction follows the requested L2 unknown scale norm", "[N4-05][N4-06]") {
    Fixture fixture{2};
    const std::array<Triplet, 2> entries{{{0, 0, 3}, {1, 1, 3}}};
    auto matrix = fixture.matrix(2, entries);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{1.0, 1.0};
    auto b = fixture.rhs(2, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    auto demand = accuracy(Norm::L2);
    demand.absolute_tolerance = 1e-30;
    demand.componentwise_backward_tolerance = 1e-30;
    auto result = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(), demand,
                                    limits(), prepared, fixture.plan));
    // Each initial correction is about 1.85e-17; their L2 norm exceeds 2e-17.
    const double correction = result.certificate().columns()[0].scaled_correction_upper;
    CHECK(correction > 2e-17);
    CHECK(correction < 3e-17);
}

TEST_CASE("Sparse scale generations cannot restamp different values or SI units",
          "[N4-03][N4-05]") {
    Fixture fixture{1};
    const std::array<Triplet, 1> entries{{{0, 0, 1}}};
    auto matrix = fixture.matrix(1, entries);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array data{2.0};
    auto b = fixture.rhs(1, 1, data);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    fixture.unknown[0] = 2;
    auto scale = solve_linear(original, rhs, fixture.scales(), fixture.precision(), accuracy(),
                              limits(), prepared, fixture.plan);
    REQUIRE_FALSE(scale);
    CHECK(scale.error().code == core::FailureCode::StaleGeneration);
    CHECK(scale.error().details.model == context().cache.model);
    CHECK(scale.error().details.trial == context().trial);
    CHECK(scale.error().details.numerical_generation == context().cache.matrix);
    fixture.unknown[0] = 1;
    fixture.residual_units[0][0] = 1;
    auto unit = solve_linear(original, rhs, fixture.scales(), fixture.precision(), accuracy(),
                             limits(), prepared, fixture.plan);
    REQUIRE_FALSE(unit);
    CHECK(unit.error().code == core::FailureCode::StaleGeneration);
    fixture.residual_units[0][0] = 0;
    auto recovered = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(),
                                       accuracy(), limits(), prepared, fixture.plan));
    CHECK(values(recovered)[0] == 2);
    fixture.residual[0] = 0;
    auto invalid = SparsePrepared::prepare(original, fixture.scales(), SparseMethod::GeneralLu,
                                           limits(), fixture.plan);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
}
