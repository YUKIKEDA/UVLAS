#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <uvlas/numerics/dense_solve.hpp>
#include <uvlas/numerics/sparse_solve.hpp>
#include <variant>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
using Exact = boost::multiprecision::cpp_rational;
constexpr auto operation = core::Operation::SolveLinear;
enum class RealBackend : std::uint8_t { Dense, SparseLu, SparseLlt };
constexpr std::array original_a{7.0, -2.0, 6.0, 5.0, -5.0, 5.0, 8.0, 8.0, 6.0};
constexpr std::array original_b{41.0, 12.0, 34.0};
constexpr std::array exact_x{1.0, 2.0, 3.0};

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
template <StorageScalar Scalar>
class NoChecks final : public ConsumerEvaluator<Scalar> {
public:
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<Scalar>&, const PrecisionBuffer<Scalar>&, const Context&,
           const ConsumerCheck& request, core::CertificateContext expected) const override {
        FAIL_CHECK("No consumer was requested by this algebraic fixture");
        return core::ErrorCertificate::unevaluated(expected, request.budget);
    }
};

class RealSystem {
public:
    RealSystem(RealBackend backend, core::Index order, std::span<const double> coefficients)
        : backend_(backend), n_(order) {
        if (n_ < 1 || n_ > 3 || coefficients.size() != static_cast<std::size_t>(n_ * n_)) {
            throw std::invalid_argument("Fixture order must be in [1,3]");
        }
        auto copied = DenseReal::copy({n_, n_, n_}, coefficients, input_plan_, operation);
        REQUIRE(copied);
        matrix_ = std::move(*copied);
        std::array<Triplet, 9> triplets{};
        for (core::Index column = 0; column < n_; ++column) {
            for (core::Index row = 0; row < n_; ++row) {
                const auto index = static_cast<std::size_t>(column * n_ + row);
                triplets[index] = {row, column, coefficients[index]};
            }
        }
        auto sparse = SparseReal::from_triplets(
            n_, n_, std::span{triplets}.first(coefficients.size()), input_plan_, operation);
        REQUIRE(sparse);
        sparse_ = std::move(*sparse);
        if (backend_ == RealBackend::Dense) {
            auto prepared = DensePrepared::prepare(DenseOperator{matrix_, context()}, scales(),
                                                   limits(5), workspace_plan_);
            REQUIRE(prepared);
            dense_.emplace(std::move(*prepared));
        } else {
            auto prepared = SparsePrepared::prepare(SparseOperator{sparse_, context()}, scales(),
                                                    backend_ == RealBackend::SparseLu
                                                        ? SparseMethod::GeneralLu
                                                        : SparseMethod::PositiveDefiniteLlt,
                                                    limits(5), workspace_plan_);
            REQUIRE(prepared);
            prepared_sparse_.emplace(std::move(*prepared));
        }
    }
    [[nodiscard]] PrecisionBuffer<double> rhs(std::span<const double> values) {
        auto copied = DenseReal::copy({n_, 1, n_}, values, input_plan_, operation);
        REQUIRE(copied);
        auto result =
            PrecisionBuffer<double>::binary64(std::move(*copied), core::NumericalGeneration{30});
        REQUIRE(result);
        return std::move(*result);
    }
    [[nodiscard]] core::Result<Solution<double>>
    solve(const PrecisionBuffer<double>& rhs, double budget, unsigned refinements, Norm norm,
          double relative = 0, bool project = false, std::span<const ConsumerCheck> checks = {},
          const ConsumerEvaluator<double>* evaluator = nullptr) {
        const Accuracy accuracy{norm,
                                budget,
                                relative,
                                std::max(budget, relative),
                                core::ObservationId{20},
                                core::ObservationId{21}};
        const auto generation = core::NumericalGeneration{next_generation_++};
        ReturnedPrecision returned = RetainedPrecision{generation};
        if (project) {
            returned =
                Binary64Projection{generation, core::NumericalGeneration{next_generation_++}};
        }
        const auto& consumer = evaluator != nullptr
                                   ? *evaluator
                                   : static_cast<const ConsumerEvaluator<double>&>(consumer_);
        const Precision<double> precision{returned, checks, consumer};
        const Rhs<double> input{rhs, context(), core::RhsGeneration{31}};
        if (backend_ == RealBackend::Dense) {
            if (!dense_) {
                throw std::logic_error("Missing dense fixture preparation");
            }
            return solve_linear(DenseOperator{matrix_, context()}, input, scales(), precision,
                                accuracy, limits(refinements), *dense_, output_plan_);
        }
        if (!prepared_sparse_) {
            throw std::logic_error("Missing sparse fixture preparation");
        }
        return solve_linear(SparseOperator{sparse_, context()}, input, scales(), precision,
                            accuracy, limits(refinements), *prepared_sparse_, output_plan_);
    }
    [[nodiscard]] PrecisionBuffer<double> rhs_words(std::span<const double> hi,
                                                    std::span<const double> lo) {
        auto high = DenseReal::copy({n_, 1, n_}, hi, input_plan_, operation);
        auto low = DenseReal::copy({n_, 1, n_}, lo, input_plan_, operation);
        REQUIRE(high);
        REQUIRE(low);
        auto result = PrecisionBuffer<double>::two_word(std::move(*high), std::move(*low),
                                                        core::NumericalGeneration{30});
        REQUIRE(result);
        return std::move(*result);
    }
    [[nodiscard]] const DenseReal& original_dense() const noexcept {
        return matrix_;
    }
    [[nodiscard]] const SparseReal& original_sparse() const noexcept {
        return sparse_;
    }

private:
    [[nodiscard]] Scales scales() const {
        const auto n = static_cast<std::size_t>(n_);
        return {std::span{scales_}.first(n), std::span{scales_}.first(n),
                std::span{units_}.first(n), std::span{units_}.first(n), context()};
    }
    [[nodiscard]] Limits limits(unsigned refinements) const {
        return {n_, 1, 0, refinements, EstimateAdmission{1 << 20}};
    }
    RealBackend backend_;
    core::Index n_;
    core::ResourcePlan input_plan_{1 << 20}, workspace_plan_{1 << 20}, output_plan_{1 << 20};
    DenseReal matrix_;
    SparseReal sparse_;
    std::optional<DensePrepared> dense_;
    std::optional<SparsePrepared> prepared_sparse_;
    std::array<double, 3> scales_{1, 1, 1};
    std::array<UnitDimensions, 3> units_{};
    NoChecks<double> consumer_;
    std::uint64_t next_generation_ = 40;
};

using Complex = std::complex<double>;
class UpdatedSystem {
public:
    UpdatedSystem(core::Index order, std::span<const double> coefficients,
                  std::span<const Complex> left, std::span<const Complex> right)
        : n_(order) {
        auto a = DenseReal::copy({n_, n_, n_}, coefficients, input_plan_, operation);
        auto u = DenseComplex::copy({n_, 1, n_}, left, input_plan_, operation);
        auto v = DenseComplex::copy({n_, 1, n_}, right, input_plan_, operation);
        REQUIRE(a);
        REQUIRE(u);
        REQUIRE(v);
        a_ = std::move(*a);
        u_ = std::move(*u);
        v_ = std::move(*v);
        auto prepared = UpdatedPrepared::prepare(op(), scales(), limits(5), workspace_plan_);
        REQUIRE(prepared);
        prepared_.emplace(std::move(*prepared));
    }
    [[nodiscard]] PrecisionBuffer<Complex> rhs(std::span<const Complex> hi,
                                               std::span<const Complex> lo = {}) {
        auto h = DenseComplex::copy({n_, 1, n_}, hi, input_plan_, operation);
        REQUIRE(h);
        if (lo.empty()) {
            auto result =
                PrecisionBuffer<Complex>::binary64(std::move(*h), core::NumericalGeneration{30});
            REQUIRE(result);
            return std::move(*result);
        }
        auto l = DenseComplex::copy({n_, 1, n_}, lo, input_plan_, operation);
        REQUIRE(l);
        auto result = PrecisionBuffer<Complex>::two_word(std::move(*h), std::move(*l),
                                                         core::NumericalGeneration{30});
        REQUIRE(result);
        return std::move(*result);
    }
    [[nodiscard]] core::Result<Solution<Complex>>
    solve(const PrecisionBuffer<Complex>& rhs, double budget, unsigned refinements, Norm norm,
          double relative = 0, bool project = false, std::span<const ConsumerCheck> checks = {},
          const ConsumerEvaluator<Complex>* evaluator = nullptr) {
        if (!prepared_) {
            throw std::logic_error("Missing updated fixture preparation");
        }
        const Accuracy accuracy{norm,
                                budget,
                                relative,
                                std::max(budget, relative),
                                core::ObservationId{20},
                                core::ObservationId{21}};
        const auto generation = core::NumericalGeneration{next_generation_++};
        ReturnedPrecision returned = RetainedPrecision{generation};
        if (project) {
            returned =
                Binary64Projection{generation, core::NumericalGeneration{next_generation_++}};
        }
        const auto& consumer = evaluator != nullptr
                                   ? *evaluator
                                   : static_cast<const ConsumerEvaluator<Complex>&>(consumer_);
        return solve_linear(op(), Rhs{rhs, context(), core::RhsGeneration{31}}, scales(),
                            Precision<Complex>{returned, checks, consumer}, accuracy,
                            limits(refinements), *prepared_, output_plan_);
    }
    [[nodiscard]] const DenseReal& original() const noexcept {
        return a_;
    }
    [[nodiscard]] const DenseComplex& left() const noexcept {
        return u_;
    }
    [[nodiscard]] const DenseComplex& right() const noexcept {
        return v_;
    }

private:
    [[nodiscard]] ColumnUpdatedOperator op() const {
        return {a_, u_, v_, context()};
    }
    [[nodiscard]] Scales scales() const {
        const auto n = static_cast<std::size_t>(n_);
        return {std::span{scales_}.first(n), std::span{scales_}.first(n),
                std::span{units_}.first(n), std::span{units_}.first(n), context()};
    }
    [[nodiscard]] Limits limits(unsigned refinements) const {
        return {n_, 1, 1, refinements, EstimateAdmission{1 << 20}};
    }
    core::Index n_;
    core::ResourcePlan input_plan_{1 << 20}, workspace_plan_{1 << 20}, output_plan_{1 << 20};
    DenseReal a_;
    DenseComplex u_, v_;
    std::optional<UpdatedPrepared> prepared_;
    std::array<double, 3> scales_{1, 1, 1};
    std::array<UnitDimensions, 3> units_{};
    NoChecks<Complex> consumer_;
    std::uint64_t next_generation_ = 40;
};
[[nodiscard]] Exact complex_actual(const PrecisionBuffer<Complex>& buffer, std::size_t index,
                                   bool imaginary) {
    const auto part = [imaginary](Complex value) {
        return imaginary ? value.imag() : value.real();
    };
    if (const auto* binary = std::get_if<Binary64<Complex>>(&buffer.storage())) {
        return Exact{part(binary->values.values()[index])};
    }
    const auto* words = std::get_if<TwoWord<Complex>>(&buffer.storage());
    if (words == nullptr) {
        throw std::runtime_error("Missing complex retained words");
    }
    return Exact{part(words->hi.values()[index])} + Exact{part(words->lo.values()[index])};
}
[[nodiscard]] Exact actual(const PrecisionBuffer<double>& buffer, std::size_t index) {
    if (const auto* binary = std::get_if<Binary64<double>>(&buffer.storage())) {
        return Exact{binary->values.values()[index]};
    }
    const auto* words = std::get_if<TwoWord<double>>(&buffer.storage());
    if (words == nullptr) {
        throw std::runtime_error("Missing retained numerical words");
    }
    return Exact{words->hi.values()[index]} + Exact{words->lo.values()[index]};
}
void verify_real(const Solution<double>& result, double budget, Norm norm,
                 std::span<const double> coefficients, std::span<const double> right,
                 std::span<const double> truth) {
    Exact squared = 0;
    for (std::size_t i = 0; i < right.size(); ++i) {
        const Exact x = actual(result.values(), i);
        CHECK(abs(x - Exact{truth[i]}) < Exact{"1/1000000000000"});
        Exact residual = Exact{right[i]}, denominator = abs(Exact{right[i]});
        for (std::size_t j = 0; j < right.size(); ++j) {
            const Exact coefficient{coefficients[j * right.size() + i]};
            const Exact value = actual(result.values(), j);
            residual -= coefficient * value;
            denominator += abs(coefficient) * abs(value);
        }
        if (budget == 0) {
            CHECK(residual == 0);
        } else if (norm == Norm::Infinity) {
            CHECK(abs(residual) <= Exact{budget});
        } else {
            const Exact normalized = residual / Exact{budget};
            squared += normalized * normalized;
        }
        CHECK(abs(residual) <= Exact{budget} * denominator);
    }
    if (budget > 0 && norm == Norm::L2) {
        CHECK(squared <= 1);
    }
    CHECK(result.certificate().columns()[0].normalized_original_residual_upper <= 1);
    CHECK(result.certificate().columns()[0].componentwise_backward_error_upper <= budget);
}
} // namespace

TEST_CASE("Direct refinement reaches the fixed integer solution after unrepresentable metrics",
          "[verification][negative][N4-05][N4-06]") {
    // Independently fixed row products are 7+10+24, -2-10+24, 6+10+18.
    // det(A)=-30, so x=(1,2,3) is the unique exact solution.
    const double tiny = std::ldexp(1.0, -1063);
    for (const auto backend : {RealBackend::Dense, RealBackend::SparseLu}) {
        for (const auto norm : {Norm::Infinity, Norm::L2}) {
            RealSystem system{backend, 3, original_a};
            const auto rhs = system.rhs(original_b);
            const auto control = system.solve(rhs, 1.e-12, 0, norm);
            REQUIRE(control);
            CHECK(control->certificate().refinements() == 0);
            verify_real(*control, 1.e-12, norm, original_a, original_b, exact_x);
            const std::array snapshot{actual(control->values(), 0), actual(control->values(), 1),
                                      actual(control->values(), 2)};
            const auto no_refinement = system.solve(rhs, 0, 0, norm);
            // Different conforming factorizations may already return the exact solution.
            // Both outcomes are checked; no platform or fixture is skipped.
            if (no_refinement) {
                CHECK(no_refinement->certificate().refinements() == 0);
                verify_real(*no_refinement, 0, norm, original_a, original_b, exact_x);
            } else {
                CHECK(no_refinement.error().code == core::FailureCode::LinearUnconverged);
            }
            for (const double budget : {0.0, tiny}) {
                const auto refined = system.solve(rhs, budget, 5, norm);
                REQUIRE(refined);
                if (!no_refinement) {
                    CHECK(refined->certificate().refinements() > 0);
                }
                CHECK(refined->certificate().refinements() <= 5);
                verify_real(*refined, budget, norm, original_a, original_b, exact_x);
                CHECK(refined->certificate().context() == context());
                CHECK(refined->certificate().rhs_buffer_generation() == rhs.generation());
            }
            for (std::size_t i = 0; i < snapshot.size(); ++i) {
                CHECK(actual(control->values(), i) == snapshot[i]);
                CHECK(actual(rhs, i) == Exact{original_b[i]});
            }
            CHECK(std::ranges::equal(system.original_dense().values(), original_a));
            CHECK(std::ranges::equal(system.original_sparse().values(), original_a));
        }
    }
}

TEST_CASE("Direct refinement preserves genuine underflow-grid enclosure failure",
          "[negative][N4-05][N4-06]") {
    // fl(minsub/1.5)=minsub, but the exact product 1.5*minsub has a bit
    // below the binary64 subnormal grid. The EFT cannot return that enclosure.
    const std::array coefficient{1.5};
    const std::array right{std::numeric_limits<double>::denorm_min()};
    for (const auto backend : {RealBackend::Dense, RealBackend::SparseLu}) {
        RealSystem system{backend, 1, coefficient};
        const auto rhs = system.rhs(right);
        const auto rejected = system.solve(rhs, 0, 5, Norm::Infinity);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::PrecisionUnresolved);
        CHECK(actual(rhs, 0) == Exact{right[0]});
        CHECK(system.original_dense().values()[0] == 1.5);
        CHECK(system.original_sparse().values()[0] == 1.5);
    }
}

TEST_CASE("Sparse LLT refines an initially unmet exact SPD solution",
          "[verification][negative][N4-05][N4-06]") {
    // Sylvester minors 3 and 5 prove SPD. The exact row products are 3+2 and 1+4.
    const std::array coefficients{3., 1., 1., 2.};
    const std::array right{5., 5.};
    const std::array truth{1., 2.};
    for (const auto norm : {Norm::Infinity, Norm::L2}) {
        RealSystem system{RealBackend::SparseLlt, 2, coefficients};
        const auto rhs = system.rhs(right);
        const auto control = system.solve(rhs, 1.e-12, 0, norm);
        REQUIRE(control);
        CHECK(control->certificate().refinements() == 0);
        verify_real(*control, 1.e-12, norm, coefficients, right, truth);
        const std::array snapshot{actual(control->values(), 0), actual(control->values(), 1)};
        const auto capped = system.solve(rhs, 0, 0, norm);
        if (capped) {
            CHECK(capped->certificate().refinements() == 0);
            verify_real(*capped, 0, norm, coefficients, right, truth);
        } else {
            CHECK(capped.error().code == core::FailureCode::LinearUnconverged);
        }
        for (const double budget : {0.0, std::ldexp(1., -1063)}) {
            const auto refined = system.solve(rhs, budget, 5, norm);
            REQUIRE(refined);
            if (!capped) {
                CHECK(refined->certificate().refinements() > 0);
            }
            CHECK(refined->certificate().refinements() <= 5);
            verify_real(*refined, budget, norm, coefficients, right, truth);
        }
        for (std::size_t i = 0; i < right.size(); ++i) {
            CHECK(actual(control->values(), i) == snapshot[i]);
            CHECK(actual(rhs, i) == Exact{right[i]});
        }
        CHECK(std::ranges::equal(system.original_sparse().values(), coefficients));
    }
}

TEST_CASE("An unrepresentable first row does not hide a later raw enclosure failure",
          "[negative][N4-05][N4-06]") {
    // Row zero has a finite nonzero residual with zero allowance. Row one instead
    // fails the actual product enclosure: 1.5*minsub needs a bit below the grid.
    // With no IR allowed, swallowing that raw failure would give LinearUnconverged.
    const std::array coefficients{3., 0., 0., 1.5};
    const std::array right{1., std::numeric_limits<double>::denorm_min()};
    for (const auto backend : {RealBackend::Dense, RealBackend::SparseLu, RealBackend::SparseLlt}) {
        RealSystem system{backend, 2, coefficients};
        const auto rhs = system.rhs(right);
        const auto rejected = system.solve(rhs, 0, 0, Norm::Infinity);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::PrecisionUnresolved);
        CHECK(actual(rhs, 0) == 1);
        CHECK(actual(rhs, 1) == Exact{right[1]});
    }
}

TEST_CASE("Complex low-rank refinement reaches an exact original integer system",
          "[verification][negative][N4-05][N4-06]") {
    // Aeff=A0+i*e0*e0^T; det(Aeff)=-30-70i. This is V^T, not V^H.
    // x=(1+i)*(1,2,3), A0*x=(1+i)*(41,12,34); the update adds -1+i in row zero.
    const std::array left{Complex{1}, Complex{}, Complex{}};
    const std::array right{Complex{0, 1}, Complex{}, Complex{}};
    const std::array load{Complex{40, 42}, Complex{12, 12}, Complex{34, 34}};
    for (const auto norm : {Norm::Infinity, Norm::L2}) {
        UpdatedSystem system{3, original_a, left, right};
        const auto rhs = system.rhs(load);
        const auto control = system.solve(rhs, 1.e-12, 0, norm);
        REQUIRE(control);
        CHECK(control->certificate().refinements() == 0);
        const std::array snapshot{complex_actual(control->values(), 0, false),
                                  complex_actual(control->values(), 1, false),
                                  complex_actual(control->values(), 2, false),
                                  complex_actual(control->values(), 0, true),
                                  complex_actual(control->values(), 1, true),
                                  complex_actual(control->values(), 2, true)};
        const auto capped = system.solve(rhs, 0, 0, norm);
        if (capped) {
            CHECK(capped->certificate().refinements() == 0);
            for (std::size_t i = 0; i < 3; ++i) {
                CHECK(complex_actual(capped->values(), i, false) == Exact{exact_x[i]});
                CHECK(complex_actual(capped->values(), i, true) == Exact{exact_x[i]});
            }
        } else {
            CHECK(capped.error().code == core::FailureCode::LinearUnconverged);
        }
        for (const double budget : {0., std::ldexp(1., -1063)}) {
            const auto refined = system.solve(rhs, budget, 5, norm);
            REQUIRE(refined);
            if (!capped) {
                CHECK(refined->certificate().refinements() > 0);
            }
            CHECK(refined->certificate().refinements() <= 5);
            // Exact rational equalities verify both retained planes and the original
            // updated operator independently of the production residual evaluator.
            for (std::size_t i = 0; i < 3; ++i) {
                CHECK(complex_actual(refined->values(), i, false) == Exact{exact_x[i]});
                CHECK(complex_actual(refined->values(), i, true) == Exact{exact_x[i]});
                for (const bool imaginary : {false, true}) {
                    Exact residual{imaginary ? load[i].imag() : load[i].real()};
                    for (std::size_t j = 0; j < 3; ++j) {
                        residual -= Exact{original_a[j * 3 + i]} *
                                    complex_actual(refined->values(), j, imaginary);
                    }
                    if (i == 0) {
                        if (imaginary)
                            residual -= complex_actual(refined->values(), 0, false);
                        else
                            residual += complex_actual(refined->values(), 0, true);
                    }
                    CHECK(residual == 0);
                }
            }
            CHECK(refined->certificate().columns()[0].normalized_original_residual_upper <= 1);
            CHECK(refined->certificate().columns()[0].componentwise_backward_error_upper <= budget);
        }
        for (std::size_t i = 0; i < 3; ++i) {
            CHECK(complex_actual(control->values(), i, false) == snapshot[i]);
            CHECK(complex_actual(control->values(), i, true) == snapshot[i + 3]);
            CHECK(complex_actual(rhs, i, false) == Exact{load[i].real()});
            CHECK(complex_actual(rhs, i, true) == Exact{load[i].imag()});
        }
        CHECK(std::ranges::equal(system.original().values(), original_a));
        CHECK(std::ranges::equal(system.left().values(), left));
        CHECK(std::ranges::equal(system.right().values(), right));
    }
}

TEST_CASE("A direct retained low word cannot be silently projected under zero allowances",
          "[verification][negative][N4-05][N4-06]") {
    const double delta = std::ldexp(1., -80);
    const std::array identity{1.}, hi{1.}, lo{delta};
    const Exact expected = Exact{1} + Exact{delta};
    for (const auto backend : {RealBackend::Dense, RealBackend::SparseLu, RealBackend::SparseLlt}) {
        RealSystem system{backend, 1, identity};
        const auto rhs = system.rhs_words(hi, lo);
        const auto retained = system.solve(rhs, 0, 0, Norm::Infinity);
        REQUIRE(retained);
        CHECK(retained->values().precision() == core::Precision::TwoWord);
        CHECK(retained->certificate().refinements() == 0);
        CHECK(actual(retained->values(), 0) == expected);
        const auto projected = system.solve(rhs, 0, 0, Norm::Infinity, 0, true);
        REQUIRE_FALSE(projected);
        CHECK(projected.error().code == core::FailureCode::PrecisionUnresolved);
        CHECK(actual(retained->values(), 0) == expected);
        CHECK(actual(rhs, 0) == expected);
    }
    const std::array zero{Complex{}}, complex_hi{Complex{1, 1}}, complex_lo{Complex{delta, delta}};
    UpdatedSystem system{1, identity, zero, zero};
    const auto rhs = system.rhs(complex_hi, complex_lo);
    const auto retained = system.solve(rhs, 0, 0, Norm::Infinity);
    REQUIRE(retained);
    CHECK(retained->values().precision() == core::Precision::TwoWord);
    CHECK(retained->certificate().refinements() == 0);
    const auto projected = system.solve(rhs, 0, 0, Norm::Infinity, 0, true);
    REQUIRE_FALSE(projected);
    CHECK(projected.error().code == core::FailureCode::PrecisionUnresolved);
    for (const bool imaginary : {false, true}) {
        CHECK(complex_actual(retained->values(), 0, imaginary) == expected);
        CHECK(complex_actual(rhs, 0, imaginary) == expected);
    }
}

TEST_CASE("Direct relative-only qualification keeps zero right-hand-side rows exact",
          "[verification][negative][N4-05][N4-06]") {
    // All products below are exact integers, fixed independently of the factorization.
    // atol=0 gives zero allowance at the zero RHS row even though rtol is positive.
    const std::array real_right{0., 21., -3.}, real_truth{1., -3., 1.};
    for (const auto backend : {RealBackend::Dense, RealBackend::SparseLu, RealBackend::SparseLlt}) {
        const bool llt = backend == RealBackend::SparseLlt;
        const std::array spd{3., 1., 1., 2.};
        const std::array spd_right{0., -5.}, spd_truth{1., -3.};
        const std::span<const double> coefficients =
            llt ? std::span<const double>{spd} : std::span<const double>{original_a};
        const std::span<const double> right =
            llt ? std::span<const double>{spd_right} : std::span<const double>{real_right};
        const std::span<const double> truth =
            llt ? std::span<const double>{spd_truth} : std::span<const double>{real_truth};
        for (const auto norm : {Norm::Infinity, Norm::L2}) {
            RealSystem system{backend, static_cast<core::Index>(right.size()), coefficients};
            const auto rhs = system.rhs(right);
            const auto result = system.solve(rhs, 0, 5, norm, 1.e-12);
            REQUIRE(result);
            Exact squared = 0;
            for (std::size_t i = 0; i < right.size(); ++i) {
                CHECK(abs(actual(result->values(), i) - Exact{truth[i]}) <
                      Exact{"1/1000000000000"});
                Exact residual{right[i]}, denominator = abs(Exact{right[i]});
                for (std::size_t j = 0; j < right.size(); ++j) {
                    const Exact coefficient{coefficients[j * right.size() + i]};
                    residual -= coefficient * actual(result->values(), j);
                    denominator += abs(coefficient) * abs(actual(result->values(), j));
                }
                const Exact allowance = Exact{1.e-12} * abs(Exact{right[i]});
                if (right[i] == 0)
                    CHECK(residual == 0);
                else {
                    const Exact normalized = residual / allowance;
                    if (norm == Norm::Infinity)
                        CHECK(abs(normalized) <= 1);
                    else
                        squared += normalized * normalized;
                }
                CHECK(abs(residual) <= Exact{1.e-12} * denominator);
                CHECK(actual(rhs, i) == Exact{right[i]});
            }
            CHECK(squared <= 1);
            CHECK(result->certificate().columns()[0].normalized_original_residual_upper <= 1);
        }
    }
    // x=(1+i)*(1,-12,9); base products are (1+i)*(19,130,0),
    // and i*x0=-1+i is added only to row zero.
    const std::array left{Complex{1}, Complex{}, Complex{}},
        right{Complex{0, 1}, Complex{}, Complex{}};
    const std::array load{Complex{18, 20}, Complex{130, 130}, Complex{}};
    const std::array truth{1., -12., 9.};
    for (const auto norm : {Norm::Infinity, Norm::L2}) {
        UpdatedSystem system{3, original_a, left, right};
        const auto rhs = system.rhs(load);
        const auto result = system.solve(rhs, 0, 5, norm, 1.e-12);
        REQUIRE(result);
        Exact squared = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            Exact residual_squared = 0;
            for (const bool imaginary : {false, true}) {
                CHECK(abs(complex_actual(result->values(), i, imaginary) - Exact{truth[i]}) <
                      Exact{"1/1000000000000"});
                Exact residual{imaginary ? load[i].imag() : load[i].real()};
                for (std::size_t j = 0; j < 3; ++j)
                    residual -= Exact{original_a[j * 3 + i]} *
                                complex_actual(result->values(), j, imaginary);
                if (i == 0) {
                    if (imaginary)
                        residual -= complex_actual(result->values(), 0, false);
                    else
                        residual += complex_actual(result->values(), 0, true);
                }
                residual_squared += residual * residual;
                CHECK(complex_actual(rhs, i, imaginary) ==
                      Exact{imaginary ? load[i].imag() : load[i].real()});
            }
            const Exact allowance_squared = Exact{1.e-12} * Exact{1.e-12} *
                                            (Exact{load[i].real()} * Exact{load[i].real()} +
                                             Exact{load[i].imag()} * Exact{load[i].imag()});
            if (i == 2)
                CHECK(residual_squared == 0);
            else if (norm == Norm::Infinity)
                CHECK(residual_squared <= allowance_squared);
            else
                squared += residual_squared / allowance_squared;
        }
        CHECK(squared <= 1);
        CHECK(result->certificate().columns()[0].normalized_original_residual_upper <= 1);
        CHECK(result->certificate().columns()[0].componentwise_backward_error_upper <= 1.e-12);
    }
}
namespace {
template <StorageScalar Scalar>
class TwoChecks final : public ConsumerEvaluator<Scalar> {
public:
    mutable std::array<std::optional<core::CertificateContext>, 2> seen{};
    mutable unsigned calls = 0;
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<Scalar>&, const PrecisionBuffer<Scalar>& candidate,
           const Context& system, const ConsumerCheck& check,
           core::CertificateContext expected) const override {
        CHECK(system == context());
        CHECK(expected.numerical_generation == candidate.generation());
        CHECK(expected.representation == candidate.precision());
        const bool first = check.observation == core::ObservationId{80};
        seen[first ? 0 : 1] = expected;
        ++calls;
        // A protocol-only gate: the first necessary consumer succeeds and the
        // second reports an exceeded bound. This is not a physical error oracle.
        return core::ErrorCertificate::assess(expected, check.budget, first ? 0. : 1.);
    }
};
template <StorageScalar Scalar>
Exact component(const PrecisionBuffer<Scalar>& values, std::size_t index, bool imaginary) {
    if constexpr (std::is_same_v<Scalar, double>)
        return imaginary ? Exact{0} : actual(values, index);
    else
        return complex_actual(values, index, imaginary);
}
template <class System, StorageScalar Scalar>
void verify_consumer_gate(System& system, const PrecisionBuffer<Scalar>& rhs) {
    const auto retained = system.solve(rhs, 0, 5, Norm::Infinity);
    REQUIRE(retained);
    const auto generation = retained->values().generation();
    std::array<Exact, 6> snapshot{};
    const auto n = static_cast<std::size_t>(rhs.shape().rows);
    for (std::size_t i = 0; i < n; ++i) {
        snapshot[i] = component(retained->values(), i, false);
        snapshot[i + 3] = component(retained->values(), i, true);
    }
    const std::array checks{
        ConsumerCheck{CheckKind::Observation, core::ObservationId{80}, core::ScopeId{81}, 1.e-8},
        ConsumerCheck{CheckKind::Reaction, core::ObservationId{82}, core::ScopeId{83}, 1.e-8}};
    const TwoChecks<Scalar> evaluator;
    const auto refused = system.solve(rhs, 0, 5, Norm::Infinity, 0, false, checks, &evaluator);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().code == core::FailureCode::PrecisionUnresolved);
    CHECK(refused.error().details.problem == core::ProblemId{core::ObservationId{82}});
    CHECK(evaluator.calls == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        if (!evaluator.seen[i]) {
            FAIL_CHECK("Every necessary consumer must evaluate the current candidate");
            return;
        }
        const auto& seen = *evaluator.seen[i];
        CHECK(seen.observation == checks[i].observation);
        CHECK(seen.scope == checks[i].scope);
        CHECK(seen.numerical_generation != generation);
        CHECK(seen.numerical_generation == refused.error().details.numerical_generation);
    }
    CHECK(retained->values().generation() == generation);
    CHECK(retained->certificate().consumer_count() == 0);
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(component(retained->values(), i, false) == snapshot[i]);
        CHECK(component(retained->values(), i, true) == snapshot[i + 3]);
    }
}
} // namespace
TEST_CASE("A newly refined direct candidate must pass every necessary consumer",
          "[negative][N4-05][N4-06]") {
    for (const auto backend : {RealBackend::Dense, RealBackend::SparseLu, RealBackend::SparseLlt}) {
        const bool llt = backend == RealBackend::SparseLlt;
        const std::array spd{3., 1., 1., 2.};
        const std::array spd_right{5., 5.};
        const std::span<const double> coefficients =
            llt ? std::span<const double>{spd} : std::span<const double>{original_a};
        const std::span<const double> right =
            llt ? std::span<const double>{spd_right} : std::span<const double>{original_b};
        RealSystem system{backend, static_cast<core::Index>(right.size()), coefficients};
        const auto rhs = system.rhs(right);
        verify_consumer_gate(system, rhs);
        for (std::size_t i = 0; i < right.size(); ++i)
            CHECK(actual(rhs, i) == Exact{right[i]});
    }
    const std::array left{Complex{1}, Complex{}, Complex{}},
        right{Complex{0, 1}, Complex{}, Complex{}};
    const std::array load{Complex{40, 42}, Complex{12, 12}, Complex{34, 34}};
    UpdatedSystem system{3, original_a, left, right};
    const auto rhs = system.rhs(load);
    verify_consumer_gate(system, rhs);
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(complex_actual(rhs, i, false) == Exact{load[i].real()});
        CHECK(complex_actual(rhs, i, true) == Exact{load[i].imag()});
    }
}
