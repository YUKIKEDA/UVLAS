#include "exact_polynomial_detail.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <uvlas/structure/exact_polynomial.hpp>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using Json = nlohmann::json;
constexpr auto budget = std::uint64_t{64} * 1024 * 1024;
ExactLimits limits() {
    return {10000,
            1000,
            10000,
            100,
            200,
            4096,
            1000000,
            ExactEstimateAdmission{budget},
            ExactTimingMode::Disabled};
}
template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw std::runtime_error("Unexpected exact preparation failure " +
                                 std::to_string(static_cast<unsigned>(result.error().code)));
    return std::move(*result);
}
ExactRational rational(std::string_view numerator, std::string_view denominator,
                       core::ResourcePlan& plan) {
    return take(ExactRational::from_integer_pair(numerator, denominator, limits(), plan));
}
std::pair<std::string, std::string> canonical(const ExactRational& value) {
    const auto& backend = detail::ExactFactory::value(value).backend();
    return {backend.num().str(0, std::ios_base::dec), backend.denom().str(0, std::ios_base::dec)};
}
void fraction_is(const ExactRational& value, std::string_view numerator,
                 std::string_view denominator) {
    const auto actual = canonical(value);
    CHECK(actual.first == numerator);
    CHECK(actual.second == denominator);
}
Json fixture() {
    std::ifstream input(UVLAS_EXACT_POLYNOMIAL_FIXTURE);
    if (!input)
        throw std::runtime_error("Missing frozen exact fixture");
    Json result;
    input >> result;
    if (result.at("schema") != "uvlas.exact-polynomial-reference/1" ||
        result.at("id") != "exact-polynomial-canonical-v1")
        throw std::runtime_error("Unsupported frozen exact fixture");
    return result;
}
ExactChartPolynomial polynomial(std::span<const std::array<int, 4>> terms,
                                core::ResourcePlan& plan) {
    std::vector<ExactRational> coefficients;
    coefficients.reserve(terms.size());
    for (const auto& term : terms)
        coefficients.push_back(rational(std::to_string(term[2]), std::to_string(term[3]), plan));
    std::vector<ExactInputTerm> input;
    input.reserve(terms.size());
    for (std::size_t i = 0; i < terms.size(); ++i)
        input.push_back({{terms[i][0], terms[i][1]}, coefficients[i]});
    return take(ExactChartPolynomial::make(input, limits(), plan));
}
void polynomial_is(const ExactChartPolynomial& value, const Json& expected) {
    const auto actual = value.terms();
    REQUIRE(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) {
        CHECK(actual[i].powers.first == expected[i][0].get<core::Index>());
        CHECK(actual[i].powers.second == expected[i][1].get<core::Index>());
        fraction_is(actual[i].coefficient, expected[i][2].get<std::string>(),
                    expected[i][3].get<std::string>());
    }
}
ExactChartPolynomial F(core::ResourcePlan& plan) {
    constexpr std::array terms{std::array{2, 1, 1, 3}, std::array{0, 0, 5, 6},
                               std::array{2, 1, 1, 6}, std::array{1, 2, -3, 4},
                               std::array{4, 4, 0, 1}};
    return polynomial(terms, plan);
}
template <class A, class B>
concept CanAdd = requires(const A& a, const B& b, const ExactLimits& l, core::ResourcePlan& p) {
    a.add(b, l, p);
};
template <class T>
concept RvalueTerms = requires(const T&& value) { std::move(value).terms(); };
static_assert(!CanAdd<ExactChartPolynomial, ExactCellPolynomial>);
static_assert(!RvalueTerms<ExactChartPolynomial>);
static_assert(!std::is_copy_constructible_v<ExactRational>);
static_assert(!std::is_copy_constructible_v<ExactChartPolynomial>);
} // namespace

TEST_CASE("Exact fractions match independent canonical numerator and denominator literals",
          "[SH-D01][SH-D03][N4-02][N4-08][verification]") {
    const auto reference = fixture();
    core::ResourcePlan plan{budget};
    for (const auto& row : reference["rational"]) {
        CAPTURE(row);
        auto a = rational(row["a"][0].get<std::string>(), row["a"][1].get<std::string>(), plan);
        const auto op = row["op"].get<std::string>();
        auto result = take(a.clone(limits(), plan));
        if (op != "normalize") {
            const auto b =
                rational(row["b"][0].get<std::string>(), row["b"][1].get<std::string>(), plan);
            if (op == "add")
                result = take(a.add(b, limits(), plan));
            else if (op == "subtract")
                result = take(a.subtract(b, limits(), plan));
            else if (op == "multiply")
                result = take(a.multiply(b, limits(), plan));
            else
                result = take(a.divide(b, limits(), plan));
        }
        fraction_is(result, row["expected"][0].get<std::string>(),
                    row["expected"][1].get<std::string>());
    }
    for (const auto& row : reference["decimal"]) {
        const auto result =
            take(ExactRational::from_decimal(row["input"].get<std::string>(), limits(), plan));
        fraction_is(result, row["expected"][0].get<std::string>(),
                    row["expected"][1].get<std::string>());
    }
    for (const auto& row : reference["binary64"]) {
        const auto text = row["hex"].get<std::string>();
        const auto input = std::strtod(text.c_str(), nullptr);
        const auto result = take(ExactRational::from_binary64(input, limits(), plan));
        fraction_is(result, row["expected"][0].get<std::string>(),
                    row["expected"][1].get<std::string>());
    }
    const auto negative_zero = take(ExactRational::from_binary64(-0., limits(), plan));
    fraction_is(negative_zero, "0", "1");
}
TEST_CASE("Exact polynomial operations match hand-expanded canonical term arrays",
          "[SH-D01][SH-D03][N4-02][verification]") {
    const auto reference = fixture();
    core::ResourcePlan plan{budget};
    const auto f = F(plan);
    polynomial_is(f, reference["F"]);
    CHECK(f.degree() == 3);
    constexpr std::array terms{std::array{1, 0, 1, 1}, std::array{0, 1, -2, 1}};
    const auto g = polynomial(terms, plan);
    const auto fg = take(f.multiply(g, limits(), plan));
    polynomial_is(fg, reference["FG"]);
    const auto fs = take(f.differentiate(ExactAxis::First, 1, limits(), plan));
    const auto ft = take(f.differentiate(ExactAxis::Second, 1, limits(), plan));
    polynomial_is(fs, reference["Fs"]);
    polynomial_is(ft, reference["Ft"]);
    const auto fst = take(fs.differentiate(ExactAxis::Second, 1, limits(), plan));
    polynomial_is(fst, reference["Fst"]);
    const auto x = rational("2", "3", plan), y = rational("-1", "2", plan);
    const auto evaluated = take(f.evaluate(x, y, limits(), plan));
    fraction_is(evaluated, "43", "72");
    const auto lower = rational("1", "3", plan);
    const auto integrated = take(f.integrate(ExactAxis::First, lower, limits(), plan));
    polynomial_is(integrated, reference["Fis"]);
    const auto at_lower = take(integrated.evaluate(lower, y, limits(), plan));
    fraction_is(at_lower, "0", "1");
    constexpr std::array x_terms{std::array{1, 0, 1, 1}, std::array{0, 0, 1, 1}};
    constexpr std::array<std::array<int, 4>, 1> y_terms{{{0, 1, 2, 1}}};
    const auto sx = polynomial(x_terms, plan), sy = polynomial(y_terms, plan);
    const auto composed = take(f.substitute(sx, sy, limits(), plan));
    polynomial_is(composed, reference["Fsub"]);
    const auto cancelled = take(f.subtract(f, limits(), plan));
    CHECK(cancelled.terms().empty());
    CHECK_FALSE(cancelled.degree());
    const auto zero_order = take(f.differentiate(ExactAxis::First, 0, limits(), plan));
    polynomial_is(zero_order, reference["F"]);
    const auto high_order = take(f.differentiate(ExactAxis::Second, 99, limits(), plan));
    CHECK(high_order.terms().empty());
    const auto doubled = take(f.add(f, limits(), plan));
    const auto half = rational("1", "2", plan);
    const auto recovered = take(doubled.scale(half, limits(), plan));
    polynomial_is(recovered, reference["F"]);
}
TEST_CASE("Exact preparation rejects invalid values and insufficient independent budgets",
          "[SH-D03][N4-02][N4-08][negative]") {
    auto l = limits();
    core::ResourcePlan input{budget};
    const auto one = rational("1", "1", input);
    const auto zero = rational("0", "1", input);
    for (const auto bad : {"", "+", "1.0", " 1", "1e2", "0x10", "NaN"}) {
        core::ResourcePlan p{budget};
        const auto r = ExactRational::from_integer_pair(bad, "1", l, p);
        REQUIRE_FALSE(r);
        CHECK(r.error().code == core::FailureCode::InvalidInput);
        CHECK(p.allocated_bytes() == 0);
    }
    for (const auto bad : {"", ".", "+", "1e", "1e+", "1 2", "1.2.3", "nan", "Infinity"}) {
        core::ResourcePlan p{budget};
        const auto r = ExactRational::from_decimal(bad, l, p);
        REQUIRE_FALSE(r);
        CHECK(r.error().code == core::FailureCode::InvalidInput);
    }
    for (const auto bad :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()}) {
        core::ResourcePlan p{budget};
        const auto r = ExactRational::from_binary64(bad, l, p);
        REQUIRE_FALSE(r);
        CHECK(r.error().code == core::FailureCode::Nonfinite);
    }
    core::ResourcePlan none{0};
    const auto division = one.divide(zero, l, none);
    REQUIRE_FALSE(division);
    CHECK(division.error().code == core::FailureCode::InvalidInput);
    auto strict = l;
    strict.external_memory = ExactRequireProvenBound{budget};
    const auto unproved = one.clone(strict, none);
    REQUIRE_FALSE(unproved);
    CHECK(unproved.error().code == core::FailureCode::Unsupported);
    auto no_external = l;
    no_external.external_memory = ExactEstimateAdmission{0};
    core::ResourcePlan external_output{budget};
    const auto no_memory = one.add(one, no_external, external_output);
    REQUIRE_FALSE(no_memory);
    CHECK(no_memory.error().code == core::FailureCode::ResourceExceeded);
    CHECK(external_output.allocated_bytes() == 0);
    CHECK(none.allocated_bytes() == 0);
    const auto no_payload = one.clone(l, none);
    REQUIRE_FALSE(no_payload);
    CHECK(no_payload.error().code == core::FailureCode::ResourceExceeded);
    const auto probe = take(one.add(one, l, input));
    auto exact_cap = l;
    exact_cap.external_memory =
        ExactEstimateAdmission{probe.preparation().estimated_external_payload_bytes};
    core::ResourcePlan p{budget};
    CHECK(one.add(one, exact_cap, p).has_value());
    exact_cap.external_memory =
        ExactEstimateAdmission{probe.preparation().estimated_external_payload_bytes - 1};
    CHECK_FALSE(one.add(one, exact_cap, p));
    fraction_is(one, "1", "1");
    fraction_is(zero, "0", "1");
}
TEST_CASE("Exact values preserve ownership frame boundaries and immutable prior results",
          "[SH-D03][N4-02][N4-08][property][negative]") {
    auto kept = []() {
        core::ResourcePlan p{budget};
        return F(p);
    }();
    const auto expected = fixture();
    polynomial_is(kept, expected["F"]);
    auto moved = std::move(kept);
    // NOLINTNEXTLINE(bugprone-use-after-move): the public moved-from contract is canonical zero.
    CHECK(kept.terms().empty());
    polynomial_is(moved, expected["F"]);
    core::ResourcePlan p{budget};
    const auto one = rational("1", "1", p);
    const std::array bad{ExactInputTerm{{-1, 0}, one}};
    const auto rejected = ExactChartPolynomial::make(bad, limits(), p);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    auto bounded = limits();
    bounded.maximum_raw_terms = 1;
    const std::array duplicate{ExactInputTerm{{0, 0}, one}, ExactInputTerm{{0, 0}, one}};
    const auto raw_limit = ExactChartPolynomial::make(duplicate, bounded, p);
    REQUIRE_FALSE(raw_limit);
    CHECK(raw_limit.error().code == core::FailureCode::ResourceExceeded);
    bounded = limits();
    bounded.maximum_work_bits = 1;
    CHECK_FALSE(one.add(one, bounded, p));
    bounded = limits();
    bounded.maximum_value_bits = 1;
    const auto value_limit = ExactRational::from_integer_pair("3", "1", bounded, p);
    REQUIRE_FALSE(value_limit);
    CHECK(value_limit.error().code == core::FailureCode::ResourceExceeded);
    const auto unknown = moved.differentiate(static_cast<ExactAxis>(255), 1, limits(), p);
    REQUIRE_FALSE(unknown);
    CHECK(unknown.error().code == core::FailureCode::InvalidInput);
    const auto negative = moved.differentiate(ExactAxis::First, -1, limits(), p);
    REQUIRE_FALSE(negative);
    CHECK(negative.error().code == core::FailureCode::InvalidInput);
    auto measured = limits();
    measured.timing = ExactTimingMode::Enabled;
    const auto timed = take(moved.clone(measured, p));
    const auto ordinary = take(moved.clone(limits(), p));
    const auto& timing = timed.preparation().timing;
    if (!timing) {
        FAIL_CHECK("Enabled preparation did not retain its timing");
        return;
    }
    CHECK(timing->calls == 1);
    CHECK_FALSE(ordinary.preparation().timing);
    polynomial_is(timed, expected["F"]);
    polynomial_is(ordinary, expected["F"]);
    std::array<bool, 2> results{};
    std::array<std::thread, 2> threads;
    for (std::size_t i = 0; i < 2; ++i)
        threads[i] = std::thread([&, i]() {
            core::ResourcePlan local{budget};
            const auto v = moved.clone(limits(), local);
            results[i] = v.has_value() && v->terms().size() == 3;
        });
    for (auto& thread : threads)
        thread.join();
    CHECK(results[0]);
    CHECK(results[1]);
    polynomial_is(moved, expected["F"]);
}

TEST_CASE("Exact substitution validates resulting degrees before allocating its payload",
          "[SH-D01][SH-D03][N4-02][N4-08][negative]") {
    core::ResourcePlan inputs{budget};
    constexpr std::array<std::array<int, 4>, 1> square{{{2, 0, 1, 1}}};
    constexpr std::array<std::array<int, 4>, 1> second{{{0, 1, 1, 1}}};
    const auto f = polynomial(square, inputs);
    const auto x = polynomial(square, inputs);
    const auto y = polynomial(second, inputs);
    auto l = limits();
    l.maximum_axis_degree = 2;
    l.maximum_total_degree = 2;
    core::ResourcePlan output{budget};
    const auto invalid = f.substitute(x, y, l, output);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::ResourceExceeded);
    CHECK(output.allocated_bytes() == 0);
    const ExactChartPolynomial zero;
    const auto annihilated = take(f.substitute(zero, y, l, output));
    CHECK(annihilated.terms().empty());
    CHECK(output.allocated_bytes() == 0);
    const auto one = rational("1", "1", inputs);
    const std::array large_input{ExactInputTerm{{std::numeric_limits<core::Index>::max(), 1}, one}};
    l.maximum_axis_degree = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    l.maximum_total_degree = l.maximum_axis_degree;
    const auto sum_overflow = ExactChartPolynomial::make(large_input, l, output);
    REQUIRE_FALSE(sum_overflow);
    CHECK(sum_overflow.error().code == core::FailureCode::ResourceExceeded);
    CHECK(output.allocated_bytes() == 0);
}

TEST_CASE("Exact work caps include backend word rounding and guards",
          "[SH-D03][N4-02][N4-08][negative][property]") {
    core::ResourcePlan inputs{budget};
    const auto one = rational("1", "1", inputs);
    const auto probe = take(one.add(one, limits(), inputs));
    const auto guarded = probe.preparation().planned_work_bits;
    CHECK(guarded > probe.preparation().algebraic_work_bits);
    auto l = limits();
    l.maximum_work_bits = guarded;
    core::ResourcePlan accepted{budget};
    CHECK(one.add(one, l, accepted).has_value());
    --l.maximum_work_bits;
    core::ResourcePlan rejected{budget};
    const auto insufficient = one.add(one, l, rejected);
    REQUIRE_FALSE(insufficient);
    CHECK(insufficient.error().code == core::FailureCode::ResourceExceeded);
    CHECK(rejected.allocated_bytes() == 0);
    const ExactChartPolynomial empty;
    l = limits();
    l.maximum_work_bits = 0;
    core::ResourcePlan none{0};
    const auto zero = take(empty.differentiate(ExactAxis::First, 0, l, none));
    CHECK(zero.terms().empty());
    CHECK(zero.preparation().planned_work_bits == 0);
    CHECK(zero.preparation().estimated_external_payload_bytes == 0);
    const auto invalid = empty.differentiate(static_cast<ExactAxis>(255), 0, l, none);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
}

TEST_CASE("Exact ordering and coordinate conversion retain literal algebra",
          "[SH-D01][SH-D03][N4-02][verification]") {
    core::ResourcePlan plan{budget};
    const auto a = rational("-7", "3", plan), b = rational("2", "5", plan);
    CHECK(take(a.compare(b, limits(), plan)) == std::strong_ordering::less);
    CHECK(take(b.compare(a, limits(), plan)) == std::strong_ordering::greater);
    const auto same = rational("-14", "6", plan);
    CHECK(take(a.compare(same, limits(), plan)) == std::strong_ordering::equal);
    CHECK(a == same);
    const auto decimal = take(ExactRational::from_decimal("0.1", limits(), plan));
    const auto binary = take(ExactRational::from_binary64(0.1, limits(), plan));
    CHECK(take(decimal.compare(binary, limits(), plan)) == std::strong_ordering::less);
    const auto zero = take(a.subtract(same, limits(), plan));
    CHECK(zero.sign() == 0);
    CHECK(zero.numerator_bits() == 0);
    CHECK(zero.denominator_bits() == 1);
    CHECK(a.sign() == -1);
    const auto one = rational("1", "1", plan), two = rational("2", "1", plan);
    const std::array u_terms{ExactInputTerm{{1, 0}, one}};
    const std::array v_terms{ExactInputTerm{{0, 1}, two}};
    const auto u = take(ExactCellPolynomial::make(u_terms, limits(), plan));
    const auto v = take(ExactCellPolynomial::make(v_terms, limits(), plan));
    const auto source = F(plan);
    const auto cell = take(source.substitute(u, v, limits(), plan));
    // F(u,2v) = 5/6 + u^2 v - 3 u v^2. Literal values are independent.
    REQUIRE(cell.terms().size() == 3);
    const auto at = take(cell.evaluate(two, one, limits(), plan));
    fraction_is(at, "-7", "6");
    const auto lower = rational("-1", "1", plan);
    const auto integrated = take(cell.integrate(ExactAxis::Second, lower, limits(), plan));
    const auto at_lower = take(integrated.evaluate(two, lower, limits(), plan));
    fraction_is(at_lower, "0", "1");
    const auto at_one = take(integrated.evaluate(two, one, limits(), plan));
    fraction_is(at_one, "-7", "3");
    const auto derivative = take(integrated.differentiate(ExactAxis::Second, 1, limits(), plan));
    const auto recovered = take(derivative.evaluate(two, one, limits(), plan));
    fraction_is(recovered, "-7", "6");
}
