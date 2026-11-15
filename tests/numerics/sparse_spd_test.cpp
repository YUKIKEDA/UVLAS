#include "sparse_spd.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <span>
#include <uvlas/numerics/precision.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using Big = boost::multiprecision::cpp_bin_float_100;
constexpr auto operation = core::Operation::SolveLinear;

SparseReal sparse(core::Index n, std::initializer_list<Triplet> entries, core::ResourcePlan& plan) {
    auto result =
        SparseReal::from_triplets(n, n, {entries.begin(), entries.size()}, plan, operation);
    if (!result) {
        FAIL("Cannot construct the fixed sparse reference");
    }
    return std::move(result).value();
}

Big reference_shift(std::span<const unsigned> counts, std::span<const double> diagonal) {
    const Big unit = ldexp(Big{1}, -53);
    const Big eta = ldexp(Big{1}, -1074);
    Big total = 0;
    double largest = 0;
    for (std::size_t i = 0; i < counts.size(); ++i) {
        const Big ku = counts[i] * unit;
        total += ku / (1 - 2 * ku) * diagonal[i];
        largest = std::max(largest, diagonal[i]);
    }
    const Big n = counts.size();
    return total + 3 * n * (2 * n + largest) * eta;
}

void encloses(const SparseReal& matrix, const Big& reference) {
    auto result = linear::detail::rump_sparse_shift(matrix);
    if (!result) {
        FAIL("Expected finite Rump shift");
        return;
    }
    CHECK(Big{*result} >= reference);
    // This small fixture has at most a few dozen outward operations. This fixed
    // margin checks useful tightness independently, not a solver acceptance tolerance.
    CHECK(Big{*result} <= reference * (1 + ldexp(Big{1}, -45)));
}
} // namespace

TEST_CASE("Rump sparse bound encloses a non-diagonally-dominant SPD reference", "[N4-05]") {
    core::ResourcePlan plan{10000};
    auto matrix = sparse(
        3, {{0, 0, 2}, {1, 0, 4}, {2, 0, -6}, {0, 1, 4}, {1, 1, 11}, {0, 2, -6}, {2, 2, 71}}, plan);
    const std::array<unsigned, 3> counts{2, 3, 4};
    const std::array diagonal{2.0, 11.0, 71.0};
    encloses(matrix, reference_shift(counts, diagonal));
}

TEST_CASE("Rump counts the permuted envelope rather than stored entries", "[N4-05]") {
    core::ResourcePlan plan{10000};
    auto matrix =
        sparse(4, {{0, 0, 4}, {3, 0, 1}, {1, 1, 4}, {2, 2, 4}, {0, 3, 1}, {3, 3, 4}}, plan);
    const std::array<unsigned, 4> counts{2, 2, 2, 5};
    const std::array diagonal{4.0, 4.0, 4.0, 4.0};
    encloses(matrix, reference_shift(counts, diagonal));
}

TEST_CASE("Rump underflow term does not form an overflowing diagonal sum", "[N4-05]") {
    core::ResourcePlan plan{10000};
    const double largest = std::numeric_limits<double>::max();
    auto matrix = sparse(1, {{0, 0, largest}}, plan);
    const std::array<unsigned, 1> counts{2};
    const std::array diagonal{largest};
    encloses(matrix, reference_shift(counts, diagonal));
}

TEST_CASE("Rump rejects nonpositive diagonals and unsupported rounding", "[N4-05]") {
    core::ResourcePlan plan{10000};
    auto matrix = sparse(2, {{0, 0, 1}, {1, 1, -1}}, plan);
    auto result = linear::detail::rump_sparse_shift(matrix);
    REQUIRE_FALSE(result);
    CHECK(result.error().code == core::FailureCode::RankUnresolved);
    auto empty = SparseReal::copy_csc(0, 0, std::array<core::Index, 1>{0}, {}, {}, plan, operation);
    if (!empty) {
        FAIL("Empty reference");
        return;
    }
    auto zero = linear::detail::rump_sparse_shift(*empty);
    if (!zero) {
        FAIL("Empty bound");
        return;
    }
    CHECK(*zero == 0);
    auto positive = sparse(1, {{0, 0, 1}}, plan);
    const int saved = std::fegetround();
    REQUIRE(std::fesetround(FE_DOWNWARD) == 0);
    auto invalid = linear::detail::rump_sparse_shift(positive);
    const int restored = std::fesetround(saved);
    CHECK(restored == 0);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::PrecisionUnsupported);
}

TEST_CASE("CSC symmetry compares absent zeros without averaging", "[N4-05]") {
    core::ResourcePlan plan{10000};
    auto symmetric = sparse(2, {{0, 0, 2}, {1, 0, 0}, {1, 1, 3}}, plan);
    auto exact = linear::detail::sparse_exact_symmetry(symmetric);
    if (!exact) {
        FAIL("Symmetry check failed");
        return;
    }
    CHECK(*exact);
    auto asymmetric = sparse(2, {{0, 0, 2}, {1, 0, 1}, {0, 1, 2}, {1, 1, 3}}, plan);
    auto mismatch = linear::detail::sparse_exact_symmetry(asymmetric);
    if (!mismatch) {
        FAIL("Symmetry check failed");
        return;
    }
    CHECK_FALSE(*mismatch);
}

TEST_CASE("SPD congruence keeps every coefficient exact or skips scaling", "[N4-05]") {
    core::ResourcePlan plan{10000};
    const double tiny = std::ldexp(1.0, -500);
    const double coupling = std::ldexp(0.5, -250);
    auto matrix = sparse(2, {{0, 0, 1}, {1, 0, coupling}, {0, 1, coupling}, {1, 1, tiny}}, plan);
    auto transformed = linear::detail::sparse_exact_congruence(matrix, plan);
    if (!transformed) {
        FAIL("Exact congruence failed");
        return;
    }
    CHECK(transformed->scaled);
    const auto data = transformed->matrix.values();
    CHECK(data[0] == 1);
    CHECK(data[1] == 0.5);
    CHECK(data[2] == 0.5);
    CHECK(data[3] == 1);
    const double huge = std::numeric_limits<double>::max();
    const double subnormal = std::numeric_limits<double>::denorm_min();
    auto extreme =
        sparse(2, {{0, 0, huge}, {1, 0, subnormal}, {0, 1, subnormal}, {1, 1, huge}}, plan);
    auto preserved = linear::detail::sparse_exact_congruence(extreme, plan);
    if (!preserved) {
        FAIL("Fallback copy failed");
        return;
    }
    CHECK_FALSE(preserved->scaled);
    CHECK(std::equal(preserved->matrix.values().begin(), preserved->matrix.values().end(),
                     extreme.values().begin(), extreme.values().end()));
    core::ResourcePlan exhausted{0};
    auto rejected = linear::detail::sparse_exact_congruence(matrix, exhausted);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::ResourceExceeded);
    CHECK(matrix.values()[1] == coupling);
}

TEST_CASE("Shifted Cholesky certifies original non-diagonal-dominant SPD matrices", "[N4-05]") {
    core::ResourcePlan plan{100000};
    auto matrix = sparse(
        3, {{0, 0, 2}, {1, 0, 4}, {2, 0, -6}, {0, 1, 4}, {1, 1, 11}, {0, 2, -6}, {2, 2, 71}}, plan);
    auto proof = linear::detail::qualify_sparse_positive_definite(matrix, plan);
    if (!proof) {
        FAIL("SPD qualification did not complete");
        return;
    }
    CHECK(proof->certified);
    REQUIRE(proof->shift.has_value());
    if (proof->shift.has_value()) {
        CHECK(proof->shift.value() > 0);
    }
    CHECK(matrix.values()[0] == 2);
    CHECK(matrix.values()[1] == 4);
    auto tiny = sparse(2, {{0, 0, 1}, {1, 1, std::ldexp(1.0, -500)}}, plan);
    auto tiny_proof = linear::detail::qualify_sparse_positive_definite(tiny, plan);
    if (!tiny_proof) {
        FAIL("Scaled SPD qualification did not complete");
        return;
    }
    CHECK(tiny_proof->certified);
}

TEST_CASE("Failed shifted Cholesky remains uncertified instead of proving negative curvature",
          "[N4-05]") {
    core::ResourcePlan plan{100000};
    auto indefinite = sparse(2, {{0, 0, 1}, {1, 0, 2}, {0, 1, 2}, {1, 1, 1}}, plan);
    auto negative = linear::detail::qualify_sparse_positive_definite(indefinite, plan);
    if (!negative) {
        FAIL("Expected an unqualified diagnostic");
        return;
    }
    CHECK_FALSE(negative->certified);
    auto rigid = sparse(2, {{0, 0, 1}, {1, 0, -1}, {0, 1, -1}, {1, 1, 1}}, plan);
    auto null = linear::detail::qualify_sparse_positive_definite(rigid, plan);
    if (!null) {
        FAIL("Expected an unqualified diagnostic");
        return;
    }
    CHECK_FALSE(null->certified);
    auto nonsymmetric = sparse(2, {{0, 0, 2}, {1, 0, 1}, {0, 1, 2}, {1, 1, 3}}, plan);
    auto invalid = linear::detail::qualify_sparse_positive_definite(nonsymmetric, plan);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
}
