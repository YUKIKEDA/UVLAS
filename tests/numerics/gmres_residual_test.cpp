#include "gmres_residual_detail.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

namespace {
using namespace uvlas;
using numerics::linear::detail::combine_gmres_residual;
using Exact = boost::multiprecision::cpp_rational;
} // namespace

TEST_CASE("GMRES residual enclosures include both signed tails and producer error endpoints",
          "[verification][N4-05][N4-06]") {
    constexpr std::size_t count = 6;
    const std::array<double, count> bh{std::ldexp(1., 60), 1., std::ldexp(1., -1022), 0., 1., 1.};
    const std::array<double, count> bl{1.,  std::ldexp(1., -80), std::ldexp(1., -1074), 0.,
                                       -1., std::ldexp(1., -80)};
    const std::array<double, count> ah{std::ldexp(1., 60), 1., std::ldexp(1., -1022), 0., 0., 0.};
    const std::array<double, count> al{0., std::ldexp(1., -81), 0., 0., 0., 0.};
    const std::array<double, count> ath{0., -std::ldexp(1., -83), 0., 0., 0., 0.};
    const std::array<double, count> atl{0., std::ldexp(1., -84), 0., 0., 0., 0.};
    const std::array<double, count> eh{0., std::ldexp(1., -90), std::ldexp(1., -1074), 0., 0., 0.};
    const std::array<double, count> el{0., std::ldexp(1., -100), 0., 0., 0., 0.};
    std::array<double, count> residual{}, error{};
    const auto combined =
        combine_gmres_residual(6, bh, bl, ah, al, ath, atl, eh, el, residual, error);
    REQUIRE(combined);
    for (std::size_t i = 0; i < count; ++i) {
        for (const int high_sign : {-1, 1}) {
            for (const int low_sign : {-1, 1}) {
                // Exact rationals are independent even at the subnormal boundary.
                const Exact truth = Exact{bh[i]} + Exact{bl[i]} - Exact{ah[i]} - Exact{ath[i]} -
                                    Exact{al[i]} - Exact{atl[i]} - high_sign * Exact{eh[i]} -
                                    low_sign * Exact{el[i]};
                CHECK(abs(truth - Exact{residual[i]}) <= Exact{error[i]});
            }
        }
    }
    CHECK(residual[0] == 1.);
    CHECK(error[0] == 0.);
    CHECK(residual[2] == std::numeric_limits<double>::denorm_min());
    CHECK(residual[3] == 0.);
    CHECK(error[3] == 0.);
    CHECK(residual[4] == 0.);
    CHECK(error[4] == 0.);
    CHECK(residual[5] == 1.);
    CHECK(error[5] >= std::ldexp(1., -80));
}

TEST_CASE("GMRES residual rejects malformed inputs before changing output scratch",
          "[negative][N4-05][N4-08]") {
    const std::array<double, 8> valid{1., 0., 1., 0., 0., 0., 0., 0.};
    std::array<double, 1> residual{42.}, error{43.};
    const auto run = [&](const std::array<double, 8>& data, core::Index rows = 1) {
        const std::span all{data};
        return combine_gmres_residual(rows, all.subspan(0, 1), all.subspan(1, 1), all.subspan(2, 1),
                                      all.subspan(3, 1), all.subspan(4, 1), all.subspan(5, 1),
                                      all.subspan(6, 1), all.subspan(7, 1), residual, error);
    };
    for (std::size_t site = 0; site < valid.size(); ++site) {
        auto bad = valid;
        bad[site] = std::numeric_limits<double>::quiet_NaN();
        const auto rejected = run(bad);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::Nonfinite);
        CHECK(residual[0] == 42.);
        CHECK(error[0] == 43.);
    }
    for (const std::size_t site : {std::size_t{6}, std::size_t{7}}) {
        auto bad = valid;
        bad[site] = -1.;
        const auto rejected = run(bad);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    }
    for (const core::Index rows : {-1, 0, 2}) {
        const auto rejected = run(valid, rows);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    }
    CHECK(combine_gmres_residual(0, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}));
    const std::array<double, 2> late_bad{1., std::numeric_limits<double>::quiet_NaN()}, zeros{};
    std::array<double, 2> out{42., 44.}, bounds{43., 45.};
    const auto late = combine_gmres_residual(2, late_bad, zeros, zeros, zeros, zeros, zeros, zeros,
                                             zeros, out, bounds);
    REQUIRE_FALSE(late);
    CHECK(late.error().code == core::FailureCode::Nonfinite);
    CHECK((out == std::array{42., 44.}));
    CHECK((bounds == std::array{43., 45.}));
    for (std::size_t site = 0; site < 10; ++site) {
        std::array<std::span<const double>, 8> inputs;
        inputs.fill(std::span<const double>{zeros}.first(1));
        auto output = std::span{out}.first(1), bound = std::span{bounds}.first(1);
        if (site < 8)
            inputs[site] = {};
        else if (site == 8)
            output = {};
        else
            bound = {};
        const auto mismatch =
            combine_gmres_residual(1, inputs[0], inputs[1], inputs[2], inputs[3], inputs[4],
                                   inputs[5], inputs[6], inputs[7], output, bound);
        REQUIRE_FALSE(mismatch);
        CHECK(mismatch.error().code == core::FailureCode::InvalidInput);
    }
}

TEST_CASE("Finite GMRES residual overflow is unresolved rather than a nonfinite input",
          "[negative][N4-05][N4-06]") {
    const std::array<double, 1> huge{std::numeric_limits<double>::max()}, zero{};
    std::array<double, 1> residual{}, error{};
    const auto arithmetic =
        combine_gmres_residual(1, huge, huge, huge, huge, zero, zero, zero, zero, residual, error);
    REQUIRE_FALSE(arithmetic);
    CHECK(arithmetic.error().code == core::FailureCode::PrecisionUnresolved);
    const auto enclosure =
        combine_gmres_residual(1, zero, zero, zero, zero, zero, zero, huge, huge, residual, error);
    REQUIRE_FALSE(enclosure);
    CHECK(enclosure.error().code == core::FailureCode::PrecisionUnresolved);
}
