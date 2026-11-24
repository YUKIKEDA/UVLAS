#pragma once

#include "exact_polynomial_detail.hpp"
#include "pg_partition_test_support.hpp"

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <array>
#include <bit>
#include <cstdint>
#include <span>
#include <uvlas/structure/pg_conditioning.hpp>
#include <vector>

namespace pg_conditioning_test {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;

using Q = detail::Rational;
using I = detail::Integer;
using Dec = boost::multiprecision::cpp_dec_float_100;
constexpr auto profile = std::uint64_t{48} << 30U;
constexpr auto steps = 10'000'000'000'000ULL;
constexpr std::uint32_t dyadic_fractional_bits = 128;

inline PgReferenceLimits query_limits() {
    auto a = exact_limits();
    a.maximum_value_bits = 8192;
    a.maximum_work_bits = 32768;
    return {a, 20'000'000, 10'000'000, {}};
}

inline PgReferenceBasis make_basis(Fixture& f, const PgReferenceFrame& frame, core::Workspace& ws,
                                   core::ResourcePlan& plan) {
    auto dictionary = take(PgCoefficientDictionary::prepare(
        f.chart_graph(), {exact_limits(), 2000, 20000000, 100000, {}}, ws, plan));
    auto partition = take(PgIntegrationPartition::prepare(f.request(), limits(), ws, plan));
    PgReferenceInputs input{std::move(f.chart_graph()), std::move(dictionary),
                            std::move(partition)};
    auto geometry =
        take(PgReferenceGeometry::prepare(std::move(input), frame, query_limits(), plan));
    const PgRankPreparationLimits rank_limits{profile, steps, 128, {}};
    const auto source = take(geometry.source());
    auto rank = take(
        PgRankPreparation::prepare(source.dictionary, source.partition, rank_limits, ws, plan));
    return take(
        PgReferenceBasis::prepare(std::move(geometry), std::move(rank), rank_limits, ws, plan));
}

inline PgReferenceBasis make_basis(Fixture& f, const ExactRational& length, core::Workspace& ws,
                                   core::ResourcePlan& plan, bool cyclic_frame = false) {
    std::array<ExactRational, 9> rotation;
    rotation[cyclic_frame ? 2 : 0] = f.r("1");
    rotation[cyclic_frame ? 3 : 4] = f.r("1");
    rotation[cyclic_frame ? 7 : 8] = f.r("1");
    const std::array<ExactRational, 3> origin{};
    return make_basis(f, {rotation, origin, length}, ws, plan);
}

inline const Q& q(const ExactRational& value) {
    return detail::ExactFactory::value(value);
}
inline ExactRational rat(Fixture& f, const I& n, const I& d) {
    return f.r(n.str(), d.str());
}
inline std::size_t slot(std::size_t i, std::size_t j) {
    return j * (j + 1) / 2 + i;
}

struct Term {
    PgRankRow row;
    Q value;
};

inline std::vector<std::vector<Term>> basis_terms(const PgRankPreparation& rank, std::size_t count,
                                                  core::ResourcePlan& plan) {
    const auto rows = take(rank.rows());
    const auto algebra = query_limits().algebra;
    std::vector<std::vector<Term>> basis(count);
    for (std::size_t k = 0; k < count; ++k)
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto x = take(rank.basis({static_cast<core::Index>(k)},
                                           static_cast<core::Index>(i), algebra, plan));
            if (x.sign() != 0)
                basis[k].push_back({rows[i], q(x)});
        }
    return basis;
}

template <class Decimal = Dec>
std::vector<ExactRational> proposal(Fixture& f, const std::vector<Q>& gram, std::size_t count,
                                    std::uint32_t bits = dyadic_fractional_bits) {
    // Approximate test-side generator, independent of the checker.
    std::vector<Decimal> factor(gram.size()), inverse(gram.size());
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i) {
            const auto& x = gram[slot(i, j)];
            Decimal residual =
                numerator(x).convert_to<Decimal>() / denominator(x).convert_to<Decimal>();
            for (std::size_t k = 0; k < i; ++k)
                residual -= factor[slot(k, i)] * factor[slot(k, j)];
            if (i == j) {
                REQUIRE(residual > 0);
                factor[slot(i, j)] = sqrt(residual);
            } else {
                factor[slot(i, j)] = residual / factor[slot(i, i)];
            }
        }
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t k = j + 1; k > 0; --k) {
            const auto i = k - 1;
            Decimal residual = i == j ? 1 : 0;
            for (std::size_t l = i + 1; l <= j; ++l)
                residual -= factor[slot(i, l)] * inverse[slot(l, j)];
            inverse[slot(i, j)] = residual / factor[slot(i, i)];
        }
    const I scale = I{1} << bits;
    const auto decimal_scale = scale.convert_to<Decimal>();
    std::vector<ExactRational> result;
    for (const auto& value : inverse) {
        const Decimal integer = floor(value * decimal_scale + Decimal{"0.5"});
        // Boost's decimal-to-integer conversion negates ilogb(0)'s sentinel.
        const auto units = integer == 0 ? I{} : integer.template convert_to<I>();
        result.push_back(rat(f, units, scale));
    }
    return result;
}

template <class T>
void check_failure(const core::Result<T>& failure, core::FailureCode expected) {
    REQUIRE_FALSE(failure);
    CHECK(failure.error().code == expected);
}
inline void check_generated_required(const PgConditioningCertificate& accepted, std::uint32_t bits,
                                     std::span<const ExactRational> independent,
                                     core::Workspace& workspace, core::ResourcePlan& plan) {
    const auto info = take(accepted.info());
    const auto count = static_cast<std::size_t>(info.rank);
    const auto leading = static_cast<std::size_t>(info.protected_rank);
    REQUIRE(independent.size() == slot(0, count));
    const auto limits = query_limits();
    const auto* basis = take(accepted.reference_basis());
    const auto* rank = take(basis->rank());
    const Q ulp = Q{1} / (I{1} << bits);
    std::vector<Q> transform(independent.size());
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i) {
            const auto k = slot(i, j);
            transform[k] = q(take(accepted.transform_coefficient({static_cast<core::Index>(i)},
                                                                 {static_cast<core::Index>(j)},
                                                                 limits.algebra, plan)));
            CHECK(abs(transform[k] - q(independent[k])) <= ulp);
        }
    for (core::Index t = 0; t < 33; ++t) {
        const auto coordinates = take(accepted.required_coordinates({t}, limits, workspace, plan));
        REQUIRE(coordinates.size() == count);
        for (std::size_t i = 0; i < count; ++i) {
            if (i >= leading)
                CHECK(coordinates.values()[i].sign() == 0);
            Q restored{};
            for (std::size_t j = i; j < leading; ++j)
                restored += transform[slot(i, j)] * q(coordinates.values()[j]);
            const auto expected =
                take(rank->required({static_cast<core::Index>(i)}, {t}, limits.algebra, plan));
            CHECK(restored == q(expected));
        }
    }
}

inline Q power(const Q& x, core::Index exponent) {
    Q result{1};
    for (core::Index i = 0; i < exponent; ++i)
        result *= x;
    return result;
}

inline Q independent_double(double value) {
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto biased = static_cast<unsigned>((bits >> 52U) & 0x7ffU);
    REQUIRE(biased != 0x7ffU);
    const auto fraction = bits & 0x000fffffffffffffULL;
    const auto significand = biased == 0 ? fraction : fraction | (1ULL << 52U);
    const auto exponent = biased == 0 ? -1074 : static_cast<int>(biased) - 1075;
    I n{significand};
    if ((bits >> 63U) != 0)
        n = -n;
    const auto shift = static_cast<unsigned>(std::abs(exponent));
    return exponent >= 0 ? Q{n * (I{1} << shift)} : Q{n} / (I{1} << shift);
}

} // namespace pg_conditioning_test
