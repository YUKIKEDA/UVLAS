#pragma once

#include "pg_coefficients_maps.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace uvlas::structure::detail {
inline std::uint64_t pg_add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        PgCoefficientGeneration::exceeded();
    return a + b;
}
inline std::uint64_t pg_mul(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        PgCoefficientGeneration::exceeded();
    return a * b;
}
inline std::uint64_t pg_log(std::uint64_t value) noexcept {
    // bit_width of an unsigned value is in [0,64]; widening preserves its value.
    return value <= 1 ? 0 : static_cast<std::uint64_t>(std::bit_width(value - 1));
}
inline std::uint64_t pg_monomials(std::uint64_t degree) {
    return pg_mul(pg_add(degree, 1), pg_add(degree, 2)) / 2;
}
/// Native denominators use only these primes. Other input factors remain in one
/// opaque common denominator; none is approximated, factored with heap scratch,
/// or assumed absent because its bit length is large.
struct PgDenominatorBound {
    std::array<std::uint64_t, 6> prime{};
    std::uint64_t opaque = 0;
};
inline PgDenominatorBound pg_den_join(const PgDenominatorBound& a, const PgDenominatorBound& b) {
    PgDenominatorBound result;
    result.opaque = std::max(a.opaque, b.opaque);
    for (std::size_t j = 0; j < 6; ++j)
        result.prime[j] = std::max(a.prime[j], b.prime[j]);
    return result;
}
inline PgDenominatorBound pg_den_product(const PgDenominatorBound& a, const PgDenominatorBound& b) {
    PgDenominatorBound result;
    result.opaque = pg_add(a.opaque, b.opaque);
    for (std::size_t j = 0; j < 6; ++j)
        result.prime[j] = pg_add(a.prime[j], b.prime[j]);
    return result;
}
inline PgDenominatorBound pg_den_native(std::uint64_t value) {
    constexpr std::array<std::uint64_t, 6> primes{2, 3, 5, 7, 11, 13};
    PgDenominatorBound result;
    for (std::size_t j = 0; j < 6; ++j)
        while (value > 1 && value % primes[j] == 0) {
            ++result.prime[j];
            value /= primes[j];
        }
    if (value > 1)
        result.opaque = 1;
    return result;
}
inline PgDenominatorBound pg_den_lcm_range(std::uint64_t last) {
    if (last > 16)
        PgCoefficientGeneration::exceeded();
    PgDenominatorBound result;
    for (std::uint64_t j = 1; j <= last; ++j)
        result = pg_den_join(result, pg_den_native(j));
    return result;
}
inline std::uint64_t pg_native_backend(const IntegerBackend& value) noexcept {
    if constexpr (IntegerBackend::limb_bits == 64)
        return value.size() == 0 ? 0 : static_cast<std::uint64_t>(value.limbs()[0]);
    else {
        std::uint64_t result = 0;
        for (std::size_t j = value.size(); j > 0; --j)
            result = (result << 32) | static_cast<std::uint64_t>(value.limbs()[j - 1]);
        return result;
    }
}
struct PgBudgetPolynomial {
    std::uint64_t terms = 0, first = 0, second = 0, degree = 0, l1 = 0;
    PgDenominatorBound denominator;
    std::uint64_t storage = 0;
};
/// Pure arithmetic forecast for this fixed degree-four/P7 construction. It uses
/// no dynamic integers or arrays. Generic leaf admission and actual storage
/// bounds are distinct: U(V+1) remains required by the existing exact API.
struct PgCoefficientBudget {
    const PgCoefficientLimits& limits;
    std::uint64_t opaque_bits = 0, payload = 0, live_external = 0, leaf_external = 0;
    std::uint64_t maximum_work = 0, maximum_raw = 0, output_terms = 0;
    std::uint64_t maximum_storage = 0, maximum_value = 0, canonical_terms = 0;
    std::uint64_t first_degree = 0, second_degree = 0, total_degree = 0;
    std::uint64_t repeat = 1;
    bool charging = true;
    explicit PgCoefficientBudget(const PgCoefficientLimits& caps) : limits(caps) {}
    static PgDenominatorBound input_denominator(const IntegerBackend& value, std::uint64_t bits) {
        if (bits <= 64)
            return pg_den_native(pg_native_backend(value));
        PgDenominatorBound result;
        result.opaque = 1;
        return result;
    }
    void opaque_input(const IntegerBackend& value, std::uint64_t bits) {
        if (bits > 64) {
            opaque_bits = pg_add(opaque_bits, bits);
            return;
        }
        constexpr std::array<std::uint64_t, 6> primes{2, 3, 5, 7, 11, 13};
        auto native = pg_native_backend(value);
        for (const auto prime : primes)
            while (native > 1 && native % prime == 0)
                native /= prime;
        if (native > 1)
            opaque_bits = pg_add(opaque_bits, static_cast<std::uint64_t>(std::bit_width(native)));
    }
    [[nodiscard]] std::uint64_t value_bits(const PgBudgetPolynomial& value) const {
        if (value.terms == 0)
            return 0;
        constexpr std::array<std::uint64_t, 6> weights{1, 2, 3, 3, 4, 4};
        auto result = pg_add(value.l1, pg_mul(value.denominator.opaque, opaque_bits));
        for (std::size_t j = 0; j < 6; ++j)
            result = pg_add(result, pg_mul(value.denominator.prime[j], weights[j]));
        return pg_add(result, 1);
    }
    static std::uint64_t capacity(std::uint64_t bits) {
        if (bits == 0)
            return 0;
        constexpr auto w = IntegerBackend::limb_bits;
        const auto words = pg_add(pg_add(bits / w, bits % w != 0 ? 1 : 0), 2);
        return words <= IntegerBackend::internal_limb_count ? 0 : pg_mul(4, words);
    }
    [[nodiscard]] std::uint64_t external_owner(const PgBudgetPolynomial& value) const {
        return pg_mul(pg_mul(pg_mul(2, value.terms), capacity(value.storage)),
                      sizeof(boost::multiprecision::limb_type));
    }
    [[nodiscard]] std::uint64_t generic_scratch(std::uint64_t work, std::uint64_t raw) const {
        if (work == 0)
            return 0;
        constexpr auto w = IntegerBackend::limb_bits;
        const auto words = pg_add(pg_add(work / w, work % w != 0 ? 1 : 0), 2);
        const auto c = words <= IntegerBackend::internal_limb_count ? 0 : pg_mul(4, words);
        const auto r = c == 0 ? 0 : words;
        const auto gcd = pg_add(pg_mul(6, pg_add(words, 1)), pg_mul(5, c));
        const auto division = pg_mul(2, c), karatsuba = words >= 60 ? pg_mul(5, words) : 0;
        return pg_mul(pg_add(pg_add(pg_mul(pg_add(16, pg_mul(6, raw)), c),
                                    std::max({gcd, division, karatsuba})),
                             r),
                      sizeof(boost::multiprecision::limb_type));
    }
    void retain(const PgBudgetPolynomial& value, std::uint64_t count, bool published = false) {
        live_external = pg_add(live_external, pg_mul(external_owner(value), count));
        if (published)
            output_terms = pg_add(output_terms, pg_mul(value.terms, count));
    }
    void charge(std::uint64_t raw, std::uint64_t canonical, std::uint64_t raw_bits,
                std::uint64_t extra_raw = 0) {
        if (!charging)
            return;
        const auto work = raw == 0 ? 0 : pg_mul(raw, pg_add(raw_bits, 1));
        maximum_work = std::max(maximum_work, work);
        maximum_raw = std::max(maximum_raw, raw);
        leaf_external = std::max(leaf_external, generic_scratch(work, raw));
        const auto owned =
            pg_add(pg_mul(pg_add(raw, extra_raw), exact_polynomial_raw_term_bytes()),
                   pg_add(pg_mul(canonical, pg_add(sizeof(ExactTerm), sizeof(RationalStorage))),
                          canonical == 0 ? 0 : sizeof(PolynomialStorage)));
        payload = pg_add(payload, pg_mul(owned, repeat));
    }
    PgBudgetPolynomial finish(PgBudgetPolynomial result, std::uint64_t raw, std::uint64_t raw_bits,
                              std::uint64_t operand_bits = 0) {
        if (result.terms != 0) {
            result.terms = std::min({result.terms, pg_monomials(result.degree),
                                     pg_mul(pg_add(result.first, 1), pg_add(result.second, 1))});
            result.storage = pg_add(pg_mul(2, std::max(value_bits(result), operand_bits)), 1);
        }
        if (charging) {
            maximum_storage = std::max(maximum_storage, result.storage);
            maximum_value = std::max(maximum_value, value_bits(result));
            canonical_terms = std::max(canonical_terms, result.terms);
            first_degree = std::max(first_degree, result.first);
            second_degree = std::max(second_degree, result.second);
            total_degree = std::max(total_degree, result.degree);
        }
        charge(raw, result.terms, raw_bits);
        return result;
    }
    PgBudgetPolynomial number(std::uint64_t numerator, std::uint64_t denominator = 1) {
        PgBudgetPolynomial result;
        if (numerator == 0)
            return result;
        result.terms = 1;
        result.l1 = pg_log(numerator);
        result.denominator = pg_den_native(denominator);
        pg_require_identity(result.denominator.opaque == 0);
        result.storage = value_bits(result);
        return result;
    }
    PgBudgetPolynomial input(const ExactRational& value, bool inverse = false) {
        if (value.sign() == 0)
            return {};
        const auto& v = ExactFactory::value(value).backend();
        PgBudgetPolynomial result;
        result.terms = 1;
        result.l1 = inverse ? value.denominator_bits() : value.numerator_bits();
        result.denominator = inverse ? input_denominator(v.num(), value.numerator_bits())
                                     : input_denominator(v.denom(), value.denominator_bits());
        result.storage = value_bits(result);
        return result;
    }
    PgBudgetPolynomial coordinate(ExactAxis axis, std::uint64_t n) {
        auto result = number(1, 4 * n);
        result.terms = 2;
        result.degree = 1;
        result.first = axis == ExactAxis::First ? 1 : 0;
        result.second = axis == ExactAxis::Second ? 1 : 0;
        result.storage = pg_add(pg_mul(2, value_bits(result)), 1);
        return result;
    }
    PgBudgetPolynomial a(const PgBudgetPolynomial& x, const PgBudgetPolynomial& y) {
        PgBudgetPolynomial r;
        r.terms = pg_add(x.terms, y.terms);
        r.first = std::max(x.first, y.first);
        r.second = std::max(x.second, y.second);
        r.degree = std::max(x.degree, y.degree);
        r.l1 = pg_add(std::max(x.l1, y.l1), x.terms != 0 && y.terms != 0 ? 1 : 0);
        r.denominator = pg_den_join(x.denominator, y.denominator);
        return finish(r, r.terms, pg_add(std::max(value_bits(x), value_bits(y)), 1),
                      std::max(value_bits(x), value_bits(y)));
    }
    PgBudgetPolynomial m(const PgBudgetPolynomial& x, const PgBudgetPolynomial& y) {
        PgBudgetPolynomial r;
        if (x.terms == 0 || y.terms == 0)
            return finish(r, 0, 0);
        r.terms = pg_mul(x.terms, y.terms);
        r.first = pg_add(x.first, y.first);
        r.second = pg_add(x.second, y.second);
        r.degree = pg_add(x.degree, y.degree);
        r.l1 = pg_add(x.l1, y.l1);
        r.denominator = pg_den_product(x.denominator, y.denominator);
        return finish(r, r.terms, pg_add(value_bits(x), value_bits(y)),
                      std::max(value_bits(x), value_bits(y)));
    }
    PgBudgetPolynomial s(const PgBudgetPolynomial& x, const PgBudgetPolynomial& scalar) {
        auto r = x;
        if (scalar.terms == 0 || x.terms == 0)
            return finish({}, 0, 0);
        r.l1 = pg_add(x.l1, scalar.l1);
        r.denominator = pg_den_product(x.denominator, scalar.denominator);
        return finish(r, x.terms, pg_add(value_bits(x), value_bits(scalar)),
                      std::max(value_bits(x), value_bits(scalar)));
    }
    PgBudgetPolynomial d(const PgBudgetPolynomial& x, ExactAxis axis) {
        auto r = x;
        auto& power = axis == ExactAxis::First ? r.first : r.second;
        if (power == 0 || x.terms == 0)
            return finish({}, 0, 0);
        const auto growth = static_cast<std::uint64_t>(std::bit_width(power));
        --power;
        --r.degree;
        r.l1 = pg_add(x.l1, growth);
        return finish(r, x.terms, pg_add(value_bits(x), growth), value_bits(x));
    }
    PgBudgetPolynomial i(const PgBudgetPolynomial& x, ExactAxis axis) {
        auto r = x;
        if (x.terms == 0)
            return finish({}, 0, 0);
        auto& power = axis == ExactAxis::First ? r.first : r.second;
        ++power;
        ++r.degree;
        r.terms = pg_mul(2, x.terms);
        r.l1 = pg_add(x.l1, 1);
        r.denominator = pg_den_product(x.denominator, pg_den_lcm_range(power));
        const auto raw_bits =
            pg_add(pg_add(value_bits(x), static_cast<std::uint64_t>(std::bit_width(power))), power);
        return finish(r, r.terms, raw_bits, value_bits(x));
    }
    PgBudgetPolynomial clone(const PgBudgetPolynomial& x) {
        auto r = x;
        charge(x.terms, x.terms, value_bits(x));
        r.storage = value_bits(x);
        return r;
    }
    PgBudgetPolynomial restriction(const PgBudgetPolynomial& x) {
        auto r = x;
        r.first = 0;
        r.degree = r.second;
        r.terms = std::min(x.terms, pg_add(r.second, 1));
        return finish(r, x.terms, value_bits(x), value_bits(x));
    }
};
inline PgBudgetPolynomial pg_budget_union(const PgBudgetPolynomial& a,
                                          const PgBudgetPolynomial& b) {
    PgBudgetPolynomial r;
    r.terms = std::max(a.terms, b.terms);
    r.first = std::max(a.first, b.first);
    r.second = std::max(a.second, b.second);
    r.degree = std::max(a.degree, b.degree);
    r.l1 = std::max(a.l1, b.l1);
    r.denominator = pg_den_join(a.denominator, b.denominator);
    r.storage = std::max(a.storage, b.storage);
    return r;
}
struct PgBudgetGeometry {
    PgBudgetPolynomial s, t, f, fx, fy, U, V, xy, g;
    std::array<PgBudgetPolynomial, 5> s_power, t_power;
};
struct PgBudgetPair {
    std::array<PgBudgetPolynomial, 3> v, p;
};
inline std::array<PgBudgetPolynomial, 5> pg_budget_powers(PgCoefficientBudget& b,
                                                          const PgBudgetPolynomial& x) {
    std::array<PgBudgetPolynomial, 5> p;
    p[0] = b.number(1);
    b.charge(1, 1, 1);
    p[1] = b.clone(x);
    for (std::size_t j = 2; j < 5; ++j)
        p[j] = b.m(p[j - 1], x);
    return p;
}
inline PgBudgetPolynomial pg_budget_curvature() {
    PgBudgetPolynomial p;
    p.terms = 4;
    p.first = 7;
    p.degree = 7;
    p.l1 = 7;
    p.denominator.prime[0] = 4;
    p.storage = 25;
    return p;
}
inline PgBudgetGeometry pg_budget_geometry(PgCoefficientBudget& b, const PgGraphView& graph) {
    const auto n = static_cast<std::uint64_t>(graph.n);
    PgBudgetGeometry g;
    g.s = b.coordinate(ExactAxis::First, n);
    g.t = b.coordinate(ExactAxis::Second, n);
    b.charge(2, 2, 5);
    b.charge(2, 2, 5);
    g.s_power = pg_budget_powers(b, g.s);
    g.t_power = pg_budget_powers(b, g.t);
    const auto active = b.charging;
    b.charging = false;
    std::uint64_t raw = 0, maximum = 1;
    for (const auto& term : graph.height.terms()) {
        const auto a = static_cast<std::size_t>(term.powers.first),
                   c = static_cast<std::size_t>(term.powers.second);
        const auto piece = b.s(b.m(g.s_power[a], g.t_power[c]), b.input(term.coefficient));
        g.f = b.a(g.f, piece);
        raw = pg_add(raw, std::uint64_t{1} << static_cast<unsigned>(a + c));
        maximum = std::max(maximum, pg_add(std::max(term.coefficient.numerator_bits(),
                                                    term.coefficient.denominator_bits()),
                                           pg_mul(a + c, 5)));
    }
    b.charging = active;
    // substitute allocates global raw U plus (2U-m) expansion raw slots.
    b.charge(raw, g.f.terms, maximum,
             raw == 0 ? 0
                      : pg_add(pg_mul(2, raw), 0) -
                            static_cast<std::uint64_t>(graph.height.terms().size()));
    b.maximum_storage = std::max(b.maximum_storage, g.f.storage);
    b.maximum_value = std::max(b.maximum_value, b.value_bits(g.f));
    const auto scale = b.number(4 * n), inverse = b.input(graph.width, true),
               skew = b.input(graph.skew);
    g.fx = b.s(b.d(g.f, ExactAxis::First), b.number(4 * n));
    g.fx = b.s(g.fx, inverse);
    g.fy = b.a(b.s(b.d(g.f, ExactAxis::Second), scale), b.s(g.fx, skew));
    const auto xx = b.m(g.fx, g.fx), yy = b.m(g.fy, g.fy);
    g.U = b.a(b.number(1), xx);
    g.V = b.a(b.number(1), yy);
    g.xy = b.m(g.fx, g.fy);
    g.g = b.a(g.U, yy);
    return g;
}
inline std::array<PgBudgetPolynomial, 5> pg_budget_moments(PgCoefficientBudget& b, std::uint64_t n,
                                                           ExactAxis axis, bool completed) {
    auto coordinate = b.coordinate(axis, n);
    const auto powers = pg_budget_powers(b, coordinate);
    auto curvature = pg_budget_curvature();
    if (axis == ExactAxis::Second) {
        curvature.second = 7;
        curvature.first = 0;
    }
    const auto h = b.number(1, 4 * n);
    std::array<PgBudgetPolynomial, 5> result;
    b.charge(2, 2, 5); // seed coordinate affine; native constants are counted separately.
    for (std::size_t m = 0; m < 5; ++m) {
        const auto integrated = b.i(b.s(b.m(powers[m], curvature), h), axis);
        result[m] = integrated;
        const auto evaluation_bits = pg_add(b.value_bits(integrated), integrated.degree);
        const auto work = pg_mul(integrated.terms, pg_add(evaluation_bits, 1));
        if (b.charging) {
            b.maximum_work = std::max(b.maximum_work, work);
            b.maximum_raw = std::max(b.maximum_raw, integrated.terms);
            b.leaf_external = std::max(b.leaf_external, b.generic_scratch(work, integrated.terms));
            b.payload = pg_add(b.payload, pg_mul(sizeof(RationalStorage), b.repeat));
        }
        auto endpoint = integrated;
        endpoint.terms = integrated.terms == 0 ? 0 : 1;
        endpoint.first = 0;
        endpoint.second = 0;
        endpoint.degree = 0;
        b.charge(1, endpoint.terms, b.value_bits(endpoint));
        endpoint.storage = b.value_bits(endpoint);
        if (completed)
            result[m] = endpoint;
    }
    return result;
}
inline PgBudgetPair pg_budget_pair(PgCoefficientBudget& b, const PgGraphView& graph,
                                   const PgBudgetGeometry& g,
                                   const std::array<PgBudgetPolynomial, 5>& moments,
                                   PgPairedAxis axis) {
    const auto& transverse = axis == PgPairedAxis::S ? g.t_power : g.s_power;
    const auto& along = axis == PgPairedAxis::S ? g.s : g.t;
    PgBudgetPolynomial J, Jt, Jst;
    for (const auto& term : graph.height.terms()) {
        const auto i = static_cast<std::size_t>(axis == PgPairedAxis::S ? term.powers.first
                                                                        : term.powers.second);
        const auto j = static_cast<std::size_t>(axis == PgPairedAxis::S ? term.powers.second
                                                                        : term.powers.first);
        const auto coefficient = b.input(term.coefficient);
        J = b.a(J, b.s(b.m(transverse[j], moments[i]), coefficient));
        if (j != 0) {
            auto differentiated = coefficient;
            differentiated.l1 = pg_add(differentiated.l1, pg_log(j));
            Jt = b.a(Jt, b.s(b.m(transverse[j - 1], moments[i]), differentiated));
            Jst = b.a(Jst, b.s(b.m(transverse[j - 1], moments[i + 1]), differentiated));
        }
    }
    auto first = b.a(J, b.m(g.f, moments[0]));
    auto second = b.a(Jst, b.m(along, Jt));
    auto third = b.a(b.m(along, moments[0]), moments[1]);
    if (axis == PgPairedAxis::T)
        std::swap(first, second);
    const auto inverse = b.input(graph.width, true), skew = b.input(graph.skew),
               scale = b.number(static_cast<std::uint64_t>(4 * graph.n));
    first = b.s(first, inverse);
    second = b.a(second, b.s(first, skew));
    PgBudgetPair result;
    result.v = {first, second, third};
    std::array<PgBudgetPolynomial, 3> dx, dy;
    for (std::size_t j = 0; j < 3; ++j) {
        auto sx = scale;
        sx.l1 = pg_add(sx.l1, inverse.l1);
        sx.denominator = pg_den_product(sx.denominator, inverse.denominator);
        dx[j] = b.s(b.d(result.v[j], ExactAxis::First), sx);
        dy[j] = b.a(b.s(b.d(result.v[j], ExactAxis::Second), scale), b.s(dx[j], skew));
    }
    const auto tx = b.a(b.a(dx[2], b.m(g.fx, dx[0])), b.m(g.fy, dx[1]));
    const auto ty = b.a(b.a(dy[2], b.m(g.fx, dy[0])), b.m(g.fy, dy[1]));
    result.p[0] = b.a(b.m(g.xy, ty), b.m(g.V, tx));
    result.p[1] = b.a(b.m(g.xy, tx), b.m(g.U, ty));
    result.p[2] = b.a(b.m(g.fx, tx), b.m(g.fy, ty));
    return result;
}
inline std::uint64_t pg_actual_external(const ExactRational& value) {
    const auto& v = ExactFactory::value(value).backend();
    const auto heap = [](const IntegerBackend& i) {
        return i.capacity() > IntegerBackend::internal_limb_count
                   ? static_cast<std::uint64_t>(i.capacity())
                   : 0;
    };
    return pg_mul(pg_add(heap(v.num()), heap(v.denom())), sizeof(boost::multiprecision::limb_type));
}
inline PgBudgetPolynomial pg_budget_scalar(PgCoefficientBudget& b, const PgBudgetPolynomial& p) {
    auto scalar = p;
    scalar.terms = p.terms == 0 ? 0 : 1;
    scalar.first = 0;
    scalar.second = 0;
    scalar.degree = 0;
    scalar.storage = b.value_bits(scalar);
    return scalar;
}
inline void pg_budget_evaluation(PgCoefficientBudget& b, const PgBudgetPolynomial& p) {
    const auto raw_bits = pg_add(b.value_bits(p), p.degree);
    const auto work = p.terms == 0 ? 0 : pg_mul(p.terms, pg_add(raw_bits, 1));
    b.maximum_work = std::max(b.maximum_work, work);
    b.maximum_raw = std::max(b.maximum_raw, p.terms);
    b.leaf_external = std::max(b.leaf_external, b.generic_scratch(work, p.terms));
    b.payload = pg_add(b.payload, pg_mul(sizeof(RationalStorage), b.repeat));
}
inline PgBudgetPolynomial pg_budget_global_curvature(PgCoefficientBudget& b, std::uint64_t n,
                                                     std::uint64_t degree) {
    auto q = b.coordinate(ExactAxis::First, n);
    q.denominator = pg_den_native(2 * n);
    std::array<PgBudgetPolynomial, 8> powers;
    powers[0] = b.number(1);
    b.charge(1, 1, 1);
    powers[1] = b.clone(q);
    b.charge(2, 2, 5);
    for (std::size_t j = 2; j < 8; ++j)
        powers[j] = b.m(powers[j - 1], q);
    PgBudgetPolynomial result;
    const auto leading = pg_budget_curvature();
    // Every P_k has at most four terms; only the parity matching k is present.
    for (std::uint64_t power = degree % 2; power <= degree; power += 2) {
        auto coefficient = leading;
        coefficient.terms = 1;
        coefficient.first = 0;
        coefficient.degree = 0;
        result = b.a(result, b.s(powers[static_cast<std::size_t>(power)], coefficient));
    }
    return result;
}
inline std::array<PgBudgetPolynomial, 8> pg_budget_decompose(PgCoefficientBudget& b,
                                                             const PgBudgetPolynomial& curvature,
                                                             std::uint64_t degree) {
    constexpr std::array<std::uint64_t, 8> lead_num{1, 1, 3, 5, 35, 63, 231, 429};
    constexpr std::array<std::uint64_t, 8> lead_den{1, 1, 2, 2, 8, 8, 16, 16};
    auto remainder = curvature;
    std::array<PgBudgetPolynomial, 8> result;
    for (std::uint64_t offset = 0; offset <= degree; ++offset) {
        const auto j = degree - offset;
        auto coefficient = pg_budget_scalar(b, remainder);
        coefficient.denominator =
            pg_den_product(coefficient.denominator, pg_den_native(lead_num[j]));
        coefficient.l1 = pg_add(coefficient.l1, pg_log(lead_den[j]));
        coefficient.storage = b.value_bits(coefficient);
        result[j] = coefficient;
        auto basis = pg_budget_curvature();
        basis.degree = j;
        basis.first = j;
        basis.terms = j / 2 + 1;
        remainder = b.a(remainder, b.s(basis, coefficient));
        // The chosen leading coefficient cancels exactly. Product verification
        // also checks this identity; the forecast does not remove a small value.
        if (j == 0)
            remainder = {};
        else {
            remainder.first = j - 1;
            remainder.degree = j - 1;
            remainder.terms = std::min(remainder.terms, j);
        }
    }
    return result;
}
inline std::uint64_t pg_payload_poly(std::uint64_t terms) {
    return pg_add(pg_mul(terms, pg_add(pg_add(exact_polynomial_raw_term_bytes(), sizeof(ExactTerm)),
                                       sizeof(RationalStorage))),
                  terms == 0 ? 0 : sizeof(PolynomialStorage));
}
inline void pg_budget_legendre(PgCoefficientBudget& b) {
    const auto x = b.coordinate(ExactAxis::First, 1);
    const auto l = pg_budget_curvature();
    b.charge(2, 2, 1);
    b.charge(1, 1, 1);
    (void)b.clone(x);
    for (std::uint64_t k = 1; k < 7; ++k) {
        auto left = b.s(b.m(x, l), b.number(2 * k + 1));
        auto right = b.s(l, b.number(k));
        (void)b.s(b.a(left, right), b.number(1, k + 1));
    }
}
inline void pg_budget_basis(PgCoefficientBudget& b, std::uint64_t n) {
    const auto axis = 2 * n;
    const auto saved = b.repeat;
    auto linear = b.coordinate(ExactAxis::First, 1);
    linear.denominator = pg_den_native(8);
    auto previous = linear;
    previous.l1 = 1;
    b.repeat = pg_mul(2 * axis, 6 * n + 8);
    const auto first = b.m(linear, previous), second = b.m(linear, previous);
    (void)b.a(first, second);
    b.charge(2, 2, 5);
    b.charge(2, 2, 5);
    b.repeat = pg_mul(2 * axis, 3 * n + 7);
    b.charge(1, 1, 1);
    b.repeat = 4 * axis;
    b.charge(2, 2, 5);
    b.repeat = saved;
}
inline void pg_budget_native_payload(PgCoefficientBudget& b, std::uint64_t calls,
                                     std::uint64_t value_bits) {
    b.payload = pg_add(b.payload, pg_mul(calls, sizeof(RationalStorage)));
    const auto work = std::max(std::uint64_t{9}, pg_add(pg_mul(2, value_bits), 1));
    b.maximum_work = std::max(b.maximum_work, work);
    b.maximum_storage = std::max(b.maximum_storage, work);
    // 16+6 here is a conservative superset of the existing rational slots 7/6;
    // GCD/division/resize are still included, never measured allocation counts.
    b.leaf_external = std::max(b.leaf_external, b.generic_scratch(work, 1));
    b.maximum_value = std::max(b.maximum_value, value_bits);
}
/// Pure whole-call admission. Every repetition is a fixed family/cardinality
/// from the contract; no Model-scale dense rank or numerical candidate is built.
inline PgCoefficientPlan pg_coefficient_plan(const PgGraphView& graph,
                                             const PgCoefficientLimits& limits,
                                             core::ResourcePlan& owned) {
    const auto& l = limits.algebra;
    if (l.external_memory.valueless_by_exception() ||
        (l.timing != ExactTimingMode::Disabled && l.timing != ExactTimingMode::Enabled))
        throw PgCoefficientAbort{{core::FailureCode::InvalidInput, core::Operation::PrepareModel}};
    constexpr auto index_max = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    for (const auto value :
         {l.maximum_input_characters, l.maximum_terms, l.maximum_raw_terms, l.maximum_axis_degree,
          l.maximum_total_degree, l.maximum_value_bits, l.maximum_work_bits, limits.maximum_columns,
          limits.maximum_stored_terms, limits.maximum_map_entries})
        if (value > index_max)
            PgCoefficientGeneration::exceeded();
    if (graph.n != 1 && graph.n != 2 && graph.n != 4)
        throw PgCoefficientAbort{{core::FailureCode::InvalidInput, core::Operation::PrepareModel}};
    if (std::holds_alternative<ExactRequireProvenBound>(l.external_memory))
        throw PgCoefficientAbort{{core::FailureCode::Unsupported, core::Operation::PrepareModel}};
    if (l.maximum_input_characters < 2)
        PgCoefficientGeneration::exceeded();
    if (limits.cancellation.stop_requested())
        throw PgCoefficientAbort{{core::FailureCode::Cancelled, core::Operation::PrepareModel}};
    const auto n = static_cast<std::uint64_t>(graph.n), axis = 2 * n, cells = axis * axis;
    const auto nr = (n + 2) * (n + 2) + 1, nd = (2 * n + 1) * (2 * n + 1) + 3, pairs = 96 * n,
               nc = (n + 1) * (n + 1);
    PgCoefficientPlan result;
    result.cells = cells;
    result.raw_columns = nr + nd + pairs + nc;
    result.full_columns = result.raw_columns + 99 * n + 7;
    result.record_capacity = pg_mul(cells, pg_add(result.full_columns, pairs));
    result.map_capacity = 5 * (n + 2) * (n + 2) + 4 + 108 * axis;
    result.c_map_capacity = n * (n + 1);
    if (result.full_columns > limits.maximum_columns ||
        pg_add(result.map_capacity, result.c_map_capacity) > limits.maximum_map_entries)
        PgCoefficientGeneration::exceeded();
    PgCoefficientBudget b{limits};
    std::uint64_t input_heap = 0, max_input = 1;
    const auto inspect = [&](const ExactRational& value) {
        const auto bits = std::max(value.numerator_bits(), value.denominator_bits());
        if (bits > l.maximum_value_bits)
            PgCoefficientGeneration::exceeded();
        max_input = std::max(max_input, bits);
        const auto& v = ExactFactory::value(value).backend();
        b.opaque_input(v.denom(), value.denominator_bits());
        input_heap = pg_add(input_heap, pg_actual_external(value));
    };
    inspect(graph.x0);
    inspect(graph.width);
    inspect(graph.skew);
    for (const auto& term : graph.height.terms())
        inspect(term.coefficient);
    b.opaque_input(ExactFactory::value(graph.width).backend().num(), graph.width.numerator_bits());
    b.live_external = input_heap;
    for (const auto* value : {&graph.x0, &graph.width, &graph.skew})
        b.retain(b.input(*value), 1);
    for (const auto& term : graph.height.terms())
        b.retain(b.input(term.coefficient), 1);
    // Materialize cell geometry once; the two local-power caches remain owned
    // until verification ends. The graph snapshot uses explicit fresh clones.
    b.repeat = cells;
    const auto geometry = pg_budget_geometry(b, graph);
    for (const auto* value : {&geometry.s, &geometry.t, &geometry.f, &geometry.fx, &geometry.fy,
                              &geometry.U, &geometry.V, &geometry.xy})
        b.retain(*value, cells);
    b.retain(geometry.g, cells, true);
    for (const auto& value : geometry.s_power)
        b.retain(value, cells);
    for (const auto& value : geometry.t_power)
        b.retain(value, cells);
    // Fixed support: at most three quadratic and two linear bases per axis/cell.
    pg_budget_basis(b, n);
    auto quadratic = geometry.s;
    quadratic.terms = 3;
    quadratic.first = 2;
    quadratic.degree = 2;
    quadratic.denominator = pg_den_native(64);
    quadratic.l1 = 2;
    quadratic.storage = pg_add(pg_mul(2, b.value_bits(quadratic)), 1);
    auto transverse = quadratic;
    transverse.first = 0;
    transverse.second = 2;
    b.repeat = 9 * cells;
    const auto r_tensor = b.m(quadratic, transverse);
    b.retain(r_tensor, 9 * cells, true);
    b.repeat = cells;
    b.retain(b.clone(geometry.f), cells, true);
    b.repeat = 9 * cells;
    const auto d_tensor = b.m(geometry.g, b.m(quadratic, transverse));
    b.retain(d_tensor, 9 * cells, true);
    b.repeat = cells;
    b.retain(b.s(b.m(geometry.g, geometry.fx), b.number(1)), cells, true);
    b.retain(b.s(b.m(geometry.g, geometry.fy), b.number(1)), cells, true);
    b.retain(b.clone(geometry.g), cells, true);
    b.repeat = 4 * cells;
    const auto c_tensor = b.m(geometry.s, geometry.t);
    b.retain(c_tensor, 4 * cells, true);
    b.repeat = 4;
    pg_budget_legendre(b);
    // A pair exists only on and after its origin. Active and completed-carry
    // cells have different polynomial degrees and therefore different costs.
    std::array<PgBudgetPolynomial, 3> maximum_v{}, maximum_p{};
    PgBudgetPolynomial anchor_v_envelope, anchor_p_envelope;
    for (const auto axis_kind : {PgPairedAxis::S, PgPairedAxis::T}) {
        const auto derivative_axis =
            axis_kind == PgPairedAxis::S ? ExactAxis::First : ExactAxis::Second;
        b.repeat = 8 * axis;
        const auto active = pg_budget_moments(b, n, derivative_axis, false);
        const auto enabled = b.charging;
        b.charging = false;
        const auto completed = pg_budget_moments(b, n, derivative_axis, true);
        b.charging = enabled;
        b.repeat = 8 * axis * axis;
        const auto direct = pg_budget_pair(b, graph, geometry, active, axis_kind);
        b.repeat = 4 * axis * axis * (axis - 1);
        const auto tail = pg_budget_pair(b, graph, geometry, completed, axis_kind);
        for (std::size_t j = 0; j < 3; ++j) {
            b.retain(direct.v[j], 8 * axis * axis, true);
            b.retain(direct.p[j], 8 * axis * axis, true);
            b.retain(tail.v[j], 4 * axis * axis * (axis - 1), true);
            b.retain(tail.p[j], 4 * axis * axis * (axis - 1), true);
            maximum_v[j] = pg_budget_union(maximum_v[j], direct.v[j]);
            maximum_p[j] = pg_budget_union(maximum_p[j], direct.p[j]);
            if (axis_kind == PgPairedAxis::T) {
                anchor_v_envelope = pg_budget_union(anchor_v_envelope, direct.v[j]);
                anchor_p_envelope = pg_budget_union(anchor_p_envelope, direct.p[j]);
            }
        }
    }
    // All raw R/P traces are checked, including the structural-zero s family.
    // Coefficients with u removed are one-dimensional before multiplying ell.
    auto anchor_v = anchor_v_envelope, anchor_p = anchor_p_envelope;
    b.repeat = (nr + nd + 2 * pairs) * axis;
    anchor_v = b.restriction(anchor_v);
    anchor_p = b.restriction(anchor_p);
    b.repeat = 1;
    b.charge(1, 1, 1); // geometry unit constant
    b.charge(1, 1, 1); // anchor unit constant
    b.repeat = axis;
    (void)b.a(b.number(1), geometry.s);
    auto ell = geometry.s;
    ell.l1 = 0;
    b.repeat = (n + 3 + 48 * n) * cells;
    const auto ar = b.m(ell, anchor_v);
    b.retain(ar, (n + 3 + 48 * n) * cells, true);
    b.repeat = (2 * n + 4 + 48 * n) * cells;
    const auto ap = b.m(ell, anchor_p);
    b.retain(ap, (2 * n + 4 + 48 * n) * cells, true);
    b.repeat = 2 * (result.full_columns - result.raw_columns) * axis;
    (void)b.restriction(ar);
    (void)b.restriction(ap);
    b.repeat = nc * axis;
    (void)b.restriction(c_tensor);
    b.repeat = nc * cells;
    (void)b.m(ell, b.restriction(c_tensor));
    b.repeat = pg_mul(result.c_map_capacity, cells);
    (void)b.a(c_tensor, b.s(c_tensor, b.number(1)));
    // Native exact maps: every interval has sum_(k=0)^7(k+1)=36 potential
    // coefficients, copied to three component fields. Original zeros remain.
    PgBudgetPolynomial map_max;
    for (std::uint64_t degree = 0; degree < 8; ++degree) {
        b.repeat = 2 * axis;
        const auto curvature = pg_budget_global_curvature(b, n, degree);
        b.repeat = axis;
        const auto coefficients = pg_budget_decompose(b, curvature, degree);
        for (std::uint64_t j = 0; j <= degree; ++j) {
            b.retain(coefficients[j], 3 * axis);
            map_max = pg_budget_union(map_max, coefficients[j]);
        }
    }
    const auto output_rows = 5 * (n + 2) * (n + 2) + 4;
    auto coordinate_map = b.input(graph.x0);
    // Greville coefficients have denominator 2n. Include width and skew even
    // when x0 is zero: every nonzero coordinate is a scalar map entry.
    coordinate_map.terms = 1;
    coordinate_map.first = coordinate_map.second = coordinate_map.degree = 0;
    coordinate_map.l1 = pg_add(max_input, 4);
    const auto greville = pg_den_native(2 * n);
    coordinate_map.denominator = pg_den_join(
        coordinate_map.denominator,
        pg_den_join(greville,
                    pg_den_join(pg_den_product(b.input(graph.width).denominator, greville),
                                pg_den_product(b.input(graph.skew).denominator, greville))));
    coordinate_map.storage = b.value_bits(coordinate_map);
    b.retain(coordinate_map, output_rows);
    b.retain(b.number(1, n), result.c_map_capacity);
    map_max = pg_budget_union(map_max, coordinate_map);
    map_max = pg_budget_union(map_max, b.number(1));
    // Independent basic expected fields own x0/unit, affine X, and normal completions.
    b.repeat = 1;
    b.charge(1, 1, b.value_bits(b.input(graph.x0)));
    b.charge(1, 1, 1);
    b.repeat = cells;
    (void)b.a(b.a(b.input(graph.x0), b.s(geometry.s, b.input(graph.width))),
              b.s(geometry.t, b.input(graph.skew)));
    (void)b.s(b.m(geometry.g, geometry.fx), b.number(1));
    (void)b.s(b.m(geometry.g, geometry.fy), b.number(1));
    (void)b.clone(geometry.g);
    // Reconstruction consumes only entries in the sparse 33-field maps. Its
    // expected global curvature is built separately with carried moments.
    auto map_scalar = map_max;
    map_scalar.terms = 1;
    map_scalar.first = 0;
    map_scalar.second = 0;
    map_scalar.degree = 0;
    PgBudgetPolynomial worst;
    for (std::size_t j = 0; j < 3; ++j) {
        worst = pg_budget_union(worst, maximum_p[j]);
        worst = pg_budget_union(worst, maximum_v[j]);
    }
    worst.terms = std::min(pg_monomials(worst.degree),
                           pg_mul(pg_add(worst.first, 1), pg_add(worst.second, 1)));
    const auto weighted = b.s(worst, map_scalar);
    auto partial_sum = weighted;
    // Partial sums include every tensor-basis or interval-curvature entry.
    partial_sum.l1 = pg_add(partial_sum.l1, pg_log(std::max((n + 2) * (n + 2), 8 * axis)));
    partial_sum.storage = std::max(partial_sum.storage, b.value_bits(partial_sum));
    b.repeat = 3 * result.map_capacity * cells;
    (void)b.a(partial_sum, b.s(worst, map_scalar));
    for (std::uint64_t degree = 0; degree < 8; ++degree) {
        const auto saved = b.charging;
        b.charging = false;
        auto curvature = pg_budget_global_curvature(b, n, degree);
        b.charging = saved;
        std::array<PgBudgetPolynomial, 5> carry;
        for (std::uint64_t cs = 0; cs < axis; ++cs) {
            b.repeat = 1;
            std::array<PgBudgetPolynomial, 5> moments;
            for (std::size_t m = 0; m < 5; ++m) {
                const auto active = b.i(
                    b.s(b.m(geometry.s_power[m], curvature), b.number(1, 4 * n)), ExactAxis::First);
                b.charge(carry[m].terms, carry[m].terms, b.value_bits(carry[m]));
                moments[m] = b.a(active, carry[m]);
                pg_budget_evaluation(b, moments[m]);
                carry[m] = pg_budget_scalar(b, moments[m]);
            }
            b.repeat = axis;
            (void)pg_budget_pair(b, graph, geometry, moments, PgPairedAxis::S);
        }
    }
    // Distinct declared native/rational operations, not a free multiplier over
    // all polynomial calls. T can negate every input term before make().
    const auto trace_calls = (nr + nd + 2 * pairs) * axis +
                             2 * (result.full_columns - result.raw_columns) * axis + nc * axis;
    const auto nodes = 2 * axis * (6 * n + 8);
    auto native_calls = pg_add(9 + 4 * cells, pg_add(4 * nodes + 2 * axis * (3 * n + 7) + 8 * axis,
                                                     std::uint64_t{4} * 19 + 80 * axis));
    native_calls = pg_add(native_calls, pg_mul(8 * axis * axis * (axis + 1), 21));
    native_calls = pg_add(native_calls, pg_mul(trace_calls, 254));
    native_calls = pg_add(native_calls, result.c_map_capacity + 9 * (n + 2) * (n + 2) + (n + 2) +
                                            5 + 144 * axis + 48 * axis + 3);
    pg_budget_native_payload(b, native_calls,
                             std::max(b.value_bits(map_max), pg_mul(2, max_input)));
    // The bounded live work set is counted separately from published values:
    // J jets, six plane derivatives, two contractions, active operands,
    // expected v/p during verification, moment/power/Legendre rows, native T
    // negations and Greville/scalar locals. No result-wide scratch multiplier.
    const auto working = 3 * pg_monomials(12) + 6 * pg_monomials(11) + 2 * pg_monomials(14) +
                         2 * pg_monomials(20) + 3 * pg_monomials(12) + 3 * pg_monomials(20) + 65 +
                         5 + 15 + 32 + 253 + 12 + (12 * n + 10);
    b.live_external =
        pg_add(b.live_external,
               pg_mul(pg_mul(pg_mul(2, working), PgCoefficientBudget::capacity(b.maximum_storage)),
                      sizeof(boost::multiprecision::limb_type)));
    const auto graph_clone = pg_add(
        pg_add(sizeof(PgGraphData), pg_mul(cells, sizeof(PgCellIndex))),
        pg_add(pg_mul(3, sizeof(RationalStorage)), pg_payload_poly(graph.height.terms().size())));
    auto baseline =
        pg_add(sizeof(PgCoefficientData), pg_mul(result.full_columns, sizeof(PgCoefficientSource)));
    baseline = pg_add(baseline, pg_mul(result.record_capacity, sizeof(PgCoefficientRecord)));
    baseline = pg_add(baseline, pg_mul(cells, sizeof(PgPoly)));
    baseline = pg_add(baseline, pg_mul(result.map_capacity, sizeof(PgRequiredMapEntry)));
    baseline = pg_add(baseline, pg_mul(result.c_map_capacity, sizeof(PgAnchorCMapEntry)));
    baseline = pg_add(baseline, pg_mul(cells, sizeof(PgCoefficientCellCache)));
    baseline = pg_add(baseline, graph_clone);
    result.report.planned_owned_payload_upper_bytes = pg_add(baseline, b.payload);
    result.report.planned_column_capacity = result.full_columns;
    result.report.arithmetic.algebraic_work_bits =
        std::max(b.maximum_work, pg_mul(graph.height.terms().size(), pg_add(max_input, 1)));
    constexpr auto limb_bits = IntegerBackend::limb_bits;
    const auto work = result.report.arithmetic.algebraic_work_bits;
    result.report.arithmetic.planned_work_bits =
        work == 0
            ? 0
            : pg_mul(pg_add(pg_add(work / limb_bits, work % limb_bits != 0 ? 1 : 0), 2), limb_bits);
    result.report.arithmetic.planned_raw_terms =
        std::max(b.maximum_raw, static_cast<std::uint64_t>(graph.height.terms().size()));
    result.report.arithmetic.estimated_external_payload_bytes =
        pg_add(b.live_external, b.leaf_external);
    if (b.canonical_terms > l.maximum_terms ||
        result.report.arithmetic.planned_raw_terms > l.maximum_raw_terms ||
        b.first_degree > l.maximum_axis_degree || b.second_degree > l.maximum_axis_degree ||
        b.total_degree > l.maximum_total_degree || b.maximum_value > l.maximum_value_bits ||
        result.report.arithmetic.planned_work_bits > l.maximum_work_bits ||
        b.output_terms > limits.maximum_stored_terms)
        PgCoefficientGeneration::exceeded();
    const auto* permission = std::get_if<ExactEstimateAdmission>(&l.external_memory);
    if (permission == nullptr)
        throw PgCoefficientAbort{{core::FailureCode::InvalidInput, core::Operation::PrepareModel}};
    if (result.report.arithmetic.estimated_external_payload_bytes > permission->bytes ||
        result.report.planned_owned_payload_upper_bytes > owned.remaining_bytes())
        PgCoefficientGeneration::exceeded();
    return result;
}
} // namespace uvlas::structure::detail
