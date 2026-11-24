#pragma once

#include "../aero/p1_interval_private.hpp"
#include "../structure/pg_rank_storage.hpp"

namespace uvlas::io::detail::nastran_directed {

namespace m = aero::detail::p1;
namespace g = structure::detail;
using g::partition_fail;
using Rational = structure::ExactRational;
using Code = core::FailureCode;
using Position = std::array<m::Interval, 3>;
using FixedInteger = aero::detail::ExactInteger;

[[nodiscard]] inline m::Interval norm_squared(const Position&, m::Work&);

[[nodiscard]] inline g::RankFraction<FixedInteger> fixed_fraction(const Rational& source,
                                                                  m::Work& work) {
    work.step();
    constexpr auto width = static_cast<std::uint64_t>(std::numeric_limits<g::RankWide>::digits);
    if (std::max(source.numerator_bits(), source.denominator_bits()) > width) {
        partition_fail(Code::ResourceExceeded);
    }
    return g::rank_fraction<FixedInteger>(source);
}

[[nodiscard]] inline m::Interval ratio(const FixedInteger& numerator,
                                       const FixedInteger& denominator, m::Work& work) {
    return m::divide(m::exact(numerator, 0), m::exact(denominator, 0), work);
}

[[nodiscard]] inline Rational retain(const m::Bound& input, const structure::ExactLimits& exact,
                                     core::ResourcePlan& plan, m::Work& work) {
    work.step();
    if (input.mantissa == 0) {
        return {};
    }
    constexpr auto width = std::numeric_limits<g::RankWide>::digits;
    if (input.power < -width + 1 ||
        static_cast<std::int64_t>(m::bits(input.mantissa)) + input.power > width) {
        partition_fail(Code::ResourceExceeded);
    }
    const g::RankWide units{input.mantissa};
    if (input.power < 0) {
        return g::rank_materialize_dyadic(units, static_cast<std::uint32_t>(-input.power), exact,
                                          plan);
    }
    auto shifted = g::rank_abs(units) << static_cast<unsigned>(input.power);
    if (units < 0) {
        shifted = -shifted;
    }
    return g::rank_exact(g::RankFraction<g::RankWide>{std::move(shifted), g::RankWide{1}}, exact,
                         plan);
}

[[nodiscard]] inline Position subtract(const Position& first, const Position& second,
                                       m::Work& work) {
    return {m::subtract(first[0], second[0], work), m::subtract(first[1], second[1], work),
            m::subtract(first[2], second[2], work)};
}

[[nodiscard]] inline Position cross(const Position& first, const Position& second, m::Work& work) {
    Position output;
    for (std::size_t i = 0; i < output.size(); ++i) {
        const auto j = (i + 1) % 3, k = (i + 2) % 3;
        output[i] = m::subtract(m::multiply(first[j], second[k], work),
                                m::multiply(first[k], second[j], work), work);
    }
    return output;
}

[[nodiscard]] inline Position normalized(const Position& direction, m::Work& work) {
    const auto squared = norm_squared(direction, work);
    if (squared.upper.mantissa == 0) {
        partition_fail(Code::GeometryInvalid);
    }
    if (squared.lower.mantissa <= 0) {
        partition_fail(Code::PrecisionUnresolved);
    }
    const auto length = m::square_root(squared, work);
    return {m::divide(direction[0], length, work), m::divide(direction[1], length, work),
            m::divide(direction[2], length, work)};
}

[[nodiscard]] inline m::Interval number(const Rational& value, m::Work& work) {
    const auto value_fraction = fixed_fraction(value, work);
    return ratio(value_fraction.numerator, value_fraction.denominator, work);
}

[[nodiscard]] inline m::Interval norm_squared(const Position& position, m::Work& work) {
    m::Interval squared;
    for (const auto& item : position) {
        squared = m::add(squared, m::multiply(item, item, work), work);
    }
    return squared;
}

} // namespace uvlas::io::detail::nastran_directed
