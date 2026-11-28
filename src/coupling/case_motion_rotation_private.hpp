#pragma once

#include "aero/certified_trigonometry_private.hpp"

#include <array>
#include <cstdint>

namespace uvlas::coupling::detail::motion {

namespace m = aero::detail::p1;

struct ScalarJet {
    m::Interval value, first, second;
};
using VectorJet = std::array<ScalarJet, 3>;
using MatrixJet = std::array<ScalarJet, 9>;
using Vector = std::array<m::Interval, 3>;
using Matrix = std::array<m::Interval, 9>;

struct RigidJet {
    VectorJet centre;
    MatrixJet rotation;
};

namespace rotation_detail {

inline ScalarJet constant(const m::Interval& value) {
    return {value, m::integer(0), m::integer(0)};
}

inline ScalarJet add(const ScalarJet& a, const ScalarJet& b, m::Work& work) {
    return {m::add(a.value, b.value, work), m::add(a.first, b.first, work),
            m::add(a.second, b.second, work)};
}

inline ScalarJet subtract(const ScalarJet& a, const ScalarJet& b, m::Work& work) {
    return {m::subtract(a.value, b.value, work), m::subtract(a.first, b.first, work),
            m::subtract(a.second, b.second, work)};
}

inline ScalarJet scale(const ScalarJet& a, const m::Interval& factor, m::Work& work) {
    return {m::multiply(a.value, factor, work), m::multiply(a.first, factor, work),
            m::multiply(a.second, factor, work)};
}

inline ScalarJet multiply(const ScalarJet& a, const ScalarJet& b, m::Work& work) {
    const auto first =
        m::add(m::multiply(a.first, b.value, work), m::multiply(a.value, b.first, work), work);
    const auto second = m::add(
        m::add(m::multiply(a.second, b.value, work), m::multiply(a.value, b.second, work), work),
        m::scale(m::multiply(a.first, b.first, work), 1), work);
    return {m::multiply(a.value, b.value, work), first, second};
}

inline ScalarJet divide(const ScalarJet& a, const ScalarJet& b, m::Work& work) {
    // Differentiate a = q*b so each derivative uses the same certified divisor.
    const auto value = m::divide(a.value, b.value, work);
    const auto first =
        m::divide(m::subtract(a.first, m::multiply(value, b.first, work), work), b.value, work);
    const auto subtract_products = m::add(m::multiply(value, b.second, work),
                                          m::scale(m::multiply(first, b.first, work), 1), work);
    const auto second = m::divide(m::subtract(a.second, subtract_products, work), b.value, work);
    return {value, first, second};
}

inline ScalarJet with_time_derivatives(const ScalarJet& spatial, const ScalarJet& argument,
                                       m::Work& work) {
    return {
        spatial.value, m::multiply(spatial.first, argument.first, work),
        m::add(m::multiply(spatial.second, m::multiply(argument.first, argument.first, work), work),
               m::multiply(spatial.first, argument.second, work), work)};
}

/// a(x) = sin(sqrt(x))/sqrt(x) uses offset 1. b(x) = (1-cos(sqrt(x)))/x
/// uses offset 2. Both series and their first two x derivatives are alternating
/// with decreasing terms after term 32 on [0, 1]; term 32 bounds the tail.
inline ScalarJet small_coefficient(const ScalarJet& argument, std::int64_t offset, m::Work& work) {
    const auto x = argument.value;
    auto coefficient = m::divide(m::integer(1), m::integer(offset), work);
    auto series = constant(m::integer(0));
    auto power = m::integer(1);
    auto previous = m::integer(1);
    auto previous2 = m::integer(1);
    for (std::int64_t n = 0; n <= 32; ++n) {
        const auto value = m::multiply(coefficient, power, work);
        const auto first =
            m::multiply(m::multiply(coefficient, m::integer(n), work), previous, work);
        const auto second =
            m::multiply(m::multiply(coefficient, m::integer(n * (n - 1)), work), previous2, work);
        if (n == 32) {
            const auto enclose = [&](const m::Interval& sum, const m::Interval& omitted) {
                const m::Interval radius{{}, m::absolute(omitted).upper};
                return m::Interval{m::subtract(sum, radius, work).lower,
                                   m::add(sum, radius, work).upper};
            };
            return with_time_derivatives({enclose(series.value, value),
                                          enclose(series.first, first),
                                          enclose(series.second, second)},
                                         argument, work);
        }
        series = add(series, {value, first, second}, work);
        const auto divisor = 2 * n + offset + 1;
        coefficient =
            m::divide(m::negate(coefficient, work), m::integer(divisor * (divisor + 1)), work);
        previous2 = previous;
        previous = power;
        power = m::multiply(power, x, work);
    }
    throw m::ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
}

struct RodriguesCoefficients {
    ScalarJet a, b;
};

inline RodriguesCoefficients coefficients(const ScalarJet& squared, const m::Interval& pi,
                                          m::Work& work) {
    if (m::compare(squared.value.upper, m::integer(1).upper) <= 0)
        return {small_coefficient(squared, 1, work), small_coefficient(squared, 2, work)};

    const auto theta = m::square_root(squared.value, work);
    const auto twice_theta = m::scale(theta, 1);
    const auto rate = m::divide(squared.first, twice_theta, work);
    const auto rate_squared = m::multiply(rate, rate, work);
    const auto acceleration =
        m::divide(m::subtract(squared.second, m::scale(rate_squared, 1), work), twice_theta, work);
    const auto trig = m::sine_cosine(theta, m::PhaseUnits::Radians, pi, work);
    const ScalarJet sine{trig.sine, m::multiply(trig.cosine, rate, work),
                         m::subtract(m::multiply(trig.cosine, acceleration, work),
                                     m::multiply(trig.sine, rate_squared, work), work)};
    const ScalarJet one_minus_cosine{m::subtract(m::integer(1), trig.cosine, work),
                                     m::multiply(trig.sine, rate, work),
                                     m::add(m::multiply(trig.cosine, rate_squared, work),
                                            m::multiply(trig.sine, acceleration, work), work)};
    return {divide(sine, {theta, rate, acceleration}, work),
            divide(one_minus_cosine, squared, work)};
}

inline MatrixJet exponential(const VectorJet& phi, const m::Interval& pi, m::Work& work) {
    const auto zero = constant(m::integer(0));
    const auto one = constant(m::integer(1));
    const auto minus_one = m::integer(-1);
    auto squared = zero;
    for (const auto& component : phi)
        squared = add(squared, multiply(component, component, work), work);
    const auto [a, b] = coefficients(squared, pi, work);
    const MatrixJet skew{zero,
                         scale(phi[2], minus_one, work),
                         phi[1],
                         phi[2],
                         zero,
                         scale(phi[0], minus_one, work),
                         scale(phi[1], minus_one, work),
                         phi[0],
                         zero};
    MatrixJet result{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            const auto slot = 3 * row + column;
            auto skew_squared = multiply(phi[row], phi[column], work);
            if (row == column)
                skew_squared = subtract(skew_squared, squared, work);
            const auto increment =
                add(multiply(a, skew[slot], work), multiply(b, skew_squared, work), work);
            result[slot] = row == column ? add(one, increment, work) : increment;
        }
    }
    return result;
}

} // namespace rotation_detail

/// The origin frame has already been qualified. The two input jets use frame-local
/// metres and radians, with time derivatives per SI second.
inline RigidJet compose_local(const Vector& centre, const Matrix& frame, const VectorJet& offset,
                              const VectorJet& rotation, const m::Interval& pi, m::Work& work) {
    using namespace rotation_detail;
    const auto relative = exponential(rotation, pi, work);
    RigidJet result{};
    for (std::size_t row = 0; row < 3; ++row) {
        result.centre[row] = constant(centre[row]);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto& factor = frame[3 * row + axis];
            result.centre[row] = add(result.centre[row], scale(offset[axis], factor, work), work);
            for (std::size_t column = 0; column < 3; ++column) {
                auto& target = result.rotation[3 * row + column];
                target = add(target, scale(relative[3 * axis + column], factor, work), work);
            }
        }
    }
    return result;
}

/// Apply Q(t) Q0^T to a vector expressed in the original world frame.
inline VectorJet map_reference_vector(const Vector& reference, const Matrix& frame,
                                      const MatrixJet& rotation, VectorJet origin, m::Work& work) {
    using namespace rotation_detail;
    Vector local{};
    for (std::size_t axis = 0; axis < 3; ++axis)
        for (std::size_t world = 0; world < 3; ++world)
            local[axis] = m::add(
                local[axis], m::multiply(frame[3 * world + axis], reference[world], work), work);
    auto result = std::move(origin);
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t axis = 0; axis < 3; ++axis)
            result[row] =
                add(result[row], scale(rotation[3 * row + axis], local[axis], work), work);
    return result;
}

/// Map a fixed reference point, expressed in the original world frame.
inline VectorJet map_reference_point(const Vector& reference, const Vector& centre,
                                     const Matrix& frame, const RigidJet& rigid, m::Work& work) {
    Vector difference{};
    for (std::size_t axis = 0; axis < 3; ++axis)
        difference[axis] = m::subtract(reference[axis], centre[axis], work);
    return map_reference_vector(difference, frame, rigid.rotation, rigid.centre, work);
}

} // namespace uvlas::coupling::detail::motion
