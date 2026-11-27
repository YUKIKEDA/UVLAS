#pragma once

#include "p1_interval_private.hpp"

#include <span>

namespace uvlas::aero::detail::p1 {
inline std::uint64_t checked_sum(std::uint64_t left, std::uint64_t right) {
    if (right > std::numeric_limits<std::uint64_t>::max() - left)
        throw ArithmeticFailure{core::FailureCode::ResourceExceeded};
    return left + right;
}
inline std::uint64_t checked_product(std::uint64_t left, std::uint64_t right) {
    if (left != 0 && right > std::numeric_limits<std::uint64_t>::max() / left)
        throw ArithmeticFailure{core::FailureCode::ResourceExceeded};
    return left * right;
}
enum class Opcode : std::uint8_t {
    Input,
    IntegerConstant,
    Pi,
    Add,
    Subtract,
    Negate,
    Multiply,
    Divide,
    Sqrt,
    Log1p,
    Atan2,
    ValueWithDerivative
};
inline constexpr std::uint64_t root_steps = UINT64_C(2) * 257;
inline constexpr std::uint64_t log_series_steps = 1 + 128 * 3 + 1 + 4 + 2;
inline constexpr std::uint64_t log_endpoint_steps = 3 + log_series_steps + 2;
inline constexpr std::uint64_t log_one_plus_steps = 1 + 2 * log_endpoint_steps;
inline constexpr std::uint64_t atan_series_steps = 1 + 128 * 3 + 1 + 2;
inline constexpr std::uint64_t atan_half_steps = 4 + root_steps + atan_series_steps + 1;
inline constexpr std::uint64_t atan_point_steps = 1 + atan_half_steps + 2;
inline constexpr std::uint64_t atan_two_steps = 4 * atan_point_steps;
static_assert(log_one_plus_steps + 14 == 809 && atan_two_steps + 55 == 3695);
struct Node {
    Opcode op = Opcode::IntegerConstant;
    std::size_t left = 0, right = 0;
    std::int64_t value = 0;
};
struct Program;
struct Expression {
    Program* program;
    std::size_t index;
};
/// The first identical pass has an empty span and performs no allocation.
struct Program {
    std::span<Node> nodes;
    std::size_t count = 0;
    std::uint64_t work = 0;
    Expression emit(Opcode op, std::size_t a = 0, std::size_t b = 0, std::int64_t v = 0) {
        // Counts are scalar interval operations and individual root candidates.
        // Smooth-gap substitution only copies already calculated components.
        constexpr std::array<std::uint64_t, 12> costs{0,
                                                      0,
                                                      0,
                                                      14,
                                                      14,
                                                      14,
                                                      40,
                                                      40,
                                                      root_steps + 14,
                                                      log_one_plus_steps + 14,
                                                      atan_two_steps + 55,
                                                      0};
        const auto cost = costs[static_cast<std::size_t>(op)];
        if (work > std::numeric_limits<std::uint64_t>::max() - cost ||
            count == std::numeric_limits<std::size_t>::max())
            throw ArithmeticFailure{core::FailureCode::ResourceExceeded};
        if (!nodes.empty()) {
            if (count >= nodes.size())
                throw ArithmeticFailure{core::FailureCode::ResourceExceeded};
            nodes[count] = {op, a, b, v};
        }
        const std::size_t index = count++;
        work += cost;
        return {this, index};
    }
    Expression input(std::size_t n) {
        return emit(Opcode::Input, n);
    }
    Expression constant(std::int64_t n) {
        return emit(Opcode::IntegerConstant, 0, 0, n);
    }
    Expression pi() {
        return emit(Opcode::Pi);
    }
};
inline Expression operator+(Expression a, Expression b) {
    return a.program->emit(Opcode::Add, a.index, b.index);
}
inline Expression operator-(Expression a, Expression b) {
    return a.program->emit(Opcode::Subtract, a.index, b.index);
}
inline Expression operator-(Expression a) {
    return a.program->emit(Opcode::Negate, a.index);
}
inline Expression operator*(Expression a, Expression b) {
    return a.program->emit(Opcode::Multiply, a.index, b.index);
}
inline Expression operator/(Expression a, Expression b) {
    return a.program->emit(Opcode::Divide, a.index, b.index);
}
inline Expression operator+(Expression a, std::int64_t b) {
    return a + a.program->constant(b);
}
inline Expression operator-(Expression a, std::int64_t b) {
    return a - a.program->constant(b);
}
inline Expression operator*(Expression a, std::int64_t b) {
    return a * a.program->constant(b);
}
inline Expression operator/(Expression a, std::int64_t b) {
    return a / a.program->constant(b);
}
inline Expression sqrt(Expression a) {
    return a.program->emit(Opcode::Sqrt, a.index);
}
inline Expression log1p(Expression a) {
    return a.program->emit(Opcode::Log1p, a.index);
}
inline Expression atan2(Expression a, Expression b) {
    return a.program->emit(Opcode::Atan2, a.index, b.index);
}
using Vector = std::array<Expression, 3>;
inline Vector operator+(const Vector& a, const Vector& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
inline Vector operator-(const Vector& a, const Vector& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline Vector operator*(const Vector& a, Expression b) {
    return {a[0] * b, a[1] * b, a[2] * b};
}
inline Vector operator/(const Vector& a, Expression b) {
    return {a[0] / b, a[1] / b, a[2] / b};
}
inline Expression dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vector cross(const Vector& a, const Vector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline Expression norm(const Vector& a) {
    return sqrt(dot(a, a));
}
struct Options {
    std::array<std::int8_t, 3> start_sign{}, end_sign{};
    bool exterior = false;
    std::array<bool, 3> boundary{true, true, true};
};
struct Output {
    std::array<std::size_t, 9> nodes{};
};
struct Inputs {
    Vector point;
    std::array<Vector, 3> vertices;
    Expression radius;
};
inline Inputs inputs(Program& p) {
    return {{p.input(0), p.input(1), p.input(2)},
            {{{p.input(3), p.input(4), p.input(5)},
              {p.input(6), p.input(7), p.input(8)},
              {p.input(9), p.input(10), p.input(11)}}},
            p.input(12)};
}
inline Expression contour_integral(Expression length, Expression s, Expression z, Expression ra,
                                   Expression rb, Expression k2, bool same) {
    if (same)
        return length * (s + z) / (ra * rb * (s * rb + z * ra));
    return (s / ra - z / rb) / k2;
}
/// Shared edge geometry makes every P1 column consume the identical RM kernel.
inline Output near_program(Program& p, const Options& options) {
    const auto in = inputs(p);
    const auto e = in.vertices[1] - in.vertices[0], f = in.vertices[2] - in.vertices[0];
    const auto c = cross(e, f);
    const auto jac = norm(c);
    const auto n = c * (p.constant(1) / jac);
    const auto d = dot(in.point - in.vertices[0], n);
    const auto ee = dot(e, e), ef = dot(e, f), ff = dot(f, f), determinant = ee * ff - ef * ef;
    std::array<Vector, 3> line{}, tangent{};
    auto sheet = Vector{p.constant(0), p.constant(0), p.constant(0)};
    std::array<Expression, 3> ra{p.constant(0), p.constant(0), p.constant(0)};
    std::array<Expression, 3> rb = ra;
    for (std::size_t edge = 0; edge < 3; ++edge) {
        const auto& a = in.vertices[edge];
        const auto& b = in.vertices[(edge + 1) % 3];
        const auto length = norm(b - a);
        tangent[edge] = (b - a) * (p.constant(1) / length);
        const auto r = in.point - a;
        const auto s = dot(r, tangent[edge]), z = s - length;
        const auto perpendicular = r - tangent[edge] * s;
        const auto k2 = dot(perpendicular, perpendicular) + in.radius * in.radius;
        ra[edge] = sqrt(s * s + k2);
        rb[edge] = sqrt(z * z + k2);
        if (options.boundary[edge]) {
            const bool same = options.start_sign[edge] * options.end_sign[edge] > 0;
            const auto integral = contour_integral(length, s, z, ra[edge], rb[edge], k2, same);
            // Negative contour reverses the positive geometric boundary.
            line[edge] = cross(tangent[edge], perpendicular) * (-integral);
        }
        // The absolute-value sum has an exact zero middle branch.
        const auto gapbase = options.start_sign[edge] < 0
                                 ? -(s * 2)
                                 : (options.end_sign[edge] > 0 ? z * 2 : p.constant(0));
        const auto abs_s = options.start_sign[edge] < 0 ? -s : s;
        const auto abs_z = options.end_sign[edge] > 0 ? z : -z;
        const auto stable_gap = gapbase + k2 / (ra[edge] + abs_s) + k2 / (rb[edge] + abs_z);
        // A distinct opcode substitutes the smooth differential of ra+rb-L below.
        const auto smooth_gap = ra[edge] + rb[edge] - length;
        const auto gap = p.emit(Opcode::ValueWithDerivative, stable_gap.index, smooth_gap.index);
        const auto logarithm = log1p(length * 2 / gap);
        // stable_gap and smooth_gap are identical physical functions; the
        // edge-sign branches only avoid cancellation in the value.

        sheet = sheet + cross(tangent[edge], n) * logarithm;
    }
    Expression normal = p.constant(0);
    if (options.exterior) {
        const auto r0 = in.point - in.vertices[0], r1 = in.point - in.vertices[1],
                   r2 = in.point - in.vertices[2];
        const auto l0 = norm(r0), l1 = norm(r1), l2 = norm(r2);
        const auto delta = l0 * l1 * l2 + dot(r0, r1) * l2 + dot(r1, r2) * l0 + dot(r2, r0) * l1;
        normal = d * jac * 2 / delta;
    } else {
        const auto h = sqrt(d * d + in.radius * in.radius);
        const auto virtual_point = in.point - n * d + n * h;
        const auto r0 = virtual_point - in.vertices[0], r1 = virtual_point - in.vertices[1],
                   r2 = virtual_point - in.vertices[2];
        const auto l0 = norm(r0), l1 = norm(r1), l2 = norm(r2);
        const auto delta = l0 * l1 * l2 + dot(r0, r1) * l2 + dot(r1, r2) * l0 + dot(r2, r0) * l1;
        normal = d / h * atan2(h * jac, delta) * 2;
    }
    sheet = sheet + n * normal;
    const auto factor = p.constant(1) / (p.pi() * 4);
    Output output;
    for (std::size_t k = 0; k < 3; ++k) {
        const auto m0 = p.constant(k == 0 ? 1 : 0), m1 = p.constant(k == 1 ? 1 : 0),
                   m2 = p.constant(k == 2 ? 1 : 0);
        const auto b1 = m1 - m0, b2 = m2 - m0;
        const auto gradient =
            e * ((ff * b1 - ef * b2) / determinant) + f * ((ee * b2 - ef * b1) / determinant);
        auto velocity = cross(cross(n, gradient), sheet);
        const std::array<Expression, 3> mu{m0, m1, m2};
        for (std::size_t edge = 0; edge < 3; ++edge) {
            if (!options.boundary[edge])
                continue;
            const auto strength =
                (mu[edge] * rb[edge] + mu[(edge + 1) % 3] * ra[edge]) / (ra[edge] + rb[edge]);
            velocity = velocity + line[edge] * strength;
        }
        for (std::size_t axis = 0; axis < 3; ++axis)
            output.nodes[3 * k + axis] = (velocity[axis] * factor).index;
    }
    return output;
}
// One unit-circulation column per facet edge. Bound rings follow the negative
// geometric boundary, the same physical convention as the P1 sheet.
inline Output bound_program(Program& p, const Options& options) {
    const auto in = inputs(p);
    const auto zero = p.constant(0);
    const auto factor = p.constant(-1) / (p.pi() * 4);
    const auto radius_squared = in.radius * in.radius;
    std::array<Vector, 3> relative{};
    std::array<Expression, 3> endpoint{};
    for (std::size_t vertex = 0; vertex < 3; ++vertex) {
        if (options.boundary[vertex] || options.boundary[(vertex + 2) % 3]) {
            relative[vertex] = in.point - in.vertices[vertex];
            endpoint[vertex] = sqrt(dot(relative[vertex], relative[vertex]) + radius_squared);
        }
    }
    Output out;
    for (std::size_t edge = 0; edge < 3; ++edge) {
        Vector velocity{zero, zero, zero};
        if (options.boundary[edge]) {
            const auto next = (edge + 1) % 3;
            const auto delta = in.vertices[next] - in.vertices[edge];
            const auto length_squared = dot(delta, delta);
            const auto from = dot(relative[edge], delta);
            const auto to = from - length_squared;
            const auto perpendicular = cross(delta, relative[edge]);
            const auto ra = endpoint[edge], rb = endpoint[next];
            const bool same = options.start_sign[edge] * options.end_sign[edge] > 0;
            const auto contour = [&] {
                if (same)
                    return (from + to) / (ra * rb * (from * rb + to * ra));
                const auto cross_squared = dot(perpendicular, perpendicular);
                return (from / ra - to / rb) / (cross_squared + radius_squared * length_squared);
            }();
            velocity = perpendicular * (factor * contour);
        }
        for (std::size_t axis = 0; axis < 3; ++axis)
            out.nodes[3 * edge + axis] = velocity[axis].index;
    }
    return out;
}
struct FarOutput {
    Output output;
    std::size_t sample_begin{};
    std::uint64_t setup_work{};
};

inline FarOutput far_program(Program& p) {
    const auto in = inputs(p);
    const auto one = p.constant(1);
    const auto three = p.constant(3);
    const auto cross_product =
        cross(in.vertices[1] - in.vertices[0], in.vertices[2] - in.vertices[0]);
    const auto jac = norm(cross_product);
    const auto n = cross_product * (one / jac);
    const auto radius_squared = in.radius * in.radius;
    const auto twice_radius_squared = radius_squared * 2;
    const auto four_pi = p.pi() * 4;
    // These jets depend only on the facet and query, not on the Gauss point.
    const auto sample_begin = p.count;
    const auto setup_work = p.work;
    const auto u = p.input(13), v = p.input(14), weight = p.input(15);
    const std::array<Expression, 3> bary{(one - u) * (one - v), u, (one - u) * v};
    const auto point =
        in.vertices[0] * bary[0] + in.vertices[1] * bary[1] + in.vertices[2] * bary[2];
    const auto r = in.point - point;
    const auto r2 = dot(r, r), den = r2 + radius_squared, rn = dot(r, n);
    const auto kernel = (n * (r2 - twice_radius_squared) - r * (rn * three)) *
                        (jac * weight * (one - u) / (den * den * sqrt(den))) / four_pi;
    Output out;
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t axis = 0; axis < 3; ++axis)
            out.nodes[3 * k + axis] = (kernel[axis] * bary[k]).index;
    return {out, sample_begin, setup_work};
}
struct Jet {
    std::array<Interval, 14> components;
};
inline void evaluate_node(const Node& node, const std::span<Jet> values, std::size_t index,
                          const std::array<Interval, 16>& inputs, int seed_power,
                          const Interval& pi, const Interval& ln2, Work& work) {
    auto& out = values[index].components;
    if (node.op == Opcode::Input || node.op == Opcode::IntegerConstant || node.op == Opcode::Pi) {
        out.fill(integer(0));
        out[0] = node.op == Opcode::Input ? inputs[node.left]
                                          : (node.op == Opcode::Pi ? pi : integer(node.value));
        if (node.op == Opcode::Input && node.left < 13)
            out[node.left + 1] = scale(integer(1), seed_power);
        return;
    }
    const auto& a = values[node.left].components;
    const auto& b = values[node.right].components;
    switch (node.op) {
    case Opcode::Add:
    case Opcode::Subtract:
    case Opcode::Negate:
        for (std::size_t j = 0; j < 14; ++j)
            out[j] = node.op == Opcode::Add
                         ? add(a[j], b[j], work)
                         : (node.op == Opcode::Subtract ? subtract(a[j], b[j], work)
                                                        : negate(a[j], work));
        break;
    case Opcode::Multiply:
        out[0] = multiply(a[0], b[0], work);
        for (std::size_t j = 1; j < 14; ++j)
            out[j] = add(multiply(a[j], b[0], work), multiply(a[0], b[j], work), work);
        break;
    case Opcode::Divide:
        out[0] = divide(a[0], b[0], work);
        for (std::size_t j = 1; j < 14; ++j)
            out[j] = divide(subtract(a[j], multiply(out[0], b[j], work), work), b[0], work);
        break;
    case Opcode::Sqrt: {
        out[0] = square_root(a[0], work);
        const auto twice = multiply(integer(2), out[0], work);
        for (std::size_t j = 1; j < 14; ++j)
            out[j] = divide(a[j], twice, work);
        break;
    }
    case Opcode::Log1p: {
        out[0] = log_one_plus(a[0], ln2, work);
        const auto denominator = add(integer(1), a[0], work);
        for (std::size_t j = 1; j < 14; ++j)
            out[j] = divide(a[j], denominator, work);
        break;
    }
    case Opcode::ValueWithDerivative:
        out = a;
        for (std::size_t j = 1; j < 14; ++j)
            out[j] = b[j];
        break;
    case Opcode::Atan2: {
        out[0] = atan_two(a[0], b[0], pi, work);
        const auto denominator = add(multiply(a[0], a[0], work), multiply(b[0], b[0], work), work);
        for (std::size_t j = 1; j < 14; ++j)
            out[j] = divide(subtract(multiply(b[0], a[j], work), multiply(a[0], b[j], work), work),
                            denominator, work);
        break;
    }
    default:
        throw ArithmeticFailure{core::FailureCode::InvalidInput};
    }
}
} // namespace uvlas::aero::detail::p1
