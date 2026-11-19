#pragma once

#include "si_interval.hpp"

#include <span>

namespace uvlas::io::detail::si {
enum class Opcode : std::uint8_t { Input, IntegerValue, Add, Subtract, Multiply, Divide };
struct Instruction {
    Opcode op = Opcode::IntegerValue;
    std::size_t left = 0, right = 0;
    std::int64_t value = 0;
};
struct Symbol {
    std::size_t index;
};
struct Builder {
    std::span<Instruction> nodes;
    std::size_t count = 0;
    std::uint64_t work = 0;
    constexpr Symbol emit(Opcode op, std::size_t left = 0, std::size_t right = 0,
                          std::int64_t value = 0) {
        if (!nodes.empty())
            nodes[count] = {op, left, right, value};
        // One dispatch visit; arithmetic opcodes also charge the interval primitive.
        work += op == Opcode::Input || op == Opcode::IntegerValue ? 1 : 2;
        return {count++};
    }
    constexpr Symbol input(std::size_t index) {
        return emit(Opcode::Input, index);
    }
    constexpr Symbol integer(std::int64_t value) {
        return emit(Opcode::IntegerValue, 0, 0, value);
    }
    constexpr Symbol add(Symbol a, Symbol b) {
        return emit(Opcode::Add, a.index, b.index);
    }
    constexpr Symbol subtract(Symbol a, Symbol b) {
        return emit(Opcode::Subtract, a.index, b.index);
    }
    constexpr Symbol multiply(Symbol a, Symbol b) {
        return emit(Opcode::Multiply, a.index, b.index);
    }
    constexpr Symbol divide(Symbol a, Symbol b) {
        return emit(Opcode::Divide, a.index, b.index);
    }
};
enum class Body : std::uint8_t {
    Digit,
    MultiplyOnly,
    SquareOnly,
    MultiplyAndSquare,
    AtanAdd,
    AtanSubtract
};
constexpr std::array<std::size_t, 2> definition(Body body, Builder& b) {
    const auto x = b.input(0), y = b.input(1);
    if (body == Body::Digit)
        return {b.add(b.multiply(x, b.integer(10)), y).index, 0};
    if (body == Body::MultiplyOnly)
        return {b.multiply(x, y).index, y.index};
    if (body == Body::SquareOnly)
        return {x.index, b.multiply(y, y).index};
    if (body == Body::MultiplyAndSquare)
        return {b.multiply(x, y).index, b.multiply(y, y).index};
    const auto square = b.input(2), denominator = b.input(3);
    const auto term = b.divide(x, denominator);
    const auto sum = body == Body::AtanAdd ? b.add(y, term) : b.subtract(y, term);
    return {b.multiply(x, square).index, sum.index};
}
constexpr std::size_t node_count(Body body) {
    Builder b{};
    (void)definition(body, b);
    return b.count;
}
template <std::size_t N>
struct StaticProgram {
    std::array<Instruction, N> nodes{};
    std::array<std::size_t, 2> outputs{};
    std::uint64_t work = 0;
};
template <Body B>
constexpr auto make_program() {
    StaticProgram<node_count(B)> result;
    Builder builder{result.nodes};
    result.outputs = definition(B, builder);
    result.work = builder.work;
    return result;
}
inline constexpr auto digit_program = make_program<Body::Digit>();
inline constexpr auto multiply_program = make_program<Body::MultiplyOnly>();
inline constexpr auto square_program = make_program<Body::SquareOnly>();
inline constexpr auto both_program = make_program<Body::MultiplyAndSquare>();
inline constexpr auto atan_add_program = make_program<Body::AtanAdd>();
inline constexpr auto atan_subtract_program = make_program<Body::AtanSubtract>();
inline constexpr std::size_t program_slots = std::max(
    {digit_program.nodes.size(), multiply_program.nodes.size(), square_program.nodes.size(),
     both_program.nodes.size(), atan_add_program.nodes.size(), atan_subtract_program.nodes.size()});
struct ProgramScratch {
    std::array<Interval, program_slots> values;
};
template <std::size_t N>
std::array<Interval, 2> execute(const StaticProgram<N>& program, std::span<const Interval> inputs,
                                ProgramScratch& scratch, Work& work) {
    for (std::size_t n = 0; n < N; ++n) {
        work.step();
        const auto& op = program.nodes[n];
        auto& output = scratch.values[n];
        switch (op.op) {
        case Opcode::Input:
            output = inputs[op.left];
            break;
        case Opcode::IntegerValue:
            output = integer(op.value);
            break;
        case Opcode::Add:
            output = add(scratch.values[op.left], scratch.values[op.right], work);
            break;
        case Opcode::Subtract:
            output = subtract(scratch.values[op.left], scratch.values[op.right], work);
            break;
        case Opcode::Multiply:
            output = multiply(scratch.values[op.left], scratch.values[op.right], work);
            break;
        case Opcode::Divide:
            output = divide(scratch.values[op.left], scratch.values[op.right], work);
            break;
        }
    }
    return {scratch.values[program.outputs[0]], scratch.values[program.outputs[1]]};
}
constexpr std::uint64_t power_work_bound(std::uint64_t magnitude) {
    std::uint64_t result = 0;
    while (magnitude) {
        const bool multiply = (magnitude & 1U) != 0;
        magnitude >>= 1;
        result += 1;
        if (multiply && magnitude)
            result += both_program.work;
        else if (multiply)
            result += multiply_program.work;
        else if (magnitude)
            result += square_program.work;
    }
    return result;
}
inline Interval programmed_power(Interval base, std::int64_t exponent, ProgramScratch& scratch,
                                 Work& w) {
    auto result = integer(1);
    auto remaining = magnitude(exponent);
    if (exponent < 0)
        base = divide(integer(1), base, w);
    while (remaining) {
        w.step();
        const bool multiply_bit = (remaining & 1U) != 0;
        remaining >>= 1;
        const std::array inputs{result, base};
        std::array<Interval, 2> values{result, base};
        if (multiply_bit && remaining)
            values = execute(both_program, inputs, scratch, w);
        else if (multiply_bit)
            values = execute(multiply_program, inputs, scratch, w);
        else if (remaining)
            values = execute(square_program, inputs, scratch, w);
        result = std::move(values[0]);
        base = std::move(values[1]);
    }
    return result;
}
inline Interval programmed_decimal(SiDecimalView input, ProgramScratch& scratch, Work& w) {
    auto result = integer(0);
    for (const char c : input.digits) {
        w.step();
        const std::array values{result, integer(c - '0')};
        result = execute(digit_program, values, scratch, w)[0];
    }
    if (!input.sign)
        return integer(0);
    if (input.sign < 0)
        result = {negate(result.upper), negate(result.lower)};
    return multiply(result, programmed_power(integer(10), input.exponent10, scratch, w), w);
}
inline constexpr std::uint64_t pi_work =
    2 * (2 + 64 * atan_add_program.work + 64 * atan_subtract_program.work + 128 + 2) + 3;
inline Interval programmed_pi(ProgramScratch& scratch, Work& w) {
    const auto atan = [&](std::int64_t denominator) {
        const auto x = divide(integer(1), integer(denominator), w), square = multiply(x, x, w);
        auto term = x, sum = integer(0);
        for (std::int64_t j = 0; j < 128; ++j) {
            w.step();
            const std::array inputs{term, sum, square, integer(2 * j + 1)};
            auto values = (j & 1) ? execute(atan_subtract_program, inputs, scratch, w)
                                  : execute(atan_add_program, inputs, scratch, w);
            term = std::move(values[0]);
            sum = std::move(values[1]);
        }
        const auto remainder = divide(term, integer(257), w);
        return Interval{sum.lower, add(sum, remainder, w).upper};
    };
    const auto a = atan(5), b = atan(239);
    return subtract(multiply(integer(16), a, w), multiply(integer(4), b, w), w);
}
} // namespace uvlas::io::detail::si
