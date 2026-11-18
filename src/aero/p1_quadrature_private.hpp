#pragma once
#include "p1_gauss_private.hpp"
#include "p1_interval_private.hpp"
namespace uvlas::aero::detail::p1 {
struct QuadratureRule {
    std::array<p1::Interval, 12> nodes, weights;
    std::size_t count = 0;
};
template <std::size_t N>
void build_rule(const std::array<p1::GaussRoot, N>& source, QuadratureRule& out, p1::Work& work) {
    static_assert(N == 8 || N == 12);
    out.count = N;
    for (std::size_t i = 0; i < N; ++i) {
        const Interval u{{Integer{source[i].lower.data()}, -256},
                         {Integer{source[i].upper.data()}, -256}};
        out.nodes[i] = u;
        const Interval x = subtract(multiply(integer(2), u, work), integer(1), work);
        Interval previous = integer(1), current = x;
        for (unsigned k = 2; k <= N; ++k) {
            const Interval next =
                divide(subtract(multiply(multiply(integer(2 * k - 1), x, work), current, work),
                                multiply(integer(k - 1), previous, work), work),
                       integer(k), work);
            previous = current;
            current = next;
        }
        const Interval derivative =
            divide(multiply(integer(N), subtract(multiply(x, current, work), previous, work), work),
                   subtract(multiply(x, x, work), integer(1), work), work);
        out.weights[i] = divide(integer(1),
                                multiply(subtract(integer(1), multiply(x, x, work), work),
                                         multiply(derivative, derivative, work), work),
                                work);
    }
}
} // namespace uvlas::aero::detail::p1
