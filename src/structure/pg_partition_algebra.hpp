#pragma once
#include "pg_partition_plan.hpp"

#include <array>
#include <charconv>

namespace uvlas::structure::detail {
struct PgPartitionAlgebra {
    const PgPartitionLimits& caps;
    PgPartitionBudget budget;
    core::ResourcePlan& payload;
    ExactLimits exact;
    std::uint64_t numeric_used = 0, topology_used = 0, observed_bits = 0;
    PgPartitionAlgebra(const PgPartitionLimits& limits, PgPartitionBudget plan,
                       core::ResourcePlan& resource)
        : caps(limits), budget(plan), payload(resource), exact(limits.algebra) {
        exact.timing = ExactTimingMode::Disabled;
    }
    void numeric(std::uint64_t amount = 1) {
        partition_cancel(caps);
        numeric_used = partition_add(numeric_used, amount);
        if (numeric_used > budget.report.numeric_operations)
            partition_fail(core::FailureCode::ResourceExceeded);
    }
    void visit(std::uint64_t amount = 1) {
        partition_cancel(caps);
        topology_used = partition_add(topology_used, amount);
        if (topology_used > budget.report.topology_visits)
            partition_fail(core::FailureCode::ResourceExceeded);
    }
    ExactRational finish(core::Result<ExactRational> result) {
        auto x = partition_take(std::move(result));
        observed_bits = std::max(observed_bits, partition_bits(x));
        if (x.preparation().planned_work_bits > budget.report.planned_work_bits)
            partition_fail(core::FailureCode::ResourceExceeded);
        return x;
    }
    ExactRational clone(const ExactRational& x) {
        numeric();
        return finish(x.clone(exact, payload));
    }
    ExactRational add(const ExactRational& a, const ExactRational& b) {
        numeric();
        return finish(a.add(b, exact, payload));
    }
    ExactRational sub(const ExactRational& a, const ExactRational& b) {
        numeric();
        return finish(a.subtract(b, exact, payload));
    }
    ExactRational mul(const ExactRational& a, const ExactRational& b) {
        numeric();
        return finish(a.multiply(b, exact, payload));
    }
    ExactRational div(const ExactRational& a, const ExactRational& b) {
        numeric();
        return finish(a.divide(b, exact, payload));
    }
    int compare(const ExactRational& a, const ExactRational& b) {
        numeric();
        const auto result = partition_take(a.compare(b, exact, payload));
        return result < 0 ? -1 : result > 0 ? 1 : 0;
    }
    ExactRational number(core::Index a, core::Index b = 1) {
        numeric();
        std::array<char, 32> na{}, nb{};
        const auto va = std::to_chars(na.data(), na.data() + na.size(), a),
                   vb = std::to_chars(nb.data(), nb.data() + nb.size(), b);
        if (va.ec != std::errc{} || vb.ec != std::errc{})
            partition_fail(core::FailureCode::ResourceExceeded);
        return finish(ExactRational::from_integer_pair(
            {na.data(), static_cast<std::size_t>(va.ptr - na.data())},
            {nb.data(), static_cast<std::size_t>(vb.ptr - nb.data())}, exact, payload));
    }
    PgPartitionPoint point_clone(const PgPartitionPoint& x) {
        return {clone(x.s), clone(x.t)};
    }
    ExactRational determinant(const ExactRational& a, const ExactRational& b,
                              const ExactRational& c, const ExactRational& d) {
        const auto x = mul(a, d), y = mul(b, c);
        return sub(x, y);
    }
    ExactRational orientation(const PgPartitionPoint& a, const PgPartitionPoint& b,
                              const PgPartitionPoint& c) {
        const auto x = sub(b.s, a.s), y = sub(b.t, a.t), u = sub(c.s, a.s), v = sub(c.t, a.t);
        return determinant(x, y, u, v);
    }
    bool between(const ExactRational& a, const ExactRational& b, const ExactRational& x) {
        const auto xa = compare(x, a), xb = compare(x, b);
        return (xa >= 0 && xb <= 0) || (xa <= 0 && xb >= 0);
    }
    bool on_segment(const PgPartitionPoint& a, const PgPartitionPoint& b,
                    const PgPartitionPoint& x) {
        return orientation(a, b, x).sign() == 0 && between(a.s, b.s, x.s) && between(a.t, b.t, x.t);
    }
    int point_order(const PgPartitionPoint& a, const PgPartitionPoint& b) {
        const auto s = compare(a.s, b.s);
        return s != 0 ? s : compare(a.t, b.t);
    }
    ExactRational line_value(const PgPartitionLine& l, const PgPartitionPoint& p) {
        const auto s = mul(l.coefficients[0], p.s), t = mul(l.coefficients[1], p.t);
        const auto sum = add(s, t);
        return add(sum, l.coefficients[2]);
    }
    ExactRational publish(Rational value) {
        numeric();
        const auto bits = [](const IntegerBackend& x) {
            if (x.size() == 1 && x.limbs()[0] == 0)
                return std::uint64_t{0};
            return static_cast<std::uint64_t>((x.size() - 1) * IntegerBackend::limb_bits) +
                   static_cast<std::uint64_t>(IntegerBackend::limb_bits) -
                   static_cast<std::uint64_t>(std::countl_zero(x.limbs()[x.size() - 1]));
        };
        if (std::max(bits(value.backend().num()), bits(value.backend().denom())) >
            exact.maximum_value_bits)
            partition_fail(core::FailureCode::ResourceExceeded);
        return finish(ExactFactory::rational(std::move(value), exact, payload, {}));
    }
    Integer gcd(Integer a, Integer b) {
        if (a < 0)
            a = -a;
        if (b < 0)
            b = -b;
        // Fixed finite Euclidean temporaries fit the pinned division/GCD forecast.
        while (b != 0) {
            partition_cancel(caps);
            Integer r = a % b;
            a = std::move(b);
            b = std::move(r);
        }
        return a;
    }
    std::array<Integer, 3> homogeneous(const PgPartitionPoint& p) {
        const auto& s = ExactFactory::value(p.s).backend();
        const auto& t = ExactFactory::value(p.t).backend();
        numeric(3);
        return {Integer{s.num()} * Integer{t.denom()}, Integer{t.num()} * Integer{s.denom()},
                Integer{s.denom()} * Integer{t.denom()}};
    }
    std::array<Integer, 3> cross(const std::array<Integer, 3>& a, const std::array<Integer, 3>& b) {
        numeric(9);
        return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    }
    PgPartitionLine line(const PgPartitionPoint& a, const PgPartitionPoint& b) {
        auto values = cross(homogeneous(a), homogeneous(b));
        if (values[0] == 0 && values[1] == 0)
            partition_fail(core::FailureCode::GeometryInvalid);
        numeric(8);
        const auto divisor = gcd(gcd(values[0], values[1]), values[2]);
        for (auto& x : values)
            x /= divisor;
        const auto sign = values[0] != 0 ? values[0] < 0 : values[1] < 0;
        if (sign)
            for (auto& x : values)
                x = -x;
        PgPartitionLine output;
        for (std::size_t i = 0; i < 3; ++i)
            output.coefficients[i] = publish(Rational{values[i]});
        return output;
    }
    PgPartitionPoint intersection(const PgPartitionLine& a, const PgPartitionLine& b) {
        const auto ints = [](const PgPartitionLine& x) {
            return std::array<Integer, 3>{
                Integer{ExactFactory::value(x.coefficients[0]).backend().num()},
                Integer{ExactFactory::value(x.coefficients[1]).backend().num()},
                Integer{ExactFactory::value(x.coefficients[2]).backend().num()}};
        };
        const auto c = cross(ints(a), ints(b));
        if (c[2] == 0)
            partition_fail(core::FailureCode::InvalidInput);
        numeric(6);
        // Both coordinates come from these original support lines, never from
        // interpolating a rounded or repeatedly clipped generated vertex.
        Rational s{c[0]}, t{c[1]};
        s /= Rational{c[2]};
        t /= Rational{c[2]};
        return {publish(std::move(s)), publish(std::move(t))};
    }
};
template <class Compare>
void partition_sort(PgPartitionAlgebra& op, std::span<core::Index> order,
                    std::span<core::Index> merge, Compare compare) {
    for (std::size_t width = 1; width < order.size();) {
        for (std::size_t first = 0; first < order.size();) {
            const auto middle = first + std::min(width, order.size() - first);
            const auto last = middle + std::min(width, order.size() - middle);
            auto a = first, b = middle, out = first;
            while (a < middle || b < last) {
                op.visit();
                if (b == last || (a < middle && !compare(order[b], order[a])))
                    merge[out++] = order[a++];
                else
                    merge[out++] = order[b++];
            }
            first = last;
        }
        std::copy_n(merge.begin(), order.size(), order.begin());
        if (width > order.size() / 2)
            break;
        width *= 2;
    }
}
} // namespace uvlas::structure::detail
