#pragma once
#include "pg_partition_detail.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <type_traits>

namespace uvlas::structure::detail {
inline constexpr auto partition_operation = core::Operation::PrepareModel;
struct PartitionAbort {
    core::Failure failure;
};
[[noreturn]] inline void partition_fail(core::FailureCode code) {
    throw PartitionAbort{{code, partition_operation}};
}
inline std::uint64_t partition_add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        partition_fail(core::FailureCode::ResourceExceeded);
    return a + b;
}
inline std::uint64_t partition_mul(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        partition_fail(core::FailureCode::ResourceExceeded);
    return a * b;
}
inline core::Index partition_index(std::uint64_t n) {
    if (n > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        partition_fail(core::FailureCode::ResourceExceeded);
    return static_cast<core::Index>(n);
}
template <class T>
T partition_take(core::Result<T> value) {
    if (!value)
        throw PartitionAbort{value.error()};
    return std::move(*value);
}
inline void partition_take(core::Result<void> value) {
    if (!value)
        throw PartitionAbort{value.error()};
}
inline void partition_cancel(const PgPartitionLimits& limits) {
    if (limits.cancellation.stop_requested())
        partition_fail(core::FailureCode::Cancelled);
}
inline std::uint64_t partition_log(std::uint64_t n) noexcept {
    return n <= 1 ? 0 : static_cast<std::uint64_t>(std::bit_width(n - 1));
}
inline std::uint64_t partition_bits(const ExactRational& x) noexcept {
    return std::max(x.numerator_bits(), x.denominator_bits());
}
inline std::uint64_t partition_retained(const ExactRational& x) {
    const auto& b = ExactFactory::value(x).backend();
    const auto count = [](const IntegerBackend& v) {
        return v.capacity() > IntegerBackend::internal_limb_count ? v.capacity() : 0;
    };
    return partition_mul(partition_add(count(b.num()), count(b.denom())),
                         sizeof(boost::multiprecision::limb_type));
}
inline bool partition_equal_point(const PgPartitionPoint& a, const PgPartitionPoint& b) noexcept {
    return a.s == b.s && a.t == b.t;
}
/// Charge a conservative allocation-free comparison workload before each read.
/// A 32-bit-word bound covers the native limb loop on both supported platforms.
template <class Charge>
inline bool partition_same_graph(const PgGraphView& a, const PgGraphView& b, Charge&& charge) {
    const auto equal = [&](const ExactRational& x, const ExactRational& y) {
        charge(1);
        const auto words = [](std::uint64_t bits) { return bits / 32 + (bits % 32 != 0 ? 1 : 0); };
        charge(words(x.numerator_bits()) + words(x.denominator_bits()) + words(y.numerator_bits()) +
               words(y.denominator_bits()));
        return x == y;
    };

    charge(1);
    if (a.n != b.n || a.convention != b.convention || a.cells.size() != b.cells.size())
        return false;
    if (!equal(a.x0, b.x0) || !equal(a.width, b.width) || !equal(a.skew, b.skew))
        return false;
    for (std::size_t i = 0; i < a.cells.size(); ++i) {
        charge(1);
        if (a.cells[i] != b.cells[i])
            return false;
    }
    const auto x = a.height.terms(), y = b.height.terms();
    if (x.size() != y.size())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i) {
        charge(1);
        if (x[i].powers.first != y[i].powers.first || x[i].powers.second != y[i].powers.second ||
            !equal(x[i].coefficient, y[i].coefficient))
            return false;
    }
    return true;
}

inline bool partition_same_graph(const PgGraphView& a, const PgGraphView& b) noexcept {
    return partition_same_graph(a, b, [](std::uint64_t) noexcept {});
}
[[nodiscard]] inline PgGraphView partition_graph(const PgPartitionData& data) {
    if (!data.graph)
        partition_fail(core::FailureCode::InvalidInput);
    return partition_take(data.graph->view());
}
template <class F>
auto partition_entry(F&& f) noexcept -> std::invoke_result_t<F> {
    using Result = std::invoke_result_t<F>;
    return core::guard_entry(partition_operation, core::FailureCode::InvalidInput, {},
                             [&]() -> Result {
                                 try {
                                     return f();
                                 } catch (const PartitionAbort& a) {
                                     return std::unexpected(a.failure);
                                 }
                             });
}
} // namespace uvlas::structure::detail
