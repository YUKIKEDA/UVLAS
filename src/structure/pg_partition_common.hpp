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
inline bool partition_same_graph(const PgGraphView& a, const PgGraphView& b) noexcept {
    if (a.n != b.n || a.convention != b.convention || !(a.x0 == b.x0) || !(a.width == b.width) ||
        !(a.skew == b.skew) || a.cells.size() != b.cells.size())
        return false;
    if (!std::equal(a.cells.begin(), a.cells.end(), b.cells.begin()))
        return false;
    const auto x = a.height.terms(), y = b.height.terms();
    if (x.size() != y.size())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (x[i].powers.first != y[i].powers.first || x[i].powers.second != y[i].powers.second ||
            !(x[i].coefficient == y[i].coefficient))
            return false;
    return true;
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
