#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <uvlas/core/failure.hpp>

namespace uvlas::core {

/// Signed 64-bit array index or extent; object identities use separate unsigned types.
using Index = std::int64_t;

/// Multiply nonnegative extents without signed overflow or allocation.
/// @param extents Array dimensions; a negative dimension is InvalidInput even if another is zero.
/// @param operation Operation to identify in any returned Failure.
/// @return Exact product, zero if any extent is zero, or one for an empty list.
/// A product above the maximum Index returns ResourceExceeded without evaluating that product.
[[nodiscard]] Result<Index> checked_extent_product(std::span<const Index> extents,
                                                   Operation operation) noexcept;

/// Compute the requested payload bytes without overflowing the host allocation size.
/// @param count Nonnegative element count; zero requires no payload storage.
/// @param element_bytes Positive size of one element, normally sizeof(T).
/// @param operation Operation to identify in any returned Failure.
/// @return Exact byte count, InvalidInput for a negative count or zero element size,
/// or ResourceExceeded when the count or product cannot fit in size_t.
[[nodiscard]] Result<std::size_t> checked_byte_count(Index count, std::size_t element_bytes,
                                                     Operation operation) noexcept;

/// Convert a nonnegative extent to a signed 32-bit external API argument.
/// @return Exact value, InvalidInput for a negative value, or ResourceExceeded above INT32_MAX.
[[nodiscard]] Result<std::int32_t> checked_api32(Index value, Operation operation) noexcept;

/// Multiply extents with checked arithmetic and convert the product to signed 32-bit.
/// @return Same empty/negative rules as checked_extent_product; ResourceExceeded above INT32_MAX.
/// This checks the product only. Check each dimension with checked_api32 as well when
/// the external API accepts dimensions separately; a zero product does not validate them.
[[nodiscard]] Result<std::int32_t> checked_api32_product(std::span<const Index> extents,
                                                         Operation operation) noexcept;

} // namespace uvlas::core
