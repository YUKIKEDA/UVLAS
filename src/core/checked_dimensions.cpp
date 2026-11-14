#include <limits>
#include <uvlas/core/checked_dimensions.hpp>

namespace uvlas::core {

Result<Index> checked_extent_product(std::span<const Index> extents, Operation operation) noexcept {
    bool empty = false;
    // Validate every axis first: zero cannot hide an invalid negative dimension,
    // and an empty array needs no intermediate multiplication (N4-02).
    for (const Index extent : extents) {
        if (extent < 0) {
            return std::unexpected(Failure{FailureCode::InvalidInput, operation});
        }
        empty = empty || extent == 0;
    }
    if (empty) {
        return Index{0};
    }
    Index product = 1;
    for (const Index extent : extents) {
        if (product > std::numeric_limits<Index>::max() / extent) {
            return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
        }
        product *= extent;
    }
    return product;
}

Result<std::size_t> checked_byte_count(Index count, std::size_t element_bytes,
                                       Operation operation) noexcept {
    if (count < 0 || element_bytes == 0) {
        return std::unexpected(Failure{FailureCode::InvalidInput, operation});
    }
    if (static_cast<std::uint64_t>(count) > std::numeric_limits<std::size_t>::max()) {
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
    }
    const auto elements = static_cast<std::size_t>(count);
    // Division precedes multiplication, so the overflow check itself cannot overflow.
    if (elements > std::numeric_limits<std::size_t>::max() / element_bytes) {
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
    }
    return elements * element_bytes;
}

Result<std::int32_t> checked_api32(Index value, Operation operation) noexcept {
    if (value < 0) {
        return std::unexpected(Failure{FailureCode::InvalidInput, operation});
    }
    if (value > std::numeric_limits<std::int32_t>::max()) {
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
    }
    return static_cast<std::int32_t>(value);
}

Result<std::int32_t> checked_api32_product(std::span<const Index> extents,
                                           Operation operation) noexcept {
    const auto product = checked_extent_product(extents, operation);
    if (!product) {
        return std::unexpected(product.error());
    }
    return checked_api32(*product, operation);
}

} // namespace uvlas::core
