#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <string>
#include <uvlas/structure/exact_polynomial.hpp>
#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
// This process-local fault injector does not ship in the product. Disarm before
// unwinding the first allocation failure so reporting cannot fail recursively.
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
void before_allocate() {
    const auto next = remaining.load(std::memory_order_relaxed);
    if (next == -1)
        return;
    attempts.fetch_add(1, std::memory_order_relaxed);
    if (next == 0) {
        remaining.store(-1, std::memory_order_relaxed);
        throw std::bad_alloc{};
    }
    if (next > 0)
        remaining.fetch_sub(1, std::memory_order_relaxed);
}
void* allocate(std::size_t size) {
    before_allocate();
    // NOLINTNEXTLINE(misc-const-correctness): operator new returns writable storage.
    void* value = std::malloc(size == 0 ? 1 : size);
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
void* allocate_aligned(std::size_t size, std::size_t alignment) {
    before_allocate();
    size = size == 0 ? 1 : size;
#if defined(_WIN32)
    // NOLINTNEXTLINE(misc-const-correctness): operator new returns writable storage.
    void* value = _aligned_malloc(size, alignment);
#else
    if (size > static_cast<std::size_t>(-1) - (alignment - 1))
        throw std::bad_alloc{};
    // NOLINTNEXTLINE(misc-const-correctness): operator new returns writable storage.
    void* value = std::aligned_alloc(alignment, ((size + alignment - 1) / alignment) * alignment);
#endif
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
void free_aligned(void* value) noexcept {
#if defined(_WIN32)
    _aligned_free(value);
#else
    std::free(value);
#endif
}
struct FaultScope {
    explicit FaultScope(std::int64_t count) noexcept {
        attempts.store(0, std::memory_order_relaxed);
        remaining.store(count, std::memory_order_relaxed);
    }
    ~FaultScope() {
        remaining.store(-1, std::memory_order_relaxed);
    }
    FaultScope(const FaultScope&) = delete;
    FaultScope& operator=(const FaultScope&) = delete;
};
} // namespace
void* operator new(std::size_t size) {
    return allocate(size);
}
void* operator new[](std::size_t size) {
    return allocate(size);
}
void operator delete(void* value) noexcept {
    std::free(value);
}
void operator delete[](void* value) noexcept {
    std::free(value);
}
void operator delete(void* value, std::size_t) noexcept {
    std::free(value);
}
void operator delete[](void* value, std::size_t) noexcept {
    std::free(value);
}
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
    try {
        return allocate(size);
    } catch (...) {
        return nullptr;
    }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
    try {
        return allocate(size);
    } catch (...) {
        return nullptr;
    }
}
void operator delete(void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void operator delete[](void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void* operator new(std::size_t size, std::align_val_t alignment) {
    return allocate_aligned(size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
    return allocate_aligned(size, static_cast<std::size_t>(alignment));
}
void operator delete(void* value, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete(void* value, std::size_t, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::size_t, std::align_val_t) noexcept {
    free_aligned(value);
}
void* operator new(std::size_t size, std::align_val_t alignment, const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(size, static_cast<std::size_t>(alignment));
    } catch (...) {
        return nullptr;
    }
}
void* operator new[](std::size_t size, std::align_val_t alignment, const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(size, static_cast<std::size_t>(alignment));
    } catch (...) {
        return nullptr;
    }
}
void operator delete(void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}

namespace {
using namespace uvlas;
using namespace uvlas::structure;
constexpr std::uint64_t budget = std::uint64_t{64} * 1024 * 1024;
ExactLimits limits() {
    return {10000,
            1000,
            10000,
            100,
            200,
            4096,
            1000000,
            ExactEstimateAdmission{budget},
            ExactTimingMode::Disabled};
}
} // namespace
TEST_CASE("Exact preparation translates both external and owned allocation failures",
          "[SH-D03][N4-02][N4-08][negative]") {
    core::ResourcePlan inputs{budget};
    const std::string numerator = "1" + std::string(200, '0');
    const auto source = ExactRational::from_integer_pair(numerator, "3", limits(), inputs);
    REQUIRE(source);
    const auto control = source->clone(limits(), inputs);
    REQUIRE(control);
    for (const std::int64_t site : {0, 1}) {
        core::ResourcePlan output{budget};
        const auto result = [&]() {
            const FaultScope fault{site};
            return source->clone(limits(), output);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(result.error().operation == core::Operation::PrepareModel);
        CHECK(output.allocated_bytes() == 0);
        CHECK(*source == *control);
    }
    const auto recovered = source->clone(limits(), inputs);
    REQUIRE(recovered);
    CHECK(*recovered == *source);
}
TEST_CASE("Exact admission rejects known bad input and budgets without allocating",
          "[SH-D03][N4-02][N4-08][negative]") {
    core::ResourcePlan output{budget};
    auto l = limits();
    l.external_memory = ExactRequireProvenBound{budget};
    const auto unproved = [&]() {
        const FaultScope counter{-2};
        return ExactRational::from_integer_pair("1", "3", l, output);
    }();
    REQUIRE_FALSE(unproved);
    CHECK(unproved.error().code == core::FailureCode::Unsupported);
    CHECK(attempts.load() == 0);
    CHECK(output.allocated_bytes() == 0);
    l = limits();
    l.external_memory = ExactEstimateAdmission{0};
    const auto refused = [&]() {
        const FaultScope counter{-2};
        return ExactRational::from_integer_pair("1", "3", l, output);
    }();
    REQUIRE_FALSE(refused);
    CHECK(refused.error().code == core::FailureCode::ResourceExceeded);
    CHECK(attempts.load() == 0);
    CHECK(output.allocated_bytes() == 0);
    const auto invalid = [&]() {
        const FaultScope counter{-2};
        return ExactRational::from_integer_pair("1", "0", l, output);
    }();
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
    CHECK(attempts.load() == 0);
}
