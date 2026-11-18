#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <string>
#include <uvlas/structure/pg_coefficients.hpp>
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

TEST_CASE("PG dictionary contains allocation faults and releases its exclusive workspace",
          "[SH-D03][SH-D05][N4-08][negative]") {
    using namespace uvlas;
    using namespace uvlas::structure;
    constexpr auto budget = std::uint64_t{64} * 1024 * 1024 * 1024;
    const PgCoefficientLimits limits{
        {64, 1000, 4096, 32, 32, 65536, 2000000,
         ExactEstimateAdmission{std::uint64_t{16} * 1024 * 1024 * 1024}, ExactTimingMode::Disabled},
        1000,
        20000000,
        2000,
        {}};
    core::ResourcePlan input{budget};
    const auto huge =
        ExactRational::from_integer_pair("1" + std::string(50, '0'), "3", limits.algebra, input);
    REQUIRE(huge);
    const auto one = ExactRational::from_integer_pair("1", "1", limits.algebra, input);
    REQUIRE(one);
    const ExactChartPolynomial zero;
    const auto graph = PgGraph::prepare({*huge, *one, *one, zero, 1, PgGraphConvention::Version1},
                                        {limits.algebra, 64, {}}, input);
    REQUIRE(graph);
    const auto before = graph->view();
    REQUIRE(before);
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    for (const std::int64_t site : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
        core::ResourcePlan output{budget};
        const auto rejected = [&] {
            const FaultScope fault{site};
            return PgCoefficientDictionary::prepare(*graph, limits, workspace, output);
        }();
        REQUIRE_FALSE(rejected);
        CHECK(attempts.load(std::memory_order_relaxed) == static_cast<std::uint64_t>(site + 1));
        CHECK(rejected.error().code == core::FailureCode::ResourceExceeded);
        const auto after = graph->view();
        REQUIRE(after);
        CHECK(after->identity == before->identity);
        const auto lease = workspace.acquire();
        REQUIRE(lease);
        if (site == 0)
            CHECK(output.allocated_bytes() == 0);
    }
    core::ResourcePlan recovered_plan{budget};
    const auto recovered =
        PgCoefficientDictionary::prepare(*graph, limits, workspace, recovered_plan);
    REQUIRE(recovered);
    CHECK(recovered->preparation().observed_columns == 161);
}
