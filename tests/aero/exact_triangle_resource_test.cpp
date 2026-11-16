#include "aero/exact_triangle_predicates.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>
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
namespace subject = aero::detail;
constexpr std::uint64_t budget = 4ULL * 1024 * 1024;
const subject::PredicateLimits limits{numerics::linear::EstimateAdmission{28672}, 100};
} // namespace

TEST_CASE("Triangle preparation preserves owners after payload and Boost allocation failures",
          "[negative][UI-D06][N4-02][N4-03][N4-08][N4-09]") {
    core::ResourcePlan preparation{budget}, old_output{budget};
    auto workspace = subject::PredicateWorkspace::create(preparation);
    REQUIRE(workspace);
    const std::array<double, 12> input{0, 0, 0, 3, 0, 0, 0, 4, 0, .75, 1, 20};
    const auto preserved = input;
    const core::FailureDetails context{core::ModelGeneration{1},
                                       core::AcceptedStateId{2},
                                       core::TrialId{3},
                                       core::NumericalGeneration{4},
                                       core::ComponentId{5},
                                       core::ProblemId{core::PanelId{6}},
                                       std::nullopt,
                                       core::AcceptedStateId{7}};
    const auto control = subject::classify_triangle(subject::Binary64Coordinates{input},
                                                    **workspace, old_output, limits, {}, context);
    REQUIRE(control);
    const auto* original = control->outcome();
    REQUIRE(original != nullptr);
    const auto* old = std::get_if<subject::NondegenerateTriangle>(original);
    REQUIRE(old != nullptr);
    const auto old_numerator = old->distance_squared.numerator;
    const auto old_denominator = old->distance_squared.denominator;
    // With workspace prepared, the only owning new before arithmetic is the result
    // array. The next C++ allocation comes from fixed Boost Karatsuba scratch;
    // integer limbs themselves are inline. Neither site is inferred from PG counts.
    for (const std::int64_t site : {0, 1}) {
        core::ResourcePlan output{budget};
        const auto result = [&] {
            const FaultScope fault{site};
            return subject::classify_triangle(subject::Binary64Coordinates{input}, **workspace,
                                              output, limits, {}, context);
        }();
        REQUIRE_FALSE(result);
        CHECK(attempts.load(std::memory_order_relaxed) == static_cast<std::uint64_t>(site + 1));
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(result.error().operation == core::Operation::EvaluateAero);
        auto expected = context;
        if (site == 0)
            expected.resources = core::ResourceUsage{
                core::ResourceKind::Bytes, subject::ExactTriangleResult::payload_bytes(), 0};
        CHECK(result.error().details == expected);
        CHECK(output.allocated_bytes() ==
              (site == 0 ? 0 : subject::ExactTriangleResult::payload_bytes()));
        CHECK(input == preserved);
        CHECK(old->distance_squared.numerator == old_numerator);
        CHECK(old->distance_squared.denominator == old_denominator);
        CHECK(control->diagnostic_context() == context);
        {
            const auto lease = (*workspace)->acquire(context);
            REQUIRE(lease);
        }
        core::ResourcePlan recovery_output{budget};
        const auto recovered = subject::classify_triangle(
            subject::Binary64Coordinates{input}, **workspace, recovery_output, limits, {}, context);
        REQUIRE(recovered);
        const auto* outcome = recovered->outcome();
        REQUIRE(outcome != nullptr);
        const auto* value = std::get_if<subject::NondegenerateTriangle>(outcome);
        REQUIRE(value != nullptr);
        CHECK(value->distance_squared.numerator == old_numerator);
        CHECK(value->distance_squared.denominator == old_denominator);
        CHECK(value->far == subject::FarRelation::Boundary);
    }
}
