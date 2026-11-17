#include "geometry_snapshot_fixture.hpp"

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
using namespace uvlas::test_geometry;
using namespace subject;
} // namespace
TEST_CASE("Geometry allocation failures preserve owners and release workspace leases",
          "[UI-D05][UI-D06][N4-02][N4-03][N4-09][negative]") {
    const auto details = context().diagnostics;
    core::ResourcePlan short_workspace{PredicateWorkspace::payload_bytes()};
    const auto denied_workspace = GeometryWorkspace::create(short_workspace, details);
    REQUIRE_FALSE(denied_workspace.has_value());
    CHECK(denied_workspace.error().code == core::FailureCode::ResourceExceeded);
    CHECK(denied_workspace.error().details.parent == details.parent);
    CHECK(denied_workspace.error().details.trial == details.trial);
    CHECK(denied_workspace.error().details.last_accepted == details.last_accepted);
    REQUIRE(denied_workspace.error().details.resources.has_value());
    TriangleFixture fixture;
    core::ResourcePlan prep{budget}, old_plan{budget};
    auto w = workspace(prep);
    auto g = take(prepare_owned_geometry(fixture.input(), limits(), *w, old_plan));
    const auto identity = take(g.identity());
    for (const std::int64_t site : {0, 1}) {
        core::ResourcePlan plan{budget};
        const auto denied = [&] {
            const FaultScope fault{site};
            return prepare_owned_geometry(fixture.input(), limits(), *w, plan);
        }();
        failure(denied, core::FailureCode::ResourceExceeded);
        CHECK(take(g.identity()) == identity);
        auto retry = prepare_owned_geometry(fixture.input(), limits(), *w, plan);
        REQUIRE(retry.has_value());
        CHECK(take(retry->identity()) == identity);
    }
    core::ResourcePlan qp{budget};
    const QueryRequest request{QueryId{1},
                               SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}};
    const auto denied_query = [&] {
        const FaultScope fault{0};
        return resolve_query(g, request, g.context(), limits(), *w, qp);
    }();
    failure(denied_query, core::FailureCode::ResourceExceeded);
    CHECK(qp.allocated_bytes() == 0);
    auto query = take(resolve_query(g, request, g.context(), limits(), *w, qp));
    core::ResourcePlan old_result{budget};
    const auto saved =
        take(classify_resolved_query(g, query, FacetId{70}, limits(), *w, old_result));
    REQUIRE(saved.outcome() != nullptr);
    for (const std::int64_t site : {0, 1}) {
        core::ResourcePlan result_plan{budget};
        const auto denied = [&] {
            const FaultScope fault{site};
            return classify_resolved_query(g, query, FacetId{70}, limits(), *w, result_plan);
        }();
        failure(denied, core::FailureCode::ResourceExceeded);
        CHECK(take(g.identity()) == identity);
        REQUIRE(saved.outcome() != nullptr);
        REQUIRE(
            classify_resolved_query(g, query, FacetId{70}, limits(), *w, result_plan).has_value());
    }
}
