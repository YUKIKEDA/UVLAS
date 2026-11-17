#include "aero/local_geometry_packet.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
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

using namespace uvlas;
using namespace uvlas::test_geometry;
using namespace subject;
TEST_CASE("Local packet workspace payload is planned and empty budget rejects",
          "[N4-02][N4-09][UI-D06][negative][verification]") {
    core::ResourcePlan empty{0};
    auto rejected = LocalPacketWorkspace::create(empty);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::ResourceExceeded);
    CHECK(empty.allocated_bytes() == 0);
    core::ResourcePlan plan{LocalPacketWorkspace::payload_bytes()};
    auto workspace = LocalPacketWorkspace::create(plan);
    REQUIRE(workspace);
    CHECK(plan.allocated_bytes() == LocalPacketWorkspace::payload_bytes());
    auto lease = (*workspace)->acquire();
    REQUIRE(lease);
    auto busy = (*workspace)->acquire();
    REQUIRE_FALSE(busy);
    CHECK(busy.error().code == core::FailureCode::InvalidInput);
}
TEST_CASE("Local packet total budgets cancellation and limits are admitted per call",
          "[UI-D06][N4-02][N4-09][negative][verification]") {
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    auto created = LocalPacketWorkspace::create(plan);
    REQUIRE(created);
    auto current = std::move(*created);
    const TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    auto query = take(resolve_query(geometry, {QueryId{5}, SourceVertex{PointId{11}}},
                                    geometry.context(), limits(), *old, plan));
    const auto bytes = ExactTriangleResult::payload_bytes() + LocalGeometryPacket::payload_bytes();
    core::ResourcePlan exact{bytes};
    auto retained = take(prepare_local_packet(
        geometry, query, FacetId{70}, {core::Precision::TwoWord, 0}, limits(), *current, exact));
    CHECK(exact.allocated_bytes() == bytes);
    core::ResourcePlan short_plan{bytes - 1};
    auto denied = prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0},
                                       limits(), *current, short_plan);
    failure(denied, core::FailureCode::ResourceExceeded);
    REQUIRE_FALSE(denied);
    const auto resource = denied.error().cause.details.resources;
    REQUIRE(resource.has_value());
    if (resource)
        CHECK(*resource == core::ResourceUsage{core::ResourceKind::Bytes, bytes,
                                               ExactTriangleResult::payload_bytes()});
    auto tight = limits();
    tight.max_steps = 0;
    auto refused = prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0},
                                        tight, *current, plan);
    failure(refused, core::FailureCode::ResourceExceeded);
    REQUIRE_FALSE(refused);
    CHECK_FALSE(refused.error().problem.has_value());
    tight.max_steps = 1;
    failure(prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0}, tight,
                                 *current, plan),
            core::FailureCode::ResourceExceeded);
    tight = limits();
    tight.predicate.max_steps = 0;
    failure(prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0}, tight,
                                 *current, plan),
            core::FailureCode::ResourceExceeded);
    tight = limits();
    tight.predicate.external_memory = numerics::linear::RequireProvenBound{28672};
    failure(prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0}, tight,
                                 *current, plan),
            core::FailureCode::Unsupported);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC14.51 request_stop is non-const.
    std::stop_source stop;
    CHECK(stop.request_stop());
    tight = limits();
    tight.max_steps = 0;
    auto cancelled =
        prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0}, tight,
                             *current, plan, stop.get_token());
    failure(cancelled, core::FailureCode::Cancelled);
    REQUIRE_FALSE(cancelled);
    CHECK_FALSE(cancelled.error().problem.has_value());
    auto lease = current->acquire();
    REQUIRE(lease);
    failure(prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0},
                                 limits(), *current, plan),
            core::FailureCode::InvalidInput);
    CHECK(retained.active());
}
TEST_CASE("Every actual local preparation allocation failure preserves owners and releases lease",
          "[UI-D06][N4-02][N4-03][N4-09][negative][verification]") {
    core::ResourcePlan setup{budget};
    auto old = workspace(setup);
    auto created = LocalPacketWorkspace::create(setup);
    REQUIRE(created);
    auto current = std::move(*created);
    const TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, setup));
    auto query = take(resolve_query(
        geometry, {QueryId{5}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}},
        geometry.context(), limits(), *old, setup));
    const auto identity = take(geometry.identity());
    core::ResourcePlan baseline_plan{budget};
    auto baseline = [&] {
        const FaultScope count{std::numeric_limits<std::int64_t>::max()};
        return prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0},
                                    limits(), *current, baseline_plan);
    }();
    const auto total = attempts.load(std::memory_order_relaxed);
    REQUIRE(total > 0);
    auto retained = take(std::move(baseline));
    for (std::uint64_t index = 0; index < total; ++index) {
        CAPTURE(index, total);
        core::ResourcePlan failure_plan{budget};
        auto failed = [&] {
            const FaultScope inject{static_cast<std::int64_t>(index)};
            return prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0},
                                        limits(), *current, failure_plan);
        }();
        failure(failed, core::FailureCode::ResourceExceeded);
        CHECK(take(geometry.identity()) == identity);
        CHECK(retained.active());
        auto lease = current->acquire();
        REQUIRE(lease);
    }
    auto recovered = prepare_local_packet(geometry, query, FacetId{70},
                                          {core::Precision::TwoWord, 0}, limits(), *current, setup);
    REQUIRE(recovered);
}
TEST_CASE("Local packet cumulative admission has an exact remainder boundary",
          "[UI-D06][N4-09][negative][verification]") {
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    auto created = LocalPacketWorkspace::create(plan);
    REQUIRE(created);
    auto current = std::move(*created);
    const TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    auto query = take(resolve_query(geometry, {QueryId{5}, SourceVertex{PointId{11}}},
                                    geometry.context(), limits(), *old, plan));
    // Vertex loading: 3 query + 9 target scalars; 17 A0 stages; A2 reserves 12*259.
    constexpr std::uint64_t admission = 3 + 9 + 17 + 3108;
    for (const auto cap : {admission - 1, admission}) {
        auto bounded = limits();
        bounded.max_steps = cap;
        auto result = prepare_local_packet(geometry, query, FacetId{70},
                                           {core::Precision::TwoWord, 0}, bounded, *current, plan);
        if (cap < admission) {
            failure(result, core::FailureCode::ResourceExceeded);
            REQUIRE_FALSE(result);
            CHECK_FALSE(result.error().problem);
        } else {
            REQUIRE(result);
            CHECK(take(result->steps()) <= cap);
        }
    }
    for (const std::uint64_t cap : {16, 17}) {
        auto bounded = limits();
        bounded.predicate.max_steps = cap;
        auto result = prepare_local_packet(geometry, query, FacetId{70},
                                           {core::Precision::TwoWord, 0}, bounded, *current, plan);
        if (cap == 16)
            failure(result, core::FailureCode::ResourceExceeded);
        else
            REQUIRE(result);
    }
}
