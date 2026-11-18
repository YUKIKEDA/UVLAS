#include "pg_partition_common.hpp"
#include "pg_partition_test_support.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <string>
#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
// This process-local fault injector does not ship in the product. Disarm before
// unwinding the first allocation failure so reporting cannot fail recursively.
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
constinit void (*entry_hook)() = nullptr;
constinit uvlas::core::ResourcePlan* watched_plan = nullptr;
constinit bool external_only = false;
void before_allocate() {
    if (const auto callback = std::exchange(entry_hook, nullptr))
        callback();
    // PlannedBuffer reserves before calling new. With no pending reservation,
    // this product path allocates only pinned Boost limbs (no std containers).
    if (external_only && watched_plan != nullptr &&
        watched_plan->byte_limit() - watched_plan->remaining_bytes() -
                watched_plan->allocated_bytes() !=
            0)
        return;
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
using namespace partition_test;
struct Reentry {
    const PgPartitionRequest& request;
    const PgPartitionLimits& limits;
    core::Workspace& workspace;
    core::ResourcePlan& nested_plan;
    std::optional<core::FailureCode> rejected;
};
constinit Reentry* reentry = nullptr;
constinit std::stop_source* stop_request = nullptr;
void invoke_reentry() {
    const auto result = PgIntegrationPartition::prepare(reentry->request, reentry->limits,
                                                        reentry->workspace, reentry->nested_plan);
    if (!result)
        reentry->rejected = result.error().code;
}
void invoke_cancel() {
    (void)stop_request->request_stop();
}
[[nodiscard]] PgPartitionLimits exact_caps(PgPartitionLimits l, const PgPartitionPlan& p,
                                           std::uint64_t layers) {
    l.maximum_input_points = p.input_points;
    l.maximum_segments = p.segments;
    l.maximum_lines = p.line_capacity;
    l.maximum_pieces = p.piece_capacity;
    l.maximum_vertices = p.intersection_capacity;
    l.maximum_edges = p.edge_capacity;
    l.maximum_layers = layers;
    l.maximum_numeric_operations = p.numeric_operations;
    l.maximum_topology_visits = p.topology_visits;
    l.maximum_text_bytes = p.text_bytes;
    l.algebra.maximum_value_bits = p.value_bits;
    l.algebra.maximum_work_bits = p.planned_work_bits;
    l.algebra.external_memory = ExactEstimateAdmission{p.estimated_external_bytes};
    return l;
}
} // namespace
TEST_CASE("PG partition admits exact caps and rejects every one-less cap before allocation",
          "[SH-D07][N4-02][N4-08][negative][property]") {
    Fixture f;
    f.named.push_back(f.point("0", "0"));
    f.named.push_back(f.point("1", "0"));
    f.names.push_back({"root", f.named});
    const auto request = f.request();
    const auto forecast = take(PgIntegrationPartition::plan(request, limits()));
    const auto bounded = exact_caps(limits(), forecast, f.layers.size());
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan admitted_plan{forecast.owned_payload_bytes};
    const auto admitted =
        PgIntegrationPartition::prepare(request, bounded, workspace, admitted_plan);
    REQUIRE(admitted);
    CHECK(admitted_plan.allocated_bytes() <= forecast.owned_payload_bytes);
    const auto old_identity = take(admitted->identity());
    const std::array members{
        &PgPartitionLimits::maximum_input_points,    &PgPartitionLimits::maximum_segments,
        &PgPartitionLimits::maximum_lines,           &PgPartitionLimits::maximum_pieces,
        &PgPartitionLimits::maximum_vertices,        &PgPartitionLimits::maximum_edges,
        &PgPartitionLimits::maximum_layers,          &PgPartitionLimits::maximum_numeric_operations,
        &PgPartitionLimits::maximum_topology_visits, &PgPartitionLimits::maximum_text_bytes};
    for (const auto member : members) {
        auto changed = bounded;
        REQUIRE(changed.*member > 0);
        --(changed.*member);
        core::ResourcePlan plan{cap};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(request, changed, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
        CHECK(plan.allocated_bytes() == 0);
        CHECK(take(admitted->identity()) == old_identity);
    }
    for (const auto member : {&ExactLimits::maximum_value_bits, &ExactLimits::maximum_work_bits}) {
        auto changed = bounded;
        --(changed.algebra.*member);
        core::ResourcePlan plan{cap};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(request, changed, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
    }
    {
        auto changed = bounded;
        changed.algebra.external_memory =
            ExactEstimateAdmission{forecast.estimated_external_bytes - 1};
        core::ResourcePlan plan{cap};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(request, changed, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
    }
    {
        core::ResourcePlan plan{forecast.owned_payload_bytes - 1};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(request, bounded, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
        CHECK(plan.allocated_bytes() == 0);
    }
    {
        auto changed = bounded;
        changed.algebra.external_memory = ExactRequireProvenBound{cap};
        core::ResourcePlan plan{cap};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(request, changed, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::Unsupported);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
    }
}
TEST_CASE("PG partition contains owned and Boost allocation faults and preserves old results",
          "[SH-D07][N4-02][N4-08][negative]") {
    Fixture f;
    const auto tiny = take(ExactRational::from_decimal("1e-400", exact_limits(), f.input_plan));
    f.domain[2].t = take(tiny.clone(exact_limits(), f.input_plan));
    f.domain[3].t = take(tiny.clone(exact_limits(), f.input_plan));
    const auto request = f.request();
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan old_plan{cap};
    const auto old = take(PgIntegrationPartition::prepare(request, limits(), workspace, old_plan));
    const auto id = take(old.identity());
    for (const auto site :
         {std::int64_t{0}, std::int64_t{1}, std::int64_t{5}, std::int64_t{15}, std::int64_t{28},
          std::int64_t{40}, std::int64_t{60}, std::int64_t{100}, std::int64_t{180}}) {
        core::ResourcePlan plan{cap};
        const auto result = [&] {
            const FaultScope fault{site};
            return PgIntegrationPartition::prepare(request, limits(), workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == static_cast<std::uint64_t>(site + 1));
        CHECK(take(old.identity()) == id);
        CHECK(old.validate_for(request));
        const auto lease = workspace.acquire();
        REQUIRE(lease);
    }
    core::ResourcePlan recovered{cap};
    const auto result = PgIntegrationPartition::prepare(request, limits(), workspace, recovered);
    REQUIRE(result);
    CHECK(take(result->identity()) == id);
}
TEST_CASE("PG partition rejects allocator reentry and cancellation without altering old output",
          "[SH-D07][N4-02][N4-08][negative][property]") {
    const Fixture f;
    const auto request = f.request();
    const auto configuration = limits();
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan nested{cap}, outer{cap};
    Reentry context{request, configuration, workspace, nested, {}};
    reentry = &context;
    entry_hook = invoke_reentry;
    const auto result = PgIntegrationPartition::prepare(request, configuration, workspace, outer);
    entry_hook = nullptr;
    reentry = nullptr;
    REQUIRE(result);
    if (!context.rejected) {
        FAIL("Allocator reentry must record its rejected failure");
        return;
    }
    const auto reentry_failure = *context.rejected;
    CHECK(reentry_failure == core::FailureCode::InvalidInput);
    CHECK(nested.allocated_bytes() == 0);
    const auto id = take(result->identity());
    // MSVC request_stop() is nonconst; libstdc++ exposes a const overload.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source stopped;
    (void)stopped.request_stop();
    auto cancelled = configuration;
    cancelled.cancellation = stopped.get_token();
    core::ResourcePlan before{cap};
    const auto rejected = PgIntegrationPartition::prepare(request, cancelled, workspace, before);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::Cancelled);
    CHECK(before.allocated_bytes() == 0);
    std::stop_source during;
    auto mid = configuration;
    mid.cancellation = during.get_token();
    stop_request = &during;
    entry_hook = invoke_cancel;
    core::ResourcePlan partial{cap};
    const auto interrupted = PgIntegrationPartition::prepare(request, mid, workspace, partial);
    entry_hook = nullptr;
    stop_request = nullptr;
    REQUIRE_FALSE(interrupted);
    CHECK(interrupted.error().code == core::FailureCode::Cancelled);
    CHECK(partial.allocated_bytes() > 0);
    CHECK(take(result->identity()) == id);
    CHECK(result->validate_for(request));
    const auto released = workspace.acquire();
    CHECK(released);
}
TEST_CASE("PG partition checked dimensions and graph clone quotas reject without allocation",
          "[SH-D07][N4-02][N4-08][negative]") {
    Fixture f;
    auto overflow = limits();
    overflow.maximum_lines = std::numeric_limits<std::uint64_t>::max();
    CHECK_FALSE(PgIntegrationPartition::plan(f.request(), overflow));
    const auto expect_overflow = [](auto callable) {
        bool rejected = false;
        try {
            callable();
        } catch (const detail::PartitionAbort& e) {
            rejected = e.failure.code == core::FailureCode::ResourceExceeded;
        }
        CHECK(rejected);
    };
    expect_overflow(
        [] { (void)detail::partition_add(std::numeric_limits<std::uint64_t>::max(), 1); });
    expect_overflow(
        [] { (void)detail::partition_mul(std::numeric_limits<std::uint64_t>::max(), 2); });
    expect_overflow(
        [] { (void)detail::partition_index(std::numeric_limits<std::uint64_t>::max()); });
    const auto one = f.r("1");
    const std::array term{ExactInputTerm{{4, 0}, one}};
    const auto height = take(ExactChartPolynomial::make(term, exact_limits(), f.input_plan));
    const auto zero = f.r("0");
    f.graph = take(PgGraph::prepare({zero, one, zero, height, 1, PgGraphConvention::Version1},
                                    {exact_limits(), 64, {}}, f.input_plan));
    for (const auto member :
         {&ExactLimits::maximum_terms, &ExactLimits::maximum_raw_terms,
          &ExactLimits::maximum_axis_degree, &ExactLimits::maximum_total_degree}) {
        auto changed = limits();
        changed.algebra.*member = 0;
        core::ResourcePlan plan{cap};
        core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(f.request(), changed, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
    }
}
TEST_CASE("PG partition contains an external Boost allocation failure after owned arrays",
          "[SH-D07][N4-02][N4-08][negative]") {
    Fixture f;
    const auto tiny = take(ExactRational::from_decimal("1e-400", exact_limits(), f.input_plan));
    f.domain[2].t = take(tiny.clone(exact_limits(), f.input_plan));
    f.domain[3].t = take(tiny.clone(exact_limits(), f.input_plan));
    core::ResourcePlan output{cap};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    watched_plan = &output;
    external_only = true;
    const auto result = [&] {
        const FaultScope fault{0};
        return PgIntegrationPartition::prepare(f.request(), limits(), workspace, output);
    }();
    external_only = false;
    watched_plan = nullptr;
    REQUIRE_FALSE(result);
    CHECK(result.error().code == core::FailureCode::ResourceExceeded);
    CHECK(attempts.load(std::memory_order_relaxed) == 1);
    CHECK(output.allocated_bytes() > 0);
    CHECK(f.domain[2].t == tiny);
    const auto lease = workspace.acquire();
    CHECK(lease);
}

TEST_CASE("PG n4 incidence plan admits equality and rejects each resource one below",
          "[SH-D07][N4-02][N4-08][negative][property]") {
    Fixture f{4};
    incidence_fixture(f);
    const auto forecast = take(PgIntegrationPartition::plan(f.request(), limits()));
    const auto bounded = exact_caps(limits(), forecast, f.layers.size());
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan equal{forecast.owned_payload_bytes};
    const auto accepted = PgIntegrationPartition::prepare(f.request(), bounded, workspace, equal);
    REQUIRE(accepted);
    const auto old = take(accepted->identity());
    for (const auto member : {&PgPartitionLimits::maximum_numeric_operations,
                              &PgPartitionLimits::maximum_topology_visits}) {
        auto below = bounded;
        --(below.*member);
        core::ResourcePlan plan{cap};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(f.request(), below, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
        CHECK(plan.allocated_bytes() == 0);
        CHECK(take(accepted->identity()) == old);
    }
    for (const auto external : {false, true}) {
        auto below = bounded;
        if (external)
            below.algebra.external_memory =
                ExactEstimateAdmission{forecast.estimated_external_bytes - 1};
        core::ResourcePlan plan{external ? cap : forecast.owned_payload_bytes - 1};
        const auto result = [&] {
            const FaultScope fault{0};
            return PgIntegrationPartition::prepare(f.request(), below, workspace, plan);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
        CHECK(plan.allocated_bytes() == 0);
        CHECK(take(accepted->identity()) == old);
    }
}
