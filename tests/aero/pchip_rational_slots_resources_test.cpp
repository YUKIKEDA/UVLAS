#include "aero/pchip_slot_private.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <optional>
#if defined(_WIN32)
#include <malloc.h>
#endif
namespace {
// This process-local fault injector does not ship in the product. Disarm before
// unwinding the first allocation failure so reporting cannot fail recursively.
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
constinit void (*allocation_callback)() = nullptr;
void before_allocate() {
    if (allocation_callback != nullptr) {
        const auto callback = allocation_callback;
        allocation_callback = nullptr;
        callback();
    }
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
namespace p = uvlas::aero::detail::pchip;
using namespace uvlas;
template <class T>
T take(core::Result<T> result) {
    if (!result) {
        FAIL("unexpected failure " << static_cast<int>(result.error().code));
    }
    return std::move(*result);
}
core::Result<p::RationalPlan> program(std::uint64_t steps = UINT64_MAX) noexcept {
    p::Counter c{{128, 1024, steps, 64}};
    auto a = c.input({64, 64}), b = c.input({64, 64});
    if (!a)
        return std::unexpected(a.error());
    if (!b)
        return std::unexpected(b.error());
    for (const auto op : {p::RationalOp::Add, p::RationalOp::Subtract, p::RationalOp::Multiply,
                          p::RationalOp::Divide}) {
        auto value = c.binary(op, *a, *b);
        if (!value)
            return std::unexpected(value.error());
    }
    for (const auto op : {p::RationalOp::Clone, p::RationalOp::Negate}) {
        auto value = c.unary(op, *a);
        if (!value)
            return std::unexpected(value.error());
    }
    auto floor = c.floor(*a);
    if (!floor)
        return std::unexpected(floor.error());
    auto compare = c.compare(*a, *b);
    if (!compare)
        return std::unexpected(compare.error());
    return c.finish();
}
core::FailureDetails details() {
    core::FailureDetails d;
    d.model = core::ModelGeneration{11};
    d.parent = core::AcceptedStateId{12};
    d.trial = core::TrialId{13};
    d.numerical_generation = core::NumericalGeneration{14};
    d.component = core::ComponentId{15};
    d.problem = core::DofId{16};
    d.last_accepted = core::AcceptedStateId{17};
    return d;
}
constinit std::stop_source* final_cancel = nullptr;
constinit unsigned cancel_after = 0;
void cancel_last() {
    if (--cancel_after == 0)
        final_cancel->request_stop();
    else
        allocation_callback = cancel_last;
}
struct Nested {
    const p::RationalPlan* plan = nullptr;
    core::ResourcePlan* memory = nullptr;
    std::optional<p::Workspace> retained;
    bool succeeded = false;
};
constinit Nested* nested = nullptr;
void nested_then_fail() {
    auto result = p::Workspace::create(*nested->plan, *nested->memory);
    nested->succeeded = bool(result);
    if (result)
        nested->retained.emplace(std::move(*result));
    remaining.store(0);
}
} // namespace
TEST_CASE("Private PCHIP workspace plans all allocations and preserves diagnostics at every fault",
          "[N4-02][N4-03][N4-08][negative]") {
    const auto plan = take(program());
    const auto d = details();
    core::ResourcePlan measured{plan.payload_bytes()};
    auto good = [&] {
        const FaultScope scope{INT64_MAX};
        return p::Workspace::create(plan, measured, d);
    }();
    REQUIRE(good);
    CHECK(attempts.load() == 3);
    CHECK(measured.allocated_bytes() == plan.payload_bytes());
    for (std::int64_t at = 0; at < 3; ++at) {
        core::ResourcePlan memory{plan.payload_bytes()};
        auto failed = [&] {
            const FaultScope scope{at};
            return p::Workspace::create(plan, memory, d);
        }();
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
        CHECK(failed.error().operation == core::Operation::EvaluateAero);
        auto actual = failed.error().details;
        actual.resources.reset();
        CHECK(actual == d);
        const auto expected =
            at == 0   ? 0
            : at == 1 ? sizeof(p::slot_detail::WorkspaceData)
                      : plan.payload_bytes() - plan.slots() * sizeof(p::slot_detail::SlotMetadata);
        CHECK(memory.allocated_bytes() == expected);
        auto session = take(good->begin());
        CHECK(session.integer(1)); // Older workspace survives every new failure.
    }
    core::ResourcePlan short_memory{plan.payload_bytes() - 1};
    auto rejected = [&] {
        const FaultScope track{INT64_MAX};
        return p::Workspace::create(plan, short_memory, d);
    }();
    REQUIRE_FALSE(rejected);
    CHECK(attempts.load() == 0);
    CHECK(short_memory.allocated_bytes() == 0);
    CHECK(rejected.error().details == d);
    const core::ResourcePlan unused{UINT64_MAX};
    auto short_work = [&] {
        const FaultScope track{INT64_MAX};
        return program(plan.reserved_steps() - 1);
    }();
    REQUIRE_FALSE(short_work);
    CHECK(attempts.load() == 0);
    CHECK(unused.allocated_bytes() == 0);
}
TEST_CASE("Private PCHIP prepared arithmetic and routine refusals use no C++ allocation",
          "[N4-02][N4-06][N4-08][negative]") {
    const auto plan = take(program());
    core::ResourcePlan memory{plan.payload_bytes()};
    auto workspace = take(p::Workspace::create(plan, memory));
    const std::array<std::uint32_t, 1> n{17}, d{5}, one{1};
    const auto diagnostic = details();
    auto session = take(workspace.begin({}, diagnostic));
    auto outcomes = [&] {
        const FaultScope track{INT64_MAX};
        std::array<bool, 10> result{};
        auto a = session.input({1, n, one}), b = session.input({1, d, one});
        if (!a || !b)
            return result;
        result[0] = bool(session.binary(p::RationalOp::Add, *a, *b));
        result[1] = bool(session.binary(p::RationalOp::Subtract, *a, *b));
        result[2] = bool(session.binary(p::RationalOp::Multiply, *a, *b));
        result[3] = bool(session.binary(p::RationalOp::Divide, *a, *b));
        result[4] = bool(session.unary(p::RationalOp::Clone, *a));
        result[5] = bool(session.unary(p::RationalOp::Negate, *a));
        result[6] = bool(session.floor(*a));
        result[7] = bool(session.compare(*a, *b));
        result[8] = bool(session.view(*a));
        auto quota = session.integer(1);
        result[9] = !quota && quota.error().code == core::FailureCode::ResourceExceeded;
        return result;
    }();
    CHECK(attempts.load() == 0);
    for (const auto success : outcomes)
        CHECK(success);
    // Only C++ new/delete are instrumented; absence of C heap calls is a separate source/call-graph
    // check. The numerical executor never throws for expected failure or cancellation.
    auto failed_context = session.integer(1);
    REQUIRE_FALSE(failed_context);
    CHECK(failed_context.error().details == diagnostic);
}
TEST_CASE(
    "Private PCHIP invalid canonical input zero division and cancellation stay allocation-free",
    "[N4-08][negative]") {
    const auto plan = take(program());
    core::ResourcePlan memory{plan.payload_bytes()};
    auto workspace = take(p::Workspace::create(plan, memory));
    const std::array<std::uint32_t, 1> one{1}, two{2}, four{4};
    {
        auto session = take(workspace.begin());
        auto baseline = take(session.input({1, one, one}));
        const auto before = take(session.view(baseline));
        auto result = [&] {
            const FaultScope track{INT64_MAX};
            return session.input({1, two, four});
        }();
        CHECK(attempts.load() == 0);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::InvalidInput);
        CHECK(before.numerator[0] == 1);
    }
    {
        auto session = take(workspace.begin());
        auto a = take(session.integer(1)), b = take(session.integer(0));
        auto result = [&] {
            const FaultScope track{INT64_MAX};
            return session.binary(p::RationalOp::Divide, a, b);
        }();
        CHECK(attempts.load() == 0);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::InvalidInput);
    }
    // MSVC request_stop() is non-const; libstdc++ declares it const.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source stop;
    {
        auto session = take(workspace.begin(stop.get_token()));
        auto a = take(session.integer(1));
        stop.request_stop();
        auto result = [&] {
            const FaultScope track{INT64_MAX};
            return session.unary(p::RationalOp::Negate, a);
        }();
        CHECK(attempts.load() == 0);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::Cancelled);
    }
    {
        auto session = take(workspace.begin());
        CHECK(session.integer(1));
    }
}
TEST_CASE(
    "Private PCHIP final allocation cancellation and nested reservation retain only owned charges",
    "[N4-02][N4-08][N4-09][negative]") {
    const auto plan = take(program());
    // The callback escapes a mutable stop_source pointer, so this declaration is not const.
    std::stop_source stop;
    final_cancel = &stop;
    cancel_after = 3;
    allocation_callback = cancel_last;
    core::ResourcePlan cancelled_memory{plan.payload_bytes()};
    auto cancelled = p::Workspace::create(plan, cancelled_memory, details(), stop.get_token());
    final_cancel = nullptr;
    allocation_callback = nullptr;
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details());
    CHECK(cancelled_memory.allocated_bytes() == plan.payload_bytes());
    core::ResourcePlan memory{plan.payload_bytes() * 2};
    Nested context{&plan, &memory, {}, false};
    nested = &context;
    allocation_callback = nested_then_fail;
    auto failed = p::Workspace::create(plan, memory);
    nested = nullptr;
    allocation_callback = nullptr;
    remaining.store(-1);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
    REQUIRE(context.succeeded);
    if (!context.retained) {
        FAIL_CHECK("Missing nested workspace");
        return;
    }
    CHECK(memory.allocated_bytes() == plan.payload_bytes());
    auto& retained = context.retained.value();
    auto session = take(retained.begin());
    CHECK(session.integer(1));
    CHECK_FALSE(p::slot_detail::next_epoch(UINT64_MAX));
}
TEST_CASE("Private PCHIP step meter accepts the exact budget then refuses before another visit",
          "[N4-02][N4-08][negative]") {
    p::slot_detail::WorkspaceData data;
    data.plan = take(program());
    data.details = details();
    p::slot_detail::Engine engine{data};
    const auto visits = data.plan.reserved_steps();
    auto outcomes = [&] {
        const FaultScope track{INT64_MAX};
        return std::array<bool, 2>{engine.step(visits), engine.step(1)};
    }();
    CHECK(attempts.load() == 0);
    CHECK(outcomes[0]);
    CHECK_FALSE(outcomes[1]);
    CHECK(data.steps == visits);
    CHECK(engine.error().code == core::FailureCode::ResourceExceeded);
    CHECK(engine.error().details == details());
}
