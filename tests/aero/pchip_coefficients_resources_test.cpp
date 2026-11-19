#include "aero/pchip_coefficients.hpp"

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
namespace c = uvlas::core;
template <class T>
T take(c::Result<T> value) {
    if (!value)
        FAIL("Unexpected failure " << static_cast<int>(value.error().code));
    return std::move(*value);
}
struct Fixture {
    std::array<std::uint32_t, 1> one{1}, two{2}, four{4};
    std::array<p::RationalView, 3> t{{{0, {}, one}, {1, one, one}, {1, two, one}}};
    std::array<p::RationalView, 3> y{{{0, {}, one}, {1, one, one}, {1, four, one}}};
    [[nodiscard]] p::SampleGrid view() const {
        return {t, y, 1};
    }
};
p::RationalLimits limits() {
    return {64, 100000, UINT64_MAX, 10000};
}
c::FailureDetails details() {
    c::FailureDetails d;
    d.model = c::ModelGeneration{11};
    d.parent = c::AcceptedStateId{12};
    d.trial = c::TrialId{13};
    d.component = c::ComponentId{14};
    d.problem = c::DofId{15};
    return d;
}
constinit std::stop_source* cancel_source = nullptr;
constinit unsigned until_cancel = 0;
void cancel_last() {
    if (--until_cancel == 0)
        cancel_source->request_stop();
    else
        allocation_callback = cancel_last;
}
struct Reentry {
    p::CoefficientWorkspace* workspace = nullptr;
    const Fixture* fixture = nullptr;
    c::ResourcePlan* memory = nullptr;
    std::optional<c::Failure> failure;
};
constinit Reentry* reentry = nullptr;
void try_reentry() {
    const auto result = p::prepare_coefficients(reentry->fixture->view(), *reentry->workspace,
                                                limits(), *reentry->memory);
    if (!result)
        reentry->failure = result.error();
}
} // namespace
TEST_CASE("Private PCHIP coefficient preparation preflights all seven allocations and every fault",
          "[N4-02][N4-03][N4-08][negative]") {
    const Fixture fixture;
    const auto grid = fixture.view();
    const auto plan = [&] {
        const FaultScope audit{INT64_MAX};
        return p::forecast_coefficients(grid, limits());
    }();
    REQUIRE(plan);
    CHECK(attempts.load() == 0);
    p::CoefficientWorkspace workspace;
    c::ResourcePlan measured{plan->payload_bytes()};
    const auto good = [&] {
        const FaultScope audit{INT64_MAX};
        return p::prepare_coefficients(grid, workspace, limits(), measured);
    }();
    REQUIRE(good);
    CHECK(attempts.load() == 7);
    CHECK(measured.allocated_bytes() == plan->payload_bytes());
    const auto previous = take(good->coefficient(0, 0, 0));
    for (std::int64_t fail = 0; fail < 7; ++fail) {
        CAPTURE(fail);
        c::ResourcePlan memory{plan->payload_bytes()};
        const auto d = details();
        const auto result = [&] {
            const FaultScope audit{fail};
            return p::prepare_coefficients(grid, workspace, limits(), memory, {}, d);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == c::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == static_cast<std::uint64_t>(fail + 1));
        CHECK(result.error().details.model == d.model);
        CHECK(result.error().details.parent == d.parent);
        CHECK(result.error().details.trial == d.trial);
        CHECK(result.error().details.component == d.component);
        CHECK(result.error().details.problem == d.problem);
        CHECK(previous.sign == 0);
        CHECK(previous.denominator[0] == 1);
    }
    c::ResourcePlan too_small{plan->payload_bytes() - 1};
    const auto refused = [&] {
        const FaultScope audit{0};
        return p::prepare_coefficients(grid, workspace, limits(), too_small);
    }();
    REQUIRE_FALSE(refused);
    CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
    CHECK(attempts.load() == 0);
    CHECK(too_small.allocated_bytes() == 0);
}
TEST_CASE("Private PCHIP coefficient bit and work boundaries reject before allocating",
          "[N4-02][N4-08][negative]") {
    const Fixture fixture;
    const auto plan = take(p::forecast_coefficients(fixture.view(), limits()));
    auto exact = limits();
    exact.work_bits = plan.arithmetic().integer_bits();
    exact.steps = plan.reserved_steps();
    exact.input_bits = 3;
    exact.slots = plan.arithmetic().slots();
    const auto accepted = take(p::forecast_coefficients(fixture.view(), exact));
    CHECK(accepted.payload_bytes() == plan.payload_bytes());
    p::CoefficientWorkspace workspace;
    c::ResourcePlan memory{plan.payload_bytes()};
    const auto block = take(p::prepare_coefficients(fixture.view(), workspace, exact, memory));
    CHECK(block.used_steps() <= exact.steps);
    for (const auto changed : {0, 1, 2, 3}) {
        auto cap = exact;
        if (changed == 0)
            --cap.steps;
        else if (changed == 1)
            --cap.work_bits;
        else if (changed == 2)
            --cap.input_bits;
        else
            --cap.slots;
        c::ResourcePlan untouched{plan.payload_bytes()};
        const auto refused = [&] {
            const FaultScope audit{0};
            return p::prepare_coefficients(fixture.view(), workspace, cap, untouched);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == 0);
        CHECK(untouched.allocated_bytes() == 0);
    }
}
TEST_CASE("Private PCHIP coefficient final allocation cancellation and reentry preserve ownership",
          "[N4-03][N4-08][negative]") {
    const Fixture fixture;
    const auto plan = take(p::forecast_coefficients(fixture.view(), limits()));
    p::CoefficientWorkspace workspace;
    c::ResourcePlan memory{plan.payload_bytes() * 3};
    std::stop_source stop;
    cancel_source = &stop;
    until_cancel = 7;
    allocation_callback = cancel_last;
    const auto cancelled = p::prepare_coefficients(fixture.view(), workspace, limits(), memory,
                                                   stop.get_token(), details());
    cancel_source = nullptr;
    allocation_callback = nullptr;
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == c::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details());
    CHECK(memory.allocated_bytes() == plan.payload_bytes());
    Reentry context{&workspace, &fixture, &memory, {}};
    reentry = &context;
    allocation_callback = try_reentry;
    const auto outer = p::prepare_coefficients(fixture.view(), workspace, limits(), memory);
    reentry = nullptr;
    allocation_callback = nullptr;
    REQUIRE(outer);
    if (!context.failure) {
        FAIL_CHECK("Nested call was not rejected");
        return;
    }
    CHECK(context.failure->code == c::FailureCode::InvalidInput);
    CHECK(memory.allocated_bytes() == 2 * plan.payload_bytes());
    const auto recovered = p::prepare_coefficients(fixture.view(), workspace, limits(), memory);
    REQUIRE(recovered);
    CHECK(memory.allocated_bytes() == 3 * plan.payload_bytes());
}
TEST_CASE("Private PCHIP canonical mathematical rejection performs only planned allocations",
          "[N4-02][N4-03][N4-08][negative]") {
    Fixture fixture;
    fixture.y[1] = {1, fixture.two, fixture.four};
    const auto plan = take(p::forecast_coefficients(fixture.view(), limits()));
    p::CoefficientWorkspace workspace;
    c::ResourcePlan memory{plan.payload_bytes()};
    const auto result = [&] {
        const FaultScope audit{INT64_MAX};
        return p::prepare_coefficients(fixture.view(), workspace, limits(), memory, {}, details());
    }();
    REQUIRE_FALSE(result);
    CHECK(result.error().code == c::FailureCode::InvalidInput);
    CHECK(attempts.load() == 7);
    CHECK(memory.allocated_bytes() == plan.payload_bytes());
    CHECK(fixture.y[1].numerator[0] == 2);
    CHECK(fixture.y[1].denominator[0] == 4);
}
