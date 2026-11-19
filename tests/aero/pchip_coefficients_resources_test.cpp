#include "aero/pchip_interval_jet.hpp"

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

namespace {

constexpr std::array<std::uint32_t, 1> half_numerator{1};
constexpr std::array<std::uint32_t, 1> half_denominator{2};
constexpr p::IntervalRequest midpoint{
    0, {1, half_numerator, half_denominator}, p::JetOrder::Second};
constexpr std::int64_t interval_allocations = 7;

p::CoefficientBlock prepare_interval_source() {
    const Fixture fixture;
    p::CoefficientWorkspace workspace;
    c::ResourcePlan memory{UINT64_MAX};
    return take(p::prepare_coefficients(fixture.view(), workspace, limits(), memory));
}

void check_midpoint_unchanged(const p::IntervalValues& value) {
    // With samples 0,1,4 and unit spacing, the first Hermite cubic is
    // (3*u^2-u^3)/2. At u=1/2 its value and first two derivatives are exact.
    const std::array orders{p::JetOrder::Value, p::JetOrder::First, p::JetOrder::Second};
    constexpr std::array<std::uint32_t, 3> numerators{5, 9, 3};
    constexpr std::array<std::uint32_t, 3> denominators{16, 8, 2};
    for (std::size_t i = 0; i < orders.size(); ++i) {
        CAPTURE(i);
        const auto rational = take(value.derivative(0, orders[i]));
        CHECK(rational.sign == 1);
        REQUIRE(rational.numerator.size() == 1);
        REQUIRE(rational.denominator.size() == 1);
        CHECK(rational.numerator.front() == numerators[i]);
        CHECK(rational.denominator.front() == denominators[i]);
    }
}

struct IntervalReentry {
    const p::CoefficientBlock* source = nullptr;
    p::IntervalWorkspace* workspace = nullptr;
    c::ResourcePlan* memory = nullptr;
    std::optional<c::Failure> failure{};
};

constinit IntervalReentry* pending_interval_reentry = nullptr;

void try_interval_reentry() {
    const auto refused = p::evaluate_interval(*pending_interval_reentry->source, midpoint,
                                              *pending_interval_reentry->workspace, limits(),
                                              *pending_interval_reentry->memory);
    if (!refused)
        pending_interval_reentry->failure = refused.error();
}

} // namespace

TEST_CASE("Private PCHIP interval evaluation admits every allocation and preserves prior values",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = prepare_interval_source();
    const auto forecast = [&] {
        const FaultScope audit{0};
        return p::forecast_interval(source, midpoint, limits());
    }();
    REQUIRE(forecast);
    CHECK(attempts.load() == 0);
    auto exact = limits();
    exact.input_bits = 3; // Second derivative imports the fixed integer 6.
    exact.work_bits = forecast->arithmetic().integer_bits();
    exact.steps = forecast->reserved_steps();
    exact.slots = forecast->arithmetic().slots();
    const auto exact_forecast = take(p::forecast_interval(source, midpoint, exact));
    CHECK(exact_forecast.payload_bytes() == forecast->payload_bytes());

    p::IntervalWorkspace workspace;
    c::ResourcePlan initial{forecast->payload_bytes()};
    const auto old = [&] {
        const FaultScope audit{INT64_MAX};
        return p::evaluate_interval(source, midpoint, workspace, exact, initial);
    }();
    REQUIRE(old);
    CHECK(attempts.load() == static_cast<std::uint64_t>(interval_allocations));
    CHECK(initial.allocated_bytes() == forecast->payload_bytes());
    CHECK(old->used_steps() <= exact.steps);
    CHECK(old->used_slots() <= exact.slots);
    check_midpoint_unchanged(*old);

    std::uint64_t prior_payload = 0;
    const auto context = details();
    for (std::int64_t fault = 0; fault < interval_allocations; ++fault) {
        CAPTURE(fault);
        c::ResourcePlan memory{forecast->payload_bytes()};
        const auto rejected = [&] {
            const FaultScope audit{fault};
            return p::evaluate_interval(source, midpoint, workspace, exact, memory, {}, context);
        }();
        REQUIRE_FALSE(rejected);
        const auto& failure = rejected.error();
        CHECK(failure.code == c::FailureCode::ResourceExceeded);
        CHECK(failure.details.model == context.model);
        CHECK(failure.details.parent == context.parent);
        CHECK(failure.details.trial == context.trial);
        CHECK(failure.details.component == context.component);
        CHECK(failure.details.problem == context.problem);
        CHECK(attempts.load() == static_cast<std::uint64_t>(fault + 1));
        if (fault == 0)
            CHECK(memory.allocated_bytes() == 0);
        else
            CHECK(memory.allocated_bytes() > prior_payload);
        CHECK(memory.allocated_bytes() < forecast->payload_bytes());
        prior_payload = memory.allocated_bytes();
        check_midpoint_unchanged(*old);
    }

    c::ResourcePlan insufficient{forecast->payload_bytes() - 1};
    const auto short_bytes = [&] {
        const FaultScope audit{0};
        return p::evaluate_interval(source, midpoint, workspace, exact, insufficient);
    }();
    REQUIRE_FALSE(short_bytes);
    CHECK(short_bytes.error().code == c::FailureCode::ResourceExceeded);
    CHECK(attempts.load() == 0);
    CHECK(insufficient.allocated_bytes() == 0);

    std::array short_caps{exact, exact, exact, exact};
    --short_caps[0].input_bits;
    --short_caps[1].work_bits;
    --short_caps[2].steps;
    --short_caps[3].slots;
    for (std::size_t i = 0; i < short_caps.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{forecast->payload_bytes()};
        const auto rejected = [&] {
            const FaultScope audit{0};
            return p::evaluate_interval(source, midpoint, workspace, short_caps[i], memory);
        }();
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == c::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == 0);
        CHECK(memory.allocated_bytes() == 0);
    }
    c::ResourcePlan final_memory{forecast->payload_bytes()};
    const auto recovered =
        take(p::evaluate_interval(source, midpoint, workspace, exact, final_memory));
    check_midpoint_unchanged(recovered);
}

TEST_CASE("Private PCHIP interval evaluation cancels at final allocation and rejects reentry",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = prepare_interval_source();
    const auto plan = take(p::forecast_interval(source, midpoint, limits()));
    p::IntervalWorkspace workspace;
    c::ResourcePlan memory{plan.payload_bytes() * 3};
    std::stop_source stop;
    cancel_source = &stop;
    until_cancel = static_cast<unsigned>(interval_allocations);
    allocation_callback = cancel_last;
    const auto cancelled = p::evaluate_interval(source, midpoint, workspace, limits(), memory,
                                                stop.get_token(), details());
    allocation_callback = nullptr;
    cancel_source = nullptr;
    REQUIRE_FALSE(cancelled);
    CHECK(stop.stop_requested());
    CHECK(cancelled.error().code == c::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details());
    CHECK(memory.allocated_bytes() == plan.payload_bytes());

    IntervalReentry context{&source, &workspace, &memory, {}};
    pending_interval_reentry = &context;
    allocation_callback = try_interval_reentry;
    const auto outer = p::evaluate_interval(source, midpoint, workspace, limits(), memory);
    allocation_callback = nullptr;
    pending_interval_reentry = nullptr;
    REQUIRE(outer);
    if (context.failure) {
        CHECK(context.failure->code == c::FailureCode::InvalidInput);
    } else {
        FAIL("The shared interval workspace accepted a nested call");
    }
    CHECK(memory.allocated_bytes() == plan.payload_bytes() * 2);
    check_midpoint_unchanged(*outer);
    const auto next = take(p::evaluate_interval(source, midpoint, workspace, limits(), memory));
    CHECK(memory.allocated_bytes() == plan.payload_bytes() * 3);
    check_midpoint_unchanged(next);
}
