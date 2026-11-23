#include "aero/pchip_boundary.hpp"
#include "aero/pchip_curve.hpp"
#include "aero/pchip_curve_plan.hpp"
#include "aero/pchip_domain_route.hpp"
#include "aero/pchip_interval_jet.hpp"
#include "aero/pchip_knot_locator.hpp"
#include "aero/pchip_preparation_input.hpp"
#include "aero/pchip_shared_preparation_input.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
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

namespace {
constexpr std::int64_t locator_allocations = 3;

void check_midpoint_location(const p::KnotLocation& location) {
    CHECK_FALSE(location.is_knot());
    CHECK(location.nodes() == 3);
    CHECK(take(location.left_interval()) == 0);
    CHECK(take(location.right_interval()) == 0);
    CHECK(take(location.unique_interval()) == 0);
}

struct LocatorReentry {
    const p::CoefficientBlock* source = nullptr;
    p::LocatorWorkspace* workspace = nullptr;
    c::ResourcePlan* memory = nullptr;
    std::optional<c::Failure> failure{};
};
constinit LocatorReentry* pending_locator_reentry = nullptr;

void try_locator_reentry() {
    const auto rejected = p::locate_time(*pending_locator_reentry->source, midpoint.time,
                                         *pending_locator_reentry->workspace, limits(),
                                         *pending_locator_reentry->memory);
    if (!rejected)
        pending_locator_reentry->failure = rejected.error();
}
} // namespace

TEST_CASE("Private PCHIP locator admits the whole arena and retains old scalar results",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = prepare_interval_source();
    const auto forecast = [&] {
        const FaultScope audit{0};
        return p::forecast_location(source, midpoint.time, limits());
    }();
    REQUIRE(forecast);
    CHECK(attempts.load() == 0);
    CHECK(forecast->maximum_comparisons() == 5);
    auto exact = limits();
    exact.input_bits = 2; // Widest stored knot (2/1) and query (1/2).
    exact.work_bits = forecast->arithmetic().integer_bits();
    exact.steps = forecast->reserved_steps();
    exact.slots = forecast->arithmetic().slots();
    const auto payload = forecast->payload_bytes();

    p::LocatorWorkspace workspace;
    c::ResourcePlan initial{payload};
    const auto old = [&] {
        const FaultScope audit{INT64_MAX};
        return p::locate_time(source, midpoint.time, workspace, exact, initial);
    }();
    REQUIRE(old);
    CHECK(initial.allocated_bytes() == payload);
    CHECK(attempts.load() == static_cast<std::uint64_t>(locator_allocations));
    CHECK(old->comparisons() <= forecast->maximum_comparisons());
    CHECK(old->used_slots() <= exact.slots);
    CHECK(old->used_steps() <= exact.steps);
    check_midpoint_location(*old);

    std::uint64_t committed = 0;
    const auto context = details();
    for (std::int64_t fault = 0; fault < locator_allocations; ++fault) {
        CAPTURE(fault);
        c::ResourcePlan memory{payload};
        const auto rejected = [&] {
            const FaultScope audit{fault};
            return p::locate_time(source, midpoint.time, workspace, exact, memory, {}, context);
        }();
        REQUIRE_FALSE(rejected);
        const auto& failure = rejected.error();
        CHECK(failure.code == c::FailureCode::ResourceExceeded);
        CHECK(failure.details.model == context.model);
        CHECK(failure.details.component == context.component);
        CHECK(failure.details.problem == context.problem);
        CHECK(attempts.load() == static_cast<std::uint64_t>(fault + 1));
        if (fault == 0)
            CHECK(memory.allocated_bytes() == 0);
        else
            CHECK(memory.allocated_bytes() > committed);
        CHECK(memory.allocated_bytes() < payload);
        committed = memory.allocated_bytes();
        check_midpoint_location(*old);
    }

    std::array caps{exact, exact, exact, exact, exact};
    --caps[0].input_bits;
    --caps[1].work_bits;
    --caps[2].steps;
    --caps[3].slots;
    for (std::size_t i = 0; i < caps.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{i == 4 ? payload - 1 : payload};
        const auto refused = [&] {
            const FaultScope audit{0};
            return p::locate_time(source, midpoint.time, workspace, caps[i], memory);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(memory.allocated_bytes() == 0);
        CHECK(attempts.load() == 0);
    }
    c::ResourcePlan retry{payload};
    check_midpoint_location(take(p::locate_time(source, midpoint.time, workspace, exact, retry)));
}

TEST_CASE("Private PCHIP locator cancels at the last allocation and rejects reentry",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = prepare_interval_source();
    const auto plan = take(p::forecast_location(source, midpoint.time, limits()));
    const auto payload = plan.payload_bytes();
    p::LocatorWorkspace workspace;
    c::ResourcePlan memory{payload * 3};
    std::stop_source stop;
    cancel_source = &stop;
    until_cancel = static_cast<unsigned>(locator_allocations);
    allocation_callback = cancel_last;
    const auto cancelled = p::locate_time(source, midpoint.time, workspace, limits(), memory,
                                          stop.get_token(), details());
    allocation_callback = nullptr;
    cancel_source = nullptr;
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == c::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details());
    CHECK(memory.allocated_bytes() == payload);

    LocatorReentry context{&source, &workspace, &memory, {}};
    pending_locator_reentry = &context;
    allocation_callback = try_locator_reentry;
    const auto outer = p::locate_time(source, midpoint.time, workspace, limits(), memory);
    allocation_callback = nullptr;
    pending_locator_reentry = nullptr;
    REQUIRE(outer);
    if (context.failure) {
        CHECK(context.failure->code == c::FailureCode::InvalidInput);
    } else {
        FAIL("The shared locator workspace accepted a nested call");
    }
    CHECK(memory.allocated_bytes() == payload * 2);
    check_midpoint_location(*outer);
    const auto next = take(p::locate_time(source, midpoint.time, workspace, limits(), memory));
    CHECK(memory.allocated_bytes() == payload * 3);
    check_midpoint_location(next);
}

namespace {
constexpr std::array<std::uint32_t, 1> hold_end_numerator{2};
constexpr p::RationalView hold_end_query{1, hold_end_numerator, half_numerator};
constexpr auto hold_policy = p::DomainPolicy::Hold;
constexpr std::int64_t domain_allocations = 3;

void check_hold_endpoint(const p::DomainRoute& route) {
    CHECK(route.relation() == p::DomainRelation::AtLast);
    CHECK(route.has_inside());
    CHECK(route.has_constant());
    CHECK(take(route.constant_node()) == 2);
}

struct DomainReentry {
    const p::CoefficientBlock* source = nullptr;
    p::DomainWorkspace* workspace = nullptr;
    c::ResourcePlan* memory = nullptr;
    std::optional<c::Failure> failure{};
};
constinit DomainReentry* pending_domain_reentry = nullptr;

void try_domain_reentry() {
    const auto refused = p::route_domain(*pending_domain_reentry->source, hold_end_query,
                                         hold_policy, *pending_domain_reentry->workspace, limits(),
                                         *pending_domain_reentry->memory);
    if (!refused)
        pending_domain_reentry->failure = refused.error();
}
} // namespace

TEST_CASE("Private PCHIP domain route preflights the arena and preserves old candidates",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = prepare_interval_source();
    const auto forecast = [&] {
        const FaultScope audit{0};
        return p::forecast_domain(source, hold_end_query, hold_policy, limits());
    }();
    REQUIRE(forecast);
    CHECK(attempts.load() == 0);

    auto exact = limits();
    exact.input_bits = 2;
    exact.work_bits = forecast->arithmetic().integer_bits();
    exact.steps = forecast->reserved_steps();
    exact.slots = forecast->arithmetic().slots();
    const auto payload = forecast->payload_bytes();
    p::DomainWorkspace workspace;
    c::ResourcePlan initial{payload};
    const auto old = [&] {
        const FaultScope audit{INT64_MAX};
        return p::route_domain(source, hold_end_query, hold_policy, workspace, exact, initial);
    }();
    REQUIRE(old);
    CHECK(attempts.load() == static_cast<std::uint64_t>(domain_allocations));
    CHECK(initial.allocated_bytes() == payload);
    CHECK(old->used_steps() <= exact.steps);
    CHECK(old->used_slots() == 3);
    check_hold_endpoint(*old);

    std::uint64_t prior_payload = 0;
    const auto context = details();
    for (std::int64_t fault = 0; fault < domain_allocations; ++fault) {
        CAPTURE(fault);
        c::ResourcePlan memory{payload};
        const auto refused = [&] {
            const FaultScope audit{fault};
            return p::route_domain(source, hold_end_query, hold_policy, workspace, exact, memory,
                                   {}, context);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(refused.error().details.model == context.model);
        CHECK(refused.error().details.component == context.component);
        CHECK(attempts.load() == static_cast<std::uint64_t>(fault + 1));
        if (fault == 0)
            CHECK(memory.allocated_bytes() == 0);
        else
            CHECK(memory.allocated_bytes() > prior_payload);
        CHECK(memory.allocated_bytes() < payload);
        prior_payload = memory.allocated_bytes();
        check_hold_endpoint(*old);
    }

    std::array caps{exact, exact, exact, exact, exact};
    --caps[0].input_bits;
    --caps[1].work_bits;
    --caps[2].steps;
    --caps[3].slots;
    for (std::size_t i = 0; i < caps.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{i == 4 ? payload - 1 : payload};
        const auto refused = [&] {
            const FaultScope audit{0};
            return p::route_domain(source, hold_end_query, hold_policy, workspace, caps[i], memory);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == 0);
        CHECK(memory.allocated_bytes() == 0);
    }
    c::ResourcePlan retry{payload};
    check_hold_endpoint(
        take(p::route_domain(source, hold_end_query, hold_policy, workspace, exact, retry)));
}

TEST_CASE("Private PCHIP domain route cancels at the last allocation and rejects reentry",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = prepare_interval_source();
    const auto plan = take(p::forecast_domain(source, hold_end_query, hold_policy, limits()));
    const auto payload = plan.payload_bytes();
    p::DomainWorkspace workspace;
    c::ResourcePlan memory{payload * 3};
    std::stop_source stop;
    cancel_source = &stop;
    until_cancel = static_cast<unsigned>(domain_allocations);
    allocation_callback = cancel_last;
    const auto cancelled = p::route_domain(source, hold_end_query, hold_policy, workspace, limits(),
                                           memory, stop.get_token(), details());
    allocation_callback = nullptr;
    cancel_source = nullptr;
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == c::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details());
    CHECK(memory.allocated_bytes() == payload);

    DomainReentry context{&source, &workspace, &memory, {}};
    pending_domain_reentry = &context;
    allocation_callback = try_domain_reentry;
    const auto outer =
        p::route_domain(source, hold_end_query, hold_policy, workspace, limits(), memory);
    allocation_callback = nullptr;
    pending_domain_reentry = nullptr;
    REQUIRE(outer);
    if (context.failure) {
        CHECK(context.failure->code == c::FailureCode::InvalidInput);
    } else {
        FAIL("The shared domain workspace accepted a nested call");
    }
    CHECK(memory.allocated_bytes() == payload * 2);
    check_hold_endpoint(*outer);
    const auto next =
        take(p::route_domain(source, hold_end_query, hold_policy, workspace, limits(), memory));
    CHECK(memory.allocated_bytes() == payload * 3);
    check_hold_endpoint(next);
}

TEST_CASE("Private PCHIP periodic time admits all four allocations before computing",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = prepare_interval_source();
    const p::RationalView cycle{1, hold_end_numerator, half_numerator};
    p::DomainWorkspace workspace;
    const auto plan = [&] {
        const FaultScope audit{0};
        return p::forecast_periodic_time(source, midpoint.time, cycle, limits());
    }();
    REQUIRE(plan);
    CHECK(attempts.load() == 0);
    auto exact = limits();
    exact.input_bits = 2;
    exact.work_bits = plan->arithmetic().integer_bits();
    exact.steps = plan->reserved_steps();
    exact.slots = plan->arithmetic().slots();
    const auto payload = plan->payload_bytes();
    const auto call = [&](c::ResourcePlan& memory, p::RationalLimits cap,
                          const std::stop_token& stop = {}) {
        return p::normalize_periodic_time(source, midpoint.time, cycle, workspace, cap, memory,
                                          stop, details());
    };
    c::ResourcePlan initial{payload};
    const auto old = [&] {
        const FaultScope audit{INT64_MAX};
        return call(initial, exact);
    }();
    REQUIRE(old);
    CHECK(attempts.load() == 4);
    CHECK(initial.allocated_bytes() == payload);

    std::uint64_t charged = 0;
    for (std::int64_t fault = 0; fault < 4; ++fault) {
        CAPTURE(fault);
        c::ResourcePlan memory{payload};
        const auto refused = [&] {
            const FaultScope audit{fault};
            return call(memory, exact);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(refused.error().details.model == details().model);
        CHECK(attempts.load() == static_cast<std::uint64_t>(fault + 1));
        if (fault == 0)
            CHECK(memory.allocated_bytes() == 0);
        else
            CHECK(memory.allocated_bytes() > charged);
        CHECK(memory.allocated_bytes() < payload);
        charged = memory.allocated_bytes();
    }
    std::array caps{exact, exact, exact, exact, exact};
    --caps[0].input_bits;
    --caps[1].work_bits;
    --caps[2].steps;
    --caps[3].slots;
    for (std::size_t i = 0; i < caps.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{i == 4 ? payload - 1 : payload};
        const auto refused = [&] {
            const FaultScope audit{0};
            return call(memory, caps[i]);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == 0);
        CHECK(memory.allocated_bytes() == 0);
    }

    c::ResourcePlan retried{payload * 2};
    std::stop_source stop;
    cancel_source = &stop;
    until_cancel = 4;
    allocation_callback = cancel_last;
    const auto cancelled = call(retried, exact, stop.get_token());
    allocation_callback = nullptr;
    cancel_source = nullptr;
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == c::FailureCode::Cancelled);
    CHECK(retried.allocated_bytes() == payload);
    const auto again = take(call(retried, exact));
    CHECK(retried.allocated_bytes() == payload * 2);
    const auto a = take(old->time());
    const auto b = take(again.time());
    CHECK(a.sign == midpoint.time.sign);
    CHECK(b.sign == a.sign);
    CHECK(std::ranges::equal(a.numerator, midpoint.time.numerator));
    CHECK(std::ranges::equal(a.denominator, midpoint.time.denominator));
    CHECK(std::ranges::equal(b.numerator, a.numerator));
    CHECK(std::ranges::equal(b.denominator, a.denominator));
}

namespace {
struct BoundaryReentry {
    const p::CoefficientBlock& source;
    p::BoundaryRequest request;
    p::BoundaryWorkspace& workspace;
    c::ResourcePlan& memory;
    std::optional<c::Failure> failure;
};
constinit BoundaryReentry* pending_boundary_reentry = nullptr;

void try_boundary_reentry() {
    auto& entry = *pending_boundary_reentry;
    const auto result = p::check_boundary(entry.source, entry.request, entry.workspace, limits(),
                                          entry.memory, {}, details());
    if (!result)
        entry.failure = result.error();
}
} // namespace

TEST_CASE("Private PCHIP boundary admits its entire arena and isolates faults",
          "[AE-07][N4-02][N4-03][N4-08][N4-09][negative]") {
    Fixture fixture;
    fixture.y.fill(fixture.y.front());
    p::CoefficientWorkspace coefficients;
    c::ResourcePlan source_memory{UINT64_MAX};
    const auto source =
        take(p::prepare_coefficients(fixture.view(), coefficients, limits(), source_memory));
    const p::BoundaryRequest request{fixture.t.front(), fixture.t.back(),
                                     p::PeriodicBoundary{fixture.t.back()}};
    const auto plan = [&] {
        const FaultScope audit{0};
        return p::forecast_boundary(source, request, limits());
    }();
    REQUIRE(plan);
    CHECK(attempts.load() == 0);
    auto exact = limits();
    exact.input_bits = 2;
    exact.work_bits = plan->arithmetic().integer_bits();
    exact.steps = plan->reserved_steps();
    exact.slots = plan->arithmetic().slots();
    const auto payload = plan->payload_bytes();
    constexpr std::int64_t allocations = 3;

    p::BoundaryWorkspace workspace;
    const auto call = [&](c::ResourcePlan& memory, p::RationalLimits cap,
                          const std::stop_token& stop = {}) {
        return p::check_boundary(source, request, workspace, cap, memory, stop, details());
    };
    c::ResourcePlan initial{payload};
    const auto result = [&] {
        const FaultScope audit{INT64_MAX};
        return call(initial, exact);
    }();
    REQUIRE(result);
    CHECK(attempts.load() == static_cast<std::uint64_t>(allocations));
    CHECK(initial.allocated_bytes() == payload);
    CHECK(result->slots == exact.slots);
    CHECK(result->steps <= exact.steps);

    std::uint64_t previous_charge = 0;
    for (std::int64_t fault = 0; fault < allocations; ++fault) {
        CAPTURE(fault);
        c::ResourcePlan memory{payload};
        const auto rejected = [&] {
            const FaultScope audit{fault};
            return call(memory, exact);
        }();
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == c::FailureCode::ResourceExceeded);
        CHECK(rejected.error().details.component == details().component);
        CHECK(attempts.load() == static_cast<std::uint64_t>(fault + 1));
        CHECK(memory.allocated_bytes() < payload);
        if (fault == 0)
            CHECK(memory.allocated_bytes() == 0);
        else
            CHECK(memory.allocated_bytes() > previous_charge);
        previous_charge = memory.allocated_bytes();
    }

    std::array caps{exact, exact, exact, exact, exact};
    --caps[0].input_bits;
    --caps[1].work_bits;
    --caps[2].steps;
    --caps[3].slots;
    for (std::size_t i = 0; i < caps.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{i == 4 ? payload - 1 : payload};
        const auto rejected = [&] {
            const FaultScope audit{0};
            return call(memory, caps[i]);
        }();
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == c::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == 0);
        CHECK(memory.allocated_bytes() == 0);
    }

    c::ResourcePlan retried{payload * 2};
    std::stop_source stop;
    cancel_source = &stop;
    until_cancel = static_cast<unsigned>(allocations);
    allocation_callback = cancel_last;
    const auto cancellation = call(retried, exact, stop.get_token());
    allocation_callback = nullptr;
    cancel_source = nullptr;
    REQUIRE_FALSE(cancellation);
    CHECK(cancellation.error().code == c::FailureCode::Cancelled);
    CHECK(cancellation.error().details == details());
    CHECK(retried.allocated_bytes() == payload);

    auto malformed = request;
    malformed.end = {};
    c::ResourcePlan unused{payload};
    const auto stopped_before_work = [&] {
        const FaultScope audit{0};
        return call(unused, exact, stop.get_token());
    }();
    REQUIRE_FALSE(stopped_before_work);
    CHECK(stopped_before_work.error().code == c::FailureCode::Cancelled);
    CHECK(attempts.load() == 0);
    const auto malformed_before_stop = [&] {
        const FaultScope audit{0};
        return p::check_boundary(source, malformed, workspace, exact, unused, stop.get_token());
    }();
    REQUIRE_FALSE(malformed_before_stop);
    CHECK(malformed_before_stop.error().code == c::FailureCode::InvalidInput);
    CHECK(attempts.load() == 0);
    CHECK(unused.allocated_bytes() == 0);

    BoundaryReentry entry{source, request, workspace, retried, {}};
    pending_boundary_reentry = &entry;
    allocation_callback = try_boundary_reentry;
    const auto after_failure = call(retried, exact);
    allocation_callback = nullptr;
    pending_boundary_reentry = nullptr;
    if (entry.failure) {
        CHECK(entry.failure->code == c::FailureCode::InvalidInput);
        CHECK(entry.failure->details == details());
    } else {
        FAIL("The boundary workspace accepted a nested call");
    }
    REQUIRE(after_failure);
    CHECK(after_failure->slots == result->slots);
    CHECK(after_failure->steps == result->steps);
    CHECK(retried.allocated_bytes() == payload * 2);
    for (std::uint64_t knot = 0; knot < source.nodes(); ++knot) {
        const auto value = take(source.value(knot, 0));
        CHECK(value.sign == 0);
        CHECK(value.numerator.empty());
        CHECK(std::ranges::equal(value.denominator, fixture.one));
    }
}

namespace {
constexpr p::PchipInputStamp preparation_stamp{p::PchipIssuerId{17}, p::PchipInputGeneration{3}};
constexpr std::array<uvlas::numerics::SiDimension, 1> preparation_dimensions{
    {{{1, 0, -1, 0, 0, 0, 0}}}};

struct PreparationReentry {
    const p::PreparationRequest& request;
    p::PreparationWorkspace& workspace;
    c::ResourcePlan& memory;
    std::optional<c::Failure> failure;
};

constinit PreparationReentry* pending_preparation = nullptr;

void try_preparation_reentry() {
    auto& item = *pending_preparation;
    const auto result = p::make_preparation_input(item.request, item.workspace, limits(),
                                                  item.memory, {}, details());
    if (!result)
        item.failure = result.error();
}
} // namespace

TEST_CASE("Private PCHIP preparation input preflights storage and retains old owners",
          "[AE-07][N4-02][N4-03][N4-08][N4-09][negative]") {
    const Fixture fixture;
    const p::PreparationRequest request{fixture.view(),
                                        preparation_dimensions,
                                        {fixture.t.front(), fixture.t.back(), p::RejectBoundary{}},
                                        {},
                                        preparation_stamp};
    p::PreparationWorkspace workspace;
    const auto plan = [&] {
        const FaultScope audit{0};
        return p::forecast_preparation_input(request, limits());
    }();
    REQUIRE(plan);
    CHECK(attempts.load() == 0);

    auto exact = limits();
    exact.input_bits = 3;
    exact.work_bits = plan->arithmetic().integer_bits();
    exact.slots = plan->arithmetic().slots();
    exact.steps = plan->reserved_steps();
    const auto payload = plan->payload_bytes();
    constexpr std::int64_t allocations = 7;
    const auto make = [&](c::ResourcePlan& memory, p::RationalLimits cap,
                          const std::stop_token& stop = {}) {
        return p::make_preparation_input(request, workspace, cap, memory, stop, details());
    };
    c::ResourcePlan charged{payload};
    const auto initial = [&] {
        const FaultScope audit{INT64_MAX};
        return make(charged, exact);
    }();
    REQUIRE(initial);
    CHECK(attempts.load() == static_cast<std::uint64_t>(allocations));
    CHECK(charged.allocated_bytes() == payload);
    CHECK(initial->used_steps() <= exact.steps);
    const auto before = take(initial->view());

    std::uint64_t previous = 0;
    for (std::int64_t fault = 0; fault < allocations; ++fault) {
        CAPTURE(fault);
        c::ResourcePlan memory{payload};
        const auto refused = [&] {
            const FaultScope audit{fault};
            return make(memory, exact);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(refused.error().details.component == details().component);
        CHECK(attempts.load() == static_cast<std::uint64_t>(fault + 1));
        CHECK(memory.allocated_bytes() < payload);
        if (fault == 0)
            CHECK(memory.allocated_bytes() == 0);
        else
            CHECK(memory.allocated_bytes() > previous);
        previous = memory.allocated_bytes();
    }

    std::array caps{exact, exact, exact, exact, exact};
    --caps[0].input_bits;
    --caps[1].work_bits;
    --caps[2].steps;
    --caps[3].slots;
    for (std::size_t i = 0; i < caps.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{i == 4 ? payload - 1 : payload};
        const auto refused = [&] {
            const FaultScope audit{0};
            return make(memory, caps[i]);
        }();
        REQUIRE_FALSE(refused);
        CHECK(refused.error().code == c::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == 0);
        CHECK(memory.allocated_bytes() == 0);
    }

    c::ResourcePlan retried{payload * 2};
    std::stop_source stop;
    cancel_source = &stop;
    until_cancel = static_cast<unsigned>(allocations);
    allocation_callback = cancel_last;
    const auto cancelled = make(retried, exact, stop.get_token());
    allocation_callback = nullptr;
    cancel_source = nullptr;
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == c::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details());
    CHECK(retried.allocated_bytes() == payload);

    PreparationReentry context{request, workspace, retried, {}};
    pending_preparation = &context;
    allocation_callback = try_preparation_reentry;
    const auto retry = make(retried, exact);
    allocation_callback = nullptr;
    pending_preparation = nullptr;
    REQUIRE(retry);
    if (context.failure) {
        CHECK(context.failure->code == c::FailureCode::InvalidInput);
        CHECK(context.failure->details == details());
    } else {
        FAIL("The preparation workspace accepted reentry");
    }
    CHECK(retried.allocated_bytes() == payload * 2);
    const auto compare = [&] {
        const FaultScope audit{0};
        return p::match_preparation_input(*initial, *retry, UINT64_MAX);
    }();
    CHECK(compare.has_value());
    CHECK(attempts.load() == 0);
    const auto after = take(initial->view());
    CHECK(before.samples.times.data() == after.samples.times.data());
    CHECK(before.stamp == after.stamp);
    CHECK(std::ranges::equal(after.samples.values.back().numerator, fixture.four));
}

TEST_CASE("PCHIP curve numerical forecast never allocates before aggregate admission",
          "[AE-07][N4-02][N4-09][negative]") {
    const Fixture fixture;
    const p::PreparationRequest request{fixture.view(),
                                        preparation_dimensions,
                                        {fixture.t.front(), fixture.t.back(), p::RejectBoundary{}},
                                        {},
                                        preparation_stamp};
    c::ResourcePlan immutable_memory{UINT64_MAX};
    p::PreparationWorkspace workspace;
    const auto input =
        take(p::make_preparation_input(request, workspace, limits(), immutable_memory));
    const auto planned = [&] {
        const FaultScope audit{0};
        return p::forecast_curve_numerics(input, limits());
    }();
    REQUIRE(planned);
    CHECK(attempts.load() == 0);
    auto small = limits();
    small.steps = planned->reserved_steps() - 1;
    CHECK(small.steps >= planned->coefficients().reserved_steps());
    CHECK(small.steps >= planned->boundary().reserved_steps());
    const auto rejected = [&] {
        const FaultScope audit{0};
        return p::forecast_curve_numerics(input, small);
    }();
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == c::FailureCode::ResourceExceeded);
    CHECK(attempts.load() == 0);
}

namespace {
namespace linear_memory = uvlas::numerics::linear;

struct ShareReentry {
    p::PchipPreparationInput& original;
    p::SharedPreparationWorkspace& workspace;
    p::SharedPreparationLimits limits;
    c::ResourcePlan& memory;
    std::optional<c::Failure> failure;
};

constinit ShareReentry* pending_share = nullptr;

void try_share_reentry() {
    auto& call = *pending_share;
    const auto result = p::share_preparation_input(std::move(call.original), call.workspace,
                                                   call.limits, call.memory, {}, details());
    if (!result)
        call.failure = result.error();
}

} // namespace

TEST_CASE("PCHIP shared preparation admits resources before moving the complete input",
          "[AE-07][N4-02][N4-03][N4-08][N4-09][negative]") {
    const Fixture fixture;
    const p::PreparationRequest request{fixture.view(),
                                        preparation_dimensions,
                                        {fixture.t.front(), fixture.t.back(), p::RejectBoundary{}},
                                        {},
                                        preparation_stamp};
    c::ResourcePlan preparation_memory{UINT64_MAX};
    p::PreparationWorkspace preparation_workspace;
    const auto new_input = [&] {
        return take(p::make_preparation_input(request, preparation_workspace, limits(),
                                              preparation_memory));
    };
    auto input = new_input();
    const auto expected = new_input();
    const auto initial_charge = preparation_memory.allocated_bytes();
    const p::SharedPreparationLimits unlimited{linear_memory::EstimateAdmission{UINT64_MAX},
                                               UINT64_MAX};
    const auto plan = [&] {
        const FaultScope audit{0};
        return p::forecast_shared_preparation(input, unlimited);
    }();
    REQUIRE(plan);
    CHECK(attempts.load() == 0);
    const auto payload = plan->payload_bytes();
    const auto meta = plan->metadata_estimate_bytes();
    const auto steps = plan->reserved_steps();
    REQUIRE(payload > 0);
    REQUIRE(meta > 0);
    REQUIRE(steps > 0);
    const p::SharedPreparationLimits exact{linear_memory::EstimateAdmission{meta}, steps};
    p::SharedPreparationWorkspace workspace;
    const auto share = [&](c::ResourcePlan& memory, p::SharedPreparationLimits cap,
                           const std::stop_token& stop = {}) {
        return p::share_preparation_input(std::move(input), workspace, cap, memory, stop,
                                          details());
    };
    const auto unchanged = [&] {
        return p::match_preparation_input(input, expected, UINT64_MAX).has_value();
    };

    const std::array rejected_limits{
        p::SharedPreparationLimits{linear_memory::RequireProvenBound{UINT64_MAX}, steps},
        p::SharedPreparationLimits{linear_memory::EstimateAdmission{meta - 1}, steps},
        p::SharedPreparationLimits{linear_memory::EstimateAdmission{meta}, steps - 1}, exact};
    for (std::size_t i = 0; i < rejected_limits.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{i == 3 ? payload - 1 : payload};
        const auto rejected = [&] {
            const FaultScope audit{0};
            return share(memory, rejected_limits[i]);
        }();
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code ==
              (i == 0 ? c::FailureCode::Unsupported : c::FailureCode::ResourceExceeded));
        CHECK(rejected.error().details.component == details().component);
        CHECK(attempts.load() == 0);
        CHECK(memory.allocated_bytes() == 0);
        CHECK(unchanged());
    }

    for (const std::int64_t allocation : {std::int64_t{0}, std::int64_t{1}}) {
        CAPTURE(allocation);
        c::ResourcePlan memory{payload};
        const auto rejected = [&] {
            const FaultScope audit{allocation};
            return share(memory, exact);
        }();
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == c::FailureCode::ResourceExceeded);
        CHECK(rejected.error().details.component == details().component);
        CHECK(attempts.load() == static_cast<std::uint64_t>(allocation + 1));
        CHECK(memory.allocated_bytes() == (allocation == 0 ? 0 : payload));
        CHECK(unchanged());
    }

    for (const unsigned allocation : {1U, 2U}) {
        CAPTURE(allocation);
        c::ResourcePlan memory{payload};
        std::stop_source stop;
        cancel_source = &stop;
        until_cancel = allocation;
        allocation_callback = cancel_last;
        const auto rejected = share(memory, exact, stop.get_token());
        allocation_callback = nullptr;
        cancel_source = nullptr;
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == c::FailureCode::Cancelled);
        CHECK(rejected.error().details == details());
        CHECK(memory.allocated_bytes() == payload);
        CHECK(unchanged());
    }

    // Verify that reentry attempted by each allocation leaves the complete input untouched.
    for (const unsigned allocation : {1U, 2U}) {
        CAPTURE(allocation);
        auto nested_input = new_input();
        c::ResourcePlan memory{payload};
        ShareReentry context{nested_input, workspace, exact, memory, {}};
        pending_share = &context;
        if (allocation == 1) {
            allocation_callback = try_share_reentry;
        } else {
            allocation_callback = [] { allocation_callback = try_share_reentry; };
        }
        const auto handle =
            p::share_preparation_input(std::move(nested_input), workspace, exact, memory);
        allocation_callback = nullptr;
        pending_share = nullptr;
        REQUIRE(handle);
        if (context.failure) {
            CHECK(context.failure->code == c::FailureCode::InvalidInput);
            CHECK(context.failure->details == details());
        } else {
            FAIL("The shared-input workspace accepted reentry");
        }
        CHECK(memory.allocated_bytes() == payload);
        CHECK(p::match_preparation_input(take(handle->input()).get(), expected, UINT64_MAX));
    }
    CHECK(preparation_memory.allocated_bytes() > initial_charge);
    CHECK(unchanged());

    c::ResourcePlan final_memory{payload};
    const auto accepted = [&] {
        const FaultScope audit{INT64_MAX};
        return share(final_memory, exact);
    }();
    REQUIRE(accepted);
    CHECK(attempts.load() == 2);
    CHECK(final_memory.allocated_bytes() == payload);
    CHECK(p::match_preparation_input(take(accepted->input()).get(), expected, UINT64_MAX));
}

namespace {
constinit unsigned until_curve_reentry = 0;
constinit const void* curve_reentry_context = nullptr;
constinit void (*curve_reentry_invoke)(const void*) = nullptr;

void reenter_curve() {
    if (--until_curve_reentry == 0)
        curve_reentry_invoke(curve_reentry_context);
    else
        allocation_callback = reenter_curve;
}

template <class Call, class Unchanged>
void verify_curve_faults(const p::CurveAdmission& plan, p::CurveLimits cap,
                         unsigned planned_allocations, const Call& call,
                         const Unchanged& unchanged) {
    cap.arithmetic.steps = plan.reserved_steps;
    cap.arithmetic.work_bits = plan.maximum_integer_bits;
    cap.metadata = linear_memory::EstimateAdmission{plan.metadata_estimate_bytes};
    const auto exact = [&] {
        c::ResourcePlan memory{plan.payload_bytes};
        auto result = [&] {
            const FaultScope audit{INT64_MAX};
            return call(memory, cap, {});
        }();
        CHECK(memory.allocated_bytes() == plan.payload_bytes);
        return result;
    }();
    REQUIRE(exact);
    const auto actual_allocations = attempts.load();
    REQUIRE(actual_allocations == planned_allocations);
    CHECK(unchanged());

    auto work_less = cap;
    --work_less.arithmetic.steps;
    auto opaque_less = cap;
    opaque_less.metadata = linear_memory::EstimateAdmission{plan.metadata_estimate_bytes - 1};
    auto proven = cap;
    proven.metadata = linear_memory::RequireProvenBound{UINT64_MAX};
    auto input_one_bit = cap;
    input_one_bit.arithmetic.input_bits = 1;
    auto bit_less = cap;
    --bit_less.arithmetic.work_bits;
    const std::array rejected{work_less, opaque_less, proven, input_one_bit, bit_less, cap};
    for (std::size_t i = 0; i < rejected.size(); ++i) {
        CAPTURE(i);
        c::ResourcePlan memory{i == rejected.size() - 1 ? plan.payload_bytes - 1
                                                        : plan.payload_bytes};
        const auto refusal = [&] {
            const FaultScope audit{0};
            return call(memory, rejected[i], {});
        }();
        REQUIRE_FALSE(refusal);
        CHECK(refusal.error().code ==
              (i == 2 ? c::FailureCode::Unsupported : c::FailureCode::ResourceExceeded));
        CHECK(refusal.error().operation == c::Operation::EvaluateAero);
        CHECK(refusal.error().details.component == details().component);
        CHECK(attempts.load() == 0);
        CHECK(memory.allocated_bytes() == 0);
        CHECK(unchanged());
    }
    for (unsigned position = 1; position <= actual_allocations; ++position) {
        CAPTURE(position);
        c::ResourcePlan fault_memory{plan.payload_bytes};
        const auto fault = [&] {
            const FaultScope audit{static_cast<std::int64_t>(position) - 1};
            return call(fault_memory, cap, {});
        }();
        REQUIRE_FALSE(fault);
        CHECK(fault.error().code == c::FailureCode::ResourceExceeded);
        CHECK(fault.error().details.component == details().component);
        CHECK(attempts.load() == position);
        CHECK(fault_memory.allocated_bytes() <= plan.payload_bytes);
        CHECK(unchanged());

        // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 declares request_stop non-const.
        std::stop_source stop;
        cancel_source = &stop;
        until_cancel = position;
        allocation_callback = cancel_last;
        c::ResourcePlan cancel_memory{plan.payload_bytes};
        const auto cancellation = call(cancel_memory, cap, stop.get_token());
        allocation_callback = nullptr;
        cancel_source = nullptr;
        REQUIRE_FALSE(cancellation);
        CHECK(cancellation.error().code == c::FailureCode::Cancelled);
        CHECK(cancellation.error().details == details());
        CHECK(cancel_memory.allocated_bytes() <= plan.payload_bytes);
        CHECK(unchanged());

        c::ResourcePlan nested_memory{plan.payload_bytes};
        std::optional<c::Failure> refusal;
        const auto reenter = [&] {
            auto denied = call(nested_memory, cap, {});
            if (!denied)
                refusal = denied.error();
        };
        const auto invoke = []<class F>(const F& operation) {
            curve_reentry_context = &operation;
            curve_reentry_invoke = [](const void* context) { (*static_cast<const F*>(context))(); };
        };
        invoke(reenter);
        until_curve_reentry = position;
        allocation_callback = reenter_curve;
        c::ResourcePlan outer_memory{plan.payload_bytes};
        const auto outer = call(outer_memory, cap, {});
        allocation_callback = nullptr;
        curve_reentry_context = nullptr;
        curve_reentry_invoke = nullptr;
        REQUIRE(outer);
        if (!refusal)
            FAIL("Same curve workspace accepted reentry");
        CHECK(refusal->code == c::FailureCode::InvalidInput);
        CHECK(refusal->details == details());
        CHECK(nested_memory.allocated_bytes() == 0);
        CHECK(unchanged());
    }
}
} // namespace

TEST_CASE("PCHIP curve entrypoints admit all resources and reject faults without publication",
          "[AE-07][N4-02][N4-03][N4-08][N4-09][negative]") {
    const Fixture fixture;
    const p::PreparationRequest request{fixture.view(),
                                        preparation_dimensions,
                                        {fixture.t.front(), fixture.t.back(), p::RejectBoundary{}},
                                        {},
                                        preparation_stamp};
    const p::CurveLimits cap{limits(), linear_memory::EstimateAdmission{UINT64_MAX}};
    p::CurveWorkspace workspace;
    p::PreparationWorkspace input_workspace;
    p::SharedPreparationWorkspace sharing;
    c::ResourcePlan parents{UINT64_MAX};
    auto independent = take(p::make_preparation_input(request, input_workspace, limits(), parents));
    const auto current = take(p::share_preparation_input(
        std::move(independent), sharing, {cap.metadata, cap.arithmetic.steps}, parents));
    const auto curve = take(p::prepare_pchip_curve(request, workspace, cap, parents));
    const auto saving = take(p::snapshot_pchip_curve(curve, workspace, cap, parents));
    const auto saved = take(saving.view());
    const p::PchipQuery query{
        {1, fixture.one, half_denominator}, p::JetOrder::Second, p::RequestedSide::Right};
    const auto older_jet =
        take(p::evaluate_pchip_curve(curve, current, query, workspace, cap, parents));
    const auto older_value = take(older_jet.derivative(0, p::JetOrder::Value));
    const auto older_coefficient = take(curve.coefficient(1, 0, 0));
    const auto ledger = parents.allocated_bytes();

    const auto unchanged = [&] {
        const auto present = older_jet.derivative(0, p::JetOrder::Value);
        const auto coefficient = curve.coefficient(1, 0, 0);
        return present && coefficient &&
               present->denominator.data() == older_value.denominator.data() &&
               coefficient->numerator.data() == older_coefficient.numerator.data() &&
               std::ranges::equal(coefficient->numerator, fixture.one) &&
               saving.view().has_value() && parents.allocated_bytes() == ledger;
    };

    SECTION("prepare from immutable borrows") {
        const auto plan = take(p::forecast_pchip_prepare(request, cap));
        verify_curve_faults(
            plan, cap, 21,
            [&](c::ResourcePlan& memory, p::CurveLimits budget, const std::stop_token& stop) {
                return p::prepare_pchip_curve(request, workspace, budget, memory, stop, details());
            },
            unchanged);

        auto corrupted = request;
        auto values = fixture.y;
        values[1] = {1, fixture.two, fixture.two}; // Shape-valid but noncanonical 2/2.
        corrupted.samples.values = values;
        const auto bound = take(p::forecast_pchip_prepare(corrupted, cap));
        c::ResourcePlan payload{bound.payload_bytes};
        const auto invalid = [&] {
            const FaultScope audit{INT64_MAX};
            return p::prepare_pchip_curve(corrupted, workspace, cap, payload, {}, details());
        }();
        REQUIRE_FALSE(invalid);
        CHECK(invalid.error().code == c::FailureCode::InvalidInput);
        CHECK(attempts.load() == 21);
        CHECK(unchanged());
    }
    SECTION("evaluate against independently owned expected input") {
        const auto plan = take(p::forecast_pchip_evaluation(curve, current, query, cap));
        verify_curve_faults(
            plan, cap, 16,
            [&](c::ResourcePlan& memory, p::CurveLimits budget, const std::stop_token& stop) {
                return p::evaluate_pchip_curve(curve, current, query, workspace, budget, memory,
                                               stop, details());
            },
            unchanged);

        auto corrupted = query;
        corrupted.time = {1, fixture.two, fixture.two}; // Canonicality follows all allocation.
        const auto bound = take(p::forecast_pchip_evaluation(curve, current, corrupted, cap));
        c::ResourcePlan payload{bound.payload_bytes};
        const auto invalid = [&] {
            const FaultScope audit{INT64_MAX};
            return p::evaluate_pchip_curve(curve, current, corrupted, workspace, cap, payload, {},
                                           details());
        }();
        REQUIRE_FALSE(invalid);
        CHECK(invalid.error().code == c::FailureCode::InvalidInput);
        CHECK(attempts.load() == 16);
        CHECK(unchanged());
    }
    SECTION("restore against independently owned expected input") {
        const auto plan = take(p::forecast_pchip_restore(saved, current, cap));
        verify_curve_faults(
            plan, cap, 14,
            [&](c::ResourcePlan& memory, p::CurveLimits budget, const std::stop_token& stop) {
                return p::restore_pchip_curve(saved, current, workspace, budget, memory, stop,
                                              details());
            },
            unchanged);

        auto corrupted = saved;
        std::array<p::RationalView, 8> coefficients;
        std::ranges::copy(saved.coefficients, coefficients.begin());
        coefficients[0] = {1, fixture.two, fixture.two};
        corrupted.coefficients = coefficients;
        const auto bound = take(p::forecast_pchip_restore(corrupted, current, cap));
        c::ResourcePlan payload{bound.payload_bytes};
        const auto invalid = [&] {
            const FaultScope audit{INT64_MAX};
            return p::restore_pchip_curve(corrupted, current, workspace, cap, payload, {},
                                          details());
        }();
        REQUIRE_FALSE(invalid);
        CHECK(invalid.error().code == c::FailureCode::InvalidInput);
        CHECK(attempts.load() == 14);
        CHECK(unchanged());
    }
}

TEST_CASE("Periodic exact PCHIP reserves both automatic seam candidates before normalization",
          "[AE-07][N4-02][N4-03][N4-08][N4-09][negative]") {
    const Fixture fixture;
    const std::array<p::RationalView, 3> samples{fixture.t[0], fixture.t[0], fixture.t[0]};
    const p::PreparationRequest request{
        {fixture.t, samples, 1},
        preparation_dimensions,
        {fixture.t.front(), fixture.t.back(), p::PeriodicBoundary{fixture.t.back()}},
        {},
        preparation_stamp};
    const p::CurveLimits cap{limits(), linear_memory::EstimateAdmission{UINT64_MAX}};
    p::CurveWorkspace workspace;
    p::PreparationWorkspace input_workspace;
    p::SharedPreparationWorkspace sharing;
    c::ResourcePlan parents{UINT64_MAX};
    auto owned = take(p::make_preparation_input(request, input_workspace, cap.arithmetic, parents));
    const auto current = take(p::share_preparation_input(
        std::move(owned), sharing, {cap.metadata, cap.arithmetic.steps}, parents));
    const auto curve = take(p::prepare_pchip_curve(request, workspace, cap, parents));
    const p::PchipQuery seam{{-1, fixture.two, fixture.one}, p::JetOrder::Second};
    const auto plan = take(p::forecast_pchip_evaluation(curve, current, seam, cap));
    const auto before = parents.allocated_bytes();
    verify_curve_faults(
        plan, cap, 24,
        [&](c::ResourcePlan& memory, p::CurveLimits allowed, const std::stop_token& stop) {
            return p::evaluate_pchip_curve(curve, current, seam, workspace, allowed, memory, stop,
                                           details());
        },
        [&] {
            const auto value = curve.coefficient(1, 0, 0);
            return value && value->sign == 0 && parents.allocated_bytes() == before;
        });
}
