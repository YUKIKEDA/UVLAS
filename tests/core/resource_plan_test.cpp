#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <utility>
#include <uvlas/core/resource_plan.hpp>

namespace {
using namespace uvlas::core;
constexpr Operation operation = Operation::PrepareModel;

// Per-executable probe: observe allocation itself, not just element construction.
struct AllocationProbe {
    static inline int allocations = 0;
    static inline int constructions = 0;
    int value = 0;
    AllocationProbe() noexcept {
        ++constructions;
    }
    static void* operator new[](std::size_t bytes) {
        ++allocations;
        return ::operator new(bytes);
    }
    static void operator delete[](void* pointer) noexcept {
        ::operator delete(pointer);
    }
};

struct ThrowsBadAlloc {
    ThrowsBadAlloc() {
        throw std::bad_alloc{};
    }
};
struct ThrowsOther {
    ThrowsOther() {
        throw std::runtime_error{"element construction failed"};
    }
};

enum class NestedOutcome : std::uint8_t { Success, BadAlloc, OtherException };

struct NestedProbe {
    static inline ResourcePlan* plan = nullptr;
    static inline NestedOutcome outcome = NestedOutcome::Success;
    static inline std::optional<Failure> child_failure;
    static inline std::uint64_t remaining_before_child = 0;
    PlannedBuffer<AllocationProbe> child;

    NestedProbe() {
        remaining_before_child = plan->remaining_bytes();
        auto result = plan->allocate<AllocationProbe>(1, operation);
        if (result) {
            child = std::move(*result);
        } else {
            child_failure = result.error();
        }
        if (outcome == NestedOutcome::BadAlloc) {
            throw std::bad_alloc{};
        }
        if (outcome == NestedOutcome::OtherException) {
            throw std::runtime_error{"outer construction failed"};
        }
    }
};

// Restore the test callback context even when a REQUIRE aborts a test case.
struct NestedProbeContext {
    explicit NestedProbeContext(ResourcePlan& plan,
                                NestedOutcome outcome = NestedOutcome::Success) noexcept {
        NestedProbe::plan = &plan;
        NestedProbe::outcome = outcome;
        NestedProbe::child_failure.reset();
        AllocationProbe::allocations = 0;
        AllocationProbe::constructions = 0;
    }
    ~NestedProbeContext() {
        NestedProbe::plan = nullptr;
    }
};
} // namespace

TEST_CASE("ResourcePlan checks the cumulative budget before any allocation", "[negative][N4-02]") {
    AllocationProbe::allocations = 0;
    AllocationProbe::constructions = 0;
    ResourcePlan plan{3 * sizeof(AllocationProbe)};
    auto first = plan.allocate<AllocationProbe>(2, operation);
    REQUIRE(first);
    CHECK(first->size() == 2);
    CHECK(AllocationProbe::allocations == 1);
    CHECK(AllocationProbe::constructions == 2);
    CHECK(plan.allocated_bytes() == 2 * sizeof(AllocationProbe));
    CHECK(plan.remaining_bytes() == sizeof(AllocationProbe));

    const auto rejected = plan.allocate<AllocationProbe>(2, operation);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == FailureCode::ResourceExceeded);
    CHECK(rejected.error().operation == operation);
    if (const auto& resources = rejected.error().details.resources; resources) {
        CHECK(resources->required == 4 * sizeof(AllocationProbe));
        CHECK(resources->observed == 2 * sizeof(AllocationProbe));
    } else {
        FAIL("Missing resource diagnostic");
    }
    CHECK(AllocationProbe::allocations == 1);
    CHECK(AllocationProbe::constructions == 2);
    CHECK(plan.allocated_bytes() == 2 * sizeof(AllocationProbe));

    const auto last = plan.allocate<AllocationProbe>(1, operation);
    REQUIRE(last);
    CHECK(plan.remaining_bytes() == 0);
    CHECK(plan.allocated_bytes() == plan.byte_limit());
    CHECK(AllocationProbe::allocations == 2);
}

TEST_CASE("Empty buffers and invalid extents do not allocate or consume budget",
          "[negative][N4-02]") {
    AllocationProbe::allocations = 0;
    ResourcePlan plan{0};
    const auto empty = plan.allocate<AllocationProbe>(0, operation);
    REQUIRE(empty);
    CHECK(empty->size() == 0);
    CHECK(empty->values().empty());
    CHECK(AllocationProbe::allocations == 0);
    const auto negative = plan.allocate<AllocationProbe>(-1, operation);
    REQUIRE_FALSE(negative);
    CHECK(negative.error().code == FailureCode::InvalidInput);
    const auto over_budget = plan.allocate<AllocationProbe>(1, operation);
    REQUIRE_FALSE(over_budget);
    CHECK(over_budget.error().code == FailureCode::ResourceExceeded);
    CHECK(AllocationProbe::allocations == 0);
    CHECK(plan.allocated_bytes() == 0);
}

TEST_CASE("Buffer ownership survives plan destruction and moves leave empty sources",
          "[unit][N4-02]") {
    PlannedBuffer<double> owned;
    {
        ResourcePlan plan{24};
        auto allocated = plan.allocate<double>(3, operation);
        REQUIRE(allocated);
        CHECK(allocated->values()[0] == 0.0);
        CHECK(allocated->values()[1] == 0.0);
        CHECK(allocated->values()[2] == 0.0);
        allocated->values()[1] = 42.0;
        owned = std::move(*allocated);
        CHECK(allocated->size() == 0);
        CHECK(allocated->values().empty());
        CHECK(plan.allocated_bytes() == 24);
    }
    CHECK(owned.values()[1] == 42.0);
    PlannedBuffer<double> moved{std::move(owned)};
    // NOLINTNEXTLINE(bugprone-use-after-move): the public move contract resets size.
    CHECK(owned.size() == 0);
    CHECK(owned.values().empty());
    CHECK(moved.values()[1] == 42.0);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<PlannedBuffer<double>>);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<ResourcePlan>);
    STATIC_REQUIRE_FALSE(std::is_move_constructible_v<ResourcePlan>);
}

TEST_CASE("Failed element construction leaves resource accounting unchanged", "[negative][N4-02]") {
    ResourcePlan plan{64};
    const auto allocation_failure = plan.allocate<ThrowsBadAlloc>(1, operation);
    REQUIRE_FALSE(allocation_failure);
    CHECK(allocation_failure.error().code == FailureCode::ResourceExceeded);
    CHECK(plan.allocated_bytes() == 0);
    const auto other_failure = plan.allocate<ThrowsOther>(1, operation);
    REQUIRE_FALSE(other_failure);
    CHECK(other_failure.error().code == FailureCode::InvalidInput);
    CHECK(plan.allocated_bytes() == 0);
    const auto success = plan.allocate<double>(1, operation);
    REQUIRE(success);
    CHECK(plan.allocated_bytes() == 8);
}

TEST_CASE("Destroying a buffer cannot refund an incomplete preparation plan", "[unit][N4-02]") {
    ResourcePlan plan{8};
    {
        const auto first = plan.allocate<double>(1, operation);
        REQUIRE(first);
    }
    CHECK(plan.remaining_bytes() == 0);
    const auto second = plan.allocate<double>(1, operation);
    REQUIRE_FALSE(second);
    CHECK(second.error().code == FailureCode::ResourceExceeded);
}

TEST_CASE("An overflowing cumulative byte demand fails without allocation", "[negative][N4-02]") {
    ResourcePlan plan{2};
    const auto small = plan.allocate<std::uint16_t>(1, operation);
    REQUIRE(small);
    const auto impossible =
        plan.allocate<std::uint16_t>(std::numeric_limits<Index>::max(), operation);
    REQUIRE_FALSE(impossible);
    CHECK(impossible.error().code == FailureCode::ResourceExceeded);
    CHECK_FALSE(impossible.error().details.resources.has_value());
    CHECK(plan.allocated_bytes() == 2);
}

TEST_CASE("Nested allocation cannot reuse the outer construction budget", "[negative][N4-02]") {
    ResourcePlan plan{sizeof(NestedProbe)};
    const NestedProbeContext context{plan};
    const auto outer = plan.allocate<NestedProbe>(1, operation);
    REQUIRE(outer);
    CHECK(plan.allocated_bytes() <= plan.byte_limit());
    CHECK(plan.remaining_bytes() <= plan.byte_limit());
    CHECK(NestedProbe::remaining_before_child == 0);
    if (const auto& failure = NestedProbe::child_failure; failure) {
        CHECK(failure->code == FailureCode::ResourceExceeded);
        if (const auto& resources = failure->details.resources; resources) {
            CHECK(resources->required == sizeof(NestedProbe) + sizeof(AllocationProbe));
            CHECK(resources->observed == sizeof(NestedProbe));
        } else {
            FAIL("Nested rejection lost its reserved-byte diagnostic");
        }
    } else {
        FAIL("Nested allocation reused the outer reservation");
    }
    CHECK(AllocationProbe::allocations == 0);
    CHECK(AllocationProbe::constructions == 0);
    CHECK(outer->values()[0].child.size() == 0);
    CHECK(plan.allocated_bytes() == plan.byte_limit());
    CHECK(plan.remaining_bytes() == 0);
    const auto following = plan.allocate<unsigned char>(1, operation);
    REQUIRE_FALSE(following);
    CHECK(following.error().code == FailureCode::ResourceExceeded);
}

TEST_CASE("Nested allocations can fill the combined budget exactly", "[unit][N4-02]") {
    ResourcePlan plan{sizeof(NestedProbe) + sizeof(AllocationProbe)};
    const NestedProbeContext context{plan};
    const auto outer = plan.allocate<NestedProbe>(1, operation);
    REQUIRE(outer);
    CHECK_FALSE(NestedProbe::child_failure.has_value());
    CHECK(NestedProbe::remaining_before_child == sizeof(AllocationProbe));
    CHECK(outer->values()[0].child.size() == 1);
    CHECK(AllocationProbe::allocations == 1);
    CHECK(plan.allocated_bytes() == plan.byte_limit());
    CHECK(plan.remaining_bytes() == 0);
}

TEST_CASE("Array construction reserves the entire outer payload before nested calls",
          "[negative][N4-02]") {
    ResourcePlan plan{2 * sizeof(NestedProbe) + sizeof(AllocationProbe)};
    const NestedProbeContext context{plan};
    const auto outer = plan.allocate<NestedProbe>(2, operation);
    REQUIRE(outer);
    CHECK(outer->values()[0].child.size() == 1);
    CHECK(outer->values()[1].child.size() == 0);
    CHECK(AllocationProbe::allocations == 1);
    CHECK(plan.allocated_bytes() == plan.byte_limit());
    CHECK(plan.remaining_bytes() == 0);
}

TEST_CASE("Outer construction failure rolls back only its own reservation", "[negative][N4-02]") {
    for (const auto outcome : {NestedOutcome::BadAlloc, NestedOutcome::OtherException}) {
        for (const bool child_fits : {false, true}) {
            const std::uint64_t child_bytes = child_fits ? sizeof(AllocationProbe) : 0;
            ResourcePlan plan{sizeof(NestedProbe) + child_bytes};
            const NestedProbeContext context{plan, outcome};
            const auto outer = plan.allocate<NestedProbe>(1, operation);
            REQUIRE_FALSE(outer);
            CHECK(outer.error().code == (outcome == NestedOutcome::BadAlloc
                                             ? FailureCode::ResourceExceeded
                                             : FailureCode::InvalidInput));
            CHECK(NestedProbe::child_failure.has_value() == !child_fits);
            // The child succeeded before the outer failure. Its cumulative charge must
            // survive even though unwinding destroyed the child buffer.
            CHECK(plan.allocated_bytes() == child_bytes);
            CHECK(plan.remaining_bytes() == sizeof(NestedProbe));
            const auto following = plan.allocate<unsigned char>(sizeof(NestedProbe), operation);
            REQUIRE(following);
            CHECK(plan.remaining_bytes() == 0);
        }
    }
}
