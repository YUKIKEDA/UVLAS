#include "aero/p1_induction.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#if defined(_WIN32)
#include <malloc.h>
#include <xmmintrin.h>
#endif

namespace {
// This process-local fault injector does not ship in the product. Disarm before
// unwinding the first allocation failure so reporting cannot fail recursively.
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
constinit void (*allocation_hook)() noexcept = nullptr;
void before_allocate() {
    const auto hook = std::exchange(allocation_hook, nullptr);
    if (hook != nullptr)
        hook();
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
TEST_CASE("P1 workspace exact planned payload and zero budget",
          "[UI-D06][N4-02][negative][verification]") {
    core::ResourcePlan empty{0};
    const auto rejected = P1Workspace::create(empty);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::ResourceExceeded);
    core::ResourcePlan exact{P1Workspace::payload_bytes()};
    auto result = P1Workspace::create(exact);
    REQUIRE(result);
    CHECK(exact.allocated_bytes() == P1Workspace::payload_bytes());
    core::ResourcePlan short_plan{P1Workspace::payload_bytes() - 1};
    CHECK_FALSE(P1Workspace::create(short_plan));
}
namespace {
LocalGeometryPacket resource_packet(core::ResourcePlan& plan, const TriangleFixture& fixture,
                                    QueryDefinition query) {
    auto geometry_workspace = workspace(plan);
    auto packet_workspace = LocalPacketWorkspace::create(plan);
    if (!packet_workspace)
        throw std::runtime_error("packet workspace");
    auto geometry =
        take(prepare_owned_geometry(fixture.input(), limits(), *geometry_workspace, plan));
    auto resolved = take(resolve_query(geometry, {QueryId{44}, query}, geometry.context(), limits(),
                                       *geometry_workspace, plan));
    return take(prepare_local_packet(geometry, resolved, FacetId{70},
                                     {core::Precision::TwoWord, 1e-10}, limits(),
                                     **packet_workspace, plan));
}
} // namespace
TEST_CASE("P1 allocation failures preserve inputs and recover without scratch allocation",
          "[UI-D06][N4-02][N4-03][N4-09][negative][verification]") {
    core::ResourcePlan setup{128ULL * 1024 * 1024};
    auto local = resource_packet(setup, TriangleFixture{}, FreePoint{FreeSingle{{.75, 1, 2}}});
    const auto identity = take(local.snapshot_identity());
    for (std::int64_t site = 0; site < 5; ++site) {
        core::ResourcePlan plan{P1Workspace::payload_bytes()};
        std::optional<core::Result<std::unique_ptr<P1Workspace>>> outcome;
        {
            const FaultScope fault{site};
            outcome.emplace(P1Workspace::create(plan));
        }
        REQUIRE(outcome.has_value());
        if (!outcome)
            return;
        if (!*outcome)
            CHECK(outcome->error().code == core::FailureCode::ResourceExceeded);
    }
    auto kernel = P1Workspace::create(setup);
    REQUIRE(kernel);
    core::ResourcePlan result_plan{P1Basis::payload_bytes()};
    std::optional<AeroResult<P1Basis>> failed_result;
    {
        const FaultScope fault{0};
        failed_result.emplace(
            evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                              {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)},
                              **kernel, result_plan));
    }
    REQUIRE(failed_result.has_value());
    if (!failed_result)
        return;
    failure(*failed_result, core::FailureCode::ResourceExceeded);
    CHECK(result_plan.allocated_bytes() == 0);
    // Exactly the owning result allocation is allowed. Any later scratch allocation fails.
    std::optional<AeroResult<P1Basis>> success;
    {
        const FaultScope fault{1};
        success.emplace(
            evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                              {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)},
                              **kernel, result_plan));
    }
    REQUIRE(success.has_value());
    if (!success)
        return;
    REQUIRE(success->has_value());
    CHECK(attempts.load() == 1);
    CHECK(result_plan.allocated_bytes() == P1Basis::payload_bytes());
    CHECK(take(local.snapshot_identity()) == identity);
}

namespace {
const LocalGeometryPacket* reentry_packet = nullptr;
P1Workspace* reentry_workspace = nullptr;
bool reentry_called = false, reentry_rejected = false;
void enter_again() noexcept {
    if (reentry_packet == nullptr || reentry_workspace == nullptr)
        return;
    reentry_called = true;
    core::ResourcePlan empty{0};
    const auto result =
        evaluate_p1_basis(*reentry_packet, {core::Precision::TwoWord, 1e-10, 1e-10},
                          {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)},
                          *reentry_workspace, empty);
    reentry_rejected = !result && result.error().cause.code == core::FailureCode::InvalidInput;
}
} // namespace
TEST_CASE("P1 allocator reentry is rejected before any nested action",
          "[UI-D06][N4-02][N4-03][N4-08][negative][verification]") {
    core::ResourcePlan plan{128ULL * 1024 * 1024};
    const auto local = resource_packet(plan, TriangleFixture{}, FreePoint{FreeSingle{{.75, 1, 2}}});
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    reentry_packet = &local;
    reentry_workspace = kernel->get();
    reentry_called = false;
    reentry_rejected = false;
    allocation_hook = &enter_again;
    auto result = evaluate_p1_basis(
        local, {core::Precision::TwoWord, 1e-10, 1e-10},
        {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)}, **kernel, plan);
    allocation_hook = nullptr;
    reentry_packet = nullptr;
    reentry_workspace = nullptr;
    REQUIRE(result.has_value());
    CHECK(reentry_called);
    CHECK(reentry_rejected);
}
TEST_CASE("P1 all evaluation families use only their owning result allocation",
          "[UI-D02][UI-D03][UI-D04][UI-D06][N4-02][N4-09][verification]") {
    core::ResourcePlan plan{256ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    for (unsigned which = 0; which < 6; ++which) {
        TriangleFixture fixture;
        QueryDefinition query = SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}};
        if (which == 0)
            query = FreePoint{FreeSingle{{.75, 1, 30}}};
        if (which == 1) {
            fixture.radius = SingleWord{0};
            query = FreePoint{FreeSingle{{6, 6, 0}}};
        }
        if (which == 2 || which == 3) {
            const int scale = which == 2 ? -900 : 900;
            for (auto& value : fixture.reference)
                value = std::scalbn(value, scale);
            fixture.current = fixture.reference;
            fixture.radius = SingleWord{std::scalbn(.25, scale)};
        }
        if (which == 4) {
            fixture.two_word = true;
            fixture.reference_low = fixture.reference;
            fixture.current_low = fixture.current;
            fixture.reference.fill(std::scalbn(1., 500));
            fixture.current = fixture.reference;
        }
        if (which == 5) {
            fixture.reference = {.25, -.125, .5, 1.25, .25, .625, .375, 1.125, .875};
            fixture.current = fixture.reference;
        }
        const auto local = resource_packet(plan, fixture, query);
        core::ResourcePlan output{P1Basis::payload_bytes()};
        std::optional<AeroResult<P1Basis>> evaluated;
        {
            const FaultScope fault{1};
            evaluated.emplace(evaluate_p1_basis(
                local, {core::Precision::TwoWord, 1e-10, 1e-10},
                {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)}, **kernel,
                output));
        }
        REQUIRE(evaluated.has_value());
        if (!evaluated)
            return;
        REQUIRE(evaluated->has_value());
        CHECK(attempts.load() == 1);
    }
}

TEST_CASE("P1 workspace allocation failures retain all supplied context and inner byte usage",
          "[UI-D06][N4-02][N4-03][N4-08][negative][verification]") {
    core::FailureDetails supplied;
    supplied.model = core::ModelGeneration{17};
    supplied.parent = core::AcceptedStateId{19};
    supplied.trial = core::TrialId{23};
    supplied.numerical_generation = core::NumericalGeneration{29};
    supplied.component = core::ComponentId{31};
    supplied.problem = core::ScopeId{37};
    supplied.last_accepted = core::AcceptedStateId{41};
    supplied.resources = core::ResourceUsage{core::ResourceKind::Nanoseconds, 7, 3};
    for (std::int64_t site = 0; site < 4; ++site) {
        core::ResourcePlan plan{P1Workspace::payload_bytes()};
        std::optional<core::Result<std::unique_ptr<P1Workspace>>> evaluated;
        {
            const FaultScope fault{site};
            evaluated.emplace(P1Workspace::create(plan, supplied));
        }
        REQUIRE(evaluated.has_value());
        if (!evaluated)
            return;
        REQUIRE_FALSE(evaluated->has_value());
        const auto& error = evaluated->error();
        CHECK(error.code == core::FailureCode::ResourceExceeded);
        CHECK(error.operation == core::Operation::EvaluateAero);
        CHECK(error.details.model == supplied.model);
        CHECK(error.details.parent == supplied.parent);
        CHECK(error.details.trial == supplied.trial);
        CHECK(error.details.numerical_generation == supplied.numerical_generation);
        CHECK(error.details.component == supplied.component);
        CHECK(error.details.problem == supplied.problem);
        CHECK(error.details.last_accepted == supplied.last_accepted);
        const auto& usage = error.details.resources;
        REQUIRE(usage.has_value());
        if (!usage)
            return;
        if (site < 3)
            CHECK(usage->kind == core::ResourceKind::Bytes);
        else
            CHECK(usage == supplied.resources);
    }
}
namespace {
struct EnvironmentRestore {
    int rounding = std::fegetround();
    ~EnvironmentRestore() {
        std::fesetround(rounding);
    }
};
void alter_rounding() noexcept {
    std::fesetround(FE_DOWNWARD);
}
} // namespace
TEST_CASE("P1 rejects unsupported rounding before and after allocation callbacks",
          "[N4-02][N4-03][negative][verification]") {
    core::ResourcePlan setup{128ULL * 1024 * 1024};
    auto local = resource_packet(setup, TriangleFixture{}, FreePoint{FreeSingle{{.75, 1, 2}}});
    auto kernel = P1Workspace::create(setup);
    REQUIRE(kernel);
    const EnvironmentRestore restore;
    core::ResourcePlan before{P1Basis::payload_bytes()};
    REQUIRE(std::fesetround(FE_DOWNWARD) == 0);
    auto rejected = evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                      {numerics::linear::EstimateAdmission{28672}, 1000000000},
                                      **kernel, before);
    std::fesetround(FE_TONEAREST);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().cause.code == core::FailureCode::PrecisionUnsupported);
    CHECK(before.remaining_bytes() == P1Basis::payload_bytes());
    core::ResourcePlan after{P1Basis::payload_bytes()};
    allocation_hook = &alter_rounding;
    rejected = evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                 {numerics::linear::EstimateAdmission{28672}, 1000000000}, **kernel,
                                 after);
    allocation_hook = nullptr;
    std::fesetround(FE_TONEAREST);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().cause.code == core::FailureCode::PrecisionUnsupported);
    CHECK(after.remaining_bytes() == 0);
    CHECK(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                            {numerics::linear::EstimateAdmission{28672}, 1000000000}, **kernel,
                            setup));
}
#if defined(_WIN32) && (defined(_M_X64) || defined(_M_IX86))
namespace {
struct MxcsrRestore {
    unsigned saved = _mm_getcsr();
    ~MxcsrRestore() {
        _mm_setcsr(saved);
    }
};
} // namespace
TEST_CASE("P1 rejects either subnormal flushing mode without reserving output",
          "[N4-02][N4-03][negative][verification]") {
    core::ResourcePlan setup{128ULL * 1024 * 1024};
    auto local = resource_packet(setup, TriangleFixture{}, FreePoint{FreeSingle{{.75, 1, 2}}});
    auto kernel = P1Workspace::create(setup);
    REQUIRE(kernel);
    const MxcsrRestore restore;
    for (const unsigned mask : {0x8000U, 0x0040U}) {
        core::ResourcePlan result{P1Basis::payload_bytes()};
        _mm_setcsr((restore.saved & ~0x8040U) | mask);
        auto rejected = evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                          {numerics::linear::EstimateAdmission{28672}, 1000000000},
                                          **kernel, result);
        _mm_setcsr(restore.saved);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().cause.code == core::FailureCode::PrecisionUnsupported);
        CHECK(result.remaining_bytes() == P1Basis::payload_bytes());
    }
}
#endif
