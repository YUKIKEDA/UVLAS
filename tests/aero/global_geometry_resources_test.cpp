#include "aero/global_geometry_private.hpp"
#include "geometry_snapshot_fixture.hpp"
#include "global_kinematics_fixture.hpp"

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
using namespace uvlas;
namespace s = aero::detail;
constexpr std::uint64_t memory_budget = 64ULL * 1024 * 1024;
s::GlobalLimits global_limits() {
    return {UINT64_MAX / 4, numerics::linear::EstimateAdmission{4096}};
}
template <class T>
T take_global(s::GlobalResult<T> value) {
    if (!value) {
        FAIL_CHECK("global failure " << static_cast<unsigned>(value.error().origin.cause.code));
        throw std::runtime_error("global failure");
    }
    return std::move(*value);
}
struct Setup {
    test_geometry::TriangleFixture fixture;
    core::ResourcePlan memory{memory_budget};
    std::unique_ptr<s::GeometryWorkspace> local = test_geometry::workspace(memory);
    std::unique_ptr<s::GlobalWorkspace> workspace;
    s::OwnedGeometry geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.input(), test_geometry::limits(), *local, memory));
    std::array<s::FaceRoleInput, 1> roles{
        s::FaceRoleInput{fixture.faces[0].id, s::AeroFaceRole::Bound}};
    std::array<s::EdgeBindingInput, 3> edges{
        s::EdgeBindingInput{fixture.edges[0].id, s::CanonicalEdge{}},
        s::EdgeBindingInput{fixture.edges[1].id, s::CanonicalEdge{}},
        s::EdgeBindingInput{fixture.edges[2].id, s::CanonicalEdge{}}};
    Setup() {
        auto result = s::GlobalWorkspace::create(memory, numerics::linear::EstimateAdmission{4096});
        REQUIRE(result);
        workspace = std::move(*result);
    }
    [[nodiscard]] s::GlobalTopologyInput input() const {
        return {roles, edges, {}};
    }
};
struct Reentry {
    Setup* setup{};
    core::ResourcePlan* plan{};
    bool rejected = false;
    bool unchanged = false;
};
constinit Reentry* reentry = nullptr;
constinit std::stop_source* final_stop = nullptr;
constinit unsigned final_allocations = 0;
void cancel_final_allocation() {
    if (--final_allocations == 0)
        final_stop->request_stop();
    else
        allocation_callback = cancel_final_allocation;
}
void attempt_reentry() {
    const auto before = reentry->plan->allocated_bytes();
    const auto result =
        s::build_global_topology(reentry->setup->geometry, reentry->setup->input(), global_limits(),
                                 *reentry->setup->workspace, *reentry->plan);
    reentry->rejected =
        !result && result.error().origin.cause.code == core::FailureCode::InvalidInput;
    reentry->unchanged = before == reentry->plan->allocated_bytes();
}
} // namespace
TEST_CASE("Global workspace blocks allocation-time reentry before charging the nested call",
          "[N4-02][N4-03][N4-09]") {
    Setup setup;
    core::ResourcePlan plan{memory_budget};
    Reentry state{&setup, &plan};
    reentry = &state;
    allocation_callback = attempt_reentry;
    const auto result = s::build_global_topology(setup.geometry, setup.input(), global_limits(),
                                                 *setup.workspace, plan);
    allocation_callback = nullptr;
    reentry = nullptr;
    REQUIRE(result);
    CHECK(state.rejected);
    CHECK(state.unchanged);
}
TEST_CASE("Cancellation from the final metadata allocator cannot publish a topology",
          "[UW-D07][N4-03][N4-09]") {
    Setup setup;
    core::ResourcePlan plan{memory_budget};
    std::stop_source stop;
    final_stop = &stop;
    final_allocations = 6;
    allocation_callback = cancel_final_allocation;
    const auto result = s::build_global_topology(setup.geometry, setup.input(), global_limits(),
                                                 *setup.workspace, plan, stop.get_token());
    allocation_callback = nullptr;
    final_stop = nullptr;
    CHECK(stop.stop_requested());
    REQUIRE_FALSE(result);
    CHECK(result.error().origin.cause.code == core::FailureCode::Cancelled);
    core::ResourcePlan recovery{memory_budget};
    CHECK(s::build_global_topology(setup.geometry, setup.input(), global_limits(), *setup.workspace,
                                   recovery));
}
TEST_CASE("Global workspace allocation faults preserve all supplied diagnostics",
          "[N4-02][N4-03][N4-09]") {
    auto details = test_geometry::context().diagnostics;
    details.component = core::ComponentId{71};
    details.problem = core::DofId{73};
    details.model = core::ModelGeneration{75};
    details.numerical_generation = core::NumericalGeneration{76};
    for (std::int64_t fail = 0; fail < 2; ++fail) {
        core::ResourcePlan memory{memory_budget};
        const auto result = [&] {
            FaultScope const fault{fail};
            return s::GlobalWorkspace::create(memory, numerics::linear::EstimateAdmission{4096},
                                              details);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        CHECK(result.error().details.model == details.model);
        CHECK(result.error().details.parent == details.parent);
        CHECK(result.error().details.trial == details.trial);
        CHECK(result.error().details.numerical_generation == details.numerical_generation);
        CHECK(result.error().details.component == details.component);
        CHECK(result.error().details.problem == details.problem);
        CHECK(result.error().details.last_accepted == details.last_accepted);
    }
    core::ResourcePlan equal{s::GlobalWorkspace::payload_bytes()};
    CHECK(s::GlobalWorkspace::create(equal, numerics::linear::EstimateAdmission{4096}));
    core::ResourcePlan short_plan{s::GlobalWorkspace::payload_bytes() - 1};
    auto short_result =
        s::GlobalWorkspace::create(short_plan, numerics::linear::EstimateAdmission{4096}, details);
    REQUIRE_FALSE(short_result);
    CHECK(short_plan.allocated_bytes() == 0);
    CHECK(short_result.error().details.parent == details.parent);
    CHECK(short_result.error().details.trial == details.trial);
    core::ResourcePlan unused{memory_budget};
    CHECK_FALSE(
        s::GlobalWorkspace::create(unused, numerics::linear::RequireProvenBound{UINT64_MAX}));
    CHECK(unused.allocated_bytes() == 0);
}
TEST_CASE(
    "Global topology and contact calls admit whole payload and recover from every allocation fault",
    "[UW-D07][N4-02][N4-03][N4-09]") {
    Setup setup;
    core::ResourcePlan measured{memory_budget};
    auto topology = [&] {
        FaultScope const track{std::numeric_limits<std::int64_t>::max()};
        return s::build_global_topology(setup.geometry, setup.input(), global_limits(),
                                        *setup.workspace, measured);
    }();
    REQUIRE(topology);
    const auto allocations = attempts.load();
    CHECK(allocations == 6); // Five nonempty planned arrays + shared control owner.
    const auto footprint = topology->work_plan()->payload_bytes;
    CHECK(measured.allocated_bytes() == footprint);
    for (std::uint64_t point = 0; point < allocations; ++point) {
        core::ResourcePlan memory{memory_budget};
        const auto result = [&] {
            FaultScope const fault{static_cast<std::int64_t>(point)};
            return s::build_global_topology(setup.geometry, setup.input(), global_limits(),
                                            *setup.workspace, memory);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().origin.cause.code == core::FailureCode::ResourceExceeded);
        CHECK(result.error().origin.cause.details.parent ==
              setup.geometry.context().diagnostics.parent);
        core::ResourcePlan recovery{memory_budget};
        CHECK(s::build_global_topology(setup.geometry, setup.input(), global_limits(),
                                       *setup.workspace, recovery));
    }
    core::ResourcePlan equal{footprint};
    CHECK(s::build_global_topology(setup.geometry, setup.input(), global_limits(), *setup.workspace,
                                   equal));
    core::ResourcePlan one_short{footprint - 1};
    CHECK_FALSE(s::build_global_topology(setup.geometry, setup.input(), global_limits(),
                                         *setup.workspace, one_short));
    CHECK(one_short.allocated_bytes() == 0);
    auto metadata = global_limits();
    metadata.metadata =
        numerics::linear::EstimateAdmission{topology->work_plan()->metadata_forecast};
    core::ResourcePlan metadata_equal{memory_budget};
    CHECK(s::build_global_topology(setup.geometry, setup.input(), metadata, *setup.workspace,
                                   metadata_equal));
    metadata.metadata =
        numerics::linear::EstimateAdmission{topology->work_plan()->metadata_forecast - 1};
    core::ResourcePlan metadata_short{memory_budget};
    CHECK_FALSE(s::build_global_topology(setup.geometry, setup.input(), metadata, *setup.workspace,
                                         metadata_short));
    CHECK(metadata_short.allocated_bytes() == 0);
    auto wl = global_limits();
    wl.max_steps = topology->work_plan()->reserved_steps;
    core::ResourcePlan exact_work{memory_budget};
    CHECK(
        s::build_global_topology(setup.geometry, setup.input(), wl, *setup.workspace, exact_work));
    --wl.max_steps;
    core::ResourcePlan short_work{memory_budget};
    CHECK_FALSE(
        s::build_global_topology(setup.geometry, setup.input(), wl, *setup.workspace, short_work));
    CHECK(short_work.allocated_bytes() == 0);
    core::ResourcePlan contact_memory{memory_budget};
    const auto contact = [&] {
        FaultScope const track{std::numeric_limits<std::int64_t>::max()};
        return s::check_global_geometry(*topology, setup.geometry, global_limits(),
                                        *setup.workspace, contact_memory);
    }();
    REQUIRE(contact);
    const auto contact_allocations = attempts.load();
    CHECK(contact_allocations == 3);
    for (std::uint64_t point = 0; point < contact_allocations; ++point) {
        core::ResourcePlan memory{memory_budget};
        const auto result = [&] {
            FaultScope const fault{static_cast<std::int64_t>(point)};
            return s::check_global_geometry(*topology, setup.geometry, global_limits(),
                                            *setup.workspace, memory);
        }();
        REQUIRE_FALSE(result);
        CHECK(result.error().origin.cause.code == core::FailureCode::ResourceExceeded);
    }
    const auto contact_plan = *contact->work_plan();
    core::ResourcePlan contact_equal{contact_plan.payload_bytes};
    CHECK(s::check_global_geometry(*topology, setup.geometry, global_limits(), *setup.workspace,
                                   contact_equal));
    core::ResourcePlan contact_short{contact_plan.payload_bytes - 1};
    CHECK_FALSE(s::check_global_geometry(*topology, setup.geometry, global_limits(),
                                         *setup.workspace, contact_short));
    CHECK(contact_short.allocated_bytes() == 0);
    auto cl = global_limits();
    cl.max_steps = contact_plan.reserved_steps;
    core::ResourcePlan exact_contact{memory_budget};
    CHECK(s::check_global_geometry(*topology, setup.geometry, cl, *setup.workspace, exact_contact));
    --cl.max_steps;
    core::ResourcePlan short_contact{memory_budget};
    CHECK_FALSE(
        s::check_global_geometry(*topology, setup.geometry, cl, *setup.workspace, short_contact));
}
TEST_CASE("Fixed dense integer products have no hidden heap arena", "[UW-D07][N4-02][N4-06]") {
    core::ResourcePlan memory{memory_budget};
    auto workspace = s::GlobalWorkspace::create(memory, numerics::linear::EstimateAdmission{4096});
    REQUIRE(workspace);
    auto prepared = memory.allocate<s::GlobalScratch>(1, core::Operation::EvaluateAero);
    REQUIRE(prepared);
    auto scratch = std::move(*prepared); // Prepared test-only primitive arena.
    for (const auto widths : std::array<std::pair<unsigned, unsigned>, 4>{
             {{3200, 3200}, {6400, 6400}, {19000, 9000}, {16000, 16000}}}) {
        s::ExactInteger left = (s::ExactInteger{1} << widths.first) - 1,
                        right = (s::ExactInteger{1} << widths.second) - 1;
        const boost::multiprecision::cpp_int expected =
            boost::multiprecision::cpp_int{left} * boost::multiprecision::cpp_int{right};
        s::global_detail::Run run{
            global_limits(), test_geometry::context(), scratch.values()[0], {}};
        run.admit(s::global_detail::multiply_cost);
        auto actual = [&] {
            FaultScope const no_heap{0};
            return s::global_detail::multiply(left, right, run);
        }();
        CHECK(attempts.load() == 0);
        CHECK(boost::multiprecision::cpp_int{actual} == expected);
        run.reserved += s::global_detail::multiply_cost;
        CHECK(s::global_detail::multiply(-left, right, run) == -actual);
    }
    s::global_detail::Run run{global_limits(), test_geometry::context(), scratch.values()[0], {}};
    run.admit(s::global_detail::multiply_cost);
    const s::ExactInteger oversized = (s::ExactInteger{1} << 19000U) - 1;
    CHECK_THROWS_AS(s::global_detail::multiply(oversized, oversized, run), s::global_detail::Abort);
    CHECK_THROWS_AS(s::global_detail::plus(UINT64_MAX, 1), s::global_detail::Abort);
    CHECK_THROWS_AS(s::global_detail::times(UINT64_MAX, 2), s::global_detail::Abort);
    CHECK_THROWS_AS(s::global_detail::tree_shape(UINT64_MAX), s::global_detail::Abort);
}
TEST_CASE("Alias CSR payload and every allocation fault preserve the source and qualification",
          "[UW-D07][N4-02][N4-03][N4-09]") {
    test_global_geometry::KinematicsFixture fixture;
    core::ResourcePlan initial{memory_budget};
    auto old = test_geometry::workspace(initial);
    auto geometry = test_geometry::take(
        s::prepare_owned_geometry(fixture.geometry(), test_geometry::limits(), *old, initial));
    auto made = s::GlobalWorkspace::create(initial, numerics::linear::EstimateAdmission{4096});
    REQUIRE(made);
    auto w = std::move(*made);
    core::ResourcePlan measured{memory_budget};
    auto qualified = [&] {
        FaultScope const track{std::numeric_limits<std::int64_t>::max()};
        return s::check_alias_kinematics(geometry, fixture.input(), global_limits(), *w, measured);
    }();
    REQUIRE(qualified);
    const auto count = attempts.load();
    CHECK(count == 5);
    const auto original = *qualified->velocity(s::PointId{0}, 0);
    const auto original_word = original.limbs[0];
    for (std::uint64_t fault = 0; fault < count; ++fault) {
        core::ResourcePlan plan{memory_budget};
        const auto rejected = [&] {
            FaultScope const injection{static_cast<std::int64_t>(fault)};
            return s::check_alias_kinematics(geometry, fixture.input(), global_limits(), *w, plan);
        }();
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().origin.cause.code == core::FailureCode::ResourceExceeded);
        CHECK(qualified->velocity(s::PointId{0}, 0)->limbs[0] == original_word);
    }
    const auto data = *qualified->work_plan();
    core::ResourcePlan equal{data.payload_bytes};
    CHECK(s::check_alias_kinematics(geometry, fixture.input(), global_limits(), *w, equal));
    core::ResourcePlan short_plan{data.payload_bytes - 1};
    CHECK_FALSE(
        s::check_alias_kinematics(geometry, fixture.input(), global_limits(), *w, short_plan));
    CHECK(short_plan.allocated_bytes() == 0);
    auto work = global_limits();
    work.max_steps = data.reserved_steps;
    core::ResourcePlan enough_work{memory_budget};
    CHECK(s::check_alias_kinematics(geometry, fixture.input(), work, *w, enough_work));
    --work.max_steps;
    core::ResourcePlan not_enough{memory_budget};
    CHECK_FALSE(s::check_alias_kinematics(geometry, fixture.input(), work, *w, not_enough));
    CHECK(not_enough.allocated_bytes() == 0);
}
