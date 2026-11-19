// Allocation replacement follows the existing typed-schema resource executable.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <stop_token>
#include <string_view>
#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
// Process-local injection belongs only to this dedicated test executable.
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
// NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop mutates the stop source.
thread_local std::stop_source* cancel_at_allocation = nullptr;
thread_local std::int64_t cancel_after = -1;
// Replacement new injects only its standard allocation failure. Throwing
// arbitrary exceptions here can violate standard-library noexcept contracts;
// core/failure_test.cpp tests guard_entry's arbitrary-exception guarantee.
void allocation_attempt() {
    if (cancel_at_allocation && cancel_after-- == 0)
        cancel_at_allocation->request_stop();
    const auto next = remaining.load(std::memory_order_relaxed);
    if (next != -1) {
        attempts.fetch_add(1, std::memory_order_relaxed);
        if (next == 0) {
            // Disable before exception construction and stack unwinding.
            remaining.store(-1, std::memory_order_relaxed);
            throw std::bad_alloc{};
        }
        if (next > 0)
            remaining.fetch_sub(1, std::memory_order_relaxed);
    }
}
void* allocate(std::size_t count) {
    allocation_attempt();
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable pointee storage.
    void* const value = std::malloc(count == 0 ? 1 : count);
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
void* allocate_aligned(std::size_t count, std::size_t alignment) {
    allocation_attempt();
    count = count == 0 ? 1 : count;
#if defined(_WIN32)
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* const value = _aligned_malloc(count, alignment);
#else
    if (count > std::numeric_limits<std::size_t>::max() - (alignment - 1))
        throw std::bad_alloc{};
    const auto rounded = ((count + alignment - 1) / alignment) * alignment;
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* const value = std::aligned_alloc(alignment, rounded);
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
void* operator new(std::size_t count) {
    return allocate(count);
}
void* operator new[](std::size_t count) {
    return allocate(count);
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
// Cover the forms used by library temporary buffers as well as user new
// expressions. Otherwise ASan's default nothrow new can pair with our free.
void* operator new(std::size_t count, const std::nothrow_t&) noexcept {
    try {
        return allocate(count);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void* operator new[](std::size_t count, const std::nothrow_t&) noexcept {
    try {
        return allocate(count);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void operator delete(void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void operator delete[](void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void* operator new(std::size_t count, std::align_val_t alignment) {
    return allocate_aligned(count, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t count, std::align_val_t alignment) {
    return allocate_aligned(count, static_cast<std::size_t>(alignment));
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
void* operator new(std::size_t count, std::align_val_t alignment, const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(count, static_cast<std::size_t>(alignment));
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void* operator new[](std::size_t count, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(count, static_cast<std::size_t>(alignment));
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void operator delete(void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}

#include "case_schema_snapshot.hpp"
#include "reference_declarations_fixture.hpp"

#include <future>
#include <stop_token>
namespace {
using namespace uvlas;
using namespace test::declaration_fixture;
std::uint64_t nodes(const JsonToken& v) {
    std::uint64_t count = 1;
    for (const auto& m : v.members)
        count += nodes(m.second);
    for (const auto& child : v.elements)
        count += nodes(child);
    return count;
}
std::uint64_t depth(const JsonToken& v, std::uint64_t edge = 0) {
    auto maximum = edge + (v.kind == JsonToken::Kind::Scalar ? 0 : 1);
    for (const auto& m : v.members)
        maximum = std::max(maximum, depth(m.second, edge + 1));
    for (const auto& child : v.elements)
        maximum = std::max(maximum, depth(child, edge + 1));
    return maximum;
}
} // namespace
TEST_CASE("Declaration caps are checked before result allocation and preserve earlier ownership",
          "[IO4-01][negative][property]") {
    const auto fixture = all_roles();
    const auto source = checked(fixture);
    core::ResourcePlan oldplan{budget};
    const auto old = take(rd::collect_case_file_references(source, limits(), oldplan, {}, nullptr));
    const auto before = snapshot(old);
    const auto parent = test::schema_snapshot::serialize(take(source.root()));
    const auto report = take(old.memory_report());
    std::uint64_t path_bytes = 0, text_bytes = 0;
    for (core::Index i = 0; i < 8; ++i) {
        const auto row = take(old.entry(i));
        path_bytes += row.original_path.size();
        text_bytes += row.original_path.size() + row.source_id_spelling.size() +
                      row.digest_spelling.size() + row.declaration_pointer.size();
    }
    const std::array boundaries{
        std::pair{rd::DeclarationQuota::SourceNodes, nodes(fixture)},
        std::pair{rd::DeclarationQuota::SourceDepth, depth(fixture)},
        std::pair{rd::DeclarationQuota::Declarations, std::uint64_t{8}},
        std::pair{rd::DeclarationQuota::SourcePathBytes, path_bytes},
        std::pair{rd::DeclarationQuota::OwnedTextBytes, text_bytes},
        std::pair{rd::DeclarationQuota::OwnedBytes, report.planned_owned_bytes},
        std::pair{rd::DeclarationQuota::ScratchBytes, report.planned_scratch_bytes},
        std::pair{rd::DeclarationQuota::CumulativeBytes, report.admission_bytes},
        std::pair{rd::DeclarationQuota::Work, report.planned_work}};
    for (const auto& [q, bound] : boundaries) {
        INFO(static_cast<unsigned>(q));
        REQUIRE(bound > 0);
        core::ResourcePlan exact{budget};
        CHECK(rd::collect_case_file_references(source, with(q, bound), exact, {}, nullptr)
                  .has_value());
        core::ResourcePlan refused{budget};
        std::optional<core::Result<rd::DeclaredReferences>> result;
        {
            const FaultScope audit{-2};
            result.emplace(
                rd::collect_case_file_references(source, with(q, bound - 1), refused, {}, nullptr));
        }
        REQUIRE_FALSE(result->has_value());
        CHECK(result->error().code == core::FailureCode::ResourceExceeded);
        CHECK(attempts.load() == 0);
        CHECK(refused.allocated_bytes() == 0);
        CHECK(snapshot(old) == before);
        CHECK(test::schema_snapshot::serialize(take(source.root())) == parent);
    }
    core::ResourcePlan exact{report.admission_bytes};
    CHECK(rd::collect_case_file_references(source, limits(), exact, {}, nullptr).has_value());
    core::ResourcePlan small{report.admission_bytes - 1};
    CHECK_FALSE(rd::collect_case_file_references(source, limits(), small, {}, nullptr));
    CHECK(small.allocated_bytes() == 0);
    core::ResourcePlan no_path{budget};
    CHECK(rd::collect_case_file_references(source, with(rd::DeclarationQuota::DiagnosticBytes, 0),
                                           no_path, {}, nullptr)
              .has_value());
}
TEST_CASE("Every declaration allocation failure preserves complete payload and successful charges",
          "[IO4-01][negative][property]") {
    const auto source = checked(all_roles());
    core::ResourcePlan oldplan{budget};
    const auto old = take(rd::collect_case_file_references(source, limits(), oldplan, {}, nullptr));
    const auto before = snapshot(old);
    const auto parent = test::schema_snapshot::serialize(take(source.root()));
    core::ResourcePlan measured{budget};
    std::optional<core::Result<rd::DeclaredReferences>> control;
    {
        const FaultScope audit{-2};
        control.emplace(rd::collect_case_file_references(source, limits(), measured, {}, nullptr));
    }
    REQUIRE(control->has_value());
    const auto count = attempts.load();
    REQUIRE(count == 3);
    for (std::uint64_t i = 0; i < count; ++i) {
        core::ResourcePlan plan{budget};
        rd::DeclarationDiagnostic slot;
        std::optional<core::Result<rd::DeclaredReferences>> output;
        {
            const FaultScope fault{static_cast<std::int64_t>(i)};
            output.emplace(rd::collect_case_file_references(source, limits(), plan, {}, &slot));
        }
        REQUIRE_FALSE(output->has_value());
        CHECK(output->error().code == core::FailureCode::ResourceExceeded);
        const auto detail = slot.detail();
        if (!detail) {
            FAIL("Expected owned value");
            return;
        }
        CHECK(detail->reason == rd::DeclarationReason::Allocation);
        CHECK_FALSE(detail->node.has_value());
        CHECK_FALSE(slot.path().has_value());
        const auto charged = plan.allocated_bytes();
        CHECK((i == 0 ? charged == 0 : charged > 0));
        slot.reset();
        CHECK(plan.allocated_bytes() == charged);
        CHECK(snapshot(old) == before);
        CHECK(test::schema_snapshot::serialize(take(source.root())) == parent);
    }
}
TEST_CASE("Declaration diagnostics keep primary failure and cancellation never invents a node",
          "[IO4-01][negative]") {
    const auto source = checked(all_roles());
    core::ResourcePlan plan{budget};
    rd::DeclarationDiagnostic slot;
    auto failed = rd::collect_case_file_references(
        source, with(rd::DeclarationQuota::Declarations, 0), plan, {}, &slot);
    REQUIRE_FALSE(failed);
    const auto first = slot.detail();
    if (!first) {
        FAIL("Expected owned value");
        return;
    }
    CHECK(first->node.has_value());
    const auto path = slot.path();
    if (!path) {
        FAIL("Expected complete diagnostic path");
        return;
    }
    CHECK(*path == "/model/geometry/0/mesh");
    const auto charged = plan.allocated_bytes();
    CHECK(charged == path->size());
    for (const auto q :
         {rd::DeclarationQuota::DiagnosticBytes, rd::DeclarationQuota::CumulativeBytes}) {
        std::array<std::uint64_t, static_cast<std::size_t>(rd::DeclarationQuota::Count)> caps{};
        caps.fill(budget);
        caps[static_cast<std::size_t>(rd::DeclarationQuota::Declarations)] = 0;
        caps[static_cast<std::size_t>(q)] = 0;
        core::ResourcePlan limited{budget};
        auto negative = rd::collect_case_file_references(source, rd::DeclarationLimits{caps},
                                                         limited, {}, &slot);
        REQUIRE_FALSE(negative);
        CHECK(negative.error().code == failed.error().code);
        const auto detail = slot.detail();
        if (!detail) {
            FAIL("Expected owned value");
            return;
        }
        CHECK(detail->reason == first->reason);
        CHECK(detail->node == first->node);
        CHECK_FALSE(slot.path());
        CHECK(limited.allocated_bytes() == 0);
    }
    core::ResourcePlan faultplan{budget};
    std::optional<core::Result<rd::DeclaredReferences>> fault;
    {
        const FaultScope inject{0};
        fault.emplace(rd::collect_case_file_references(
            source, with(rd::DeclarationQuota::Declarations, 0), faultplan, {}, &slot));
    }
    REQUIRE_FALSE(fault->has_value());
    CHECK(fault->error().code == failed.error().code);
    CHECK_FALSE(slot.path());
    CHECK(faultplan.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop is non-const; libstdc++ differs.
    std::stop_source stop;
    stop.request_stop();
    core::ResourcePlan cancelplan{budget};
    auto cancelled =
        rd::collect_case_file_references(source, limits(), cancelplan, stop.get_token(), &slot);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    const auto detail = slot.detail();
    if (!detail) {
        FAIL("Expected owned value");
        return;
    }
    CHECK_FALSE(detail->parent);
    CHECK_FALSE(detail->node);
    CHECK_FALSE(slot.path());
    CHECK(cancelplan.allocated_bytes() == 0);
}
TEST_CASE("Independent declaration calls can use the same immutable case concurrently",
          "[IO4-01][property]") {
    const auto source = checked(all_roles());
    auto worker = [&source] {
        core::ResourcePlan plan{budget};
        rd::DeclarationDiagnostic slot;
        return snapshot(take(rd::collect_case_file_references(source, limits(), plan, {}, &slot)));
    };
    auto first = std::async(std::launch::async, worker);
    auto second = std::async(std::launch::async, worker);
    CHECK(first.get() == second.get());
}

TEST_CASE("Cancellation during the final declaration allocation cannot publish a partial result",
          "[IO4-01][negative]") {
    const auto source = checked(all_roles());
    core::ResourcePlan plan{budget};
    rd::DeclarationDiagnostic slot;
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop is non-const.
    std::stop_source stop;
    cancel_at_allocation = &stop;
    cancel_after = 2;
    auto result = rd::collect_case_file_references(source, limits(), plan, stop.get_token(), &slot);
    cancel_at_allocation = nullptr;
    cancel_after = -1;
    REQUIRE_FALSE(result);
    CHECK(result.error().code == core::FailureCode::Cancelled);
    CHECK(plan.allocated_bytes() > 0);
    const auto charged = plan.allocated_bytes();
    slot.reset();
    CHECK(plan.allocated_bytes() == charged);
}
