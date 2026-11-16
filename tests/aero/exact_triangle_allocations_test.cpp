#include "../support/allocation_audit.hpp"
#include "aero/exact_triangle_predicates.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <chrono>
#include <iostream>
#include <optional>
#include <uvlas/core/instrumentation.hpp>

namespace {
namespace subject = uvlas::aero::detail;
using namespace uvlas;
std::unique_ptr<subject::PredicateWorkspace> prepared(core::ResourcePlan& plan) {
    auto result = subject::PredicateWorkspace::create(plan);
    if (!result) {
        FAIL_CHECK("Workspace preparation failed");
        return nullptr;
    }
    return std::move(*result);
}
} // namespace
TEST_CASE("Exact predicates report actual preparation allocations separately from forecasts",
          "[unit][N4-09][UI-D06]") {
    core::ResourcePlan preparation{4ULL * 1024 * 1024}, output{4ULL * 1024 * 1024};
    auto workspace = prepared(preparation);
    REQUIRE(workspace != nullptr);
    const std::array<double, 12> coordinates{0, 0, 0, 3, 0, 0, 0, 4, 0, .75, 1, 20};
    std::optional<core::Result<subject::ExactTriangleResult>> actual;
    const auto start = std::chrono::steady_clock::now();
    test_support::AllocationAudit audit{test_support::cpp_backend,
                                        test_support::AuditThreadDomain::CallingThread};
    actual.emplace(subject::classify_triangle(subject::Binary64Coordinates{coordinates}, *workspace,
                                              output,
                                              {numerics::linear::EstimateAdmission{28672}, 100}));
    const auto observed = audit.finish();
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                             std::chrono::steady_clock::now() - start)
                             .count();
    if (!actual.has_value()) {
        FAIL_CHECK("No response captured");
        return;
    }
    const auto& response = *actual;
    REQUIRE(response.has_value());
    CHECK(observed.status == test_support::AllocationAuditStatus::AllocationObserved);
    CHECK(observed.cpp_allocation_attempts >=
          1); // At least the independently owned result payload.
    CHECK(output.allocated_bytes() == subject::ExactTriangleResult::payload_bytes());
    CHECK(preparation.allocated_bytes() == subject::PredicateWorkspace::payload_bytes());
    CHECK_FALSE(response->memory_plan().proven_peak_payload_upper_bound_bytes());
    const auto memory = core::sample_process_memory(core::Operation::EvaluateAero);
    std::cout << "aero exact preparation allocations=" << observed.cpp_allocation_attempts
              << " elapsed_ns=" << elapsed
              << " planned_workspace_bytes=" << preparation.allocated_bytes()
              << " planned_result_bytes=" << output.allocated_bytes();
    if (memory)
        std::cout << " process_peak_bytes=" << memory->peak_resident_bytes;
    std::cout << '\n';
}
TEST_CASE("An unsupported strict predicate request performs no hidden backend allocation",
          "[negative][N4-09][UI-D06]") {
    core::ResourcePlan preparation{4ULL * 1024 * 1024}, output{4ULL * 1024 * 1024};
    auto workspace = prepared(preparation);
    REQUIRE(workspace != nullptr);
    const std::array<double, 12> coordinates{0, 0, 0, 3, 0, 0, 0, 4, 0, .75, 1, 20};
    std::optional<core::Result<subject::ExactTriangleResult>> actual;
    test_support::AllocationAudit audit{test_support::cpp_backend,
                                        test_support::AuditThreadDomain::CallingThread};
    actual.emplace(subject::classify_triangle(subject::Binary64Coordinates{coordinates}, *workspace,
                                              output,
                                              {numerics::linear::RequireProvenBound{28672}, 100}));
    const auto observed = audit.finish();
    if (!actual.has_value()) {
        FAIL_CHECK("No response captured");
        return;
    }
    const auto& response = *actual;
    REQUIRE_FALSE(response.has_value());
    CHECK(response.error().code == core::FailureCode::Unsupported);
    CHECK(observed.status == test_support::AllocationAuditStatus::Clean);
    CHECK(observed.cpp_allocation_attempts == 0);
    CHECK(output.allocated_bytes() == 0);
}
