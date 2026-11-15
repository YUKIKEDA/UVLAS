#include "../support/allocation_audit.hpp"
#include "backend_memory_detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <limits>
#include <type_traits>
#include <uvlas/numerics/linear_solve.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics::linear;
using detail::BackendPlanFactory;
} // namespace

TEST_CASE("Backend forecast admission does not manufacture proof or observed bytes",
          "[unit][negative][N4-09]") {
    const auto plan = BackendPlanFactory::make(10);
    const auto exact = admit_backend_memory(plan, EstimateAdmission{10});
    REQUIRE(exact);
    CHECK(exact->admitted_payload_bytes() == 10);
    CHECK_FALSE(exact->plan().proven_peak_payload_upper_bound_bytes());
    const auto exceeded = admit_backend_memory(plan, EstimateAdmission{9});
    REQUIRE_FALSE(exceeded);
    CHECK(exceeded.error().code == core::FailureCode::ResourceExceeded);
    CHECK_FALSE(exceeded.error().details.resources);
    for (const auto limit :
         std::array<std::uint64_t, 3>{0, 10, std::numeric_limits<std::uint64_t>::max()}) {
        const auto strict = admit_backend_memory(plan, RequireProvenBound{limit});
        REQUIRE_FALSE(strict);
        CHECK(strict.error().code == core::FailureCode::Unsupported);
        CHECK_FALSE(strict.error().details.resources);
    }
}
TEST_CASE("Backend strict admission uses its proven amount including zero and uint64 endpoint",
          "[unit][N4-09]") {
    const auto plan = BackendPlanFactory::make(10, 20);
    const auto strict = admit_backend_memory(plan, RequireProvenBound{20});
    REQUIRE(strict);
    CHECK(strict->admitted_payload_bytes() == 20);
    CHECK_FALSE(admit_backend_memory(plan, RequireProvenBound{19}));
    REQUIRE(admit_backend_memory(plan, EstimateAdmission{10}));
    const auto empty = BackendPlanFactory::make(0, 0);
    REQUIRE(admit_backend_memory(empty, RequireProvenBound{0}));
    REQUIRE(admit_backend_memory(empty, EstimateAdmission{0}));
    const auto largest = std::numeric_limits<std::uint64_t>::max();
    const auto endpoint = BackendPlanFactory::make(largest, largest);
    REQUIRE(admit_backend_memory(endpoint, RequireProvenBound{largest}));
    CHECK_FALSE(admit_backend_memory(endpoint, EstimateAdmission{largest - 1}));
}
TEST_CASE("Backend admission is allocation free and leaves owned array accounting separate",
          "[property][N4-09]") {
    const auto plan = BackendPlanFactory::make(12);
    test_support::AllocationAudit audit{test_support::cpp_backend,
                                        test_support::AuditThreadDomain::AllThreads};
    const auto pass = admit_backend_memory(plan, EstimateAdmission{12});
    const auto fail = admit_backend_memory(plan, RequireProvenBound{12});
    const auto result = audit.finish();
    REQUIRE(pass);
    REQUIRE_FALSE(fail);
    CHECK(result.status == test_support::AllocationAuditStatus::Clean);
    CHECK(result.cpp_allocation_attempts == 0);
}
static_assert(!std::is_default_constructible_v<BackendMemoryPlan>);
static_assert(
    !std::is_constructible_v<BackendMemoryPlan, std::uint64_t, std::optional<std::uint64_t>>);

static_assert(!std::is_default_constructible_v<BackendMemoryLimit>);

static_assert(!std::is_default_constructible_v<Limits>);
