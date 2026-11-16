#include "../support/allocation_audit.hpp"
#include "eigen_full_phase_private.hpp"

#include <catch2/catch_test_macros.hpp>

#include <limits>
TEST_CASE("Full exact phase comparisons allocate no heap even at extreme exponents",
          "[property][N4-07][N4-09]") {
    namespace d = uvlas::numerics::eigen::detail;
    namespace t = uvlas::test_support;
    using uvlas::numerics::eigen::PairSide;
    const double maximum = std::numeric_limits<double>::max(),
                 tiny = std::numeric_limits<double>::denorm_min();
    t::AllocationAudit audit{t::strict_kernel_backends, t::AuditThreadDomain::CallingThread};
    const bool a = d::full_phase_coordinate_larger({maximum, maximum}, maximum, {maximum, 0},
                                                   maximum, PairSide::Left);
    const bool b =
        d::full_phase_coordinate_larger({maximum, tiny}, tiny, {maximum, 0}, tiny, PairSide::Right);
    const auto result = audit.finish();
    CHECK(a);
    CHECK(b);
    CHECK(result.status == t::AllocationAuditStatus::Clean);
    CHECK(result.cpp_allocation_attempts == 0);
}
