#include "../support/allocation_audit.hpp"
#include "sparse_test_support.hpp"
#include "sparse_triangular.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <uvlas/core/resource_plan.hpp>

TEST_CASE("Prepared sparse triangular kernel allocates no C++ or Eigen storage", "[N4-05][N4-09]") {
    using namespace uvlas;
    using namespace numerics;
    using namespace numerics::linear::detail;
    using namespace test_support;
    core::ResourcePlan plan{std::uint64_t{1024} * 1024};
    constexpr auto operation = core::Operation::SolveLinear;
    const std::array<Triplet, 4> entries{{{0, 0, 2}, {1, 0, 1}, {0, 1, 1}, {1, 1, -1}}};
    auto original = SparseReal::from_triplets(2, 2, entries, plan, operation);
    REQUIRE(original);
    const Eigen::Map<const Eigen::SparseMatrix<double, Eigen::ColMajor, core::Index>> mapped{
        2,
        2,
        4,
        original->column_offsets().data(),
        original->row_indices().data(),
        original->values().data()};
    SparseLuFactor factor;
    factor.compute(mapped);
    REQUIRE(factor.info() == Eigen::Success);
    auto factors = extract_sparse_lu(factor, plan);
    REQUIRE(factors);
    const std::array rhs{3.0, 0.0};
    std::array<double, 2> work{}, output{};
    AllocationAudit audit{strict_kernel_backends, AuditThreadDomain::CallingThread};
    auto solved = solve_sparse_triangular(*factors, rhs, 1, work, output);
    const auto observed = audit.finish();
    REQUIRE(solved);
    CHECK(observed.status == AllocationAuditStatus::Clean);
    CHECK(observed.cpp_allocation_attempts == 0);
    CHECK(output[0] == 1);
    CHECK(output[1] == 1);
}

TEST_CASE("Sparse refinement adds no allocations beyond fixed owning result preparation",
          "[N4-05][N4-06][N4-09]") {
    using namespace sparse_test;
    using namespace test_support;
    Fixture fixture{1};
    const std::array<Triplet, 1> entries{{{0, 0, 3}}};
    auto matrix = fixture.matrix(1, entries);
    const SparseOperator original{matrix, context()};
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, limits(), fixture.plan));
    const std::array right{1.0};
    auto b = fixture.rhs(1, 1, right);
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    AllocationAudit initial{strict_kernel_backends, AuditThreadDomain::CallingThread};
    auto binary = solve_linear(original, rhs, fixture.scales(), fixture.precision(), accuracy(),
                               limits(), prepared, fixture.plan);
    const auto binary_allocations = initial.finish();
    REQUIRE(binary);
    auto demand = accuracy();
    demand.absolute_tolerance = 1e-30;
    demand.componentwise_backward_tolerance = 1e-30;
    AllocationAudit refined{strict_kernel_backends, AuditThreadDomain::CallingThread};
    auto pair = solve_linear(original, rhs, fixture.scales(), fixture.precision(), demand, limits(),
                             prepared, fixture.plan);
    const auto pair_allocations = refined.finish();
    REQUIRE(pair);
    CHECK(pair->values().precision() == core::Precision::TwoWord);
    CHECK(binary_allocations.status == AllocationAuditStatus::AllocationObserved);
    CHECK(binary_allocations.cpp_allocation_attempts == pair_allocations.cpp_allocation_attempts);
    // Six metadata arrays and two owning solution planes are prepared before either kernel.
    CHECK(pair_allocations.cpp_allocation_attempts == 8);
}
