#include "sparse_test_support.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <uvlas/core/instrumentation.hpp>

using namespace sparse_test;

TEST_CASE("Sparse synthetic 100k tridiagonal is a capacity check only", "[N4-05][R-02]") {
    // Closed-form x_i=i+1 gives b_i=2(i+1), with the last b=3n+1.
    // This input has benign fill and does not certify a product graph or RSS cap.
    constexpr core::Index n = 100000;
    Fixture fixture{n};
    std::vector<Triplet> entries;
    entries.reserve(static_cast<std::size_t>(3 * n - 2));
    std::vector<double> right(static_cast<std::size_t>(n));
    for (core::Index i = 0; i < n; ++i) {
        entries.push_back({i, i, 4});
        if (i != 0) {
            entries.push_back({i - 1, i, -1});
        }
        if (i + 1 != n) {
            entries.push_back({i + 1, i, -1});
        }
        right[static_cast<std::size_t>(i)] = static_cast<double>(2 * (i + 1));
    }
    right.back() = static_cast<double>(3 * n + 1);
    auto matrix = fixture.matrix(n, entries);
    auto b = fixture.rhs(n, 1, right);
    const SparseOperator original{matrix, context()};
    const Rhs<double> rhs{b, context(), core::RhsGeneration{41}};
    const auto budget = limits(n, 1);
    const auto start = std::chrono::steady_clock::now();
    auto prepared = take(SparsePrepared::prepare(original, fixture.scales(),
                                                 SparseMethod::GeneralLu, budget, fixture.plan));
    const auto ready = std::chrono::steady_clock::now();
    auto demand = accuracy();
    demand.absolute_tolerance = 1e-9;
    demand.relative_tolerance = 1e-12;
    auto result = take(solve_linear(original, rhs, fixture.scales(), fixture.precision(), demand,
                                    budget, prepared, fixture.plan));
    const auto done = std::chrono::steady_clock::now();
    double maximum_relative = 0;
    for (core::Index i = 0; i < n; ++i) {
        const double expected = static_cast<double>(i + 1);
        maximum_relative =
            std::max(maximum_relative,
                     std::abs(values(result)[static_cast<std::size_t>(i)] - expected) / expected);
    }
    CHECK(maximum_relative <= 1e-12);
    auto memory = take(core::sample_process_memory(operation));
    std::cout << "capacity_only n=" << n << " nnz=" << entries.size()
              << " prepare_seconds=" << std::chrono::duration<double>(ready - start).count()
              << " solve_seconds=" << std::chrono::duration<double>(done - ready).count()
              << " lifetime_peak_rss_bytes=" << memory.peak_resident_bytes
              << " maximum_relative_error=" << maximum_relative << '\n';
}
