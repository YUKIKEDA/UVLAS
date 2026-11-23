#include "pg_partition_test_support.hpp"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <uvlas/structure/pg_rank_preparation.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using partition_test::Fixture;
using partition_test::take;
constexpr auto cap = std::uint64_t{48} << 30U;

constinit std::atomic<std::int64_t> arrays_before_failure{-1};

// This executable replaces only the ordinary array path used by the product's
// natural-aligned PlannedBuffer element arrays. Disarm before unwinding.
void maybe_fail_array() {
    const auto remaining = arrays_before_failure.load(std::memory_order_relaxed);
    if (remaining < 0)
        return;
    if (remaining == 0) {
        arrays_before_failure.store(-1, std::memory_order_relaxed);
        throw std::bad_alloc{};
    }
    arrays_before_failure.fetch_sub(1, std::memory_order_relaxed);
}

struct ArrayFault {
    explicit ArrayFault(std::int64_t count) noexcept {
        arrays_before_failure.store(count, std::memory_order_relaxed);
    }
    ~ArrayFault() {
        arrays_before_failure.store(-1, std::memory_order_relaxed);
    }
    ArrayFault(const ArrayFault&) = delete;
    ArrayFault& operator=(const ArrayFault&) = delete;
};
} // namespace

void* operator new[](std::size_t bytes) {
    maybe_fail_array();
    // NOLINTNEXTLINE(misc-const-correctness): operator new returns writable storage.
    void* result = std::malloc(bytes == 0 ? 1 : bytes);
    if (result != nullptr)
        return result;
    throw std::bad_alloc{};
}
void operator delete[](void* pointer) noexcept {
    std::free(pointer);
}
void operator delete[](void* pointer, std::size_t) noexcept {
    std::free(pointer);
}

TEST_CASE("PG local rank abandons a partial result after allocation failure",
          "[SH-D03][N4-08][negative][verification]") {
    Fixture f;
    core::ResourcePlan producers{cap};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    const auto dictionary = take(PgCoefficientDictionary::prepare(
        *f.graph, {partition_test::exact_limits(), 2000, 20000000, 100000, {}}, workspace,
        producers));
    const auto partition = take(PgIntegrationPartition::prepare(
        f.request(), partition_test::limits(), workspace, producers));
    const auto saved = take(partition.identity());
    const PgRankPreparationLimits limits{cap, 10000000000000ULL, 128, {}};

    core::ResourcePlan control{cap};
    const auto old =
        take(PgRankPreparation::prepare(dictionary, partition, limits, workspace, control));
    const auto previous = take(old.preparation());
    for (const std::int64_t after : {0, 43}) {
        CAPTURE(after);
        core::ResourcePlan trial{cap};
        const auto attempt = [&] {
            const ArrayFault fault(after);
            return PgRankPreparation::prepare(dictionary, partition, limits, workspace, trial);
        }();
        REQUIRE_FALSE(attempt);
        CHECK(attempt.error().code == core::FailureCode::ResourceExceeded);
        if (after == 0)
            CHECK(trial.allocated_bytes() == 0);
        else
            CHECK(trial.allocated_bytes() > 0);
        CHECK(take(partition.identity()) == saved);
        CHECK(take(old.preparation()).steps == previous.steps);
        core::ResourcePlan retry_plan{cap};
        const auto retry =
            take(PgRankPreparation::prepare(dictionary, partition, limits, workspace, retry_plan));
        CHECK(take(retry.preparation()).steps == previous.steps);
    }
}
