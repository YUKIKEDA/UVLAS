#include "reference_dependency_walk.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace uvlas;
namespace rd = io::detail;
using Index = core::Index;
using Code = core::FailureCode;

rd::DependencyEdge edge(std::optional<Index> from, Index to) noexcept {
    return {from, to, io::SourceId{from ? static_cast<std::uint64_t>(*from + 1) * 10 : 1},
            rd::ReferenceRole::ObservationReference, 102};
}

template <class T>
T take(core::Result<T> result) {
    REQUIRE(result.has_value());
    return std::move(*result);
}

rd::DependencyWalkView observe(const rd::DependencyWalk& graph) {
    if (const auto view = graph.view()) {
        return *view;
    }
    throw std::runtime_error{"Missing graph view"};
}

template <class T>
void check_failure(const core::Result<T>& result, Code expected) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == expected);
    CHECK(result.error().operation == core::Operation::PrepareModel);
}

} // namespace

TEST_CASE("Incremental dependency walk visits and finishes the frozen graph in order",
          "[IO4-01][N4-02][property][verification]") {
    // Frozen reference-snapshot-design-cases graph, with IDs 10, 20, 30 mapped to 0, 1, 2.
    constexpr std::array<Index, 3> expected_visits{0, 1, 2};
    constexpr std::array<Index, 3> expected_finishes{1, 2, 0};
    auto walk = [] {
        constexpr auto bytes = 12 * sizeof(Index);
        core::ResourcePlan plan{bytes};
        auto owned = take(rd::DependencyWalk::create({3, 5, 3, 20, 0}, plan, {}));
        CHECK(plan.allocated_bytes() == bytes);
        return owned;
    }();
    auto moved = std::move(walk);
    // NOLINTNEXTLINE(bugprone-use-after-move): the previously owning view is absent after move.
    CHECK_FALSE(walk.view());
    core::ResourcePlan no_allocations{0};
    std::vector<Index> first_visits;
    auto enter = [&](std::optional<Index> from, Index to) {
        const auto next = take(moved.enter(edge(from, to), no_allocations, {}, nullptr));
        if (next == rd::DependencyStep::Visit) {
            first_visits.push_back(to);
        }
        return next;
    };
    CHECK(enter({}, 0) == rd::DependencyStep::Visit);
    CHECK(enter(0, 1) == rd::DependencyStep::Visit);
    REQUIRE(moved.finish(1, {}));
    CHECK(enter(0, 2) == rd::DependencyStep::Visit);
    CHECK(enter(2, 1) == rd::DependencyStep::Shared);
    REQUIRE(moved.finish(2, {}));
    REQUIRE(moved.finish(0, {}));
    CHECK(enter({}, 2) == rd::DependencyStep::Shared);
    const auto view = observe(moved);
    CHECK(view.active.empty());
    CHECK(std::ranges::equal(first_visits, expected_visits));
    CHECK(std::ranges::equal(view.dependency_first, expected_finishes));
    CHECK(view.edges == 5);
    CHECK(view.work == 20);
    CHECK(no_allocations.allocated_bytes() == 0);
}

TEST_CASE("Cycle diagnostics own the closing declaration and bounded active suffix",
          "[IO4-01][N4-08][negative][verification]") {
    constexpr std::array<Index, 2> stack{0, 2};
    constexpr std::array<Index, 3> chain{0, 2, 0};
    constexpr std::array<Index, 2> self{2, 2};
    constexpr auto diagnostic_bytes = 3 * sizeof(Index);
    core::ResourcePlan plan{1024};
    auto walk = take(rd::DependencyWalk::create({3, 8, 3, 40, diagnostic_bytes}, plan, {}));
    REQUIRE(walk.enter(edge({}, 0), plan, {}, nullptr));
    REQUIRE(walk.enter(edge(0, 2), plan, {}, nullptr));
    const auto original = observe(walk);
    rd::DependencyWalkDiagnostic diagnostic;
    const auto closing = edge(2, 0);
    core::ResourcePlan enough{diagnostic_bytes};
    check_failure(walk.enter(closing, enough, {}, &diagnostic), Code::InvalidInput);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::Cycle);
    if (const auto& location = diagnostic.edge) {
        CHECK(location->declaring_source == io::SourceId{30});
        CHECK(location->role == rd::ReferenceRole::ObservationReference);
        CHECK(location->declaration_node == 102);
    } else {
        FAIL("Missing closing edge provenance");
    }
    CHECK(std::ranges::equal(diagnostic.cycle.values(), chain));
    CHECK(enough.remaining_bytes() == 0);
    check_failure(walk.enter(closing, plan, {}, nullptr), Code::InvalidInput);
    core::ResourcePlan short_plan{diagnostic_bytes - 1};
    check_failure(walk.enter(closing, short_plan, {}, &diagnostic), Code::InvalidInput);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::Cycle);
    CHECK(diagnostic.cycle.size() == 0);
    CHECK(short_plan.allocated_bytes() == 0);
    check_failure(walk.enter(edge(2, 2), plan, {}, &diagnostic), Code::InvalidInput);
    CHECK(std::ranges::equal(diagnostic.cycle.values(), self));
    CHECK(std::ranges::equal(observe(walk).active, stack));
    CHECK(observe(walk).dependency_first.empty());
    CHECK(observe(walk).edges == original.edges);
    CHECK(observe(walk).work == original.work);

    auto limited = take(rd::DependencyWalk::create({3, 8, 3, 40, diagnostic_bytes - 1}, plan, {}));
    REQUIRE(limited.enter(edge({}, 0), plan, {}, nullptr));
    REQUIRE(limited.enter(edge(0, 2), plan, {}, nullptr));
    check_failure(limited.enter(closing, plan, {}, &diagnostic), Code::InvalidInput);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::Cycle);
    CHECK(diagnostic.cycle.size() == 0);

    auto tight_work = take(rd::DependencyWalk::create({3, 8, 3, 15, diagnostic_bytes}, plan, {}));
    REQUIRE(tight_work.enter(edge({}, 0), plan, {}, nullptr));
    REQUIRE(tight_work.enter(edge(0, 2), plan, {}, nullptr));
    const auto before = plan.allocated_bytes();
    check_failure(tight_work.enter(closing, plan, {}, &diagnostic), Code::InvalidInput);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::Cycle);
    CHECK(diagnostic.cycle.size() == 0);
    CHECK(plan.allocated_bytes() == before);
    CHECK(observe(tight_work).work == 14);
}

TEST_CASE("Dependency bounds and rejected transitions preserve admitted traversal state",
          "[IO4-01][N4-02][N4-08][negative]") {
    core::ResourcePlan plan{1024};
    constexpr auto bytes = 6 * sizeof(Index);
    core::ResourcePlan one_short{bytes - 1};
    check_failure(rd::DependencyWalk::create({2, 2, 1, 7, 0}, one_short, {}),
                  Code::ResourceExceeded);
    CHECK(one_short.allocated_bytes() == 0);
    check_failure(rd::DependencyWalk::create({2, 2, 1, 5, 0}, plan, {}), Code::ResourceExceeded);
    check_failure(rd::DependencyWalk::create(
                      {std::numeric_limits<std::uint64_t>::max(), 0, 0, 0, 0}, plan, {}),
                  Code::ResourceExceeded);
    CHECK(plan.allocated_bytes() == 0);

    auto walk = take(rd::DependencyWalk::create({2, 2, 1, 9, 0}, plan, {}));
    rd::DependencyWalkDiagnostic diagnostic;
    check_failure(walk.enter(edge({}, 2), plan, {}, &diagnostic), Code::InvalidInput);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::InvalidEdge);
    REQUIRE(walk.enter(edge({}, 0), plan, {}, &diagnostic));
    CHECK_FALSE(diagnostic.reason);
    check_failure(walk.enter(edge(0, 1), plan, {}, &diagnostic), Code::ResourceExceeded);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::DepthLimit);
    CHECK(diagnostic.required == 2);
    CHECK(diagnostic.limit == 1);
    check_failure(walk.enter(edge({}, 1), plan, {}, &diagnostic), Code::InvalidInput);
    check_failure(walk.finish(1, {}), Code::InvalidInput);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop requires mutable source.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    check_failure(walk.enter(edge(0, 1), plan, stop.get_token(), &diagnostic), Code::Cancelled);
    CHECK_FALSE(diagnostic.reason);
    check_failure(walk.finish(0, stop.get_token()), Code::Cancelled);
    CHECK(observe(walk).active.size() == 1);
    CHECK(observe(walk).edges == 1);
    CHECK(observe(walk).work == 7);
    REQUIRE(walk.finish(0, {}));
    REQUIRE(walk.enter(edge({}, 0), plan, {}, &diagnostic));
    CHECK(observe(walk).edges == 2);
    CHECK(observe(walk).work == 9);
    check_failure(walk.enter(edge({}, 1), plan, {}, &diagnostic), Code::ResourceExceeded);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::EdgeLimit);

    auto short_work = take(rd::DependencyWalk::create({2, 2, 1, 7, 0}, plan, {}));
    REQUIRE(short_work.enter(edge({}, 0), plan, {}, nullptr));
    check_failure(short_work.finish(0, {}), Code::ResourceExceeded);
    CHECK(observe(short_work).active.size() == 1);
    CHECK(observe(short_work).dependency_first.empty());

    core::ResourcePlan none{0};
    const auto empty = take(rd::DependencyWalk::create({}, none, {}));
    CHECK(observe(empty).dependency_first.empty());
    CHECK(observe(empty).work == 0);
}

TEST_CASE("Shared reference depth is independent of root declaration order",
          "[IO4-01][N4-02][N4-08][negative][verification]") {
    // The longest path is 0 -> 1 -> 2 -> 3. The manifest also declares 2 directly.
    constexpr std::array<Index, 4> expected_order{3, 2, 1, 0};
    for (const bool shared_root_first : {false, true}) {
        for (const std::uint64_t depth : {3, 4}) {
            CAPTURE(shared_root_first, depth);
            core::ResourcePlan plan{(8 + 2 * depth) * sizeof(Index)};
            auto walk = take(rd::DependencyWalk::create({4, 5, depth, 25, 0}, plan, {}));
            core::ResourcePlan none{0};
            rd::DependencyWalkDiagnostic diagnostic;
            const auto enter = [&](std::optional<Index> from, Index to) {
                return walk.enter(edge(from, to), none, {}, &diagnostic);
            };
            if (shared_root_first) {
                REQUIRE(take(enter({}, 2)) == rd::DependencyStep::Visit);
                REQUIRE(take(enter(2, 3)) == rd::DependencyStep::Visit);
                REQUIRE(walk.finish(3, {}));
                REQUIRE(walk.finish(2, {}));
            }
            REQUIRE(take(enter({}, 0)) == rd::DependencyStep::Visit);
            REQUIRE(take(enter(0, 1)) == rd::DependencyStep::Visit);
            if (!shared_root_first) {
                REQUIRE(take(enter(1, 2)) == rd::DependencyStep::Visit);
            }
            if (depth == 3) {
                const auto before = observe(walk);
                const auto result = shared_root_first ? enter(1, 2) : enter(2, 3);
                check_failure(result, Code::ResourceExceeded);
                CHECK(diagnostic.reason == rd::DependencyWalkFailure::DepthLimit);
                CHECK(diagnostic.required == 4);
                CHECK(diagnostic.limit == 3);
                const auto after = observe(walk);
                constexpr std::array<Index, 3> path{0, 1, 2};
                constexpr std::array<Index, 2> finished{3, 2};
                const auto path_size = shared_root_first ? 2U : 3U;
                const auto finished_size = shared_root_first ? 2U : 0U;
                CHECK(std::ranges::equal(after.active, std::span{path}.first(path_size)));
                CHECK(std::ranges::equal(after.dependency_first,
                                         std::span{finished}.first(finished_size)));
                CHECK(after.edges == before.edges);
                CHECK(after.work == before.work);
            } else {
                if (shared_root_first) {
                    REQUIRE(take(enter(1, 2)) == rd::DependencyStep::Shared);
                } else {
                    REQUIRE(take(enter(2, 3)) == rd::DependencyStep::Visit);
                    REQUIRE(walk.finish(3, {}));
                    REQUIRE(walk.finish(2, {}));
                }
                REQUIRE(walk.finish(1, {}));
                REQUIRE(walk.finish(0, {}));
                if (!shared_root_first) {
                    REQUIRE(take(enter({}, 2)) == rd::DependencyStep::Shared);
                }
                const auto after = observe(walk);
                CHECK(after.active.empty());
                CHECK(std::ranges::equal(after.dependency_first, expected_order));
                CHECK(after.edges == 5);
                CHECK(after.work == 25);
            }
            CHECK(none.allocated_bytes() == 0);
        }
    }
}

TEST_CASE("A finished reference retains the longest of its shared dependency paths",
          "[IO4-01][N4-02][N4-08][negative]") {
    // 3 -> 4 finishes first. Source 2 then shares 3 and visits a shorter leaf.
    core::ResourcePlan plan{20 * sizeof(Index)};
    auto walk = take(rd::DependencyWalk::create({7, 8, 3, 40, 0}, plan, {}));
    core::ResourcePlan none{0};
    const auto visit = [&](std::optional<Index> from, Index to) {
        REQUIRE(take(walk.enter(edge(from, to), none, {}, nullptr)) == rd::DependencyStep::Visit);
    };
    visit({}, 3);
    visit(3, 4);
    REQUIRE(walk.finish(4, {}));
    REQUIRE(walk.finish(3, {}));
    visit({}, 2);
    REQUIRE(take(walk.enter(edge(2, 3), none, {}, nullptr)) == rd::DependencyStep::Shared);
    visit(2, 5);
    REQUIRE(walk.finish(5, {}));
    REQUIRE(walk.finish(2, {}));

    visit({}, 0);
    visit(0, 1);
    const auto before = observe(walk);
    rd::DependencyWalkDiagnostic diagnostic;
    check_failure(walk.enter(edge(1, 2), none, {}, &diagnostic), Code::ResourceExceeded);
    CHECK(diagnostic.reason == rd::DependencyWalkFailure::DepthLimit);
    CHECK(diagnostic.required == 5);
    CHECK(diagnostic.limit == 3);
    const auto after = observe(walk);
    constexpr std::array<Index, 2> expected_active{0, 1};
    constexpr std::array<Index, 4> expected_finished{4, 3, 5, 2};
    CHECK(std::ranges::equal(after.active, expected_active));
    CHECK(std::ranges::equal(after.dependency_first, expected_finished));
    CHECK(after.edges == before.edges);
    CHECK(after.work == before.work);
    CHECK(none.allocated_bytes() == 0);
}

TEST_CASE("Deep reference dependencies use an explicit stack and reverse completion order",
          "[IO4-01][N4-02][property]") {
    constexpr Index n = 2048;
    core::ResourcePlan plan{4 * n * sizeof(Index)};
    auto walk = take(rd::DependencyWalk::create({n, n, n, 6 * n, 0}, plan, {}));
    core::ResourcePlan none{0};
    REQUIRE(walk.enter(edge({}, 0), none, {}, nullptr));
    for (Index to = 1; to < n; ++to) {
        REQUIRE(walk.enter(edge(to - 1, to), none, {}, nullptr));
    }
    CHECK(observe(walk).active.size() == n);
    for (Index to = n; to > 0; --to) {
        REQUIRE(walk.finish(to - 1, {}));
    }
    const auto output = observe(walk);
    REQUIRE(output.dependency_first.size() == n);
    for (Index i = 0; i < n; ++i) {
        CHECK(output.dependency_first[static_cast<std::size_t>(i)] == n - i - 1);
    }
    CHECK(output.work == 6 * n);
    CHECK(none.allocated_bytes() == 0);
}
