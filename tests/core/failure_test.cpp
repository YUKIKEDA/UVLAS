#include <catch2/catch_test_macros.hpp>

#include <array>
#include <concepts>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <uvlas/core/failure.hpp>

using namespace uvlas::core;

namespace {
[[nodiscard]] FailureDetails context() {
    FailureDetails result;
    result.model = ModelGeneration{11};
    result.parent = AcceptedStateId{12};
    result.trial = TrialId{13};
    result.numerical_generation = NumericalGeneration{14};
    result.component = ComponentId{2};
    result.problem = PanelId{42};
    result.resources = ResourceUsage{ResourceKind::Bytes, 1000, 900};
    result.last_accepted = AcceptedStateId{12};
    return result;
}

template <class T, class... Others>
consteval bool separate_from_all() {
    return ((std::same_as<T, Others> || !std::is_constructible_v<T, Others>) && ...);
}

template <class... T>
consteval bool separate_identifiers(std::tuple<T...>*) {
    return (separate_from_all<T, T...>() && ...);
}
} // namespace

TEST_CASE("failure codes retain their distinct identity and owned context", "[unit][N4-08]") {
    constexpr std::array codes{
        FailureCode::InvalidInput,          FailureCode::Unsupported,
        FailureCode::GeometryInvalid,       FailureCode::RankUnresolved,
        FailureCode::ChartSingular,         FailureCode::Nonfinite,
        FailureCode::FactorizationFailed,   FailureCode::LinearUnconverged,
        FailureCode::NonlinearUnconverged,  FailureCode::EigenCoverageUnresolved,
        FailureCode::PrecisionUnsupported,  FailureCode::PrecisionUnresolved,
        FailureCode::ResourceExceeded,      FailureCode::Cancelled,
        FailureCode::StaleGeneration,       FailureCode::CheckpointCorrupt,
        FailureCode::CheckpointIncompatible};
    STATIC_REQUIRE(codes.size() == 17);
    for (std::size_t first = 0; first < codes.size(); ++first) {
        for (std::size_t second = first + 1; second < codes.size(); ++second) {
            CHECK(codes[first] != codes[second]);
        }
    }
    const auto original = context();
    for (const auto code : codes) {
        const Result<int> result = std::unexpected(Failure{code, Operation::SolveLinear, original});
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error().code == code);
        CHECK(result.error().operation == Operation::SolveLinear);
        CHECK(result.error().details == original);
    }
    auto copy = Failure{FailureCode::Cancelled, Operation::AdvanceTrial, original};
    copy.details.problem = DofId{99};
    CHECK(original.problem == std::optional<ProblemId>{PanelId{42}});
    CHECK(copy.details.problem != original.problem);
}

TEST_CASE("diagnostic identifier kinds cannot be implicitly or explicitly mixed", "[unit][N4-02]") {
    using Ids =
        std::tuple<DofId, PanelId, WakeVertexId, WakeFaceId, LayerId, TrialId, ModelGeneration,
                   NumericalGeneration, AcceptedStateId, ObservationId, ScopeId, ComponentId>;
    STATIC_REQUIRE(separate_identifiers(static_cast<Ids*>(nullptr)));
    STATIC_REQUIRE_FALSE(std::is_convertible_v<std::uint64_t, ModelGeneration>);
    STATIC_REQUIRE_FALSE(std::is_convertible_v<ModelGeneration, std::uint64_t>);
    CHECK(ModelGeneration{0}.value() == 0);
    CHECK(ModelGeneration{std::numeric_limits<std::uint64_t>::max()}.value() ==
          std::numeric_limits<std::uint64_t>::max());
    STATIC_REQUIRE(std::is_nothrow_copy_constructible_v<Failure>);
    STATIC_REQUIRE_FALSE(std::is_default_constructible_v<Failure>);
    STATIC_REQUIRE_FALSE(std::is_default_constructible_v<ResourceUsage>);
}

TEST_CASE("entry guard preserves values and existing expected failures", "[unit][N4-08]") {
    const auto details = context();
    const auto success = guard_entry(Operation::PrepareModel, FailureCode::InvalidInput, details,
                                     []() -> Result<int> { return 23; });
    REQUIRE(success.has_value());
    CHECK(*success == 23);
    const Failure original{FailureCode::StaleGeneration, Operation::Commit, details};
    const auto failure =
        guard_entry(Operation::PrepareModel, FailureCode::InvalidInput, {},
                    [&original]() -> Result<int> { return std::unexpected(original); });
    REQUIRE_FALSE(failure.has_value());
    CHECK(failure.error() == original);
    const auto empty = guard_entry(Operation::Commit, FailureCode::InvalidInput, {},
                                   []() -> Result<void> { return {}; });
    CHECK(empty.has_value());
}

TEST_CASE("entry guard translates allocation standard and nonstandard exceptions",
          "[negative][N4-08]") {
    const auto details = context();
    const auto allocation = guard_entry(Operation::SolveLinear, FailureCode::FactorizationFailed,
                                        details, []() -> Result<int> { throw std::bad_alloc{}; });
    REQUIRE_FALSE(allocation.has_value());
    CHECK(allocation.error() ==
          Failure{FailureCode::ResourceExceeded, Operation::SolveLinear, details});
    const auto standard =
        guard_entry(Operation::SolveLinear, FailureCode::FactorizationFailed, details,
                    []() -> Result<void> { throw std::runtime_error("library failure"); });
    REQUIRE_FALSE(standard.has_value());
    CHECK(standard.error() ==
          Failure{FailureCode::FactorizationFailed, Operation::SolveLinear, details});
    const auto nonstandard = guard_entry(Operation::Recover, FailureCode::Unsupported, details,
                                         []() -> Result<int> { throw 17; });
    REQUIRE_FALSE(nonstandard.has_value());
    CHECK(nonstandard.error() == Failure{FailureCode::Unsupported, Operation::Recover, details});
    CHECK(details == context());
    STATIC_REQUIRE(noexcept(guard_entry(Operation::Commit, FailureCode::InvalidInput, details,
                                        []() -> Result<void> { return {}; })));
}
