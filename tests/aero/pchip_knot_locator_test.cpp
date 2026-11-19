#include "aero/pchip_interval_jet.hpp"
#include "aero/pchip_knot_locator.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <fstream>
#include <future>
#include <map>
#include <vector>
namespace {
namespace p = uvlas::aero::detail::pchip;
namespace core = uvlas::core;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
template <class T>
T take(core::Result<T> result) {
    if (!result)
        FAIL("Unexpected failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}
struct Input {
    int sign = 0;
    std::vector<std::uint32_t> n, d;
    explicit Input(const Rat& value) {
        auto numerator = boost::multiprecision::numerator(value);
        sign = numerator == 0 ? 0 : numerator < 0 ? -1 : 1;
        if (numerator < 0)
            numerator = -numerator;
        boost::multiprecision::export_bits(numerator, std::back_inserter(n), 32, false);
        if (numerator == 0)
            n.clear();
        const auto denominator = boost::multiprecision::denominator(value);
        boost::multiprecision::export_bits(denominator, std::back_inserter(d), 32, false);
    }
    [[nodiscard]] p::RationalView view() const {
        return {sign, n, d};
    }
};
p::RationalLimits limits() {
    return {4096, 1000000, UINT64_MAX, 100000};
}
Rat decoded(p::RationalView value) {
    Big n = 0, d = 0;
    for (auto i = value.numerator.size(); i > 0; --i) {
        n <<= 32;
        n += value.numerator[i - 1];
    }
    for (auto i = value.denominator.size(); i > 0; --i) {
        d <<= 32;
        d += value.denominator[i - 1];
    }
    if (value.denominator.empty() || d == 0) {
        FAIL_CHECK("Missing canonical positive denominator");
        return 0;
    }
    CHECK(value.denominator.back() != 0);
    if (value.numerator.empty()) {
        CHECK(value.sign == 0);
        CHECK(d == 1);
    } else {
        CHECK(value.numerator.back() != 0);
        CHECK((value.sign == -1 || value.sign == 1));
    }
    Big a = n, b = d;
    while (b != 0) {
        const Big remainder = a % b;
        a = b;
        b = remainder;
    }
    CHECK(a == 1);
    if (value.sign < 0)
        n = -n;
    return Rat{n} / d;
}
p::CoefficientBlock coefficients(const std::vector<Rat>& times, const std::vector<Rat>& values,
                                 std::uint64_t m = 1) {
    std::vector<Input> t, y;
    std::vector<p::RationalView> tv, yv;
    t.reserve(times.size());
    y.reserve(values.size());
    for (const auto& v : times)
        t.emplace_back(v);
    for (const auto& v : values)
        y.emplace_back(v);
    for (const auto& v : t)
        tv.push_back(v.view());
    for (const auto& v : y)
        yv.push_back(v.view());
    const p::SampleGrid samples{tv, yv, m};
    const auto plan = take(p::forecast_coefficients(samples, limits()));
    core::ResourcePlan memory{plan.payload_bytes()};
    p::CoefficientWorkspace workspace;
    return take(p::prepare_coefficients(samples, workspace, limits(), memory));
}
} // namespace
TEST_CASE("Private PCHIP locator matches independent exact linear search",
          "[AE-07][N4-02][N4-03][N4-06][verification]") {
    std::ifstream source_stream{UVLAS_PCHIP_COEFFICIENT_REFERENCE},
        query_stream{UVLAS_PCHIP_INTERVAL_REFERENCE};
    REQUIRE(source_stream.good());
    REQUIRE(query_stream.good());
    const auto raw = nlohmann::json::parse(source_stream),
               queries = nlohmann::json::parse(query_stream);
    std::map<std::string, p::CoefficientBlock> sources;
    std::map<std::string, std::vector<Rat>> expected_times;
    for (const auto& record : raw.at("cases")) {
        std::vector<Rat> times, values;
        for (const auto& value : record.at("times"))
            times.emplace_back(value.get<std::string>());
        for (const auto& value : record.at("values"))
            values.emplace_back(value.get<std::string>());
        const auto name = record.at("name").get<std::string>();
        sources.try_emplace(
            name, coefficients(times, values, record.at("components").get<std::uint64_t>()));
        expected_times.emplace(name, std::move(times));
    }
    for (const auto& record : queries.at("cases")) {
        CAPTURE(record.at("source_name"), record.at("time"));
        const auto name = record.at("source_name").get<std::string>();
        const auto& source = sources.at(name);
        const auto& times = expected_times.at(name);
        const Rat query{record.at("time").get<std::string>()};
        const Input input{query};
        // Linear scan over independent cpp_rational values, never the product's time comparisons.
        std::size_t position = 0;
        while (position < times.size() && times[position] < query)
            ++position;
        REQUIRE(position < times.size());
        const bool exact = times[position] == query;
        const auto plan = take(p::forecast_location(source, input.view(), limits()));
        core::ResourcePlan memory{plan.payload_bytes()};
        p::LocatorWorkspace workspace;
        const auto found = take(p::locate_time(source, input.view(), workspace, limits(), memory));
        CHECK(found.is_knot() == exact);
        CHECK(found.nodes() == times.size());
        CHECK(found.comparisons() <= plan.maximum_comparisons());
        CHECK(found.used_steps() <= plan.reserved_steps());
        CHECK(found.used_slots() <= plan.arithmetic().slots());
        CHECK(memory.allocated_bytes() == plan.payload_bytes());
        if (exact)
            CHECK(take(found.knot_index()) == position);
        else
            CHECK_FALSE(found.knot_index());
        const bool left = !exact || position > 0, right = !exact || position + 1 < times.size();
        CHECK(found.has_left() == left);
        CHECK(found.has_right() == right);
        for (const auto side : {p::KnotSide::Left, p::KnotSide::Right}) {
            const auto actual = found.choose(side);
            const bool present = side == p::KnotSide::Left ? left : right;
            if (present) {
                REQUIRE(actual);
                const auto expected = exact && side == p::KnotSide::Right ? position : position - 1;
                CHECK(*actual == expected);
            } else {
                REQUIRE_FALSE(actual);
                CHECK(actual.error().code == core::FailureCode::InvalidInput);
            }
        }
        const auto unique = found.unique_interval();
        if (!exact || !left || !right) {
            REQUIRE(unique);
            CHECK(*unique == (!exact ? position - 1 : position == 0 ? 0 : position - 1));
        } else {
            REQUIRE_FALSE(unique);
            CHECK(unique.error().code == core::FailureCode::InvalidInput);
        }
    }
}
TEST_CASE("Private PCHIP locator bounds binary probes across powers of two",
          "[N4-02][N4-06][verification]") {
    for (const auto n : std::array<std::uint64_t, 9>{2, 3, 4, 7, 8, 9, 15, 16, 17}) {
        std::vector<Rat> times, values;
        for (std::uint64_t i = 0; i < n; ++i) {
            times.emplace_back(i);
            values.emplace_back(0);
        }
        const auto source = coefficients(times, values);
        p::LocatorWorkspace workspace;
        for (std::uint64_t j = 0; j < 2 * n - 1; ++j) {
            const Input query{Rat{j} / 2};
            const auto plan = take(p::forecast_location(source, query.view(), limits()));
            core::ResourcePlan memory{plan.payload_bytes()};
            const auto result =
                take(p::locate_time(source, query.view(), workspace, limits(), memory));
            CHECK(plan.maximum_comparisons() == static_cast<std::uint64_t>(std::bit_width(n)) + 3);
            CHECK(result.comparisons() <= plan.maximum_comparisons());
            CHECK(result.is_knot() == (j % 2 == 0));
            if (j % 2 == 0)
                CHECK(take(result.knot_index()) == j / 2);
            else
                CHECK(take(result.unique_interval()) == j / 2);
        }
    }
}
TEST_CASE("Private PCHIP locator connects explicit knot sides to existing interval evaluation",
          "[AE-07][N4-03][N4-06][verification]") {
    const auto source = coefficients({0, 1, 2}, {0, 1, 4});
    const Input time{Rat{1}};
    const auto plan = take(p::forecast_location(source, time.view(), limits()));
    core::ResourcePlan memory{plan.payload_bytes()};
    p::LocatorWorkspace workspace;
    const auto at = take(p::locate_time(source, time.view(), workspace, limits(), memory));
    CHECK(at.is_knot());
    CHECK_FALSE(at.unique_interval());
    const auto evaluate = [&](p::KnotSide side) {
        const p::IntervalRequest request{take(at.choose(side)), time.view(), p::JetOrder::Second};
        const auto interval_plan = take(p::forecast_interval(source, request, limits()));
        core::ResourcePlan result_memory{interval_plan.payload_bytes()};
        p::IntervalWorkspace interval_workspace;
        return take(
            p::evaluate_interval(source, request, interval_workspace, limits(), result_memory));
    };
    const auto left = evaluate(p::KnotSide::Left), right = evaluate(p::KnotSide::Right);
    CHECK(take(left.interval()) == 0);
    CHECK(take(right.interval()) == 1);
    CHECK(decoded(take(left.derivative(0, p::JetOrder::Value))) == 1);
    CHECK(decoded(take(right.derivative(0, p::JetOrder::Value))) == 1);
    CHECK(decoded(take(left.derivative(0, p::JetOrder::First))) == Rat{3} / 2);
    CHECK(decoded(take(right.derivative(0, p::JetOrder::First))) == Rat{3} / 2);
    CHECK(decoded(take(left.derivative(0, p::JetOrder::Second))) == 0);
    CHECK(decoded(take(right.derivative(0, p::JetOrder::Second))) == 4);
}
TEST_CASE("Private PCHIP locator rejects malformed inputs and keeps scalar results",
          "[N4-03][N4-08][negative]") {
    const auto source = coefficients({0, 1, 2}, {0, 1, 4});
    const Input input{Rat{1} / 2}, below{Rat{-1} / 2}, above{Rat{5} / 2};
    const auto plan = take(p::forecast_location(source, input.view(), limits()));
    p::LocatorWorkspace workspace;
    core::ResourcePlan memory{plan.payload_bytes()};
    const auto old = take(p::locate_time(source, input.view(), workspace, limits(), memory));
    const std::array<std::uint32_t, 1> one{1}, two{2}, four{4}, zero{0};
    const std::array<std::uint32_t, 2> leading{1, 0};
    const std::array<p::RationalView, 8> bad{{below.view(),
                                              above.view(),
                                              {1, two, four},
                                              {1, one, {}},
                                              {1, {}, one},
                                              {0, one, one},
                                              {1, leading, one},
                                              {1, one, zero}}};
    core::FailureDetails details;
    details.component = core::ComponentId{17};
    details.parent = core::AcceptedStateId{18};
    for (const auto query : bad) {
        core::ResourcePlan candidate{UINT64_MAX};
        const auto failed =
            p::locate_time(source, query, workspace, limits(), candidate, {}, details);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == core::FailureCode::InvalidInput);
        CHECK(failed.error().details.component == details.component);
        CHECK(failed.error().details.parent == details.parent);
        CHECK(take(old.unique_interval()) == 0);
    }
    CHECK_FALSE(old.choose(static_cast<p::KnotSide>(255)));
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() is non-const.
    std::stop_source stop;
    stop.request_stop();
    core::ResourcePlan cancelled_memory{plan.payload_bytes()};
    const auto cancelled = p::locate_time(source, input.view(), workspace, limits(),
                                          cancelled_memory, stop.get_token(), details);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details);
    CHECK(cancelled_memory.allocated_bytes() == 0);
    const auto kept = [] {
        const auto parent = coefficients({0, 1, 2}, {0, 1, 4});
        const Input time{Rat{1}};
        const auto counted = take(p::forecast_location(parent, time.view(), limits()));
        p::LocatorWorkspace scratch;
        core::ResourcePlan arena{counted.payload_bytes()};
        return take(p::locate_time(parent, time.view(), scratch, limits(), arena));
    }();
    CHECK(take(kept.knot_index()) == 1);
    CHECK(take(kept.left_interval()) == 0);
    CHECK(take(kept.right_interval()) == 1);
    auto moved_source = coefficients({0, 1}, {0, 1});
    auto retained = std::move(moved_source);
    // NOLINTBEGIN(bugprone-use-after-move): inactive source wrappers are explicitly rejected.
    const auto inactive = p::forecast_location(moved_source, input.view(), limits());
    // NOLINTEND(bugprone-use-after-move)
    REQUIRE_FALSE(inactive);
    CHECK(inactive.error().code == core::FailureCode::InvalidInput);
    CHECK(retained.active());
}
TEST_CASE("Private PCHIP locator callers return identical scalar locations and counters",
          "[N4-03][N4-09][property]") {
    const auto source = coefficients({0, 1, 3, 4}, {0, 2, -1, 4});
    const Input query{Rat{4} / 3};
    using Encoded = std::array<std::uint64_t, 7>;
    const auto evaluate = [&]() -> std::optional<Encoded> {
        const auto plan = p::forecast_location(source, query.view(), limits());
        if (!plan)
            return std::nullopt;
        core::ResourcePlan memory{plan->payload_bytes()};
        p::LocatorWorkspace workspace;
        const auto found = p::locate_time(source, query.view(), workspace, limits(), memory);
        if (!found)
            return std::nullopt;
        const auto interval = found->unique_interval();
        if (!interval)
            return std::nullopt;
        return Encoded{found->nodes(),          *interval,           found->is_knot() ? 1U : 0U,
                       found->comparisons(),    found->used_steps(), found->used_slots(),
                       memory.allocated_bytes()};
    };
    const auto baseline = evaluate();
    if (!baseline) {
        FAIL_CHECK("Missing locator baseline");
        return;
    }
    const auto& expected = baseline.value();
    for (const auto workers : std::array<std::size_t, 3>{1, 2, 4}) {
        std::vector<std::future<std::optional<Encoded>>> pending;
        pending.reserve(workers);
        for (std::size_t i = 0; i < workers; ++i)
            pending.push_back(std::async(std::launch::async, evaluate));
        for (auto& task : pending) {
            const auto result = task.get();
            if (!result) {
                FAIL_CHECK("Missing independent locator");
                continue;
            }
            const auto& actual = result.value();
            CHECK(actual == expected);
        }
    }
}
