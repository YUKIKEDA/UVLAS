#include "aero/pchip_interval_jet.hpp"

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
TEST_CASE("Private PCHIP interval values match independent Hermite time derivatives",
          "[AE-07][N4-02][N4-03][N4-06][verification]") {
    std::ifstream source_stream{UVLAS_PCHIP_COEFFICIENT_REFERENCE},
        query_stream{UVLAS_PCHIP_INTERVAL_REFERENCE};
    REQUIRE(source_stream.good());
    REQUIRE(query_stream.good());
    const auto source_reference = nlohmann::json::parse(source_stream),
               reference = nlohmann::json::parse(query_stream);
    std::map<std::string, p::CoefficientBlock> sources;
    for (const auto& record : source_reference.at("cases")) {
        std::vector<Rat> times, values;
        for (const auto& v : record.at("times"))
            times.emplace_back(v.get<std::string>());
        for (const auto& v : record.at("values"))
            values.emplace_back(v.get<std::string>());
        sources.try_emplace(
            record.at("name").get<std::string>(),
            coefficients(times, values, record.at("components").get<std::uint64_t>()));
    }
    const std::array<const char*, 3> names{"value", "dt", "dtt"};
    for (const auto& record : reference.at("cases")) {
        CAPTURE(record.at("source_name"), record.at("interval"), record.at("u"));
        const auto& source = sources.at(record.at("source_name").get<std::string>());
        const Input time{Rat{record.at("time").get<std::string>()}};
        for (const auto order : {p::JetOrder::Value, p::JetOrder::First, p::JetOrder::Second}) {
            CAPTURE(static_cast<int>(order));
            const p::IntervalRequest query{record.at("interval").get<std::uint64_t>(), time.view(),
                                           order};
            const auto plan = take(p::forecast_interval(source, query, limits()));
            core::ResourcePlan memory{plan.payload_bytes()};
            p::IntervalWorkspace workspace;
            const auto value =
                take(p::evaluate_interval(source, query, workspace, limits(), memory));
            CHECK(value.components() == source.components());
            CHECK(take(value.interval()) == query.interval);
            CHECK(take(value.requested_order()) == order);
            CHECK(memory.allocated_bytes() == plan.payload_bytes());
            CHECK(value.used_steps() <= plan.reserved_steps());
            CHECK(value.used_slots() == plan.arithmetic().slots());
            CHECK(decoded(take(value.time())) == Rat{record.at("time").get<std::string>()});
            CHECK(decoded(take(value.coordinate_u())) == Rat{record.at("u").get<std::string>()});
            for (std::uint64_t c = 0; c < value.components(); ++c)
                for (std::uint64_t d = 0; d < 3; ++d) {
                    const auto actual = value.derivative(c, static_cast<p::JetOrder>(d));
                    if (d <= static_cast<std::uint64_t>(order)) {
                        REQUIRE(actual);
                        CHECK(decoded(*actual) == Rat{record.at("values")
                                                          .at(static_cast<std::size_t>(c))
                                                          .at(names[static_cast<std::size_t>(d)])
                                                          .get<std::string>()});
                    } else {
                        REQUIRE_FALSE(actual);
                        CHECK(actual.error().code == core::FailureCode::InvalidInput);
                    }
                }
        }
    }
}
TEST_CASE("Private PCHIP interval requests reject invalid shape canonicality and range",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const auto source = coefficients({0, 1, 2}, {0, 1, 4});
    const Input good{Rat{1} / 2}, below{Rat{-1} / 2}, beyond{Rat{3} / 2};
    const p::IntervalRequest query{0, good.view(), p::JetOrder::Second};
    const auto plan = take(p::forecast_interval(source, query, limits()));
    p::IntervalWorkspace workspace;
    core::ResourcePlan memory{plan.payload_bytes()};
    const auto old = take(p::evaluate_interval(source, query, workspace, limits(), memory));
    const auto borrowed = take(old.derivative(0, p::JetOrder::Value));
    const std::array<std::uint32_t, 1> one{1}, two{2}, four{4}, zero{0};
    const std::array<std::uint32_t, 2> leading{1, 0};
    const std::array<p::IntervalRequest, 10> invalid{{{0, below.view(), p::JetOrder::Second},
                                                      {0, beyond.view(), p::JetOrder::Value},
                                                      {2, good.view(), p::JetOrder::Value},
                                                      {UINT64_MAX, good.view(), p::JetOrder::Value},
                                                      {0, good.view(), static_cast<p::JetOrder>(3)},
                                                      {0, {1, two, four}, p::JetOrder::First},
                                                      {0, {1, one, {}}, p::JetOrder::Value},
                                                      {0, {1, {}, one}, p::JetOrder::Value},
                                                      {0, {1, leading, one}, p::JetOrder::Value},
                                                      {0, {1, one, zero}, p::JetOrder::Value}}};
    core::FailureDetails details;
    details.component = core::ComponentId{17};
    details.parent = core::AcceptedStateId{18};
    for (const auto& bad : invalid) {
        core::ResourcePlan candidate{UINT64_MAX};
        const auto result =
            p::evaluate_interval(source, bad, workspace, limits(), candidate, {}, details);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::InvalidInput);
        CHECK(result.error().details.component == details.component);
        CHECK(result.error().details.parent == details.parent);
        CHECK(decoded(borrowed) == Rat{5} / 16);
    }
    CHECK_FALSE(old.derivative(1, p::JetOrder::Value));
    CHECK_FALSE(old.derivative(0, static_cast<p::JetOrder>(255)));
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() is non-const.
    std::stop_source stop;
    stop.request_stop();
    core::ResourcePlan cancelled_memory{plan.payload_bytes()};
    const auto cancelled = p::evaluate_interval(source, query, workspace, limits(),
                                                cancelled_memory, stop.get_token(), details);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK(cancelled.error().details == details);
    CHECK(cancelled_memory.allocated_bytes() == 0);
    core::ResourcePlan recovered_memory{plan.payload_bytes()};
    const auto recovered =
        take(p::evaluate_interval(source, query, workspace, limits(), recovered_memory));
    CHECK(decoded(take(recovered.derivative(0, p::JetOrder::First))) == Rat{9} / 8);
}
TEST_CASE("Private PCHIP interval owners remain valid after source destruction and moves",
          "[N4-03][N4-08][property]") {
    auto result = [] {
        const auto source = coefficients({0, 1, 2}, {0, 1, 4});
        const Input time{Rat{1} / 2};
        const p::IntervalRequest query{0, time.view(), p::JetOrder::Second};
        const auto plan = take(p::forecast_interval(source, query, limits()));
        core::ResourcePlan memory{plan.payload_bytes()};
        p::IntervalWorkspace workspace;
        return take(p::evaluate_interval(source, query, workspace, limits(), memory));
    }();
    const auto borrowed = take(result.derivative(0, p::JetOrder::Second));
    auto moved = std::move(result);
    // NOLINTBEGIN(bugprone-use-after-move): explicitly supported inactive wrapper queries.
    CHECK_FALSE(result.active());
    CHECK(result.components() == 0);
    CHECK(result.used_steps() == 0);
    CHECK(result.used_slots() == 0);
    CHECK_FALSE(result.interval());
    CHECK_FALSE(result.requested_order());
    CHECK_FALSE(result.time());
    CHECK_FALSE(result.coordinate_u());
    CHECK_FALSE(result.derivative(0, p::JetOrder::Value));
    // NOLINTEND(bugprone-use-after-move)
    CHECK(decoded(borrowed) == Rat{3} / 2);
    CHECK(decoded(take(moved.time())) == Rat{1} / 2);
    CHECK(decoded(take(moved.coordinate_u())) == Rat{1} / 2);
    auto source = coefficients({0, 1}, {0, 1});
    auto source_owner = std::move(source);
    const Input time{Rat{1} / 2};
    // NOLINTBEGIN(bugprone-use-after-move): evaluation rejects inactive coefficient wrappers.
    const auto inactive =
        p::forecast_interval(source, {0, time.view(), p::JetOrder::Value}, limits());
    // NOLINTEND(bugprone-use-after-move)
    REQUIRE_FALSE(inactive);
    CHECK(inactive.error().code == core::FailureCode::InvalidInput);
    CHECK(source_owner.active());
}
TEST_CASE("Private PCHIP interval callers produce identical canonical words and resource counters",
          "[N4-03][N4-09][property]") {
    const auto source = coefficients({0, 1, 3}, {0, 2, -1});
    const Input time{Rat{4} / 3};
    const p::IntervalRequest query{1, time.view(), p::JetOrder::Second};
    struct Encoded {
        std::vector<std::uint32_t> words;
        std::uint64_t work, slots, payload;
        bool operator==(const Encoded&) const = default;
    };
    const auto evaluate = [&]() -> std::optional<Encoded> {
        const auto plan = p::forecast_interval(source, query, limits());
        if (!plan)
            return std::nullopt;
        core::ResourcePlan memory{plan->payload_bytes()};
        p::IntervalWorkspace workspace;
        const auto value = p::evaluate_interval(source, query, workspace, limits(), memory);
        if (!value)
            return std::nullopt;
        Encoded result{{}, value->used_steps(), value->used_slots(), memory.allocated_bytes()};
        for (const auto order : {p::JetOrder::Value, p::JetOrder::First, p::JetOrder::Second}) {
            const auto view = value->derivative(0, order);
            if (!view)
                return std::nullopt;
            result.words.push_back(static_cast<std::uint32_t>(view->sign + 1));
            result.words.push_back(static_cast<std::uint32_t>(view->numerator.size()));
            result.words.push_back(static_cast<std::uint32_t>(view->denominator.size()));
            result.words.insert(result.words.end(), view->numerator.begin(), view->numerator.end());
            result.words.insert(result.words.end(), view->denominator.begin(),
                                view->denominator.end());
        }
        return result;
    };
    const auto baseline = evaluate();
    if (!baseline) {
        FAIL_CHECK("Missing independent baseline");
        return;
    }
    const auto& expected = baseline.value();
    for (const auto workers : std::array<std::size_t, 3>{1, 2, 4}) {
        std::vector<std::future<std::optional<Encoded>>> futures;
        futures.reserve(workers);
        for (std::size_t i = 0; i < workers; ++i)
            futures.push_back(std::async(std::launch::async, evaluate));
        for (auto& pending : futures) {
            const auto result = pending.get();
            if (!result) {
                FAIL_CHECK("Missing independent worker result");
                continue;
            }
            const auto& actual = result.value();
            CHECK(actual == expected);
        }
    }
}
