#include "aero/pchip_domain_route.hpp"
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
TEST_CASE("Private PCHIP domain matches independent rational classification",
          "[AE-07][N4-02][N4-03][N4-06][verification]") {
    std::ifstream stream{UVLAS_PCHIP_COEFFICIENT_REFERENCE};
    REQUIRE(stream.good());
    const auto reference = nlohmann::json::parse(stream);
    for (const auto& record : reference.at("cases")) {
        std::vector<Rat> times, values;
        for (const auto& v : record.at("times"))
            times.emplace_back(v.get<std::string>());
        for (const auto& v : record.at("values"))
            values.emplace_back(v.get<std::string>());
        const auto source =
            coefficients(times, values, record.at("components").get<std::uint64_t>());
        const Rat small = Rat{1} / (Big{1} << 400);
        std::vector<Rat> queries{times.front() - small, times.front(), times.front() + small,
                                 times.back() - small,  times.back(),  times.back() + small};
        for (const auto& t : times)
            queries.push_back(t);
        for (const auto& q : queries)
            for (const auto policy : {p::DomainPolicy::Reject, p::DomainPolicy::Hold}) {
                const Input input{q};
                const bool outside = q < times.front() || q > times.back();
                const bool endpoint = q == times.front() || q == times.back();
                const auto plan = take(p::forecast_domain(source, input.view(), policy, limits()));
                core::ResourcePlan memory{plan.payload_bytes()};
                p::DomainWorkspace workspace;
                auto result =
                    p::route_domain(source, input.view(), policy, workspace, limits(), memory);
                if (policy == p::DomainPolicy::Reject && outside) {
                    REQUIRE_FALSE(result);
                    CHECK(result.error().code == core::FailureCode::InvalidInput);
                    continue;
                }
                const auto route = take(std::move(result));
                const auto expected = q < times.front()    ? p::DomainRelation::Before
                                      : q == times.front() ? p::DomainRelation::AtFirst
                                      : q > times.back()   ? p::DomainRelation::After
                                      : q == times.back()  ? p::DomainRelation::AtLast
                                                           : p::DomainRelation::Inside;
                CHECK(route.relation() == expected);
                CHECK(route.has_inside() == !outside);
                CHECK(route.has_constant() ==
                      (policy == p::DomainPolicy::Hold && (outside || endpoint)));
                CHECK(route.used_steps() <= plan.reserved_steps());
                CHECK(route.used_slots() == 3);
                if (route.has_constant()) {
                    const auto node = take(route.constant_node());
                    CHECK(node == (q <= times.front() ? std::uint64_t{0} : source.nodes() - 1));
                    for (std::uint64_t c = 0; c < source.components(); ++c)
                        CHECK(decoded(take(source.value(node, c))) ==
                              values[static_cast<std::size_t>(node * source.components() + c)]);
                } else
                    CHECK_FALSE(route.constant_node());
                if (route.has_inside()) {
                    const auto lp = take(p::forecast_location(source, input.view(), limits()));
                    core::ResourcePlan lm{lp.payload_bytes()};
                    p::LocatorWorkspace lw;
                    const auto location =
                        take(p::locate_time(source, input.view(), lw, limits(), lm));
                    const auto interval = location.has_right() ? take(location.right_interval())
                                                               : take(location.left_interval());
                    const p::IntervalRequest request{interval, input.view(), p::JetOrder::Value};
                    const auto ep = take(p::forecast_interval(source, request, limits()));
                    core::ResourcePlan em{ep.payload_bytes()};
                    p::IntervalWorkspace ew;
                    const auto evaluated =
                        take(p::evaluate_interval(source, request, ew, limits(), em));
                    CHECK(decoded(take(evaluated.time())) == q);
                    if (endpoint)
                        for (std::uint64_t c = 0; c < source.components(); ++c) {
                            const auto node =
                                q == times.front() ? std::uint64_t{0} : source.nodes() - 1;
                            CHECK(decoded(take(evaluated.derivative(c, p::JetOrder::Value))) ==
                                  decoded(take(source.value(node, c))));
                        }
                }
            }
    }
}
TEST_CASE("Private PCHIP domain rejects invalid inputs and preserves scalar",
          "[AE-07][N4-08][N4-09][negative]") {
    auto source = coefficients({0, 1, 2}, {0, 1, 4});
    const Input input{Rat{1} / 2};
    const auto unknown =
        p::forecast_domain(source, input.view(), static_cast<p::DomainPolicy>(255), limits());
    REQUIRE_FALSE(unknown);
    CHECK(unknown.error().code == core::FailureCode::Unsupported);
    const std::uint32_t two = 2;
    const p::RationalView invalid{1, {&two, 1}, {&two, 1}};
    const auto plan = take(p::forecast_domain(source, invalid, p::DomainPolicy::Hold, limits()));
    core::ResourcePlan memory{plan.payload_bytes()};
    p::DomainWorkspace workspace;
    const auto rejected =
        p::route_domain(source, invalid, p::DomainPolicy::Hold, workspace, limits(), memory);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    const Input endpoint{0};
    const auto valid =
        take(p::forecast_domain(source, endpoint.view(), p::DomainPolicy::Hold, limits()));
    core::ResourcePlan good{valid.payload_bytes()};
    const auto route = take(
        p::route_domain(source, endpoint.view(), p::DomainPolicy::Hold, workspace, limits(), good));
    auto moved = std::move(source);
    CHECK(route.has_inside());
    CHECK(route.has_constant());
    CHECK(take(route.constant_node()) == 0);
    CHECK(decoded(take(moved.value(0, 0))) == 0);
    // NOLINTNEXTLINE(bugprone-use-after-move): the documented inactive owner is rejected.
    const auto inactive = p::forecast_domain(source, input.view(), p::DomainPolicy::Hold, limits());
    REQUIRE_FALSE(inactive);
    CHECK(inactive.error().code == core::FailureCode::InvalidInput);
}
