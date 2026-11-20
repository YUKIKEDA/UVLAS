#include "aero/pchip_boundary.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <future>
#include <iterator>
#include <limits>
#include <vector>

namespace {
namespace p = uvlas::aero::detail::pchip;
namespace c = uvlas::core;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
constexpr p::RationalLimits limits{4096, 10000, std::numeric_limits<std::uint64_t>::max(), 10000};

template <class T>
T take(c::Result<T> result) {
    if (!result)
        FAIL("Unexpected failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}

struct ExactInput {
    int sign;
    std::vector<std::uint32_t> numerator, denominator;

    explicit ExactInput(const Rat& value) : sign(value < 0 ? -1 : value == 0 ? 0 : 1) {
        Big n = boost::multiprecision::numerator(value);
        if (sign < 0)
            n = -n;
        if (sign != 0)
            boost::multiprecision::export_bits(n, std::back_inserter(numerator), 32, false);
        boost::multiprecision::export_bits(boost::multiprecision::denominator(value),
                                           std::back_inserter(denominator), 32, false);
    }

    [[nodiscard]] p::RationalView view() const& noexcept {
        return {sign, numerator, denominator};
    }
    [[nodiscard]] p::RationalView view() const&& = delete;
};

std::vector<ExactInput> inputs(const std::vector<Rat>& values) {
    std::vector<ExactInput> result;
    result.reserve(values.size());
    for (const auto& value : values)
        result.emplace_back(value);
    return result;
}

std::vector<p::RationalView> views(const std::vector<ExactInput>& values) {
    std::vector<p::RationalView> result;
    result.reserve(values.size());
    for (const auto& value : values)
        result.push_back(value.view());
    return result;
}

p::CoefficientBlock prepare(const std::vector<Rat>& times, const std::vector<Rat>& values,
                            std::uint64_t components = 1) {
    const auto t = inputs(times), y = inputs(values);
    const auto tv = views(t), yv = views(y);
    c::ResourcePlan memory{UINT64_MAX};
    p::CoefficientWorkspace workspace;
    return take(p::prepare_coefficients({tv, yv, components}, workspace, limits, memory));
}

using Policy = decltype(p::BoundaryRequest::policy);

c::Result<p::BoundaryUsage> check(const p::CoefficientBlock& source, const Rat& begin,
                                  const Rat& end, Policy policy) {
    const ExactInput a{begin}, b{end};
    const p::BoundaryRequest request{a.view(), b.view(), policy};
    const auto plan = take(p::forecast_boundary(source, request, limits));
    c::ResourcePlan memory{plan.payload_bytes()};
    p::BoundaryWorkspace workspace;
    const auto result = p::check_boundary(source, request, workspace, limits, memory);
    if (result) {
        CHECK(result->steps <= plan.reserved_steps());
        CHECK(result->slots <= plan.arithmetic().slots());
    }
    return result;
}

void matches(const c::Result<p::BoundaryUsage>& result, bool accepted) {
    if (accepted) {
        REQUIRE(result);
    } else {
        REQUIRE_FALSE(result);
        CHECK(result.error().code == c::FailureCode::InvalidInput);
        CHECK(result.error().operation == c::Operation::EvaluateAero);
    }
}

// Independent physical slope from the endpoint secants, not from the u-coefficients.
Rat endpoint(const Rat& h, const Rat& next_h, const Rat& d, const Rat& next_d) {
    Rat candidate = d + h * (d - next_d) / (h + next_h);
    const auto sign = [](const Rat& value) { return value < 0 ? -1 : value == 0 ? 0 : 1; };
    if (sign(candidate) != sign(d))
        return 0;
    if (sign(d) != sign(next_d) &&
        boost::multiprecision::abs(candidate) > 3 * boost::multiprecision::abs(d))
        return 3 * d;
    return candidate;
}
} // namespace

TEST_CASE("PCHIP boundary qualification uses the inclusive demanded interval",
          "[AE-07][N4-02][N4-03][N4-06][V-AERO][verification]") {
    const auto linear = prepare({0, 2}, {0, 2});
    const auto constant = prepare({0, 2}, {1, 1});
    const auto left_flat = prepare({0, 1, 2}, {0, 0, 1});
    const auto right_flat = prepare({0, 1, 2}, {1, 0, 0});
    const auto two = prepare({0, 1, 2}, {1, 0, 0, 0, 0, 1}, 2);
    const Rat low_bit = Rat{1} / (Big{1} << 200);
    const std::array<Rat, 7> points{-low_bit, 0, low_bit, 1, 2 - low_bit, 2, 2 + low_bit};
    for (const auto& begin : points) {
        for (const auto& end : points) {
            CAPTURE(begin, end);
            const bool ordered = begin <= end;
            const bool first_contact = begin <= 0 && end >= 0;
            const bool last_contact = begin <= 2 && end >= 2;
            matches(check(linear, begin, end, p::RejectBoundary{}),
                    ordered && begin >= 0 && end <= 2);
            matches(check(constant, begin, end, p::HoldBoundary{}), ordered);
            matches(check(linear, begin, end, p::HoldBoundary{}),
                    ordered && !first_contact && !last_contact);
            matches(check(left_flat, begin, end, p::HoldBoundary{}), ordered && !last_contact);
            matches(check(right_flat, begin, end, p::HoldBoundary{}), ordered && !first_contact);
            matches(check(two, begin, end, p::HoldBoundary{}),
                    ordered && !first_contact && !last_contact);
        }
    }
}

TEST_CASE("PCHIP periodic qualification compares all seam values and physical derivatives",
          "[AE-07][N4-03][N4-06][N4-08][V-AERO][verification]") {
    // Left width 1, right width 2. The first valid case has physical slopes (2,2)
    // and u derivatives (2,4). The next has physical slopes (2,3/2).
    const std::vector<Rat> times{0, 1, 2, 3, 5};
    const Rat low_bit = Rat{1} / (Big{1} << 200);
    const std::vector<std::vector<Rat>> series{
        {0, 1, 0, Rat{-4} / 3, 0},           {0, 1, 0, -1, 0}, {0, 1, 0, Rat{-4} / 3, low_bit},
        {0, 1, 0, Rat{-4} / 3 + low_bit, 0}, {0, 0, 1, 0, 0},  {0, 0, 0, 0, 0}};
    const ExactInput period{5};
    const p::PeriodicBoundary policy{period.view()};
    for (std::size_t i = 0; i < series.size(); ++i) {
        CAPTURE(i);
        const auto& y = series[i];
        const Rat left = endpoint(1, 1, y[1] - y[0], y[2] - y[1]);
        const Rat right = endpoint(2, 1, (y[4] - y[3]) / 2, y[3] - y[2]);
        if (i == 0) {
            CHECK(left == 2);
            CHECK(right == 2);
        }
        if (i == 1)
            CHECK(right == Rat{3} / 2);
        const bool qualifies = y.front() == y.back() && left == right;
        matches(check(prepare(times, y), -17, 28, policy), qualifies);

        std::vector<Rat> vector_values;
        for (std::size_t row = 0; row < times.size(); ++row) {
            vector_values.push_back(series.front()[row]);
            vector_values.push_back(y[row]);
        }
        matches(check(prepare(times, vector_values, 2), -17, 28, policy), qualifies);
    }
    const auto source = prepare(times, series.front());
    for (const Rat& wrong : std::vector<Rat>{0, -5, Rat{5} - low_bit, Rat{5} + low_bit}) {
        const ExactInput supplied{wrong};
        matches(check(source, 1, 1, p::PeriodicBoundary{supplied.view()}), false);
    }
    matches(check(source, 1, 0, policy), false);
    const ExactInput start{0}, finish{5};
    const auto missing = p::forecast_boundary(
        source, {start.view(), finish.view(), p::PeriodicBoundary{{}}}, limits);
    REQUIRE_FALSE(missing);
    CHECK(missing.error().code == c::FailureCode::InvalidInput);

    const std::array<std::uint32_t, 1> ten{10}, two{2};
    matches(check(source, 1, 1, p::PeriodicBoundary{{1, ten, two}}), false);
}

TEST_CASE("PCHIP boundary work is deterministic across independent workspaces",
          "[AE-07][N4-03][N4-09][property]") {
    const auto source = prepare({0, 1, 2, 3, 5}, {0, 1, 0, Rat{-4} / 3, 0});
    const ExactInput begin{-7}, end{8}, period{5};
    const p::BoundaryRequest request{begin.view(), end.view(), p::PeriodicBoundary{period.view()}};
    const auto plan = take(p::forecast_boundary(source, request, limits));
    const auto run = [&] {
        c::ResourcePlan memory{plan.payload_bytes()};
        p::BoundaryWorkspace workspace;
        return p::check_boundary(source, request, workspace, limits, memory);
    };
    const auto baseline = take(run());
    for (const std::size_t count : {1UZ, 2UZ, 4UZ}) {
        std::vector<std::future<c::Result<p::BoundaryUsage>>> results;
        results.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            results.push_back(std::async(std::launch::async, run));
        for (auto& pending : results) {
            const auto result = take(pending.get());
            CHECK(result.steps == baseline.steps);
            CHECK(result.slots == baseline.slots);
        }
    }
}
