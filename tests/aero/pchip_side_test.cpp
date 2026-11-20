#include "../support/allocation_audit.hpp"
#include "aero/pchip_boundary.hpp"
#include "aero/pchip_domain_route.hpp"
#include "aero/pchip_knot_locator.hpp"
#include "aero/pchip_side.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <future>
#include <iterator>
#include <vector>

namespace {
namespace c = uvlas::core;
namespace p = uvlas::aero::detail::pchip;
namespace audit = uvlas::test_support;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
using Side = p::RequestedSide;
using Order = p::JetOrder;
using Choice = p::JetCandidate;
using Location = p::SideLocation;

constexpr p::RationalLimits limits{4096, 10000, UINT64_MAX, 10000};
constexpr std::array sides{Side::Automatic, Side::Left, Side::Right};
constexpr std::array orders{Order::Value, Order::First, Order::Second};

template <class T>
T take(c::Result<T> result) {
    if (!result)
        FAIL("Unexpected failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}

struct Input {
    int sign;
    std::vector<std::uint32_t> n, d;

    explicit Input(const Rat& q) : sign(q < 0 ? -1 : q == 0 ? 0 : 1) {
        const Big magnitude = boost::multiprecision::abs(boost::multiprecision::numerator(q));
        if (sign != 0)
            boost::multiprecision::export_bits(magnitude, std::back_inserter(n), 32, false);
        boost::multiprecision::export_bits(boost::multiprecision::denominator(q),
                                           std::back_inserter(d), 32, false);
    }
    [[nodiscard]] p::RationalView view() const& noexcept {
        return {sign, n, d};
    }
    [[nodiscard]] p::RationalView view() const&& = delete;
};

Rat decode(p::RationalView value) {
    Big n = 0, d = 0;
    boost::multiprecision::import_bits(d, value.denominator.begin(), value.denominator.end(), 32,
                                       false);
    if (!value.numerator.empty())
        boost::multiprecision::import_bits(n, value.numerator.begin(), value.numerator.end(), 32,
                                           false);
    if (value.sign < 0)
        n = -n;
    return Rat{n} / d;
}

std::vector<Input> inputs(const std::vector<Rat>& values) {
    std::vector<Input> result;
    result.reserve(values.size());
    for (const auto& value : values)
        result.emplace_back(value);
    return result;
}

std::vector<p::RationalView> views(const std::vector<Input>& values) {
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

p::IntervalValues evaluate(const p::CoefficientBlock& block, std::uint64_t interval,
                           p::RationalView time, Order order) {
    c::ResourcePlan memory{UINT64_MAX};
    p::IntervalWorkspace workspace;
    return take(p::evaluate_interval(block, {interval, time, order}, workspace, limits, memory));
}

p::DomainRoute route(const p::CoefficientBlock& block, p::RationalView time, p::DomainPolicy mode) {
    c::ResourcePlan memory{UINT64_MAX};
    p::DomainWorkspace workspace;
    return take(p::route_domain(block, time, mode, workspace, limits, memory));
}

p::KnotLocation locate(const p::CoefficientBlock& block, p::RationalView time) {
    c::ResourcePlan memory{UINT64_MAX};
    p::LocatorWorkspace workspace;
    return take(p::locate_time(block, time, workspace, limits, memory));
}

struct End {
    Rat t, y, slope;
};

// Independent endpoint-basis Hermite values, with physical-time derivatives.
std::array<Rat, 3> hermite(const End& a, const End& b, const Rat& time) {
    const Rat h = b.t - a.t;
    const Rat u = (time - a.t) / h;
    const Rat u2 = u * u;
    const Rat u3 = u2 * u;
    return {(2 * u3 - 3 * u2 + 1) * a.y + (u3 - 2 * u2 + u) * h * a.slope +
                (-2 * u3 + 3 * u2) * b.y + (u3 - u2) * h * b.slope,
            ((6 * u2 - 6 * u) * a.y + (3 * u2 - 4 * u + 1) * h * a.slope + (-6 * u2 + 6 * u) * b.y +
             (3 * u2 - 2 * u) * h * b.slope) /
                h,
            ((12 * u - 6) * a.y + (6 * u - 4) * h * a.slope + (6 - 12 * u) * b.y +
             (6 * u - 2) * h * b.slope) /
                (h * h)};
}

void compare_jet(const p::IntervalValues& actual, std::uint64_t component, const End& a,
                 const End& b, const Rat& time, Order requested) {
    const auto exact = hermite(a, b, time);
    CHECK(decode(take(actual.time())) == time);
    for (std::size_t d = 0; d < orders.size(); ++d) {
        const auto value = actual.derivative(component, orders[d]);
        if (orders[d] <= requested) {
            if (value) {
                CHECK(decode(*value) == exact[d]);
            } else {
                FAIL("A requested derivative is unavailable");
            }
        } else {
            REQUIRE_FALSE(value);
            CHECK(value.error().code == c::FailureCode::InvalidInput);
        }
    }
}

void expect_match(const c::Result<p::SideAgreement>& result, bool equal, Choice common) {
    if (equal) {
        REQUIRE(result);
        CHECK(result->selected == common);
    } else {
        REQUIRE_FALSE(result);
        CHECK(result.error().code == c::FailureCode::InvalidInput);
        CHECK(result.error().operation == c::Operation::EvaluateAero);
    }
}

void check_join(const p::CoefficientBlock& source, std::uint64_t left, std::uint64_t right,
                p::RationalView left_time, p::RationalView right_time, const End& la, const End& lb,
                const End& ra, const End& rb) {
    const auto lt = decode(left_time), rt = decode(right_time);
    for (const auto side : sides) {
        for (const auto order : orders) {
            const auto decision = take(p::decide_side(Location::C1Join, side, order));
            if (const auto* direct = std::get_if<Choice>(&decision); direct != nullptr) {
                const bool choose_left = side == Side::Left;
                CHECK(*direct == (choose_left ? Choice::Left : Choice::Right));
                CHECK_FALSE((side == Side::Automatic && order == Order::Second));
                const auto jet = evaluate(source, choose_left ? left : right,
                                          choose_left ? left_time : right_time, order);
                compare_jet(jet, 0, choose_left ? la : ra, choose_left ? lb : rb,
                            choose_left ? lt : rt, order);
            } else {
                const auto* compare = std::get_if<p::CompareSecondDerivatives>(&decision);
                if (compare == nullptr) {
                    FAIL("Expected derivative agreement");
                    continue;
                }
                REQUIRE(side == Side::Automatic);
                REQUIRE(order == Order::Second);
                const auto l = evaluate(source, left, left_time, order);
                const auto r = evaluate(source, right, right_time, order);
                compare_jet(l, 0, la, lb, lt, order);
                compare_jet(r, 0, ra, rb, rt, order);
                const auto lhs = hermite(la, lb, lt), rhs = hermite(ra, rb, rt);
                REQUIRE(lhs[0] == rhs[0]);
                REQUIRE(lhs[1] == rhs[1]);
                expect_match(p::agree_side(*compare, l, r, UINT64_MAX), lhs[2] == rhs[2],
                             Choice::Right);
            }
        }
    }
}
} // namespace

TEST_CASE("PCHIP selects a single candidate throughout Reject and exterior Hold routes",
          "[AE-07][N4-06][V-AERO][verification]") {
    const auto source = prepare({0, 2}, {1, 3});
    const End a{0, 1, 1}, b{2, 3, 1};
    for (const Rat& time : std::array<Rat, 3>{0, 1, 2}) {
        const Input q{time};
        const auto domain = route(source, q.view(), p::DomainPolicy::Reject);
        REQUIRE(domain.has_inside());
        CHECK_FALSE(domain.has_constant());
        const auto interval = take(locate(source, q.view()).unique_interval());
        for (const auto side : sides)
            for (const auto order : orders) {
                const auto decision = take(p::decide_side(Location::Single, side, order));
                const auto* direct = std::get_if<Choice>(&decision);
                if (direct != nullptr) {
                    CHECK(*direct == Choice::Only);
                } else {
                    FAIL("Expected the single available candidate");
                }
                compare_jet(evaluate(source, interval, q.view(), order), 0, a, b, time, order);
            }
    }
    const Rat nearby = Rat{1} / (Big{1} << 200);
    const std::array<Rat, 2> exterior{-nearby, 2 + nearby};
    for (const auto& time : exterior) {
        const Input q{time};
        const auto domain = route(source, q.view(), p::DomainPolicy::Hold);
        REQUIRE(domain.has_constant());
        CHECK_FALSE(domain.has_inside());
        const auto node = take(domain.constant_node());
        for (const auto side : sides)
            for (const auto order : orders) {
                const auto decision = take(p::decide_side(Location::Single, side, order));
                const auto* direct = std::get_if<Choice>(&decision);
                if (direct != nullptr) {
                    CHECK(*direct == Choice::Only);
                } else {
                    FAIL("Expected the single available candidate");
                }
                CHECK(decode(take(source.value(node, 0))) == (time < 0 ? 1 : 3));
            }
    }
}

TEST_CASE("PCHIP knot and periodic seam Automatic compare second derivatives exactly",
          "[AE-07][N4-03][N4-06][V-AERO][verification]") {
    const Rat small = Rat{1} / (Big{1} << 200);
    for (const Rat& perturbation : std::array<Rat, 2>{0, small}) {
        const auto source = prepare({-1, 0, 1}, {1, 0, 1 + perturbation});
        const Input at_knot{0};
        const auto domain = route(source, at_knot.view(), p::DomainPolicy::Reject);
        REQUIRE(domain.has_inside());
        const auto knot = locate(source, at_knot.view());
        REQUIRE(knot.is_knot());
        const End a{-1, 1, -2 - perturbation / 2}, b{0, 0, 0},
            d{1, 1 + perturbation, 2 + 3 * perturbation / 2};
        CHECK(hermite(a, b, 0)[2] == 2 - perturbation);
        CHECK(hermite(b, d, 0)[2] == 2 + 3 * perturbation);
        check_join(source, take(knot.left_interval()), take(knot.right_interval()), at_knot.view(),
                   at_knot.view(), a, b, b, d);
    }

    const Input query{-6}, start{0}, period{3};
    for (const bool constant : {false, true}) {
        const auto source = prepare({0, 1, 2, 3}, constant ? std::vector<Rat>{0, 0, 0, 0}
                                                           : std::vector<Rat>{0, 1, -1, 0});
        c::ResourcePlan boundary_memory{UINT64_MAX};
        p::BoundaryWorkspace boundary_workspace;
        REQUIRE(p::check_boundary(source,
                                  {start.view(), period.view(), p::PeriodicBoundary{period.view()}},
                                  boundary_workspace, limits, boundary_memory));
        c::ResourcePlan time_memory{UINT64_MAX};
        p::DomainWorkspace time_workspace;
        const auto normalized = take(p::normalize_periodic_time(
            source, query.view(), period.view(), time_workspace, limits, time_memory));
        REQUIRE(normalized.is_seam());
        const auto right_time = take(normalized.time()), left_time = take(source.time(3));
        CHECK(decode(right_time) == 0);
        CHECK(decode(left_time) == 3);
        const auto a = End{0, 0, constant ? Rat{0} : Rat{5} / 2};
        const auto b = End{1, constant ? 0 : 1, 0};
        const auto d = End{2, constant ? 0 : -1, 0};
        const auto e = End{3, 0, a.slope};
        CHECK(hermite(d, e, 3)[2] == (constant ? 0 : 4));
        CHECK(hermite(a, b, 0)[2] == (constant ? 0 : -4));
        check_join(source, take(locate(source, left_time).left_interval()),
                   take(locate(source, right_time).right_interval()), left_time, right_time, d, e,
                   a, b);
    }
}

TEST_CASE("PCHIP Hold endpoint compares every requested derivative against zero",
          "[AE-07][N4-03][N4-06][V-AERO][verification]") {
    const Input first_time{0}, last_time{2};
    const std::array<std::vector<Rat>, 4> inputs{
        std::vector<Rat>{0, 1, 4}, {4, 1, 0}, {0, 1, 2}, {2, 2, 2}};
    // Endpoint slopes and exact interior slopes are known from the sample secants.
    const std::array<std::array<Rat, 3>, 4> slopes{
        {{0, Rat{3} / 2, 4}, {-4, Rat{-3} / 2, 0}, {1, 1, 1}, {0, 0, 0}}};
    for (std::size_t series = 0; series < inputs.size(); ++series) {
        CAPTURE(series);
        const auto& samples = inputs[series];
        const auto source = prepare({0, 1, 2}, samples);
        for (const bool at_first : {true, false}) {
            CAPTURE(at_first);
            const auto time = at_first ? first_time.view() : last_time.view();
            const auto domain = route(source, time, p::DomainPolicy::Hold);
            REQUIRE(domain.has_inside());
            REQUIRE(domain.has_constant());
            const auto knot = locate(source, time);
            const auto interval = take(knot.unique_interval());
            const auto condition = at_first ? Location::HoldAtFirst : Location::HoldAtLast;
            const auto inside = at_first ? Choice::Right : Choice::Left;
            const std::size_t a = at_first ? 0U : 1U;
            const End begin{a, samples[a], slopes[series][a]};
            const End end{a + 1, samples[a + 1], slopes[series][a + 1]};
            const Rat exact_time = at_first ? 0 : 2;
            const auto reference = hermite(begin, end, exact_time);
            REQUIRE(decode(take(source.value(take(domain.constant_node()), 0))) == reference[0]);
            for (const auto side : sides) {
                for (const auto order : orders) {
                    const auto decision = take(p::decide_side(condition, side, order));
                    if (const auto* direct = std::get_if<Choice>(&decision); direct != nullptr) {
                        const auto expected = side == Side::Automatic ? inside
                                              : side == Side::Left    ? Choice::Left
                                                                      : Choice::Right;
                        CHECK(*direct == expected);
                        CHECK_FALSE((side == Side::Automatic && order != Order::Value));
                        if (*direct == inside)
                            compare_jet(evaluate(source, interval, time, order), 0, begin, end,
                                        exact_time, order);
                    } else {
                        const auto* compare = std::get_if<p::CompareHeldDerivatives>(&decision);
                        if (compare == nullptr) {
                            FAIL("Expected derivative agreement");
                            continue;
                        }
                        CHECK(compare->through == order);
                        CHECK(compare->inside == inside);
                        REQUIRE(side == Side::Automatic);
                        const auto jet = evaluate(source, interval, time, order);
                        compare_jet(jet, 0, begin, end, exact_time, order);
                        const bool zero =
                            reference[1] == 0 && (order == Order::First || reference[2] == 0);
                        expect_match(p::agree_side(*compare, jet, UINT64_MAX), zero, inside);
                    }
                }
            }
        }
    }
}

TEST_CASE("PCHIP side agreement bounds work and compares more than the first component",
          "[AE-07][N4-02][N4-08][N4-09][negative]") {
    const auto continuous = prepare({-1, 0, 1}, {1, 8, 0, 6, 1, 8}, 2);
    const Rat small = Rat{1} / (Big{1} << 200);
    const auto unequal = prepare({-1, 0, 1}, {1, 8, 0, 6, 1, 8 + small}, 2);
    const auto scalar = prepare({-1, 0, 1}, {1, 0, 1});
    const Input t{0};
    const auto left = evaluate(continuous, 0, t.view(), Order::Second);
    const auto right = evaluate(continuous, 1, t.view(), Order::Second);
    const auto bad_left = evaluate(unequal, 0, t.view(), Order::Second);
    const auto bad_right = evaluate(unequal, 1, t.view(), Order::Second);
    const auto first_only = evaluate(continuous, 0, t.view(), Order::First);
    const auto scalar_jet = evaluate(scalar, 1, t.view(), Order::Second);
    const p::CompareSecondDerivatives comparison;
    expect_match(p::agree_side(comparison, bad_left, bad_right, UINT64_MAX), false, Choice::Right);
    expect_match(p::agree_side(comparison, left, first_only, UINT64_MAX), false, Choice::Right);
    expect_match(p::agree_side(comparison, left, scalar_jet, UINT64_MAX), false, Choice::Right);
    CHECK(decode(take(bad_left.derivative(0, Order::Second))) ==
          decode(take(bad_right.derivative(0, Order::Second))));

    const auto baseline = take(p::agree_side(comparison, left, right, UINT64_MAX));
    CHECK(baseline.selected == Choice::Right);
    audit::AllocationAudit observer{audit::cpp_backend, audit::AuditThreadDomain::CallingThread};
    const auto equal = p::agree_side(comparison, left, right, baseline.used_steps);
    const auto short_budget = p::agree_side(comparison, left, right, baseline.used_steps - 1);
    const auto observed = observer.finish();
    CHECK(take(equal).used_steps == baseline.used_steps);
    REQUIRE_FALSE(short_budget);
    CHECK(short_budget.error().code == c::FailureCode::ResourceExceeded);
    CHECK(observed.status == audit::AllocationAuditStatus::Clean);
    CHECK(observed.cpp_allocation_attempts == 0);

    for (const std::size_t count : {1UZ, 2UZ, 4UZ}) {
        std::vector<std::future<c::Result<p::SideAgreement>>> jobs;
        jobs.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            jobs.push_back(std::async(std::launch::async, [&] {
                return p::agree_side(comparison, left, right, baseline.used_steps);
            }));
        }
        for (auto& job : jobs) {
            const auto result = take(job.get());
            CHECK(result.selected == baseline.selected);
            CHECK(result.used_steps == baseline.used_steps);
        }
    }

    const auto constant = prepare({0, 2}, {1, 4, 1, 4}, 2);
    const auto varying = prepare({0, 2}, {1, 4, 1, 5}, 2);
    const auto flat = evaluate(constant, 0, t.view(), Order::Second);
    const auto not_flat = evaluate(varying, 0, t.view(), Order::Second);
    const p::CompareHeldDerivatives hold{Choice::Right, Order::Second};
    expect_match(p::agree_side(hold, not_flat, UINT64_MAX), false, Choice::Right);
    const auto hold_result = take(p::agree_side(hold, flat, UINT64_MAX));
    CHECK(hold_result.used_steps == 5);
    CHECK(take(p::agree_side(hold, flat, hold_result.used_steps)).selected == Choice::Right);
    const auto hold_short = p::agree_side(hold, flat, hold_result.used_steps - 1);
    REQUIRE_FALSE(hold_short);
    CHECK(hold_short.error().code == c::FailureCode::ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() is non-const.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    c::FailureDetails details;
    details.model = c::ModelGeneration{225};
    const auto stopped =
        p::agree_side(comparison, left, right, UINT64_MAX, stop.get_token(), details);
    REQUIRE_FALSE(stopped);
    CHECK(stopped.error().code == c::FailureCode::Cancelled);
    CHECK(stopped.error().details == details);
    const auto stopped_hold = p::agree_side(hold, flat, UINT64_MAX, stop.get_token(), details);
    REQUIRE_FALSE(stopped_hold);
    CHECK(stopped_hold.error().code == c::FailureCode::Cancelled);
    expect_match(p::agree_side(comparison, left, first_only, UINT64_MAX, stop.get_token()), false,
                 Choice::Right);
    CHECK(decode(take(left.derivative(0, Order::Second))) == 2);
    CHECK(decode(take(flat.derivative(1, Order::Second))) == 0);
}

TEST_CASE("PCHIP rejects unknown side decisions before establishing a candidate",
          "[AE-07][N4-08][negative]") {
    const auto unknown_side =
        p::decide_side(Location::Single, static_cast<Side>(255), Order::Value);
    REQUIRE_FALSE(unknown_side);
    CHECK(unknown_side.error().code == c::FailureCode::InvalidInput);
    const auto unknown_order =
        p::decide_side(Location::Single, Side::Left, static_cast<Order>(255));
    REQUIRE_FALSE(unknown_order);
    CHECK(unknown_order.error().code == c::FailureCode::InvalidInput);
    const auto unknown_location =
        p::decide_side(static_cast<Location>(255), Side::Right, Order::Second);
    REQUIRE_FALSE(unknown_location);
    CHECK(unknown_location.error().code == c::FailureCode::InvalidInput);
}
