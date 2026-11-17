#include "aero/exact_triangle_predicates.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <stop_token>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using namespace uvlas;
namespace subject = aero::detail;
using Rat = boost::multiprecision::cpp_rational;
using DynamicInteger = boost::multiprecision::cpp_int;
using RV = std::array<Rat, 3>;
using Words = std::array<double, 12>;
constexpr std::uint64_t owned_budget = 4ULL * 1024 * 1024;
const subject::PredicateLimits limits{numerics::linear::EstimateAdmission{28672}, 100};

template <class T>
T take(core::Result<T> value) {
    if (!value) {
        FAIL_CHECK("Unexpected failure " << static_cast<unsigned>(value.error().code));
        throw std::runtime_error("Unexpected predicate failure");
    }
    return std::move(*value);
}
Rat canonical(const Rat& value) {
    return boost::multiprecision::numerator(value) == 0 ? Rat{0} : value;
}
int sign(const Rat& value) {
    const auto numerator = boost::multiprecision::numerator(value);
    return numerator == 0 ? 0 : numerator < 0 ? -1 : 1;
}
RV minus(const RV& a, const RV& b) {
    return {canonical(a[0] - b[0]), canonical(a[1] - b[1]), canonical(a[2] - b[2])};
}
RV plus(const RV& a, const RV& b) {
    return {canonical(a[0] + b[0]), canonical(a[1] + b[1]), canonical(a[2] + b[2])};
}
RV times(const RV& a, const Rat& t) {
    return {canonical(a[0] * t), canonical(a[1] * t), canonical(a[2] * t)};
}
Rat inner(const RV& a, const RV& b) {
    return canonical(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}
RV vector_product(const RV& a, const RV& b) {
    return {canonical(a[1] * b[2] - a[2] * b[1]), canonical(a[2] * b[0] - a[0] * b[2]),
            canonical(a[0] * b[1] - a[1] * b[0])};
}
Rat divide(const Rat& a, const Rat& b) {
    return sign(a) == 0 ? Rat{0} : canonical(a / b);
}
Rat distance(const RV& a, const RV& b) {
    const auto d = minus(a, b);
    return inner(d, d);
}
bool equal_point(const RV& a, const RV& b) {
    return sign(distance(a, b)) == 0;
}
Rat physical(const subject::ExactInteger& integer, int exponent) {
    Rat value{DynamicInteger{integer}};
    if (sign(value) == 0)
        return Rat{0};
    if (exponent < 0)
        value /= DynamicInteger{1} << static_cast<unsigned>(-exponent);
    else
        value *= DynamicInteger{1} << static_cast<unsigned>(exponent);
    return value;
}
const subject::NondegenerateTriangle& nondegenerate(const subject::ExactTriangleResult& value) {
    const auto* outcome = value.outcome();
    if (outcome == nullptr) {
        FAIL_CHECK("Missing owned outcome");
        throw std::runtime_error("Moved outcome");
    }
    const auto* actual = std::get_if<subject::NondegenerateTriangle>(outcome);
    if (actual == nullptr) {
        FAIL_CHECK("Expected nondegenerate geometry");
        throw std::runtime_error("Unexpected degeneracy");
    }
    return *actual;
}
std::unique_ptr<subject::PredicateWorkspace> workspace() {
    core::ResourcePlan plan{owned_budget};
    return take(subject::PredicateWorkspace::create(plan));
}
subject::ExactTriangleResult evaluate(const Words& high, const Words* low = nullptr) {
    auto scratch = workspace();
    core::ResourcePlan result{owned_budget};
    subject::CoordinateView view = subject::Binary64Coordinates{high};
    if (low != nullptr)
        view = subject::TwoWordCoordinates{high, *low};
    return take(subject::classify_triangle(view, *scratch, result, limits));
}
Words fixture(double px, double py, double pz) {
    return {0, 0, 0, 3, 0, 0, 0, 4, 0, px, py, pz};
}
void compare_oracle(const Words& high, const Words* low = nullptr) {
    const auto actual = evaluate(high, low);
    std::array<RV, 4> p;
    for (std::size_t i = 0; i < 12; ++i)
        p[i / 3][i % 3] = canonical(Rat{high[i]} + (low == nullptr ? Rat{0} : Rat{(*low)[i]}));
    const auto normal = vector_product(minus(p[1], p[0]), minus(p[2], p[0]));
    const auto norm_squared = inner(normal, normal);
    if (sign(norm_squared) == 0) {
        REQUIRE(actual.outcome() != nullptr);
        CHECK(std::holds_alternative<subject::DegenerateTriangle>(*actual.outcome()));
        return;
    }
    const auto& got = nondegenerate(actual);
    const auto signed_height = inner(normal, minus(p[3], p[0]));
    CHECK(got.plane == (sign(signed_height) < 0   ? subject::PlaneSide::Below
                        : sign(signed_height) > 0 ? subject::PlaneSide::Above
                                                  : subject::PlaneSide::On));
    // Independently form the closest physical points. Do not reuse the product
    // Gram numerators or squared-distance candidate expressions.
    const auto projected = minus(p[3], times(normal, divide(signed_height, norm_squared)));
    bool inside = true;
    std::array<Rat, 3> edge_side;
    for (std::size_t i = 0; i < 3; ++i) {
        edge_side[i] =
            inner(vector_product(minus(p[(i + 1) % 3], p[i]), minus(projected, p[i])), normal);
        inside = inside && sign(edge_side[i]) >= 0;
    }
    auto projection = subject::ProjectionFeature::Outside;
    if (inside) {
        if (equal_point(projected, p[0]))
            projection = subject::ProjectionFeature::Vertex0;
        else if (equal_point(projected, p[1]))
            projection = subject::ProjectionFeature::Vertex1;
        else if (equal_point(projected, p[2]))
            projection = subject::ProjectionFeature::Vertex2;
        else if (sign(edge_side[0]) == 0)
            projection = subject::ProjectionFeature::Edge01;
        else if (sign(edge_side[1]) == 0)
            projection = subject::ProjectionFeature::Edge12;
        else if (sign(edge_side[2]) == 0)
            projection = subject::ProjectionFeature::Edge20;
        else
            projection = subject::ProjectionFeature::Face;
    }
    CHECK(got.projection == projection);
    std::vector<RV> candidates{p[0], p[1], p[2]};
    Rat longest = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto edge = minus(p[(i + 1) % 3], p[i]);
        const auto length = inner(edge, edge);
        if (length > longest)
            longest = length;
        auto t = divide(inner(minus(p[3], p[i]), edge), length);
        if (sign(t) < 0)
            t = 0;
        if (t > 1)
            t = 1;
        candidates.push_back(plus(p[i], times(edge, t)));
    }
    if (inside)
        candidates.push_back(projected);
    Rat minimum = distance(p[3], candidates[0]);
    RV closest = candidates[0];
    std::uint8_t mask = 0;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const auto d = distance(p[3], candidates[i]);
        if (d < minimum) {
            minimum = d;
            closest = candidates[i];
            mask = 0;
        }
        if (d == minimum)
            mask = static_cast<std::uint8_t>(mask | static_cast<std::uint8_t>(1U << i));
    }
    REQUIRE(got.distance_squared.denominator > 0);
    CHECK(got.distance_squared.unit_exponent == -4296);
    const Rat measured = divide(physical(got.distance_squared.numerator, -4296),
                                Rat{DynamicInteger{got.distance_squared.denominator}});
    // The minimum must be <= EVERY valid candidate, and attained by at least one.
    bool attained = false;
    for (const auto& candidate : candidates) {
        const auto d = distance(p[3], candidate);
        CHECK(measured <= d);
        attained = attained || measured == d;
    }
    CHECK(attained);
    CHECK(measured == distance(p[3], closest));
    CHECK(measured == minimum);
    CHECK(got.minimizing_features == mask);
    CHECK(physical(got.longest_edge_squared, -4296) == longest);
    CHECK(got.far == (minimum > 16 * longest    ? subject::FarRelation::Far
                      : minimum == 16 * longest ? subject::FarRelation::Boundary
                                                : subject::FarRelation::Near));
    if (sign(minimum) == 0)
        CHECK(got.distance_squared.denominator == 1);
}
void failure(const core::Result<subject::ExactTriangleResult>& value, core::FailureCode code) {
    REQUIRE_FALSE(value.has_value());
    CHECK(value.error().code == code);
    CHECK(value.error().operation == core::Operation::EvaluateAero);
}
} // namespace

TEST_CASE("Exact triangle literals retain finite face distance and every support",
          "[verification][UI-D01][UI-D03][UI-D06][N4-06]") {
    for (const auto point :
         std::array{std::array{.75, 1., 0.}, std::array{.75, 0., 0.}, std::array{1.5, 2., 0.},
                    std::array{0., 1., 0.}, std::array{0., 0., 0.}, std::array{3., 0., 0.},
                    std::array{0., 4., 0.}, std::array{-.5, -.5, 0.}, std::array{.75, 1., 2.},
                    std::array{.75, 1., -2.}})
        compare_oracle(fixture(point[0], point[1], point[2]));
    const auto vertex = evaluate(fixture(0, 0, 0));
    CHECK(nondegenerate(vertex).minimizing_features == 105); // 1 + 8 + 32 + 64.
    const Words exterior{0, 0, 0, 2, 0, 0, 0, 2, 0, 3, 3, 0};
    compare_oracle(exterior);
    const auto value = evaluate(exterior);
    const auto& result = nondegenerate(value);
    CHECK(divide(physical(result.distance_squared.numerator, -4296),
                 Rat{DynamicInteger{result.distance_squared.denominator}}) == 8);
    for (const double z : std::array{std::nextafter(20., 0.), 20., std::nextafter(20., 21.)})
        compare_oracle(fixture(.75, 1, z));
}
TEST_CASE("Exact triangle values survive all vertex orders and physical unit changes",
          "[property][verification][UI-D05][UI-D06][N4-06]") {
    const auto original = fixture(.75, 1, 2);
    std::array<std::size_t, 3> order{0, 1, 2};
    do {
        Words permuted = original;
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t k = 0; k < 3; ++k)
                permuted[3 * i + k] = original[3 * order[i] + k];
        compare_oracle(permuted);
    } while (std::next_permutation(order.begin(), order.end()));
    for (const double scale : std::array{1e-3, 1e3, 0x1p-900, 0x1p900}) {
        Words high{}, low{};
        for (std::size_t i = 0; i < 12; ++i) {
            high[i] = original[i] * scale;
            low[i] = std::fma(original[i], scale, -high[i]);
        }
        compare_oracle(high, &low);
    }
    auto rotated = original;
    for (std::size_t i = 0; i < 4; ++i) {
        rotated[3 * i + 1] = -original[3 * i + 1];
        rotated[3 * i + 2] = -original[3 * i + 2];
    }
    compare_oracle(rotated);
}
TEST_CASE("Every finite binary64 exponent and the retained low word remains exact",
          "[verification][negative][UI-D06][N4-06]") {
    const double tiny = std::numeric_limits<double>::denorm_min();
    compare_oracle(Words{0, 0, 0, tiny, 0, 0, 0, tiny, 0, tiny, tiny, tiny});
    const double large = std::numeric_limits<double>::max();
    compare_oracle(Words{large, large, large, -large, large, large, large, -large, large, large,
                         large, -large});
    compare_oracle(Words{0, 0, 0, large, 0, 0, 0, tiny, 0, 0, tiny, tiny});
    compare_oracle(Words{0, 0, 0, 1, 0, 0, 1, 0x1p-60, 0, .5, 0x1p-62, 1});
    compare_oracle(Words{0, 0, 0, 1, 0, 0, 1, 0x1p-60, 0, .5, 0x1p-62, 0});
    compare_oracle(Words{0, 0, 0, 1, 0, 0, 2, 0, 0, .5, 0, 1});
    Words high{};
    high.fill(0x1p500);
    const Words low{0, 0, 0, 3, 0, 0, 0, 4, 0, .75, 1, 0};
    compare_oracle(high, &low);
    compare_oracle(high);
    auto negative_zero = fixture(.75, 1, 0);
    for (auto& value : negative_zero)
        if (value == 0)
            value = -0.0;
    compare_oracle(negative_zero);
    auto sum_high = fixture(.75, 1, 0), sum_low = Words{};
    sum_high[0] = large;
    sum_low[0] = large;
    compare_oracle(sum_high, &sum_low);
}
TEST_CASE("Every binary64 power of two is classified without floating point zero thresholds",
          "[verification][UI-D06][N4-06]") {
    auto scratch = workspace();
    for (int exponent = -1074; exponent <= 1023; ++exponent) {
        const auto input = fixture(.75, 1, std::ldexp(1., exponent));
        core::ResourcePlan output{owned_budget};
        const auto result = take(subject::classify_triangle(subject::Binary64Coordinates{input},
                                                            *scratch, output, limits));
        const auto& got = nondegenerate(result);
        CHECK(got.plane == subject::PlaneSide::Above);
        Rat expected{1};
        if (exponent < 0)
            expected /= DynamicInteger{1} << static_cast<unsigned>(-2 * exponent);
        else
            expected *= DynamicInteger{1} << static_cast<unsigned>(2 * exponent);
        CHECK(divide(physical(got.distance_squared.numerator, -4296),
                     Rat{DynamicInteger{got.distance_squared.denominator}}) == expected);
    }
}
TEST_CASE("Predicate limits and failures preserve borrowed words and existing results",
          "[unit][negative][N4-02][N4-08][N4-09]") {
    auto scratch = workspace();
    const auto input = fixture(.75, 1, 2);
    core::ResourcePlan exact_plan{subject::ExactTriangleResult::payload_bytes()};
    auto kept = take(subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch,
                                                exact_plan, limits));
    const auto old = nondegenerate(kept).distance_squared.numerator;
    core::ResourcePlan none{0};
    failure(subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, none, limits),
            core::FailureCode::ResourceExceeded);
    core::ResourcePlan sufficient{owned_budget};
    failure(subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, sufficient,
                                       {numerics::linear::EstimateAdmission{28671}, 100}),
            core::FailureCode::ResourceExceeded);
    for (const auto cap :
         std::array<std::uint64_t, 3>{0, 28672, std::numeric_limits<std::uint64_t>::max()})
        failure(subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch,
                                           sufficient,
                                           {numerics::linear::RequireProvenBound{cap}, 100}),
                core::FailureCode::Unsupported);
    failure(subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, sufficient,
                                       {numerics::linear::EstimateAdmission{28672}, 0}),
            core::FailureCode::ResourceExceeded);
    // MSVC 14.51 request_stop() is non-const; libstdc++ declares it const.
    // NOLINTNEXTLINE(misc-const-correctness): keep the MSVC call well-formed.
    std::stop_source stop;
    stop.request_stop();
    failure(subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, sufficient,
                                       limits, stop.get_token()),
            core::FailureCode::Cancelled);
    for (const auto value : std::array{std::numeric_limits<double>::infinity(),
                                       -std::numeric_limits<double>::infinity(),
                                       std::numeric_limits<double>::quiet_NaN()}) {
        for (std::size_t i = 0; i < 12; ++i) {
            auto invalid = input;
            invalid[i] = value;
            failure(subject::classify_triangle(subject::Binary64Coordinates{invalid}, *scratch,
                                               sufficient, limits),
                    core::FailureCode::Nonfinite);
            failure(subject::classify_triangle(subject::TwoWordCoordinates{input, invalid},
                                               *scratch, sufficient, limits),
                    core::FailureCode::Nonfinite);
        }
    }
    failure(subject::classify_triangle(subject::Binary64Coordinates{std::span{input}.first(11)},
                                       *scratch, sufficient, limits),
            core::FailureCode::InvalidInput);
    failure(subject::classify_triangle(subject::TwoWordCoordinates{input, {}}, *scratch, sufficient,
                                       limits),
            core::FailureCode::InvalidInput);
    CHECK(nondegenerate(kept).distance_squared.numerator == old);
    CHECK(kept.memory_plan().peak_payload_estimate_bytes() == 28672);
    CHECK_FALSE(kept.memory_plan().proven_peak_payload_upper_bound_bytes());
    const auto exact_steps = kept.steps();
    auto repeated =
        subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, sufficient,
                                   {numerics::linear::EstimateAdmission{28672}, exact_steps});
    REQUIRE(repeated);
    failure(
        subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, sufficient,
                                   {numerics::linear::EstimateAdmission{28672}, exact_steps - 1}),
        core::FailureCode::ResourceExceeded);
    auto moved = std::move(kept);
    // NOLINTNEXTLINE(bugprone-use-after-move): inactive query is the documented move contract.
    CHECK(kept.outcome() == nullptr);
    CHECK(moved.outcome() != nullptr);
}
TEST_CASE("Predicate workspace leases reject concurrent sharing and release on failures",
          "[property][negative][N4-03][N4-09]") {
    auto scratch = workspace();
    const auto input = fixture(.75, 1, 20);
    {
        auto held = take(scratch->acquire());
        core::ResourcePlan plan{owned_budget};
        failure(
            subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, plan, limits),
            core::FailureCode::InvalidInput);
        bool rejected = false;
        std::jthread worker([&] {
            core::ResourcePlan own{owned_budget};
            auto value = subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch,
                                                    own, limits);
            rejected = !value && value.error().code == core::FailureCode::InvalidInput;
        });
        worker.join();
        CHECK(rejected);
    }
    core::ResourcePlan plan{owned_budget};
    auto value =
        subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch, plan, limits);
    REQUIRE(value);
    std::array<bool, 2> agreed{};
    std::array<std::jthread, 2> workers;
    for (std::size_t i = 0; i < 2; ++i)
        workers[i] = std::jthread([&, i] {
            core::ResourcePlan preparation{owned_budget}, output{owned_budget};
            auto own = subject::PredicateWorkspace::create(preparation);
            if (!own)
                return;
            auto result = subject::classify_triangle(subject::Binary64Coordinates{input}, **own,
                                                     output, limits);
            const auto* a = result && result->outcome() != nullptr
                                ? std::get_if<subject::NondegenerateTriangle>(result->outcome())
                                : nullptr;
            const auto* b = value->outcome() != nullptr
                                ? std::get_if<subject::NondegenerateTriangle>(value->outcome())
                                : nullptr;
            if (a != nullptr && b != nullptr) {
                agreed[i] =
                    a->distance_squared.numerator == b->distance_squared.numerator &&
                    a->distance_squared.denominator == b->distance_squared.denominator &&
                    a->distance_squared.unit_exponent == b->distance_squared.unit_exponent &&
                    a->plane == b->plane && a->projection == b->projection &&
                    a->longest_edge_squared == b->longest_edge_squared &&
                    a->minimizing_features == b->minimizing_features && a->far == b->far &&
                    result->steps() == value->steps() &&
                    result->diagnostic_context() == value->diagnostic_context();
            }
        });
    for (auto& worker : workers)
        worker.join();
    CHECK(agreed[0]);
    CHECK(agreed[1]);
}
TEST_CASE("Predicate owners and known context survive every preparation scope",
          "[property][negative][N4-03][N4-08][N4-09]") {
    const core::FailureDetails context{core::ModelGeneration{1},
                                       core::AcceptedStateId{2},
                                       core::TrialId{3},
                                       core::NumericalGeneration{4},
                                       core::ComponentId{5},
                                       core::ProblemId{core::PanelId{6}},
                                       core::ResourceUsage{core::ResourceKind::Iterations, 7, 8},
                                       core::AcceptedStateId{9}};
    auto result = [&] {
        auto input = fixture(.75, 1, 20);
        const auto preserved = input;
        core::ResourcePlan preparation{subject::PredicateWorkspace::payload_bytes()};
        auto scratch = take(subject::PredicateWorkspace::create(preparation, context));
        core::ResourcePlan output{subject::ExactTriangleResult::payload_bytes()};
        auto value = take(subject::classify_triangle(subject::Binary64Coordinates{input}, *scratch,
                                                     output, limits, {}, context));
        CHECK(input == preserved);
        CHECK(preparation.allocated_bytes() == subject::PredicateWorkspace::payload_bytes());
        CHECK(output.allocated_bytes() == subject::ExactTriangleResult::payload_bytes());
        CHECK(value.diagnostic_context() == context);
        core::ResourcePlan no_payload{0};
        const auto denied = subject::classify_triangle(subject::Binary64Coordinates{input},
                                                       *scratch, no_payload, limits, {}, context);
        REQUIRE_FALSE(denied);
        auto expected = context;
        expected.resources = core::ResourceUsage{core::ResourceKind::Bytes,
                                                 subject::ExactTriangleResult::payload_bytes(), 0};
        CHECK(denied.error().details == expected); // Inner actual demand is preserved.
        const auto strict = subject::classify_triangle(
            subject::Binary64Coordinates{input}, *scratch, output,
            {numerics::linear::RequireProvenBound{28672}, 100}, {}, context);
        REQUIRE_FALSE(strict);
        CHECK(strict.error().details == context);
        input.fill(std::numeric_limits<double>::quiet_NaN());
        return value;
    }();
    CHECK(nondegenerate(result).far == subject::FarRelation::Boundary);
    CHECK(result.diagnostic_context() == context);
    core::ResourcePlan short_workspace{subject::PredicateWorkspace::payload_bytes() - 1};
    const auto no_workspace = subject::PredicateWorkspace::create(short_workspace, context);
    REQUIRE_FALSE(no_workspace);
    CHECK(no_workspace.error().code == core::FailureCode::ResourceExceeded);
    CHECK(no_workspace.error().details.model == context.model);
    CHECK(short_workspace.allocated_bytes() == 0);
}
static_assert(!std::is_copy_constructible_v<subject::ExactTriangleResult>);
template <class Workspace>
concept ForgedWorkspace = requires(core::PlannedBuffer<subject::PredicateScratch> storage) {
    Workspace({}, std::move(storage));
};
static_assert(!ForgedWorkspace<subject::PredicateWorkspace>);
static_assert(!std::is_default_constructible_v<subject::PredicateMemoryPlan>);
static_assert(!std::is_constructible_v<subject::PredicateMemoryPlan, std::uint64_t>);
