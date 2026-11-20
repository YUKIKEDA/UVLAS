#include "aero/pchip_preparation_input.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <future>
#include <iterator>
#include <type_traits>
#include <vector>

namespace {
namespace c = uvlas::core;
namespace p = uvlas::aero::detail::pchip;
namespace n = uvlas::numerics;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
constexpr p::RationalLimits limits{4096, 10000, UINT64_MAX, 10000};
constexpr p::PchipInputStamp zero{p::PchipIssuerId{0}, p::PchipInputGeneration{0}};
constexpr n::SiDimension velocity{{1, 0, -1, 0, 0, 0, 0}};
static_assert(!std::is_convertible_v<p::PchipInputGeneration, p::PchipIssuerId>);
static_assert(!std::is_convertible_v<p::PchipIssuerId, p::PchipInputGeneration>);

template <class T>
T take(c::Result<T> result) {
    if (!result)
        FAIL("Unexpected failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}

struct Input {
    int sign;
    std::vector<std::uint32_t> numerator, denominator;

    explicit Input(const Rat& value) : sign(value < 0 ? -1 : value == 0 ? 0 : 1) {
        const Big magnitude = boost::multiprecision::abs(boost::multiprecision::numerator(value));
        if (sign != 0)
            boost::multiprecision::export_bits(magnitude, std::back_inserter(numerator), 32, false);
        boost::multiprecision::export_bits(boost::multiprecision::denominator(value),
                                           std::back_inserter(denominator), 32, false);
    }
    [[nodiscard]] p::RationalView view() const& noexcept {
        return {sign, numerator, denominator};
    }
    [[nodiscard]] p::RationalView view() const&& = delete;
};

Rat decode(p::RationalView view) {
    Big n = 0, d = 0;
    boost::multiprecision::import_bits(d, view.denominator.begin(), view.denominator.end(), 32,
                                       false);
    if (!view.numerator.empty())
        boost::multiprecision::import_bits(n, view.numerator.begin(), view.numerator.end(), 32,
                                           false);
    if (view.sign < 0)
        n = -n;
    return Rat{n} / d;
}

std::vector<Input> input_words(const std::vector<Rat>& data) {
    std::vector<Input> words;
    words.reserve(data.size());
    for (const auto& item : data)
        words.emplace_back(item);
    return words;
}

std::vector<p::RationalView> views(const std::vector<Input>& input) {
    std::vector<p::RationalView> result;
    result.reserve(input.size());
    for (const auto& value : input)
        result.push_back(value.view());
    return result;
}

using Policy = decltype(p::BoundaryRequest::policy);

struct Definition {
    std::vector<Rat> times{0, 1, 2, 3};
    std::vector<Rat> values{0, 3, 1, 4, -1, 2, 0, 3};
    std::vector<Rat> events{-7, Rat{1} / 2, 9};
    std::vector<n::SiDimension> dimensions{velocity, velocity};
    Rat begin = -2, end = 5, period = 3;
};

c::Result<p::PchipPreparationInput> make(const Definition& definition,
                                         p::PchipInputStamp stamp = zero, bool periodic = true) {
    const auto t = input_words(definition.times), y = input_words(definition.values);
    const auto e = input_words(definition.events);
    const auto tv = views(t), yv = views(y), ev = views(e);
    const Input begin{definition.begin}, end{definition.end}, period{definition.period};
    const auto outside =
        periodic ? Policy{p::PeriodicBoundary{period.view()}} : Policy{p::HoldBoundary{}};
    const p::PreparationRequest request{{tv, yv, definition.dimensions.size()},
                                        definition.dimensions,
                                        {begin.view(), end.view(), outside},
                                        ev,
                                        stamp};
    p::PreparationWorkspace workspace;
    c::ResourcePlan memory{UINT64_MAX};
    return p::make_preparation_input(request, workspace, limits, memory);
}

template <class T>
void refuses(const c::Result<T>& result, c::FailureCode code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    CHECK(result.error().operation == c::Operation::EvaluateAero);
}

} // namespace

TEST_CASE("PCHIP owns the complete preparation input and connects it to curve checks",
          "[AE-07][N4-03][N4-06][V-AERO][verification]") {
    auto original = take(make(Definition{})); // All borrows and preparation resources are dead.
    const auto saved = take(original.view());
    auto input = std::move(original);
    // NOLINTNEXTLINE(bugprone-use-after-move): verify documented inactive owner.
    CHECK_FALSE(original.active());
    // NOLINTNEXTLINE(bugprone-use-after-move): verify documented inactive owner.
    refuses(original.view(), c::FailureCode::InvalidInput);
    // NOLINTNEXTLINE(bugprone-use-after-move): verify documented inactive owner.
    refuses(original.formula_version(), c::FailureCode::InvalidInput);
    const auto owned = take(input.view());
    const Definition expected;
    REQUIRE(owned.samples.times.size() == expected.times.size());
    REQUIRE(owned.samples.values.size() == expected.values.size());
    REQUIRE(owned.events.size() == expected.events.size());
    REQUIRE(owned.dimensions.size() == expected.dimensions.size());
    CHECK(owned.samples.times.data() == saved.samples.times.data());
    CHECK(owned.stamp == zero);
    CHECK(take(input.formula_version()) == p::CoefficientBlock::formula_version());
    for (std::size_t i = 0; i < expected.times.size(); ++i)
        CHECK(decode(owned.samples.times[i]) == expected.times[i]);
    for (std::size_t i = 0; i < expected.values.size(); ++i)
        CHECK(decode(owned.samples.values[i]) == expected.values[i]);
    for (std::size_t i = 0; i < expected.events.size(); ++i)
        CHECK(decode(owned.events[i]) == expected.events[i]);
    for (const auto& unit : owned.dimensions)
        CHECK(unit.powers == velocity.powers);
    CHECK(decode(owned.boundary.begin) == expected.begin);
    CHECK(decode(owned.boundary.end) == expected.end);
    if (const auto* periodic = std::get_if<p::PeriodicBoundary>(&owned.boundary.policy)) {
        CHECK(decode(periodic->period) == expected.period);
    } else {
        FAIL("Expected an owned explicit period");
    }

    // Use only owned views. A successful preparation input can drive both existing checks.
    c::ResourcePlan remaining{UINT64_MAX};
    p::CoefficientWorkspace coefficients;
    const auto block =
        take(p::prepare_coefficients(owned.samples, coefficients, limits, remaining));
    p::BoundaryWorkspace boundary;
    CHECK(p::check_boundary(block, owned.boundary, boundary, limits, remaining).has_value());
}

TEST_CASE("PCHIP compares independent owners by both stamp and all input fields",
          "[AE-07][N4-03][N4-06][N4-08][verification]") {
    const auto at = take(make(Definition{}));
    const auto same = take(make(Definition{}));
    CHECK(take(at.view()).samples.times.data() != take(same.view()).samples.times.data());
    const auto work = take(p::match_preparation_input(at, same, UINT64_MAX));
    CHECK(take(p::match_preparation_input(at, same, work)) == work);
    refuses(p::match_preparation_input(at, same, work - 1), c::FailureCode::ResourceExceeded);
    for (const auto stamp : {p::PchipInputStamp{p::PchipIssuerId{1}, p::PchipInputGeneration{0}},
                             p::PchipInputStamp{p::PchipIssuerId{0}, p::PchipInputGeneration{1}}}) {
        const auto different = take(make(Definition{}, stamp));
        refuses(p::match_preparation_input(at, different, UINT64_MAX),
                c::FailureCode::StaleGeneration);
    }
    const auto other_policy = take(make(Definition{}, zero, false));
    refuses(p::match_preparation_input(at, other_policy, UINT64_MAX),
            c::FailureCode::StaleGeneration);

    const Rat low_word = Rat{1} / (Big{1} << 200);
    for (const unsigned field : {0U, 1U, 2U, 3U, 4U, 5U, 6U, 7U}) {
        CAPTURE(field);
        Definition other;
        switch (field) {
        case 0:
            other.times[1] += low_word;
            break;
        case 1:
            other.values.back() += low_word;
            break;
        case 2:
            other.dimensions.back().powers[2] = -2;
            break;
        case 3:
            other.begin += low_word;
            break;
        case 4:
            other.end += low_word;
            break;
        case 5:
            other.events.back() += low_word;
            break;
        case 6:
            other.events.pop_back();
            break;
        case 7:
            other.times.back() += low_word;
            other.period += low_word;
            break;
        default:
            FAIL("Unexpected field");
            break;
        }
        const auto different = take(make(other));
        refuses(p::match_preparation_input(at, different, UINT64_MAX),
                c::FailureCode::StaleGeneration);
    }
    for (const std::size_t count : {1UZ, 2UZ, 4UZ}) {
        std::vector<std::future<c::Result<std::uint64_t>>> jobs;
        jobs.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            jobs.push_back(std::async(std::launch::async,
                                      [&] { return p::match_preparation_input(at, same, work); }));
        for (auto& job : jobs)
            CHECK(take(job.get()) == work);
    }
}

TEST_CASE("PCHIP rejects malformed input before publishing an owner",
          "[AE-07][N4-02][N4-08][negative]") {
    for (const unsigned defect : {0U, 1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U}) {
        CAPTURE(defect);
        Definition data;
        switch (defect) {
        case 0:
            data.times[1] = data.times[0];
            break;
        case 1:
            data.times[1] = data.times[0] - 1;
            break;
        case 2:
            data.events[1] = data.events[0];
            break;
        case 3:
            data.events[1] = data.events[0] - 1;
            break;
        case 4:
            data.begin = data.end + 1;
            break;
        case 5:
            data.period = 0;
            break;
        case 6:
            data.period += Rat{1} / (Big{1} << 200);
            break;
        case 7:
            data.values.pop_back();
            break;
        case 8:
            data.dimensions.pop_back();
            break;
        default:
            FAIL("Unexpected defect");
            break;
        }
        refuses(make(data), c::FailureCode::InvalidInput);
    }

    const Definition original;
    const auto t = input_words(original.times), y = input_words(original.values);
    auto tv = views(t);
    const auto yv = views(y);
    const Input begin{-2}, end{5}, period{3};
    const std::array<std::uint32_t, 1> two{2}, four{4};
    p::PreparationWorkspace workspace;
    for (const bool empty_period : {false, true}) {
        CAPTURE(empty_period);
        tv[1] = empty_period ? t[1].view() : p::RationalView{1, two, four};
        const p::PreparationRequest bad{
            {tv, yv, 2},
            original.dimensions,
            {begin.view(), end.view(),
             p::PeriodicBoundary{empty_period ? p::RationalView{} : period.view()}},
            {},
            zero};
        c::ResourcePlan memory{UINT64_MAX};
        refuses(p::make_preparation_input(bad, workspace, limits, memory),
                c::FailureCode::InvalidInput);
        if (empty_period)
            CHECK(memory.allocated_bytes() == 0);
    }
}
