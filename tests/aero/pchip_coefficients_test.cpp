#include "aero/pchip_coefficients.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <fstream>
#include <future>
#include <iterator>
#include <vector>
namespace {
namespace p = uvlas::aero::detail::pchip;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
template <class T>
T take(uvlas::core::Result<T> value) {
    if (!value)
        FAIL("Unexpected failure " << static_cast<int>(value.error().code));
    return std::move(*value);
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
    if (value.sign < 0)
        n = -n;
    return Rat{n} / d;
}
p::RationalLimits limits() {
    return {4096, 1000000, UINT64_MAX, 100000};
}
struct Samples {
    std::vector<Input> times, values;
    std::vector<p::RationalView> tv, vv;
    Samples(const std::vector<Rat>& t, const std::vector<Rat>& y) {
        times.reserve(t.size());
        values.reserve(y.size());
        for (const auto& v : t)
            times.emplace_back(v);
        for (const auto& v : y)
            values.emplace_back(v);
        for (const auto& v : times)
            tv.push_back(v.view());
        for (const auto& v : values)
            vv.push_back(v.view());
    }
    [[nodiscard]] p::SampleGrid view(std::uint64_t m = 1) const {
        return {tv, vv, m};
    }
};
} // namespace
TEST_CASE("Private PCHIP coefficients own exact two point polynomials",
          "[AE-07][N4-03][N4-06][verification]") {
    const Samples samples{{0, 2}, {1, 3}};
    const auto grid = samples.view();
    const auto plan = take(p::forecast_coefficients(grid, limits()));
    uvlas::core::ResourcePlan memory{plan.payload_bytes()};
    p::CoefficientWorkspace workspace;
    const auto block = take(p::prepare_coefficients(grid, workspace, limits(), memory));
    CHECK(block.nodes() == 2);
    CHECK(block.components() == 1);
    CHECK(memory.allocated_bytes() == plan.payload_bytes());
    CHECK(block.used_steps() <= plan.reserved_steps());
    for (std::uint64_t power = 0; power < 4; ++power)
        CHECK(decoded(take(block.coefficient(0, 0, power))) ==
              std::array<int, 4>{1, 2, 0, 0}[static_cast<std::size_t>(power)]);
    CHECK(decoded(take(block.time(1))) == 2);
    CHECK(decoded(take(block.value(1, 0))) == 3);
}

TEST_CASE(
    "Private PCHIP coefficients match frozen exact witnesses and independent Hermite equations",
    "[AE-07][N4-02][N4-03][N4-06][verification]") {
    std::ifstream stream{UVLAS_PCHIP_COEFFICIENT_REFERENCE};
    REQUIRE(stream.good());
    const auto reference = nlohmann::json::parse(stream);
    CHECK(reference.at("source_sha256") ==
          "497287324f6084d162a41aa98bfad9993adb91287fb78e6e5467adef9a2e22c5");
    for (const auto& record : reference.at("cases")) {
        CAPTURE(record.at("name"));
        std::vector<Rat> times, values, slopes;
        for (const auto& value : record.at("times"))
            times.emplace_back(value.get<std::string>());
        for (const auto& value : record.at("values"))
            values.emplace_back(value.get<std::string>());
        for (const auto& value : record.at("slopes"))
            slopes.emplace_back(value.get<std::string>());
        const auto components = record.at("components").get<std::uint64_t>();
        const auto block = [&] {
            const Samples samples{times, values};
            const auto grid = samples.view(components);
            const auto plan = take(p::forecast_coefficients(grid, limits()));
            uvlas::core::ResourcePlan memory{plan.payload_bytes()};
            p::CoefficientWorkspace workspace;
            auto owner = take(p::prepare_coefficients(grid, workspace, limits(), memory));
            CHECK(memory.allocated_bytes() == plan.payload_bytes());
            CHECK(owner.used_steps() <= plan.reserved_steps());
            CHECK(owner.used_slots() <= plan.arithmetic().slots());
            return owner;
        }();
        // Source borrows, scratch, and ResourcePlan have all ended before inspecting the owner.
        CHECK(block.nodes() == times.size());
        CHECK(block.components() == components);
        for (std::uint64_t i = 0; i < times.size(); ++i) {
            CHECK(decoded(take(block.time(i))) == times[static_cast<std::size_t>(i)]);
            for (std::uint64_t c = 0; c < components; ++c)
                CHECK(decoded(take(block.value(i, c))) ==
                      values[static_cast<std::size_t>(i * components + c)]);
        }
        for (std::uint64_t i = 0; i + 1 < times.size(); ++i)
            for (std::uint64_t c = 0; c < components; ++c) {
                std::array<Rat, 4> coefficient;
                for (std::uint64_t power = 0; power < 4; ++power) {
                    coefficient[static_cast<std::size_t>(power)] =
                        decoded(take(block.coefficient(i, c, power)));
                    CHECK(coefficient[static_cast<std::size_t>(power)] ==
                          Rat{record.at("coefficients")
                                  .at(static_cast<std::size_t>(i * components + c))
                                  .at(static_cast<std::size_t>(power))
                                  .get<std::string>()});
                }
                const auto at = static_cast<std::size_t>(i * components + c),
                           next = static_cast<std::size_t>((i + 1) * components + c);
                const Rat h =
                    times[static_cast<std::size_t>(i + 1)] - times[static_cast<std::size_t>(i)];
                CHECK(coefficient[0] == values[at]);
                CHECK(coefficient[0] + coefficient[1] + coefficient[2] + coefficient[3] ==
                      values[next]);
                CHECK(coefficient[1] / h == slopes[at]);
                CHECK((coefficient[1] + 2 * coefficient[2] + 3 * coefficient[3]) / h ==
                      slopes[next]);
            }
    }
}
TEST_CASE(
    "Private PCHIP coefficient forecast covers alternative branches without duplicating prefix",
    "[N4-02][N4-06][negative]") {
    p::Counter left{{1, 1, 67, 3}};
    const auto original = take(left.integer(1));
    auto right = left;
    take(left.unary(p::RationalOp::Negate, original));
    take(right.unary(p::RationalOp::Clone, original));
    REQUIRE(left.cover_alternative(right));
    const auto plan = take(left.finish());
    CHECK(plan.slots() == 3);
    CHECK(plan.reserved_steps() == 67);
    const p::Counter incompatible{{1, 2, 67, 3}};
    const auto mismatch = left.cover_alternative(incompatible);
    REQUIRE_FALSE(mismatch);
    CHECK(mismatch.error().code == uvlas::core::FailureCode::InvalidInput);
    CHECK(take(left.finish()).reserved_steps() == 67);
}
TEST_CASE("Private PCHIP coefficient preparation rejects malformed grids and keeps existing owners",
          "[AE-07][N4-02][N4-03][N4-08][negative]") {
    const Samples good{{0, 1, 2}, {0, 1, 4}}, duplicate{{0, 0, 1}, {0, 1, 2}},
        descending{{0, 2, 1}, {0, 1, 2}}, short_grid{{0}, {1}};
    const auto plan = take(p::forecast_coefficients(good.view(), limits()));
    uvlas::core::ResourcePlan memory{plan.payload_bytes() * 8};
    p::CoefficientWorkspace workspace;
    auto owner = take(p::prepare_coefficients(good.view(), workspace, limits(), memory));
    const auto old = take(owner.coefficient(0, 0, 0));
    for (const auto grid :
         {duplicate.view(), descending.view(), short_grid.view(), good.view(0), good.view(2)}) {
        const auto failed = p::prepare_coefficients(grid, workspace, limits(), memory);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == uvlas::core::FailureCode::InvalidInput);
        CHECK(decoded(old) == 0);
    }
    const auto overflow = p::forecast_coefficients(good.view(UINT64_MAX), limits());
    REQUIRE_FALSE(overflow);
    CHECK(overflow.error().code == uvlas::core::FailureCode::ResourceExceeded);
    const auto invalid = owner.coefficient(0, 0, 4);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == uvlas::core::FailureCode::InvalidInput);
    CHECK_FALSE(owner.time(3));
    CHECK_FALSE(owner.value(0, 1));
    std::stop_source stop;
    stop.request_stop();
    const auto cancelled =
        p::prepare_coefficients(good.view(), workspace, limits(), memory, stop.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == uvlas::core::FailureCode::Cancelled);
    const auto recovered = take(p::prepare_coefficients(good.view(), workspace, limits(), memory));
    CHECK(decoded(take(recovered.coefficient(0, 0, 0))) == 0);
}

TEST_CASE(
    "Private PCHIP coefficient word representations and moved owners retain explicit contracts",
    "[N4-03][N4-08][negative]") {
    const Samples samples{{0, 1, 2}, {0, 1, 4}};
    const std::array<std::uint32_t, 1> one{1}, zero{0};
    const std::array<std::uint32_t, 2> leading{1, 0};
    for (const auto bad : std::array<p::RationalView, 5>{
             {{1, one, {}}, {1, {}, one}, {0, one, one}, {1, leading, one}, {1, one, zero}}}) {
        auto values = samples.vv;
        values[0] = bad;
        const p::SampleGrid grid{samples.tv, values, 1};
        const auto result = p::forecast_coefficients(grid, limits());
        REQUIRE_FALSE(result);
        CHECK(result.error().code == uvlas::core::FailureCode::InvalidInput);
    }
    const auto plan = take(p::forecast_coefficients(samples.view(), limits()));
    uvlas::core::ResourcePlan memory{plan.payload_bytes()};
    p::CoefficientWorkspace workspace;
    uvlas::core::FailureDetails details;
    details.component = uvlas::core::ComponentId{17};
    auto owner =
        take(p::prepare_coefficients(samples.view(), workspace, limits(), memory, {}, details));
    const auto borrowed = take(owner.value(2, 0));
    auto moved = std::move(owner);
    // NOLINTBEGIN(bugprone-use-after-move): explicitly supported inactive owner queries.
    CHECK_FALSE(owner.active());
    CHECK(owner.nodes() == 0);
    CHECK(owner.components() == 0);
    CHECK(owner.used_steps() == 0);
    CHECK(owner.used_slots() == 0);
    CHECK_FALSE(owner.time(0));
    // NOLINTEND(bugprone-use-after-move)
    CHECK(decoded(borrowed) == 4);
    CHECK(decoded(take(moved.value(2, 0))) == 4);
    const auto invalid = moved.coefficient(0, 0, 4);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().details.component == details.component);
}
TEST_CASE("Private PCHIP coefficient owners are identical for one two and four callers",
          "[N4-03][N4-09][property]") {
    struct Encoded {
        std::vector<std::uint32_t> words;
        std::uint64_t work, slots, payload;
        bool operator==(const Encoded&) const = default;
    };
    const auto evaluate = []() -> std::optional<Encoded> {
        const Samples samples{{0, 1, 3, 4}, {0, 3, 2, 5}};
        const auto plan = p::forecast_coefficients(samples.view(), limits());
        if (!plan)
            return std::nullopt;
        uvlas::core::ResourcePlan memory{plan->payload_bytes()};
        p::CoefficientWorkspace workspace;
        const auto owner = p::prepare_coefficients(samples.view(), workspace, limits(), memory);
        if (!owner)
            return std::nullopt;
        Encoded result{{}, owner->used_steps(), owner->used_slots(), memory.allocated_bytes()};
        for (std::uint64_t i = 0; i + 1 < owner->nodes(); ++i)
            for (std::uint64_t power = 0; power < 4; ++power) {
                const auto view = owner->coefficient(i, 0, power);
                if (!view)
                    return std::nullopt;
                result.words.push_back(static_cast<std::uint32_t>(view->sign + 1));
                result.words.push_back(static_cast<std::uint32_t>(view->numerator.size()));
                result.words.push_back(static_cast<std::uint32_t>(view->denominator.size()));
                result.words.insert(result.words.end(), view->numerator.begin(),
                                    view->numerator.end());
                result.words.insert(result.words.end(), view->denominator.begin(),
                                    view->denominator.end());
            }
        return result;
    };
    const auto baseline = evaluate();
    if (!baseline) {
        FAIL_CHECK("Missing baseline owner");
        return;
    }
    const auto& expected = baseline.value();
    for (const auto workers : std::array<std::size_t, 3>{1, 2, 4}) {
        std::vector<std::future<std::optional<Encoded>>> pending;
        pending.reserve(workers);
        for (std::size_t i = 0; i < workers; ++i)
            pending.push_back(std::async(std::launch::async, evaluate));
        for (auto& future : pending) {
            const auto result = future.get();
            if (!result) {
                FAIL_CHECK("Missing independent owner");
                continue;
            }
            const auto& actual = result.value();
            CHECK(actual == expected);
        }
    }
}
