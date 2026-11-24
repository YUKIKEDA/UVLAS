#include "aero/geometry_snapshot_fixture.hpp"
#include "coupling/case_flow_clock.hpp"
#include "coupling/case_known_flow.hpp"
#include "coupling/case_uniform_flow.hpp"
#include "coupling/case_uniform_flow_projection.hpp"
#include "coupling/discrete_nonpenetration.hpp"
#include "io/case_schema_test_adapter.hpp"
#include "io/si_rational_projection.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
namespace a = aero::detail;
namespace p = a::pchip;
namespace cf = coupling::detail;

constexpr std::uint64_t budget = std::uint64_t{1} << 48U;
const p::PchipInputStamp stamp{p::PchipIssuerId{71}, p::PchipInputGeneration{8}};

structure::ExactLimits exact_limits() noexcept {
    return {100000,
            100000,
            100000,
            64,
            64,
            100000,
            200000,
            structure::ExactEstimateAdmission{budget},
            structure::ExactTimingMode::Disabled};
}
const cf::CaseUniformFlowLimits limits{
    exact_limits(),
    {{128, 100000, budget, 10000}, numerics::linear::EstimateAdmission{budget}},
    100};

template <class T>
T take(core::Result<T> value) {
    if (!value)
        FAIL("Unexpected case-flow failure " << static_cast<int>(value.error().code));
    return std::move(*value);
}

template <class T, class Diagnostic>
T take_io(core::Result<T> value, const Diagnostic& diag) {
    if (!value) {
        const auto path = diag.path();
        if (path)
            FAIL("Unexpected case-flow input at " << *path << " code "
                                                  << static_cast<int>(value.error().code));
        FAIL("Unexpected case-flow preflight " << static_cast<int>(value.error().code));
    }
    return std::move(*value);
}

template <class T>
void reject(const core::Result<T>& value, core::FailureCode code) {
    if (value)
        FAIL("Unexpected case-flow success");
    CHECK(value.error().code == code);
}

io::SiCase materialize(const JsonToken& source, io::SourceId id) {
    core::ResourcePlan memory{budget};
    const auto wire = source.wire();
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    std::optional<io::InputDiagnostic> parsing;
    const auto doc = take(io::read_input_document(
        std::as_bytes(std::span{wire.data(), wire.size()}), id, syntax, memory, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(budget);
    io::SchemaDiagnostic schema_diagnostic;
    const auto checked = take_io(
        io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, memory, &schema_diagnostic),
        schema_diagnostic);

    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(budget);
    const io::SiLimits si_limits{si_caps, exact_limits()};
    io::SiDiagnostic si_diagnostic;
    const auto catalog =
        take_io(io::make_si_catalog(checked, si_limits, memory, &si_diagnostic), si_diagnostic);
    std::vector<io::SiRequirementRow> requirements;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto quantity = take(catalog.entry(i));
        requirements.push_back({quantity.key,
                                {io::SiDecimalInput{"1"}, quantity.recipe.dimensions},
                                io::SiDecimalInput{"1e-25"},
                                io::SiDecimalInput{"1e-25"}});
    }
    const auto request =
        take_io(io::make_si_request(catalog, requirements, si_limits, memory, &si_diagnostic),
                si_diagnostic);
    return take_io(
        io::convert_case_to_si(checked, catalog, request, si_limits, memory, &si_diagnostic),
        si_diagnostic);
}

// Keep every source number as its original JSON token before entering the native parser.
JsonToken source(std::string_view series, std::string_view intervals) {
    auto input = root("coupled_transient");
    input.field("units").set("length", string("ft"));
    input.field("units").set("time", string("min"));
    input.field("flow").set("uniform", parse(series));
    input.field("analysis").set("intervals", parse(intervals));
    return input;
}

constexpr std::string_view two_step =
    R"([{"start":0,"end":1,"step":1,"restart":"initial"},{"start":1,"end":2,"step":1,"restart":"event"}])";
constexpr std::string_view one_step = R"([{"start":0,"end":1,"step":1,"restart":"initial"}])";

cf::CaseUniformFlowReady import(const JsonToken& input, io::SourceId id,
                                a::UniformFlowWorkspace& workspace,
                                p::PchipInputStamp identity = stamp) {
    const auto si = materialize(input, id);
    core::ResourcePlan memory{budget};
    return take(cf::prepare_case_uniform_flow(si, identity, {}, workspace, limits, memory));
}

struct Number {
    Number(std::uint32_t numerator, std::uint32_t denominator = 1, int sign = 1) noexcept
        : n{numerator}, d{denominator}, sign_(numerator == 0 ? 0 : sign) {}
    [[nodiscard]] p::RationalView view() const& noexcept {
        return {sign_, sign_ == 0 ? std::span<const std::uint32_t>{} : n, d};
    }
    [[nodiscard]] p::RationalView view() const&& = delete;
    std::array<std::uint32_t, 1> n, d;
    int sign_;
};

void exact(p::RationalView a, p::RationalView b) {
    CHECK(a.sign == b.sign);
    CHECK(std::ranges::equal(a.numerator, b.numerator));
    CHECK(std::ranges::equal(a.denominator, b.denominator));
}
void exact(p::RationalView actual, const Number& expected) {
    exact(actual, expected.view());
}

cf::CaseUniformFlowEvaluation evaluate(const cf::CaseUniformFlowReady& field,
                                       const cf::CaseUniformFlowReady& current,
                                       a::UniformFlowWorkspace& workspace, const Number& time,
                                       p::RequestedSide side = p::RequestedSide::Automatic) {
    const cf::CaseUniformFlowQuery query{stamp, time.view(), side};
    const auto forecast = take(cf::forecast_case_uniform_flow(field, current, query, limits.flow));
    auto exact_limits = limits.flow;
    exact_limits.arithmetic.steps = forecast.reserved_steps;
    exact_limits.metadata = numerics::linear::EstimateAdmission{forecast.metadata_estimate_bytes};
    core::ResourcePlan memory{forecast.payload_bytes};
    auto result = take(
        cf::evaluate_case_uniform_flow(field, current, query, workspace, exact_limits, memory));
    CHECK(result.used_steps <= forecast.reserved_steps);
    CHECK(memory.allocated_bytes() == forecast.payload_bytes);
    --exact_limits.arithmetic.steps;
    core::ResourcePlan refused{forecast.payload_bytes};
    reject(cf::evaluate_case_uniform_flow(field, current, query, workspace, exact_limits, refused),
           core::FailureCode::ResourceExceeded);
    CHECK(refused.allocated_bytes() == 0);
    return result;
}

void jet(const a::ExactUniformFlowJet& result, const std::array<Number, 3>& velocity,
         const std::array<Number, 3>& derivative, const Number& time) {
    CHECK(take(result.stamp()) == stamp);
    exact(take(result.time()), time);
    for (std::size_t i = 0; i < 3; ++i) {
        exact(take(result.velocity(i)), velocity[i]);
        exact(take(result.time_partial(i)), derivative[i]);
        for (std::size_t j = 0; j < 3; ++j)
            exact(take(result.spatial_partial(i, j)), Number{0});
    }
}
} // namespace

TEST_CASE("Typed SI constant flow owns full-width signed magnitudes and restores exact velocity",
          "[AE-07][IO4-01][N4-02][N4-03][property]") {
    const auto input =
        source(R"({"kind":"constant","value":[18446744073709551617,-2.5,0]})", one_step);
    a::UniformFlowWorkspace work;
    const auto field = import(input, io::SourceId{1}, work);
    const auto current = import(input, io::SourceId{2}, work);
    const Number time{60};
    core::ResourcePlan save_memory{budget};
    const auto saved = take(cf::snapshot_case_uniform_flow(field, work, limits.flow, save_memory));
    const auto replay = take(
        cf::restore_case_uniform_flow(take(saved.view()), current, work, limits.flow, save_memory));
    const auto value = evaluate(replay, current, work, time);
    // ft/min = (381/1250)/60 m/s = 127/25000 m/s. First coefficient is 2^64+1.
    constexpr std::array<std::uint32_t, 3> large_words{127, 0, 127};
    constexpr std::array<std::uint32_t, 1> denominator{25000};
    exact(take(value.jet.velocity(0)), {1, large_words, denominator});
    exact(take(value.jet.velocity(1)), Number{127, 10000, -1});
    exact(take(value.jet.velocity(2)), Number{0});
    for (std::size_t axis = 0; axis < 3; ++axis)
        exact(take(value.jet.time_partial(axis)), Number{0});
    CHECK(value.used_steps > value.jet.used_steps());

    auto altered = input;
    altered.field("flow").field("uniform").field("value").elements[0] =
        parse("18446744073709551616");
    const auto changed = import(altered, io::SourceId{3}, work);
    const cf::CaseUniformFlowQuery query{stamp, time.view()};
    core::ResourcePlan refused{budget};
    reject(cf::forecast_case_uniform_flow(field, changed, query, limits.flow),
           core::FailureCode::StaleGeneration);
    reject(cf::evaluate_case_uniform_flow(field, changed, query, work, limits.flow, refused),
           core::FailureCode::StaleGeneration);
    reject(cf::restore_case_uniform_flow(take(saved.view()), changed, work, limits.flow, refused),
           core::FailureCode::StaleGeneration);
    CHECK(refused.allocated_bytes() == 0);
}

TEST_CASE("Typed SI polynomial flow converts time origin scale coefficients and events exactly",
          "[AE-07][IO4-01][N4-03][N4-09][property][negative]") {
    auto input = source(
        R"({"kind":"polynomial","origin":1,"scale":2,"domain":[0,3],"coefficients":[[1,-2,0],[2,0,3],[-1,4,-2],[3,0,1]]})",
        two_step);
    input.field("flow").set("events", parse("[1]"));
    input.field("analysis").field("intervals").elements[1].set("end", parse("3"));
    a::UniformFlowWorkspace work;
    const auto field = import(input, io::SourceId{4}, work);
    const auto current = import(input, io::SourceId{5}, work);
    const auto converted = take(field.uniform().input());
    if (converted.events.size() != 1)
        FAIL("Expected one retained event");
    exact(converted.events[0], Number{60});
    const auto* polynomial = std::get_if<a::PolynomialUniformFlow>(&converted.series);
    if (!polynomial)
        FAIL("Expected imported polynomial");
    exact(polynomial->origin, Number{60});
    exact(polynomial->scale, Number{120});
    exact(polynomial->end, Number{180});

    core::ResourcePlan memory{budget};
    const auto saved = take(cf::snapshot_case_uniform_flow(field, work, limits.flow, memory));
    const auto restored =
        take(cf::restore_case_uniform_flow(take(saved.view()), current, work, limits.flow, memory));
    const Number time{120};
    const auto result = evaluate(restored, current, work, time, p::RequestedSide::Right);
    // z=(120-60)/120=1/2. Differentiate each source coefficient before ft/min to SI.
    jet(result.jet, {Number{2159, 200000}, Number{127, 25000, -1}, Number{1143, 200000}},
        {Number{1651, 12000000}, Number{127, 750000}, Number{889, 12000000}}, time);
    const Number outside{181};
    reject(cf::evaluate_case_uniform_flow(restored, current, {stamp, outside.view()}, work,
                                          limits.flow, memory),
           core::FailureCode::InvalidInput);
    const auto next = p::PchipInputStamp{stamp.issuer, p::PchipInputGeneration{9}};
    const auto new_generation = import(input, io::SourceId{6}, work, next);
    reject(
        cf::forecast_case_uniform_flow(restored, new_generation, {stamp, time.view()}, limits.flow),
        core::FailureCode::StaleGeneration);
}

TEST_CASE("Typed SI polynomial covers both ends of the transient analysis before publication",
          "[AE-07][IO4-01][N4-08][negative]") {
    const auto input = source(
        R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,2],"coefficients":[[1,2,3]]})",
        R"([{"start":0,"end":1,"step":1,"restart":"initial"},{"start":1,"end":2,"step":1,"restart":"step_change"}])");
    a::UniformFlowWorkspace work;
    const auto check_bad = [&](std::string_view domain, io::SourceId id) {
        auto edited = input;
        edited.field("flow").field("uniform").set("domain", parse(domain));
        const auto si = materialize(edited, id);
        core::ResourcePlan memory{budget};
        const auto refusal = cf::prepare_case_uniform_flow(si, stamp, {}, work, limits, memory);
        reject(refusal, core::FailureCode::InvalidInput);
        if (!refusal)
            CHECK(refusal.error().operation == core::Operation::PrepareModel);
    };
    check_bad("[0.5,2]", io::SourceId{18});
    check_bad("[0,1]", io::SourceId{19});

    const auto at_limits = materialize(input, io::SourceId{20});
    auto short_limits = limits;
    short_limits.maximum_exact_numbers = 8;
    core::ResourcePlan memory{budget};
    reject(cf::prepare_case_uniform_flow(at_limits, stamp, {}, work, short_limits, memory),
           core::FailureCode::ResourceExceeded);
    CHECK(memory.allocated_bytes() == 0);
    ++short_limits.maximum_exact_numbers;
    const auto control =
        take(cf::prepare_case_uniform_flow(at_limits, stamp, {}, work, short_limits, memory));
    const auto current =
        take(cf::prepare_case_uniform_flow(at_limits, stamp, {}, work, short_limits, memory));
    const Number begin{0}, end{120};
    CHECK(evaluate(control, current, work, begin).jet.active());
    CHECK(evaluate(control, current, work, end).jet.active());

    auto steady = root("aerodynamic_steady");
    steady.field("flow").set(
        "uniform",
        parse(
            R"({"kind":"polynomial","origin":1,"scale":1,"domain":[1,2],"coefficients":[[1,2,3]]})"));
    const auto steady_flow = import(steady, io::SourceId{21}, work);
    CHECK(steady_flow.uniform().active());
}

TEST_CASE("Typed SI PCHIP uses the declared closed time window and independent current case",
          "[AE-07][IO4-01][N4-02][N4-03][N4-08][N4-09][property][negative]") {
    auto input = source(
        R"({"kind":"pchip","times":[1,2,3],"values":[[0,-2,1],[1,-2,-1],[4,-2,-7]],"outside":{"kind":"reject"}})",
        R"([{"start":1,"end":2,"step":1,"restart":"initial"},{"start":2,"end":3,"step":1,"restart":"event"}])");
    input.field("flow").set("events", parse("[2]"));
    a::UniformFlowWorkspace work;
    const auto field = import(input, io::SourceId{7}, work);
    const auto current = import(input, io::SourceId{8}, work);
    const auto converted = take(field.uniform().input());
    const auto* pchip = std::get_if<a::PchipUniformFlow>(&converted.series);
    if (!pchip || converted.events.size() != 1)
        FAIL("Expected typed PCHIP with one retained event");
    exact(pchip->boundary.begin, Number{60});
    exact(pchip->boundary.end, Number{180});
    exact(converted.events[0], Number{120});

    core::ResourcePlan memory{budget};
    const auto saved = take(cf::snapshot_case_uniform_flow(field, work, limits.flow, memory));
    const auto saved_view = take(saved.view());
    const auto replay =
        take(cf::restore_case_uniform_flow(saved_view, current, work, limits.flow, memory));
    // In local source minutes x(u)=3u²/2-u³/2 on [0,1]. The third entry is 1-2x.
    const Number time{90};
    const auto result = evaluate(replay, current, work, time);
    jet(result.jet, {Number{127, 80000}, Number{127, 12500, -1}, Number{381, 200000}},
        {Number{381, 4000000}, Number{0}, Number{381, 2000000, -1}}, time);

    const Number knot{120};
    for (const auto side : {p::RequestedSide::Left, p::RequestedSide::Right}) {
        const auto at_knot = evaluate(replay, current, work, knot, side);
        jet(at_knot.jet, {Number{127, 25000}, Number{127, 12500, -1}, Number{127, 25000, -1}},
            {Number{127, 1000000}, Number{0}, Number{127, 500000, -1}}, knot);
    }

    auto changed_events = input;
    changed_events.field("flow").set("events", parse("[]"));
    changed_events.field("analysis")
        .field("intervals")
        .elements[1]
        .set("restart", string("step_change"));
    auto changed_sample = input;
    changed_sample.field("flow").field("uniform").field("values").elements[1].elements[0] =
        parse("2");
    for (const auto& changed : {std::cref(changed_sample), std::cref(changed_events)}) {
        const auto other = import(changed.get(), io::SourceId{9}, work);
        const cf::CaseUniformFlowQuery query{stamp, time.view()};
        core::ResourcePlan refused{budget};
        reject(cf::forecast_case_uniform_flow(replay, other, query, limits.flow),
               core::FailureCode::StaleGeneration);
        reject(cf::evaluate_case_uniform_flow(replay, other, query, work, limits.flow, refused),
               core::FailureCode::StaleGeneration);
        reject(cf::restore_case_uniform_flow(saved_view, other, work, limits.flow, refused),
               core::FailureCode::StaleGeneration);
        CHECK(refused.allocated_bytes() == 0);
    }
}

TEST_CASE("Typed SI PCHIP imports held endpoints periodic length and explicit nontransient demand",
          "[AE-07][IO4-01][N4-03][property][negative]") {
    auto input = source(
        R"({"kind":"pchip","times":[0,1,2,3,4],"values":[[0,1,-1],[0,1,-1],[2,1,-3],[0,1,-1],[0,1,-1]],"outside":{"kind":"hold"}})",
        R"([{"start":0,"end":4,"step":1,"restart":"initial"}])");
    a::UniformFlowWorkspace work;
    const Number outside{120, 1, -1};
    const std::array zeros{Number{0}, Number{0}, Number{0}};
    const auto hold = import(input, io::SourceId{13}, work);
    const auto hold_current = import(input, io::SourceId{14}, work);
    const auto held = evaluate(hold, hold_current, work, outside);
    jet(held.jet, {Number{0}, Number{127, 25000}, Number{127, 25000, -1}}, zeros, outside);

    input.field("flow").field("uniform").set("outside", parse(R"({"kind":"periodic","period":4})"));
    const auto periodic = import(input, io::SourceId{15}, work);
    const auto periodic_current = import(input, io::SourceId{16}, work);
    const auto normalized = evaluate(periodic, periodic_current, work, outside);
    jet(normalized.jet, {Number{127, 12500}, Number{127, 25000}, Number{381, 25000, -1}}, zeros,
        outside);

    auto steady = root("aerodynamic_steady");
    steady.field("flow").set(
        "uniform",
        parse(
            R"({"kind":"pchip","times":[0,1],"values":[[1,2,3],[2,3,4]],"outside":{"kind":"reject"}})"));
    const auto si = materialize(steady, io::SourceId{17});
    core::ResourcePlan memory{budget};
    reject(cf::prepare_case_uniform_flow(si, stamp, {}, work, limits, memory),
           core::FailureCode::Unsupported);
    CHECK(memory.allocated_bytes() == 0);
    const Number zero{0}, one{1}, middle{1, 2};
    const cf::CaseUniformTimeWindow required{zero.view(), one.view()};
    const auto explicit_coverage =
        take(cf::prepare_case_uniform_flow(si, stamp, required, work, limits, memory));
    const auto current =
        take(cf::prepare_case_uniform_flow(si, stamp, required, work, limits, memory));
    const auto inside = evaluate(explicit_coverage, current, work, middle);
    jet(inside.jet, {Number{3, 2}, Number{5, 2}, Number{7, 2}}, {Number{1}, Number{1}, Number{1}},
        middle);
}

TEST_CASE("Typed case flow rejects incompatible coverage opaque variants quotas and cancellation",
          "[AE-07][IO4-01][N4-08][negative]") {
    a::UniformFlowWorkspace work;
    const auto constant = source(R"({"kind":"constant","value":[1,2,3]})", one_step);
    const auto si = materialize(constant, io::SourceId{10});
    auto short_limits = limits;
    short_limits.maximum_exact_numbers = 2;
    core::ResourcePlan refused{budget};
    reject(cf::prepare_case_uniform_flow(si, stamp, {}, work, short_limits, refused),
           core::FailureCode::ResourceExceeded);
    CHECK(refused.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 declares request_stop non-const.
    std::stop_source cancellation;
    cancellation.request_stop();
    reject(cf::prepare_case_uniform_flow(si, stamp, {}, work, limits, refused,
                                         cancellation.get_token()),
           core::FailureCode::Cancelled);
    CHECK(refused.allocated_bytes() == 0);

    auto gust = constant;
    gust.field("flow").set("gust",
                           parse(R"({"kind":"sine","amplitude":{"kind":"constant","value":[1]},
                           "direction":[0,0,1],"wavevector":[1,0,0],"angular_frequency":1,
                           "phase":0,"origin":[0,0,0],"time_origin":0})"));
    const auto gust_case = materialize(gust, io::SourceId{11});
    reject(cf::prepare_case_uniform_flow(gust_case, stamp, {}, work, limits, refused),
           core::FailureCode::Unsupported);
    CHECK(refused.allocated_bytes() == 0);

    const auto incomplete = source(
        R"({"kind":"pchip","times":[0,0.5],"values":[[1,2,3],[1,2,3]],"outside":{"kind":"reject"}})",
        one_step);
    const auto uncovered = materialize(incomplete, io::SourceId{12});
    reject(cf::prepare_case_uniform_flow(uncovered, stamp, {}, work, limits, refused),
           core::FailureCode::InvalidInput);
}

namespace {
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;

Rat two_power(std::int64_t exponent) {
    if (exponent >= 0)
        return Rat{Big{1} << static_cast<unsigned>(exponent)};
    return Rat{1} / (Big{1} << static_cast<unsigned>(-exponent));
}

/// Independent exact reconstruction of a binary word using its binary decomposition.
Rat rational_word(double word) {
    REQUIRE(std::isfinite(word));
    int exponent = 0;
    const auto fraction = std::frexp(word, &exponent);
    const auto significand = static_cast<std::int64_t>(std::ldexp(fraction, 53));
    return Rat{significand} * two_power(exponent - 53);
}

void within(const io::SiEnclosure& enclosure, const Rat& actual, const Rat& allowed,
            const Rat& reference_uncertainty = Rat{0}) {
    std::visit(
        [&](const auto& result) {
            REQUIRE(result.exponent2 > -4096);
            REQUIRE(result.exponent2 < 4096);
            const auto c = [&]() -> Rat {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(result)>,
                                             io::ScaledBinary64>)
                    return rational_word(result.center);
                else
                    return rational_word(result.high) + rational_word(result.low);
            }();
            const auto radius = rational_word(result.radius) * two_power(result.exponent2);
            CHECK(radius >= 0);
            CHECK(abs(c * two_power(result.exponent2) - actual) <= radius + reference_uncertainty);
            CHECK(radius <= allowed);
        },
        enclosure);
}
} // namespace

namespace {
JsonToken gust_case(std::string_view units, std::string_view amplitude) {
    auto input = source(R"({"kind":"constant","value":[1,2,3]})", one_step);
    input.field("units").set("length", string("m"));
    input.field("units").set("time", string("s"));
    input.field("units").set("angle", string(units));
    auto gust = parse(
        R"({"kind":"sine","amplitude":{"kind":"constant","value":[0]},"direction":[0,3,4],"wavevector":[90,0,0],"angular_frequency":90,"phase":45,"origin":[1,0,0],"time_origin":0.5})");
    gust.set("amplitude", parse(amplitude));
    input.field("flow").set("gust", std::move(gust));
    return input;
}

const cf::CaseKnownFlowLimits gust_limits{limits, 2048, 10000};

cf::CaseKnownFlowReady prepare_gust(const JsonToken& input, io::SourceId id,
                                    a::UniformFlowWorkspace& work, core::ResourcePlan& memory) {
    const auto si = materialize(input, id);
    return take(cf::prepare_case_known_flow(si, stamp, {}, work, gust_limits, memory));
}

cf::CaseKnownFlowProjection gust_projection(core::Precision precision, const Number& absolute) {
    const static Number unit{1}, zero{0};
    const cf::CaseFlowComponentBudget per_component{unit.view(), absolute.view(), zero.view()};
    std::array<cf::CaseFlowComponentBudget, 3> vector{per_component, per_component, per_component};
    std::array<cf::CaseFlowComponentBudget, 9> matrix{per_component, per_component, per_component,
                                                      per_component, per_component, per_component,
                                                      per_component, per_component, per_component};
    return {precision, vector, vector, matrix};
}

Rat oracle(const char* signed_digits) {
    Big power{1};
    for (int i = 0; i < 64; ++i)
        power *= 10;
    return Rat{Big{signed_digits}} / power;
}
} // namespace

TEST_CASE("Typed analytic gust encloses the SI velocity and both physical partials",
          "[AE-07][IO4-01][N4-02][N4-03][property][verification][negative]") {
    constexpr std::string_view polynomial =
        R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[2],[3]]})";
    const auto input = gust_case("deg", polynomial);
    a::UniformFlowWorkspace work;
    core::ResourcePlan memory{budget};
    const auto field = prepare_gust(input, io::SourceId{51}, work, memory);
    const auto current = prepare_gust(input, io::SourceId{52}, work, memory);
    const auto saved = take(cf::snapshot_case_known_flow(field, work, gust_limits, memory));
    const auto restored =
        take(cf::restore_case_known_flow(saved, current, work, gust_limits, memory));
    const Number t{1, 2}, x{1}, z{0}, x90{3, 2}, tolerance{1, 1000000000};
    const cf::CaseKnownFlowQuery query{{stamp, t.view()}, {x.view(), z.view(), z.view()}};
    const auto projection = gust_projection(core::Precision::TwoWord, tolerance);
    const auto result = take(cf::evaluate_case_known_flow(restored, current, query, projection,
                                                          work, gust_limits, memory));
    CHECK(take(result.stamp()) == stamp);
    exact(take(result.time()), t);
    exact(take(result.position(0)), x);
    exact(take(result.position(1)), z);
    exact(take(result.position(2)), z);
    CHECK(result.requested_side() == p::RequestedSide::Automatic);
    CHECK(result.interval_steps() <= gust_limits.interval_steps);

    // Independent closed values from sin(45 deg)=cos(45 deg)=sqrt(2)/2 and e=(0,3/5,4/5).
    const std::array<std::array<const char*, 3>, 2> digits{{
        {"34849242404917498012417731604201829824981554691457954768355137249",
         "-10597213363973567357618176679931359210599716990829841012608209295",
         "23325135425331422796833375197818641917726763869222373671198326937"},
        {"49798989873223330683223642138935773099975406255277273024473516332",
         "-14129617818631423143490902239908478947466289321106454683477612393",
         "31100180567108563729111166930424855890302351825629831561597769249"},
    }};
    const Rat oracle_error = two_power(-205);
    const Rat allowed = two_power(-98);
    within(take(result.velocity(0)), Rat{1}, Rat{0});
    within(take(result.time_partial(0)), Rat{0}, Rat{0});
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            if (i == 0 || j != 0)
                within(take(result.spatial_partial(i, j)), Rat{0}, Rat{0});
        }
    }
    for (std::size_t i = 0; i < digits.size(); ++i) {
        const auto component = i + 1;
        const auto& values = digits[i];
        within(take(result.velocity(component)), oracle(values[0]), allowed, oracle_error);
        within(take(result.time_partial(component)), oracle(values[1]), allowed, oracle_error);
        within(take(result.spatial_partial(component, 0)), oracle(values[2]), allowed,
               oracle_error);
    }

    const cf::CaseKnownFlowQuery cardinal{{stamp, t.view()}, {x90.view(), z.view(), z.view()}};
    const auto binary = gust_projection(core::Precision::Binary64, tolerance);
    const auto exact_cardinal = take(cf::evaluate_case_known_flow(
        restored, current, cardinal, binary, work, gust_limits, memory));
    CHECK(std::holds_alternative<io::ScaledBinary64>(take(exact_cardinal.velocity(1))));
    within(take(exact_cardinal.velocity(1)), Rat{41} / 10, two_power(-49));
    within(take(exact_cardinal.velocity(2)), Rat{29} / 5, two_power(-49));
    within(take(exact_cardinal.spatial_partial(1, 0)), Rat{0}, Rat{0});
    within(take(exact_cardinal.time_partial(1)), Rat{9} / 5, two_power(-49));

    const auto exact_demand = gust_projection(core::Precision::TwoWord, z);
    reject(cf::evaluate_case_known_flow(field, current, query, exact_demand, work, gust_limits,
                                        memory),
           core::FailureCode::PrecisionUnresolved);
    auto small = gust_limits;
    small.interval_steps = result.interval_steps() - 1;
    reject(cf::evaluate_case_known_flow(restored, current, query, projection, work, small, memory),
           core::FailureCode::ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 declares request_stop non-const.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    reject(cf::evaluate_case_known_flow(field, current, query, projection, work, gust_limits,
                                        memory, stop.get_token()),
           core::FailureCode::Cancelled);

    auto changed = input;
    changed.field("flow").field("gust").set("phase", parse("45.000000000000000000001"));
    const auto changed_phase = prepare_gust(changed, io::SourceId{53}, work, memory);
    core::ResourcePlan rejected{budget};
    reject(cf::restore_case_known_flow(saved, changed_phase, work, gust_limits, rejected),
           core::FailureCode::StaleGeneration);
    reject(cf::evaluate_case_known_flow(field, changed_phase, query, projection, work, gust_limits,
                                        rejected),
           core::FailureCode::StaleGeneration);
    CHECK(rejected.allocated_bytes() == 0);
    changed.field("flow").field("gust").set("direction", parse("[1,3,4]"));
    const auto skew = materialize(changed, io::SourceId{54});
    reject(cf::prepare_case_known_flow(skew, stamp, {}, work, gust_limits, rejected),
           core::FailureCode::InvalidInput);
}

TEST_CASE("A typed gust restores PCHIP amplitude and checks the current knot content",
          "[AE-07][IO4-01][N4-03][N4-09][verification][negative]") {
    auto input = gust_case(
        "rad",
        R"({"kind":"pchip","times":[0,0.5,1],"values":[[1],[2],[5]],"outside":{"kind":"reject"}})");
    input.field("analysis")
        .set(
            "intervals",
            parse(
                R"([{"start":0,"end":0.5,"step":0.5,"restart":"initial"},{"start":0.5,"end":1,"step":0.5,"restart":"event"}])"));
    input.field("flow").set("events", parse("[0.5]"));
    auto& gust = input.field("flow").field("gust");
    gust.set("wavevector", parse("[2,0,0]"));
    gust.set("angular_frequency", parse("3"));
    gust.set("phase", parse("0"));
    a::UniformFlowWorkspace work;
    core::ResourcePlan memory{budget};
    const auto field = prepare_gust(input, io::SourceId{55}, work, memory);
    const auto current = prepare_gust(input, io::SourceId{56}, work, memory);
    const auto saved = take(cf::snapshot_case_known_flow(field, work, gust_limits, memory));
    const auto restored =
        take(cf::restore_case_known_flow(saved, current, work, gust_limits, memory));
    const Number t{1, 2}, x{1}, z{0}, tolerance{1, 1000000000};
    const auto projection = gust_projection(core::Precision::TwoWord, tolerance);
    for (const auto side : {p::RequestedSide::Left, p::RequestedSide::Right}) {
        const cf::CaseKnownFlowQuery query{{stamp, t.view(), side}, {x.view(), z.view(), z.view()}};
        const auto result = take(cf::evaluate_case_known_flow(restored, current, query, projection,
                                                              work, gust_limits, memory));
        CHECK(result.requested_side() == side);
        const auto radius = two_power(-98);
        within(take(result.velocity(1)), Rat{2}, Rat{0});
        within(take(result.velocity(2)), Rat{3}, Rat{0});
        within(take(result.spatial_partial(1, 0)), Rat{12} / 5, radius);
        within(take(result.spatial_partial(2, 0)), Rat{16} / 5, radius);
        within(take(result.time_partial(1)), Rat{-18} / 5, radius);
        within(take(result.time_partial(2)), Rat{-24} / 5, radius);
    }
    auto changed = input;
    changed.field("flow").field("gust").field("amplitude").field("values").elements[2].elements[0] =
        parse("6");
    const auto changed_samples = prepare_gust(changed, io::SourceId{57}, work, memory);
    const cf::CaseKnownFlowQuery query{{stamp, t.view()}, {x.view(), z.view(), z.view()}};
    reject(cf::evaluate_case_known_flow(restored, changed_samples, query, projection, work,
                                        gust_limits, memory),
           core::FailureCode::StaleGeneration);
    reject(cf::restore_case_known_flow(saved, changed_samples, work, gust_limits, memory),
           core::FailureCode::StaleGeneration);
}

TEST_CASE("Known flow uses SI offsets and returns zero spatial partials without a gust",
          "[AE-07][IO4-01][N4-03][verification][negative]") {
    auto input = gust_case("deg", R"({"kind":"constant","value":[1]})");
    input.field("units").set("length", string("ft"));
    input.field("units").set("time", string("min"));
    auto& gust = input.field("flow").field("gust");
    gust.set("phase", parse("0"));
    gust.set("angular_frequency", parse("0"));
    gust.set("direction", parse("[0,0,-2]"));
    a::UniformFlowWorkspace work;
    core::ResourcePlan memory{budget};
    const auto field = prepare_gust(input, io::SourceId{58}, work, memory);
    const Number zero{0}, tolerance{1, 1000000000};
    const cf::CaseKnownFlowQuery query{{stamp, zero.view()},
                                       {zero.view(), zero.view(), zero.view()}};
    const auto projection = gust_projection(core::Precision::TwoWord, tolerance);
    // At x=0 the spatial phase is -90 degrees relative to the one-foot origin.
    const auto sample = take(
        cf::evaluate_case_known_flow(field, field, query, projection, work, gust_limits, memory));
    const auto radius = two_power(-98);
    within(take(sample.velocity(0)), Rat{127} / 25000, radius);
    within(take(sample.velocity(1)), Rat{127} / 12500, radius);
    within(take(sample.velocity(2)), Rat{127} / 6250, radius);
    within(take(sample.time_partial(2)), Rat{0}, Rat{0});
    // The mixed-unit reciprocal conversion can surround zero by a narrow interval.
    within(take(sample.spatial_partial(2, 0)), Rat{0}, radius);

    input.field("flow").set("gust", parse(R"({"kind":"none"})"));
    const auto ambient = prepare_gust(input, io::SourceId{59}, work, memory);
    const auto baseline = take(cf::evaluate_case_known_flow(ambient, ambient, query, projection,
                                                            work, gust_limits, memory));
    within(take(baseline.velocity(2)), Rat{381} / 25000, radius);
    for (std::size_t i = 0; i < 3; ++i) {
        within(take(baseline.time_partial(i)), Rat{0}, Rat{0});
        for (std::size_t j = 0; j < 3; ++j)
            within(take(baseline.spatial_partial(i, j)), Rat{0}, Rat{0});
    }
    reject(
        cf::evaluate_case_known_flow(ambient, field, query, projection, work, gust_limits, memory),
        core::FailureCode::StaleGeneration);
}

TEST_CASE("Typed SI prescribed wake owns its velocity independently of ambient flow",
          "[AE-07][UW-D02][UW-D08][IO4-01][N4-02][N4-03][N4-06][verification][negative]") {
    static_assert(!std::is_convertible_v<cf::CasePrescribedWakeReady, cf::CaseUniformFlowReady>);
    static_assert(
        !std::is_convertible_v<cf::ProjectedCasePrescribedWake, cf::ProjectedCaseUniformFlow>);
    a::UniformFlowWorkspace work;
    const auto case_with = [](std::string_view velocity) {
        auto input = source(R"({"kind":"constant","value":[100,200,300]})", two_step);
        input.field("flow").set("events", parse("[1]"));
        input.field("analysis")
            .field("wake")
            .set("convection", parse(std::string{R"({"kind":"prescribed","velocity":)"} +
                                     std::string{velocity} + "}"));
        return input;
    };
    const auto prepare_transport = [&](const io::SiCase& si) {
        core::ResourcePlan memory{budget};
        return take(cf::prepare_case_prescribed_wake(si, stamp, {}, work, limits, memory));
    };
    const auto import_transport = [&](const JsonToken& input, io::SourceId id) {
        const auto si = materialize(input, id);
        return prepare_transport(si);
    };
    const auto evaluate_transport = [&](const cf::CasePrescribedWakeReady& field,
                                        const cf::CasePrescribedWakeReady& current, const Number& t,
                                        p::RequestedSide side = p::RequestedSide::Automatic) {
        const cf::CaseUniformFlowQuery query{stamp, t.view(), side};
        const auto forecast =
            take(cf::forecast_case_prescribed_wake(field, current, query, limits.flow));
        auto admitted = limits.flow;
        admitted.arithmetic.steps = forecast.reserved_steps;
        admitted.metadata = numerics::linear::EstimateAdmission{forecast.metadata_estimate_bytes};
        core::ResourcePlan memory{forecast.payload_bytes};
        auto result =
            take(cf::evaluate_case_prescribed_wake(field, current, query, work, admitted, memory));
        CHECK(result.used_steps() <= forecast.reserved_steps);
        CHECK(result.requested_side() == side);
        CHECK(memory.allocated_bytes() == forecast.payload_bytes);
        return result;
    };

    const Number t{90}, mid{30}, zero{0}, one{1}, accuracy{1, 1000000000}, knot{60}, late{121};
    SECTION("Constant input uses transport units after the SI case owner is gone") {
        const auto input = case_with(R"({"kind":"constant","value":[2,-5,0]})");
        const auto field = import_transport(input, io::SourceId{301});
        auto changed_ambient = input;
        changed_ambient.field("flow").set("uniform",
                                          parse(R"({"kind":"constant","value":[9,8,7]})"));
        const auto current = import_transport(changed_ambient, io::SourceId{302});
        const auto value = evaluate_transport(field, current, t);
        jet(value.jet(), {Number{127, 12500}, Number{127, 5000, -1}, zero}, {zero, zero, zero}, t);
        const cf::CaseUniformFlowQuery query{stamp, t.view()};
        const auto forecast =
            take(cf::forecast_case_prescribed_wake(field, current, query, limits.flow));
        auto short_limit = limits.flow;
        REQUIRE(forecast.reserved_steps > 0);
        short_limit.arithmetic.steps = forecast.reserved_steps - 1;
        core::ResourcePlan refused{forecast.payload_bytes};
        reject(cf::evaluate_case_prescribed_wake(field, current, query, work, short_limit, refused),
               core::FailureCode::ResourceExceeded);
        CHECK(refused.allocated_bytes() == 0);
    }

    SECTION("Polynomial conversion uses the case start and end, origin and scale") {
        const auto input = case_with(
            R"({"kind":"polynomial","origin":1,"scale":2,"domain":[0,2],"coefficients":[[3,-2,5],[2,4,0]]})");
        const auto si_case = materialize(input, io::SourceId{303});
        const auto field = prepare_transport(si_case);
        const auto current = prepare_transport(si_case);
        const auto value = evaluate_transport(field, current, t);
        jet(value.jet(), {Number{889, 50000}, Number{127, 25000, -1}, Number{127, 5000}},
            {Number{127, 1500000}, Number{127, 750000}, zero}, t);
        core::ResourcePlan memory{budget};
        reject(cf::evaluate_case_prescribed_wake(field, current, {stamp, late.view()}, work,
                                                 limits.flow, memory),
               core::FailureCode::InvalidInput);

        auto short_window = input;
        short_window.field("analysis")
            .field("wake")
            .field("convection")
            .field("velocity")
            .set("domain", parse("[0,1]"));
        const auto si = materialize(short_window, io::SourceId{305});
        reject(cf::prepare_case_prescribed_wake(si, stamp, {}, work, limits, memory),
               core::FailureCode::InvalidInput);
    }

    SECTION("The PCHIP result restores by content and can be projected without reuse as ambient") {
        const auto input = case_with(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,-2,1],[1,-2,-1],[4,-2,-7]],"outside":{"kind":"reject"}})");
        const auto field = import_transport(input, io::SourceId{306});
        const auto current = import_transport(input, io::SourceId{307});
        core::ResourcePlan memory{budget};
        const auto saved =
            take(cf::snapshot_case_prescribed_wake(field, work, limits.flow, memory));
        const auto replay =
            take(cf::restore_case_prescribed_wake(saved, current, work, limits.flow, memory));
        auto value = evaluate_transport(replay, current, mid, p::RequestedSide::Right);
        jet(value.jet(), {Number{127, 80000}, Number{127, 12500, -1}, Number{381, 200000}},
            {Number{381, 4000000}, zero, Number{381, 2000000, -1}}, mid);
        for (const auto side :
             {p::RequestedSide::Automatic, p::RequestedSide::Left, p::RequestedSide::Right}) {
            const auto at_knot = evaluate_transport(replay, current, knot, side);
            jet(at_knot.jet(), {Number{127, 25000}, Number{127, 12500, -1}, Number{127, 25000, -1}},
                {Number{127, 1000000}, zero, Number{127, 500000, -1}}, knot);
        }

        const cf::CaseFlowComponentBudget strict{one.view(), zero.view(), zero.view()};
        const cf::CaseFlowComponentBudget allowed{one.view(), accuracy.view(), zero.view()};
        const std::array strict_components{strict, strict, strict};
        const std::array components{allowed, allowed, allowed};
        const auto project = [](cf::CasePrescribedWakeEvaluation& source,
                                const cf::CaseFlowProjectionRequest& request) {
            return cf::project_case_prescribed_wake(std::move(source), request);
        };
        reject(
            project(
                value,
                {core::Precision::Binary64, strict_components, strict_components, {4096, budget}}),
            core::FailureCode::PrecisionUnresolved);
        CHECK(value.jet().active());
        for (const auto representation : {core::Precision::Binary64, core::Precision::TwoWord}) {
            auto trial = evaluate_transport(replay, current, mid, p::RequestedSide::Right);
            const auto projected =
                take(project(trial, {representation, components, components, {4096, budget}}));
            CHECK_FALSE(trial.jet().active());
            CHECK(projected.jet().active());
            CHECK(projected.requested_side() == p::RequestedSide::Right);
            CHECK(projected.representation() == representation);
            CHECK(projected.used_steps() > 0);
            exact(take(projected.jet().time()), mid);
            const std::array<Rat, 3> velocity{Rat{127} / 80000, Rat{-127} / 12500,
                                              Rat{381} / 200000};
            const std::array<Rat, 3> derivative{Rat{381} / 4000000, Rat{0}, Rat{-381} / 2000000};
            for (std::size_t i = 0; i < velocity.size(); ++i) {
                within(take(projected.velocity(i)), velocity[i], Rat{1} / 1000000000);
                within(take(projected.time_partial(i)), derivative[i], Rat{1} / 1000000000);
                for (std::size_t j = 0; j < velocity.size(); ++j)
                    within(take(projected.spatial_partial(i, j)), Rat{0}, Rat{0});
            }
        }

        auto changed = input;
        changed.field("analysis")
            .field("wake")
            .field("convection")
            .field("velocity")
            .field("values")
            .elements[1]
            .elements[0] = parse("2");
        const auto stale = import_transport(changed, io::SourceId{308});
        const cf::CaseUniformFlowQuery query{stamp, mid.view()};
        core::ResourcePlan refused{budget};
        reject(cf::forecast_case_prescribed_wake(replay, stale, query, limits.flow),
               core::FailureCode::StaleGeneration);
        reject(cf::evaluate_case_prescribed_wake(replay, stale, query, work, limits.flow, refused),
               core::FailureCode::StaleGeneration);
        reject(cf::restore_case_prescribed_wake(saved, stale, work, limits.flow, refused),
               core::FailureCode::StaleGeneration);
        CHECK(refused.allocated_bytes() == 0);
    }

    SECTION("Outside policy distinguishes explicit endpoint hold from a signed periodic query") {
        auto input = case_with(
            R"({"kind":"pchip","times":[0,1,2,3,4],"values":[[0,1,-1],[0,1,-1],[2,1,-3],[0,1,-1],[0,1,-1]],"outside":{"kind":"hold"}})");
        input.field("flow").set("events", parse("[]"));
        input.field("analysis")
            .set("intervals", parse(R"([{"start":0,"end":4,"step":1,"restart":"initial"}])"));
        const Number before{120, 1, -1};
        const auto hold_case = materialize(input, io::SourceId{312});
        const auto held = prepare_transport(hold_case);
        const auto held_current = prepare_transport(hold_case);
        const auto held_value = evaluate_transport(held, held_current, before);
        jet(held_value.jet(), {zero, Number{127, 25000}, Number{127, 25000, -1}},
            {zero, zero, zero}, before);

        input.field("analysis")
            .field("wake")
            .field("convection")
            .field("velocity")
            .set("outside", parse(R"({"kind":"periodic","period":4})"));
        const auto periodic_case = materialize(input, io::SourceId{314});
        const auto cyclic = prepare_transport(periodic_case);
        const auto cyclic_current = prepare_transport(periodic_case);
        const auto wrapped = evaluate_transport(cyclic, cyclic_current, before);
        jet(wrapped.jet(), {Number{127, 12500}, Number{127, 25000}, Number{381, 25000, -1}},
            {zero, zero, zero}, before);
    }

    SECTION("Free and absent wake modes are rejected before allocating a curve") {
        const auto free = materialize(source(R"({"kind":"constant","value":[2,-5,0]})", one_step),
                                      io::SourceId{309});
        const auto none = materialize(root("structural_static"), io::SourceId{310});
        const auto valid =
            materialize(case_with(R"({"kind":"constant","value":[2,-5,0]})"), io::SourceId{311});
        core::ResourcePlan memory{budget};
        for (const auto* si : std::array{&free, &none})
            reject(cf::prepare_case_prescribed_wake(*si, stamp, {}, work, limits, memory),
                   core::FailureCode::Unsupported);
        auto restricted = limits;
        restricted.maximum_exact_numbers = 2;
        reject(cf::prepare_case_prescribed_wake(valid, stamp, {}, work, restricted, memory),
               core::FailureCode::ResourceExceeded);
        // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
        std::stop_source cancel;
        REQUIRE(cancel.request_stop());
        reject(cf::prepare_case_prescribed_wake(valid, stamp, {}, work, limits, memory,
                                                cancel.get_token()),
               core::FailureCode::Cancelled);
        CHECK(memory.allocated_bytes() == 0);
    }
}

TEST_CASE("A typed transient clock binds both flow fields to one exact stage",
          "[AE-07][UW-D02][IO4-01][N4-03][N4-06][verification][negative]") {
    a::UniformFlowWorkspace work;
    const cf::CaseFlowClockLimits clock_limits{limits, 2, 4};
    const auto transport = [](JsonToken input, std::string_view velocity) {
        input.field("analysis")
            .field("wake")
            .set("convection", parse(std::string{R"({"kind":"prescribed","velocity":)"} +
                                     std::string{velocity} + "}"));
        return input;
    };
    const Number one{1}, zero{0}, absolute{1, 1000000000};
    const cf::CaseFlowComponentBudget allowed{one.view(), absolute.view(), zero.view()};
    const std::array components{allowed, allowed, allowed};
    const cf::CaseFlowProjectionRequest projection{
        core::Precision::TwoWord, components, components, {4096, budget}};
    const auto prepare_clock = [&](const JsonToken& input, io::SourceId id) {
        const auto si = materialize(input, id);
        core::ResourcePlan memory{budget};
        return take(cf::prepare_case_clocked_flows(si, stamp, work, clock_limits, memory));
    };

    SECTION("Both projections share a non-dyadic instant and each endpoint side keeps its own h") {
        auto input = transport(
            source(
                R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[1,0,0],[1000,0,0]]})",
                R"([{"start":0.001,"end":0.002,"step":0.0005,"restart":"initial"},{"start":0.002,"end":0.003,"step":0.00025,"restart":"event"}])"),
            R"({"kind":"pchip","times":[0.001,0.002,0.003],"values":[[2,0,0],[4,0,0],[6,0,0]],"outside":{"kind":"reject"}})");
        input.field("flow").set("events", parse("[0.002]"));
        const auto field = prepare_clock(input, io::SourceId{401});
        const auto current = prepare_clock(input, io::SourceId{402});
        CHECK(field.interval_count() == 2);
        CHECK(take(field.interval_steps(0)) == 2);
        CHECK(take(field.interval_steps(1)) == 4);
        reject(field.interval_steps(2), core::FailureCode::InvalidInput);
        const auto run = [&](cf::CaseFlowClockQuery address) {
            core::ResourcePlan memory{budget};
            return take(cf::evaluate_case_clocked_flows(field, current, address, projection,
                                                        projection, work, clock_limits, memory));
        };
        const auto compare = [&](const cf::CaseClockedFlowSample& result, const Number& time,
                                 const Number& step, p::RequestedSide side, const Rat& ambient,
                                 const Rat& wake) {
            exact(take(result.time()), time);
            exact(take(result.step_seconds()), step);
            exact(take(result.ambient().source().jet.time()), time);
            exact(take(result.prescribed().jet().time()), time);
            CHECK(result.ambient().source().requested_side == side);
            CHECK(result.prescribed().requested_side() == side);
            CHECK(result.address().stamp == stamp);
            CHECK(result.address().side == side);
            CHECK(result.source().source.source == io::SourceId{401});
            within(take(result.ambient().velocity(0)), ambient, Rat{1} / 1000000000);
            within(take(result.prescribed().velocity(0)), wake, Rat{1} / 1000000000);
            within(take(result.ambient().time_partial(0)), Rat{127} / 1500, Rat{1} / 1000000000);
            within(take(result.prescribed().time_partial(0)), Rat{127} / 750, Rat{1} / 1000000000);
        };
        const auto middle = run({stamp, 0, 1, p::RequestedSide::Automatic});
        compare(middle, Number{9, 100}, Number{3, 100}, p::RequestedSide::Automatic,
                Rat{127} / 10000, Rat{381} / 25000);
        CHECK(middle.restart() == cf::CaseClockRestart::ItemInitial);
        const auto left = run({stamp, 0, 2, p::RequestedSide::Left});
        const auto right = run({stamp, 1, 0, p::RequestedSide::Right});
        compare(left, Number{3, 25}, Number{3, 100}, p::RequestedSide::Left, Rat{381} / 25000,
                Rat{127} / 6250);
        compare(right, Number{3, 25}, Number{3, 200}, p::RequestedSide::Right, Rat{381} / 25000,
                Rat{127} / 6250);
        CHECK(left.restart() == cf::CaseClockRestart::ItemInitial);
        CHECK(right.restart() == cf::CaseClockRestart::ItemEvent);
        CHECK(right.address().interval == 1);
        CHECK(right.address().step == 0);

        core::ResourcePlan refused{budget};
        const auto refuse = [&](const cf::CaseClockedFlows& supplied,
                                cf::CaseFlowClockQuery query) {
            return cf::evaluate_case_clocked_flows(field, supplied, query, projection, projection,
                                                   work, clock_limits, refused);
        };
        const auto next = p::PchipInputStamp{stamp.issuer, p::PchipInputGeneration{9}};
        reject(refuse(current, {next, 0, 1, p::RequestedSide::Left}),
               core::FailureCode::StaleGeneration);
        for (const auto& bad :
             {cf::CaseFlowClockQuery{stamp, 2, 0, p::RequestedSide::Right},
              cf::CaseFlowClockQuery{stamp, 0, 3, p::RequestedSide::Left},
              cf::CaseFlowClockQuery{stamp, 0, 2, p::RequestedSide::Right},
              cf::CaseFlowClockQuery{stamp, 1, 0, p::RequestedSide::Left},
              cf::CaseFlowClockQuery{stamp, 1, 0, p::RequestedSide::Automatic},
              cf::CaseFlowClockQuery{stamp, 0, 1, static_cast<p::RequestedSide>(255)}})
            reject(refuse(current, bad), core::FailureCode::InvalidInput);
        auto changed_grid = input;
        changed_grid.field("analysis").field("intervals").elements[1].set("step", parse("0.0005"));
        const auto other_grid = prepare_clock(changed_grid, io::SourceId{403});
        reject(refuse(other_grid, {stamp, 0, 1, p::RequestedSide::Left}),
               core::FailureCode::StaleGeneration);
        CHECK(refused.allocated_bytes() == 0);

        auto altered_wake = input;
        altered_wake.field("analysis")
            .field("wake")
            .field("convection")
            .field("velocity")
            .field("values")
            .elements[1]
            .elements[0] = parse("5");
        const auto other_wake = prepare_clock(altered_wake, io::SourceId{404});
        reject(refuse(other_wake, {stamp, 0, 1, p::RequestedSide::Left}),
               core::FailureCode::StaleGeneration);
        const cf::CaseFlowComponentBudget exact_bound{one.view(), zero.view(), zero.view()};
        const std::array exact_bounds{exact_bound, exact_bound, exact_bound};
        const cf::CaseFlowProjectionRequest exact_only{
            core::Precision::Binary64, exact_bounds, exact_bounds, {4096, budget}};
        reject(cf::evaluate_case_clocked_flows(field, current,
                                               {stamp, 0, 1, p::RequestedSide::Left}, projection,
                                               exact_only, work, clock_limits, refused),
               core::FailureCode::PrecisionUnresolved);
        const auto retry = run({stamp, 0, 1, p::RequestedSide::Left});
        CHECK(retry.prescribed().jet().active());
    }

    SECTION("Exact time identity distinguishes schedules with the same binary64 time words") {
        constexpr auto velocity = R"({"kind":"constant","value":[2,0,0]})";
        constexpr auto ambient = R"({"kind":"constant","value":[1,0,0]})";
        static_assert(.1 == .100000000000000001);
        const auto first = transport(
            source(ambient, R"([{"start":0.1,"end":0.2,"step":0.1,"restart":"initial"}])"),
            velocity);
        const auto second = transport(
            source(
                ambient,
                R"([{"start":0.100000000000000001,"end":0.200000000000000001,"step":0.1,"restart":"initial"}])"),
            velocity);
        const auto old = prepare_clock(first, io::SourceId{405});
        const auto different = prepare_clock(second, io::SourceId{406});
        core::ResourcePlan memory{budget};
        reject(cf::evaluate_case_clocked_flows(old, different,
                                               {stamp, 0, 0, p::RequestedSide::Right}, projection,
                                               projection, work, clock_limits, memory),
               core::FailureCode::StaleGeneration);
        CHECK(memory.allocated_bytes() == 0);
        const auto si = materialize(first, io::SourceId{407});
        auto narrow = clock_limits;
        narrow.maximum_intervals = 0;
        reject(cf::prepare_case_clocked_flows(si, stamp, work, narrow, memory),
               core::FailureCode::ResourceExceeded);
        CHECK(memory.allocated_bytes() == 0);
        narrow = clock_limits;
        narrow.maximum_steps_per_interval = 0;
        reject(cf::prepare_case_clocked_flows(si, stamp, work, narrow, memory),
               core::FailureCode::ResourceExceeded);
    }
}

TEST_CASE("Typed SI uniform jet projects all components and retains the exact source",
          "[AE-07][N4-06][N4-08][N4-09][property][negative]") {
    const auto input = source(
        R"({"kind":"pchip","times":[0,1,2],"values":[[0,-2,1],[1,-2,-1],[4,-2,-7]],"outside":{"kind":"reject"}})",
        R"([{"start":0,"end":2,"step":1,"restart":"initial"}])");
    a::UniformFlowWorkspace work;
    const auto field = import(input, io::SourceId{31}, work);
    const auto current = import(input, io::SourceId{32}, work);
    const Number time{30}, zero{0}, scale{2};
    constexpr std::array<std::uint32_t, 1> one_word{1};
    constexpr std::array<std::uint32_t, 3> den90{0, 0, 1U << 26U};
    constexpr std::array<std::uint32_t, 2> den40{0, 1U << 8U};
    const cf::CaseFlowComponentBudget strict{scale.view(), zero.view(),
                                             p::RationalView{1, one_word, den90}};
    const cf::CaseFlowComponentBudget wide{scale.view(), p::RationalView{1, one_word, den40},
                                           zero.view()};
    const std::array host_requirements{strict, strict, strict};
    auto exact_jet = evaluate(field, current, work, time, p::RequestedSide::Right);
    cf::CaseFlowProjectionRequest request{
        core::Precision::Binary64, host_requirements, host_requirements, {4096, budget}};
    // The projection consumes its input only on success.
    const auto project = [&](cf::CaseUniformFlowEvaluation& evaluation,
                             const std::stop_token& stop = {}) {
        return cf::project_case_uniform_flow(std::move(evaluation), request, stop);
    };
    reject(project(exact_jet), core::FailureCode::PrecisionUnresolved);
    CHECK(exact_jet.jet.active());
    const std::array wide_values{wide, wide, wide};
    request.velocity = wide_values;
    request.time_partial = wide_values;
    request.time_partial[2] = {scale.view(), zero.view(), zero.view()};
    reject(project(exact_jet), core::FailureCode::PrecisionUnresolved);
    CHECK(exact_jet.jet.active());
    request.time_partial = wide_values;
    auto first = evaluate(field, current, work, time, p::RequestedSide::Right);
    const auto binary = take(project(first));
    CHECK(binary.representation() == core::Precision::Binary64);

    request.representation = core::Precision::TwoWord;
    request.velocity = host_requirements;
    request.time_partial = host_requirements;
    const auto projected = take(project(exact_jet));
    CHECK_FALSE(exact_jet.jet.active());
    CHECK(projected.representation() == core::Precision::TwoWord);
    CHECK(projected.source().requested_side == p::RequestedSide::Right);
    CHECK(take(projected.source().jet.stamp()) == stamp);
    exact(take(projected.source().jet.time()), time);

    const std::array<Rat, 3> truth_u{Rat{127} / 80000, Rat{-127} / 12500, Rat{381} / 200000};
    const std::array<Rat, 3> truth_dt{Rat{381} / 4000000, Rat{0}, Rat{-381} / 2000000};
    for (std::size_t c = 0; c < 3; ++c) {
        within(take(projected.velocity(c)), truth_u[c], abs(truth_u[c]) * two_power(-90));
        within(take(projected.time_partial(c)), truth_dt[c], abs(truth_dt[c]) * two_power(-90));
        within(take(binary.velocity(c)), truth_u[c], two_power(-39));
        within(take(binary.time_partial(c)), truth_dt[c], two_power(-39));
        for (std::size_t x = 0; x < 3; ++x)
            within(take(projected.spatial_partial(c, x)), Rat{0}, Rat{0});
    }
    reject(projected.spatial_partial(0, 3), core::FailureCode::InvalidInput);

    auto trial = evaluate(field, current, work, time);
    REQUIRE(projected.used_steps() > 1);
    request.limits.steps = projected.used_steps() - 1;
    reject(project(trial), core::FailureCode::ResourceExceeded);
    CHECK(trial.jet.active());
    request.limits.steps = projected.used_steps();
    const auto at_limit = take(project(trial));
    CHECK(at_limit.used_steps() == projected.used_steps());

    auto refusals = evaluate(field, current, work, time);
    request.velocity[0].scale = zero.view();
    reject(project(refusals), core::FailureCode::InvalidInput);
    request.velocity = host_requirements;
    request.representation = static_cast<core::Precision>(255);
    reject(project(refusals), core::FailureCode::PrecisionUnsupported);
    request.representation = core::Precision::TwoWord;
    request.limits.input_integer_bits = 64;
    reject(project(refusals), core::FailureCode::ResourceExceeded);
    request.limits.input_integer_bits = 4096;
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 declares request_stop non-const.
    std::stop_source cancellation;
    cancellation.request_stop();
    reject(project(refusals, cancellation.get_token()), core::FailureCode::Cancelled);
    CHECK(refusals.jet.active());
}

TEST_CASE("Scaled rational projection preserves distant exact words and their exponent",
          "[AE-07][N4-06][N4-08][property][negative]") {
    using io::detail::SiRationalWordsView;
    constexpr std::array<std::uint32_t, 1> one{1};
    const SiRationalWordsView zero{0, {}, one}, unit{1, one, one};
    const io::detail::SiRationalProjectionBudget exact{unit, zero, zero};
    const io::detail::SiRationalProjectionLimits work{4096, budget};
    const auto check_two = [&](SiRationalWordsView source, const Rat& expected, double high,
                               double low) {
        const auto result = take(io::detail::project_rational_words(
            source, exact, core::Precision::TwoWord, work, {}, {}));
        const auto* words = std::get_if<io::ScaledTwoWord>(&result.value);
        if (!words)
            FAIL("Requested two-word output");
        CHECK(words->high == high);
        CHECK(words->low == low);
        CHECK(words->radius == 0.0);
        within(result.value, expected, Rat{0});
    };
    std::array<std::uint32_t, 29> numerator{}, denominator{};
    numerator.front() = 1;
    numerator.back() = denominator.back() = 1U << 4U;
    check_two({1, numerator, denominator}, 1 + two_power(-900), 1.0, std::ldexp(1.0, -900));

    std::array<std::uint32_t, 34> at_two{}, grid{};
    at_two.fill(std::numeric_limits<std::uint32_t>::max());
    at_two.back() = (1U << 19U) - 1;
    grid.back() = 1U << 18U;
    check_two({1, at_two, grid}, 2 - two_power(-1074), 2.0,
              -std::numeric_limits<double>::denorm_min());

    std::array<std::uint32_t, 35> boundary{};
    boundary.back() = 1U << 12U;
    for (const auto& [source, power] : {std::pair{SiRationalWordsView{1, boundary, one}, 1100},
                                        std::pair{SiRationalWordsView{-1, one, boundary}, -1100}}) {
        const auto scaled = take(io::detail::project_rational_words(
            source, exact, core::Precision::Binary64, work, {}, {}));
        within(scaled.value, source.sign * two_power(power), Rat{0});
    }
    // 1+2^-53 is halfway between consecutive binary64 words. Exactly 2^-53 is admissible.
    constexpr std::array<std::uint32_t, 2> half_n{1, 1U << 21U}, half_d{0, 1U << 21U};
    constexpr std::array<std::uint32_t, 2> below_n{0xffffffffU, (1U << 21U) - 1};
    constexpr std::array<std::uint32_t, 4> below_d{0, 0, 0, 1U << 10U};
    const SiRationalWordsView halfway{1, half_n, half_d};
    const io::detail::SiRationalProjectionBudget enough{unit, {1, one, half_d}, zero};
    const auto tie = take(io::detail::project_rational_words(
        halfway, enough, core::Precision::Binary64, work, {}, {}));
    const auto* binary = std::get_if<io::ScaledBinary64>(&tie.value);
    if (!binary)
        FAIL("Requested one-word output");
    CHECK(binary->center == 1.0);
    within(tie.value, 1 + two_power(-53), two_power(-53));
    const io::detail::SiRationalProjectionBudget too_small{unit, {1, below_n, below_d}, zero};
    reject(io::detail::project_rational_words(halfway, too_small, core::Precision::Binary64, work,
                                              {}, {}),
           core::FailureCode::PrecisionUnresolved);

    const std::array<std::uint32_t, 1> three{3};
    reject(io::detail::project_rational_words({1, one, three}, exact, core::Precision::TwoWord,
                                              work, {}, {}),
           core::FailureCode::PrecisionUnresolved);
}

namespace {
cf::NonpenetrationVector coordinates(double x, double y, double z) noexcept {
    return {a::SingleWord{x}, a::SingleWord{y}, a::SingleWord{z}};
}

template <class T>
T take_global(a::GlobalResult<T> result) {
    if (!result)
        FAIL("Unexpected global qualification failure "
             << static_cast<int>(result.error().origin.cause.code));
    return std::move(*result);
}

struct ResidualFixture {
    std::array<a::PointId, 6> id{a::PointId{1}, a::PointId{2}, a::PointId{3},
                                 a::PointId{4}, a::PointId{5}, a::PointId{6}};
    std::array<a::PointBinding, 6> bindings{test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding()};
    std::array<double, 18> points{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 2, 1, 0, 2, 0, 0};
    std::array<a::BandId, 2> bands{a::BandId{50}, a::BandId{51}};
    std::array<a::FaceInput, 2> faces{a::FaceInput{a::PhysicalFaceId{100}, bands[0]},
                                      a::FaceInput{a::PhysicalFaceId{101}, bands[1]}};
    std::array<a::FaceRoleInput, 2> roles{a::FaceRoleInput{faces[0].id, a::AeroFaceRole::Bound},
                                          a::FaceRoleInput{faces[1].id, a::AeroFaceRole::Wake}};
    std::array<a::EdgeInput, 9> edges{
        a::EdgeInput{a::EdgeId{10}, {id[0], id[1]}}, a::EdgeInput{a::EdgeId{11}, {id[1], id[2]}},
        a::EdgeInput{a::EdgeId{12}, {id[2], id[3]}}, a::EdgeInput{a::EdgeId{13}, {id[3], id[0]}},
        a::EdgeInput{a::EdgeId{14}, {id[0], id[2]}}, a::EdgeInput{a::EdgeId{15}, {id[2], id[4]}},
        a::EdgeInput{a::EdgeId{16}, {id[4], id[5]}}, a::EdgeInput{a::EdgeId{17}, {id[5], id[3]}},
        a::EdgeInput{a::EdgeId{18}, {id[3], id[4]}}};
    std::array<a::EdgeBindingInput, 9> edge_owners{
        a::EdgeBindingInput{edges[0].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[1].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[2].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[3].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[4].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[5].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[6].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[7].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[8].id, a::CanonicalEdge{}}};
    std::array<a::TrailingEdgeInput, 1> trailing{{{edges[2].id, {id[3], id[2]}}}};
    std::array<a::QuadInput, 2> quads{
        a::QuadInput{faces[0].id,
                     {id[0], id[1], id[2], id[3]},
                     {edges[0].id, edges[1].id, edges[2].id, edges[3].id},
                     edges[4].id,
                     {a::FacetId{200}, a::FacetId{201}}},
        a::QuadInput{faces[1].id,
                     {id[3], id[2], id[4], id[5]},
                     {edges[2].id, edges[5].id, edges[6].id, edges[7].id},
                     edges[8].id,
                     {a::FacetId{202}, a::FacetId{203}}}};
    void reverse_wake() noexcept {
        edges[8].points = {id[2], id[5]};
        trailing[0].points = {id[2], id[3]};
        quads[1].points = {id[2], id[3], id[5], id[4]};
        quads[1].edges = {edges[2].id, edges[7].id, edges[6].id, edges[5].id};
    }
    [[nodiscard]] a::GeometryInput input() const noexcept {
        const a::SinglePositions xyz{points};
        return {id,
                bindings,
                xyz,
                xyz,
                edges,
                faces,
                bands,
                {},
                quads,
                a::SingleWord{.125},
                test_geometry::context()};
    }
    [[nodiscard]] a::GlobalTopologyInput topology() const noexcept {
        return {roles, edge_owners, trailing};
    }
};

void residual_expected(const a::ScaledEnclosure& value, const Rat& oracle,
                       const Rat& comparison_uncertainty = Rat{Big{1},
                                                               Big{"1000000000000000000000000"}}) {
    const auto center =
        std::visit(
            [](const auto& word) -> Rat {
                if constexpr (std::is_same_v<std::decay_t<decltype(word)>, a::SingleWord>)
                    return rational_word(word.value);
                else
                    return rational_word(word.high) + rational_word(word.low);
            },
            value.center()) *
        two_power(value.exponent());
    const auto radius = rational_word(value.radius()) * two_power(value.exponent());
    CHECK(abs(center - oracle) <= radius + comparison_uncertainty);
    CHECK(radius < comparison_uncertainty);
}

namespace ls = numerics::linear;
ls::Context circulation_context(const a::GeometryContext& geometry) {
    const auto trial = geometry.diagnostics.trial;
    if (!trial)
        FAIL("The circulation fixture requires a trial");
    return {{geometry.fingerprint, geometry.model, geometry.current, core::ChartGeneration{1},
             core::SparsityGeneration{2}, core::NumericalGeneration{3}, core::ScaleGeneration{4},
             geometry.precision, core::Precision::Binary64, core::WakePolicyGeneration{5},
             core::TimeIntervalGeneration{6}, core::BackendGeneration{7}},
            *trial,
            core::ScopeId{8},
            core::ScopeId{9},
            geometry.time};
}

cf::DiscreteCirculationRequest circulation_request(const ls::Scales& scales, core::Index order,
                                                   ls::Norm norm) {
    return {scales,
            {order, 1, 0, 5, ls::EstimateAdmission{1U << 20U}},
            {norm, 1e-28, 0., 1e-26, core::ObservationId{20}, core::ObservationId{21}},
            core::RhsGeneration{3},
            core::NumericalGeneration{4},
            1e-14,
            1e-14,
            1e-13,
            0.};
}

Rat solution_word(const ls::Solution<double>& value, std::size_t i) {
    return std::visit(
        [&](const auto& words) -> Rat {
            if constexpr (std::is_same_v<std::decay_t<decltype(words)>, numerics::Binary64<double>>)
                return rational_word(words.values.values()[i]);
            else
                return rational_word(words.hi.values()[i]) + rational_word(words.lo.values()[i]);
        },
        value.values().storage());
}

void residual_small(const a::ScaledEnclosure& value, const Rat& limit) {
    const auto center = std::visit(
        [](const auto& word) -> Rat {
            if constexpr (std::is_same_v<std::decay_t<decltype(word)>, a::SingleWord>)
                return rational_word(word.value);
            else
                return rational_word(word.high) + rational_word(word.low);
        },
        value.center());
    const auto scale = two_power(value.exponent());
    CHECK((abs(center) + rational_word(value.radius())) * scale <= limit);
}

Rat sum_words(const a::ScalarWords& words) {
    return std::visit(
        [](const auto& word) -> Rat {
            if constexpr (std::is_same_v<std::decay_t<decltype(word)>, a::SingleWord>)
                return rational_word(word.value);
            else
                return rational_word(word.high) + rational_word(word.low);
        },
        words);
}

struct WorkEnclosure {
    Rat value{};
    Rat radius{};

    void add(const a::ScaledEnclosure& dual, const a::ScalarWords& direction) {
        const auto scale = two_power(dual.exponent());
        const auto displacement = sum_words(direction);
        value += sum_words(dual.center()) * scale * displacement;
        radius += rational_word(dual.radius()) * scale * abs(displacement);
    }

    void add(const cf::CirculationVectorCovector& dual, const cf::NonpenetrationVector& direction) {
        for (std::size_t axis = 0; axis < dual.size(); ++axis)
            add(dual[axis], direction[axis]);
    }
};

WorkEnclosure pullback_work(const cf::DiscreteCirculationPullback& value,
                            const cf::DiscreteCirculationDirection& direction) {
    REQUIRE(value.point_count() == direction.points.size());
    REQUIRE(value.row_count() == direction.rows.size());
    WorkEnclosure work;
    work.add(value.uniform_flow(), direction.uniform_flow);
    for (std::size_t i = 0; i < value.point_count(); ++i) {
        const auto point = take(value.point(i));
        CHECK(point.point == direction.points[i]);
        work.add(point.position, direction.displacements[i]);
    }
    for (std::size_t i = 0; i < value.row_count(); ++i) {
        const auto row = take(value.row(i));
        CHECK(row.face == direction.rows[i].face);
        const auto& input = direction.rows[i];
        REQUIRE(row.free_query.has_value() == input.free_query.has_value());
        if (row.free_query && input.free_query)
            work.add(*row.free_query, *input.free_query);
        work.add(row.tangents[0], input.tangents[0]);
        work.add(row.tangents[1], input.tangents[1]);
        work.add(row.velocity, input.velocity);
    }
    return work;
}

} // namespace

TEST_CASE("Discrete nonpenetration joins the original residual and implicit circulation action",
          "[AE-07][UW-D01][UW-D03][UW-D06][N4-06][verification]") {
    const auto source_case = source(R"({"kind":"constant","value":[0,0,100]})", one_step);
    a::UniformFlowWorkspace flow_workspace;
    const auto provider = import(source_case, io::SourceId{101}, flow_workspace);
    const auto current = import(source_case, io::SourceId{102}, flow_workspace);
    const Number one{1}, zero{0}, absolute{1, 1000000000}, time{5, 4}, wrong_time{127, 100};
    const cf::CaseFlowComponentBudget budget_component{one.view(), absolute.view(), zero.view()};
    const std::array budgets{budget_component, budget_component, budget_component};
    const cf::CaseFlowProjectionRequest projection{
        core::Precision::TwoWord, budgets, budgets, {4096, budget}};
    const auto projected_at = [&](const Number& instant) {
        auto value = evaluate(provider, current, flow_workspace, instant);
        return take(cf::project_case_uniform_flow(std::move(value), projection));
    };
    const a::GlobalLimits cap{std::numeric_limits<std::uint64_t>::max() / 4,
                              numerics::linear::EstimateAdmission{4096}};
    const cf::DiscreteNonpenetrationLimits residual_limits{
        {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10},
        {test_geometry::limits(),
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(6000000000)};
    const Rat denominator{Big{"10000000000000000000000000000000000000000"}};
    const Rat residual = -Rat{Big{"5671057487729860677895669210567132495399"}} / denominator;
    const Rat action = Rat{Big{"3730450057081755951141985320248458351172"}} / denominator;
    for (const bool reversed : {false, true}) {
        CAPTURE(reversed);
        ResidualFixture fixture;
        if (reversed)
            fixture.reverse_wake();
        core::ResourcePlan memory{128ULL * 1024 * 1024};
        auto geometry_work = test_geometry::workspace(memory);
        const auto geometry = test_geometry::take(a::prepare_owned_geometry(
            fixture.input(), test_geometry::limits(), *geometry_work, memory));
        auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
        const auto topology = take_global(
            a::build_global_topology(geometry, fixture.topology(), cap, *global, memory));
        const auto contact =
            take_global(a::check_global_geometry(topology, geometry, cap, *global, memory));
        const auto query = test_geometry::take(
            a::resolve_query(geometry, {a::QueryId{29}, a::FreePoint{a::FreeSingle{{.75, .5, 0}}}},
                             geometry.context(), test_geometry::limits(), *geometry_work, memory));
        const auto wake_geometry = take_global(contact.uvlm_quad(a::FacetId{202}, 1000));
        const auto& attachment = wake_geometry.attachment;
        if (!attachment)
            FAIL("The wake has no directed attachment");
        CHECK(static_cast<int>(attachment->coefficient) == (reversed ? -1 : 1));
        std::array collocation{cf::DiscreteCollocation{a::FacetId{200},
                                                       query,
                                                       {coordinates(1, 0, 0), coordinates(0, 2, 0)},
                                                       coordinates(.25, .25, .5),
                                                       a::DoubleWord{.75, 0x1p-60},
                                                       a::SingleWord{-.5}}};
        // The positive front shares the bound type even when an independent history ID
        // has the same numeric face value. The negative front must skip wake ID zero.
        const auto bound_id = a::P1StrengthId::bound(fixture.faces[0].id);
        const auto historical_id = a::P1StrengthId{reversed ? 0 : fixture.faces[0].id.value()};
        std::array wake{
            cf::DiscreteWakePanel{a::FacetId{202}, cf::AttachedWakeFront{}, historical_id}};
        const std::array history{
            cf::FrozenWakeTrace{historical_id, a::SingleWord{reversed ? .25 : -.25}}};
        cf::DiscreteNonpenetrationInput view{contact,
                                             collocation,
                                             wake,
                                             {core::AcceptedStateId{2}, history},
                                             {7, 11},
                                             stamp,
                                             p::RequestedSide::Automatic};
        a::SharedP1Workspace shared;
        auto local = take(a::LocalPacketWorkspace::create(memory));
        auto kernel = take(a::P1Workspace::create(memory));
        const auto run = [&](cf::ProjectedCaseUniformFlow& flow,
                             const cf::DiscreteNonpenetrationLimits& selected,
                             core::ResourcePlan& plan) {
            return cf::evaluate_discrete_nonpenetration(view, std::move(flow), selected, shared,
                                                        *local, *kernel, plan);
        };
        auto old_time = projected_at(wrong_time);
        reject(run(old_time, residual_limits, memory), core::FailureCode::StaleGeneration);
        CHECK(old_time.source().jet.active());
        auto flow = projected_at(time);
        view.history.parent = core::AcceptedStateId{99};
        reject(run(flow, residual_limits, memory), core::FailureCode::StaleGeneration);
        view.history.parent = core::AcceptedStateId{2};
        view.flow_side = p::RequestedSide::Left;
        reject(run(flow, residual_limits, memory), core::FailureCode::StaleGeneration);
        view.flow_side = p::RequestedSide::Automatic;
        core::ResourcePlan rejected_memory{0};
        reject(run(flow, residual_limits, rejected_memory), core::FailureCode::ResourceExceeded);
        CHECK(rejected_memory.allocated_bytes() == 0);
        wake[0].leading = historical_id;
        reject(run(flow, residual_limits, memory), core::FailureCode::InvalidInput);
        wake[0].leading = cf::AttachedWakeFront{};
        wake[0].trailing = bound_id;
        reject(run(flow, residual_limits, memory), core::FailureCode::InvalidInput);
        wake[0].trailing = historical_id;
        const std::array invalid_history{cf::FrozenWakeTrace{bound_id, history[0].strength}};
        view.history.columns = invalid_history;
        reject(run(flow, residual_limits, memory), core::FailureCode::InvalidInput);
        view.history.columns = history;
        collocation[0].tangents[1] = coordinates(0, -2, 0);
        reject(run(flow, residual_limits, memory), core::FailureCode::GeometryInvalid);
        collocation[0].tangents[1] = coordinates(0, 2, 0);
        collocation[0].velocity[2] = a::SingleWord{std::numeric_limits<double>::quiet_NaN()};
        reject(run(flow, residual_limits, memory), core::FailureCode::Nonfinite);
        collocation[0].velocity[2] = a::SingleWord{.5};
        auto exhausted = residual_limits;
        exhausted.steps = 0;
        reject(run(flow, exhausted, memory), core::FailureCode::ResourceExceeded);
        CHECK(flow.source().jet.active());

        const auto result = take(run(flow, residual_limits, memory));
        CHECK_FALSE(flow.source().jet.active());
        CHECK(result.size() == 1);
        CHECK(result.geometry_context() == geometry.context());
        CHECK(result.trial_stamp() == view.trial_stamp);
        const auto row = take(result.row(0));
        CHECK(row.face == fixture.faces[0].id);
        CHECK_FALSE(row.includes_wake_quadrature);
        CHECK(result.used_steps() <= residual_limits.steps);
        residual_expected(row.residual, residual);
        residual_expected(row.circulation_action, action);
        reject(result.row(1), core::FailureCode::InvalidInput);

        auto solve_limits = residual_limits;
        solve_limits.steps = UINT64_C(12000000000);
        const std::array unit_scale{1.};
        const std::array<ls::UnitDimensions, 1> gamma_units{{{2, 0, -1, 0, 0, 0, 0}}};
        const std::array<ls::UnitDimensions, 1> normal_units{{{1, 0, -1, 0, 0, 0, 0}}};
        const ls::Scales solve_scales{unit_scale, unit_scale, gamma_units, normal_units,
                                      circulation_context(geometry.context())};
        const auto request = circulation_request(solve_scales, 1, ls::Norm::Infinity);
        const auto solve = [&](cf::ProjectedCaseUniformFlow& projected,
                               const cf::DiscreteCirculationRequest& selected) {
            return cf::solve_discrete_circulation(view, std::move(projected), solve_limits,
                                                  selected, shared, *local, *kernel, memory);
        };
        auto solve_flow = projected_at(time);
        auto too_tight = request;
        too_tight.matrix_absolute = 0;
        reject(solve(solve_flow, too_tight), core::FailureCode::PrecisionUnresolved);
        CHECK(solve_flow.source().jet.active());
        too_tight = request;
        too_tight.original_absolute = 0;
        reject(solve(solve_flow, too_tight), core::FailureCode::LinearUnconverged);
        CHECK(solve_flow.source().jet.active());
        auto stale = request;
        stale.scales.context.trial = core::TrialId{99};
        reject(solve(solve_flow, stale), core::FailureCode::StaleGeneration);

        const auto solved = take(solve(solve_flow, request));
        CHECK_FALSE(solve_flow.source().jet.active());
        CHECK(solved.solution().values().precision() == core::Precision::TwoWord);
        const Rat expected = -Rat{Big{"10103660544137754684298736617717105223494942188"}} /
                             Rat{Big{"1000000000000000000000000000000000000000000000000"}};
        CHECK(abs(solution_word(solved.solution(), 0) - expected) <= rational_word(1e-15));
        const auto original = take(solved.row(0));
        CHECK(original.face == fixture.faces[0].id);
        residual_small(original.original_residual, rational_word(1e-13));
        reject(solved.row(1), core::FailureCode::InvalidInput);
        const std::array tangent_rhs{cf::CirculationLinearRhs{row.face, a::SingleWord{.125}}};
        const auto tangent = take(cf::solve_discrete_circulation_linearization(
            solved, cf::CirculationSense::Tangent, tangent_rhs, request, 100000, memory));
        const Rat expected_tangent{Big{"167540106538491744036640792288003354"},
                                   Big{"1000000000000000000000000000000000000"}};
        CHECK(abs(solution_word(tangent.solution(), 0) - expected_tangent) <= rational_word(1e-15));
        residual_small(take(tangent.row(0)).original_residual, rational_word(1e-13));

        const auto displacement = [&](std::size_t i) {
            const auto x = fixture.points[3 * i];
            const auto y = fixture.points[3 * i + 1];
            return coordinates(x / 8 + y / 16, -x / 32 + y / 8, x / 16 + y / 32);
        };
        const std::array displacements{displacement(0), displacement(1), displacement(2),
                                       displacement(3), displacement(4), displacement(5)};
        std::array row_directions{cf::CirculationCollocationDirection{
            row.face,
            coordinates(-.03125, .0625, .125),
            {coordinates(.0625, -.0625, .125), coordinates(-.125, .125, -.0625)},
            coordinates(0, .0625, .03125)}};
        const cf::DiscreteCirculationDirection direction{fixture.id, displacements, row_directions,
                                                         coordinates(.125, 0, .09375)};
        auto tangent_request = request;
        tangent_request.rhs_generation = core::RhsGeneration{78};
        tangent_request.retained_generation = core::NumericalGeneration{79};
        const auto run_direction = [&](const cf::DiscreteNonpenetrationInput& current_input,
                                       const cf::DiscreteCirculationDirection& selected_direction,
                                       cf::ProjectedCaseUniformFlow& projected,
                                       const cf::DiscreteCirculationRequest& selected_request) {
            return cf::solve_discrete_circulation_direction(
                current_input, selected_direction, std::move(projected), solve_limits, request,
                selected_request, shared, *local, *kernel, memory);
        };

        const auto check_collision = [&](const cf::DiscreteCirculationRequest& collision) {
            auto collision_flow = projected_at(time);
            const auto allocated = memory.allocated_bytes();
            const auto collision_result = run_direction(view, direction, collision_flow, collision);
            CHECK((!collision_result &&
                   collision_result.error().code == core::FailureCode::StaleGeneration));
            CHECK(collision_flow.source().jet.active());
            CHECK(memory.allocated_bytes() == allocated);
        };
        auto rhs_collision = tangent_request;
        rhs_collision.rhs_generation = request.rhs_generation;
        check_collision(rhs_collision);
        auto retained_collision = tangent_request;
        retained_collision.retained_generation = request.retained_generation;
        check_collision(retained_collision);

        auto direction_flow = projected_at(time);
        auto zero_projection = tangent_request;
        zero_projection.rhs_absolute = 0;
        reject(run_direction(view, direction, direction_flow, zero_projection),
               core::FailureCode::PrecisionUnresolved);
        CHECK(direction_flow.source().jet.active());
        auto stale_tangent = tangent_request;
        stale_tangent.scales.context.cache.matrix = core::NumericalGeneration{80};
        reject(run_direction(view, direction, direction_flow, stale_tangent),
               core::FailureCode::StaleGeneration);
        const std::array changed_scale{2 * tangent_request.scales.unknown.front()};
        stale_tangent = tangent_request;
        stale_tangent.scales.unknown = changed_scale;
        reject(run_direction(view, direction, direction_flow, stale_tangent),
               core::FailureCode::StaleGeneration);
        CHECK(direction_flow.source().jet.active());
        auto no_query_direction = direction;
        auto other_rows = row_directions;
        other_rows[0].free_query.reset();
        no_query_direction.rows = other_rows;
        reject(run_direction(view, no_query_direction, direction_flow, tangent_request),
               core::FailureCode::InvalidInput);
        if (!reversed) {
            auto bad_order = direction;
            auto ids = fixture.id;
            std::swap(ids[0], ids[1]);
            bad_order.points = ids;
            reject(run_direction(view, bad_order, direction_flow, tangent_request),
                   core::FailureCode::InvalidInput);
            auto bad_flow_direction = direction;
            bad_flow_direction.uniform_flow[2] =
                a::SingleWord{std::numeric_limits<double>::quiet_NaN()};
            reject(run_direction(view, bad_flow_direction, direction_flow, tangent_request),
                   core::FailureCode::Nonfinite);
        }
        CHECK(direction_flow.source().jet.active());

        const Rat reference_scale{Big{"1000000000000000000000000000000000000"}};
        const Rat expected_forcing =
            Rat{Big{"94584593959108917494617025902198312"}} / reference_scale;
        const Rat expected_direction =
            Rat{Big{"126773703590472725301051768620706643"}} / reference_scale;
        const auto checked = take(run_direction(view, direction, direction_flow, tangent_request));
        CHECK_FALSE(direction_flow.source().jet.active());
        CHECK(checked.used_steps() <= solve_limits.steps);
        CHECK(checked.tangent().used_steps() < checked.used_steps());
        CHECK(checked.primal().flow().source().jet.active());
        CHECK(checked.primal().solution().certificate().rhs_generation() == request.rhs_generation);
        CHECK(checked.tangent().solution().certificate().rhs_generation() ==
              tangent_request.rhs_generation);
        CHECK(checked.primal().solution().values().generation() == request.retained_generation);
        CHECK(checked.tangent().solution().values().generation() ==
              tangent_request.retained_generation);
        CHECK(abs(solution_word(checked.primal().solution(), 0) - expected) <=
              rational_word(1e-15));
        CHECK(abs(solution_word(checked.tangent().solution(), 0) - expected_direction) <=
              rational_word(1e-15));
        const auto forcing = take(checked.forcing(0));
        CHECK(forcing.face == row.face);
        residual_expected(forcing.residual_direction, expected_forcing, rational_word(1e-16));
        residual_small(take(checked.tangent().row(0)).original_residual, rational_word(1e-13));
        reject(checked.forcing(1), core::FailureCode::InvalidInput);

        const std::array<ls::UnitDimensions, 1> dual_speed{{{-1, 0, 1, 0, 0, 0, 0}}};
        const std::array<ls::UnitDimensions, 1> dual_gamma{{{-2, 0, 1, 0, 0, 0, 0}}};
        auto dual_context = solve_scales.context;
        dual_context.cache.matrix = core::NumericalGeneration{81};
        dual_context.cache.scales = core::ScaleGeneration{82};
        dual_context.unknown_layout = core::ScopeId{83};
        dual_context.residual_layout = core::ScopeId{84};
        const ls::Scales dual_scales{unit_scale, unit_scale, dual_speed, dual_gamma, dual_context};
        auto dual_request = circulation_request(dual_scales, 1, ls::Norm::Infinity);
        dual_request.rhs_generation = core::RhsGeneration{85};
        dual_request.retained_generation = core::NumericalGeneration{86};
        const auto objective = a::DoubleWord{.375, 0x1p-58};
        const std::array seed{cf::CirculationLinearRhs{row.face, objective}};
        const auto run_pullback = [&](const cf::DiscreteNonpenetrationInput& selected,
                                      std::span<const cf::CirculationLinearRhs> objective_seed,
                                      cf::ProjectedCaseUniformFlow& projected,
                                      const cf::DiscreteCirculationRequest& transpose_request) {
            return cf::solve_discrete_circulation_pullback(
                selected, objective_seed, std::move(projected), solve_limits, request,
                transpose_request, shared, *local, *kernel, memory);
        };
        auto pullback_flow = projected_at(time);
        auto exact_seed = dual_request;
        exact_seed.rhs_absolute = 0;
        reject(run_pullback(view, seed, pullback_flow, exact_seed),
               core::FailureCode::PrecisionUnresolved);
        CHECK(pullback_flow.source().jet.active());
        if (!reversed) {
            auto stale_dual = dual_request;
            stale_dual.scales.context.cache.matrix = solve_scales.context.cache.matrix;
            reject(run_pullback(view, seed, pullback_flow, stale_dual),
                   core::FailureCode::StaleGeneration);
            const std::array wrong_face{cf::CirculationLinearRhs{fixture.faces[1].id, objective}};
            reject(run_pullback(view, wrong_face, pullback_flow, dual_request),
                   core::FailureCode::InvalidInput);
            const std::array nonfinite{cf::CirculationLinearRhs{
                row.face, a::SingleWord{std::numeric_limits<double>::quiet_NaN()}}};
            reject(run_pullback(view, nonfinite, pullback_flow, dual_request),
                   core::FailureCode::Nonfinite);
            reject(run_pullback(view, {}, pullback_flow, dual_request),
                   core::FailureCode::InvalidInput);
            CHECK(pullback_flow.source().jet.active());
        }
        const auto pullback = take(run_pullback(view, seed, pullback_flow, dual_request));
        CHECK_FALSE(pullback_flow.source().jet.active());
        CHECK(pullback.used_steps() <= solve_limits.steps);
        CHECK(pullback.transpose().used_steps() < pullback.used_steps());
        CHECK(pullback.primal().flow().source().jet.active());
        CHECK(pullback.transpose().sense() == cf::CirculationSense::Transpose);
        CHECK(pullback.transpose().solution().values().precision() == core::Precision::TwoWord);
        const auto lambda = solution_word(pullback.transpose().solution(), 0);
        const auto pullback_row = take(pullback.row(0));
        const auto uniform_dual = pullback.uniform_flow();
        residual_expected(uniform_dual[2], -lambda);
        residual_expected(pullback_row.velocity[2], lambda);
        residual_expected(pullback_row.tangents[0][2], -lambda / 4);
        residual_expected(pullback_row.tangents[1][2], -lambda / 8);
        residual_small(take(pullback.transpose().row(0)).original_residual, rational_word(1e-13));
        const auto free_work = pullback_work(pullback, direction);
        const auto seed_value = sum_words(seed[0].value);
        CHECK(free_work.radius < rational_word(1e-24));
        CHECK(abs(free_work.value - seed_value * expected_direction) <=
              free_work.radius + rational_word(1e-16));
        CHECK(abs(free_work.value - seed_value * solution_word(checked.tangent().solution(), 0)) <=
              free_work.radius + rational_word(1e-16));
        reject(pullback.row(1), core::FailureCode::InvalidInput);
        reject(pullback.point(fixture.id.size()), core::FailureCode::InvalidInput);

        if (!reversed) {
            const auto material = test_geometry::take(a::resolve_query(
                geometry,
                {a::QueryId{30},
                 a::SourceFacePoint{a::FacetId{201}, a::SingleWord{.25}, a::SingleWord{.5}}},
                geometry.context(), test_geometry::limits(), *geometry_work, memory));
            const auto& original_row = collocation[0];
            const std::array material_rows{cf::DiscreteCollocation{
                original_row.facet, material, original_row.tangents, original_row.velocity,
                original_row.circulation, original_row.circulation_direction}};
            auto material_input = view;
            material_input.collocations = material_rows;
            auto material_flow = projected_at(time);
            reject(run_direction(material_input, direction, material_flow, tangent_request),
                   core::FailureCode::InvalidInput);
            const auto material_result = take(
                run_direction(material_input, no_query_direction, material_flow, tangent_request));
            const Rat material_forcing =
                Rat{Big{"86343522569361124446608159993825338"}} / reference_scale;
            const Rat material_tangent =
                Rat{Big{"115728023761435433573650003286112081"}} / reference_scale;
            residual_expected(take(material_result.forcing(0)).residual_direction, material_forcing,
                              rational_word(1e-16));
            CHECK(abs(solution_word(material_result.tangent().solution(), 0) - material_tangent) <=
                  rational_word(1e-15));

            auto material_pullback_flow = projected_at(time);
            const auto material_dual =
                take(run_pullback(material_input, seed, material_pullback_flow, dual_request));
            const auto material_work = pullback_work(material_dual, no_query_direction);
            CHECK(material_work.radius < rational_word(1e-24));
            CHECK(abs(material_work.value - seed_value * material_tangent) <=
                  material_work.radius + rational_word(1e-16));
            CHECK(abs(material_work.value -
                      seed_value * solution_word(material_result.tangent().solution(), 0)) <=
                  material_work.radius + rational_word(1e-16));
        }
    }
}

TEST_CASE("Discrete circulation loads keep midpoint and triangle force and moment",
          "[UW-D01][UW-D03][UW-D05][UW-D06][N4-06][verification][negative]") {
    const auto source_case = source(R"({"kind":"constant","value":[0,0,100]})", one_step);
    a::UniformFlowWorkspace flow_workspace;
    const auto provider = import(source_case, io::SourceId{120}, flow_workspace);
    const auto current = import(source_case, io::SourceId{121}, flow_workspace);
    const Number one{1}, zero{0}, absolute{1, 1000000000}, time{5, 4};
    const cf::CaseFlowComponentBudget component{one.view(), absolute.view(), zero.view()};
    const std::array budgets{component, component, component};
    const cf::CaseFlowProjectionRequest projection{
        core::Precision::TwoWord, budgets, budgets, {4096, budget}};
    const auto flow_at = [&] {
        auto exact = evaluate(provider, current, flow_workspace, time);
        return take(cf::project_case_uniform_flow(std::move(exact), projection));
    };
    const auto expected = [](const cf::DiscreteLoadVector& value, std::array<double, 3> reference) {
        for (std::size_t axis = 0; axis < value.size(); ++axis)
            residual_expected(value[axis], rational_word(reference[axis]), rational_word(4e-15));
    };
    const a::GlobalLimits cap{std::numeric_limits<std::uint64_t>::max() / 4,
                              numerics::linear::EstimateAdmission{4096}};
    const cf::DiscreteNonpenetrationLimits load_limits{
        {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10},
        {test_geometry::limits(),
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(24000000000)};
    const std::array<ls::UnitDimensions, 2> gamma_units{
        {{2, 0, -1, 0, 0, 0, 0}, {2, 0, -1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 2> speed_units{
        {{1, 0, -1, 0, 0, 0, 0}, {1, 0, -1, 0, 0, 0, 0}}};
    const std::array scales{1., 1.};
    const auto moment_origin = coordinates(.25, -.125, .5);

    enum class LoadCase : std::uint8_t { Wake, ReversedWake, Shared, AliasedShared, ZeroRadius };
    for (const auto mode : {LoadCase::Wake, LoadCase::ReversedWake, LoadCase::Shared,
                            LoadCase::AliasedShared, LoadCase::ZeroRadius}) {
        CAPTURE(mode);
        const bool has_wake = mode != LoadCase::Shared && mode != LoadCase::AliasedShared;
        const bool reversed = mode == LoadCase::ReversedWake;
        ResidualFixture fixture;
        if (reversed)
            fixture.reverse_wake();
        if (!has_wake)
            fixture.roles[1].role = a::AeroFaceRole::Bound;
        auto topology_input = fixture.topology();
        auto geometry_input = fixture.input();
        if (!has_wake)
            topology_input.trailing_edges = {};
        if (mode == LoadCase::ZeroRadius)
            geometry_input.core_radius = a::SingleWord{0};
        const a::EdgeId alias{90};
        const std::array with_alias{
            fixture.edges[0], fixture.edges[1],
            fixture.edges[2], fixture.edges[3],
            fixture.edges[4], fixture.edges[5],
            fixture.edges[6], fixture.edges[7],
            fixture.edges[8], a::EdgeInput{alias, {fixture.id[3], fixture.id[2]}}};
        const std::array owners_with_alias{
            fixture.edge_owners[0], fixture.edge_owners[1],
            fixture.edge_owners[2], fixture.edge_owners[3],
            fixture.edge_owners[4], fixture.edge_owners[5],
            fixture.edge_owners[6], fixture.edge_owners[7],
            fixture.edge_owners[8], a::EdgeBindingInput{alias, a::EdgeId{12}}};
        if (mode == LoadCase::AliasedShared) {
            fixture.quads[1].edges[0] = alias;
            geometry_input.edges = with_alias;
            topology_input.edges = owners_with_alias;
        }
        core::ResourcePlan memory{128ULL * 1024 * 1024};
        auto geometry_work = test_geometry::workspace(memory);
        const auto geometry = test_geometry::take(a::prepare_owned_geometry(
            geometry_input, test_geometry::limits(), *geometry_work, memory));
        auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
        const auto topology =
            take_global(a::build_global_topology(geometry, topology_input, cap, *global, memory));
        const auto contact =
            take_global(a::check_global_geometry(topology, geometry, cap, *global, memory));
        const auto collocate = [&](a::QueryId id, std::array<double, 3> point) {
            return test_geometry::take(a::resolve_query(
                geometry, {id, a::FreePoint{a::FreeSingle{point}}}, geometry.context(),
                test_geometry::limits(), *geometry_work, memory));
        };
        const auto first = collocate(a::QueryId{90}, {.75, .5, 0});
        const auto second = collocate(a::QueryId{91}, {1.625, .375, 0});
        const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
        const std::array rows{
            cf::DiscreteCollocation{a::FacetId{200}, first, frame, coordinates(0, 0, .5),
                                    a::SingleWord{3}, a::SingleWord{0}},
            cf::DiscreteCollocation{a::FacetId{202}, second, frame, coordinates(0, 0, .75),
                                    a::SingleWord{-4}, a::SingleWord{0}}};
        const std::array wake{
            cf::DiscreteWakePanel{a::FacetId{202}, cf::AttachedWakeFront{}, a::P1StrengthId{0}}};
        const std::array old{
            cf::FrozenWakeTrace{a::P1StrengthId{0}, a::SingleWord{reversed ? .25 : -.25}}};
        const std::size_t count = has_wake ? 1U : 2U;
        const cf::DiscreteNonpenetrationInput input{
            contact,
            std::span{rows}.first(count),
            has_wake ? std::span<const cf::DiscreteWakePanel>{wake}
                     : std::span<const cf::DiscreteWakePanel>{},
            {core::AcceptedStateId{2}, has_wake ? std::span<const cf::FrozenWakeTrace>{old}
                                                : std::span<const cf::FrozenWakeTrace>{}},
            {70, 71},
            stamp,
            p::RequestedSide::Automatic};
        const ls::Scales model{std::span{scales}.first(count), std::span{scales}.first(count),
                               std::span{gamma_units}.first(count),
                               std::span{speed_units}.first(count),
                               circulation_context(geometry.context())};
        const auto request =
            circulation_request(model, static_cast<core::Index>(count), ls::Norm::Infinity);
        const std::array rates{
            cf::BoundCirculationRate{fixture.faces[0].id, a::DoubleWord{.375, 0x1p-57}},
            cf::BoundCirculationRate{fixture.faces[1].id, a::SingleWord{-.125}}};
        const std::array midpoint_velocities{
            cf::EdgeMidpointVelocity{a::EdgeId{10}, coordinates(.125, -.25, .375)},
            cf::EdgeMidpointVelocity{a::EdgeId{11}, coordinates(.25, .125, .125)},
            cf::EdgeMidpointVelocity{a::EdgeId{12}, coordinates(.25, 0, .125)},
            cf::EdgeMidpointVelocity{a::EdgeId{13}, coordinates(-.125, .25, -.25)},
            cf::EdgeMidpointVelocity{a::EdgeId{15}, coordinates(.125, -.25, .125)},
            cf::EdgeMidpointVelocity{a::EdgeId{16}, coordinates(0, .5, .25)},
            cf::EdgeMidpointVelocity{a::EdgeId{17}, coordinates(.125, .25, -.125)}};
        const std::array wake_force_edges{midpoint_velocities[0], midpoint_velocities[1],
                                          midpoint_velocities[3]};
        const cf::DiscreteVectorLoadInput load_input{
            input.trial_stamp,
            request.retained_generation,
            has_wake ? std::span<const cf::EdgeMidpointVelocity>{wake_force_edges}
                     : std::span<const cf::EdgeMidpointVelocity>{midpoint_velocities},
            std::span{rates}.first(count),
            a::SingleWord{1.25},
            moment_origin};
        a::SharedP1Workspace shared;
        auto local = take(a::LocalPacketWorkspace::create(memory));
        auto kernel = take(a::P1Workspace::create(memory));
        auto flow = flow_at();
        const auto run = [&](const cf::DiscreteVectorLoadInput& load,
                             const std::stop_token& stop = {}) {
            return cf::solve_discrete_circulation_loads(input, load, std::move(flow), load_limits,
                                                        request, *geometry_work, shared, *local,
                                                        *kernel, memory, stop);
        };
        if (mode == LoadCase::ZeroRadius) {
            reject(run(load_input), core::FailureCode::GeometryInvalid);
            CHECK(flow.source().jet.active());
            continue;
        }
        if (mode == LoadCase::Wake) {
            // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
            std::stop_source cancellation;
            REQUIRE(cancellation.request_stop());
            const auto before = memory.allocated_bytes();
            reject(run(load_input, cancellation.get_token()), core::FailureCode::Cancelled);
            CHECK(memory.allocated_bytes() == before);
            auto invalid = load_input;
            invalid.circulation = core::NumericalGeneration{901};
            reject(run(invalid), core::FailureCode::StaleGeneration);
            invalid = load_input;
            invalid.rates = {};
            reject(run(invalid), core::FailureCode::InvalidInput);
            invalid = load_input;
            invalid.density = a::SingleWord{0};
            reject(run(invalid), core::FailureCode::InvalidInput);
            invalid = load_input;
            invalid.edges = std::span{midpoint_velocities}.first(3);
            reject(run(invalid), core::FailureCode::InvalidInput); // TE is not a force edge.
            invalid = load_input;
            auto bad_rate = rates;
            bad_rate[0].rate = a::SingleWord{std::numeric_limits<double>::quiet_NaN()};
            invalid.rates = std::span{bad_rate}.first(1);
            reject(run(invalid), core::FailureCode::Nonfinite);
            CHECK(flow.source().jet.active());
        }
        const auto result = take(run(load_input));
        CHECK_FALSE(flow.source().jet.active());
        CHECK(result.primal().flow().source().jet.active());
        CHECK(result.used_steps() <= load_limits.steps);
        CHECK(result.edge_count() == load_input.edges.size());
        CHECK(result.facet_count() == 2 * count);
        CHECK(result.point_count() == (has_wake ? 4U : 6U));
        expected(result.moment_origin(), {.25, -.125, .5});
        for (std::size_t i = 0; i < result.edge_count(); ++i) {
            const auto edge = take(result.edge(i));
            CHECK(edge.edge == load_input.edges[i].edge);
            CHECK(edge.excluded_bound_occurrences == (edge.edge == a::EdgeId{12} ? 2U : 1U));
        }
        const auto first_edge = take(result.edge(0));
        const auto first_facet = take(result.facet(0));
        CHECK(first_facet.face == fixture.faces[0].id);
        CHECK(first_facet.facet == a::FacetId{200});
        expected(first_edge.position, {0, .5, 0});
        expected(first_facet.position, {1. / 3, 2. / 3, 0});
        const auto exact_force = (Rat{3, 8} + two_power(-57)) * Rat{5, 8};
        residual_expected(first_facet.force[2], exact_force);
        for (const auto axis : {0U, 1U})
            residual_expected(first_facet.force[axis], Rat{0});

        if (has_wake) {
            expected(first_edge.force, {.00169150653666965886, 0, .00157869696002152417});
            expected(take(result.edge(1)).force,
                     {0, -.00480275125434878304, -.00157869696002152417});
            expected(take(result.edge(2)).force, {0, .00953884213441335554, .00315739392004304834});
            expected(result.total_force(),
                     {.00169150653666965886, .00473609088006457251, .471907393920043057});
            expected(result.total_moment(),
                     {.294942121200026911, -.118033253268334832, .000126831134597606338});
            expected(take(result.point(0)).force,
                     {.000845753268334829431, .00476942106720667777, .158618045440032289});
        } else {
            expected(take(result.edge(2)).force, {.132938445560404041, 0, .0665551441246223153});
            expected(result.total_force(),
                     {.0228578348623865182, .118653830124245513, .574683855573414792});
            expected(result.total_moment(),
                     {.402874718995620042, -.201706821835838983, .102771929266541399});
        }

        std::array<Rat, 3> nodal_force{}, nodal_moment{};
        const std::array origin{Rat{1, 4}, -Rat{1, 8}, Rat{1, 2}};
        const auto unpack = [](const a::ScaledEnclosure& enclosure) {
            return sum_words(enclosure.center()) * two_power(enclosure.exponent());
        };
        for (std::size_t i = 0; i < result.point_count(); ++i) {
            const auto point = take(result.point(i));
            CHECK(point.point == fixture.id[i]);
            std::array<Rat, 3> local_force{}, lever{};
            for (std::size_t axis = 0; axis < 3; ++axis) {
                local_force[axis] = unpack(point.force[axis]);
                lever[axis] = rational_word(fixture.points[3 * i + axis]) - origin[axis];
                nodal_force[axis] += local_force[axis];
            }
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const auto j = (axis + 1) % 3;
                const auto k = (axis + 2) % 3;
                nodal_moment[axis] += lever[j] * local_force[k] - lever[k] * local_force[j];
            }
        }
        const auto force = result.total_force();
        const auto moment = result.total_moment();
        for (std::size_t axis = 0; axis < 3; ++axis) {
            CHECK(abs(nodal_force[axis] - unpack(force[axis])) <= rational_word(1e-24));
            CHECK(abs(nodal_moment[axis] - unpack(moment[axis])) <= rational_word(1e-24));
        }
        reject(result.edge(result.edge_count()), core::FailureCode::InvalidInput);
        reject(result.point(result.point_count()), core::FailureCode::InvalidInput);
        reject(result.facet(result.facet_count()), core::FailureCode::InvalidInput);

        const auto displacement = [&](std::size_t i) {
            const auto x = fixture.points[3 * i];
            const auto y = fixture.points[3 * i + 1];
            return coordinates(x / 8 + y / 16, -x / 32 + y / 8, x / 16 + y / 32);
        };
        const std::array motion{displacement(0), displacement(1), displacement(2),
                                displacement(3), displacement(4), displacement(5)};
        const auto velocity_direction = [](const cf::EdgeMidpointVelocity& entry) {
            return cf::EdgeMidpointVelocity{
                entry.edge,
                coordinates(.03125, -.0625, static_cast<double>(entry.edge.value()) / 256)};
        };
        const std::array all_midpoint_variations{
            velocity_direction(midpoint_velocities[0]), velocity_direction(midpoint_velocities[1]),
            velocity_direction(midpoint_velocities[2]), velocity_direction(midpoint_velocities[3]),
            velocity_direction(midpoint_velocities[4]), velocity_direction(midpoint_velocities[5]),
            velocity_direction(midpoint_velocities[6])};
        const std::array wake_midpoint_variations{
            all_midpoint_variations[0], all_midpoint_variations[1], all_midpoint_variations[3]};
        const std::array gamma_directions{
            cf::BoundLoadDirection{fixture.faces[0].id, a::DoubleWord{.1875, 0x1p-56},
                                   a::DoubleWord{-.25, 0x1p-58}},
            cf::BoundLoadDirection{fixture.faces[1].id, a::SingleWord{-.15625},
                                   a::SingleWord{.0625}}};
        const cf::DiscreteVectorLoadDirection variation{
            input.trial_stamp,
            request.retained_generation,
            fixture.id,
            motion,
            has_wake ? std::span<const cf::EdgeMidpointVelocity>{wake_midpoint_variations}
                     : std::span<const cf::EdgeMidpointVelocity>{all_midpoint_variations},
            std::span{gamma_directions}.first(count),
            coordinates(.125, -.0625, .09375)};
        auto jvp_flow = flow_at();
        const auto run_direction = [&](const cf::DiscreteVectorLoadDirection& delta) {
            return cf::solve_discrete_circulation_load_direction(
                input, load_input, delta, std::move(jvp_flow), load_limits, request, *geometry_work,
                shared, *local, *kernel, memory);
        };
        if (mode == LoadCase::Wake) {
            auto invalid = variation;
            invalid.circulation = core::NumericalGeneration{902};
            const auto before = memory.allocated_bytes();
            reject(run_direction(invalid), core::FailureCode::StaleGeneration);
            CHECK(memory.allocated_bytes() == before);
            invalid = variation;
            invalid.displacements = std::span{motion}.first(5);
            reject(run_direction(invalid), core::FailureCode::InvalidInput);
            invalid = variation;
            auto bad_rows = gamma_directions;
            bad_rows[0].circulation = a::SingleWord{std::numeric_limits<double>::quiet_NaN()};
            invalid.bound = std::span{bad_rows}.first(1);
            reject(run_direction(invalid), core::FailureCode::Nonfinite);
            invalid = variation;
            auto wrong_edge_directions = wake_midpoint_variations;
            std::swap(wrong_edge_directions[0], wrong_edge_directions[1]);
            invalid.edges = wrong_edge_directions;
            reject(run_direction(invalid), core::FailureCode::InvalidInput);
            CHECK(jvp_flow.source().jet.active());
        }
        const auto jvp = take(run_direction(variation));
        CHECK_FALSE(jvp_flow.source().jet.active());
        CHECK(jvp.loads().primal().flow().source().jet.active());
        CHECK(jvp.loads().used_steps() <= load_limits.steps);
        const auto edge_direction = take(jvp.edge(0));
        const auto facet_direction = take(jvp.facet(0));
        expected(edge_direction.position, {.03125, .0625, .015625});
        expected(facet_direction.position, {1. / 12, 7. / 96, 1. / 24});
        expected(facet_direction.force, {-.0146484375, -.00732421875, -.09765625});
        const auto exact_direction = Rat{5, 4} * ((-Rat{1, 4} + two_power(-58)) * Rat{1, 2} +
                                                  (Rat{3, 8} + two_power(-57)) * Rat{1, 8});
        residual_expected(facet_direction.force[2], exact_direction);
        if (has_wake) {
            expected(edge_direction.force,
                     {-.0313421017332066007, -.000155053438542526309, -.0300854876693103419});
            expected(take(jvp.edge(1)).force,
                     {-.0000514174166970542093, .0883611191373407818, .0290019435772196279});
            expected(take(jvp.edge(2)).force,
                     {.00010075169669772684, -.176050182858272690, -.0582494842972189574});
            expected(jvp.total_force_direction(),
                     {-.0605896424532059286, -.102492554659474434, -.254645528389309665});
            expected(jvp.total_moment_direction(),
                     {-.144922614949141765, .0349452424989099995, .0122891137124191336});
            expected(take(jvp.point(0)).force,
                     {-.0253863000182544371, -.0929854306484076081, -.109271652649931314});
        } else {
            expected(take(jvp.edge(2)).force,
                     {.237429517694973129, -.0103885011014196999, .0927190867664370914});
            expected(jvp.total_force_direction(),
                     {.0584661631586267516, -.0431664920836314646, .0342573010497334486});
            expected(jvp.total_moment_direction(),
                     {.0598305790447179028, -.387825417382047371, .0294373849541907755});
        }
        std::array<Rat, 3> sum_df{}, sum_dm{};
        for (std::size_t i = 0; i < jvp.loads().point_count(); ++i) {
            const auto value = take(jvp.loads().point(i));
            const auto delta = take(jvp.point(i));
            CHECK(value.point == delta.point);
            std::array<Rat, 3> r{}, dx{}, f{}, df{};
            for (std::size_t axis = 0; axis < r.size(); ++axis) {
                r[axis] = rational_word(fixture.points[3 * i + axis]) - origin[axis];
                dx[axis] = sum_words(motion[i][axis]);
                f[axis] = unpack(value.force[axis]);
                df[axis] = unpack(delta.force[axis]);
                sum_df[axis] += df[axis];
            }
            for (std::size_t axis = 0; axis < r.size(); ++axis) {
                const auto j = (axis + 1) % 3;
                const auto k = (axis + 2) % 3;
                sum_dm[axis] += dx[j] * f[k] - dx[k] * f[j] + r[j] * df[k] - r[k] * df[j];
            }
        }
        for (std::size_t axis = 0; axis < sum_df.size(); ++axis) {
            CHECK(abs(sum_df[axis] - unpack(jvp.total_force_direction()[axis])) <=
                  rational_word(1e-24));
            CHECK(abs(sum_dm[axis] - unpack(jvp.total_moment_direction()[axis])) <=
                  rational_word(1e-24));
        }
        reject(jvp.edge(result.edge_count()), core::FailureCode::InvalidInput);
        reject(jvp.point(result.point_count()), core::FailureCode::InvalidInput);
        reject(jvp.facet(result.facet_count()), core::FailureCode::InvalidInput);

        const auto zeros = coordinates(0, 0, 0);
        const cf::DiscreteLoadSampleWeights no_weight{zeros, zeros};
        const std::array edge_weights{cf::DiscreteLoadSampleWeights{coordinates(.125, .0625, -.125),
                                                                    coordinates(.25, -.125, .5)},
                                      no_weight,
                                      no_weight,
                                      no_weight,
                                      no_weight,
                                      no_weight,
                                      no_weight};
        const std::array facet_weights{
            cf::DiscreteLoadSampleWeights{coordinates(-.25, .125, .0625),
                                          coordinates(.0625, .25, -.125)},
            no_weight, no_weight, no_weight};
        const std::array node_weights{
            coordinates(-.375, .125, .25), zeros, zeros, zeros, zeros, zeros};
        const cf::NonpenetrationVector total_weight{a::SingleWord{.125}, a::SingleWord{-.25},
                                                    a::DoubleWord{.375, 0x1p-57}};
        const cf::DiscreteVectorLoadWeights weights{
            input.trial_stamp,
            request.retained_generation,
            std::span{edge_weights}.first(result.edge_count()),
            std::span{facet_weights}.first(result.facet_count()),
            std::span{node_weights}.first(result.point_count()),
            total_weight,
            coordinates(-.125, .25, .0625)};
        auto vjp_flow = flow_at();
        const auto run_pullback = [&](const cf::DiscreteVectorLoadWeights& dual) {
            return cf::solve_discrete_circulation_load_pullback(
                input, load_input, dual, std::move(vjp_flow), load_limits, request, *geometry_work,
                shared, *local, *kernel, memory);
        };
        if (mode == LoadCase::Wake) {
            auto invalid = weights;
            invalid.circulation = core::NumericalGeneration{903};
            const auto before = memory.allocated_bytes();
            reject(run_pullback(invalid), core::FailureCode::StaleGeneration);
            CHECK(memory.allocated_bytes() == before);
            invalid = weights;
            invalid.facets = {};
            reject(run_pullback(invalid), core::FailureCode::InvalidInput);
            invalid = weights;
            invalid.points = std::span{node_weights}.first(3);
            reject(run_pullback(invalid), core::FailureCode::InvalidInput);
            invalid = weights;
            invalid.moment[0] = a::SingleWord{std::numeric_limits<double>::quiet_NaN()};
            reject(run_pullback(invalid), core::FailureCode::Nonfinite);
            CHECK(vjp_flow.source().jet.active());
        }
        const auto vjp = take(run_pullback(weights));
        CHECK_FALSE(vjp_flow.source().jet.active());
        CHECK(vjp.loads().primal().flow().source().jet.active());
        CHECK(vjp.loads().used_steps() <= load_limits.steps);
        const auto load_work = [&](const cf::DiscreteCirculationLoadJvp& direction) {
            WorkEnclosure work;
            work.add(direction.total_force_direction(), weights.force);
            work.add(direction.total_moment_direction(), weights.moment);
            for (std::size_t i = 0; i < weights.edges.size(); ++i) {
                const auto delta = take(direction.edge(i));
                work.add(delta.position, weights.edges[i].position);
                work.add(delta.force, weights.edges[i].force);
            }
            for (std::size_t i = 0; i < weights.facets.size(); ++i) {
                const auto delta = take(direction.facet(i));
                work.add(delta.position, weights.facets[i].position);
                work.add(delta.force, weights.facets[i].force);
            }
            for (std::size_t i = 0; i < weights.points.size(); ++i)
                work.add(take(direction.point(i)).force, weights.points[i]);
            return work;
        };
        const auto forward = load_work(jvp);

        WorkEnclosure reverse;
        reverse.add(vjp.uniform_flow(), variation.uniform_flow);
        REQUIRE(vjp.point_count() == variation.points.size());
        for (std::size_t i = 0; i < vjp.point_count(); ++i) {
            const auto dual = take(vjp.point(i));
            CHECK(dual.point == variation.points[i]);
            reverse.add(dual.position, variation.displacements[i]);
        }
        for (std::size_t i = 0; i < variation.edges.size(); ++i) {
            const auto dual = take(vjp.midpoint_velocity(i));
            CHECK(dual.edge == variation.edges[i].edge);
            reverse.add(dual.velocity, variation.edges[i].velocity);
        }
        for (std::size_t i = 0; i < variation.bound.size(); ++i) {
            const auto dual = take(vjp.bound(i));
            CHECK(dual.face == variation.bound[i].face);
            reverse.add(dual.circulation, variation.bound[i].circulation);
            reverse.add(dual.circulation_rate, variation.bound[i].circulation_rate);
        }
        CHECK(abs(forward.value - reverse.value) <= forward.radius + reverse.radius);
        CHECK(forward.radius < rational_word(1e-24));
        CHECK(reverse.radius < rational_word(1e-24));
        const auto independent =
            rational_word(has_wake ? -.0958978088767875470302 : -.120307262620136262701);
        CHECK(abs(reverse.value - independent) <= reverse.radius + rational_word(1e-16));
        reject(vjp.point(vjp.point_count()), core::FailureCode::InvalidInput);
        reject(vjp.midpoint_velocity(weights.edges.size()), core::FailureCode::InvalidInput);
        reject(vjp.bound(count), core::FailureCode::InvalidInput);

        if (mode != LoadCase::Wake && mode != LoadCase::AliasedShared)
            continue;
        const std::array<ls::UnitDimensions, 2> dual_speed{
            {{-1, 0, 1, 0, 0, 0, 0}, {-1, 0, 1, 0, 0, 0, 0}}};
        const std::array<ls::UnitDimensions, 2> dual_gamma{
            {{-2, 0, 1, 0, 0, 0, 0}, {-2, 0, 1, 0, 0, 0}}};
        auto dual_context = model.context;
        dual_context.cache.matrix = core::NumericalGeneration{1000};
        dual_context.cache.scales = core::ScaleGeneration{1001};
        dual_context.unknown_layout = core::ScopeId{1002};
        dual_context.residual_layout = core::ScopeId{1003};
        const ls::Scales dual_model{model.residual, model.unknown,
                                    std::span{dual_speed}.first(count),
                                    std::span{dual_gamma}.first(count), dual_context};
        auto transpose_request =
            circulation_request(dual_model, static_cast<core::Index>(count), ls::Norm::Infinity);
        transpose_request.rhs_generation = core::RhsGeneration{1004};
        transpose_request.retained_generation = core::NumericalGeneration{1005};
        auto implicit_flow = flow_at();
        const auto run_implicit = [&](const cf::DiscreteCirculationRequest& transpose,
                                      const std::stop_token& stop = {}) {
            return cf::solve_discrete_circulation_implicit_load_pullback(
                input, load_input, weights, std::move(implicit_flow), load_limits, request,
                transpose, *geometry_work, shared, *local, *kernel, memory, stop);
        };
        if (has_wake) {
            const auto before = memory.allocated_bytes();
            auto stale = transpose_request;
            stale.scales.context.cache.matrix = model.context.cache.matrix;
            reject(run_implicit(stale), core::FailureCode::StaleGeneration);
            stale = transpose_request;
            stale.scales.context.cache.scales = model.context.cache.scales;
            reject(run_implicit(stale), core::FailureCode::StaleGeneration);
            stale = transpose_request;
            stale.scales.context.unknown_layout = model.context.residual_layout;
            reject(run_implicit(stale), core::FailureCode::StaleGeneration);
            stale = transpose_request;
            stale.rhs_generation = request.rhs_generation;
            reject(run_implicit(stale), core::FailureCode::StaleGeneration);
            stale = transpose_request;
            stale.retained_generation = request.retained_generation;
            reject(run_implicit(stale), core::FailureCode::StaleGeneration);
            // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
            std::stop_source cancel;
            REQUIRE(cancel.request_stop());
            reject(run_implicit(transpose_request, cancel.get_token()),
                   core::FailureCode::Cancelled);
            CHECK(memory.allocated_bytes() == before);
            const auto original_gamma = solution_word(result.primal().solution(), 0);
            auto tight = transpose_request;
            tight.rhs_absolute = 0;
            reject(run_implicit(tight), core::FailureCode::PrecisionUnresolved);
            tight = transpose_request;
            tight.original_absolute = 0;
            reject(run_implicit(tight), core::FailureCode::LinearUnconverged);
            CHECK(memory.allocated_bytes() >
                  before); // The plan counts successful allocations cumulatively.
            CHECK(solution_word(result.primal().solution(), 0) == original_gamma);
            CHECK(contact.active());
            CHECK(implicit_flow.source().jet.active());
        }
        const auto combined = take(run_implicit(transpose_request));
        CHECK_FALSE(implicit_flow.source().jet.active());
        CHECK(combined.loads().primal().flow().source().jet.active());
        CHECK(combined.loads().used_steps() <= load_limits.steps);
        CHECK(combined.transpose().sense() == cf::CirculationSense::Transpose);
        const auto& certificate = combined.transpose().solution().certificate();
        CHECK(certificate.context() == dual_context);
        CHECK(certificate.rhs_generation() == transpose_request.rhs_generation);
        CHECK(combined.transpose().solution().values().generation() ==
              transpose_request.retained_generation);

        const std::array collocation_directions{
            cf::CirculationCollocationDirection{
                fixture.faces[0].id,
                coordinates(.03125, -.046875, .0234375),
                {coordinates(.03125, .0625, .125), coordinates(-.09375, .015625, -.0625)},
                coordinates(.0625, -.03125, .125)},
            cf::CirculationCollocationDirection{
                fixture.faces[1].id,
                coordinates(-.046875, .0625, -.0390625),
                {coordinates(-.0625, .03125, -.125), coordinates(.125, -.03125, .0625)},
                coordinates(-.03125, .0625, -.09375)}};
        const cf::DiscreteCirculationDirection circulation_direction{
            variation.points, variation.displacements,
            std::span{collocation_directions}.first(count), variation.uniform_flow};
        auto tangent_request = request;
        tangent_request.rhs_generation = core::RhsGeneration{1100};
        tangent_request.retained_generation = core::NumericalGeneration{1101};
        const std::array rate_directions{
            cf::BoundCirculationRate{gamma_directions[0].face,
                                     gamma_directions[0].circulation_rate},
            cf::BoundCirculationRate{gamma_directions[1].face,
                                     gamma_directions[1].circulation_rate}};
        const cf::DiscreteImplicitVectorLoadDirection implicit_variation{
            variation.stamp, variation.circulation, circulation_direction, variation.edges,
            std::span{rate_directions}.first(count)};
        auto tangent_flow = flow_at();
        const auto run_tangent = [&](const cf::DiscreteCirculationRequest& settings,
                                     const cf::DiscreteImplicitVectorLoadDirection& direction,
                                     const std::stop_token& stop = {}) {
            return cf::solve_discrete_circulation_implicit_load_direction(
                input, load_input, direction, std::move(tangent_flow), load_limits, request,
                settings, *geometry_work, shared, *local, *kernel, memory, stop);
        };
        if (has_wake) {
            const auto before = memory.allocated_bytes();
            auto stale = tangent_request;
            stale.rhs_generation = request.rhs_generation;
            reject(run_tangent(stale, implicit_variation), core::FailureCode::StaleGeneration);
            stale = tangent_request;
            stale.retained_generation = request.retained_generation;
            reject(run_tangent(stale, implicit_variation), core::FailureCode::StaleGeneration);
            stale = tangent_request;
            stale.scales.context.cache.matrix = transpose_request.scales.context.cache.matrix;
            reject(run_tangent(stale, implicit_variation), core::FailureCode::StaleGeneration);
            auto missing = implicit_variation;
            missing.rates = {};
            reject(run_tangent(tangent_request, missing), core::FailureCode::InvalidInput);
            // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
            std::stop_source cancel;
            REQUIRE(cancel.request_stop());
            reject(run_tangent(tangent_request, implicit_variation, cancel.get_token()),
                   core::FailureCode::Cancelled);
            CHECK(memory.allocated_bytes() == before);

            const std::array bad_scale{2.};
            stale = tangent_request;
            stale.scales.unknown = bad_scale;
            reject(run_tangent(stale, implicit_variation), core::FailureCode::StaleGeneration);
            auto strict = tangent_request;
            strict.rhs_absolute = 0;
            reject(run_tangent(strict, implicit_variation), core::FailureCode::PrecisionUnresolved);
            strict = tangent_request;
            strict.original_absolute = 0;
            reject(run_tangent(strict, implicit_variation), core::FailureCode::LinearUnconverged);
            CHECK(tangent_flow.source().jet.active());
            CHECK(result.primal().flow().source().jet.active());
        }
        const auto tangent = take(run_tangent(tangent_request, implicit_variation));
        CHECK_FALSE(tangent_flow.source().jet.active());
        CHECK(tangent.direction().loads().primal().flow().source().jet.active());
        CHECK(tangent.direction().loads().used_steps() <= load_limits.steps);
        CHECK(tangent.tangent().sense() == cf::CirculationSense::Tangent);
        CHECK(tangent.tangent().solution().certificate().context() == model.context);
        CHECK(tangent.tangent().solution().certificate().rhs_generation() ==
              tangent_request.rhs_generation);
        CHECK(tangent.tangent().solution().values().generation() ==
              tangent_request.retained_generation);
        const auto implicit_forward = load_work(tangent.direction());
        if (has_wake) {
            const auto& storage = tangent.tangent().solution().values().storage();
            const auto* words = std::get_if<numerics::TwoWord<double>>(&storage);
            if (words == nullptr)
                FAIL("Expected a retained two-word tangent");
            CHECK(words->lo.values()[0] != 0.);
            auto exact_directions = gamma_directions;
            exact_directions[0].circulation =
                a::DoubleWord{words->hi.values()[0], words->lo.values()[0]};
            auto explicit_variation = variation;
            explicit_variation.bound = std::span{exact_directions}.first(count);
            auto explicit_flow = flow_at();
            const auto explicit_path = take(cf::solve_discrete_circulation_load_direction(
                input, load_input, explicit_variation, std::move(explicit_flow), load_limits,
                request, *geometry_work, shared, *local, *kernel, memory));
            const auto explicit_work = load_work(explicit_path);
            CHECK(abs(explicit_work.value - implicit_forward.value) <=
                  explicit_work.radius + implicit_forward.radius);
            CHECK(implicit_forward.radius < rational_word(1e-24));
        }
        WorkEnclosure implicit_reverse;
        implicit_reverse.add(combined.uniform_flow(), variation.uniform_flow);
        REQUIRE(combined.point_count() == variation.points.size());
        for (std::size_t i = 0; i < combined.point_count(); ++i) {
            const auto dual = take(combined.point(i));
            CHECK(dual.point == variation.points[i]);
            implicit_reverse.add(dual.position, variation.displacements[i]);
        }
        for (std::size_t i = 0; i < variation.edges.size(); ++i) {
            const auto dual = take(combined.midpoint_velocity(i));
            CHECK(dual.edge == variation.edges[i].edge);
            implicit_reverse.add(dual.velocity, variation.edges[i].velocity);
        }
        for (std::size_t i = 0; i < count; ++i) {
            const auto tangent_row = take(tangent.tangent().row(i));
            CHECK(tangent_row.face == variation.bound[i].face);
            residual_expected(tangent_row.original_residual, Rat{0}, rational_word(5e-15));
            const auto check = take(combined.transpose().row(i));
            CHECK(check.face == variation.bound[i].face);
            residual_expected(check.original_residual, Rat{0}, rational_word(5e-15));
            const auto rate = take(combined.circulation_rate(i));
            CHECK(rate.face == variation.bound[i].face);
            implicit_reverse.add(rate.rate, variation.bound[i].circulation_rate);
            const auto dual = take(combined.collocation(i));
            const auto& delta = collocation_directions[i];
            CHECK(dual.face == delta.face);
            if (!dual.free_query || !delta.free_query)
                FAIL("Free collocation coordinates must have covectors");
            implicit_reverse.add(*dual.free_query, *delta.free_query);
            implicit_reverse.add(dual.tangents[0], delta.tangents[0]);
            implicit_reverse.add(dual.tangents[1], delta.tangents[1]);
            implicit_reverse.add(dual.velocity, delta.velocity);
        }
        const auto work_difference = abs(implicit_forward.value - implicit_reverse.value);
        CHECK(work_difference <=
              implicit_forward.radius + implicit_reverse.radius + rational_word(5e-15));
        CHECK(implicit_reverse.radius < rational_word(1e-24));
        const auto independently_solved =
            rational_word(has_wake ? -.0552599357049885265356 : -.0753320247028545422098);
        CHECK(abs(implicit_forward.value - independently_solved) <=
              implicit_forward.radius + rational_word(5e-15));
        CHECK(abs(implicit_reverse.value - independently_solved) <=
              implicit_reverse.radius + rational_word(5e-15));
        CHECK(abs(implicit_reverse.value - reverse.value) > rational_word(1e-4));
        reject(combined.collocation(count), core::FailureCode::InvalidInput);
        reject(combined.circulation_rate(count), core::FailureCode::InvalidInput);
    }
}

TEST_CASE("Discrete nonpenetration matches an owned decimal event on its requested side",
          "[AE-07][UW-D03][N4-06][verification]") {
    auto input_case = source(
        R"({"kind":"pchip","times":[0,0.1,0.2],"values":[[0,0,0.508],[0,0,0.508],[0,0,0.508]],"outside":{"kind":"reject"}})",
        R"([{"start":0,"end":0.1,"step":0.1,"restart":"initial"},{"start":0.1,"end":0.2,"step":0.1,"restart":"event"}])");
    input_case.field("units").set("time", string("s"));
    input_case.field("units").set("length", string("m"));
    a::UniformFlowWorkspace flow_workspace;
    const auto provider = import(input_case, io::SourceId{103}, flow_workspace);
    const auto current = import(input_case, io::SourceId{104}, flow_workspace);
    const Number zero{0}, one{1}, absolute{1, 1000000000}, event{1, 10}, remote{1, 5};
    const cf::CaseFlowComponentBudget component{one.view(), absolute.view(), zero.view()};
    const std::array budgets{component, component, component};
    const cf::CaseFlowProjectionRequest projection{
        core::Precision::TwoWord, budgets, budgets, {4096, budget}};
    const auto project = [&](p::RationalView time, p::RequestedSide side) {
        core::ResourcePlan resources{budget};
        auto evaluation = take(cf::evaluate_case_uniform_flow(
            provider, current, {stamp, time, side}, flow_workspace, limits.flow, resources));
        return take(cf::project_case_uniform_flow(std::move(evaluation), projection));
    };
    constexpr auto words = [](std::uint64_t value) noexcept {
        return std::array{static_cast<std::uint32_t>(value),
                          static_cast<std::uint32_t>(value >> 32U)};
    };
    constexpr auto next_numerator = words(UINT64_C(100000000000000001));
    constexpr auto next_denominator = words(UINT64_C(1000000000000000000));
    const p::RationalView nearby{1, next_numerator, next_denominator};
    const ResidualFixture fixture;
    core::ResourcePlan memory{128ULL * 1024 * 1024};
    auto geometry_work = test_geometry::workspace(memory);
    const auto prepare = [&](p::RationalView exact_time) {
        auto input = fixture.input();
        input.context.time = 0.1;
        input.exact_time = exact_time;
        return a::prepare_owned_geometry(input, test_geometry::limits(), *geometry_work, memory);
    };
    const auto geometry = [&] {
        Number original{1, 10};
        auto captured = test_geometry::take(prepare(original.view()));
        original.n[0] = 7;
        original.d[0] = 9;
        const auto copied = test_geometry::take(captured.exact_time());
        if (!copied)
            FAIL("The geometry did not copy the original time");
        exact(*copied, event);
        return captured;
    }();
    const auto retained = test_geometry::take(geometry.exact_time());
    if (!retained)
        FAIL("The geometry lost its original event time");
    exact(*retained, event);

    const auto alternative = test_geometry::take(prepare(nearby));
    CHECK(geometry.context() == alternative.context());
    CHECK(test_geometry::take(geometry.identity()) != test_geometry::take(alternative.identity()));
    const auto low_work = [&] {
        auto context = test_geometry::context();
        context.time = 0.1;
        a::GeometryInput empty{{},
                               {},
                               a::SinglePositions{{}},
                               a::SinglePositions{{}},
                               {},
                               {},
                               {},
                               {},
                               {},
                               a::SingleWord{.125},
                               context};
        auto cap = test_geometry::limits();
        cap.max_steps = 200;
        CHECK(test_geometry::take(a::prepare_owned_geometry(empty, cap, *geometry_work, memory))
                  .active());
        empty.exact_time = event.view();
        return a::prepare_owned_geometry(empty, cap, *geometry_work, memory);
    }();
    if (low_work)
        FAIL("The geometry exceeded its exact-time work cap");
    CHECK(low_work.error().cause.code == core::FailureCode::ResourceExceeded);
    const auto inconsistent = prepare(remote.view());
    if (inconsistent)
        FAIL("The geometry accepted a disagreeing binary64 time");
    CHECK(inconsistent.error().cause.code == core::FailureCode::InvalidInput);

    const a::GlobalLimits cap{std::numeric_limits<std::uint64_t>::max() / 4,
                              numerics::linear::EstimateAdmission{4096}};
    auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
    // The reference model can be reused at another exact current event.
    const auto topology = take_global(
        a::build_global_topology(alternative, fixture.topology(), cap, *global, memory));
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, cap, *global, memory));
    const auto query = test_geometry::take(
        a::resolve_query(geometry, {a::QueryId{29}, a::FreePoint{a::FreeSingle{{.75, .5, 0}}}},
                         geometry.context(), test_geometry::limits(), *geometry_work, memory));
    const std::array collocation{
        cf::DiscreteCollocation{a::FacetId{200},
                                query,
                                {coordinates(1, 0, 0), coordinates(0, 2, 0)},
                                coordinates(.25, .25, .5),
                                a::DoubleWord{.75, 0x1p-60},
                                a::SingleWord{-.5}}};
    const auto old_id = a::P1StrengthId{0};
    const std::array wake{cf::DiscreteWakePanel{a::FacetId{202}, cf::AttachedWakeFront{}, old_id}};
    const std::array history{cf::FrozenWakeTrace{old_id, a::SingleWord{-.25}}};
    const cf::DiscreteNonpenetrationInput view{contact,
                                               collocation,
                                               wake,
                                               {core::AcceptedStateId{2}, history},
                                               {7, 11},
                                               stamp,
                                               p::RequestedSide::Right};
    const cf::DiscreteNonpenetrationLimits residual_limits{
        {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10},
        {test_geometry::limits(),
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(6000000000)};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const auto run = [&](cf::ProjectedCaseUniformFlow& flow) {
        return cf::evaluate_discrete_nonpenetration(view, std::move(flow), residual_limits, shared,
                                                    *local, *kernel, memory);
    };
    auto wrong_time = project(nearby, p::RequestedSide::Right);
    reject(run(wrong_time), core::FailureCode::StaleGeneration);
    CHECK(wrong_time.source().jet.active());
    auto wrong_side = project(event.view(), p::RequestedSide::Left);
    reject(run(wrong_side), core::FailureCode::StaleGeneration);
    CHECK(wrong_side.source().jet.active());

    auto right = project(event.view(), p::RequestedSide::Right);
    const auto result = take(run(right));
    exact(take(result.flow().source().jet.time()), event);
    CHECK(result.flow().source().requested_side == p::RequestedSide::Right);
    const auto row = take(result.row(0));
    const Rat denominator{Big{"10000000000000000000000000000000000000000"}};
    residual_expected(row.residual,
                      -Rat{Big{"5671057487729860677895669210567132495399"}} / denominator);
    residual_expected(row.circulation_action,
                      Rat{Big{"3730450057081755951141985320248458351172"}} / denominator);
}

namespace {
struct DisjointCirculationPanels {
    std::array<a::PointId, 8> points{a::PointId{1}, a::PointId{2}, a::PointId{3}, a::PointId{4},
                                     a::PointId{5}, a::PointId{6}, a::PointId{7}, a::PointId{8}};
    std::array<a::PointBinding, 8> bindings{test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding()};
    std::array<double, 24> xyz{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0,
                               0, 2, 0, 0, 3, 0, 1, 3, 0, 1, 2, 0};
    std::array<a::BandId, 2> bands{a::BandId{50}, a::BandId{51}};
    std::array<a::FaceInput, 2> faces{a::FaceInput{a::PhysicalFaceId{100}, bands[0]},
                                      a::FaceInput{a::PhysicalFaceId{101}, bands[1]}};
    std::array<a::FaceRoleInput, 2> roles{a::FaceRoleInput{faces[0].id, a::AeroFaceRole::Bound},
                                          a::FaceRoleInput{faces[1].id, a::AeroFaceRole::Bound}};
    std::array<a::EdgeInput, 10> edges{a::EdgeInput{a::EdgeId{10}, {points[0], points[1]}},
                                       a::EdgeInput{a::EdgeId{11}, {points[1], points[2]}},
                                       a::EdgeInput{a::EdgeId{12}, {points[2], points[3]}},
                                       a::EdgeInput{a::EdgeId{13}, {points[3], points[0]}},
                                       a::EdgeInput{a::EdgeId{14}, {points[0], points[2]}},
                                       a::EdgeInput{a::EdgeId{15}, {points[4], points[5]}},
                                       a::EdgeInput{a::EdgeId{16}, {points[5], points[6]}},
                                       a::EdgeInput{a::EdgeId{17}, {points[6], points[7]}},
                                       a::EdgeInput{a::EdgeId{18}, {points[7], points[4]}},
                                       a::EdgeInput{a::EdgeId{19}, {points[4], points[6]}}};
    std::array<a::EdgeBindingInput, 10> owners{
        a::EdgeBindingInput{edges[0].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[1].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[2].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[3].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[4].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[5].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[6].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[7].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[8].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[9].id, a::CanonicalEdge{}}};
    std::array<a::QuadInput, 2> quads{
        a::QuadInput{faces[0].id,
                     {points[0], points[1], points[2], points[3]},
                     {edges[0].id, edges[1].id, edges[2].id, edges[3].id},
                     edges[4].id,
                     {a::FacetId{200}, a::FacetId{201}}},
        a::QuadInput{faces[1].id,
                     {points[4], points[5], points[6], points[7]},
                     {edges[5].id, edges[6].id, edges[7].id, edges[8].id},
                     edges[9].id,
                     {a::FacetId{202}, a::FacetId{203}}}};

    [[nodiscard]] a::GeometryInput input() const noexcept {
        const a::SinglePositions positions{xyz};
        return {points,
                bindings,
                positions,
                positions,
                edges,
                faces,
                bands,
                {},
                quads,
                a::SingleWord{.125},
                test_geometry::context()};
    }
};
} // namespace

TEST_CASE("Discrete circulation and its linearization solve nonsymmetric mutual influences",
          "[UW-D01][UW-D06][N4-05][N4-06][verification][negative]") {
    const auto input_case = source(R"({"kind":"constant","value":[0,0,100]})", one_step);
    a::UniformFlowWorkspace flow_workspace;
    const auto provider = import(input_case, io::SourceId{110}, flow_workspace);
    const auto current = import(input_case, io::SourceId{111}, flow_workspace);
    const Number one{1}, zero{0}, absolute{1, 1000000000}, time{5, 4};
    const cf::CaseFlowComponentBudget budget_component{one.view(), absolute.view(), zero.view()};
    const std::array budgets{budget_component, budget_component, budget_component};
    const cf::CaseFlowProjectionRequest projection{
        core::Precision::TwoWord, budgets, budgets, {4096, budget}};
    auto exact_flow = evaluate(provider, current, flow_workspace, time);
    auto flow = take(cf::project_case_uniform_flow(std::move(exact_flow), projection));
    const DisjointCirculationPanels fixture;
    const a::GlobalLimits cap{std::numeric_limits<std::uint64_t>::max() / 4,
                              numerics::linear::EstimateAdmission{4096}};
    core::ResourcePlan memory{128ULL * 1024 * 1024};
    auto geometry_work = test_geometry::workspace(memory);
    const auto geometry = test_geometry::take(a::prepare_owned_geometry(
        fixture.input(), test_geometry::limits(), *geometry_work, memory));
    auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
    const auto topology = take_global(a::build_global_topology(
        geometry, {fixture.roles, fixture.owners, {}}, cap, *global, memory));
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, cap, *global, memory));
    const auto collocate = [&](a::QueryId id, std::array<double, 3> point) {
        return test_geometry::take(
            a::resolve_query(geometry, {id, a::FreePoint{a::FreeSingle{point}}}, geometry.context(),
                             test_geometry::limits(), *geometry_work, memory));
    };
    const auto first = collocate(a::QueryId{1}, {.75, .5, 0.});
    const auto second = collocate(a::QueryId{2}, {.625, 2.375, 0.});
    const auto rows = [&](const a::ResolvedQuery& second_query) {
        const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
        // The input may arrive in a different order than the qualified physical face order.
        return std::array{
            cf::DiscreteCollocation{a::FacetId{202}, second_query, frame, coordinates(0, 0, .75),
                                    a::SingleWord{9}, a::SingleWord{0}},
            cf::DiscreteCollocation{a::FacetId{200}, first, frame, coordinates(0, 0, .5),
                                    a::SingleWord{7}, a::SingleWord{0}}};
    };
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const cf::DiscreteNonpenetrationLimits solve_limits{
        {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10},
        {test_geometry::limits(),
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(24000000000)};
    const std::array unknown_scales{.5, 2.}, residual_scales{1., .125};
    const std::array<ls::UnitDimensions, 2> gamma_units{
        {{2, 0, -1, 0, 0, 0, 0}, {2, 0, -1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 2> velocity_units{
        {{1, 0, -1, 0, 0, 0, 0}, {1, 0, -1, 0, 0, 0, 0}}};
    const ls::Scales scales{unknown_scales, residual_scales, gamma_units, velocity_units,
                            circulation_context(geometry.context())};
    const auto request = circulation_request(scales, 2, ls::Norm::L2);
    const auto solve = [&](const a::ResolvedQuery& query) {
        const auto collocations = rows(query);
        const cf::DiscreteNonpenetrationInput input{contact,
                                                    collocations,
                                                    {},
                                                    {core::AcceptedStateId{2}, {}},
                                                    {7, 11},
                                                    stamp,
                                                    p::RequestedSide::Automatic};
        return cf::solve_discrete_circulation(input, std::move(flow), solve_limits, request, shared,
                                              *local, *kernel, memory);
    };
    reject(solve(first), core::FailureCode::RankUnresolved);
    CHECK(flow.source().jet.active());
    auto solved = take(solve(second));
    const Rat denominator{Big{"1000000000000000000000000000000000000"}};
    const std::array reference{Rat{Big{"5487032447771171614219353426944195"}} / denominator,
                               -Rat{Big{"274509732948309252476622637390296804"}} / denominator};
    for (std::size_t i = 0; i < reference.size(); ++i) {
        CAPTURE(i);
        CHECK(abs(solution_word(solved.solution(), i) - reference[i]) <= rational_word(1e-14));
        const auto original = take(solved.row(i));
        CHECK(original.face == fixture.faces[i].id);
        residual_small(original.original_residual, rational_word(1e-13) * residual_scales[i]);
    }

    constexpr std::uint64_t linearization_steps = 100000;
    std::array forcing{
        cf::CirculationLinearRhs{fixture.faces[0].id, a::DoubleWord{-.125, 0x1p-57}},
        cf::CirculationLinearRhs{fixture.faces[1].id, a::DoubleWord{.1875, -0x1p-58}}};
    const std::array seed{
        cf::CirculationLinearRhs{fixture.faces[0].id, a::DoubleWord{.25, 0x1p-57}},
        cf::CirculationLinearRhs{fixture.faces[1].id, a::DoubleWord{-.15625, -0x1p-59}}};
    const auto tangent_sense = cf::CirculationSense::Tangent;
    const auto transpose_sense = cf::CirculationSense::Transpose;
    auto tangent_request = request;
    tangent_request.rhs_generation = core::RhsGeneration{30};
    tangent_request.retained_generation = core::NumericalGeneration{31};

    const std::array<ls::UnitDimensions, 2> dual_speed{
        {{-1, 0, 1, 0, 0, 0, 0}, {-1, 0, 1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 2> dual_gamma{
        {{-2, 0, 1, 0, 0, 0, 0}, {-2, 0, 1, 0, 0, 0, 0}}};
    auto dual_context = scales.context;
    dual_context.cache.matrix = core::NumericalGeneration{32};
    dual_context.cache.scales = core::ScaleGeneration{33};
    dual_context.unknown_layout = core::ScopeId{34};
    dual_context.residual_layout = core::ScopeId{35};
    const ls::Scales dual_scales{residual_scales, unknown_scales, dual_speed, dual_gamma,
                                 dual_context};
    auto transpose_request = circulation_request(dual_scales, 2, ls::Norm::L2);
    transpose_request.rhs_generation = core::RhsGeneration{36};
    transpose_request.retained_generation = core::NumericalGeneration{37};

    const auto linearize = [&](cf::CirculationSense sense,
                               std::span<const cf::CirculationLinearRhs> rhs,
                               const cf::DiscreteCirculationRequest& selected) {
        return cf::solve_discrete_circulation_linearization(solved, sense, rhs, selected,
                                                            linearization_steps, memory);
    };
    reject(linearize(tangent_sense, std::span{forcing}.first(1), tangent_request),
           core::FailureCode::InvalidInput);
    const std::array reversed_forcing{forcing[1], forcing[0]};
    reject(linearize(tangent_sense, reversed_forcing, tangent_request),
           core::FailureCode::InvalidInput);
    auto bad_forcing = forcing;
    bad_forcing[0].value = a::SingleWord{std::numeric_limits<double>::quiet_NaN()};
    reject(linearize(tangent_sense, bad_forcing, tangent_request), core::FailureCode::Nonfinite);
    auto stale = transpose_request;
    stale.scales.context.cache.matrix = scales.context.cache.matrix;
    reject(linearize(transpose_sense, seed, stale), core::FailureCode::StaleGeneration);
    stale = transpose_request;
    stale.scales.context.cache.scales = scales.context.cache.scales;
    reject(linearize(transpose_sense, seed, stale), core::FailureCode::StaleGeneration);
    auto changed_unknown = unknown_scales;
    changed_unknown.front() *= 2;
    auto stale_tangent = tangent_request;
    stale_tangent.scales.unknown = changed_unknown;
    reject(linearize(tangent_sense, forcing, stale_tangent), core::FailureCode::StaleGeneration);
    auto changed_residual = residual_scales;
    changed_residual.back() *= 2;
    stale_tangent = tangent_request;
    stale_tangent.scales.residual = changed_residual;
    reject(linearize(tangent_sense, forcing, stale_tangent), core::FailureCode::StaleGeneration);
    stale_tangent = tangent_request;
    stale_tangent.scales.unknown_units = velocity_units;
    reject(linearize(tangent_sense, forcing, stale_tangent), core::FailureCode::StaleGeneration);

    auto fresh_tangent = tangent_request;
    fresh_tangent.scales.unknown = changed_unknown;
    fresh_tangent.scales.context.cache.scales = core::ScaleGeneration{38};
    const auto fresh = take(linearize(tangent_sense, forcing, fresh_tangent));
    CHECK(fresh.solution().certificate().context().cache.scales == core::ScaleGeneration{38});
    CHECK(std::ranges::equal(fresh.solution().certificate().unknown_scales(), changed_unknown));
    stale = transpose_request;
    stale.scales.context.trial = core::TrialId{99};
    reject(linearize(transpose_sense, seed, stale), core::FailureCode::StaleGeneration);
    for (const auto primal_layout :
         {scales.context.unknown_layout, scales.context.residual_layout}) {
        stale = transpose_request;
        stale.scales.context.unknown_layout = primal_layout;
        reject(linearize(transpose_sense, seed, stale), core::FailureCode::StaleGeneration);
        stale = transpose_request;
        stale.scales.context.residual_layout = primal_layout;
        reject(linearize(transpose_sense, seed, stale), core::FailureCode::StaleGeneration);
    }
    stale = transpose_request;
    stale.scales.context.residual_layout = stale.scales.context.unknown_layout;
    reject(linearize(transpose_sense, seed, stale), core::FailureCode::StaleGeneration);
    auto bad_units = transpose_request;
    bad_units.scales.unknown_units = velocity_units;
    reject(linearize(transpose_sense, seed, bad_units), core::FailureCode::InvalidInput);
    auto exact_rhs_projection = transpose_request;
    exact_rhs_projection.rhs_absolute = 0;
    reject(linearize(transpose_sense, seed, exact_rhs_projection),
           core::FailureCode::PrecisionUnresolved);
    auto exact_original = transpose_request;
    exact_original.original_absolute = 0;
    reject(linearize(transpose_sense, seed, exact_original), core::FailureCode::LinearUnconverged);
    // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
    std::stop_source cancellation;
    REQUIRE(cancellation.request_stop());
    core::ResourcePlan no_space{0};
    reject(cf::solve_discrete_circulation_linearization(solved, tangent_sense, forcing,
                                                        tangent_request, linearization_steps,
                                                        no_space, cancellation.get_token()),
           core::FailureCode::Cancelled);
    CHECK(no_space.allocated_bytes() == 0);
    reject(cf::solve_discrete_circulation_linearization(solved, tangent_sense, forcing,
                                                        tangent_request, 0, no_space),
           core::FailureCode::ResourceExceeded);
    reject(cf::solve_discrete_circulation_linearization(
               solved, tangent_sense, forcing, tangent_request, linearization_steps, no_space),
           core::FailureCode::ResourceExceeded);

    const auto tangent = take(linearize(tangent_sense, forcing, tangent_request));
    const auto transpose = take(linearize(transpose_sense, seed, transpose_request));
    CHECK(tangent.sense() == tangent_sense);
    CHECK(transpose.sense() == transpose_sense);
    CHECK(tangent.solution().values().precision() == core::Precision::TwoWord);
    CHECK(transpose.solution().values().precision() == core::Precision::TwoWord);
    CHECK(transpose.solution().certificate().context() == dual_context);
    CHECK(transpose.solution().certificate().rhs_generation() == core::RhsGeneration{36});
    CHECK(tangent.used_steps() <= linearization_steps);
    CHECK(transpose.used_steps() <= linearization_steps);
    forcing[0].value = a::SingleWord{500};
    const std::array tangent_reference{
        -Rat{Big{"131244355849310971713061064609083227"}} / denominator,
        Rat{Big{"210813958349156791484623892498037848"}} / denominator};
    const std::array transpose_reference{
        -Rat{Big{"264772325234835981543571156442296025"}} / denominator,
        Rat{Big{"174155889600154653987400703840364351"}} / denominator};
    const std::array forcing_exact{-Rat{1, 8} + two_power(-57), Rat{3, 16} - two_power(-58)};
    const std::array seed_exact{Rat{1, 4} + two_power(-57), -Rat{5, 32} - two_power(-59)};
    Rat dual_error{0};
    for (std::size_t i = 0; i < forcing.size(); ++i) {
        CAPTURE(i);
        const auto t = solution_word(tangent.solution(), i);
        const auto l = solution_word(transpose.solution(), i);
        CHECK(abs(t - tangent_reference[i]) <= rational_word(1e-14));
        CHECK(abs(l - transpose_reference[i]) <= rational_word(1e-14));
        dual_error += seed_exact[i] * t + forcing_exact[i] * l;
        const auto tangent_row = take(tangent.row(i));
        const auto transpose_row = take(transpose.row(i));
        CHECK(tangent_row.face == fixture.faces[i].id);
        CHECK(transpose_row.face == fixture.faces[i].id);
        residual_small(tangent_row.original_residual, rational_word(1e-13) * residual_scales[i]);
        residual_small(transpose_row.original_residual, rational_word(1e-13) * unknown_scales[i]);
        CHECK(abs(solution_word(solved.solution(), i) - reference[i]) <= rational_word(1e-14));
    }
    CHECK(abs(dual_error) <= rational_word(1e-16));
    reject(tangent.row(forcing.size()), core::FailureCode::InvalidInput);
    CHECK(solved.flow().source().jet.active());
    const std::array zero_rhs{cf::CirculationLinearRhs{fixture.faces[0].id, a::SingleWord{0}},
                              cf::CirculationLinearRhs{fixture.faces[1].id, a::SingleWord{0}}};
    auto exact_zero = tangent_request;
    exact_zero.rhs_absolute = 0;
    exact_zero.original_absolute = 0;
    const auto zero_tangent = take(linearize(tangent_sense, zero_rhs, exact_zero));
    for (std::size_t i = 0; i < zero_rhs.size(); ++i) {
        CHECK(solution_word(zero_tangent.solution(), i) == Rat{0});
        residual_small(take(zero_tangent.row(i)).original_residual, Rat{0});
    }

    // A uniform-velocity direction is shared across these two nonsymmetric rows.
    const auto two_collocations = rows(second);
    const cf::DiscreteNonpenetrationInput pullback_input{
        contact, two_collocations,           {}, {core::AcceptedStateId{2}, {}}, {7, 11},
        stamp,   p::RequestedSide::Automatic};
    auto exact_objective_flow = evaluate(provider, current, flow_workspace, time);
    auto objective_flow =
        take(cf::project_case_uniform_flow(std::move(exact_objective_flow), projection));
    // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    core::ResourcePlan no_pullback_bytes{0};
    const auto pullback = [&](core::ResourcePlan& arena, const std::stop_token& cancel = {}) {
        return cf::solve_discrete_circulation_pullback(
            pullback_input, seed, std::move(objective_flow), solve_limits, request,
            transpose_request, shared, *local, *kernel, arena, cancel);
    };
    reject(pullback(no_pullback_bytes, stop.get_token()), core::FailureCode::Cancelled);
    CHECK(no_pullback_bytes.allocated_bytes() == 0);
    CHECK(objective_flow.source().jet.active());
    const auto dual_pullback = take(pullback(memory));
    const auto lambda_a = solution_word(transpose.solution(), 0);
    const auto lambda_b = solution_word(transpose.solution(), 1);
    const auto flow_dual = dual_pullback.uniform_flow();
    residual_expected(flow_dual[2], -lambda_a - lambda_b);
    residual_expected(take(dual_pullback.row(0)).velocity[2], lambda_a);
    residual_expected(take(dual_pullback.row(1)).velocity[2], lambda_b);

    const auto relocated = std::move(solved);
    reject(linearize(tangent_sense, {}, tangent_request), core::FailureCode::InvalidInput);
    CHECK(relocated.flow().source().jet.active());
    CHECK(abs(solution_word(tangent.solution(), 0) - tangent_reference[0]) <= rational_word(1e-14));
}
