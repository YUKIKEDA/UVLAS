#include "aero/geometry_snapshot_fixture.hpp"
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

void within(const io::SiEnclosure& enclosure, const Rat& actual, const Rat& allowed) {
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
            CHECK(abs(c * two_power(result.exponent2) - actual) <= radius);
            CHECK(radius <= allowed);
        },
        enclosure);
}
} // namespace

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

void residual_expected(const a::ScaledEnclosure& value, const Rat& oracle) {
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
    const Rat comparison_uncertainty{Big{1}, Big{"1000000000000000000000000"}};
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
    stale.scales.context.trial = core::TrialId{99};
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

    auto relocated = std::move(solved);
    reject(linearize(tangent_sense, {}, tangent_request), core::FailureCode::InvalidInput);
    CHECK(relocated.flow().source().jet.active());
    CHECK(abs(solution_word(tangent.solution(), 0) - tangent_reference[0]) <= rational_word(1e-14));
}
