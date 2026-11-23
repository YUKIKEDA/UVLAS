#include "coupling/case_uniform_flow.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <limits>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
namespace a = aero::detail;
namespace p = a::pchip;
namespace cf = coupling::detail;

constexpr std::uint64_t budget = std::uint64_t{1} << 48U;
const p::PchipInputStamp stamp{p::PchipIssuerId{71}, p::PchipInputGeneration{8}};

structure::ExactLimits exact_limits() {
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

TEST_CASE("Typed SI PCHIP uses the declared closed time window and independent current case",
          "[AE-07][IO4-01][N4-02][N4-03][N4-08][N4-09][property][negative]") {
    auto input = source(
        R"({"kind":"pchip","times":[0,1,2],"values":[[0,-2,1],[1,-2,-1],[4,-2,-7]],"outside":{"kind":"reject"}})",
        two_step);
    input.field("flow").set("events", parse("[1]"));
    a::UniformFlowWorkspace work;
    const auto field = import(input, io::SourceId{7}, work);
    const auto current = import(input, io::SourceId{8}, work);
    const auto converted = take(field.uniform().input());
    const auto* pchip = std::get_if<a::PchipUniformFlow>(&converted.series);
    if (!pchip || converted.events.size() != 1)
        FAIL("Expected typed PCHIP with one retained event");
    exact(pchip->boundary.begin, Number{0});
    exact(pchip->boundary.end, Number{120});
    exact(converted.events[0], Number{60});

    core::ResourcePlan memory{budget};
    const auto saved = take(cf::snapshot_case_uniform_flow(field, work, limits.flow, memory));
    const auto saved_view = take(saved.view());
    const auto replay =
        take(cf::restore_case_uniform_flow(saved_view, current, work, limits.flow, memory));
    // In source minutes x(t)=3t²/2-t³/2 on [0,1]. The third entry is 1-2x.
    const Number time{30};
    const auto result = evaluate(replay, current, work, time);
    jet(result.jet, {Number{127, 80000}, Number{127, 12500, -1}, Number{381, 200000}},
        {Number{381, 4000000}, Number{0}, Number{381, 2000000, -1}}, time);

    const Number knot{60};
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
