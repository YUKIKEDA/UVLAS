#include "aero/geometry_snapshot_fixture.hpp"
#include "aero/query_coordinates_private.hpp"
#include "coupling/case_aero_attachment.hpp"
#include "coupling/case_flow_clock.hpp"
#include "coupling/case_known_flow.hpp"
#include "coupling/case_reference_graph.hpp"
#include "coupling/case_uniform_flow.hpp"
#include "coupling/case_uniform_flow_projection.hpp"
#include "coupling/discrete_nonpenetration.hpp"
#include "coupling/fixed_bdf2_circulation.hpp"
#include "coupling/fixed_wake_transport.hpp"
#include "coupling/pg_aero_attachment.hpp"
#include "coupling/wake_strength_transition.hpp"
#include "io/case_schema_test_adapter.hpp"
#include "io/si_rational_projection.hpp"
#include "structure/pg_conditioning_test_support.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <source_location>
#include <uvlas/core/instrumentation.hpp>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <vector>

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
T take(core::Result<T> value, const std::source_location at = std::source_location::current()) {
    if (!value) {
        INFO(at.file_name() << ':' << at.line());
        if (value.error().details.resources) {
            const auto& usage = *value.error().details.resources;
            INFO("Resource kind " << static_cast<int>(usage.kind) << ", required " << usage.required
                                  << ", observed " << usage.observed);
            FAIL("Unexpected case-flow failure " << static_cast<int>(value.error().code));
        }
        FAIL("Unexpected case-flow failure " << static_cast<int>(value.error().code));
    }
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

TEST_CASE("Typed thin surface reference owns exact Cartesian position and tangents",
          "[verification][negative][AE-02][GEO-01][FD-06]") {
    auto make = [] {
        auto value = root("aerodynamic_steady");
        value.field("units").set("length", string("mm"));
        auto& geometries = value.field("model").field("geometry").elements;
        geometries.at(0) = parse(R"({
            "kind":"surface_graph","id":"1",
            "frame":{"origin":[4000,-2000,3000],"rotation":[[0,-1,0],[1,0,0],[0,0,1]]},
            "reference_length":2000,
            "x0":{"numerator":"1","denominator":"3"},
            "width":{"numerator":"3","denominator":"2"},
            "skew":{"numerator":"-1","denominator":"4"},
            "degree":[2,3],
            "coefficients":[
                {"powers":[2,1],"coefficient":{"numerator":"1","denominator":"2"}},
                {"powers":[0,3],"coefficient":{"numerator":"-1","denominator":"3"}}],
            "domain":[
                [{"numerator":"0","denominator":"1"},{"numerator":"0","denominator":"1"}],
                [{"numerator":"1","denominator":"1"},{"numerator":"0","denominator":"1"}],
                [{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"1"}],
                [{"numerator":"1","denominator":"2"},{"numerator":"1","denominator":"2"}],
                [{"numerator":"0","denominator":"1"},{"numerator":"1","denominator":"1"}]],
            "boundaries":[],"continuity":[]
        })");
        return value;
    };
    auto input = make();
    auto& graph = input.field("model").field("geometry").elements.at(0);
    // The aerodynamic surface references leading/trailing rows of its own graph.
    auto& boundaries = graph.field("boundaries").elements;
    const auto example = root("aerodynamic_steady");
    boundaries =
        example.field("model").field("geometry").elements.at(0).field("boundaries").elements;
    constexpr io::SourceId source_id{166};
    const structure::PgReferenceLimits bounds{exact_limits(), 10000, 1000, {}};
    core::ResourcePlan memory{budget};
    const auto ready = [&] {
        const auto si = materialize(input, source_id);
        reject(cf::CaseReferenceSurface::prepare(si, 90, bounds, memory),
               core::FailureCode::InvalidInput);
        return take(cf::CaseReferenceSurface::prepare(si, 1, bounds, memory));
    }();
    const auto identity = take(ready.identity());
    CHECK(identity.geometry_id == 1);
    CHECK(identity.source.source.source == source_id);
    const auto scalar = [&](std::string_view n, std::string_view d = "1") {
        return take(structure::ExactRational::from_integer_pair(n, d, bounds.algebra, memory));
    };
    const auto chord = scalar("1", "4");
    const auto span = scalar("1", "2");
    const auto result = take(ready.evaluate(chord, span, bounds, memory));
    // Independent substitution into X=(4-2η, -2+2(1/3+3ξ/2-η/4),
    // 3+ξ²η-2η³/3), with both analytic derivatives in the same coordinates.
    CHECK(result.position[0] == scalar("3"));
    CHECK(result.position[1] == scalar("-5", "6"));
    CHECK(result.position[2] == scalar("283", "96"));
    CHECK(result.chord_tangent[0] == scalar("0"));
    CHECK(result.chord_tangent[1] == scalar("3"));
    CHECK(result.chord_tangent[2] == scalar("1", "4"));
    CHECK(result.span_tangent[0] == scalar("-2"));
    CHECK(result.span_tangent[1] == scalar("-1", "2"));
    CHECK(result.span_tangent[2] == scalar("-7", "16"));

    const auto upper = scalar("3", "4");
    const auto edge = take(ready.evaluate(chord, upper, bounds, memory));
    CHECK(edge.position[0] == scalar("5", "2"));
    reject(ready.evaluate(span, upper, bounds, memory), core::FailureCode::GeometryInvalid);
    const structure::PgReferenceLimits no_visits{exact_limits(), 10000, 0, {}};
    reject(ready.evaluate(chord, span, no_visits, memory), core::FailureCode::ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() requires a mutable source.
    std::stop_source stop;
    static_cast<void>(stop.request_stop());
    const structure::PgReferenceLimits cancelled{exact_limits(), 10000, 1000, stop.get_token()};
    reject(ready.evaluate(chord, span, cancelled, memory), core::FailureCode::Cancelled);

    graph.field("frame").set("rotation", parse("[[0,1,0],[1,0,0],[0,0,1]]"));
    const auto mirrored = materialize(input, source_id);
    reject(cf::CaseReferenceSurface::prepare(mirrored, 1, bounds, memory),
           core::FailureCode::GeometryInvalid);
    CHECK(take(ready.evaluate(chord, span, bounds, memory)).position[2] == scalar("283", "96"));
}

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
    const std::array<cf::CaseFlowComponentBudget, 3> vector{per_component, per_component,
                                                            per_component};
    const std::array<cf::CaseFlowComponentBudget, 9> matrix{
        per_component, per_component, per_component, per_component, per_component,
        per_component, per_component, per_component, per_component};
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

TEST_CASE("Scalar gust preparation caps imported SI recipes independently of provider padding",
          "[AE-07][IO4-01][N4-02][verification][negative]") {
    constexpr std::string_view pchip =
        R"({"kind":"pchip","times":[0,0.5,1],"values":[[1],[2],[5]],"outside":{"kind":"reject"}})";
    constexpr std::string_view polynomial =
        R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[1]]})";
    // PCHIP imports 3 times, 3 amplitudes and 2 analysis endpoints. Polynomial
    // imports origin, scale, domain, 2 coverage endpoints and 1 coefficient.
    for (const auto& [series, recipes] :
         std::array{std::pair{pchip, 8U}, std::pair{polynomial, 7U}}) {
        CAPTURE(recipes);
        auto input = gust_case("rad", series);
        input.field("flow").set("events", parse("[]"));
        const auto si = materialize(input, io::SourceId{60});
        a::UniformFlowWorkspace work;
        core::ResourcePlan memory{budget};
        auto exact = gust_limits;
        exact.time_curves.maximum_exact_numbers = recipes;
        CHECK(cf::prepare_case_known_flow(si, stamp, {}, work, exact, memory));
        --exact.time_curves.maximum_exact_numbers;
        reject(cf::prepare_case_known_flow(si, stamp, {}, work, exact, memory),
               core::FailureCode::ResourceExceeded);
    }
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

TEST_CASE("The exact transient stage evaluates known flow and prescribed wake together",
          "[AE-07][UW-D02][IO4-01][N4-03][N4-06][verification][negative]") {
    auto input = gust_case(
        "deg",
        R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[2],[3]]})");
    input.field("analysis")
        .set(
            "intervals",
            parse(
                R"([{"start":0,"end":0.5,"step":0.25,"restart":"initial"},{"start":0.5,"end":1,"step":0.5,"restart":"event"}])"));
    input.field("flow").set("events", parse("[0.5]"));
    input.field("flow").field("gust").set("phase", parse("90"));
    input.field("analysis")
        .field("wake")
        .set(
            "convection",
            parse(
                R"({"kind":"prescribed","velocity":{"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[4,1,0],[2,0,0]]}})"));
    a::UniformFlowWorkspace work;
    const cf::CaseClockedKnownFlowLimits clock_limits{{limits, 2, 4}, gust_limits};
    const auto make = [&](const JsonToken& fixture, io::SourceId id) {
        const auto si = materialize(fixture, id);
        core::ResourcePlan memory{budget};
        return take(cf::prepare_case_clocked_known_flows(si, stamp, work, clock_limits, memory));
    };
    const auto field = make(input, io::SourceId{420});
    const auto current = make(input, io::SourceId{421});
    CHECK(field.interval_count() == 2);
    CHECK(take(field.interval_steps(0)) == 2);
    CHECK(take(field.interval_steps(1)) == 1);

    const Number one{1}, zero{0}, tolerance{1, 1000000000}, middle_x{3, 4};
    const cf::CaseFlowComponentBudget allowed{one.view(), tolerance.view(), zero.view()};
    const std::array components{allowed, allowed, allowed};
    const cf::CaseFlowProjectionRequest wake{
        core::Precision::TwoWord, components, components, {4096, budget}};
    const auto ambient = gust_projection(core::Precision::TwoWord, tolerance);
    const auto run = [&](const cf::CaseClockedKnownFlows& host, cf::CaseFlowClockQuery address,
                         const Number& x) {
        core::ResourcePlan memory{budget};
        const std::array position{x.view(), zero.view(), zero.view()};
        return cf::evaluate_case_clocked_known_flows(field, host, address, position, ambient, wake,
                                                     work, clock_limits, memory);
    };
    const auto check = [&](const cf::CaseClockedKnownFlowSample& sample, const Number& time,
                           const Number& step, const Number& x, p::RequestedSide side,
                           const Rat& expected_u, const Rat& expected_wake) {
        exact(take(sample.time()), time);
        exact(take(sample.ambient().time()), time);
        exact(take(sample.prescribed().jet().time()), time);
        exact(take(sample.step_seconds()), step);
        exact(take(sample.ambient().position(0)), x);
        CHECK(sample.source().source.source == io::SourceId{420});
        CHECK(sample.address().side == side);
        CHECK(sample.ambient().requested_side() == side);
        CHECK(sample.prescribed().requested_side() == side);
        const Rat radius = two_power(-98);
        within(take(sample.ambient().velocity(1)), expected_u, radius);
        within(take(sample.ambient().time_partial(1)), Rat{9} / 5, radius);
        within(take(sample.ambient().spatial_partial(1, 0)), Rat{0}, Rat{0});
        within(take(sample.prescribed().velocity(0)), expected_wake, Rat{0});
        within(take(sample.prescribed().time_partial(0)), Rat{2}, Rat{0});
    };
    using enum p::RequestedSide;
    // The spatial and temporal phase corrections are -22.5 and +22.5 degrees.
    const auto middle = take(run(current, {stamp, 0, 1, Automatic}, middle_x));
    check(middle, Number{1, 4}, Number{1, 4}, middle_x, Automatic, Rat{73} / 20, Rat{9} / 2);
    const auto left = take(run(current, {stamp, 0, 2, Left}, one));
    check(left, Number{1, 2}, Number{1, 4}, one, Left, Rat{41} / 10, Rat{5});
    const auto right = take(run(current, {stamp, 1, 0, Right}, one));
    check(right, Number{1, 2}, Number{1, 2}, one, Right, Rat{41} / 10, Rat{5});
    CHECK(left.restart() == cf::CaseClockRestart::ItemInitial);
    CHECK(right.restart() == cf::CaseClockRestart::ItemEvent);

    reject(run(current, {stamp, 1, 0, Left}, one), core::FailureCode::InvalidInput);
    auto changed_grid = input;
    changed_grid.field("analysis").field("intervals").elements[1].set("step", parse("0.25"));
    const auto other_grid = make(changed_grid, io::SourceId{422});
    reject(run(other_grid, {stamp, 0, 1, Automatic}, middle_x), core::FailureCode::StaleGeneration);
    auto changed_gust = input;
    changed_gust.field("flow").field("gust").set("phase", parse("90.000000000000000001"));
    const auto other_gust = make(changed_gust, io::SourceId{423});
    reject(run(other_gust, {stamp, 0, 1, Automatic}, middle_x), core::FailureCode::StaleGeneration);
    auto changed_wake = input;
    changed_wake.field("analysis")
        .field("wake")
        .field("convection")
        .field("velocity")
        .field("coefficients")
        .elements[1]
        .elements[0] = parse("3");
    const auto other_wake = make(changed_wake, io::SourceId{424});
    reject(run(other_wake, {stamp, 0, 1, Automatic}, middle_x), core::FailureCode::StaleGeneration);
}

TEST_CASE("The known-flow clock checks periodic model events without a prescribed attachment",
          "[AE-07][BC-01][IO4-01][N4-06][verification][negative]") {
    a::UniformFlowWorkspace work;
    const cf::CaseClockedKnownFlowLimits clock_limits{{limits, 2, 4}, gust_limits};
    const auto make = [&](const JsonToken& fixture) {
        const auto si = materialize(fixture, io::SourceId{425});
        core::ResourcePlan memory{budget};
        return cf::prepare_case_clocked_known_flows(si, stamp, work, clock_limits, memory);
    };
    const auto intervals = [](bool split) {
        return parse(
            split
                ? R"([{"start":-3.5,"end":-3,"step":0.5,"restart":"initial"},{"start":-3,"end":-2.5,"step":0.5,"restart":"event"}])"
                : R"([{"start":-3.5,"end":-2.5,"step":1,"restart":"initial"}])");
    };
    const auto base = [&] {
        auto fixture = source(R"({"kind":"constant","value":[0,0,0]})", "[]");
        fixture.field("units").set("time", string("s"));
        fixture.field("analysis")
            .field("wake")
            .set("convection",
                 parse(R"({"kind":"prescribed","velocity":{"kind":"constant","value":[0,0,0]}})"));
        fixture.field("analysis").set("intervals", intervals(false));
        return fixture;
    }();
    CHECK(take(make(base)).interval_count() == 1);
    constexpr std::string_view scalar =
        R"({"kind":"pchip","times":[0,1,4],"values":[[0],[1],[0]],"outside":{"kind":"periodic","period":4}})";
    constexpr std::string_view vector =
        R"({"kind":"pchip","times":[0,1,4],"values":[[0,0,0],[1,0,0],[0,0,0]],"outside":{"kind":"periodic","period":4}})";
    auto fixture = base;
    SECTION("A structural support supplies the only periodic series") {
        auto& support = fixture.field("model").field("supports").elements.front();
        support.field("motion").set("translation", parse(vector));
        support.field("motion").set("events", parse("[0,1,4]"));
    }
    SECTION("A shell load supplies the only periodic series") {
        auto load = fragment("Load.pressure");
        load.set("multiplier", parse(scalar));
        fixture.field("model").field("shells").elements.front().field("loads").elements.push_back(
            std::move(load));
    }
    SECTION("A beam load supplies the only periodic series") {
        auto load = fragment("Load.point_force");
        load.field("target").field("target").set("kind", string("beam"));
        load.set("series", parse(vector));
        auto beam = fragment("Beam");
        beam.field("loads").elements.push_back(std::move(load));
        fixture.field("model").field("beams").elements.push_back(std::move(beam));
        auto& scales = fixture.field("numerics").field("scales").elements;
        for (const auto field : {"beam_translation", "beam_rotation", "moment_balance",
                                 "structural_angular_velocity"}) {
            auto scale = fragment("Scale.uniform_components");
            scale.set("field", string(field));
            scale.set("targets", parse(R"([{"kind":"beam","id":"1"}])"));
            scales.push_back(std::move(scale));
        }
    }
    // The source knot at 1 s recurs at -3 s. A boundary there admits both sides.
    reject(make(fixture), core::FailureCode::InvalidInput);
    fixture.field("analysis").set("intervals", intervals(true));
    CHECK(take(make(fixture)).interval_count() == 2);
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
        std::get<cf::ExplicitCollocationTangents>(collocation[0].tangents)[1] =
            coordinates(0, -2, 0);
        reject(run(flow, residual_limits, memory), core::FailureCode::GeometryInvalid);
        std::get<cf::ExplicitCollocationTangents>(collocation[0].tangents)[1] =
            coordinates(0, 2, 0);
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
            solved, cf::CirculationSense::Tangent, tangent_rhs, request, solve_limits.steps,
            memory));
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
                original_row.facet, material,
                std::get<cf::ExplicitCollocationTangents>(original_row.tangents),
                original_row.velocity, original_row.circulation,
                original_row.circulation_direction}};
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
        cf::DiscreteVectorLoadInput load_input{
            input.trial_stamp,
            request.retained_generation,
            has_wake ? std::span<const cf::EdgeMidpointVelocity>{wake_force_edges}
                     : std::span<const cf::EdgeMidpointVelocity>{midpoint_velocities},
            std::span{rates}.first(count),
            a::SingleWord{1.25},
            moment_origin};
        if (mode == LoadCase::AliasedShared)
            load_input.retention = cf::DiscreteLoadRetention::OriginalRows;
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
        if (mode == LoadCase::AliasedShared) {
            REQUIRE(result.original_count() == 28);
            std::array<Rat, 3> original_sum{};
            std::array<std::array<Rat, 3>, 2> shared_edge_face{};
            std::array<std::size_t, 2> shared_edge_ends{};
            for (std::size_t i = 0; i < result.original_count(); ++i) {
                const auto contribution = take(result.original(i));
                REQUIRE(contribution.local_vertex < 3);
                const auto points = test_geometry::take(geometry.facet_points(contribution.facet));
                CHECK(test_geometry::take(geometry.owner(points[contribution.local_vertex])) ==
                      contribution.canonical);
                const auto panel =
                    contribution.face == fixture.faces[0].id ? std::size_t{0} : std::size_t{1};
                CHECK(contribution.face == fixture.faces[panel].id);
                if (contribution.edge == a::EdgeId{12})
                    ++shared_edge_ends[panel];
                for (std::size_t axis = 0; axis < original_sum.size(); ++axis) {
                    const auto actual = unpack(contribution.force[axis]);
                    original_sum[axis] += actual;
                    if (contribution.edge == a::EdgeId{12})
                        shared_edge_face[panel][axis] += actual;
                }
                if (!contribution.edge && contribution.facet == a::FacetId{200})
                    residual_expected(contribution.force[2], exact_force / Rat{3});
            }
            CHECK(shared_edge_ends == std::array<std::size_t, 2>{2, 2});
            // Each side has its own signed contribution before cancellation. The
            // sum is the independently frozen force on the physical shared edge.
            CHECK(shared_edge_face[0][0] * shared_edge_face[1][0] < 0);
            const auto shared_edge = take(result.edge(2)).force;
            for (std::size_t axis = 0; axis < original_sum.size(); ++axis) {
                CHECK(abs(original_sum[axis] - unpack(force[axis])) <= rational_word(1e-24));
                CHECK(abs(shared_edge_face[0][axis] + shared_edge_face[1][axis] -
                          unpack(shared_edge[axis])) <= rational_word(1e-24));
            }
        } else
            CHECK(result.original_count() == 0);
        reject(result.original(result.original_count()), core::FailureCode::InvalidInput);
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
    const auto rows = [&](const a::ResolvedQuery& first_query,
                          const a::ResolvedQuery& second_query) {
        const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
        // The input may arrive in a different order than the qualified physical face order.
        return std::array{
            cf::DiscreteCollocation{a::FacetId{202}, second_query, frame, coordinates(0, 0, .75),
                                    a::SingleWord{9}, a::SingleWord{0}},
            cf::DiscreteCollocation{a::FacetId{200}, first_query, frame, coordinates(0, 0, .5),
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
        const auto qualified = take_global(contact.retain_checked(memory));
        const auto fresh_query = [&](const a::ResolvedQuery& original) {
            return test_geometry::take(a::resolve_query(
                geometry, test_geometry::take(original.definition()), geometry.context(),
                test_geometry::limits(), *geometry_work, memory));
        };
        const auto owned_first = fresh_query(first);
        const auto owned_second = fresh_query(query);
        const auto collocations = rows(owned_first, owned_second);
        const cf::DiscreteNonpenetrationInput input{qualified,
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

    const auto linearization_steps = solve_limits.steps;
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

    // There are no live input row/query objects from solve() here. The old
    // algebra-only cap does not admit exact row replay from the owned input.
    reject(cf::solve_discrete_circulation_linearization(solved, tangent_sense, forcing,
                                                        tangent_request, 100000, memory),
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
    const auto two_collocations = rows(first, second);
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

TEST_CASE("Two discrete circulation rows consume their own known-flow values and derivatives",
          "[AE-07][UW-D01][UW-D06][N4-03][N4-05][N4-06][verification][negative]") {
    auto input_case = gust_case("deg", R"({"kind":"constant","value":[2]})");
    input_case.field("analysis").field("intervals").elements[0].set("end", parse("2"));
    input_case.field("flow").set("uniform", parse(R"({"kind":"constant","value":[0,0,1]})"));
    auto& gust = input_case.field("flow").field("gust");
    gust.set("direction", parse("[0,0,1]"));
    gust.set("wavevector", parse("[360,0,0]"));
    gust.set("phase", parse("0"));
    gust.set("origin", parse("[0.625,0,0]"));
    gust.set("time_origin", parse("1.25"));

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
    const auto resolve = [&](a::QueryId id, a::QueryDefinition expression) {
        return test_geometry::take(a::resolve_query(geometry, {id, expression}, geometry.context(),
                                                    test_geometry::limits(), *geometry_work,
                                                    memory));
    };
    const auto first = resolve(a::QueryId{1}, a::FreePoint{a::FreeSingle{{.75, .5, 0.}}});
    const auto second =
        resolve(a::QueryId{2},
                a::SourceFacePoint{a::FacetId{203}, a::SingleWord{.25}, a::SingleWord{.375}});
    const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
    const std::array collocations{
        cf::DiscreteCollocation{a::FacetId{202}, second, frame, coordinates(0, 0, .75),
                                a::SingleWord{9}, a::SingleWord{0}},
        cf::DiscreteCollocation{a::FacetId{200}, first, frame, coordinates(0, 0, .5),
                                a::SingleWord{7}, a::SingleWord{0}}};
    const cf::DiscreteNonpenetrationInput input{contact,
                                                collocations,
                                                {},
                                                {core::AcceptedStateId{2}, {}},
                                                {7, 11},
                                                stamp,
                                                p::RequestedSide::Automatic};

    a::UniformFlowWorkspace flow_workspace;
    const auto provider = prepare_gust(input_case, io::SourceId{441}, flow_workspace, memory);
    const auto current = prepare_gust(input_case, io::SourceId{442}, flow_workspace, memory);
    const Number time{5, 4}, absolute{1, 1000000000};
    cf::DiscreteKnownFlowInput flow{provider,
                                    current,
                                    {stamp, time.view()},
                                    gust_projection(core::Precision::TwoWord, absolute),
                                    gust_limits,
                                    test_geometry::limits(),
                                    *geometry_work,
                                    flow_workspace};
    const cf::DiscreteNonpenetrationLimits solve_limits{
        {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10},
        {test_geometry::limits(),
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(24000000000)};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const std::array unknown{.5, 2.}, residual{1., .125};
    const std::array<ls::UnitDimensions, 2> gamma{{{2, 0, -1, 0, 0, 0, 0}, {2, 0, -1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 2> speed{{{1, 0, -1, 0, 0, 0, 0}, {1, 0, -1, 0, 0, 0, 0}}};
    const ls::Scales scales{unknown, residual, gamma, speed,
                            circulation_context(geometry.context())};
    const auto primal_request = circulation_request(scales, 2, ls::Norm::L2);
    auto tangent_request = primal_request;
    tangent_request.rhs_generation = core::RhsGeneration{30};
    tangent_request.retained_generation = core::NumericalGeneration{31};
    const std::array<ls::UnitDimensions, 2> dual_gamma{
        {{-2, 0, 1, 0, 0, 0, 0}, {-2, 0, 1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 2> dual_speed{
        {{-1, 0, 1, 0, 0, 0, 0}, {-1, 0, 1, 0, 0, 0, 0}}};
    auto dual_context = scales.context;
    dual_context.cache.matrix = core::NumericalGeneration{32};
    dual_context.cache.scales = core::ScaleGeneration{33};
    dual_context.unknown_layout = core::ScopeId{34};
    dual_context.residual_layout = core::ScopeId{35};
    const ls::Scales dual_scales{residual, unknown, dual_speed, dual_gamma, dual_context};
    auto transpose_request = circulation_request(dual_scales, 2, ls::Norm::L2);
    transpose_request.rhs_generation = core::RhsGeneration{36};
    transpose_request.retained_generation = core::NumericalGeneration{37};

    const auto run = [&](const cf::KnownCirculationDerivative& derivative) {
        core::ResourcePlan trial_memory{128ULL * 1024 * 1024};
        return cf::solve_discrete_known_circulation(input, flow, derivative, solve_limits,
                                                    primal_request, shared, *local, *kernel,
                                                    trial_memory);
    };
    const auto dx = coordinates(.25, 0, 0);
    const auto zero = coordinates(0, 0, 0);
    const std::array displacements{dx, dx, dx, dx, dx, dx, dx, dx};
    const std::array d_rows{
        cf::CirculationCollocationDirection{fixture.faces[0].id, dx, {zero, zero}, zero},
        cf::CirculationCollocationDirection{fixture.faces[1].id, std::nullopt, {zero, zero}, zero}};
    const cf::DiscreteCirculationDirection state{fixture.points, displacements, d_rows, zero};
    const auto tangent =
        take(run(cf::KnownCirculationTangent{state, a::SingleWord{.5}, tangent_request}));
    const std::array seed{cf::CirculationLinearRhs{fixture.faces[0].id, a::SingleWord{.25}},
                          cf::CirculationLinearRhs{fixture.faces[1].id, a::SingleWord{-.15625}}};
    const auto dual = take(run(cf::KnownCirculationAdjoint{seed, transpose_request}));
    const auto just_primal = take(run(std::monostate{}));
    reject(just_primal.linearization(), core::FailureCode::InvalidInput);
    const auto& values = take(tangent.solution()).get();
    const auto& dv = take(tangent.linearization()).get();
    const auto& lambda = take(dual.linearization()).get();
    const auto approx = [](const Rat& actual, double expected) {
        CHECK(abs(actual - rational_word(expected)) < rational_word(1e-12));
    };
    const auto center = [](const a::ScaledEnclosure& value) {
        return sum_words(value.center()) * two_power(value.exponent());
    };
    const std::array primal_truth{2.049364665787067014, .313915922823822824};
    const std::array tangent_truth{1.207196834651639385, 1.800156984102001041};
    const std::array dual_truth{-.264772325234835974, .174155889600154652};
    for (std::size_t i = 0; i < primal_truth.size(); ++i) {
        approx(solution_word(values, i), primal_truth[i]);
        approx(solution_word(take(dual.solution()).get(), i), primal_truth[i]);
        approx(solution_word(take(just_primal.solution()).get(), i), primal_truth[i]);
        approx(solution_word(dv.solution(), i), tangent_truth[i]);
        approx(solution_word(lambda.solution(), i), dual_truth[i]);
        CHECK(take(tangent.row(i)).face == fixture.faces[i].id);
        CHECK(take(tangent.forcing(i)).face == fixture.faces[i].id);
        residual_small(take(tangent.row(i)).original_residual, rational_word(1e-13));
        residual_small(take(dv.row(i)).original_residual, rational_word(1e-13));
        residual_small(take(lambda.row(i)).original_residual, rational_word(1e-13));
    }
    exact(take(take(tangent.flow(0)).get().position(0)), Number{3, 4});
    exact(take(take(tangent.flow(1)).get().position(0)), Number{5, 8});
    exact(take(take(tangent.flow(1)).get().position(1)), Number{19, 8});
    const auto free = take(dual.collocation(0)).free_query;
    if (!free)
        FAIL("Expected a covector for the free collocation point");
    CHECK_FALSE(take(dual.collocation(1)).free_query.has_value());
    const auto point = take(dual.point(6));
    CHECK(point.point == a::PointId{7});
    const auto time_dual = center(take(dual.flow_time_dual()));
    approx(time_dual, -.041049359793944367);
    approx(center((*free)[0]), 1.895651483847442373);
    approx(center(point.position[0]), -.826374462193735972);

    Rat action = center((*free)[0]) / 4 + time_dual / 2;
    for (std::size_t i = 0; i < fixture.points.size(); ++i)
        action += center(take(dual.point(i)).position[0]) / 4;
    const Rat transpose_work =
        solution_word(dv.solution(), 0) / 4 - 5 * solution_word(dv.solution(), 1) / 32;
    CHECK(abs(action - transpose_work) < rational_word(1e-12));

    // The tilted direction uses nonzero tangent covectors. Match the exact
    // coordinate difference route against explicit deltas for the same nominal
    // chart, without folding the tangent axes into the point directions.
    const auto scenario = [&](std::span<const cf::DiscreteCollocation> sites,
                              const cf::KnownCirculationDerivative& derivative) {
        auto view = input;
        view.collocations = sites;
        core::ResourcePlan trial_memory{128ULL * 1024 * 1024};
        return cf::solve_discrete_known_circulation(view, flow, derivative, solve_limits,
                                                    primal_request, shared, *local, *kernel,
                                                    trial_memory);
    };
    const auto from_geometry = [&](a::FacetId facet, const a::ResolvedQuery& query) {
        const auto ids = test_geometry::take(geometry.facet_points(facet));
        const auto position = [&](std::size_t i) {
            return test_geometry::take(geometry.point_words(ids[i], a::Configuration::Current));
        };
        const cf::RoundedFacetTangents nominal{{position(0), position(1), position(2)}};
        return cf::DiscreteCollocation{
            nominal, facet, query, coordinates(.25, 0, .5), a::SingleWord{0}, a::SingleWord{0}};
    };
    const std::array exact_sites{from_geometry(a::FacetId{202}, second),
                                 from_geometry(a::FacetId{200}, first)};
    auto explicit_sites = exact_sites;
    for (auto& site : explicit_sites)
        site.tangents = cf::ExplicitCollocationTangents{coordinates(1, 1, 0), coordinates(0, 1, 0)};
    const auto dz = coordinates(0, 0, .125);
    const std::array zero_points{zero, zero, zero, zero, zero, zero, zero, zero};
    const std::array tilt{
        cf::CirculationCollocationDirection{fixture.faces[0].id, zero, {dz, zero}, zero},
        cf::CirculationCollocationDirection{fixture.faces[1].id, std::nullopt, {zero, zero}, zero}};
    const cf::DiscreteCirculationDirection tilt_direction{fixture.points, zero_points, tilt, zero};
    const cf::KnownCirculationTangent tilt_tangent{tilt_direction, a::SingleWord{0},
                                                   tangent_request};
    const cf::KnownCirculationAdjoint tilt_dual{seed, transpose_request};
    const auto exact_forward = take(scenario(exact_sites, tilt_tangent));
    const auto explicit_forward = take(scenario(explicit_sites, tilt_tangent));
    const auto exact_reverse = take(scenario(exact_sites, tilt_dual));
    const auto explicit_reverse = take(scenario(explicit_sites, tilt_dual));
    const auto& exact_tangent = take(exact_forward.linearization()).get();
    const auto& explicit_tangent = take(explicit_forward.linearization()).get();
    const auto midpoint_covector = center(take(exact_reverse.collocation(0)).tangents[0][2]);
    CHECK(abs(midpoint_covector) > rational_word(1e-6));
    CHECK(midpoint_covector == center(take(explicit_reverse.collocation(0)).tangents[0][2]));
    for (std::size_t i = 0; i < exact_sites.size(); ++i) {
        CHECK(solution_word(take(exact_forward.solution()).get(), i) ==
              solution_word(take(explicit_forward.solution()).get(), i));
        CHECK(solution_word(exact_tangent.solution(), i) ==
              solution_word(explicit_tangent.solution(), i));
    }
    CHECK(abs(solution_word(exact_tangent.solution(), 0) / 4 -
              5 * solution_word(exact_tangent.solution(), 1) / 32 - midpoint_covector / 8) <
          rational_word(1e-12));

    auto stale = tangent_request;
    stale.retained_generation = primal_request.retained_generation;
    reject(run(cf::KnownCirculationTangent{state, a::SingleWord{.5}, stale}),
           core::FailureCode::StaleGeneration);
    auto stale_dual = transpose_request;
    stale_dual.rhs_generation = primal_request.rhs_generation;
    reject(run(cf::KnownCirculationAdjoint{seed, stale_dual}), core::FailureCode::StaleGeneration);
    stale_dual = transpose_request;
    stale_dual.retained_generation = primal_request.retained_generation;
    reject(run(cf::KnownCirculationAdjoint{seed, stale_dual}), core::FailureCode::StaleGeneration);
    const Number wrong_time{1};
    flow.time.time = wrong_time.view();
    reject(run(std::monostate{}), core::FailureCode::StaleGeneration);
}

TEST_CASE("Known flow drives discrete midpoint loads and their implicit actions",
          "[AE-07][UW-D05][UW-D06][N4-03][N4-06][verification][negative]") {
    auto source_case = gust_case("deg", R"({"kind":"constant","value":[1]})");
    source_case.field("analysis").field("intervals").elements[0].set("end", parse("2"));
    source_case.field("flow").set("uniform", parse(R"({"kind":"constant","value":[0,0,1]})"));
    auto& gust = source_case.field("flow").field("gust");
    gust.set("direction", parse("[0,0,1]"));
    gust.set("wavevector", parse("[180,0,0]"));
    gust.set("phase", parse("0"));
    gust.set("origin", parse("[0,0,0]"));
    gust.set("time_origin", parse("1.25"));
    const DisjointCirculationPanels fixture;
    auto geometry_input = fixture.input();
    geometry_input.points = geometry_input.points.first(4);
    geometry_input.bindings = geometry_input.bindings.first(4);
    const a::SinglePositions positions{std::span{fixture.xyz}.first(12)};
    geometry_input.reference = positions;
    geometry_input.current = positions;
    geometry_input.edges = geometry_input.edges.first(5);
    geometry_input.faces = geometry_input.faces.first(1);
    geometry_input.bands = geometry_input.bands.first(1);
    geometry_input.quads = geometry_input.quads.first(1);
    const a::GlobalLimits cap{std::numeric_limits<std::uint64_t>::max() / 4,
                              numerics::linear::EstimateAdmission{4096}};
    core::ResourcePlan memory{128ULL * 1024 * 1024};
    auto geometry_work = test_geometry::workspace(memory);
    const auto geometry = test_geometry::take(
        a::prepare_owned_geometry(geometry_input, test_geometry::limits(), *geometry_work, memory));
    auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
    const auto topology = take_global(a::build_global_topology(
        geometry, {std::span{fixture.roles}.first(1), std::span{fixture.owners}.first(5), {}}, cap,
        *global, memory));
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, cap, *global, memory));
    const auto query = test_geometry::take(
        a::resolve_query(geometry, {a::QueryId{1}, a::FreePoint{a::FreeSingle{{.75, .5, 0}}}},
                         geometry.context(), test_geometry::limits(), *geometry_work, memory));
    const auto zero = coordinates(0, 0, 0);
    const auto dx = coordinates(.25, 0, 0);
    const auto id = fixture.faces[0].id;
    const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
    const std::array rows{cf::DiscreteCollocation{
        a::FacetId{200}, query, frame, coordinates(0, 0, .25), a::SingleWord{8}, a::SingleWord{0}}};
    const cf::DiscreteNonpenetrationInput input{contact,
                                                rows,
                                                {},
                                                {core::AcceptedStateId{2}, {}},
                                                {7, 11},
                                                stamp,
                                                p::RequestedSide::Automatic};
    a::UniformFlowWorkspace flow_workspace;
    const auto provider = prepare_gust(source_case, io::SourceId{461}, flow_workspace, memory);
    const auto current = prepare_gust(source_case, io::SourceId{462}, flow_workspace, memory);
    const Number time{5, 4}, absolute{1, 1000000000};
    const cf::DiscreteKnownFlowInput flow{provider,
                                          current,
                                          {stamp, time.view()},
                                          gust_projection(core::Precision::TwoWord, absolute),
                                          gust_limits,
                                          test_geometry::limits(),
                                          *geometry_work,
                                          flow_workspace};
    const cf::DiscreteNonpenetrationLimits solve_limits{
        {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10},
        {test_geometry::limits(),
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(24000000000)};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const std::array factors{1.};
    const std::array<ls::UnitDimensions, 1> gamma{{{2, 0, -1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 1> velocity{{{1, 0, -1, 0, 0, 0, 0}}};
    const ls::Scales scales{factors, factors, gamma, velocity,
                            circulation_context(geometry.context())};
    const auto primal_request = circulation_request(scales, 1, ls::Norm::Infinity);
    auto tangent_request = primal_request;
    tangent_request.rhs_generation = core::RhsGeneration{30};
    tangent_request.retained_generation = core::NumericalGeneration{31};
    const std::array<ls::UnitDimensions, 1> dual_gamma{{{-2, 0, 1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 1> dual_speed{{{-1, 0, 1, 0, 0, 0, 0}}};
    auto dual_context = scales.context;
    dual_context.cache.matrix = core::NumericalGeneration{32};
    dual_context.cache.scales = core::ScaleGeneration{33};
    dual_context.unknown_layout = core::ScopeId{34};
    dual_context.residual_layout = core::ScopeId{35};
    const ls::Scales dual_scales{factors, factors, dual_speed, dual_gamma, dual_context};
    auto transpose_request = circulation_request(dual_scales, 1, ls::Norm::Infinity);
    transpose_request.rhs_generation = core::RhsGeneration{36};
    transpose_request.retained_generation = core::NumericalGeneration{37};

    const std::array edges{cf::EdgeMidpointVelocity{a::EdgeId{10}, zero},
                           cf::EdgeMidpointVelocity{a::EdgeId{11}, zero},
                           cf::EdgeMidpointVelocity{a::EdgeId{12}, zero},
                           cf::EdgeMidpointVelocity{a::EdgeId{13}, zero}};
    const std::array rate{cf::BoundCirculationRate{id, a::SingleWord{.375}}};
    const cf::DiscreteVectorLoadInput loads{
        input.trial_stamp,   primal_request.retained_generation, edges, rate,
        a::SingleWord{1.25}, coordinates(.25, -.125, .5)};
    const auto run = [&](const cf::KnownLoadDerivative& derivative) {
        core::ResourcePlan trial_memory{128ULL * 1024 * 1024};
        return cf::solve_discrete_known_loads(input, loads, flow, derivative, solve_limits,
                                              primal_request, shared, *local, *kernel,
                                              trial_memory);
    };
    const std::array displacements{dx, dx, dx, dx};
    const std::array d_rows{cf::CirculationCollocationDirection{id, dx, {zero, zero}, zero}};
    const cf::DiscreteCirculationDirection state{geometry_input.points, displacements, d_rows,
                                                 zero};
    const std::array drate{cf::BoundCirculationRate{id, a::SingleWord{.125}}};
    const cf::DiscreteImplicitVectorLoadDirection dload{
        input.trial_stamp, primal_request.retained_generation, state, edges, drate};
    const auto tangent =
        take(run(cf::KnownLoadTangent{dload, a::SingleWord{.25}, tangent_request}));
    const auto plain = take(run(std::monostate{}));
    reject(plain.force_direction(), core::FailureCode::InvalidInput);
    const cf::DiscreteLoadSampleWeights none{zero, zero};
    const std::array edge_seeds{cf::DiscreteLoadSampleWeights{zero, coordinates(1, 0, 0)}, none,
                                none, none};
    const std::array facet_seeds{none, none};
    const std::array point_seeds{zero, zero, zero, zero};
    const cf::DiscreteVectorLoadWeights weights{input.trial_stamp,
                                                primal_request.retained_generation,
                                                edge_seeds,
                                                facet_seeds,
                                                point_seeds,
                                                coordinates(1, 2, 3),
                                                coordinates(.125, -.25, .5)};
    const auto dual = take(run(cf::KnownLoadAdjoint{weights, transpose_request}));

    const auto center = [](const a::ScaledEnclosure& value) {
        return sum_words(value.center()) * two_power(value.exponent());
    };
    const auto compare = [&](const Rat& actual, double expected) {
        CHECK(abs(actual - rational_word(expected)) < rational_word(1e-12));
    };
    const auto vector = [&](const cf::DiscreteLoadVector& actual,
                            const std::array<double, 3>& expected) {
        for (std::size_t i = 0; i < 3; ++i)
            compare(center(actual[i]), expected[i]);
    };
    const auto& circ = take(tangent.circulation()).get();
    compare(solution_word(take(circ.solution()).get(), 0), 1.557317918884005427);
    compare(solution_word(take(circ.linearization()).get().solution(), 0), -.296777375053794743);
    residual_small(take(circ.row(0)).original_residual, rational_word(1e-13));

    const std::array edge_x{Number{0}, Number{1, 2}, Number{1}, Number{1, 2}};
    const std::array<Rat, 4> edge_speed{Rat{1}, Rat{2}, Rat{1}, Rat{2}};
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& jet = take(tangent.midpoint_flow(i)).get();
        exact(take(jet.position(0)), edge_x[i]);
        exact(take(jet.time()), time);
        within(take(jet.velocity(2)), edge_speed[i], Rat{0});
    }
    vector(take(plain.total_force()), {0, 0, .46875});
    vector(take(tangent.total_moment()), {.29296875, -.1171875, 0});
    vector(take(tangent.edge(0)).force, {-.928211413046143142, 0, 0});
    vector(take(tangent.edge(1)).force, {0, 2.874858811651149926, 0});
    vector(take(tangent.edge_direction(0)).force, {-.781640702008702757, 0, 0});
    vector(take(tangent.force_direction()), {-1.528893291646792808, 0, .15625});
    vector(take(tangent.moment_direction()), {.09765625, .608196645823396404, .955558307279245505});
    const auto time_dual = center(take(dual.flow_time_dual()));
    const cf::LoadCirculationRateCovector rate_covector = take(dual.circulation_dual(0));
    CHECK(rate_covector.face == take(circ.row(0)).face);
    residual_expected(rate_covector.rate, Rat{1005, 256});
    reject(dual.circulation_dual(1), core::FailureCode::InvalidInput);
    const auto rate_dual = center(rate_covector.rate);
    compare(time_dual, 8.095466005886887656);
    const auto free = take(dual.collocation_dual(0)).free_query;
    if (!free)
        FAIL("Expected free collocation coordinate dual");
    Rat pullback = (time_dual + center((*free)[0])) / 4 + rate_dual / 8;
    for (std::size_t i = 0; i < geometry_input.points.size(); ++i)
        pullback += center(take(dual.point_dual(i)).position[0]) / 4;
    compare(pullback, -1.503846970221721914);
    const auto df = take(tangent.force_direction());
    const auto dm = take(tangent.moment_direction());
    const Rat forward = center(df[0]) + 2 * center(df[1]) + 3 * center(df[2]) + center(dm[0]) / 8 -
                        center(dm[1]) / 4 + center(dm[2]) / 2 +
                        center(take(tangent.edge_direction(0)).force[0]);
    CHECK(abs(pullback - forward) < rational_word(1e-12));
    auto bad = transpose_request;
    bad.rhs_generation = primal_request.rhs_generation;
    reject(run(cf::KnownLoadAdjoint{weights, bad}), core::FailureCode::StaleGeneration);
}

TEST_CASE("Free wake motion uses all sources at attached and inherited material points",
          "[AE-07][UW-D01][UW-D02][UW-D03][UW-D04][UW-D07][N4-06][J1-D01][FD-07][verification]["
          "negative]") {
    auto source_case = gust_case("deg", R"({"kind":"constant","value":[1]})");
    const auto schedule = parse(R"([{"start":1,"end":1.5,"step":0.125,"restart":"initial"}])");
    source_case.field("analysis").set("intervals", schedule);
    source_case.field("analysis").field("wake").set("convection", parse(R"({"kind":"free"})"));
    source_case.field("flow").set("uniform", parse(R"({"kind":"constant","value":[0,0,1]})"));
    auto& gust = source_case.field("flow").field("gust");
    gust.set("direction", parse("[0,0,1]"));
    gust.set("phase", parse("0"));
    gust.set("time_origin", parse("1.25"));

    ResidualFixture fixture;
    const a::PointId alias{7}, coincident{8};
    const auto& id = fixture.id;
    const std::array ids{id[0], id[1], id[2], id[3], id[4], id[5], alias, coincident};
    auto inherited_alias = test_geometry::binding();
    inherited_alias.ownership = a::AliasOf{id[3]};
    const auto& bind = fixture.bindings;
    const std::array bindings{bind[0], bind[1], bind[2],         bind[3],
                              bind[4], bind[5], inherited_alias, test_geometry::binding()};
    const std::array<double, 24> xyz{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0,
                                     2, 1, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0};
    const a::EdgeId alias_edge{20};
    const auto& e = fixture.edges;
    const std::array edges{e[0],
                           e[1],
                           e[2],
                           e[3],
                           e[4],
                           e[5],
                           e[6],
                           a::EdgeInput{e[7].id, {id[5], alias}},
                           a::EdgeInput{e[8].id, {alias, id[4]}},
                           a::EdgeInput{alias_edge, {alias, id[2]}}};
    const auto& own = fixture.edge_owners;
    const std::array ownership{
        own[0], own[1], own[2], own[3], own[4],
        own[5], own[6], own[7], own[8], a::EdgeBindingInput{alias_edge, e[2].id}};
    fixture.quads[1].points[0] = alias;
    fixture.quads[1].edges[0] = alias_edge;
    const std::array<a::TrailingEdgeInput, 1> trailing{{{alias_edge, {alias, id[2]}}}};
    auto input_geometry = fixture.input();
    input_geometry.points = ids;
    input_geometry.bindings = bindings;
    input_geometry.edges = edges;
    input_geometry.reference = a::SinglePositions{xyz};
    input_geometry.current = a::SinglePositions{xyz};
    const Number transport_time{5, 4};
    input_geometry.exact_time = transport_time.view();
    const a::GlobalLimits cap{std::numeric_limits<std::uint64_t>::max() / 4,
                              numerics::linear::EstimateAdmission{4096}};
    core::ResourcePlan memory{128ULL * 1024 * 1024};
    auto geometry_work = test_geometry::workspace(memory);
    const auto geometry = test_geometry::take(
        a::prepare_owned_geometry(input_geometry, test_geometry::limits(), *geometry_work, memory));
    auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
    const auto topology = take_global(a::build_global_topology(
        geometry, {fixture.roles, ownership, trailing}, cap, *global, memory));
    REQUIRE(topology.trailing_edges().size() == 1);
    CHECK(topology.trailing_edges().front().coefficient == 1);
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, cap, *global, memory));
    const auto q = test_geometry::take(
        a::resolve_query(geometry, {a::QueryId{1}, a::FreePoint{a::FreeSingle{{.75, .5, 0}}}},
                         geometry.context(), test_geometry::limits(), *geometry_work, memory));
    const auto zero = coordinates(0, 0, 0);
    const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
    const std::array bound{cf::DiscreteCollocation{a::FacetId{200}, q, frame, zero,
                                                   a::SingleWord{.5}, a::SingleWord{.125}}};
    const a::P1StrengthId old_id{55};
    const std::array wake{cf::DiscreteWakePanel{a::FacetId{202}, cf::AttachedWakeFront{}, old_id}};
    const std::array old{cf::FrozenWakeTrace{old_id, a::SingleWord{-.25}}};
    const cf::DiscreteNonpenetrationInput input{contact,
                                                bound,
                                                wake,
                                                {core::AcceptedStateId{2}, old},
                                                {7, 11},
                                                stamp,
                                                p::RequestedSide::Automatic};
    a::UniformFlowWorkspace flow_workspace;
    const auto provider = prepare_gust(source_case, io::SourceId{481}, flow_workspace, memory);
    const auto current = prepare_gust(source_case, io::SourceId{482}, flow_workspace, memory);
    const Number time{5, 4}, absolute{1, 1000000000};
    const cf::DiscreteKnownFlowInput flow{provider,
                                          current,
                                          {stamp, time.view()},
                                          gust_projection(core::Precision::TwoWord, absolute),
                                          gust_limits,
                                          test_geometry::limits(),
                                          *geometry_work,
                                          flow_workspace};
    const cf::DiscreteFreeWakeLimits cap_velocity{
        {{core::Precision::TwoWord, 1e-10},
         {core::Precision::TwoWord, 1e-10, 1e-10},
         {test_geometry::limits(),
          {numerics::linear::EstimateAdmission{65536}, 1000000000},
          UINT64_C(5000000000)},
         100000,
         UINT64_C(24000000000)},
        2};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const std::array queries{alias, id[4]};
    const std::array points{id[0], id[1], id[2], id[3], id[4], id[5], coincident};
    const auto dx = coordinates(.125, 0, 0);
    const std::array displacements{dx, dx, dx, dx, dx, dx, dx};
    const std::array dhistory{cf::FrozenWakeTrace{old_id, a::SingleWord{.3125}}};
    const auto dambient = coordinates(.125, -.0625, .125);
    const cf::DiscreteFreeWakeDirection d{points,   displacements,      dhistory,
                                          dambient, a::SingleWord{.25}, a::SingleWord{.03125}};
    const std::array seeds{coordinates(.25, -.125, .5), coordinates(-.125, .25, -.375)};
    const auto run = [&](std::span<const a::PointId> particles,
                         const cf::DiscreteFreeWakeDirection* derivative,
                         std::span<const cf::NonpenetrationVector> weights,
                         const cf::DiscreteFreeWakeLimits& caps) {
        return cf::evaluate_discrete_free_wake_velocity(input, flow, particles, derivative, weights,
                                                        caps, shared, *local, *kernel, memory);
    };
    const auto result = take(run(queries, &d, seeds, cap_velocity));
    CHECK(result.size() == 2);
    const auto origin = take(result.source());
    CHECK(origin.parent == input.history.parent);
    CHECK(origin.stamp == input.trial_stamp);
    CHECK(origin.geometry.model == geometry.context().model);
    CHECK(result.point_count() == points.size());
    CHECK(result.strength_count() == 2);
    CHECK(result.used_steps() <= cap_velocity.action.steps);
    const auto center = [](const a::ScaledEnclosure& value) {
        return sum_words(value.center()) * two_power(value.exponent());
    };
    const auto near = [](const Rat& actual, double expected) {
        CHECK(abs(actual - rational_word(expected)) < rational_word(1e-12));
    };
    const std::array vz{.834765747523246768, 2.093805200736806681};
    const std::array dvz{-.070067978056303579, .065055705248403593};
    Rat forward{};
    for (std::size_t i = 0; i < queries.size(); ++i) {
        const auto row = take(result.row(i));
        CHECK(row.particle == queries[i]);
        CHECK(row.owner == (i == 0 ? id[3] : id[4]));
        near(center(row.velocity[0]), 0);
        near(center(row.velocity[1]), 0);
        near(center(row.velocity[2]), vz[i]);
        if (!row.direction)
            FAIL("Expected a particle velocity direction");
        const auto& tangent = row.direction;
        if (tangent) {
            near(center((*tangent)[0]), .125);
            near(center((*tangent)[1]), -.0625);
            near(center((*tangent)[2]), dvz[i]);
            for (std::size_t axis = 0; axis < 3; ++axis)
                forward += center((*tangent)[axis]) * sum_words(seeds[i][axis]);
        }
        const auto& ambient = take(result.ambient(i)).get();
        exact(take(ambient.time()), time);
        exact(take(ambient.position(0)), Number{static_cast<std::uint32_t>(i + 1)});
        within(take(ambient.velocity(2)), Rat{i + 1}, Rat{0});
    }
    const auto old_dual = take(result.strength_dual(0));
    const auto gamma_dual = take(result.strength_dual(1));
    CHECK(old_dual.id == old_id);
    CHECK(gamma_dual.id == a::P1StrengthId::bound(fixture.faces[0].id));
    near(center(old_dual.value), .105390839944920723);
    near(center(gamma_dual.value), -.182892733056897881);
    const auto fields = take(result.field_dual());
    near(center(fields.flow_time), -.785398163397448310);
    near(center(fields.core_radius), .417499074486477147);
    Rat reverse = center(old_dual.value) * 5 / 16 + center(gamma_dual.value) / 8 +
                  center(fields.flow_time) / 4 + center(fields.core_radius) / 32;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        near(center(fields.uniform_flow[axis]), .125);
        reverse += center(fields.uniform_flow[axis]) * sum_words(dambient[axis]);
    }
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto point = take(result.point_dual(i));
        CHECK(point.point == points[i]);
        reverse += center(point.position[0]) / 8;
        if (point.point == id[3])
            near(center(point.position[0]), .754642014437188603);
        if (point.point == coincident)
            for (const auto& component : point.position)
                near(center(component), 0);
    }
    CHECK(abs(forward - reverse) < rational_word(1e-12));

    const std::array canonical{id[3]};
    const auto nominal = take(run(canonical, nullptr, {}, cap_velocity));
    near(center(take(nominal.row(0)).velocity[2]), vz[0]);
    CHECK_FALSE(take(nominal.row(0)).direction.has_value());
    auto generous_coordinates = flow;
    generous_coordinates.coordinate_limits.max_steps = cap_velocity.action.steps;
    const auto with_generous_cap = take(cf::evaluate_discrete_free_wake_velocity(
        input, generous_coordinates, canonical, nullptr, {}, cap_velocity, shared, *local, *kernel,
        memory));
    const auto nominal_velocity = take(nominal.row(0)).velocity;
    const auto generous_velocity = take(with_generous_cap.row(0)).velocity;
    for (std::size_t axis = 0; axis < nominal_velocity.size(); ++axis)
        CHECK(center(generous_velocity[axis]) == center(nominal_velocity[axis]));
    reject(nominal.field_dual(), core::FailureCode::InvalidInput);
    const std::array duplicate{alias, id[3]};
    reject(run(duplicate, nullptr, {}, cap_velocity), core::FailureCode::InvalidInput);
    const std::array other{coincident};
    reject(run(other, nullptr, {}, cap_velocity), core::FailureCode::InvalidInput);
    auto reduced = cap_velocity;
    reduced.max_particles = 1;
    const auto before = memory.allocated_bytes();
    reject(run(queries, &d, seeds, reduced), core::FailureCode::ResourceExceeded);
    CHECK(memory.allocated_bytes() == before);
    reduced = cap_velocity;
    reduced.action.steps = 0;
    reject(run(canonical, nullptr, {}, reduced), core::FailureCode::ResourceExceeded);
    core::ResourcePlan no_payload{0};
    reject(cf::evaluate_discrete_free_wake_velocity(input, flow, canonical, nullptr, {},
                                                    cap_velocity, shared, *local, *kernel,
                                                    no_payload),
           core::FailureCode::ResourceExceeded);
    const cf::CaseFlowClockLimits transport_clock_limits{limits, 1, 4};
    const auto wake_si = materialize(source_case, io::SourceId{481});
    const auto wake_grid = take(
        cf::prepare_case_reference_clock(wake_si, stamp, transport_clock_limits, 1000000, memory));
    const cf::CaseFlowClockQuery destination{stamp, 0, 2, p::RequestedSide::Automatic};
    const cf::WakeTransportClock transport_clock{wake_grid, wake_grid, destination};
    const cf::WakeTransportLimits transport_limits{transport_clock_limits, 1, 10000, 4096,
                                                   core::Precision::TwoWord};
    const auto history_contact = [&](const Number& stamp_time, double seconds, double downstream,
                                     std::optional<core::GeometryGeneration> generation = {}) {
        auto shape = xyz;
        shape[12] = downstream;
        auto history_geometry = input_geometry;
        history_geometry.current = a::SinglePositions{shape};
        history_geometry.context.time = seconds;
        if (generation)
            history_geometry.context.current = *generation;
        history_geometry.exact_time = stamp_time.view();
        const auto snapshot = test_geometry::take(a::prepare_owned_geometry(
            history_geometry, test_geometry::limits(), *geometry_work, memory));
        return take_global(a::check_global_geometry(topology, snapshot, cap, *global, memory));
    };
    const auto previous = history_contact(Number{9, 8}, 1.125, 1.875);
    // The position and direction have nonzero second differences on the same time grid.
    const auto older = history_contact(Number{1}, 1., 1.625);
    const auto release = take_global(topology.release_key(id[2], a::ReleaseEpoch{8}));
    const std::array material{cf::WakeMaterialPoint{release, id[4], a::FacetId{202}}};
    const std::array stages{cf::WakeMaterialSnapshot{contact, material},
                            cf::WakeMaterialSnapshot{previous, material},
                            cf::WakeMaterialSnapshot{older, material}};
    const std::array<cf::NonpenetrationVector, 1> dprevious{coordinates(.25, 0, 0)},
        dolder{coordinates(.5, 0, 0)};
    const cf::WakeTransportHistoryDirection history_direction{dprevious, dolder};
    const std::array<cf::NonpenetrationVector, 1> transport_seeds{seeds[1]};
    const std::array selected{id[4]};
    const auto run_material = [&](const auto& history) {
        return cf::evaluate_fixed_free_wake_transport(
            transport_clock, history, input, flow, selected, &d, history_direction, transport_seeds,
            transport_limits, cap_velocity, shared, *local, *kernel, memory);
    };
    auto reversed_time = stages;
    const std::array wrong_time{stages[0], stages[2], stages[1]};
    reject(run_material(wrong_time), core::FailureCode::StaleGeneration);
    auto wrong_release = material;
    wrong_release[0].release.epoch = a::ReleaseEpoch{9};
    reversed_time[2].points = wrong_release;
    reject(run_material(reversed_time), core::FailureCode::StaleGeneration);
    auto attached = material;
    attached[0].point = id[2];
    const std::array attached_history{cf::WakeMaterialSnapshot{contact, attached}, stages[1],
                                      stages[2]};
    reject(run_material(attached_history), core::FailureCode::InvalidInput);

    const auto advected = take(run_material(stages));
    const auto transported = take(advected.row(0));
    CHECK(transported.release == release);
    residual_expected(transported.residual[0], Rat{1, 2});
    near(center(transported.residual[1]), 0.);
    near(center(transported.residual[2]), -vz[1]);
    if (!transported.direction) {
        FAIL("Expected the free transport direction");
        return;
    }
    const auto& tangent = transported.direction.value();
    residual_expected(tangent[0], Rat{-5, 8});
    near(center(tangent[1]), .0625);
    near(center(tangent[2]), -dvz[1]);
    const auto& actual_velocity = take(advected.free_velocity()).get();
    const auto time_dual = take(advected.position_dual(0));
    residual_expected(time_dual.stencil[0][0], Rat{-3, 2});
    residual_expected(time_dual.stencil[1][0], Rat{2});
    residual_expected(time_dual.stencil[2][0], Rat{-1, 2});
    WorkEnclosure left, right;
    left.add(tangent, transport_seeds[0]);
    right.add(time_dual.stencil[0], dx);
    right.add(time_dual.stencil[1], dprevious[0]);
    right.add(time_dual.stencil[2], dolder[0]);
    for (std::size_t i = 0; i < points.size(); ++i)
        right.add(take(actual_velocity.point_dual(i)).position, displacements[i]);
    right.add(take(actual_velocity.strength_dual(0)).value, dhistory[0].strength);
    right.add(take(actual_velocity.strength_dual(1)).value, bound[0].circulation_direction);
    const auto field_dual = take(actual_velocity.field_dual());
    right.add(field_dual.uniform_flow, dambient);
    right.add(field_dual.flow_time, d.flow_time);
    right.add(field_dual.core_radius, d.core_radius);
    near(left.value, .118145889468151347);
    CHECK(abs(left.value - right.value) <= left.radius + right.radius + rational_word(1e-15));

    const Number previous_time{9, 8}, start_time{1};
    const auto& previous_geometry = take_global(previous.geometry()).get();
    const auto& start_geometry = take_global(older.geometry()).get();
    const auto earlier_collocation = [&](const auto& shape) {
        return test_geometry::take(
            a::resolve_query(shape, {a::QueryId{1}, a::FreePoint{a::FreeSingle{{.75, .5, 0}}}},
                             shape.context(), test_geometry::limits(), *geometry_work, memory));
    };
    const auto earlier_query = earlier_collocation(previous_geometry);
    const auto start_query = earlier_collocation(start_geometry);
    const std::array previous_bound{cf::DiscreteCollocation{
        a::FacetId{200}, earlier_query, frame, zero, a::SingleWord{.5}, a::SingleWord{.125}}};
    const std::array start_bound{cf::DiscreteCollocation{a::FacetId{200}, start_query, frame, zero,
                                                         a::SingleWord{.5}, a::SingleWord{.125}}};
    const cf::DiscreteNonpenetrationInput previous_source{
        previous, previous_bound, wake, input.history, {7, 10}, stamp, p::RequestedSide::Automatic};
    const cf::DiscreteNonpenetrationInput start_source{
        older, start_bound, wake, input.history, {7, 9}, stamp, p::RequestedSide::Right};
    const cf::DiscreteKnownFlowInput previous_flow{
        provider,        current,          {stamp, previous_time.view()},
        flow.projection, flow.flow_limits, flow.coordinate_limits,
        *geometry_work,  flow_workspace};
    const cf::DiscreteKnownFlowInput start_flow{
        provider,        current,          {stamp, start_time.view(), p::RequestedSide::Right},
        flow.projection, flow.flow_limits, flow.coordinate_limits,
        *geometry_work,  flow_workspace};
    const std::array trailing_point{id[2]};
    const cf::FirstFreeWakeEndpoint newest{input, flow, selected, &d};
    const cf::FirstFreeWakeEndpoint at_trailing_edge{previous_source, previous_flow, trailing_point,
                                                     &d};
    const cf::FirstFreeWakeEndpoint retained_previous{previous_source, previous_flow, selected, &d};
    const cf::FirstFreeWakeEndpoint retained_start{start_source, start_flow, selected, &d};
    const std::array born{cf::FirstFreeWakePoint{material[0], std::nullopt}};
    const std::array restarted{cf::FirstFreeWakePoint{material[0], material[0]}};
    const auto first = [&](const auto& clock, const auto& trial, const auto& parent, auto epoch,
                           const auto& material_points) {
        return cf::evaluate_first_free_wake_transport(
            clock, trial, parent, epoch, material_points, transport_seeds, transport_limits,
            cap_velocity, shared, *local, *kernel, memory);
    };
    reject(first(transport_clock, newest, at_trailing_edge, a::ReleaseEpoch{9}, born),
           core::FailureCode::StaleGeneration);
    reject(first(transport_clock, newest, retained_previous, release.epoch, born),
           core::FailureCode::StaleGeneration);
    reject(first(transport_clock, newest, retained_previous, std::nullopt, restarted),
           core::FailureCode::InvalidInput);
    reject(first(transport_clock, newest, newest, release.epoch, born),
           core::FailureCode::StaleGeneration);
    const auto birth = take(first(transport_clock, newest, at_trailing_edge, release.epoch, born));
    const auto birth_row = take(birth.row(0));
    REQUIRE(birth_row.direction);
    CHECK(birth_row.release == release);
    residual_expected(birth_row.residual[0], Rat{8});
    residual_expected(birth_row.residual[1], Rat{0});
    residual_expected((*birth_row.direction)[0], Rat{-1, 8});
    const auto& parent_velocity = take(birth.parent_free_velocity()).get();
    const auto& trial_velocity = take(birth.free_velocity()).get();
    const auto parent_row = take(parent_velocity.row(0));
    const auto trial_row = take(trial_velocity.row(0));
    CHECK(parent_row.owner == release.point);
    near(center(parent_row.velocity[2]), 1.024661130622405660);
    REQUIRE(parent_row.direction);
    near(center((*parent_row.direction)[2]), -.064175653096249599);
    near(center(birth_row.residual[2]), -1.559233165679606171);
    near(center((*birth_row.direction)[2]), -.000440026076076997);
    CHECK(trial_row.owner == material[0].point);
    exact(take(take(parent_velocity.ambient(0)).get().time()), previous_time);
    CHECK(abs(center(birth_row.residual[2]) +
              (center(parent_row.velocity[2]) + center(trial_row.velocity[2])) / 2) <
          rational_word(1e-15));
    const auto birth_dual = take(birth.position_dual(0));
    residual_expected(birth_dual.stencil[0][0], Rat{-1});
    residual_expected(birth_dual.stencil[1][0], Rat{1});
    residual_expected(birth_dual.stencil[2][0], Rat{0});
    WorkEnclosure birth_left, birth_right;
    birth_left.add(*birth_row.direction, transport_seeds[0]);
    birth_right.add(birth_dual.stencil[0], dx);
    birth_right.add(birth_dual.stencil[1], dx);
    for (const auto* endpoint : {&trial_velocity, &parent_velocity}) {
        REQUIRE(endpoint->point_count() == points.size());
        for (std::size_t i = 0; i < points.size(); ++i)
            birth_right.add(take(endpoint->point_dual(i)).position, displacements[i]);
        birth_right.add(take(endpoint->strength_dual(0)).value, dhistory[0].strength);
        birth_right.add(take(endpoint->strength_dual(1)).value, bound[0].circulation_direction);
        const auto dual_field = take(endpoint->field_dual());
        birth_right.add(dual_field.uniform_flow, dambient);
        birth_right.add(dual_field.flow_time, d.flow_time);
        birth_right.add(dual_field.core_radius, d.core_radius);
    }
    CHECK(abs(birth_left.value - birth_right.value) <=
          birth_left.radius + birth_right.radius + rational_word(1e-15));

    const cf::WakeTransportClock first_interval{
        wake_grid, wake_grid, {stamp, 0, 1, p::RequestedSide::Automatic}};
    const auto restarted_result =
        take(first(first_interval, retained_previous, retained_start, std::nullopt, restarted));
    const auto restarted_row = take(restarted_result.row(0));
    residual_expected(restarted_row.residual[0], Rat{2});
    const auto& restart_parent = take(restarted_result.parent_free_velocity()).get();
    CHECK(take(restart_parent.row(0)).owner == material[0].point);
    exact(take(take(restart_parent.ambient(0)).get().time()), start_time);
    reject(first(first_interval, retained_previous, retained_start, release.epoch, restarted),
           core::FailureCode::StaleGeneration);

    const cf::FirstFreeWakeEndpoint trailing_start{start_source, start_flow, trailing_point, &d};
    const auto interval_birth =
        take(first(first_interval, retained_previous, trailing_start, release.epoch, born));
    const auto second_release = take_global(topology.release_key(id[3], release.epoch));
    const std::array growth_points{material[0],
                                   cf::WakeMaterialPoint{second_release, id[5], a::FacetId{202}}};
    const cf::WakeMaterialSnapshot growth_parent{older, {}};
    const cf::WakeMaterialSnapshot growth_trial{previous, growth_points};
    auto growth_limits = transport_limits;
    growth_limits.maximum_particles = 2;
    const auto growth = [&](const auto& parent_sample, const auto& trial_sample, auto epoch) {
        return cf::prepare_wake_material_transition(wake_si, first_interval, parent_sample,
                                                    trial_sample, epoch, growth_limits, memory);
    };
    const auto transition = take(growth(growth_parent, growth_trial, release.epoch));
    REQUIRE(transition.size() == 2);
    const std::array classified{take(transition.row(0)), take(transition.row(1))};
    CHECK(classified[0].current.release == release);
    CHECK(classified[1].current.release == second_release);
    CHECK_FALSE(classified[0].parent);
    CHECK_FALSE(classified[1].parent);
    const std::array departing{id[4], id[5]};
    const std::array origins{id[2], id[3]};
    const cf::FirstFreeWakeEndpoint growing_trial{previous_source, previous_flow, departing, {}};
    const cf::FirstFreeWakeEndpoint growing_parent{start_source, start_flow, origins, {}};
    const auto growth_residual = take(cf::evaluate_first_free_wake_transport(
        first_interval, growing_trial, growing_parent, release.epoch, classified, {}, growth_limits,
        cap_velocity, shared, *local, *kernel, memory));
    residual_expected(take(growth_residual.row(0)).residual[0], Rat{7});
    residual_expected(take(growth_residual.row(1)).residual[0], Rat{8});
    const cf::WakeMaterialSnapshot missing_release{previous, material};
    reject(growth(growth_parent, missing_release, release.epoch),
           core::FailureCode::StaleGeneration);
    auto aliased = growth_points;
    aliased[1].point = id[4];
    const cf::WakeMaterialSnapshot repeated_owner{previous, aliased};
    reject(growth(growth_parent, repeated_owner, release.epoch),
           core::FailureCode::StaleGeneration);
    auto crossed = growth_points;
    std::swap(crossed[0].point, crossed[1].point);
    const cf::WakeMaterialSnapshot crossed_releases{previous, crossed};
    reject(growth(growth_parent, crossed_releases, release.epoch),
           core::FailureCode::StaleGeneration);
    const std::array crossed_birth{cf::FirstFreeWakePoint{crossed[0], std::nullopt}};
    const std::array swapped_departure{id[5]};
    const cf::FirstFreeWakeEndpoint wrong_span{
        previous_source, previous_flow, swapped_departure, {}};
    reject(first(first_interval, wrong_span, trailing_start, release.epoch, crossed_birth),
           core::FailureCode::StaleGeneration);
    const cf::WakeMaterialSnapshot already_released{older, material};
    reject(growth(already_released, growth_trial, release.epoch),
           core::FailureCode::StaleGeneration);
    reject(growth(growth_parent, growth_trial, a::ReleaseEpoch{9}),
           core::FailureCode::StaleGeneration);
    reject(cf::prepare_wake_material_transition(wake_si, first_interval, growth_parent,
                                                growth_trial, release.epoch, transport_limits,
                                                memory),
           core::FailureCode::ResourceExceeded);

    const a::PointId far_up{9}, far_down{10};
    const a::EdgeId far_upper{21}, far_end{22}, far_lower{23}, far_diagonal{24};
    const std::array extended_ids{id[0], id[1], id[2],      id[3],  id[4],
                                  id[5], alias, coincident, far_up, far_down};
    const std::array extended_bindings{bindings[0],
                                       bindings[1],
                                       bindings[2],
                                       bindings[3],
                                       bindings[4],
                                       bindings[5],
                                       bindings[6],
                                       bindings[7],
                                       test_geometry::binding(),
                                       test_geometry::binding()};
    const std::array<double, 30> extended_xyz{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 2, 1, 0,
                                              2, 0, 0, 1, 0, 0, 1, 0, 0, 3, 1, 0, 3, 0, 0};
    const std::array extended_edges{edges[0],
                                    edges[1],
                                    edges[2],
                                    edges[3],
                                    edges[4],
                                    edges[5],
                                    edges[6],
                                    edges[7],
                                    edges[8],
                                    edges[9],
                                    a::EdgeInput{far_upper, {id[4], far_up}},
                                    a::EdgeInput{far_end, {far_up, far_down}},
                                    a::EdgeInput{far_lower, {far_down, id[5]}},
                                    a::EdgeInput{far_diagonal, {id[5], far_up}}};
    const a::PhysicalFaceId far_face{102};
    const a::BandId far_band{52};
    const a::FacetId far_facet{204};
    const std::array extended_bands{fixture.bands[0], fixture.bands[1], far_band};
    const std::array extended_faces{fixture.faces[0], fixture.faces[1],
                                    a::FaceInput{far_face, far_band}};
    const std::array extended_roles{fixture.roles[0], fixture.roles[1],
                                    a::FaceRoleInput{far_face, a::AeroFaceRole::Wake}};
    const std::array extended_quads{fixture.quads[0], fixture.quads[1],
                                    a::QuadInput{far_face,
                                                 {id[5], id[4], far_up, far_down},
                                                 {edges[6].id, far_upper, far_end, far_lower},
                                                 far_diagonal,
                                                 {far_facet, a::FacetId{205}}}};
    const std::array extended_ownership{ownership[0],
                                        ownership[1],
                                        ownership[2],
                                        ownership[3],
                                        ownership[4],
                                        ownership[5],
                                        ownership[6],
                                        ownership[7],
                                        ownership[8],
                                        ownership[9],
                                        a::EdgeBindingInput{far_upper, a::CanonicalEdge{}},
                                        a::EdgeBindingInput{far_end, a::CanonicalEdge{}},
                                        a::EdgeBindingInput{far_lower, a::CanonicalEdge{}},
                                        a::EdgeBindingInput{far_diagonal, a::CanonicalEdge{}}};
    auto extended_input = input_geometry;
    extended_input.points = extended_ids;
    extended_input.bindings = extended_bindings;
    extended_input.reference = a::SinglePositions{extended_xyz};
    extended_input.current = a::SinglePositions{extended_xyz};
    extended_input.edges = extended_edges;
    extended_input.faces = extended_faces;
    extended_input.bands = extended_bands;
    extended_input.quads = extended_quads;
    const auto extended_geometry = test_geometry::take(
        a::prepare_owned_geometry(extended_input, test_geometry::limits(), *geometry_work, memory));
    const auto extended_topology = take_global(a::build_global_topology(
        extended_geometry, {extended_roles, extended_ownership, trailing}, cap, *global, memory));
    const auto extended_contact = take_global(
        a::check_global_geometry(extended_topology, extended_geometry, cap, *global, memory));
    const a::ReleaseEpoch next_epoch{9};
    const std::array extended_points{
        cf::WakeMaterialPoint{release, far_up, far_facet},
        cf::WakeMaterialPoint{{release.point, next_epoch}, id[4], a::FacetId{202}},
        cf::WakeMaterialPoint{second_release, far_down, far_facet},
        cf::WakeMaterialPoint{{second_release.point, next_epoch}, id[5], a::FacetId{202}}};
    const cf::WakeMaterialSnapshot extended_trial{extended_contact, extended_points};
    auto extended_limits = growth_limits;
    extended_limits.maximum_particles = 4;
    const auto extended = take(
        cf::prepare_wake_material_transition(wake_si, transport_clock, growth_trial, extended_trial,
                                             next_epoch, extended_limits, memory));
    REQUIRE(extended.size() == 4);
    CHECK_FALSE(extended.window());
    CHECK_FALSE(extended.retirement_source());
    const auto inherited = take(extended.row(0));
    if (!inherited.parent) {
        FAIL("A retained row must keep its parent material point");
        return;
    }
    CHECK(inherited.parent->point == id[4]);
    CHECK(inherited.current.point == far_up);
    CHECK_FALSE(take(extended.row(1)).parent);
    CHECK(take(extended.row(2)).current.point == far_down);

    const auto low_word = std::ldexp(1., -90);
    auto parent_gamma = previous_bound;
    parent_gamma[0].circulation = a::DoubleWord{.5, low_word};
    const std::array parent_trace{cf::FrozenWakeTrace{old_id, a::DoubleWord{-.25, -low_word}}};
    const cf::DiscreteNonpenetrationInput parent_strength{previous,
                                                          parent_gamma,
                                                          wake,
                                                          {input.history.parent, parent_trace},
                                                          previous_source.trial_stamp,
                                                          stamp,
                                                          p::RequestedSide::Automatic};
    const auto extended_query = earlier_collocation(extended_geometry);
    const std::array trial_gamma{cf::DiscreteCollocation{
        a::FacetId{200}, extended_query, frame, zero, a::SingleWord{.75}, a::SingleWord{0}}};
    const a::P1StrengthId newest_id{56}, mapped_old_id{57};
    const std::array extended_panels{
        cf::DiscreteWakePanel{a::FacetId{202}, cf::AttachedWakeFront{}, newest_id},
        cf::DiscreteWakePanel{far_facet, newest_id, mapped_old_id}};
    const std::array extended_trace{
        cf::FrozenWakeTrace{newest_id, a::DoubleWord{.5, low_word}},
        cf::FrozenWakeTrace{mapped_old_id, a::DoubleWord{-.25, -low_word}}};
    const auto strengths = [&](const auto& columns, const auto& panels) {
        const cf::DiscreteNonpenetrationInput next_source{extended_contact,
                                                          trial_gamma,
                                                          panels,
                                                          {input.history.parent, columns},
                                                          input.trial_stamp,
                                                          stamp,
                                                          p::RequestedSide::Automatic};
        return cf::prepare_wake_strength_transition(
            wake_si, transport_clock, {parent_strength, growth_points},
            {next_source, extended_points}, next_epoch, extended_limits, memory);
    };
    auto mutated_trace = extended_trace;
    const auto qualified_strengths = take(strengths(mutated_trace, extended_panels));
    REQUIRE(qualified_strengths.history().size() == 2);
    CHECK(qualified_strengths.material().size() == extended.size());
    CHECK(qualified_strengths.history()[1].id == mapped_old_id);
    REQUIRE(qualified_strengths.history_material().size() == 2);
    CHECK(qualified_strengths.history_material()[1].from == second_release.point);
    CHECK(qualified_strengths.history_material()[1].to == release.point);
    CHECK(qualified_strengths.history_material()[1].epoch == release.epoch);
    CHECK(qualified_strengths.history_material()[0].epoch == next_epoch);
    CHECK(qualified_strengths.retired_history().empty());

    const auto held = std::get<a::DoubleWord>(qualified_strengths.history()[1].strength);
    CHECK(held.high == -.25);
    CHECK(held.low == -low_word);
    mutated_trace[0].strength = a::SingleWord{.75};
    reject(strengths(mutated_trace, extended_panels), core::FailureCode::StaleGeneration);
    mutated_trace = extended_trace;
    mutated_trace[0].strength = a::SingleWord{.5};
    reject(strengths(mutated_trace, extended_panels), core::FailureCode::StaleGeneration);
    mutated_trace = extended_trace;
    mutated_trace[1].strength = a::DoubleWord{-.25, low_word};
    reject(strengths(mutated_trace, extended_panels), core::FailureCode::StaleGeneration);
    CHECK(std::get<a::DoubleWord>(qualified_strengths.history()[1].strength).low == -low_word);
    const cf::DiscreteNonpenetrationInput with_owned_history{
        extended_contact,
        trial_gamma,
        extended_panels,
        {input.history.parent, qualified_strengths.history()},
        input.trial_stamp,
        stamp,
        p::RequestedSide::Automatic};
    const std::array historical_probe{far_up};
    auto extended_velocity_limits = cap_velocity;
    extended_velocity_limits.action.induction.steps *= 3;
    extended_velocity_limits.action.steps *= 3;
    extended_velocity_limits.action.action_steps *= 2;
    core::ResourcePlan extended_velocity_memory{512ULL * 1024 * 1024};
    auto evaluated_velocity = cf::evaluate_discrete_free_wake_velocity(
        with_owned_history, flow, historical_probe, nullptr, {}, extended_velocity_limits, shared,
        *local, *kernel, extended_velocity_memory);
    if (!evaluated_velocity && evaluated_velocity.error().details.resources) {
        const auto& observed = *evaluated_velocity.error().details.resources;
        INFO("Resource demand " << observed.required << ", usage " << observed.observed);
        REQUIRE(evaluated_velocity);
    }
    const auto growing_velocity = take(std::move(evaluated_velocity));
    CHECK(take(growing_velocity.row(0)).owner == far_up);

    const std::array birth_origins{cf::WakeMaterialPoint{release, release.point, {}},
                                   cf::WakeMaterialPoint{second_release, second_release.point, {}}};
    const cf::WakeMaterialSnapshot preceding{older, birth_origins, &growth_residual};
    const std::array extended_direction_points{id[0], id[1],      id[2],  id[3],   id[4],
                                               id[5], coincident, far_up, far_down};
    const std::array extended_displacements{dx, dx, dx, dx, dx, dx, dx, dx, dx};
    const cf::DiscreteFreeWakeDirection trial_direction{extended_direction_points,
                                                        extended_displacements,
                                                        {},
                                                        dambient,
                                                        d.flow_time,
                                                        d.core_radius};
    const cf::FreeWakeTrialEndpoint parent_trial{
        {parent_strength, growth_points}, previous_flow, &d};
    const cf::FreeWakeTrialEndpoint current_trial{
        {with_owned_history, extended_points}, flow, &trial_direction};
    const std::array retained_previous_direction{dprevious[0], dolder[0]};
    const std::array retained_older_direction{dolder[0], dprevious[0]};
    const cf::FreeWakeTrialLimits trial_limits{extended_limits, extended_limits,
                                               extended_velocity_limits};
    const std::array trial_seeds{coordinates(1, 0, .125), coordinates(2, 0, -.25),
                                 coordinates(3, 0, .5), coordinates(4, 0, -.125)};
    core::ResourcePlan trial_memory{1024ULL * 1024 * 1024};
    const auto step = take(cf::evaluate_free_wake_trial(
        wake_si, transport_clock, parent_trial, current_trial,
        {&preceding, {retained_previous_direction, retained_older_direction}}, next_epoch,
        trial_seeds, trial_limits, shared, *local, *kernel, trial_memory));
    if (!step.trapezoid || !step.bdf2) {
        FAIL("An interior trial must include both the released and retained rows");
        return;
    }
    REQUIRE(step.trapezoid->size() == 2);
    REQUIRE(step.bdf2->size() == 2);
    CHECK(step.strength.history_material()[0].epoch == next_epoch);
    CHECK(take(step.trapezoid->row(0)).release == extended_points[1].release);
    CHECK(take(step.trapezoid->row(1)).release == extended_points[3].release);
    CHECK(take(step.bdf2->row(0)).release == release);
    CHECK(take(step.bdf2->row(1)).release == second_release);
    residual_expected(take(step.trapezoid->row(0)).residual[0], Rat{8});
    residual_expected(take(step.trapezoid->row(1)).residual[0], Rat{8});
    residual_expected(take(step.bdf2->row(0)).residual[0], Rat{10});
    residual_expected(take(step.bdf2->row(1)).residual[0], Rat{8});
    residual_expected(take(step.trapezoid->position_dual(0)).stencil[0][0], Rat{16});
    residual_expected(take(step.trapezoid->position_dual(1)).stencil[0][0], Rat{32});
    residual_expected(take(step.bdf2->position_dual(0)).stencil[0][0], Rat{12});
    residual_expected(take(step.bdf2->position_dual(1)).stencil[0][0], Rat{36});

    const auto detached = take(step.trapezoid->row(0));
    const auto retained_top = take(step.bdf2->row(0));
    const auto retained_bottom = take(step.bdf2->row(1));
    if (!detached.direction || !retained_top.direction || !retained_bottom.direction) {
        FAIL("Expected directions for the released and retained rows");
        return;
    }
    residual_expected((*detached.direction)[0], Rat{-1, 8});
    residual_expected((*retained_top.direction)[0], Rat{-5, 8});
    residual_expected((*retained_bottom.direction)[0], Rat{-45, 8});

    const auto parent_bound = a::P1StrengthId::bound(fixture.faces[0].id);
    const std::array parent_directions{
        dhistory[0], cf::FrozenWakeTrace{parent_bound, parent_gamma[0].circulation_direction}};
    const auto transferred = take(cf::apply_wake_strength_transfer(
        step.strength, parent_directions, extended_limits.maximum_steps, trial_memory));
    REQUIRE(transferred.size() == 2);
    CHECK(transferred.values()[0].id == newest_id);
    CHECK(sum_words(transferred.values()[0].strength) == Rat{1, 8});
    CHECK(transferred.values()[1].id == mapped_old_id);
    CHECK(sum_words(transferred.values()[1].strength) == Rat{5, 16});
    for (const auto* transport : std::array{&*step.trapezoid, &*step.bdf2}) {
        const auto& velocity = take(transport->free_velocity()).get();
        auto covectors =
            std::array{take(velocity.strength_dual(0)), take(velocity.strength_dual(1))};
        const auto reversed = take(
            cf::pullback_wake_strength_transfer(step.strength, covectors, core::Precision::TwoWord,
                                                extended_limits.maximum_steps, trial_memory));
        REQUIRE(reversed.size() == 2);
        WorkEnclosure trial_work, parent_work;
        for (std::size_t i = 0; i < covectors.size(); ++i) {
            const auto seed = take(reversed.row(i));
            CHECK(seed.id == parent_directions[i].id);
            trial_work.add(covectors[i].value, transferred.values()[i].strength);
            parent_work.add(seed.value, parent_directions[i].strength);
        }
        CHECK(trial_work.value != 0);
        CHECK(abs(trial_work.value - parent_work.value) <= trial_work.radius + parent_work.radius);
        std::swap(covectors[0], covectors[1]);
        reject(cf::pullback_wake_strength_transfer(step.strength, covectors,
                                                   core::Precision::TwoWord,
                                                   extended_limits.maximum_steps, trial_memory),
               core::FailureCode::StaleGeneration);
    }

    reject(cf::evaluate_free_wake_trial(wake_si, transport_clock, parent_trial, current_trial, {},
                                        next_epoch, trial_seeds, trial_limits, shared, *local,
                                        *kernel, trial_memory),
           core::FailureCode::InvalidInput);

    auto wrong_link = extended_panels;
    wrong_link[1].leading = mapped_old_id;
    reject(strengths(extended_trace, wrong_link), core::FailureCode::StaleGeneration);

    auto limited_case = source_case;
    limited_case.field("numerics")
        .set("wake_rows", parse(R"({"kind":"limited","maximum_rows":1})"));
    const auto limited_si = materialize(limited_case, io::SourceId{486});
    const auto limited_grid = take(cf::prepare_case_reference_clock(
        limited_si, stamp, transport_clock_limits, 1000000, memory));
    const cf::WakeTransportClock limited_first{limited_grid, limited_grid,
                                               first_interval.destination};
    const cf::WakeTransportClock limited_second{limited_grid, limited_grid, destination};
    CHECK(take(cf::prepare_wake_material_transition(limited_si, limited_first, growth_parent,
                                                    growth_trial, release.epoch, extended_limits,
                                                    memory))
              .size() == 2);
    reject(cf::prepare_wake_material_transition(limited_si, limited_second, growth_trial,
                                                extended_trial, next_epoch, extended_limits,
                                                memory),
           core::FailureCode::ResourceExceeded);

    auto window_case = limited_case;
    window_case.field("analysis")
        .field("wake")
        .set(
            "retention",
            parse(R"({"kind":"closed_window","duration":0.125,"distance":1,"error_budget":0.01})"));
    const auto window_si = materialize(window_case, io::SourceId{487});
    const auto window_grid = take(cf::prepare_case_reference_clock(
        window_si, stamp, transport_clock_limits, 1000000, memory));
    const cf::WakeTransportClock window_first{window_grid, window_grid, first_interval.destination};
    const cf::WakeTransportClock window_second{window_grid, window_grid, destination};
    const auto first_window = take(
        cf::prepare_wake_material_transition(window_si, window_first, growth_parent, growth_trial,
                                             release.epoch, extended_limits, memory));
    CHECK(first_window.retired_count() == 0);
    CHECK(first_window.window().has_value());
    CHECK_FALSE(first_window.retirement_source());
    const std::array window_points{extended_points[1], extended_points[3]};
    const cf::WakeMaterialSnapshot window_trial{contact, window_points};
    const auto window = take(cf::prepare_wake_material_transition(
        window_si, window_second, growth_trial, window_trial, next_epoch, extended_limits, memory));
    CHECK(window.size() == 2);
    REQUIRE(window.retired_count() == 2);
    CHECK(take(window.retired(0)).release == release);
    CHECK(take(window.retired(1)).release == second_release);
    CHECK(take(window.retired(0)).point == id[4]);
    CHECK(take(window.row(0)).current.release.epoch == next_epoch);
    const auto& requirements = window.window();
    const auto& retirement_source = window.retirement_source();
    if (!requirements || !retirement_source) {
        FAIL("The finite window must keep its Case limits and retirement source");
        return;
    }
    CHECK(requirements->rows == 1);
    CHECK(pg_conditioning_test::q(requirements->duration_seconds) ==
          pg_conditioning_test::Q{"1/8"});
    CHECK(pg_conditioning_test::q(requirements->distance_metres) == pg_conditioning_test::Q{1});
    CHECK(pg_conditioning_test::q(requirements->error_budget) == pg_conditioning_test::Q{"1/100"});
    CHECK(retirement_source->geometry == test_geometry::take(previous_geometry.identity()));
    const auto retired_at = retirement_source->time_seconds.view();
    exact({retired_at.sign, retired_at.numerator, retired_at.denominator}, previous_time);
    const auto upper_retired = take(window.retired_sample(0));
    const auto lower_retired = take(window.retired_sample(1));
    CHECK(upper_retired.material.release == release);
    CHECK(sum_words(upper_retired.parent_position_metres[0]) == Rat{15, 8});
    CHECK(sum_words(upper_retired.parent_position_metres[1]) == 1);
    CHECK(sum_words(lower_retired.parent_position_metres[0]) == 2);
    CHECK(sum_words(lower_retired.parent_position_metres[1]) == 0);
    reject(window.retired_sample(2), core::FailureCode::InvalidInput);
    const std::array window_panels{
        cf::DiscreteWakePanel{a::FacetId{202}, cf::AttachedWakeFront{}, newest_id}};
    const std::array window_history{extended_trace.front()};
    const cf::DiscreteNonpenetrationInput next_window_source{contact,
                                                             bound,
                                                             window_panels,
                                                             {input.history.parent, window_history},
                                                             input.trial_stamp,
                                                             stamp,
                                                             p::RequestedSide::Automatic};
    const auto window_strengths = take(cf::prepare_wake_strength_transition(
        window_si, window_second, {parent_strength, growth_points},
        {next_window_source, window_points}, next_epoch, extended_limits, memory));
    CHECK(window_strengths.material().retired_count() == 2);
    REQUIRE(window_strengths.history().size() == 1);
    CHECK(window_strengths.history().front().id == newest_id);
    REQUIRE(window_strengths.retired_history().size() == 1);
    const auto& retired_strength = window_strengths.retired_history().front();
    CHECK(retired_strength.material.epoch == release.epoch);
    CHECK(retired_strength.material.from == second_release.point);
    CHECK(retired_strength.trace.id == old_id);
    CHECK(std::get<a::DoubleWord>(retired_strength.trace.strength).low == -low_word);

    ResidualFixture reversed_fixture;
    reversed_fixture.reverse_wake();
    auto reversed_geometry_input = reversed_fixture.input();
    reversed_geometry_input.context.time = 1.25;
    reversed_geometry_input.exact_time = time.view();
    const auto reverse_now = test_geometry::take(a::prepare_owned_geometry(
        reversed_geometry_input, test_geometry::limits(), *geometry_work, memory));
    const auto reverse_topology = take_global(
        a::build_global_topology(reverse_now, reversed_fixture.topology(), cap, *global, memory));
    REQUIRE(reverse_topology.trailing_edges().size() == 1);
    REQUIRE(reverse_topology.trailing_edges().front().coefficient == -1);
    const auto reverse_trial =
        take_global(a::check_global_geometry(reverse_topology, reverse_now, cap, *global, memory));
    reversed_geometry_input.context.time = 1.125;
    reversed_geometry_input.exact_time = previous_time.view();
    const auto reverse_previous = test_geometry::take(a::prepare_owned_geometry(
        reversed_geometry_input, test_geometry::limits(), *geometry_work, memory));
    const auto reverse_parent = take_global(
        a::check_global_geometry(reverse_topology, reverse_previous, cap, *global, memory));
    const auto reverse_query = earlier_collocation(reverse_now);
    const auto reverse_previous_query = earlier_collocation(reverse_previous);
    const std::array reverse_parent_gamma{
        cf::DiscreteCollocation{a::FacetId{200}, reverse_previous_query, frame, zero,
                                a::DoubleWord{.5, low_word}, a::SingleWord{0}}};
    const std::array reverse_trial_gamma{cf::DiscreteCollocation{
        a::FacetId{200}, reverse_query, frame, zero, a::SingleWord{.75}, a::SingleWord{0}}};
    const cf::DiscreteNonpenetrationInput reverse_parent_source{
        reverse_parent,
        reverse_parent_gamma,
        wake,
        {input.history.parent, parent_trace},
        previous_source.trial_stamp,
        stamp,
        p::RequestedSide::Automatic};
    const auto reverse_strengths = [&](const auto& columns) {
        const cf::DiscreteNonpenetrationInput reverse_trial_source{reverse_trial,
                                                                   reverse_trial_gamma,
                                                                   window_panels,
                                                                   {input.history.parent, columns},
                                                                   input.trial_stamp,
                                                                   stamp,
                                                                   p::RequestedSide::Automatic};
        return cf::prepare_wake_strength_transition(
            window_si, window_second, {reverse_parent_source, growth_points},
            {reverse_trial_source, window_points}, next_epoch, extended_limits, memory);
    };
    const std::array reversed_trace{cf::FrozenWakeTrace{newest_id, a::DoubleWord{-.5, -low_word}}};
    const auto signed_strengths = take(reverse_strengths(reversed_trace));
    CHECK(std::get<a::DoubleWord>(signed_strengths.history().front().strength).low == -low_word);
    const std::array bound_direction{
        cf::FrozenWakeTrace{parent_bound, a::DoubleWord{.5, low_word}}};
    const auto signed_direction = take(cf::apply_wake_strength_transfer(
        signed_strengths, bound_direction, extended_limits.maximum_steps, trial_memory));
    REQUIRE(signed_direction.size() == 1);
    const auto signed_value = std::get<a::DoubleWord>(signed_direction.values()[0].strength);
    CHECK(signed_value.high == -.5);
    CHECK(signed_value.low == -low_word);
    const auto& newest_velocity = take(step.trapezoid->free_velocity()).get();
    const std::array newest_covector{take(newest_velocity.strength_dual(0))};
    const auto signed_dual = take(cf::pullback_wake_strength_transfer(
        signed_strengths, newest_covector, core::Precision::TwoWord, extended_limits.maximum_steps,
        trial_memory));
    REQUIRE(signed_dual.size() == 1);
    const auto parent_covector = take(signed_dual.row(0));
    CHECK(parent_covector.id == parent_bound);
    WorkEnclosure newest_work, signed_parent_work;
    newest_work.add(newest_covector[0].value, signed_direction.values()[0].strength);
    signed_parent_work.add(parent_covector.value, bound_direction[0].strength);
    CHECK(newest_work.value != 0);
    CHECK(abs(newest_work.value - signed_parent_work.value) <=
          newest_work.radius + signed_parent_work.radius);
    reject(reverse_strengths(window_history), core::FailureCode::StaleGeneration);

    const cf::WakeMaterialSnapshot omitted_window{previous, {}};
    reject(cf::prepare_wake_material_transition(window_si, window_second, omitted_window,
                                                window_trial, next_epoch, extended_limits, memory),
           core::FailureCode::StaleGeneration);
    reject(cf::prepare_wake_material_transition(wake_si, transport_clock, omitted_window,
                                                window_trial, next_epoch, extended_limits, memory),
           core::FailureCode::StaleGeneration);

    auto alias_bindings = bindings;
    alias_bindings[7].ownership = a::AliasOf{id[4]};
    auto alias_reference = xyz;
    std::ranges::copy(std::span{xyz}.subspan(12, 3), alias_reference.begin() + 21);
    auto alias_current = alias_reference;
    alias_current[12] = alias_current[21] = 1.875;
    auto alias_input = input_geometry;
    alias_input.bindings = alias_bindings;
    alias_input.reference = a::SinglePositions{alias_reference};
    alias_input.current = a::SinglePositions{alias_current};
    alias_input.context.time = 1.125;
    alias_input.exact_time = previous_time.view();
    const auto alias_geometry = test_geometry::take(
        a::prepare_owned_geometry(alias_input, test_geometry::limits(), *geometry_work, memory));
    const auto alias_topology = take_global(a::build_global_topology(
        alias_geometry, {fixture.roles, ownership, trailing}, cap, *global, memory));
    const auto alias_contact =
        take_global(a::check_global_geometry(alias_topology, alias_geometry, cap, *global, memory));
    auto alias_parent_points = growth_points;
    alias_parent_points[0].point = coincident;
    const cf::WakeMaterialSnapshot alias_window_parent{alias_contact, alias_parent_points};
    const auto alias_window = take(
        cf::prepare_wake_material_transition(window_si, window_second, alias_window_parent,
                                             window_trial, next_epoch, extended_limits, memory));
    REQUIRE(alias_window.retired_count() == 2);
    CHECK(take(alias_window.retired(0)).point == id[4]);

    const cf::WakeMaterialSnapshot partial_window{previous, material};
    reject(cf::prepare_wake_material_transition(window_si, window_second, partial_window,
                                                window_trial, next_epoch, extended_limits, memory),
           core::FailureCode::StaleGeneration);

    auto fractional_case = window_case;
    fractional_case.field("analysis")
        .field("wake")
        .field("retention")
        .set("duration", parse("0.1875"));
    const auto fractional_si = materialize(fractional_case, io::SourceId{488});
    const auto fractional_grid = take(cf::prepare_case_reference_clock(
        fractional_si, stamp, transport_clock_limits, 1000000, memory));
    const cf::WakeTransportClock fractional_clock{fractional_grid, fractional_grid, destination};
    reject(cf::prepare_wake_material_transition(fractional_si, fractional_clock, growth_trial,
                                                window_trial, next_epoch, extended_limits, memory),
           core::FailureCode::InvalidInput);

    reject(cf::prepare_wake_material_transition(window_si, window_second, growth_trial,
                                                extended_trial, next_epoch, extended_limits,
                                                memory),
           core::FailureCode::StaleGeneration);
    reject(cf::prepare_wake_material_transition(limited_si, limited_second, growth_trial,
                                                window_trial, next_epoch, extended_limits, memory),
           core::FailureCode::ResourceExceeded);

    auto two_rows_case = window_case;
    two_rows_case.field("analysis").field("wake").field("retention").set("duration", parse("0.25"));
    two_rows_case.field("numerics").field("wake_rows").set("maximum_rows", parse("2"));
    const auto two_rows_si = materialize(two_rows_case, io::SourceId{490});
    const auto two_rows_grid = take(cf::prepare_case_reference_clock(
        two_rows_si, stamp, transport_clock_limits, 1000000, memory));
    const cf::WakeTransportClock two_rows_clock{
        two_rows_grid, two_rows_grid, {stamp, 0, 3, p::RequestedSide::Automatic}};
    const Number third_time{11, 8};
    auto later_input = extended_input;
    later_input.context.time = 1.375;
    later_input.exact_time = third_time.view();
    const auto later_geometry = test_geometry::take(
        a::prepare_owned_geometry(later_input, test_geometry::limits(), *geometry_work, memory));
    const auto later_contact = take_global(
        a::check_global_geometry(extended_topology, later_geometry, cap, *global, memory));
    const a::ReleaseEpoch latest_epoch{10};
    auto shifted_points = extended_points;
    for (const auto i : {std::size_t{0}, std::size_t{2}}) {
        shifted_points[i].release.epoch = next_epoch;
        shifted_points[i + 1].release.epoch = latest_epoch;
    }
    const cf::WakeMaterialSnapshot shifted_trial{later_contact, shifted_points};
    const auto shifted = take(
        cf::prepare_wake_material_transition(two_rows_si, two_rows_clock, extended_trial,
                                             shifted_trial, latest_epoch, extended_limits, memory));
    REQUIRE(shifted.retired_count() == 2);
    CHECK(take(shifted.retired(0)).release == release);
    CHECK(take(shifted.retired(1)).release == second_release);
    const auto retained = take(shifted.row(0));
    if (!retained.parent) {
        FAIL("A retained window row must keep its parent material point");
        return;
    }
    CHECK(retained.parent->release.epoch == next_epoch);
    CHECK(retained.parent->point == id[4]);
    CHECK(retained.current.point == far_up);
    auto retired_wrong_row = shifted_points;
    retired_wrong_row[0].release = release;
    const cf::WakeMaterialSnapshot wrong_tail{later_contact, retired_wrong_row};
    reject(cf::prepare_wake_material_transition(two_rows_si, two_rows_clock, extended_trial,
                                                wrong_tail, latest_epoch, extended_limits, memory),
           core::FailureCode::StaleGeneration);

    const std::array birth_origin{cf::WakeMaterialPoint{release, release.point, std::nullopt}};
    auto from_origin = stages;
    from_origin[2].points = birth_origin;
    reject(run_material(from_origin), core::FailureCode::InvalidInput);
    from_origin[2].preceding_births = &birth;
    reject(run_material(from_origin), core::FailureCode::StaleGeneration);
    from_origin[2].preceding_births = &interval_birth;
    const auto after_birth = take(run_material(from_origin));
    const auto after_birth_row = take(after_birth.row(0));
    // Xn=2, Xn-1=15/8, and the only Xn-2 is the authenticated TE origin at x=1.
    residual_expected(after_birth_row.residual[0], Rat{-2});
    near(center(after_birth_row.residual[2]), -vz[1]);
    CHECK(take(after_birth.position_dual(0)).release == release);
    const std::array different_detached{cf::WakeMaterialPoint{release, id[5], a::FacetId{202}}};
    auto wrong_detached = from_origin;
    wrong_detached[1].points = different_detached;
    reject(run_material(wrong_detached), core::FailureCode::StaleGeneration);
    auto replaced_origin = from_origin;
    replaced_origin[2].points = stages[2].points;
    reject(run_material(replaced_origin), core::FailureCode::StaleGeneration);
    const std::array origin_as_current{from_origin[2], from_origin[1], from_origin[2]};
    reject(run_material(origin_as_current), core::FailureCode::InvalidInput);

    const auto stale_origin = history_contact(start_time, 1., 1.625, core::GeometryGeneration{7});
    const std::array wrong_origin_context{
        from_origin[0], from_origin[1],
        cf::WakeMaterialSnapshot{stale_origin, birth_origin, &interval_birth}};
    reject(run_material(wrong_origin_context), core::FailureCode::StaleGeneration);

    auto other_case = source_case;
    other_case.field("flow").set("uniform", parse(R"({"kind":"constant","value":[0,0,2]})"));
    const auto other_si = materialize(other_case, io::SourceId{485});
    const auto other_grid = take(
        cf::prepare_case_reference_clock(other_si, stamp, transport_clock_limits, 1000000, memory));
    const auto other_provider = prepare_gust(other_case, io::SourceId{485}, flow_workspace, memory);
    const cf::DiscreteKnownFlowInput other_flow{
        other_provider,   other_provider,         flow.time,      flow.projection,
        flow.flow_limits, flow.coordinate_limits, *geometry_work, flow_workspace};
    const cf::WakeTransportClock other_clock{other_grid, other_grid, destination};
    reject(cf::evaluate_fixed_free_wake_transport(
               other_clock, from_origin, input, other_flow, selected, &d, history_direction,
               transport_seeds, transport_limits, cap_velocity, shared, *local, *kernel, memory),
           core::FailureCode::StaleGeneration);

    auto event_case = source_case;
    event_case.field("analysis")
        .set(
            "intervals",
            parse(
                R"([{"start":0.875,"end":1.125,"step":0.25,"restart":"initial"},{"start":1.125,"end":1.375,"step":0.125,"restart":"event"}])"));
    const auto event_si = materialize(event_case, io::SourceId{484});
    const cf::CaseFlowClockLimits event_clock_limits{limits, 2, 4};
    const auto event_grid = take(
        cf::prepare_case_reference_clock(event_si, stamp, event_clock_limits, 1000000, memory));
    const auto event_provider = prepare_gust(event_case, io::SourceId{484}, flow_workspace, memory);
    const cf::WakeTransportClock event_clock{
        event_grid, event_grid, {stamp, 1, 1, p::RequestedSide::Automatic}};
    auto event_limits = transport_limits;
    event_limits.clock = event_clock_limits;
    const cf::DiscreteKnownFlowInput event_trial_flow{
        event_provider,   event_provider,         flow.time,      flow.projection,
        flow.flow_limits, flow.coordinate_limits, *geometry_work, flow_workspace};
    const cf::FirstFreeWakeEndpoint event_trial{input, event_trial_flow, selected, &d};
    const auto after_event = [&](p::RequestedSide side) {
        const cf::DiscreteNonpenetrationInput parent_source{
            previous, previous_bound, wake, input.history, {7, 10}, stamp, side};
        const cf::DiscreteKnownFlowInput parent_flow{
            event_provider,  event_provider,   {stamp, previous_time.view(), side},
            flow.projection, flow.flow_limits, flow.coordinate_limits,
            *geometry_work,  flow_workspace};
        const cf::FirstFreeWakeEndpoint parent{parent_source, parent_flow, selected, &d};
        return cf::evaluate_first_free_wake_transport(
            event_clock, event_trial, parent, std::nullopt, restarted, transport_seeds,
            event_limits, cap_velocity, shared, *local, *kernel, memory);
    };
    reject(after_event(p::RequestedSide::Left), core::FailureCode::StaleGeneration);
    const auto post_event = take(after_event(p::RequestedSide::Right));
    residual_expected(take(post_event.row(0)).residual[0], Rat{1});
    const auto& event_parent = take(post_event.parent_free_velocity()).get();
    exact(take(take(event_parent.ambient(0)).get().time()), previous_time);
    event_case.field("analysis")
        .field("wake")
        .set("retention",
             parse(R"({"kind":"closed_window","duration":0.25,"distance":1,"error_budget":0.01})"));
    const auto event_window_si = materialize(event_case, io::SourceId{489});
    const auto event_window_grid = take(cf::prepare_case_reference_clock(
        event_window_si, stamp, event_clock_limits, 1000000, memory));
    const cf::WakeTransportClock event_window_clock{event_window_grid, event_window_grid,
                                                    event_clock.destination};
    auto event_window_limits = event_limits;
    event_window_limits.maximum_particles = 2;
    reject(cf::prepare_wake_material_transition(event_window_si, event_window_clock, growth_trial,
                                                window_trial, next_epoch, event_window_limits,
                                                memory),
           core::FailureCode::Unsupported);

    auto prescribed_case = source(R"({"kind":"constant","value":[0,0,9]})",
                                  R"([{"start":1,"end":1.5,"step":0.125,"restart":"initial"}])");
    prescribed_case.field("units").set("length", string("m"));
    prescribed_case.field("units").set("time", string("s"));
    prescribed_case.field("analysis")
        .field("wake")
        .set(
            "convection",
            parse(
                R"({"kind":"prescribed","velocity":{"kind":"polynomial","origin":1,"scale":1,"domain":[1,1.5],"coefficients":[[1,0,0],[2,0,0]]}})"));
    const auto prescribed_si = materialize(prescribed_case, io::SourceId{483});
    const auto prescribed_grid = take(cf::prepare_case_reference_clock(
        prescribed_si, stamp, transport_clock_limits, 1000000, memory));
    const auto prescribed_source = take(cf::prepare_case_clocked_flows(
        prescribed_si, stamp, flow_workspace, transport_clock_limits, memory));
    const Number projection_one{1}, projection_zero{0};
    const cf::CaseFlowComponentBudget component{projection_one.view(), absolute.view(),
                                                projection_zero.view()};
    const std::array budgets{component, component, component};
    const cf::CaseFlowProjectionRequest projection{
        core::Precision::TwoWord, budgets, budgets, {4096, budget}};
    const auto prescribed_node = take(cf::evaluate_case_clocked_flows(
        prescribed_source, prescribed_source, destination, projection, projection, flow_workspace,
        transport_clock_limits, memory));
    const cf::WakeTransportClock prescribed_clock{prescribed_grid, prescribed_grid, destination};
    reject(cf::evaluate_fixed_free_wake_transport(
               prescribed_clock, stages, input, flow, selected, &d, history_direction,
               transport_seeds, transport_limits, cap_velocity, shared, *local, *kernel, memory),
           core::FailureCode::StaleGeneration);
    const auto prescribed_ambient =
        prepare_gust(prescribed_case, io::SourceId{483}, flow_workspace, memory);
    const cf::DiscreteKnownFlowInput wrong_mode{
        prescribed_ambient, prescribed_ambient,     flow.time,        flow.projection,
        flow.flow_limits,   flow.coordinate_limits, flow.coordinates, flow.workspace};
    reject(cf::evaluate_fixed_free_wake_transport(
               prescribed_clock, stages, input, wrong_mode, selected, &d, history_direction,
               transport_seeds, transport_limits, cap_velocity, shared, *local, *kernel, memory),
           core::FailureCode::Unsupported);
    auto short_transport = transport_limits;
    short_transport.maximum_steps = 1;
    reject(cf::evaluate_fixed_free_wake_transport(
               transport_clock, stages, input, flow, selected, &d, history_direction,
               transport_seeds, short_transport, cap_velocity, shared, *local, *kernel, memory),
           core::FailureCode::ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
    std::stop_source cancellation;
    REQUIRE(cancellation.request_stop());
    reject(cf::evaluate_fixed_free_wake_transport(transport_clock, stages, input, flow, selected,
                                                  &d, history_direction, transport_seeds,
                                                  transport_limits, cap_velocity, shared, *local,
                                                  *kernel, memory, cancellation.get_token()),
           core::FailureCode::Cancelled);
    const std::array<cf::NonpenetrationVector, 1> dcurrent{dx},
        prescribed_seeds{coordinates(2, 0, 0)};
    const cf::PrescribedWakeTransportDirection direction{dcurrent, history_direction,
                                                         a::SingleWord{.25}};
    const cf::PrescribedWakeClockSample clock_sample{std::cref(prescribed_node)};
    reject(cf::evaluate_fixed_prescribed_wake_transport(transport_clock, stages, clock_sample,
                                                        &direction, prescribed_seeds,
                                                        transport_limits, memory),
           core::FailureCode::StaleGeneration);
    const auto prescribed = take(cf::evaluate_fixed_prescribed_wake_transport(
        prescribed_clock, stages, clock_sample, &direction, prescribed_seeds, transport_limits,
        memory));
    const auto row = take(prescribed.row(0));
    CHECK(row.release == release);
    residual_expected(row.residual[0], Rat{-1});
    residual_expected(row.residual[2], Rat{0});
    if (!row.prescribed_position || !row.direction || !row.prescribed_position_direction) {
        FAIL("Expected the prescribed position and both transport directions");
        return;
    }
    residual_expected(row.prescribed_position.value()[0], Rat{25, 12});
    residual_expected(row.direction.value()[0], Rat{-1});
    residual_expected(row.prescribed_position_direction.value()[0], Rat{5, 24});
    const auto prescribed_dual = take(prescribed.position_dual(0));
    residual_expected(prescribed_dual.stencil[0][0], Rat{24});
    residual_expected(prescribed_dual.stencil[1][0], Rat{-32});
    residual_expected(prescribed_dual.stencil[2][0], Rat{8});
    residual_expected(take(prescribed.prescribed_time_dual()), Rat{-4});
    WorkEnclosure prescribed_lhs, prescribed_rhs;
    prescribed_lhs.add(row.direction.value(), prescribed_seeds[0]);
    prescribed_rhs.add(prescribed_dual.stencil[0], dcurrent[0]);
    prescribed_rhs.add(prescribed_dual.stencil[1], dprevious[0]);
    prescribed_rhs.add(prescribed_dual.stencil[2], dolder[0]);
    prescribed_rhs.add(take(prescribed.prescribed_time_dual()), direction.flow_time);
    CHECK(abs(prescribed_lhs.value - prescribed_rhs.value) <=
          prescribed_lhs.radius + prescribed_rhs.radius);
}

TEST_CASE("Prescribed aero attachment feeds one pose and velocity into discrete nonpenetration",
          "[AE-07][IO4-01][UW-D06][UW-D07][N4-06][verification][negative]") {
    auto input = root("aerodynamic_prescribed_transient");
    input.field("units").set("time", string("s"));
    input.field("units").set("length", string("m"));
    input.field("units").set("angle", string("deg"));
    input.field("analysis")
        .set(
            "intervals",
            parse(
                R"([{"start":0,"end":0.1,"step":0.1,"restart":"initial"},{"start":0.1,"end":1,"step":0.1,"restart":"event"}])"));
    auto& motion =
        input.field("model").field("aero_surfaces").elements[0].field("attachment").field("motion");
    motion.set("origin", parse(R"({"origin":[0,0,0],"rotation":[[1,0,0],[0,1,0],[0,0,1]]})"));
    motion.set(
        "translation",
        parse(
            R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[0,0,0],[1,0,0]]})"));
    motion.set(
        "rotation_vector",
        parse(
            R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[0,0,0],[600,0,0]]})"));
    motion.set("events", parse("[0,0.1,1]"));

    core::ResourcePlan memory{128ULL * 1024 * 1024};
    p::CurveWorkspace curves;
    const auto different = take(
        structure::ExactRational::from_decimal("0.100000000000000001", exact_limits(), memory));
    const auto different_words = take(io::detail::copy_rational_words(different, memory, {}));
    const auto different_time = different_words.view();
    const auto geometry_caps = test_geometry::limits();
    const a::GlobalLimits contact_caps{std::numeric_limits<std::uint64_t>::max() / 4,
                                       numerics::linear::EstimateAdmission{4096}};
    const cf::PrescribedAeroLimits caps{
        {{exact_limits(), 10000, cf::case_motion_frame_interval_scratch_bytes},
         {exact_limits(),
          {{128, 100000, budget, 10000}, numerics::linear::EstimateAdmission{budget}},
          100,
          10000},
         4096,
         100000},
        geometry_caps,
        contact_caps,
        UINT64_C(10000000000)};
    const auto selected_attachment = [](const io::SiCase& si) {
        const auto top = take(take(si.source_case()).get().root());
        const auto attached =
            take(top.model.record().aero_surfaces.at(0)).record().attachment.record();
        const auto* selected = std::get_if<io::schema::AeroAttachmentPrescribed>(&attached);
        if (selected == nullptr)
            FAIL("Missing checked prescribed attachment");
        return *selected;
    };
    const Number unit{1}, exact_instant{1, 10}, abs_budget{1, 1000000000}, zero{0};
    std::array<cf::CaseFlowComponentBudget, 24> budget_rows{};
    budget_rows.fill({unit.view(), abs_budget.view(), zero.view()});
    const cf::PrescribedAeroProjection projection{
        core::Precision::Binary64, budget_rows,
        std::span<const cf::CaseFlowComponentBudget>{budget_rows}.first<6>()};

    ResidualFixture fixture;
    const auto& id = fixture.id;
    const a::PointId alias{7};
    const std::array ids{id[0], id[1], id[2], id[3], id[4], alias};
    auto alias_binding = test_geometry::binding();
    alias_binding.ownership = a::AliasOf{id[3]};
    const auto& binding = fixture.bindings;
    const std::array bindings{binding[0], binding[1], binding[2],
                              binding[3], binding[4], alias_binding};
    const std::array<double, 18> coordinates0{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 2, 1, 0, 1, 0, 0};
    const std::array moving{id[0], id[1], id[2], id[3]};
    fixture.quads[0].points[3] = alias;
    fixture.edges[2].points[1] = alias;
    fixture.edges[3].points[0] = alias;
    auto original = fixture.input();
    original.points = ids;
    original.bindings = bindings;
    original.reference = a::SinglePositions{coordinates0};
    original.current = a::SinglePositions{coordinates0};
    original.edges = std::span{fixture.edges}.first<5>();
    original.faces = std::span{fixture.faces}.first<1>();
    original.bands = std::span{fixture.bands}.first<1>();
    original.quads = std::span{fixture.quads}.first<1>();
    original.context.time = 0.1;
    original.exact_time = exact_instant.view();
    auto geometry_work = test_geometry::workspace(memory);
    auto global = take(a::GlobalWorkspace::create(memory, contact_caps.metadata));
    auto old_words = coordinates0;
    for (const auto row : std::array<std::size_t, 5>{0, 1, 2, 3, 5})
        old_words[3 * row] += 4;
    auto previously_moved = original;
    previously_moved.current = a::SinglePositions{old_words};
    const auto reference_geometry = test_geometry::take(
        a::prepare_owned_geometry(previously_moved, geometry_caps, *geometry_work, memory));
    const auto topology = take_global(a::build_global_topology(
        reference_geometry,
        {std::span{fixture.roles}.first<1>(), std::span{fixture.edge_owners}.first<5>(), {}},
        contact_caps, *global, memory));
    const std::array rows{cf::PrescribedAeroRow{
        a::QueryId{22}, a::SourceFacePoint{a::FacetId{201}, a::SingleWord{.25}, a::SingleWord{.5}},
        a::SingleWord{0}, a::SingleWord{0}}};
    const auto side = p::RequestedSide::Left;
    const cf::PrescribedAeroReference model{topology, original, moving, rows};
    a::UniformFlowWorkspace flow_workspace;
    const cf::CaseClockedPrescribedLimits clock_caps{
        {{limits, 2, 10}, gust_limits}, caps.motion, budget};
    const auto prepare_clocked = [&](const JsonToken& json) {
        const auto si = materialize(json, io::SourceId{79});
        return take(cf::prepare_case_clocked_prescribed_inputs(
            si, selected_attachment(si), stamp, flow_workspace, curves, clock_caps, memory));
    };
    const auto clock = prepare_clocked(input);
    const auto current_clock = prepare_clocked(input);
    const cf::CaseFlowClockQuery arrival{stamp, 0, 1, side};
    const auto run_clocked = [&](const cf::CaseClockedPrescribedInputs& live,
                                 cf::CaseFlowClockQuery query,
                                 const cf::PrescribedAeroReference& reference,
                                 const cf::PrescribedAeroProjection& selected,
                                 const cf::PrescribedAeroLimits& cap, core::ResourcePlan& plan) {
        return cf::evaluate_case_clocked_prescribed_attachment(
            clock, live, query, reference, selected, curves, *geometry_work, *global,
            clock_caps.flow.clock_and_wake, cap, plan);
    };
    const auto state = take(run_clocked(current_clock, arrival, model, projection, caps, memory));
    reject(run_clocked(current_clock, {stamp, 0, 1, p::RequestedSide::Right}, model, projection,
                       caps, memory),
           core::FailureCode::InvalidInput);
    const auto& sampled = take(state.motion()).get();
    const auto& instant = take(sampled.time()).get();
    exact(instant.time(), exact_instant);
    CHECK(instant.requested_side() == side);
    CHECK(state.collocation_count() == 1);
    const auto row = take(state.collocation(0));
    const auto& contact = take(state.contact()).get();
    const auto& geometry = take_global(contact.geometry()).get();
    auto larger_geometry_cap = caps;
    larger_geometry_cap.geometry.max_steps = caps.work;
    const auto same_trial =
        take(run_clocked(current_clock, arrival, model, projection, larger_geometry_cap, memory));
    const auto& same_geometry = take_global(take(same_trial.contact()).get().geometry()).get();
    CHECK(test_geometry::take(same_geometry.identity()) ==
          test_geometry::take(geometry.identity()));
    const auto same_row = take(same_trial.collocation(0));
    for (std::size_t axis = 0; axis < row.velocity.size(); ++axis)
        CHECK(sum_words(same_row.velocity[axis]) == sum_words(row.velocity[axis]));
    CHECK(test_geometry::take(geometry.owner(alias)) == id[3]);
    const auto alias_position =
        test_geometry::take(geometry.point_words(alias, a::Configuration::Current));
    const auto owner_position =
        test_geometry::take(geometry.point_words(id[3], a::Configuration::Current));
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(sum_words(alias_position[k]) == sum_words(owner_position[k]));
        const auto aliased = take(state.material_point(alias, k, p::JetOrder::First));
        const auto owned = take(state.material_point(id[3], k, p::JetOrder::First));
        const auto* alias_word = std::get_if<io::ScaledBinary64>(&aliased);
        const auto* owner_word = std::get_if<io::ScaledBinary64>(&owned);
        if (alias_word == nullptr || owner_word == nullptr)
            FAIL("Expected binary64 motion words");
        CHECK(alias_word->center == owner_word->center);
        CHECK(alias_word->exponent2 == owner_word->exponent2);
    }
    reject(state.material_point(id[4], 0, p::JetOrder::First), core::FailureCode::InvalidInput);
    const auto inherited =
        test_geometry::take(geometry.point_words(id[4], a::Configuration::Current));
    CHECK(sum_words(inherited[0]) == 2);
    CHECK(sum_words(inherited[1]) == 1);
    const auto tol = rational_word(1e-12);
    const auto near = [&](const Rat& actual, double expected) {
        CHECK(abs(actual - rational_word(expected)) < tol);
    };
    const auto irrational = 1.7320508075688772935;
    const auto angular_rate = 10.471975511965977462;
    near(sum_words(row.velocity[0]), 1);
    near(sum_words(row.velocity[1]), -angular_rate * irrational / 4);
    near(sum_words(row.velocity[2]), angular_rate / 4);
    const auto xy = test_geometry::take(
        a::exact_query_coordinates(geometry, row.query, geometry_caps, *geometry_work));
    const auto denom = two_power(static_cast<int>(a::ExactQueryCoordinates::denominator_shift));
    near(Rat{xy.numerator[0]} / denom, .85);
    near(Rat{xy.numerator[1]} / denom, .25);
    near(Rat{xy.numerator[2]} / denom, irrational / 4);

    const auto zero_flow = source(R"({"kind":"constant","value":[0,0,0]})", one_step);
    const auto ambient = import(zero_flow, io::SourceId{81}, flow_workspace);
    const auto ambient_current = import(zero_flow, io::SourceId{82}, flow_workspace);
    const std::array budgets{budget_rows[0], budget_rows[1], budget_rows[2]};
    const cf::CaseFlowProjectionRequest flow_projection{
        core::Precision::TwoWord, budgets, budgets, {4096, budget}};
    const auto collocations = std::array{row};
    const cf::DiscreteNonpenetrationInput view{
        contact, collocations, {}, {core::AcceptedStateId{2}, {}}, {7, 11}, stamp, side};
    const cf::DiscreteNonpenetrationLimits residual_limits{
        {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10},
        {geometry_caps,
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(6000000000)};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const auto residual = [&](p::RequestedSide requested) {
        auto flow = take(cf::project_case_uniform_flow(
            evaluate(ambient, ambient_current, flow_workspace, exact_instant, requested),
            flow_projection));
        return cf::evaluate_discrete_nonpenetration(view, std::move(flow), residual_limits, shared,
                                                    *local, *kernel, memory);
    };
    reject(residual(p::RequestedSide::Right), core::FailureCode::StaleGeneration);
    const auto result = take(residual(side));
    const auto value = take(result.row(0)).residual;
    const auto center = sum_words(value.center()) * two_power(value.exponent());
    near(center, -angular_rate / 2);

    const std::array<double, 18> zero_low{};
    auto wide_original = original;
    wide_original.current = a::DoublePositions{coordinates0, zero_low};
    const auto wide = take(
        run_clocked(current_clock, arrival, {topology, wide_original, moving, rows},
                    {core::Precision::TwoWord, projection.pose, projection.point}, caps, memory));
    const auto& wide_motion = take(wide.motion()).get();
    const auto& wide_time = take(wide_motion.time()).get();
    exact(wide_time.time(), exact_instant);
    CHECK(wide_time.requested_side() == side);
    const auto wide_row = take(wide.collocation(0));
    CHECK(std::holds_alternative<a::DoubleWord>(wide_row.velocity[2]));
    near(sum_words(wide_row.velocity[2]), angular_rate / 4);
    auto wrong_time = original;
    wrong_time.exact_time =
        p::RationalView{different_time.sign, different_time.numerator, different_time.denominator};
    reject(run_clocked(current_clock, arrival, {topology, wrong_time, moving, rows}, projection,
                       caps, memory),
           core::FailureCode::StaleGeneration);
    const std::array invalid_points{id[0], id[1], id[2], alias};
    reject(run_clocked(current_clock, arrival, {topology, original, invalid_points, rows},
                       projection, caps, memory),
           core::FailureCode::InvalidInput);
    const std::array extra_point{id[0], id[1], id[2], id[3], id[4]};
    reject(run_clocked(current_clock, arrival, {topology, original, extra_point, rows}, projection,
                       caps, memory),
           core::FailureCode::InvalidInput);
    auto changed_input = input;
    changed_input.field("model")
        .field("aero_surfaces")
        .elements[0]
        .field("attachment")
        .field("motion")
        .field("translation")
        .set("coefficients", parse("[[0,0,0],[2,0,0]]"));
    const auto changed_clock = prepare_clocked(changed_input);
    reject(run_clocked(changed_clock, arrival, model, projection, caps, memory),
           core::FailureCode::StaleGeneration);
    auto no_work = caps;
    no_work.work = 0;
    core::ResourcePlan refused{0};
    reject(run_clocked(current_clock, arrival, model, projection, no_work, refused),
           core::FailureCode::ResourceExceeded);
    CHECK(refused.allocated_bytes() == 0);
    CHECK(test_geometry::take(geometry.owner(alias)) == id[3]);
}

TEST_CASE("PG structural attachment retains a curved material point and transfers work",
          "[AE-02][UW-D01][UW-D06][SH-D02][N4-06][verification][negative]") {
    namespace st = uvlas::structure;
    namespace pg = pg_conditioning_test;
    pg::Fixture f;
    for (auto& thickness : f.heights)
        thickness = f.r("2");
    const auto one_metre = f.r("1");
    core::Workspace reference_scratch{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan memory{pg::profile};
    const auto pg_limits = pg::query_limits();
    auto basis = pg::make_basis(f, one_metre, reference_scratch, memory);
    auto gram = take(st::PgFlatGram::prepare(std::move(basis), 256, {pg_limits, pg::steps},
                                             reference_scratch, memory));
    auto certificate = take(std::move(gram).generate_and_certify(
        {320, 246}, {pg::profile, pg::steps, {}}, reference_scratch, memory));
    const auto patch = take(
        st::PgPreparedPatch::prepare(std::move(certificate), pg_limits, reference_scratch, memory));
    const auto curve =
        take(patch.conditioning().required_coordinates({19}, pg_limits, reference_scratch, memory));
    const auto& partition = take(take(patch.conditioning().reference_basis())->source()).partition;
    const auto nodes = take(partition.vertices());
    const auto triangles = take(partition.triangles());
    const auto faces = take(partition.horizontal_faces());
    const auto layers = take(partition.layers());
    const auto host = [&](const pg::Q& s, const pg::Q& t) {
        for (std::size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            if (face.side != st::PgPartitionFaceSide::Upper ||
                pg::q(layers[static_cast<std::size_t>(face.layer)].upper) != 1)
                continue;
            const auto& triangle = triangles[static_cast<std::size_t>(face.triangle)];
            const auto& o = nodes[static_cast<std::size_t>(triangle.vertices[0])];
            const auto& b = nodes[static_cast<std::size_t>(triangle.vertices[1])];
            const auto& c = nodes[static_cast<std::size_t>(triangle.vertices[2])];
            const pg::Q bx = pg::q(b.s) - pg::q(o.s), by = pg::q(b.t) - pg::q(o.t);
            const pg::Q cx = pg::q(c.s) - pg::q(o.s), cy = pg::q(c.t) - pg::q(o.t);
            const pg::Q x = s - pg::q(o.s), y = t - pg::q(o.t);
            const pg::Q det = bx * cy - cx * by;
            REQUIRE(det > 0);
            const pg::Q u = (x * cy - cx * y) / det, v = (bx * y - x * by) / det;
            if (u >= 0 && v >= 0 && u + v <= 1)
                return std::tuple{static_cast<core::Index>(i),
                                  pg::rat(f, numerator(u), denominator(u)),
                                  pg::rat(f, numerator(v), denominator(v))};
        }
        throw std::runtime_error("fixture material point missing");
    };
    const auto first = host(0, 0), second = host(0, 1), third = host(1, 1), fourth = host(1, 0),
               control = host(pg::Q{3, 4}, pg::Q{3, 8});
    const auto point = [](const auto& x) -> st::PgBoundaryPositionPoint {
        return st::PgReferenceHorizontalPoint{std::get<0>(x), std::get<1>(x), std::get<2>(x)};
    };
    const std::array material{point(first), point(second), point(third), point(fourth),
                              point(control)};
    const auto columns = patch.reference_coordinates();
    std::vector<double> q_high(3 * columns.size()), q_low(q_high.size()), v_high(q_high.size()),
        v_low(q_high.size());
    const auto words = [](const Rat& exact, double& high, double& low) {
        high = exact.convert_to<double>();
        low = (exact - rational_word(high)).convert_to<double>();
    };
    const auto as_rat = [](const pg::Q& q) {
        return Rat{Big{numerator(q).str()}, Big{denominator(q).str()}};
    };
    REQUIRE(curve.size() == columns.size());
    // Constant z curvature on this flat patch integrates to s²/2.
    // The current top surface has X_z=1+s²/2 and V_z=s².
    for (std::size_t i = 0; i < columns.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            const Rat change = k == 2 ? as_rat(pg::q(curve.values()[i])) : Rat{0};
            words(as_rat(pg::q(columns[i][k])) + change, q_high[3 * i + k], q_low[3 * i + k]);
            words(2 * change, v_high[3 * i + k], v_low[3 * i + k]);
        }

    ResidualFixture fixture;
    const std::array<double, 15> flat{0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1};
    const std::array<double, 15> low{};
    fixture.bindings[4].ownership = a::AliasOf{fixture.id[2]};
    fixture.quads[0].points[2] = fixture.id[4];
    auto reference = fixture.input();
    reference.points = reference.points.first(5);
    reference.bindings = reference.bindings.first(5);
    reference.reference = a::DoublePositions{flat, low};
    reference.current = reference.reference;
    reference.edges = reference.edges.first(5);
    reference.faces = reference.faces.first(1);
    reference.bands = reference.bands.first(1);
    reference.quads = reference.quads.first(1);
    const a::GlobalLimits global_limits{UINT64_C(5000000000),
                                        numerics::linear::EstimateAdmission{65536}};
    const auto topo = fixture.topology();
    const a::GlobalTopologyInput topology_input{topo.roles.first(1), topo.edges.first(5), {}};
    auto geometry_work = test_geometry::workspace(memory);
    auto global = take(a::GlobalWorkspace::create(memory, global_limits.metadata));
    const auto geometry = test_geometry::take(
        a::prepare_owned_geometry(reference, test_geometry::limits(), *geometry_work, memory));
    const auto topology = take_global(
        a::build_global_topology(geometry, topology_input, global_limits, *global, memory));
    const std::array roles{cf::PgAeroCollocation{a::FacetId{200}, a::QueryId{501}}};
    constexpr auto format = core::Precision::TwoWord;
    constexpr st::PgPreparedReferenceProjection projection{format, 1e-26, 1e-26,
                                                           1e-24,  1e-24, 1e-24};
    const auto prepare = [&](const auto& requests) {
        return cf::PgAeroMap::prepare(patch, reference.points.first(4), requests, material,
                                      projection, 1e-24, pg_limits, memory);
    };
    const auto map = take(prepare(roles));
    const cf::PgAeroState state{take(patch.info()).identity,
                                reference.context,
                                {format, q_high, q_low},
                                {format, v_high, v_low},
                                {501, 502}};
    const std::array edges{a::EdgeId{10}, a::EdgeId{11}, a::EdgeId{12}, a::EdgeId{13}};
    const cf::PgAeroChart chart{topology, reference, edges};
    const cf::PgAeroLimits kinematic_limits{{100000, 1e-21, {}},
                                            {100000, 1e-21, {}},
                                            {100000, 1e-20, {}},
                                            test_geometry::limits(),
                                            global_limits,
                                            UINT64_C(6000000000),
                                            1e-20,
                                            1e-20};
    auto storage = take(memory.allocate<std::byte>(4096, core::Operation::EvaluateAero));
    core::Workspace scratch{std::move(storage)};
    const auto run = [&](const cf::PgAeroMap& mapping, const cf::PgAeroState& current) {
        return cf::evaluate_pg_aero_trial(mapping, current, chart, kinematic_limits, scratch,
                                          *geometry_work, *global, memory);
    };
    auto wrong_patch = state;
    wrong_patch.patch.bytes.front() ^= std::byte{1};
    reject(run(map, wrong_patch), core::FailureCode::StaleGeneration);
    auto wrong_state = state;
    wrong_state.context.current = core::GeometryGeneration{999};
    reject(run(map, wrong_state), core::FailureCode::StaleGeneration);
    auto wider = kinematic_limits;
    wider.position.maximum_steps = wider.maximum_steps;
    wider.velocity.maximum_steps = wider.maximum_steps;
    wider.tangents.maximum_steps = wider.maximum_steps;
    auto wider_geometry = kinematic_limits;
    wider_geometry.geometry.max_steps = wider_geometry.maximum_steps;
    const auto check_with_limits = [&](const cf::PgAeroMap& mapping,
                                       const cf::PgAeroChart& bound_chart,
                                       const cf::PgAeroTrial& initial,
                                       const cf::PgAeroLimits& adjusted) {
        const auto again = take(cf::evaluate_pg_aero_trial(
            mapping, state, bound_chart, adjusted, scratch, *geometry_work, *global, memory));
        const auto& initial_geometry = take_global(take(initial.contact()).get().geometry()).get();
        const auto& again_geometry = take_global(take(again.contact()).get().geometry()).get();
        CHECK(test_geometry::take(again_geometry.identity()) ==
              test_geometry::take(initial_geometry.identity()));
        const auto baseline_row = take(initial.collocation(0));
        const auto wider_row = take(again.collocation(0));
        const auto first_position = test_geometry::take(a::exact_query_coordinates(
            initial_geometry, baseline_row.query, test_geometry::limits(), *geometry_work));
        const auto second_position = test_geometry::take(a::exact_query_coordinates(
            again_geometry, wider_row.query, test_geometry::limits(), *geometry_work));
        CHECK(first_position.numerator == second_position.numerator);
        const auto& first_frame = std::get<cf::ExplicitCollocationTangents>(baseline_row.tangents);
        const auto& second_frame = std::get<cf::ExplicitCollocationTangents>(wider_row.tangents);
        for (std::size_t k = 0; k < 3; ++k) {
            CHECK(sum_words(baseline_row.velocity[k]) == sum_words(wider_row.velocity[k]));
            CHECK(sum_words(first_frame[0][k]) == sum_words(second_frame[0][k]));
            CHECK(sum_words(first_frame[1][k]) == sum_words(second_frame[1][k]));
        }
    };
    auto short_position = kinematic_limits;
    short_position.position.maximum_steps = 1;
    reject(cf::evaluate_pg_aero_trial(map, state, chart, short_position, scratch, *geometry_work,
                                      *global, memory),
           core::FailureCode::ResourceExceeded);
    const auto trial = take(run(map, state));
    check_with_limits(map, chart, trial, wider);
    check_with_limits(map, chart, trial, wider_geometry);
    const auto& contact = take(trial.contact()).get();
    const auto& current = take_global(contact.geometry()).get();
    REQUIRE(trial.collocation_count() == 1);
    CHECK(test_geometry::take(current.owner(fixture.id[4])) == fixture.id[2]);
    const auto alias_xyz =
        test_geometry::take(current.point_words(fixture.id[4], a::Configuration::Current));
    const auto owner_xyz =
        test_geometry::take(current.point_words(fixture.id[2], a::Configuration::Current));
    for (std::size_t axis = 0; axis < 3; ++axis)
        CHECK(sum_words(alias_xyz[axis]) == sum_words(owner_xyz[axis]));
    const std::array mismatched_points{material[1], material[0], material[2], material[3],
                                       material[4]};
    const auto wrong_reference =
        take(cf::PgAeroMap::prepare(patch, reference.points.first(4), roles, mismatched_points,
                                    projection, 1e-24, pg_limits, memory));
    reject(run(wrong_reference, state), core::FailureCode::GeometryInvalid);
    auto no_work = kinematic_limits;
    no_work.maximum_steps = 0;
    reject(cf::evaluate_pg_aero_trial(map, state, chart, no_work, scratch, *geometry_work, *global,
                                      memory),
           core::FailureCode::ResourceExceeded);
    const std::array rows{take(trial.collocation(0))};
    const auto exact = test_geometry::take(a::exact_query_coordinates(
        current, rows[0].query, test_geometry::limits(), *geometry_work));
    REQUIRE(exact.point_count == 0);
    const auto position = [&](std::size_t k) {
        return Rat{exact.numerator[k].convert_to<Big>()} *
               two_power(-static_cast<int>(a::ExactQueryCoordinates::denominator_shift));
    };
    CHECK(abs(position(0) - Rat{3, 4}) < rational_word(1e-21));
    CHECK(abs(position(1) - Rat{3, 8}) < rational_word(1e-21));
    CHECK(abs(position(2) - Rat{41, 32}) < rational_word(1e-21));
    CHECK(abs(sum_words(rows[0].velocity[2]) - Rat{9, 16}) < rational_word(1e-21));
    const auto& tangent = std::get<cf::ExplicitCollocationTangents>(rows[0].tangents);
    const auto normal = [&](std::size_t xyz) -> Rat {
        const auto j = (xyz + 1) % 3, k = (xyz + 2) % 3;
        return sum_words(tangent[0][j]) * sum_words(tangent[1][k]) -
               sum_words(tangent[0][k]) * sum_words(tangent[1][j]);
    };
    REQUIRE(normal(2) > 0);
    CHECK(abs(normal(0) / normal(2) + Rat{3, 4}) < rational_word(1e-20));
    CHECK(abs(normal(1) / normal(2)) < rational_word(1e-20));
    const std::array midpoint{take(trial.midpoint(0)), take(trial.midpoint(1)),
                              take(trial.midpoint(2)), take(trial.midpoint(3))};
    for (std::size_t i = 0; i < midpoint.size(); ++i) {
        constexpr std::array expected{0., .5, 1., .5};
        CHECK(midpoint[i].edge == edges[i]);
        CHECK(abs(sum_words(midpoint[i].velocity[2]) - rational_word(expected[i])) <
              rational_word(1e-20));
    }

    const auto source_case = source(R"({"kind":"constant","value":[0,0,2]})", one_step);
    a::UniformFlowWorkspace flow_workspace;
    const auto provider = import(source_case, io::SourceId{501}, flow_workspace);
    const Number one{1}, zero{0}, tolerance{1, 1000000000}, time{5, 4};
    const cf::CaseFlowComponentBudget component{one.view(), tolerance.view(), zero.view()};
    const std::array budgets{component, component, component};
    auto flow =
        take(cf::project_case_uniform_flow(evaluate(provider, provider, flow_workspace, time),
                                           {format, budgets, budgets, {4096, budget}}));
    const cf::DiscreteNonpenetrationInput input{contact,
                                                rows,
                                                {},
                                                {core::AcceptedStateId{2}, {}},
                                                {501, 502},
                                                stamp,
                                                p::RequestedSide::Automatic};
    const std::array<ls::UnitDimensions, 1> gamma_units{ls::UnitDimensions{2, 0, -1, 0, 0, 0, 0}};
    const std::array<ls::UnitDimensions, 1> speed_units{ls::UnitDimensions{1, 0, -1, 0, 0, 0, 0}};
    const std::array scales{1.};
    const auto request = circulation_request(
        {scales, scales, gamma_units, speed_units, circulation_context(current.context())}, 1,
        ls::Norm::Infinity);
    const std::array rates{cf::BoundCirculationRate{fixture.faces[0].id, a::SingleWord{.25}}};
    const cf::DiscreteVectorLoadInput load_input{
        input.trial_stamp,   request.retained_generation, midpoint, rates,
        a::SingleWord{1.25}, coordinates(0, 0, 0)};
    const cf::DiscreteNonpenetrationLimits load_limits{
        {format, 1e-10},
        {format, 1e-10, 1e-10},
        {test_geometry::limits(),
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(5000000000)},
        100000,
        UINT64_C(24000000000)};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const auto loads = take(cf::solve_discrete_circulation_loads(
        input, load_input, std::move(flow), load_limits, request, *geometry_work, shared, *local,
        *kernel, memory));
    REQUIRE(loads.point_count() == 4);
    const std::array forces{take(loads.point(0)), take(loads.point(1)), take(loads.point(2)),
                            take(loads.point(3))};
    auto different_stamp = state;
    different_stamp.stamp.generation++;
    const auto other_trial = take(run(map, different_stamp));
    reject(other_trial.pullback_nodal_forces(loads, format, 1e-16, 1000000, memory),
           core::FailureCode::StaleGeneration);
    auto shifted_q = q_high;
    for (std::size_t i = 2; i < shifted_q.size(); i += 3)
        shifted_q[i] += v_high[i] * .125;
    auto same_stamp_new_geometry = state;
    same_stamp_new_geometry.position.high = shifted_q;
    const auto shifted_trial = take(run(map, same_stamp_new_geometry));
    reject(shifted_trial.pullback_nodal_forces(loads, format, 1e-16, 1000000, memory),
           core::FailureCode::StaleGeneration);
    const auto structural_forces =
        take(trial.pullback_nodal_forces(loads, format, 1e-16, 1000000, memory));
    REQUIRE(structural_forces.size() == columns.size());
    WorkEnclosure structural, fluid;
    for (std::size_t i = 0; i < columns.size(); ++i) {
        const auto& result = structural_forces.values()[i];
        if (!result)
            FAIL("missing column force");
        const auto& force = result.value();
        for (std::size_t k = 0; k < 3; ++k)
            structural.add(force[k], a::DoubleWord{v_high[3 * i + k], v_low[3 * i + k]});
    }
    for (std::size_t i = 0; i < forces.size(); ++i)
        fluid.add(forces[i].force[2], a::SingleWord{i >= 2 ? 1. : 0.});
    CHECK(abs(fluid.value) > rational_word(.001));
    CHECK(abs(structural.value - fluid.value) <=
          structural.radius + fluid.radius + rational_word(1e-20));

    // The composed material chart crosses the host's reference midsurface at
    // the collocation. Its zeta derivative adds 1/4 to the surface slope.
    const auto volumes = take(partition.volumes());
    const auto inner = [&](const auto& location,
                           const st::ExactRational& depth) -> st::PgReferencePoint {
        const auto& selected = faces[static_cast<std::size_t>(std::get<0>(location))];
        const auto volume = std::ranges::find_if(volumes, [&](const auto& v) {
            return v.triangle == selected.triangle && v.layer == selected.layer;
        });
        REQUIRE(volume != volumes.end());
        return st::PgReferenceVolumePoint{static_cast<core::Index>(volume - volumes.begin()),
                                          std::get<1>(location), std::get<2>(location), depth};
    };
    const auto near_depth = f.r("7", "64"), far_depth = f.r("19", "64"),
               control_depth = f.r("1", "4");
    const std::array interior_points{inner(first, near_depth), inner(second, near_depth),
                                     inner(third, far_depth), inner(fourth, far_depth),
                                     inner(control, control_depth)};
    st::PgMaterialDirections directions{};
    directions[0][0] = f.r("1");
    directions[0][2] = f.r("1", "4");
    directions[1][1] = f.r("1");
    const auto material_map = take(cf::PgAeroMap::prepare_material(
        patch, reference.points.first(4), roles, interior_points, std::span{&directions, 1},
        projection, 1e-24, pg_limits, memory));
    constexpr double near_z = -3. / 16.;
    constexpr double far_z = 1. / 16.;
    const std::array<double, 15> mid_grid{0,     0, near_z, 0,     1, near_z, 1,    1,
                                          far_z, 1, 0,      far_z, 1, 1,      far_z};
    auto material_geometry = reference;
    material_geometry.reference = a::DoublePositions{mid_grid, low};
    material_geometry.current = material_geometry.reference;
    const auto mid_base = test_geometry::take(a::prepare_owned_geometry(
        material_geometry, test_geometry::limits(), *geometry_work, memory));
    const auto mid_topology = take_global(
        a::build_global_topology(mid_base, topology_input, global_limits, *global, memory));
    const cf::PgAeroChart mid_chart{mid_topology, material_geometry, edges};
    const auto mid_trial =
        take(cf::evaluate_pg_aero_trial(material_map, state, mid_chart, kinematic_limits, scratch,
                                        *geometry_work, *global, memory));
    check_with_limits(material_map, mid_chart, mid_trial, wider);
    check_with_limits(material_map, mid_chart, mid_trial, wider_geometry);
    const std::array mid_rows{take(mid_trial.collocation(0))};
    const auto& mid_contact = take(mid_trial.contact()).get();
    const auto& mid_geometry = take_global(mid_contact.geometry()).get();
    const auto mid_query = test_geometry::take(a::exact_query_coordinates(
        mid_geometry, mid_rows[0].query, test_geometry::limits(), *geometry_work));
    const Rat z{mid_query.numerator[2].convert_to<Big>()};
    CHECK(abs(z * two_power(-static_cast<int>(a::ExactQueryCoordinates::denominator_shift)) -
              Rat{9, 32}) < rational_word(1e-21));
    CHECK(abs(sum_words(mid_rows[0].velocity[2]) - Rat{9, 16}) < rational_word(1e-21));
    const auto& mid_frame = std::get<cf::ExplicitCollocationTangents>(mid_rows[0].tangents);
    CHECK(abs(sum_words(mid_frame[0][0]) - 1) < rational_word(1e-20));
    CHECK(abs(sum_words(mid_frame[0][2]) - 1) < rational_word(1e-20));
    CHECK(abs(sum_words(mid_frame[1][1]) - 1) < rational_word(1e-20));
    CHECK(abs(sum_words(mid_frame[1][2])) < rational_word(1e-20));

    const cf::DiscreteNonpenetrationInput mid_input{mid_contact,
                                                    mid_rows,
                                                    {},
                                                    {core::AcceptedStateId{2}, {}},
                                                    state.stamp,
                                                    stamp,
                                                    p::RequestedSide::Automatic};
    auto mid_flow =
        take(cf::project_case_uniform_flow(evaluate(provider, provider, flow_workspace, time),
                                           {format, budgets, budgets, {4096, budget}}));
    const auto residual = take(cf::evaluate_discrete_nonpenetration(
        mid_input, std::move(mid_flow), load_limits, shared, *local, *kernel, memory));
    REQUIRE(residual.size() == 1);
    // The typed case supplies 2 ft/min, so the physical wind is 127/12500 m/s.
    const double source_mps = 127. / 12500.;
    residual_expected(take(residual.row(0)).residual,
                      rational_word((source_mps - 9. / 16.) / std::sqrt(2.)), rational_word(1e-15));
    // The same provider remains usable for exact planar points, without
    // claiming the sampled PG point has a facet interpolation stencil.
    auto planar_high = q_high;
    auto planar_low = q_low;
    for (std::size_t i = 2; i < planar_high.size(); i += 3) {
        planar_high[i] = 0.;
        planar_low[i] = 0.;
    }
    auto planar_state = state;
    planar_state.position = {format, planar_high, planar_low};
    const auto planar = take(run(map, planar_state));
    const std::array planar_rows{take(planar.collocation(0))};
    CHECK(std::holds_alternative<a::FreePoint>(
        test_geometry::take(planar_rows[0].query.definition()).definition));
    const cf::DiscreteNonpenetrationInput planar_input{
        take(planar.contact()).get(),   planar_rows, {},
        {core::AcceptedStateId{2}, {}}, state.stamp, stamp,
        p::RequestedSide::Automatic};
    auto planar_flow =
        take(cf::project_case_uniform_flow(evaluate(provider, provider, flow_workspace, time),
                                           {format, budgets, budgets, {4096, budget}}));
    const auto planar_residual = take(cf::evaluate_discrete_nonpenetration(
        planar_input, std::move(planar_flow), load_limits, shared, *local, *kernel, memory));
    residual_expected(take(planar_residual.row(0)).residual, Rat{127, 12500} - Rat{9, 16},
                      rational_word(1e-15));
}

TEST_CASE("Fixed case BDF2 circulation feeds loads and the live-history dual pairing",
          "[AE-02][UW-D06][J1-D01][FD-07][N4-06][verification][negative]") {
    constexpr auto format = core::Precision::TwoWord;
    auto model_case = source(R"({"kind":"constant","value":[0,0,1]})",
                             R"([{"start":0,"end":0.3,"step":0.1,"restart":"initial"}])");
    model_case.field("units").set("length", string("m"));
    model_case.field("units").set("time", string("s"));
    model_case.field("analysis")
        .field("wake")
        .set("convection",
             parse(R"({"kind":"prescribed","velocity":{"kind":"constant","value":[1,0,0]}})"));
    core::ResourcePlan memory{256ULL * 1024 * 1024};
    a::UniformFlowWorkspace flow_work;
    const auto si = materialize(model_case, io::SourceId{1190});
    const cf::CaseFlowClockLimits clock_limits{limits, 1, 3};
    const auto clock =
        take(cf::prepare_case_clocked_flows(si, stamp, flow_work, clock_limits, memory));
    const auto provider =
        take(cf::prepare_case_uniform_flow(si, stamp, {}, flow_work, limits, memory));
    const Number one{1}, zero{0}, tolerance{1, 1000000000};
    const cf::CaseFlowComponentBudget c{one.view(), tolerance.view(), zero.view()};
    const std::array components{c, c, c};
    const cf::CaseFlowProjectionRequest projection{format, components, components, {4096, budget}};
    const auto side = [](std::size_t i) {
        return i == 0 ? p::RequestedSide::Right : p::RequestedSide::Automatic;
    };
    const auto node = [&](std::size_t i) {
        return take(cf::evaluate_case_clocked_flows(clock, clock, {stamp, 0, i, side(i)},
                                                    projection, projection, flow_work, clock_limits,
                                                    memory));
    };
    const auto instant = [](std::size_t i) {
        return i == 0 ? Number{0} : i == 1 ? Number{1, 10} : Number{1, 5};
    };
    const auto flow_at = [&](std::size_t i) {
        const auto t = instant(i);
        return take(cf::project_case_uniform_flow(
            evaluate(provider, provider, flow_work, t, side(i)), projection));
    };

    ResidualFixture fixture;
    auto base = fixture.input();
    base.points = base.points.first(4);
    base.bindings = base.bindings.first(4);
    const a::SinglePositions xyz{std::span{fixture.points}.first(12)};
    base.reference = xyz;
    base.current = xyz;
    base.edges = base.edges.first(5);
    base.faces = base.faces.first(1);
    base.bands = base.bands.first(1);
    base.quads = base.quads.first(1);
    const a::GlobalLimits cap{UINT64_C(5000000000), ls::EstimateAdmission{65536}};
    auto geometry_work = test_geometry::workspace(memory);
    auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
    const auto reference = test_geometry::take(
        a::prepare_owned_geometry(base, test_geometry::limits(), *geometry_work, memory));
    const auto raw_topology = fixture.topology();
    const a::GlobalTopologyInput roles{
        raw_topology.roles.first(1), raw_topology.edges.first(5), {}};
    const auto topology =
        take_global(a::build_global_topology(reference, roles, cap, *global, memory));
    const auto contact_at = [&](std::size_t i, const a::GlobalTopology& owner) {
        auto positions = base;
        constexpr std::array seconds{0., .1, .2};
        const auto t = instant(i);
        positions.context.time = seconds[i];
        positions.context.current = core::GeometryGeneration{3 + i};
        positions.exact_time = t.view();
        const auto geometry = test_geometry::take(
            a::prepare_owned_geometry(positions, test_geometry::limits(), *geometry_work, memory));
        return take_global(a::check_global_geometry(owner, geometry, cap, *global, memory));
    };
    const auto geometry_at = [&](const a::GlobalContactState& contact) -> const a::OwnedGeometry& {
        return take_global(contact.geometry()).get();
    };
    const auto query_at = [&](const a::GlobalContactState& contact) {
        const auto& g = geometry_at(contact);
        return test_geometry::take(
            a::resolve_query(g, {a::QueryId{880}, a::FreePoint{a::FreeSingle{{.75, .375, 0.}}}},
                             g.context(), test_geometry::limits(), *geometry_work, memory));
    };
    const cf::DiscreteNonpenetrationLimits load_limits{
        {format, 1e-10},
        {format, 1e-10, 1e-10},
        {test_geometry::limits(), {ls::EstimateAdmission{65536}, 1000000000}, UINT64_C(5000000000)},
        100000,
        UINT64_C(24000000000)};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const std::array scales{1.};
    const std::array<ls::UnitDimensions, 1> gamma{{{2, 0, -1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 1> speed{{{1, 0, -1, 0, 0, 0, 0}}};
    const auto request_at = [&](std::size_t i, const a::GlobalContactState& contact) {
        auto setup = circulation_request(
            {scales, scales, gamma, speed, circulation_context(geometry_at(contact).context())}, 1,
            ls::Norm::Infinity);
        setup.rhs_generation = core::RhsGeneration{31 + i};
        setup.retained_generation = core::NumericalGeneration{41 + i};
        return setup;
    };
    const auto zeros = coordinates(0, 0, 0);
    const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
    const auto face = fixture.faces[0].id;
    const auto make_row = [&](const a::ResolvedQuery& query, double normal_velocity) {
        return cf::DiscreteCollocation{fixture.quads[0].facets[0],
                                       query,
                                       frame,
                                       coordinates(0, 0, normal_velocity),
                                       a::SingleWord{82},
                                       a::SingleWord{91}};
    };
    const auto history = [&](std::size_t i, const a::GlobalTopology& owner) {
        const auto contact = contact_at(i, owner);
        const auto query = query_at(contact);
        const std::array rows{make_row(query, i == 0 ? .125 : .375)};
        const cf::DiscreteNonpenetrationInput input{
            contact, rows, {}, {core::AcceptedStateId{2}, {}}, {50, i}, stamp, side(i)};
        return take(cf::solve_discrete_circulation(input, flow_at(i), load_limits,
                                                   request_at(i, contact), shared, *local, *kernel,
                                                   memory));
    };
    const auto earlier = history(0, topology);
    const auto previous = history(1, topology);
    const auto t0 = node(0);
    const auto t1 = node(1);
    auto skipped = node(2);
    reject(cf::FixedBdf2Circulation::prepare(std::move(skipped), t0, t1, previous, earlier, 4096,
                                             10000, memory),
           core::FailureCode::StaleGeneration);
    // NOLINTNEXTLINE(bugprone-use-after-move) Rejected preparation keeps the node owner active.
    CHECK(skipped.ambient().source().jet.active());
    reject(
        cf::FixedBdf2Circulation::prepare(node(2), t1, t0, earlier, previous, 4096, 10000, memory),
        core::FailureCode::StaleGeneration);
    const auto other_topology =
        take_global(a::build_global_topology(reference, roles, cap, *global, memory));
    const auto foreign_history = history(0, other_topology);
    reject(cf::FixedBdf2Circulation::prepare(node(2), t1, t0, previous, foreign_history, 4096,
                                             10000, memory),
           core::FailureCode::StaleGeneration);
    auto altered_case = model_case;
    altered_case.field("analysis").field("intervals").elements[0].set("end", parse("0.4"));
    const auto edited_si = materialize(altered_case, io::SourceId{1190});
    auto edited_limits = clock_limits;
    edited_limits.maximum_steps_per_interval = 4;
    const auto foreign_clock =
        take(cf::prepare_case_clocked_flows(edited_si, stamp, flow_work, edited_limits, memory));
    const auto wrong_case = take(
        cf::evaluate_case_clocked_flows(foreign_clock, foreign_clock, {stamp, 0, 0, side(0)},
                                        projection, projection, flow_work, edited_limits, memory));
    reject(cf::FixedBdf2Circulation::prepare(node(2), t1, wrong_case, previous, earlier, 4096,
                                             10000, memory),
           core::FailureCode::StaleGeneration);
    const auto stencil = take(
        cf::FixedBdf2Circulation::prepare(node(2), t1, t0, previous, earlier, 4096, 10000, memory));

    const auto contact = contact_at(2, topology);
    const auto query = query_at(contact);
    const std::array rows{make_row(query, .75)};
    const cf::DiscreteNonpenetrationInput input{
        contact, rows, {}, {core::AcceptedStateId{2}, {}}, {50, 2}, stamp, side(2)};
    const auto request = request_at(2, contact);
    const std::array midpoints{cf::EdgeMidpointVelocity{fixture.edges[0].id, zeros},
                               cf::EdgeMidpointVelocity{fixture.edges[1].id, zeros},
                               cf::EdgeMidpointVelocity{fixture.edges[2].id, zeros},
                               cf::EdgeMidpointVelocity{fixture.edges[3].id, zeros}};
    cf::DiscreteVectorLoadInput load{
        input.trial_stamp, request.retained_generation, midpoints, {}, a::SingleWord{1.25}, zeros};
    load.fixed_bdf2 = &stencil;
    const auto result =
        take(cf::solve_discrete_circulation_loads(input, load, flow_at(2), load_limits, request,
                                                  *geometry_work, shared, *local, *kernel, memory));
    const Rat bdf2 =
        (3 * solution_word(result.primal().solution(), 0) -
         4 * solution_word(previous.solution(), 0) + solution_word(earlier.solution(), 0)) *
        5;
    REQUIRE(abs(bdf2) > rational_word(.01));
    REQUIRE(result.bdf2_rate_count() == 1);
    const auto rate = take(result.bdf2_rate(0));
    CHECK(rate.face == face);
    residual_expected(rate.rate, bdf2, rational_word(1e-20));
    const Rat half_density{5, 8};
    for (std::size_t i = 0; i < 2; ++i)
        residual_expected(take(result.facet(i)).force[2], bdf2 * half_density,
                          rational_word(1e-20));

    const std::array displacements{zeros, zeros, zeros, zeros};
    const auto unit_z = coordinates(0, 0, 1);
    const std::array collocation_direction{
        cf::CirculationCollocationDirection{face, zeros, {zeros, zeros}, unit_z}};
    const cf::DiscreteCirculationDirection geometry_direction{base.points, displacements,
                                                              collocation_direction, zeros};
    const cf::DiscreteImplicitVectorLoadDirection variation{
        load.stamp, load.circulation, geometry_direction, midpoints, {}};
    auto tangent = request;
    tangent.rhs_generation = core::RhsGeneration{52};
    tangent.retained_generation = core::NumericalGeneration{62};
    const auto derivative = take(cf::solve_discrete_circulation_implicit_load_direction(
        input, load, variation, flow_at(2), load_limits, request, tangent, *geometry_work, shared,
        *local, *kernel, memory));
    const auto dgamma = solution_word(derivative.tangent().solution(), 0);
    REQUIRE(abs(dgamma) > rational_word(.01));
    for (std::size_t i = 0; i < 2; ++i)
        residual_expected(take(derivative.direction().facet(i)).force[2],
                          15 * dgamma * half_density, rational_word(1e-20));
    const auto sample = cf::DiscreteLoadSampleWeights{zeros, zeros};
    const std::array edge_seeds{sample, sample, sample, sample};
    const std::array facet_seeds{sample, sample};
    const cf::DiscreteVectorLoadWeights seeds{
        load.stamp, load.circulation, edge_seeds, facet_seeds, displacements, unit_z, zeros};
    const std::array<ls::UnitDimensions, 1> dual_gamma{{{-2, 0, 1, 0, 0, 0, 0}}};
    const std::array<ls::UnitDimensions, 1> dual_speed{{{-1, 0, 1, 0, 0, 0, 0}}};
    auto dual_context = request.scales.context;
    dual_context.cache.matrix = core::NumericalGeneration{70};
    dual_context.cache.scales = core::ScaleGeneration{71};
    dual_context.unknown_layout = core::ScopeId{72};
    dual_context.residual_layout = core::ScopeId{73};
    auto transpose = circulation_request({scales, scales, dual_speed, dual_gamma, dual_context}, 1,
                                         ls::Norm::Infinity);
    transpose.rhs_generation = core::RhsGeneration{74};
    transpose.retained_generation = core::NumericalGeneration{75};
    const auto dual = take(cf::solve_discrete_circulation_implicit_load_pullback(
        input, load, seeds, flow_at(2), load_limits, request, transpose, *geometry_work, shared,
        *local, *kernel, memory));
    const auto row_dual = take(dual.collocation(0));
    WorkEnclosure left, right;
    left.add(derivative.direction().total_force_direction()[2], unit_z[2]);
    right.add(row_dual.velocity[2], unit_z[2]);
    CHECK(abs(left.value) > rational_word(.01));
    CHECK(abs(left.value - right.value) <= left.radius + right.radius + rational_word(1e-15));
}

namespace {
struct WakeCostMesh {
    explicit WakeCostMesh(std::size_t wake) {
        const auto panels = wake + 1;
        const auto point = [](std::size_t line, std::size_t side) {
            return a::PointId{2 * line + side + 1};
        };
        const auto edge = [&](a::PointId from, a::PointId to) {
            const a::EdgeId id{1000 + edges.size()};
            edges.push_back({id, {from, to}});
            edge_owners.emplace_back(id, a::CanonicalEdge{});
            return id;
        };
        std::vector<a::EdgeId> transverse;
        for (std::size_t row = 0; row <= panels; ++row) {
            for (std::size_t side = 0; side < 2; ++side) {
                ids.push_back(point(row, side));
                bindings.push_back(test_geometry::binding());
                xyz.insert(xyz.end(), {static_cast<double>(row), static_cast<double>(side), 0.});
            }
            transverse.push_back(edge(point(row, 1), point(row, 0)));
        }
        trailing.push_back({transverse[1], {point(1, 0), point(1, 1)}});
        for (std::size_t row = 0; row < panels; ++row) {
            const std::array vertices{point(row, 0), point(row, 1), point(row + 1, 1),
                                      point(row + 1, 0)};
            const a::PhysicalFaceId face{100 + row};
            const a::BandId band{50 + row};
            bands.push_back(band);
            faces.push_back({face, band});
            roles.push_back({face, row == 0 ? a::AeroFaceRole::Bound : a::AeroFaceRole::Wake});
            quads.push_back(a::QuadInput{face,
                                         vertices,
                                         {transverse[row], edge(vertices[1], vertices[2]),
                                          transverse[row + 1], edge(vertices[3], vertices[0])},
                                         edge(vertices[0], vertices[2]),
                                         {a::FacetId{200 + 2 * row}, a::FacetId{201 + 2 * row}}});
        }
    }

    [[nodiscard]] a::GeometryInput input(std::span<const double> current, const Number& time,
                                         double time_word) const {
        auto context = test_geometry::context();
        context.time = time_word;
        return {ids,
                bindings,
                a::SinglePositions{xyz},
                a::SinglePositions{current},
                edges,
                faces,
                bands,
                {},
                quads,
                a::SingleWord{.125},
                context,
                time.view()};
    }

    std::vector<a::PointId> ids;
    std::vector<a::PointBinding> bindings;
    std::vector<double> xyz;
    std::vector<a::EdgeInput> edges;
    std::vector<a::EdgeBindingInput> edge_owners;
    std::vector<a::BandId> bands;
    std::vector<a::FaceInput> faces;
    std::vector<a::FaceRoleInput> roles;
    std::vector<a::QuadInput> quads;
    std::vector<a::TrailingEdgeInput> trailing;
};

struct WakeCostWords {
    void add(const a::ScaledEnclosure& x) {
        const auto word = [&](double value) {
            words.push_back(std::bit_cast<std::uint64_t>(value));
        };
        const auto center = x.center();
        words.push_back(center.index());
        if (const auto* single = std::get_if<a::SingleWord>(&center))
            word(single->value);
        else {
            const auto& pair = std::get<a::DoubleWord>(center);
            word(pair.high);
            word(pair.low);
        }
        words.push_back(std::bit_cast<std::uint64_t>(static_cast<std::int64_t>(x.exponent())));
        word(x.radius());
    }
    void add(const cf::DiscreteLoadVector& values) {
        for (const auto& value : values)
            add(value);
    }
    void add(const std::optional<cf::DiscreteLoadVector>& value) {
        words.push_back(value.has_value());
        if (value)
            add(*value);
    }

    void compare_reference(const cf::DiscreteFreeWakeVelocity& result,
                           const std::array<double, 2>& expected) const {
        for (std::size_t end = 0; end < expected.size(); ++end) {
            const auto row = take(result.row(end == 0 ? 0 : result.size() - 1));
            const auto& vertical = row.velocity[2];
            const auto quantum = two_power(vertical.exponent());
            const Rat center = sum_words(vertical.center()) * quantum;
            const Rat radius = rational_word(vertical.radius()) * quantum;
            CHECK(abs(center - rational_word(expected[end])) <= radius + rational_word(1e-12));
        }
    }

    void consume(const cf::DiscreteFreeWakeVelocity& result) {
        for (std::size_t i = 0; i < result.size(); ++i) {
            const auto row = take(result.row(i));
            words.push_back(row.owner.value());
            add(row.velocity);
            add(row.direction);
            source_facets += row.source_facets;
            far_facets += row.far_facets;
            cancelled_pairs += row.cancelled_analytic_pairs;
        }
        for (std::size_t i = 0; i < result.point_count(); ++i) {
            const auto point = take(result.point_dual(i));
            words.push_back(point.point.value());
            add(point.position);
        }
        for (std::size_t i = 0; i < result.strength_count(); ++i) {
            const auto strength = take(result.strength_dual(i));
            add(strength.value);
        }
        if (result.strength_count() != 0) {
            const auto field = take(result.field_dual());
            add(field.uniform_flow);
            add(field.flow_time);
            add(field.core_radius);
        }
        steps = result.used_steps();
    }

    void consume(const cf::WakeTransportResult& result, bool dual) {
        consume(take(result.free_velocity()).get());
        for (std::size_t i = 0; i < result.size(); ++i) {
            const auto row = take(result.row(i));
            words.push_back(row.release.point.value());
            words.push_back(row.release.epoch.value());
            add(row.residual);
            add(row.direction);
            if (dual) {
                const auto transpose = take(result.position_dual(i));
                for (const auto& term : transpose.stencil)
                    add(term);
            }
        }
    }

    [[nodiscard]] std::uint64_t fingerprint() const {
        std::uint64_t value = 14695981039346656037ULL;
        for (const auto word : words) {
            value ^= word;
            value *= 1099511628211ULL;
        }
        return value;
    }
    std::uint64_t source_facets{}, far_facets{}, cancelled_pairs{}, steps{};
    std::vector<std::uint64_t> words;
};

void observe_wake_cost(std::size_t rows) {
    using Clock = std::chrono::steady_clock;
    constexpr auto operation = core::Operation::EvaluateAero;
    const auto elapsed = [](Clock::time_point start) {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
    };
    const auto prepare_begin = Clock::now();
    const WakeCostMesh mesh{rows};
    const Number now{5, 4}, then{9, 8}, older{1};
    const a::GlobalLimits contact_caps{std::numeric_limits<std::uint64_t>::max() / 4,
                                       numerics::linear::EstimateAdmission{4096}};
    auto geometry_caps = test_geometry::limits();
    geometry_caps.max_steps = 1000000;
    geometry_caps.predicate.max_steps = 10000;
    constexpr std::uint64_t bytes = 512ULL * 1024 * 1024;
    core::ResourcePlan preparation{bytes};
    auto geometry_work = test_geometry::workspace(preparation);
    auto global = take(a::GlobalWorkspace::create(preparation, contact_caps.metadata));
    const auto geometry = test_geometry::take(a::prepare_owned_geometry(
        mesh.input(mesh.xyz, now, 1.25), geometry_caps, *geometry_work, preparation));
    const auto topology = take_global(
        a::build_global_topology(geometry, {mesh.roles, mesh.edge_owners, mesh.trailing},
                                 contact_caps, *global, preparation));
    REQUIRE(topology.trailing_edges().size() == 1);
    const auto contact = take_global(
        a::check_global_geometry(topology, geometry, contact_caps, *global, preparation));
    const auto history = [&](const Number& time, double elapsed_seconds, double shift) {
        auto xyz = mesh.xyz;
        for (std::size_t line = 2; line < mesh.ids.size() / 2; ++line) {
            const auto dx = static_cast<double>(line - 1) * shift;
            xyz[6 * line] -= dx;
            xyz[6 * line + 3] -= dx;
        }
        const auto sample = test_geometry::take(a::prepare_owned_geometry(
            mesh.input(xyz, time, elapsed_seconds), geometry_caps, *geometry_work, preparation));
        return take_global(
            a::check_global_geometry(topology, sample, contact_caps, *global, preparation));
    };
    const auto previous_contact = history(then, 1.125, .125);
    const auto older_contact = history(older, 1., .25);
    const auto geometry_ns = elapsed(prepare_begin);
    const auto geometry_payload = preparation.allocated_bytes();
    const auto prepared_memory = take(core::sample_process_memory(operation));
    const auto case_begin = Clock::now();
    auto recipe = source(R"({"kind":"constant","value":[1,0,0.125]})",
                         R"([{"start":1,"end":1.5,"step":0.125,"restart":"initial"}])");
    recipe.field("units").set("length", string("m"));
    recipe.field("units").set("time", string("s"));
    recipe.field("analysis").field("wake").set("convection", parse(R"({"kind":"free"})"));
    const auto si = materialize(recipe, io::SourceId{490});
    const cf::CaseFlowClockLimits clock_limits{limits, 1, 4};
    const auto clock =
        take(cf::prepare_case_reference_clock(si, stamp, clock_limits, 1000000, preparation));
    const auto destination = cf::CaseFlowClockQuery{stamp, 0, 2, p::RequestedSide::Automatic};
    const cf::WakeTransportClock transport_clock{clock, clock, destination};
    a::UniformFlowWorkspace uniform;
    const auto provider =
        take(cf::prepare_case_known_flow(si, stamp, {}, uniform, gust_limits, preparation));
    const auto q = test_geometry::take(
        a::resolve_query(geometry, {a::QueryId{42}, a::FreePoint{a::FreeSingle{{.75, .5, 0}}}},
                         geometry.context(), geometry_caps, *geometry_work, preparation));
    const auto zero = coordinates(0, 0, 0);
    const std::array tangents{coordinates(1, 0, 0), coordinates(0, 1, 0)};
    const std::array bound{cf::DiscreteCollocation{mesh.quads.front().facets[0], q, tangents, zero,
                                                   a::SingleWord{.5}, a::SingleWord{.0625}}};
    std::vector<cf::FrozenWakeTrace> traces, trace_direction;
    std::vector<cf::DiscreteWakePanel> wake;
    for (std::size_t i = 0; i < rows; ++i) {
        const a::P1StrengthId right{1000 + i};
        traces.push_back({right, a::DoubleWord{-.25 + static_cast<double>(i) / 64, 0x1p-58}});
        trace_direction.push_back({right, a::SingleWord{.03125}});
        wake.emplace_back(
            mesh.quads[i + 1].facets[0],
            i == 0 ? std::variant<cf::AttachedWakeFront, a::P1StrengthId>{cf::AttachedWakeFront{}}
                   : std::variant<cf::AttachedWakeFront, a::P1StrengthId>{traces[i - 1].id},
            right);
    }
    const cf::DiscreteNonpenetrationInput input{
        contact, bound, wake, {core::AcceptedStateId{2}, traces}, {7, 11}, stamp, destination.side};
    std::vector<cf::WakeMaterialPoint> material;
    std::vector<a::PointId> particles;
    for (std::size_t side = 0; side < 2; ++side) {
        for (std::size_t i = 0; i < rows; ++i) {
            const auto release =
                take_global(topology.release_key(mesh.ids[2 + side], a::ReleaseEpoch{8 + i}));
            const auto point = mesh.ids[2 * (i + 2) + side];
            material.push_back({release, point, wake[i].facet});
            particles.push_back(point);
        }
    }
    const std::array snapshots{cf::WakeMaterialSnapshot{contact, material},
                               cf::WakeMaterialSnapshot{previous_contact, material},
                               cf::WakeMaterialSnapshot{older_contact, material}};
    const std::vector variations(mesh.ids.size(), coordinates(.125, 0, 0));
    const std::vector previous_direction(particles.size(), coordinates(.25, 0, 0));
    const std::vector older_direction(particles.size(), coordinates(.375, 0, 0));
    const std::vector seeds(particles.size(), coordinates(.125, -.25, .375));
    const cf::DiscreteFreeWakeDirection direction{mesh.ids,           variations,
                                                  trace_direction,    coordinates(.125, 0, 0),
                                                  a::SingleWord{.25}, a::SingleWord{.03125}};
    const cf::WakeTransportHistoryDirection history_direction{previous_direction, older_direction};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(preparation));
    auto kernel = take(a::P1Workspace::create(preparation));
    const auto case_ns = elapsed(case_begin);
    const auto provider_memory = take(core::sample_process_memory(operation));
    std::cout << "WAKE_PREP," << rows << ',' << geometry_ns << ',' << case_ns << ','
              << geometry_payload << ',' << preparation.allocated_bytes() << ','
              << prepared_memory.resident_bytes << ',' << provider_memory.resident_bytes << ','
              << provider_memory.peak_resident_bytes << '\n';

    struct Role {
        std::string_view name;
        bool transport;
        bool tangent;
        bool adjoint;
    };
    constexpr std::array roles{
        Role{"velocity_value", false, false, false}, Role{"velocity_jvp", false, true, false},
        Role{"velocity_vjp", false, false, true}, Role{"transport_value", true, false, false},
        Role{"transport_full", true, true, true}};
    constexpr std::array end_velocity{
        std::array{.247542502068428649, .152597465749261047},
        std::array{.241277910608254584, .140543900713236010},
        std::array{.240636854147210874, .134174028149262049},
    };
    const auto& expected = end_velocity[rows == 2 ? 0 : rows == 4 ? 1 : 2];
    const Number enclosure{1, 1000000000};
    for (const auto precision : {core::Precision::Binary64, core::Precision::TwoWord}) {
        const cf::DiscreteKnownFlowInput known{
            provider,    provider,      {stamp, now.view()}, gust_projection(precision, enclosure),
            gust_limits, geometry_caps, *geometry_work,      uniform};
        const cf::DiscreteFreeWakeLimits velocity_limits{
            {{precision, 1e-10},
             {precision, 1e-10, 1e-10},
             {geometry_caps,
              {numerics::linear::EstimateAdmission{65536}, 1000000000},
              UINT64_C(24000000000)},
             10000000,
             UINT64_C(1000000000000)},
            particles.size()};
        const cf::WakeTransportLimits transport_limits{clock_limits, particles.size(), 1000000,
                                                       4096, precision};
        for (const auto role : roles) {
            std::vector<std::uint64_t> first;
            for (std::uint64_t iteration = 0; iteration < 3; ++iteration) {
                const auto before = take(core::sample_process_memory(operation));
                core::ResourcePlan plan{bytes};
                const auto* tangent = role.tangent ? &direction : nullptr;
                const auto dual = role.adjoint ? std::span<const cf::NonpenetrationVector>{seeds}
                                               : std::span<const cf::NonpenetrationVector>{};
                const auto start = Clock::now();
                WakeCostWords observed;
                std::int64_t ns{};
                core::ProcessMemory after{};
                if (role.transport) {
                    const auto value = take(cf::evaluate_fixed_free_wake_transport(
                        transport_clock, snapshots, input, known, particles, tangent,
                        role.tangent ? history_direction : cf::WakeTransportHistoryDirection{},
                        dual, transport_limits, velocity_limits, shared, *local, *kernel, plan));
                    ns = elapsed(start);
                    after = take(core::sample_process_memory(operation));
                    observed.compare_reference(take(value.free_velocity()).get(), expected);
                    observed.consume(value, role.adjoint);
                } else {
                    const auto value = take(cf::evaluate_discrete_free_wake_velocity(
                        input, known, particles, tangent, dual, velocity_limits, shared, *local,
                        *kernel, plan));
                    ns = elapsed(start);
                    after = take(core::sample_process_memory(operation));
                    observed.compare_reference(value, expected);
                    observed.consume(value);
                }
                REQUIRE(observed.source_facets == 2 * (rows + 1) * particles.size());
                if (iteration == 0)
                    first = observed.words;
                else
                    CHECK(first == observed.words);
                std::cout << "WAKE_COST," << rows << ','
                          << (precision == core::Precision::Binary64 ? "binary64" : "two_word")
                          << ',' << role.name << ',' << iteration << ',' << particles.size() << ','
                          << observed.source_facets << ',' << observed.far_facets << ','
                          << observed.cancelled_pairs << ',' << observed.steps << ',' << ns << ','
                          << plan.allocated_bytes() << ',' << before.resident_bytes << ','
                          << after.resident_bytes << ',' << after.peak_resident_bytes << ','
                          << observed.fingerprint() << '\n'
                          << std::flush;
            }
        }
    }
}
} // namespace

TEST_CASE("Wake cost with two retained panels", "[.wake-cost-2]") {
    observe_wake_cost(2);
}
TEST_CASE("Wake cost with four retained panels", "[.wake-cost-4]") {
    observe_wake_cost(4);
}
TEST_CASE("Wake cost with eight retained panels", "[.wake-cost-8]") {
    observe_wake_cost(8);
}

namespace {
struct FiniteWingMesh {
    explicit FiniteWingMesh(const nlohmann::json& coordinates) {
        const auto point = [](std::size_t chord, std::size_t span) {
            return a::PointId{1 + 5 * chord + span};
        };
        for (std::size_t i = 0; i < coordinates.size(); ++i) {
            ids.emplace_back(i + 1);
            bindings.push_back(test_geometry::binding());
            const auto x = coordinates[i].get<std::array<double, 3>>();
            xyz.insert(xyz.end(), x.begin(), x.end());
        }
        std::map<std::array<a::PointId, 2>, a::EdgeId> by_ends;
        const auto edge = [&](a::PointId from, a::PointId to) {
            const std::array key{std::min(from, to), std::max(from, to)};
            const a::EdgeId new_id{1000 + edges.size()};
            const auto [found, inserted] = by_ends.try_emplace(key, new_id);
            if (inserted) {
                edges.push_back({new_id, {from, to}});
                owners.emplace_back(new_id, a::CanonicalEdge{});
            }
            return found->second;
        };
        for (std::size_t chord = 0; chord < 3; ++chord) {
            for (std::size_t span = 0; span < 4; ++span) {
                const auto index = 4 * chord + span;
                const a::PhysicalFaceId face{100 + index};
                const a::BandId band{50 + span};
                const std::array vertex{point(chord, span), point(chord, span + 1),
                                        point(chord + 1, span + 1), point(chord + 1, span)};
                const std::array outline{edge(vertex[0], vertex[1]), edge(vertex[1], vertex[2]),
                                         edge(vertex[2], vertex[3]), edge(vertex[3], vertex[0])};
                faces.push_back({face, band});
                roles.push_back(
                    {face, chord == 2 ? a::AeroFaceRole::Wake : a::AeroFaceRole::Bound});
                quads.emplace_back(
                    face, vertex, outline, edge(vertex[0], vertex[2]),
                    std::array{a::FacetId{200 + 2 * index}, a::FacetId{201 + 2 * index}});
                if (chord == 2)
                    trailing.push_back({outline[0], {vertex[0], vertex[1]}});
            }
        }
    }

    std::vector<a::PointId> ids;
    std::vector<a::PointBinding> bindings;
    std::vector<double> xyz;
    std::vector<a::EdgeInput> edges;
    std::vector<a::EdgeBindingInput> owners;
    std::vector<a::FaceInput> faces;
    std::vector<a::FaceRoleInput> roles;
    std::vector<a::QuadInput> quads;
    std::vector<a::TrailingEdgeInput> trailing;
};

JsonToken finite_wing_case(const nlohmann::json& fixture, const nlohmann::json& wing) {
    auto recipe = root("aerodynamic_steady");
    recipe.field("model").set("geometry", parse(wing.at("case_geometry").dump()));
    recipe.field("model").set("aero_surfaces", parse(wing.at("case_surfaces").dump()));
    recipe.field("flow").set("density", parse("1.2"));
    recipe.field("flow").set(
        "reference",
        parse(
            R"({"density":1.2,"speed":20,"area":8.4,"chord":1.114285714285714,"span":8,"moment_origin":[0.375,0,0]})"));
    recipe.field("flow").field("uniform").set("value", parse(wing.at("flow_m_s").dump()));
    recipe.field("analysis").field("speed").set("value", parse("20"));
    recipe.field("observations").elements[0].set("surfaces", parse(R"(["1","2"])"));
    auto& wake = recipe.field("analysis").field("wake");
    wake.set("core_radius", parse(fixture.at("core_radius_m").dump()));
    wake.set("convection",
             parse(R"({"kind":"prescribed","velocity":{"kind":"constant","value":[20,0,0]}})"));
    return recipe;
}

bool entire_wing_enclosure_fits(const Rat& center, const Rat& radius, const Rat& oracle,
                                const Rat& tolerance) {
    return abs(center - oracle) + radius <= tolerance;
}

void check_frozen_wing_enclosure(const a::ScaledEnclosure& value, const Rat& oracle,
                                 const Rat& tolerance) {
    const auto scale = two_power(value.exponent());
    const Rat center = sum_words(value.center()) * scale;
    const Rat radius = rational_word(value.radius()) * scale;
    CHECK(entire_wing_enclosure_fits(center, radius, oracle, tolerance));
}

void verify_finite_wing(const nlohmann::json& fixture, const nlohmann::json& wing) {
    const auto recipe = finite_wing_case(fixture, wing);
    const auto si = materialize(recipe, io::SourceId{532});
    const auto typed = take(take(si.source_case()).get().root());
    const auto analysis = std::get<io::schema::AnalysisStatic>(typed.analysis.record());
    const auto wake_case = std::get<io::schema::WakeActive>(analysis.wake.record());
    const auto radius = take(si.quantity(wake_case.core_radius));
    CHECK(std::holds_alternative<io::schema::ConvectionPrescribed>(wake_case.convection.record()));
    const auto model = typed.model.record();
    REQUIRE(model.aero_surfaces.size() == 2);
    for (core::Index half = 0; half < 2; ++half) {
        const auto surface = take(model.aero_surfaces.at(half)).record();
        CHECK(std::holds_alternative<io::schema::AeroAttachmentFixed>(surface.attachment.record()));
        CHECK(surface.chord.size() == 3);
        CHECK(surface.span.size() == 3);
        CHECK(surface.sharing.size() == 1);
    }
    const double core_radius = std::visit(
        [](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, io::ScaledBinary64>)
                return std::scalbn(value.center, static_cast<int>(value.exponent2));
            else
                return std::scalbn(value.high, static_cast<int>(value.exponent2));
        },
        radius.enclosure);
    REQUIRE(core_radius == fixture.at("core_radius_m").get<double>());

    const FiniteWingMesh mesh{wing.at("vertices_m")};
    const std::array bands{a::BandId{50}, a::BandId{51}, a::BandId{52}, a::BandId{53}};
    auto geometry_caps = test_geometry::limits();
    geometry_caps.max_steps = 1000000;
    geometry_caps.predicate.max_steps = 10000;
    constexpr std::uint64_t bytes = 256ULL * 1024 * 1024;
    core::ResourcePlan memory{bytes};
    auto geometry_work = test_geometry::workspace(memory);
    const Number instant{5, 4};
    const a::SinglePositions xyz{mesh.xyz};
    const auto geometry =
        test_geometry::take(a::prepare_owned_geometry({mesh.ids,
                                                       mesh.bindings,
                                                       xyz,
                                                       xyz,
                                                       mesh.edges,
                                                       mesh.faces,
                                                       bands,
                                                       {},
                                                       mesh.quads,
                                                       a::SingleWord{core_radius},
                                                       test_geometry::context(),
                                                       instant.view()},
                                                      geometry_caps, *geometry_work, memory));
    const a::GlobalLimits contact_caps{std::numeric_limits<std::uint64_t>::max() / 4,
                                       numerics::linear::EstimateAdmission{65536}};
    auto global = take(a::GlobalWorkspace::create(memory, contact_caps.metadata));
    const auto topology = take_global(a::build_global_topology(
        geometry, {mesh.roles, mesh.owners, mesh.trailing}, contact_caps, *global, memory));
    REQUIRE(topology.trailing_edges().size() == 4);
    const auto contact =
        take_global(a::check_global_geometry(topology, geometry, contact_caps, *global, memory));

    std::vector<a::ResolvedQuery> queries;
    constexpr std::size_t bound_count = 8;
    queries.reserve(bound_count);
    for (std::size_t i = 0; i < bound_count; ++i) {
        const auto selected = wing.at("collocations_m").at(i).get<std::array<double, 3>>();
        queries.push_back(test_geometry::take(
            a::resolve_query(geometry, {a::QueryId{i + 1}, a::FreePoint{a::FreeSingle{selected}}},
                             geometry.context(), geometry_caps, *geometry_work, memory)));
    }
    const auto gamma = wing.at("circulation_m2_s").get<std::array<double, bound_count>>();
    const auto reference = wing.at("first_bound_column").get<std::array<double, bound_count>>();
    std::vector<cf::DiscreteCollocation> rows;
    for (std::size_t i = 0; i < bound_count; ++i) {
        const auto axis = wing.at("tangents_y").at(i).get<std::array<double, 3>>();
        const std::array tangent{coordinates(1, 0, 0), coordinates(axis[0], axis[1], axis[2])};
        rows.emplace_back(mesh.quads[i].facets[0], queries[i], tangent, coordinates(0, 0, 0),
                          a::SingleWord{gamma[i]}, a::SingleWord{i == 0 ? 1. : 0.});
    }
    std::vector<cf::FrozenWakeTrace> traces;
    std::vector<cf::DiscreteWakePanel> wake;
    for (std::size_t span = 0; span < 4; ++span) {
        const a::P1StrengthId history{700 + span};
        traces.push_back({history, a::SingleWord{gamma[4 + span]}});
        wake.emplace_back(mesh.quads[8 + span].facets[0], cf::AttachedWakeFront{}, history);
    }
    const cf::DiscreteNonpenetrationInput input{contact,
                                                rows,
                                                wake,
                                                {core::AcceptedStateId{2}, traces},
                                                {815, 6},
                                                stamp,
                                                p::RequestedSide::Automatic};
    const cf::DiscreteNonpenetrationLimits caps{
        {core::Precision::TwoWord, 1e-12},
        {core::Precision::TwoWord, 1e-12, 1e-12},
        {geometry_caps,
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(32000000000)},
        1000000,
        UINT64_C(1000000000000)};
    a::UniformFlowWorkspace flow_work;
    const auto flow = take(cf::prepare_case_uniform_flow(si, stamp, {}, flow_work, limits, memory));
    const auto prescribed =
        take(cf::prepare_case_prescribed_wake(si, stamp, {}, flow_work, limits, memory));
    const auto transport = take(cf::evaluate_case_prescribed_wake(
        prescribed, prescribed, {stamp, instant.view()}, flow_work, limits.flow, memory));
    const Number speed{20};
    exact(take(transport.jet().velocity(0)), speed);
    exact(take(transport.jet().velocity(1)), Number{0});
    exact(take(transport.jet().velocity(2)), Number{0});
    auto sampled = evaluate(flow, flow, flow_work, instant);
    const Number one{1}, absolute{1, 1000000000}, zero{0};
    const cf::CaseFlowComponentBudget component{one.view(), absolute.view(), zero.view()};
    const std::array bounds{component, component, component};
    auto projected = take(cf::project_case_uniform_flow(
        std::move(sampled), {core::Precision::TwoWord, bounds, bounds, {4096, budget}}));
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const auto result = take(cf::evaluate_discrete_nonpenetration(input, std::move(projected), caps,
                                                                  shared, *local, *kernel, memory));
    REQUIRE(result.size() == bound_count);
    const auto tolerance = rational_word(2e-10);
    for (std::size_t i = 0; i < bound_count; ++i) {
        CAPTURE(i);
        const auto row = take(result.row(i));
        CHECK(row.face == mesh.faces[i].id);
        check_frozen_wing_enclosure(row.residual, Rat{0}, tolerance);
        check_frozen_wing_enclosure(row.circulation_action, rational_word(reference[i]), tolerance);
    }
}
} // namespace

TEST_CASE("Frozen finite wing closes all discrete residual rows and a bound direction",
          "[V-AERO][AE-02][UW-D03][N4-06][verification]") {
    std::ifstream input{UVLAS_FINITE_WING_CASE_FIXTURE};
    REQUIRE(input);
    const auto fixture = nlohmann::json::parse(input);
    CHECK(fixture.at("wake_length_m").get<int>() == 32);
    CHECK(fixture.at("chord_panels").get<int>() == 2);
    CHECK(fixture.at("span_panels").get<int>() == 4);
    REQUIRE(fixture.at("cases").size() == 2);
    for (const auto& wing : fixture.at("cases")) {
        const auto angle = wing.at("dihedral_degrees").get<int>();
        CAPTURE(angle);
        verify_finite_wing(fixture, wing);
    }
}

TEST_CASE("A frozen-wing reference window bounds both enclosure endpoints",
          "[V-AERO][AE-02][UW-D03][N4-06][verification][negative]") {
    const Rat oracle{5, 4};
    const Rat tolerance{1, 8};
    const Rat radius{1, 16};
    for (const int sign : {-1, 1}) {
        CAPTURE(sign);
        const Rat center = oracle + Rat{3 * sign, 32};
        REQUIRE(abs(center - oracle) < tolerance);
        REQUIRE(radius < tolerance);
        CHECK_FALSE(entire_wing_enclosure_fits(center, radius, oracle, tolerance));
        CHECK(entire_wing_enclosure_fits(center, Rat{1, 32}, oracle, tolerance));
    }
}
