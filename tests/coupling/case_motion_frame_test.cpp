#include "coupling/case_motion_frame.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <thread>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
namespace cf = coupling::detail;
namespace m = cf::motion::m;
using Failure = core::FailureCode;
using Fraction = boost::multiprecision::cpp_rational;
using Integer = boost::multiprecision::cpp_int;

constexpr std::uint64_t budget = std::uint64_t{1} << 40U;
constexpr std::string_view origin = "[10000,20000,30000]";
constexpr std::string_view rotation = "[[0.6,-0.8,0],[0.8,0.6,0],[0,0,1]]";

structure::ExactLimits exact_limits() {
    return {10000,
            1000,
            1000,
            64,
            64,
            32768,
            65536,
            structure::ExactEstimateAdmission{budget},
            structure::ExactTimingMode::Disabled};
}

cf::CaseMotionFrameLimits frame_limits() {
    return {exact_limits(), 5000, cf::case_motion_frame_interval_scratch_bytes};
}

template <class T>
T take(core::Result<T> result) {
    if (!result)
        FAIL("Unexpected motion frame failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}

template <class T>
void rejects(const core::Result<T>& result, Failure code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
}

template <class T, class Diagnostic>
T checked(core::Result<T> result, const Diagnostic& diagnostic) {
    if (!result) {
        const auto path = diagnostic.path();
        const auto detail = diagnostic.detail();
        if (path && detail)
            FAIL("Motion fixture rejected at "
                 << *path << ", code " << static_cast<int>(result.error().code) << ", reason "
                 << static_cast<int>(detail->reason));
    }
    return take(std::move(result));
}

JsonToken fixture(std::string_view matrix = rotation, std::string_view position = origin,
                  bool second = false) {
    auto input = root("coupled_transient");
    input.field("units").set("length", string("mm"));
    input.field("units").set("angle", string("deg"));
    // Keep the structural aero attachment required by this coupled-physics fixture.
    auto& supports = input.field("model").field("supports").elements;
    auto& support = supports.at(0);
    auto& frame = support.field("motion").field("origin");
    frame.set("origin", parse(position));
    frame.set("rotation", parse(matrix));
    if (second) {
        auto duplicate = support;
        duplicate.set("id", string("999"));
        duplicate.field("region").set("boundary", parse("1"));
        supports.push_back(std::move(duplicate));
    }
    return input;
}

io::SiCase materialize(const JsonToken& source, io::SourceId id) {
    core::ResourcePlan memory{budget};
    const auto wire = source.wire();
    const io::InputLimits syntax{
        100000, 128,    100000, 100000,
        100000, 100000, 100000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    std::optional<io::InputDiagnostic> parsing;
    const auto doc = take(io::read_input_document(
        std::as_bytes(std::span{wire.data(), wire.size()}), id, syntax, memory, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(budget);
    io::SchemaDiagnostic schema_diagnostic;
    const auto schema = checked(
        io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, memory, &schema_diagnostic),
        schema_diagnostic);

    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(budget);
    const io::SiLimits si_limits{si_caps, exact_limits()};
    io::SiDiagnostic si_diagnostic;
    const auto catalog =
        checked(io::make_si_catalog(schema, si_limits, memory, &si_diagnostic), si_diagnostic);
    const auto count = take(catalog.size());
    std::vector<io::SiRequirementRow> rows;
    rows.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t i = 0; i < count; ++i) {
        const auto entry = take(catalog.entry(i));
        rows.push_back({entry.key,
                        {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                        io::SiDecimalInput{"1e-27"},
                        io::SiDecimalInput{"0"}});
    }
    const auto request = checked(
        io::make_si_request(catalog, rows, si_limits, memory, &si_diagnostic), si_diagnostic);
    return checked(
        io::convert_case_to_si(schema, catalog, request, si_limits, memory, &si_diagnostic),
        si_diagnostic);
}

cf::CaseMotionField field(const io::SiCase& source, core::Index index = 0) {
    const auto& typed = take(source.source_case()).get();
    const auto supports = take(typed.root()).model.record().supports;
    const auto candidate = take(supports.at(index)).record();
    const auto* trace = std::get_if<io::schema::SupportTraceConstraint>(&candidate);
    REQUIRE(trace != nullptr);
    return trace->motion;
}

Fraction endpoint(const m::Bound& value) {
    const auto shift = value.power;
    const Integer numerator{value.mantissa};
    if (shift >= 0)
        return Fraction{numerator * (Integer{1} << static_cast<unsigned>(shift))};
    return Fraction{numerator} / (Integer{1} << static_cast<unsigned>(-shift));
}

void encloses(const m::Interval& value, const Fraction& expected) {
    const auto lower = endpoint(value.lower);
    const auto upper = endpoint(value.upper);
    CHECK(lower <= expected);
    CHECK(upper >= expected);
    CHECK(upper - lower <= Fraction{1} / (Integer{1} << 240U));
}
} // namespace

TEST_CASE("Typed motion frame owns SI origin and right-handed row-major rotation",
          "[AE-07][FD-02][CN-D03][SI-D07][verification]") {
    const io::SourceId id{38};
    const auto input = fixture(rotation, origin, true);
    const auto ready = [&]() {
        const auto source = materialize(input, id);
        core::ResourcePlan plan{budget};
        return take(cf::CaseMotionFrame::prepare(source, field(source), frame_limits(), plan));
    }();
    const auto frame = take(ready.view());
    const std::array<Fraction, 3> position{10, 20, 30};
    const Fraction c = Fraction{3} / 5;
    const Fraction s = Fraction{4} / 5;
    const std::array<Fraction, 9> q{c, -s, 0, s, c, 0, 0, 0, 1};
    for (std::size_t i = 0; i < position.size(); ++i)
        encloses(frame.origin[i], position[i]);
    for (std::size_t i = 0; i < q.size(); ++i)
        encloses(frame.rotation[i], q[i]);

    // A unit local x translation follows the first physical column of Q0.
    m::Work work{0, 10000, {}};
    cf::motion::VectorJet d{}, phi{};
    d[0].value = m::integer(1);
    const auto composed =
        cf::motion::compose_local(frame.origin, frame.rotation, d, phi, m::pi_constant(work), work);
    encloses(composed.centre[0].value, Fraction{53} / 5);
    encloses(composed.centre[1].value, Fraction{104} / 5);
    encloses(composed.centre[2].value, 30);

    const auto same = materialize(input, id);
    core::ResourcePlan plan{budget};
    REQUIRE(ready.validate_for(same, field(same), frame_limits(), plan));
    const auto same_frame =
        take(cf::CaseMotionFrame::prepare(same, field(same), frame_limits(), plan));
    const auto other_node =
        take(cf::CaseMotionFrame::prepare(same, field(same, 1), frame_limits(), plan));
    REQUIRE(ready.match(same_frame));
    REQUIRE(ready.match(ready));
    rejects(ready.match(other_node), Failure::StaleGeneration);
    rejects(ready.validate_for(same, field(same, 1), frame_limits(), plan),
            Failure::StaleGeneration);
    const auto another = materialize(fixture(rotation, "[10001,20000,30000]", true), id);
    rejects(ready.validate_for(another, field(another), frame_limits(), plan),
            Failure::StaleGeneration);
    rejects(cf::CaseMotionFrame::prepare(same, field(another), frame_limits(), plan),
            Failure::StaleGeneration);
    encloses(take(ready.view()).origin[0], 10);
}

TEST_CASE("Motion frame rejects reflections, exact shears and insufficient work transactionally",
          "[AE-07][FD-02][SI-D07][N4-02][N4-08][negative]") {
    const io::SourceId id{41};
    for (const auto bad : {std::string_view{"[[-1,0,0],[0,1,0],[0,0,1]]"},
                           std::string_view{"[[1,1e-30,0],[0,1,0],[0,0,1]]"}}) {
        const auto source = materialize(fixture(bad), id);
        core::ResourcePlan plan{budget};
        rejects(cf::CaseMotionFrame::prepare(source, field(source), frame_limits(), plan),
                Failure::GeometryInvalid);
        CHECK(source.source_case().has_value());
    }

    const auto source = materialize(fixture(), id);
    const auto selected = field(source);
    core::ResourcePlan original{budget};
    auto frame = take(cf::CaseMotionFrame::prepare(source, selected, frame_limits(), original));
    const auto measured = take(frame.info()).work;
    REQUIRE(measured > 0);
    auto bounds = frame_limits();
    bounds.maximum_work = measured;
    core::ResourcePlan exact_bytes{original.allocated_bytes()};
    REQUIRE(cf::CaseMotionFrame::prepare(source, selected, bounds, exact_bytes));
    CHECK(exact_bytes.remaining_bytes() == 0);
    core::ResourcePlan one_byte_short{original.allocated_bytes() - 1};
    rejects(cf::CaseMotionFrame::prepare(source, selected, bounds, one_byte_short),
            Failure::ResourceExceeded);
    core::ResourcePlan retries{budget};
    --bounds.maximum_work;
    rejects(cf::CaseMotionFrame::prepare(source, selected, bounds, retries),
            Failure::ResourceExceeded);
    auto scratch = frame_limits();
    --scratch.maximum_scratch_bytes;
    core::ResourcePlan none{0};
    rejects(cf::CaseMotionFrame::prepare(source, selected, scratch, none),
            Failure::ResourceExceeded);
    CHECK(none.allocated_bytes() == 0);
    // A token from a different jthread is stopped before every entry starts.
    std::jthread stopped{[](const std::stop_token&) {}};
    REQUIRE(stopped.request_stop());
    rejects(cf::CaseMotionFrame::prepare(source, selected, frame_limits(), none,
                                         stopped.get_stop_token()),
            Failure::Cancelled);
    rejects(frame.validate_for(source, selected, frame_limits(), none, stopped.get_stop_token()),
            Failure::Cancelled);
    rejects(frame.match(frame, stopped.get_stop_token()), Failure::Cancelled);
    REQUIRE(frame.validate_for(source, selected, frame_limits(), retries));
    encloses(take(frame.view()).rotation[3], Fraction{4} / 5);
}
