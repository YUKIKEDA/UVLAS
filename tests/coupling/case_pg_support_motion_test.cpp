#include "coupling/case_pg_reference.hpp"
#include "coupling/case_pg_support_motion.hpp"
#include "io/case_schema_test_adapter.hpp"
#include "structure/pg_conditioning_test_support.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <uvlas/structure/pg_full_trace.hpp>
#include <uvlas/structure/pg_global_constraints.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;

PgPreparedPatch patch(Fixture& f, const PgReferenceFrame& frame, core::Workspace& ws,
                      core::ResourcePlan& plan) {
    const auto limits = query_limits();
    auto basis = make_basis(f, frame, ws, plan);
    auto gram = take(PgFlatGram::prepare(std::move(basis), 256, {limits, steps}, ws, plan));
    auto accepted =
        take(std::move(gram).generate_and_certify({320, 246}, {profile, steps, {}}, ws, plan));
    return take(PgPreparedPatch::prepare(std::move(accepted), limits, ws, plan));
}

const PgIntegrationPartition& partition(const PgPreparedPatch& patch) {
    return take(take(patch.conditioning().reference_basis())->source()).partition;
}

struct EdgeRange {
    core::Index edge{};
    ExactRational first, last;
    bool increasing = true;
};

std::vector<EdgeRange> vertical_edges(Fixture& fixture, const PgPreparedPatch& patch, const Q& s,
                                      const Q& scale) {
    const auto& part = partition(patch);
    const auto vertices = take(part.vertices());
    const auto edges = take(part.edges());
    std::vector<EdgeRange> result;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& edge = edges[i];
        if (!edge.external)
            continue;
        const auto& a = vertices[static_cast<std::size_t>(edge.vertices[0])];
        const auto& b = vertices[static_cast<std::size_t>(edge.vertices[1])];
        if (q(a.s) != s || q(b.s) != s)
            continue;
        const Q x = q(a.t) * scale, y = q(b.t) * scale;
        const auto lo = std::min(x, y), hi = std::max(x, y);
        result.push_back({static_cast<core::Index>(i), rat(fixture, numerator(lo), denominator(lo)),
                          rat(fixture, numerator(hi), denominator(hi)), x < y});
    }
    std::ranges::sort(result, [](const auto& a, const auto& b) { return q(a.first) < q(b.first); });
    return result;
}

std::vector<PgTraceEdgeSpan> edge_spans(std::span<const EdgeRange> ranges,
                                        const ExactRational& zero, const ExactRational& one,
                                        std::optional<core::Index> layer = {}) {
    std::vector<PgTraceEdgeSpan> result;
    result.reserve(ranges.size());
    for (const auto& edge : ranges)
        result.push_back({edge.edge, edge.first, edge.last, edge.increasing ? zero : one,
                          edge.increasing ? one : zero, layer});
    return result;
}

} // namespace

namespace {
using namespace uvlas;
namespace c = coupling::detail;
namespace p = aero::detail::pchip;
namespace a = uvlas::test::schema_adapter;
using enum p::JetOrder;
using enum p::RequestedSide;

constexpr std::uint64_t case_budget = std::uint64_t{1} << 48U;
constexpr p::PchipInputStamp stamp{p::PchipIssuerId{71}, p::PchipInputGeneration{8}};

structure::ExactLimits motion_exact_limits() noexcept {
    return {100000,
            100000,
            100000,
            64,
            64,
            100000,
            200000,
            structure::ExactEstimateAdmission{case_budget},
            structure::ExactTimingMode::Disabled};
}

const c::CaseMotionTimeLimits motion_limits{
    motion_exact_limits(),
    {{128, 100000, case_budget, 10000}, numerics::linear::EstimateAdmission{case_budget}},
    100,
    1000};

io::SiCase checked_case_wire(std::string_view wire) {
    core::ResourcePlan plan{case_budget};
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, case_budget}};
    std::optional<io::InputDiagnostic> parsing;
    const auto doc =
        take(io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                     io::SourceId{19}, syntax, plan, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(case_budget);
    io::SchemaDiagnostic schema_diagnostic;
    auto checked =
        io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, plan, &schema_diagnostic);
    if (!checked)
        FAIL("Motion schema fixture at " << schema_diagnostic.path().value_or(""));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(case_budget);
    const io::SiLimits bounds{si_caps, motion_exact_limits()};
    const auto catalog = take(io::make_si_catalog(*checked, bounds, plan, nullptr));
    std::vector<io::SiRequirementRow> requirements;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto entry = take(catalog.entry(i));
        requirements.push_back({entry.key,
                                {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                                io::SiDecimalInput{"1e-25"},
                                io::SiDecimalInput{"1e-25"}});
    }
    const auto request = take(io::make_si_request(catalog, requirements, bounds, plan, nullptr));
    return take(io::convert_case_to_si(*checked, catalog, request, bounds, plan, nullptr));
}

io::SiCase checked_case(const a::JsonToken& input) {
    return checked_case_wire(input.wire());
}

const c::CasePgReferenceLimits& case_pg_reference_limits() {
    static const auto exact = partition_test::exact_limits();
    static const c::CasePgReferenceLimits limits{{exact, 64, {}},
                                                 {exact, 2000, 20000000, 100000, {}},
                                                 partition_test::limits(),
                                                 query_limits(),
                                                 100000,
                                                 {}};
    return limits;
}

io::schema::AeroAttachmentPrescribed attachment_field(const io::SiCase& source) {
    const auto root = take(take(source.source_case()).get().root());
    const auto model = root.model.record();
    const auto object = take(model.aero_surfaces.at(0)).record().attachment.record();
    if (const auto* motion = std::get_if<io::schema::AeroAttachmentPrescribed>(&object))
        return *motion;
    FAIL("Missing prescribed attachment");
}

a::JsonToken& motion_json(a::JsonToken& root) {
    return root.field("model")
        .field("aero_surfaces")
        .elements[0]
        .field("attachment")
        .field("motion");
}

} // namespace

TEST_CASE("An explicit Case PG shell owns its exact reference geometry in SI",
          "[V-PG][SH-D02][verification][negative]") {
    auto input = a::root();
    input.field("units").set("length", a::string("mm"));
    auto& model = input.field("model");
    auto& geometry = model.field("geometry").elements.front();
    geometry.set("reference_length", a::parse("2000"));
    geometry.set("frame",
                 a::parse(R"({"origin":[3000,-4000,5000],"rotation":[[0,0,1],[1,0,0],[0,1,0]]})"));
    geometry.set(
        "coefficients",
        a::parse(R"([{"powers":[0,0],"coefficient":{"numerator":"1","denominator":"4"}}])"));
    geometry.field("thickness").set("values", a::parse("[50,50,50,50]"));
    auto& layers = model.field("shells").elements.front().field("layers").elements;
    auto upper = layers.front();
    layers.front().set(
        "interval",
        a::parse(R"([{"numerator":"-1","denominator":"1"},{"numerator":"0","denominator":"1"}])"));
    upper.set("id", a::string("9"));
    upper.set(
        "interval",
        a::parse(R"([{"numerator":"0","denominator":"1"},{"numerator":"1","denominator":"1"}])"));
    layers.push_back(std::move(upper));

    core::ResourcePlan plan{profile};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    const auto& limits = case_pg_reference_limits();
    const auto ready = [&] {
        const auto si = checked_case(input);
        return take(c::prepare_case_pg_reference(si, 1, limits, workspace, plan));
    }();
    const auto prepared = take(ready.source());
    const auto layer_data = take(prepared.partition.layers());
    REQUIRE(layer_data.size() == 2);
    CHECK(layer_data[0].id == core::LayerId{1});
    CHECK(layer_data[1].id == core::LayerId{9});
    CHECK(q(take(ready.frame()).length_metres) == 2);
    const auto volume = take(prepared.partition.volumes());
    REQUIRE_FALSE(volume.empty());
    Fixture values;
    const auto zero = values.r("0");
    const auto depth = values.r("1", "2");
    const auto at = take(ready.evaluate(PgReferenceVolumePoint{0, zero, zero, depth}, std::nullopt,
                                        limits.geometry, plan));
    CHECK(q(at.full_thickness_metres) == Q{1} / 20);
    // x = 3 + 2*(1/4) + zeta/40, y = -4 + 2*s, z = 5 + 2*t.
    CHECK(q(at.position_metres[0].rational) + q(at.position_metres[0].radical) ==
          Q{7} / 2 + q(at.zeta) / 40);
    CHECK(q(at.position_metres[1].rational) == -4 + 2 * q(at.s));
    CHECK(q(at.position_metres[2].rational) == 5 + 2 * q(at.t));

    const auto identity = take(ready.identity());
    const auto si = checked_case(input);
    auto reduced = limits;
    reduced.maximum_input_work = 1;
    const auto exhausted = c::prepare_case_pg_reference(si, 1, reduced, workspace, plan);
    REQUIRE_FALSE(exhausted);
    CHECK(exhausted.error().code == core::FailureCode::ResourceExceeded);
    const auto missing = c::prepare_case_pg_reference(si, 77, limits, workspace, plan);
    REQUIRE_FALSE(missing);
    CHECK(missing.error().code == core::FailureCode::InvalidInput);
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    auto cancelled = limits;
    cancelled.cancellation = stop.get_token();
    const auto entry_stop = c::prepare_case_pg_reference(si, 1, cancelled, workspace, plan);
    REQUIRE_FALSE(entry_stop);
    CHECK(entry_stop.error().code == core::FailureCode::Cancelled);
    cancelled.cancellation = {};
    cancelled.partition.cancellation = stop.get_token();
    const auto stage_stop = c::prepare_case_pg_reference(si, 1, cancelled, workspace, plan);
    REQUIRE_FALSE(stage_stop);
    CHECK(stage_stop.error().code == core::FailureCode::Cancelled);
    CHECK(take(ready.identity()) == identity);
}

TEST_CASE("The frozen V-PG Case series prepares its physical PG references",
          "[.frozen-vpg][V-PG][FD-03][SH-D02][verification]") {
    std::ifstream manifest_file{"manifest.json", std::ios::binary};
    REQUIRE(manifest_file);
    const auto manifest = nlohmann::json::parse(manifest_file);
    const auto& cases = manifest.at("cases");
    REQUIRE(cases.size() == 34);
    const auto& limits = case_pg_reference_limits();
    Fixture sample;
    const auto third = sample.r("1", "3");
    const auto half = sample.r("1", "2");
    std::size_t prepared_shells = 0;
    for (const auto& entry : cases) {
        const auto path = entry.at("case").get<std::string>();
        CAPTURE(path);
        std::ifstream case_file{path, std::ios::binary};
        REQUIRE(case_file);
        const std::string wire{std::istreambuf_iterator<char>{case_file},
                               std::istreambuf_iterator<char>{}};
        const auto si = checked_case_wire(wire);
        const auto& grid = entry.at("n");
        const auto count = grid.is_array() ? grid.size() : 1;
        const auto& thickness = entry.at("h");
        const Q h = q(sample.r(thickness.at("numerator").get<std::string>(),
                               thickness.at("denominator").get<std::string>()));
        for (std::size_t shell = 0; shell < count; ++shell) {
            core::ResourcePlan plan{profile};
            core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
            const auto result =
                c::prepare_case_pg_reference(si, shell + 1, limits, workspace, plan);
            if (!result)
                UNSCOPED_INFO("Preparation failure " << static_cast<unsigned>(result.error().code));
            REQUIRE(result);
            const auto& ready = *result;
            const auto source = take(ready.source());
            CHECK(take(source.graph.view()).n ==
                  (grid.is_array() ? grid.at(shell) : grid).get<core::Index>());
            const auto volumes = take(source.partition.volumes());
            const auto triangles = take(source.partition.triangles());
            REQUIRE_FALSE(volumes.empty());
            Q area = 0;
            for (const auto& triangle : triangles)
                area += q(triangle.twice_area) / 2;
            CHECK(area == (path.starts_with("r5-") ? Q{4} / 5 : Q{1}));
            const auto at = take(
                ready.evaluate(PgReferenceVolumePoint{static_cast<core::Index>(volumes.size() / 2),
                                                      third, third, half},
                               std::nullopt, limits.geometry, plan));
            const auto variable = entry.at("variable").get<bool>();
            const Q expected = variable ? h * (1 + q(at.s) / 10 + q(at.t) / 20) : h;
            CHECK(q(at.full_thickness_metres) == expected);
            CHECK(at.layer == core::LayerId{shell + 1});
            ++prepared_shells;
        }
    }
    CHECK(prepared_shells == 36);
}

TEST_CASE("Clocked rigid PG supports preserve original rows and reaction power",
          "[SH-D05][SH-D06][N4-04][AE-07][verification][negative]") {
    auto input = a::root("aerodynamic_prescribed_transient");
    input.field("units").set("time", a::string("s"));
    input.field("units").set("length", a::string("m"));
    input.field("units").set("angle", a::string("deg"));
    input.field("analysis")
        .set(
            "intervals",
            a::parse(
                R"([{"start":0,"end":1,"step":1,"restart":"initial"},{"start":1,"end":2,"step":1,"restart":"event"}])"));
    auto& motion = motion_json(input);
    motion.set("origin", a::parse(R"({"origin":[1,2,3],"rotation":[[0,-1,0],[1,0,0],[0,0,1]]})"));
    motion.set(
        "translation",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[1,0,0],[3,0,0]],"outside":{"kind":"reject"}})"));
    motion.set(
        "rotation_vector",
        a::parse(
            R"({"kind":"pchip","times":[0,2],"values":[[0,0,0],[0,0,180]],"outside":{"kind":"reject"}})"));
    motion.set("events", a::parse("[0,1,2]"));

    core::ResourcePlan memory{profile};
    aero::detail::UniformFlowWorkspace flow_workspace;
    p::CurveWorkspace motion_workspace;
    const c::CaseUniformFlowLimits uniform{motion_exact_limits(), motion_limits.curve, 100};
    const c::CaseClockedPrescribedLimits clocks{
        {{uniform, 8, 8}, {uniform, 2048, 50'000}},
        {{motion_exact_limits(), 10'000, c::case_motion_frame_interval_scratch_bytes},
         motion_limits,
         2048,
         50'000},
        case_budget};
    const auto prepare_case = [&](const a::JsonToken& json) {
        const auto source = checked_case(json);
        return take(c::prepare_case_clocked_prescribed_inputs(source, attachment_field(source),
                                                              stamp, flow_workspace,
                                                              motion_workspace, clocks, memory));
    };
    const auto clock = prepare_case(input), current = prepare_case(input);
    auto other = input;
    motion_json(other).set(
        "rotation_vector",
        a::parse(
            R"({"kind":"pchip","times":[0,2],"values":[[0,0,0],[0,0,178]],"outside":{"kind":"reject"}})"));
    const auto changed = prepare_case(other);

    Fixture f;
    const auto zero = f.r("0"), one = f.r("1");
    std::array<ExactRational, 9> rotation;
    for (const auto i : {0U, 4U, 8U})
        rotation[i] = f.r("1");
    const std::array<ExactRational, 3> origin{f.r("2"), f.r("-1"), f.r("0")};
    core::Workspace work{core::PlannedBuffer<std::byte>{}};
    const auto model = patch(f, {rotation, origin, one}, work, memory);
    const auto ranges = vertical_edges(f, model, Q{0}, Q{1});
    const auto spans = edge_spans(ranges, zero, one);
    auto limits = query_limits();
    limits.maximum_scalar_operations = 300'000'000;
    limits.maximum_partition_visits = 100'000'000;
    const auto support = take(
        PgFullTrace::prepare_support(model, spans, take(model.info()).identity,
                                     model.reference_coordinates(), {limits, 30000}, work, memory));
    constexpr PgConstraintId support_id{71}, duplicate_id{17};
    const std::array patches{PgConstraintPatchBinding{numerics::PatchId{12}, model}};
    const std::array sources{PgConstraintTraceBinding{support_id, support, 0, {}},
                             PgConstraintTraceBinding{duplicate_id, support, 0, {}}};
    auto exact = limits;
    exact.algebra.maximum_value_bits = 65536;
    exact.algebra.maximum_work_bits = 262144;
    const auto new_chart = [&] {
        auto original = take(PgGlobalConstraints::prepare(
            patches, sources, {limits, 1, 2, 4096, 30000}, work, memory));
        return take(
            PgSparseChart::prepare(std::move(original), {exact, 140, 140, 12000}, work, memory));
    };
    auto chart = new_chart();
    const auto identity = take(chart.identity());
    const auto count = chart.original().original_rows().size();
    const PgRigidSupportLimits support_limits{exact, exact, count};
    const auto qualify = [&](PgSparseChart& parent, std::span<const PgConstraintId> selected,
                             const PgRigidSupportLimits& bound) {
        return PgRigidSupportRows::prepare(std::move(parent), patches, selected, bound, work,
                                           memory);
    };
    // At the rest configuration both declarations agree. Moving just one
    // cannot be certified for the same whole rotation as its fixed duplicate.
    check_failure(qualify(chart, std::array{support_id}, support_limits),
                  core::FailureCode::InvalidInput);
    auto short_rows = support_limits;
    --short_rows.maximum_original_rows;
    check_failure(qualify(chart, std::array{support_id, duplicate_id}, short_rows),
                  core::FailureCode::ResourceExceeded);
    CHECK(take(chart.identity()) == identity);
    const auto prepared =
        take(qualify(chart, std::array{support_id, duplicate_id}, support_limits));

    const std::array<std::uint32_t, 1> word{1};
    const std::array<std::uint32_t, 2> denom{0, 1U << 12U};
    const p::RationalView unit{1, word, word}, null{0, {}, word};
    const c::CaseFlowComponentBudget tolerance{unit, {1, word, denom}, null};
    std::array<c::CaseFlowComponentBudget, 36> pose{};
    std::array<c::CaseFlowComponentBudget, 9> components{};
    pose.fill(tolerance);
    components.fill(tolerance);
    const std::array vector{tolerance, tolerance, tolerance};
    const c::CaseKnownFlowProjection flow{core::Precision::TwoWord, vector, vector, components};
    const c::CasePgSupportProjection projection{core::Precision::TwoWord, pose, components,
                                                components};
    const c::CasePgSupportLimits row_limits{clocks, 2048, 2'000'000, count};
    const auto evaluate = [&](const c::CaseClockedPrescribedInputs& live,
                              c::CaseFlowClockQuery node, const c::CasePgSupportLimits& caps) {
        return c::evaluate_case_clocked_pg_support(prepared, clock, live, node, {null, null, null},
                                                   flow, {}, Second, projection, flow_workspace,
                                                   motion_workspace, caps, memory);
    };
    const auto left = take(evaluate(current, {stamp, 0, 1, Left}, row_limits));
    const auto right = take(evaluate(current, {stamp, 1, 0, Right}, row_limits));
    check_failure(evaluate(changed, {stamp, 1, 0, Right}, row_limits),
                  core::FailureCode::StaleGeneration);
    check_failure(evaluate(current, {stamp, 0, 1, Right}, row_limits),
                  core::FailureCode::InvalidInput);
    auto work_short = row_limits;
    work_short.maximum_row_work = 1;
    check_failure(evaluate(current, {stamp, 1, 0, Right}, work_short),
                  core::FailureCode::ResourceExceeded);
    REQUIRE(left.rows == count);
    REQUIRE(right.rows == count);
    REQUIRE(left.orders == 3);
    CHECK(right.chart == identity);
    const auto t = take(left.inputs.motion.time()).get().time();
    const auto tr = take(right.inputs.motion.time()).get().time();
    REQUIRE(t.numerator.size() == 1);
    REQUIRE(tr.numerator.size() == 1);
    CHECK(t.numerator.front() == 1);
    CHECK(tr.numerator.front() == 1);

    const auto dec = [](const ExactRational& value) {
        const auto x = q(value);
        return Dec(numerator(x).convert_to<std::string>()) /
               Dec(denominator(x).convert_to<std::string>());
    };
    const auto decode = [](const io::SiEnclosure& x) {
        const auto& w = std::get<io::ScaledTwoWord>(x);
        const auto exponent = static_cast<int>(w.exponent2);
        return std::array<Dec, 2>{ldexp(Dec(w.high) + Dec(w.low), exponent),
                                  ldexp(Dec(w.radius), exponent)};
    };
    const Dec omega = acos(Dec(-1)) / 2;
    const auto ref = model.reference_coordinates();
    const auto unit_field = model.translation();
    const auto compare = [&](const c::CasePgSupportSample& result, const Dec& acceleration) {
        using Triple = std::array<Dec, 3>;
        std::vector<std::array<Triple, 3>> field(ref.size());
        for (std::size_t j = 0; j < ref.size(); ++j) {
            const auto b = dec(unit_field[j]);
            const Dec x = dec(ref[j][0]) - b, y = dec(ref[j][1]) - 2 * b;
            field[j][0] = {b - y, 3 * b + x, dec(ref[j][2])};
            field[j][1] = {-omega * x, (Dec(4) / 3) * b - omega * y, Dec{0}};
            field[j][2] = {omega * omega * y, acceleration * b - omega * omega * x, Dec{0}};
        }
        std::vector<Triple> force(ref.size());
        Dec power = 0, certified_power = 0, power_radius = 0;
        bool any = false;
        for (std::size_t row = 0; row < count; ++row) {
            std::array<Triple, 3> expected{};
            const auto terms = take(prepared.chart().original().row(static_cast<core::Index>(row)));
            const auto multiplier = [&](std::size_t axis) {
                return static_cast<int>((row + 2 * axis) % 5) - 2;
            };
            for (const auto& term : terms) {
                const auto column = static_cast<std::size_t>(term.column);
                const auto coefficient = dec(term.coefficient);
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    force[column][axis] += coefficient * multiplier(axis);
                    for (std::size_t d = 0; d < 3; ++d)
                        expected[d][axis] += coefficient * field[column][d][axis];
                }
            }
            for (std::size_t d = 0; d < 3; ++d)
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    const auto encl =
                        decode(result.original_rhs.values()[(d * count + row) * 3 + axis]);
                    CHECK(abs(encl[0] - expected[d][axis]) <= encl[1] + Dec("1e-70"));
                    any = any || abs(expected[d][axis]) > Dec("1e-10");
                    if (d == 1) {
                        power += expected[d][axis] * multiplier(axis);
                        certified_power += encl[0] * multiplier(axis);
                        power_radius += encl[1] * std::abs(multiplier(axis));
                    }
                }
        }
        Dec direct = 0;
        for (std::size_t j = 0; j < field.size(); ++j)
            for (std::size_t axis = 0; axis < 3; ++axis)
                direct += force[j][axis] * field[j][1][axis];
        CHECK(any);
        CHECK(abs(direct) > Dec("1e-10"));
        CHECK(abs(power - direct) <= Dec("1e-70"));
        CHECK(abs(certified_power - direct) <= power_radius + Dec("1e-70"));
    };
    compare(left, Dec(1) / 3);
    compare(right, Dec(5) / 3);
    auto fixed_chart = new_chart();
    const auto fixed_rows = take(qualify(fixed_chart, {}, support_limits));
    const auto static_result = take(c::evaluate_case_clocked_pg_support(
        fixed_rows, clock, current, {stamp, 1, 0, Right}, {null, null, null}, flow, {}, Second,
        projection, flow_workspace, motion_workspace, row_limits, memory));
    REQUIRE(static_result.rows == count);
    REQUIRE(static_result.orders == 3);
    for (std::size_t row = 0; row < count; ++row) {
        const auto& old = fixed_rows.chart().original().original_rows()[row].prescribed;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto position = decode(static_result.original_rhs.values()[row * 3 + axis]);
            CHECK(abs(position[0] - dec(old[axis])) <= position[1] + Dec("1e-70"));
            for (std::size_t d = 1; d < 3; ++d) {
                const auto derivative =
                    decode(static_result.original_rhs.values()[(d * count + row) * 3 + axis]);
                CHECK(derivative[0] == 0);
                CHECK(derivative[1] == 0);
            }
        }
    }
}
