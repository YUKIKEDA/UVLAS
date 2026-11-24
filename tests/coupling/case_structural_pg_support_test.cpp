#include "coupling/case_pg_support_motion.hpp"
#include "io/case_schema_test_adapter.hpp"
#include "numerics/sparse_test_support.hpp"
#include "structure/pg_conditioning_test_support.hpp"

#include <algorithm>
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

io::SiCase checked_case(const a::JsonToken& input) {
    core::ResourcePlan plan{case_budget};
    const auto wire = input.wire();
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
    if (!checked) {
        const auto diagnostic = schema_diagnostic.detail();
        FAIL("Motion schema fixture "
             << static_cast<unsigned>(checked.error().code) << " "
             << (diagnostic ? static_cast<unsigned>(diagnostic->reason) : 999U) << " "
             << schema_diagnostic.path().value_or("") << " "
             << (diagnostic ? diagnostic->expected_field_name().value_or("") : ""));
    }
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

} // namespace

TEST_CASE("Actual structural Support declarations use the case clock without aerodynamic inputs",
          "[SH-D05][SH-D06][N4-04][verification][negative]") {
    auto input = a::root("structural_static");
    input.field("units").set("angle", a::string("deg"));
    input.set("initial",
              a::parse(R"({"kind":"equilibrium","branch":"1","perturbation":{"kind":"none"}})"));
    auto& analysis = input.field("analysis");
    analysis.set("kind", a::string("transient"));
    analysis.erase("speed");
    analysis.set("physics", a::string("structural"));
    analysis.set(
        "intervals",
        a::parse(
            R"([{"start":0,"end":1,"step":1,"restart":"initial"},{"start":1,"end":2,"step":1,"restart":"event"}])"));
    analysis.set("wake", a::parse(R"({"kind":"none"})"));
    analysis.set("loads", a::array());
    input.set("flow", a::parse(R"({"kind":"none"})"));
    auto& model = input.field("model");
    model.set("aero_surfaces", a::array());
    model.field("geometry").elements[0].field("frame").set("origin", a::parse("[2,-1,0]"));
    auto& scales = input.field("numerics").field("scales");
    scales.elements.push_back(a::scale("structural_velocity"));
    auto declaration = a::fragment("Support.trace_constraint");
    declaration.set("id", a::string("17"));
    auto& motion = declaration.field("motion");
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
    auto duplicate = declaration;
    duplicate.set("id", a::string("71"));
    model.set("supports", a::array());
    model.field("supports").elements.push_back(std::move(duplicate));
    model.field("supports").elements.push_back(std::move(declaration));

    core::ResourcePlan memory{profile};
    p::CurveWorkspace motion_workspace;
    const c::CaseUniformFlowLimits uniform{motion_exact_limits(), motion_limits.curve, 100};
    const c::CaseClockedPrescribedLimits clocks{
        {{uniform, 8, 8}, {uniform, 2048, 50'000}},
        {{motion_exact_limits(), 10'000, c::case_motion_frame_interval_scratch_bytes},
         motion_limits,
         2048,
         50'000},
        case_budget};
    Fixture f;
    const std::array boundary{f.point("0", "0"), f.point("0", "1")};
    f.names.push_back({"leading", boundary});
    f.layers.clear();
    f.layers.push_back({core::LayerId{1}, f.zeta[0], f.zeta[2]});
    const auto zero = f.r("0"), one = f.r("1");
    std::array<ExactRational, 9> rotation;
    for (const auto i : {0U, 4U, 8U})
        rotation[i] = f.r("1");
    const std::array<ExactRational, 3> origin{f.r("2"), f.r("-1"), f.r("0")};
    core::Workspace ws{core::PlannedBuffer<std::byte>{}};
    const auto patch_ready = patch(f, {rotation, origin, one}, ws, memory);
    const auto ranges = vertical_edges(f, patch_ready, Q{0}, Q{1});
    const auto spans = edge_spans(ranges, zero, one);
    auto arithmetic = query_limits();
    arithmetic.maximum_scalar_operations = 300'000'000;
    arithmetic.maximum_partition_visits = 100'000'000;
    const PgTraceLimits trace_limits{arithmetic, 30000};
    const auto trace = take(PgFullTrace::prepare_support(
        patch_ready, spans, take(patch_ready.info()).identity, patch_ready.reference_coordinates(),
        trace_limits, ws, memory));
    const auto clipped = take(PgFullTrace::prepare_support(
        patch_ready, std::span{spans}.first(1), take(patch_ready.info()).identity,
        patch_ready.reference_coordinates(), trace_limits, ws, memory));
    const auto middle = f.r("1", "2");
    auto shorter = spans;
    shorter.clear();
    for (std::size_t i = 0; i < spans.size(); ++i) {
        const auto& span = spans[i];
        shorter.push_back({span.edge, span.first, span.last, span.along_first,
                           i == 0 ? middle : span.along_last, span.layer});
    }
    const auto shortened = take(PgFullTrace::prepare_support(
        patch_ready, shorter, take(patch_ready.info()).identity,
        patch_ready.reference_coordinates(), trace_limits, ws, memory));
    constexpr PgConstraintId support_id{71}, duplicate_id{17};
    const std::array patches{PgConstraintPatchBinding{numerics::PatchId{1}, patch_ready}};
    auto exact = arithmetic;
    exact.algebra.maximum_value_bits = 65536;
    exact.algebra.maximum_work_bits = 262144;
    const auto make_chart = [&](const PgFullTrace& second) {
        const std::array sources{PgConstraintTraceBinding{support_id, trace, 0, {}},
                                 PgConstraintTraceBinding{duplicate_id, second, 0, {}}};
        auto result = take(PgGlobalConstraints::prepare(
            patches, sources, {arithmetic, 1, 2, 4096, 30000}, ws, memory));
        return take(
            PgSparseChart::prepare(std::move(result), {exact, 140, 140, 12000}, ws, memory));
    };
    auto chart = make_chart(trace);
    const auto id = take(chart.identity());
    const auto count = chart.original().original_rows().size();
    const c::CasePgSupportLimits evaluation{clocks, 2048, 2'000'000, count};
    const c::CaseStructuralSupportLimits caps{
        evaluation, {exact, exact, count}, trace_limits, 2, 100'000};
    constexpr std::array first{support_id, duplicate_id}, reverse{duplicate_id, support_id};
    const auto ready = [&](const a::JsonToken& json, PgSparseChart& source,
                           std::span<const PgConstraintId> ids) {
        const auto typed = checked_case(json);
        return c::prepare_case_structural_pg_support(typed, std::move(source), patches, ids, stamp,
                                                     ws, motion_workspace, caps, memory);
    };
    auto other = input;
    auto& other_supports = other.field("model").field("supports").elements;
    other_supports[0].field("motion").set(
        "rotation_vector",
        a::parse(
            R"({"kind":"pchip","times":[0,2],"values":[[0,0,0],[0,0,178]],"outside":{"kind":"reject"}})"));
    check_failure(ready(other, chart, first), core::FailureCode::Unsupported);
    CHECK(take(chart.identity()) == id);
    auto wrong_region = input;
    wrong_region.field("model")
        .field("supports")
        .elements[0]
        .field("region")
        .set("boundary", a::parse("1"));
    check_failure(ready(wrong_region, chart, first), core::FailureCode::StaleGeneration);
    auto partial_chart = make_chart(clipped);
    check_failure(ready(input, partial_chart, first), core::FailureCode::StaleGeneration);
    auto shortened_chart = make_chart(shortened);
    check_failure(ready(input, shortened_chart, first), core::FailureCode::StaleGeneration);
    check_failure(ready(input, chart, std::array{support_id}), core::FailureCode::Unsupported);
    const auto prepared = take(ready(input, chart, first));
    auto duplicate_chart = make_chart(trace);
    const auto current = take(ready(input, duplicate_chart, reverse));

    const std::array<std::uint32_t, 1> word{1};
    const std::array<std::uint32_t, 2> denom{0, 1U << 12U};
    const p::RationalView unit{1, word, word}, null{0, {}, word};
    const c::CaseFlowComponentBudget tolerance{unit, {1, word, denom}, null};
    std::array<c::CaseFlowComponentBudget, 36> pose{};
    std::array<c::CaseFlowComponentBudget, 9> row_budgets{};
    pose.fill(tolerance);
    row_budgets.fill(tolerance);
    const c::CasePgSupportProjection projection{core::Precision::TwoWord, pose, row_budgets,
                                                row_budgets};
    const auto query = [&](c::CaseFlowClockQuery node) {
        return c::evaluate_case_structural_pg_support(prepared, current, node, Second, projection,
                                                      motion_workspace, evaluation, memory);
    };
    const auto left = take(query({stamp, 0, 1, Left}));
    const auto right = take(query({stamp, 1, 0, Right}));
    CHECK(q(left.node().seconds) == 1);
    CHECK(q(right.node().seconds) == 1);
    CHECK(right.chart() == id);
    REQUIRE(left.rows() == count);
    REQUIRE(right.rows() == count);
    bool unequal_accelerations = false;
    const auto centre = [](const io::SiEnclosure& x) -> Dec {
        const auto& w = std::get<io::ScaledTwoWord>(x);
        return ldexp(Dec(w.high) + Dec(w.low), static_cast<int>(w.exponent2));
    };
    for (std::size_t i = 0; i < count * 3; ++i) {
        CHECK(centre(left.original_rhs()[i]) == centre(right.original_rhs()[i]));
        CHECK(centre(left.original_rhs()[i + 3 * count]) ==
              centre(right.original_rhs()[i + 3 * count]));
        unequal_accelerations =
            unequal_accelerations || centre(left.original_rhs()[i + 6 * count]) !=
                                         centre(right.original_rhs()[i + 6 * count]);
    }
    CHECK(unequal_accelerations);
    auto numerical_chart = make_chart(trace);
    const auto rank = static_cast<core::Index>(numerical_chart.pivots().size());
    sparse_test::Fixture forward_scale{rank}, transpose_scale{rank};
    auto backward = sparse_test::context();
    backward.unknown_layout = core::ScopeId{501};
    backward.residual_layout = core::ScopeId{502};
    const auto forward = forward_scale.scales(), transpose = transpose_scale.scales(backward);
    const auto solver = sparse_test::limits(rank, 3);
    auto candidate = take(PgChartCandidate::prepare(std::move(numerical_chart), forward, transpose,
                                                    {arithmetic, 0x1p-30, solver}, memory));
    namespace lin = numerics::linear;
    const auto request = [&](const lin::Scales& scale) {
        return PgChartCandidateSolve{scale,
                                     forward_scale.precision(),
                                     sparse_test::accuracy(),
                                     solver,
                                     core::RhsGeneration{91},
                                     core::NumericalGeneration{92},
                                     core::NumericalGeneration{93},
                                     1e-5,
                                     2'000'000,
                                     {}};
    };
    const auto from_real = [&](std::size_t rows, std::span<const double> data) {
        const auto n = static_cast<core::Index>(rows);
        auto dense = take(
            numerics::DenseReal::copy({n, 3, n}, data, memory, core::Operation::EvaluateStructure));
        return take(numerics::PrecisionBuffer<double>::binary64(std::move(dense),
                                                                core::NumericalGeneration{78}));
    };
    const auto exact_dec = [](const ExactRational& value) -> Dec {
        const auto x = q(value);
        return Dec(numerator(x).convert_to<std::string>()) /
               Dec(denominator(x).convert_to<std::string>());
    };
    const auto r_columns = candidate.exact().free_columns();
    const auto q_columns = patch_ready.reference_coordinates();
    const auto b_columns = patch_ready.translation();
    const Dec omega = acos(Dec(-1)) / 2;
    using Triple = std::array<Dec, 3>;
    std::vector<std::array<Triple, 3>> full(q_columns.size());
    for (std::size_t j = 0; j < full.size(); ++j) {
        const auto b = exact_dec(b_columns[j]);
        const Dec x = exact_dec(q_columns[j][0]) - b;
        const Dec y = exact_dec(q_columns[j][1]) - 2 * b;
        full[j][0] = {b - y, 3 * b + x, exact_dec(q_columns[j][2])};
        full[j][1] = {-omega * x, (Dec(4) / 3) * b - omega * y, Dec{0}};
        full[j][2] = {omega * omega * y, (Dec(5) / 3) * b - omega * omega * x, Dec{0}};
    }
    std::vector<double> force(3 * full.size());
    force[static_cast<std::size_t>(candidate.exact().pivots().front().dependent_column)] = 1;
    const auto forcing = from_real(full.size(), force);
    const auto reaction = take(candidate.transpose({id, forcing}, request(transpose), ws, memory));
    std::vector<PgCartesianExactPosition> exact_force(full.size());
    const auto selected_force =
        static_cast<std::size_t>(candidate.exact().pivots().front().dependent_column);
    exact_force[selected_force][0] = f.r("1");
    const auto exact_dual =
        take(prepared.support().chart().transpose({id, exact_force}, exact, ws, memory));

    for (std::size_t d = 0; d < 3; ++d) {
        const auto rhs = take(c::retain_case_structural_pg_rhs(
            right, static_cast<p::JetOrder>(d), core::NumericalGeneration{109}, 100'000, memory));
        CHECK(rhs.centre().precision() == core::Precision::TwoWord);
        std::vector<double> weights(3 * r_columns.size());
        for (std::size_t axis = 0; axis < 3; ++axis)
            for (std::size_t j = 0; j < r_columns.size(); ++j)
                weights[axis * r_columns.size() + j] =
                    full[static_cast<std::size_t>(r_columns[j])][d][axis].convert_to<double>();
        const auto r = from_real(r_columns.size(), weights);
        std::vector<PgCartesianExactPosition> exact_independent(r_columns.size());
        for (std::size_t j = 0; j < r_columns.size(); ++j)
            for (std::size_t axis = 0; axis < 3; ++axis)
                exact_independent[j][axis] = take(ExactRational::from_binary64(
                    weights[axis * r_columns.size() + j], exact.algebra, memory));
        const auto compatible = take(
            c::lift_case_structural_pg_support_exact(prepared, right, static_cast<p::JetOrder>(d),
                                                     exact_independent, caps.exact, ws, memory));
        const auto& raw = prepared.support().chart().original();
        const auto recomputed =
            take(raw.apply({raw.layout(), compatible.coordinates.values()}, exact, ws, memory));
        CHECK(std::ranges::equal(recomputed.values(), compatible.original_rhs.values()));
        CHECK(exact_dec(compatible.maximum_rhs_error) <= Dec("1e-10"));
        Q exact_balance = 0;
        for (std::size_t j = 0; j < full.size(); ++j)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                CHECK(abs(exact_dec(compatible.coordinates.values()[j][axis]) - full[j][d][axis]) <=
                      Dec("1e-6"));
                exact_balance +=
                    q(compatible.coordinates.values()[j][axis]) * q(exact_force[j][axis]);
            }
        for (std::size_t j = 0; j < exact_independent.size(); ++j)
            for (std::size_t axis = 0; axis < 3; ++axis)
                exact_balance -= q(exact_independent[j][axis]) *
                                 q(exact_dual.independent_force.values()[j][axis]);
        for (std::size_t i = 0; i < count; ++i)
            for (std::size_t axis = 0; axis < 3; ++axis)
                exact_balance -= q(compatible.original_rhs.values()[i][axis]) *
                                 q(exact_dual.original_multipliers.values()[i][axis]);
        CHECK(exact_balance == 0);
        const auto primal = request(forward);
        const auto result =
            take(c::lift_case_structural_pg_support(candidate, {id, r}, rhs, primal, ws, memory));
        REQUIRE(result.original_residual_upper <= primal.maximum_absolute_original_defect);
        const auto retained_coordinate = [&](std::size_t i) {
            const auto& values = result.coordinates.storage();
            if (const auto* binary = std::get_if<numerics::Binary64<double>>(&values))
                return Dec(binary->values.values()[i]);
            const auto& two = std::get<numerics::TwoWord<double>>(values);
            return Dec(two.hi.values()[i]) + Dec(two.lo.values()[i]);
        };
        for (std::size_t axis = 0; axis < 3; ++axis)
            for (std::size_t j = 0; j < full.size(); ++j)
                CHECK(abs(retained_coordinate(axis * full.size() + j) - full[j][d][axis]) <=
                      Dec("1e-6"));
        const auto defect = take(c::work_case_structural_pg_support(
            candidate, {id, forcing}, {id, r}, rhs, result, reaction, 1e-6, 2'000'000));
        CHECK(defect <= 1e-6);
    }
}
