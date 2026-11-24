#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_kernel.hpp"
#include "coupling/case_material_chart.hpp"
#include "coupling/case_material_reference.hpp"
#include "coupling/case_material_rows.hpp"
#include "coupling/case_pg_reference.hpp"
#include "coupling/case_pg_support_motion.hpp"
#include "io/case_schema_test_adapter.hpp"
#include "structure/pg_conditioning_test_support.hpp"

#include <algorithm>
#include <cmath>
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
    const auto stop = partition_test::requested_stop_token();
    auto cancelled = limits;
    cancelled.cancellation = stop;
    const auto entry_stop = c::prepare_case_pg_reference(si, 1, cancelled, workspace, plan);
    REQUIRE_FALSE(entry_stop);
    CHECK(entry_stop.error().code == core::FailureCode::Cancelled);
    cancelled.cancellation = {};
    cancelled.partition.cancellation = stop;
    const auto stage_stop = c::prepare_case_pg_reference(si, 1, cancelled, workspace, plan);
    REQUIRE_FALSE(stage_stop);
    CHECK(stage_stop.error().code == core::FailureCode::Cancelled);
    CHECK(take(ready.identity()) == identity);
}

TEST_CASE("Case material chart owns the selected PG and preserves both sides at a shared vertex",
          "[verification][negative][AE-02][SH-D02][FD-06][V-MULTI]") {
    auto wire = a::root("coupled_flutter");
    wire.field("units").set("length", a::string("m"));
    auto& model = wire.field("model");
    auto& geometry = model.field("geometry").elements.at(0);
    geometry.set("reference_length", a::parse("2"));
    geometry.set("frame", a::parse(R"({"origin":[3,-4,5],"rotation":[[1,0,0],[0,1,0],[0,0,1]]})"));
    geometry.set("degree", a::parse("[0,0]"));
    geometry.set("coefficients", a::parse("[]"));
    geometry.field("thickness").set("values", a::parse("[0.25,0.25,0.25,0.25]"));
    auto thin_geometry = a::parse(R"({
        "kind":"surface_graph","id":"2","reference_length":2,
        "frame":{"origin":[3.5,-4,4.9375],"rotation":[[1,0,0],[0,1,0],[0,0,1]]},
        "x0":{"numerator":"0","denominator":"1"},
        "width":{"numerator":"1","denominator":"2"},
        "skew":{"numerator":"0","denominator":"1"},
        "degree":[0,1],"coefficients":[
          {"powers":[0,1],"coefficient":{"numerator":"1","denominator":"16"}}],
        "domain":[],"boundaries":[],"continuity":[]})");
    thin_geometry.set("domain", geometry.field("domain"));
    thin_geometry.set("boundaries", geometry.field("boundaries"));
    model.field("geometry").elements.push_back(std::move(thin_geometry));
    wire.field("numerics").field("patch").set("geometry_error", a::parse("0.000001"));
    auto& layers = model.field("shells").elements.at(0).field("layers").elements;
    auto upper_layer = layers.at(0);
    layers.at(0).set(
        "interval",
        a::parse(R"([{"numerator":"-1","denominator":"1"},{"numerator":"0","denominator":"1"}])"));
    upper_layer.set("id", a::string("9"));
    upper_layer.set(
        "interval",
        a::parse(R"([{"numerator":"0","denominator":"1"},{"numerator":"1","denominator":"1"}])"));
    layers.push_back(std::move(upper_layer));
    auto& surface = model.field("aero_surfaces").elements.at(0);
    surface.set("span", a::parse("[0,0.5,1]"));
    surface.set("geometry", a::string("2"));
    auto& attachment = surface.field("attachment");
    auto piece = a::parse(R"({
        "chord_cells":[0,1],"span_cells":[0,1],"layer":"1","side":"interior","terms":[
          {"powers":[0,0],"coefficient":[{"numerator":"1","denominator":"4"},
            {"numerator":"0","denominator":"1"},{"numerator":"-1","denominator":"2"}]},
          {"powers":[1,0],"coefficient":[{"numerator":"1","denominator":"2"},
            {"numerator":"0","denominator":"1"},{"numerator":"0","denominator":"1"}]},
          {"powers":[0,1],"coefficient":[{"numerator":"0","denominator":"1"},
            {"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"1"}]}]})");
    auto other = piece;
    other.set("layer", a::string("9"));
    other.set("span_cells", a::parse("[1,2]"));
    auto map = a::parse(R"({"pieces":[]})");
    map.field("pieces").elements = {piece, other};
    attachment.set("material_map", map);
    core::ResourcePlan plan{profile};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    const c::CaseMaterialChartLimits limits{case_pg_reference_limits(), query_limits(), 8};
    const auto [ready, thin] = [&] {
        const auto si = checked_case(wire);
        return std::pair{take(c::CaseMaterialChart::prepare(si, 1, limits, workspace, plan)),
                         take(c::CaseReferenceSurface::prepare(si, 2, limits.chart, plan))};
    }();
    CHECK(take(ready.identity()).shell_id == 1);
    CHECK(take(ready.identity()).surface_id == 1);
    const auto left = take(ready.vertex(0, 1, 0, limits.chart, plan));
    const auto right = take(ready.vertex(0, 1, 1, limits.chart, plan));
    CHECK(left.material == right.material);
    CHECK(q(left.material[0]) == Q{3} / 8);
    CHECK(q(left.material[1]) == Q{1} / 2);
    CHECK(q(left.material[2]) == 0);
    CHECK(left.layer == core::LayerId{1});
    CHECK(right.layer == core::LayerId{9});
    CHECK(left.side == c::CaseMaterialSide::Upper);
    CHECK(right.side == c::CaseMaterialSide::Lower);
    const auto collocation = take(ready.collocation(0, 1, limits.chart, plan));
    CHECK(q(collocation.aerodynamic[0]) == Q{3} / 4);
    CHECK(q(collocation.aerodynamic[1]) == Q{3} / 4);
    CHECK(q(collocation.material[0]) == Q{5} / 8);
    CHECK(q(collocation.material[2]) == Q{1} / 4);
    CHECK(q(collocation.chord_tangent[0]) == Q{1} / 2);
    CHECK(q(collocation.span_tangent[1]) == 1);
    CHECK(q(collocation.span_tangent[2]) == 1);
    CHECK(collocation.side == c::CaseMaterialSide::Interior);
    const auto located =
        take(c::qualify_case_material_collocation_reference(ready, thin, 0, 1, limits.chart, plan));
    const auto scalar = [](const PgReferenceScalar& value) {
        return q(value.rational) + q(value.radical);
    };
    CHECK(q(located.reference.g) == 1);
    CHECK(q(located.reference.s) == Q{5} / 8);
    CHECK(q(located.reference.t) == Q{3} / 4);
    CHECK(q(located.reference.zeta) == Q{1} / 4);
    // The flat support gives X = (3+2s, -4+2t, 5+zeta/8).
    CHECK(scalar(located.reference.position_metres[0]) == Q{17} / 4);
    CHECK(scalar(located.reference.position_metres[1]) == Q{-5} / 2);
    CHECK(scalar(located.reference.position_metres[2]) == Q{161} / 32);
    if (!located.tangents) {
        FAIL("The PG collocation must carry its exact tangent pair");
        return;
    }
    CHECK(scalar((*located.tangents)[0][0]) == 1);
    CHECK(scalar((*located.tangents)[0][1]) == 0);
    CHECK(scalar((*located.tangents)[0][2]) == 0);
    CHECK(scalar((*located.tangents)[1][0]) == 0);
    CHECK(scalar((*located.tangents)[1][1]) == 2);
    CHECK(scalar((*located.tangents)[1][2]) == Q{1} / 8);
    const auto left_reference =
        take(c::qualify_case_material_vertex_reference(ready, thin, 0, 1, 0, limits.chart, plan));
    const auto right_reference =
        take(c::qualify_case_material_vertex_reference(ready, thin, 0, 1, 1, limits.chart, plan));
    CHECK_FALSE(left_reference.tangents);
    CHECK_FALSE(right_reference.tangents);
    CHECK(std::get<PgReferenceHorizontalMeasure>(left_reference.reference.measure).side ==
          PgPartitionFaceSide::Upper);
    CHECK(std::get<PgReferenceHorizontalMeasure>(right_reference.reference.measure).side ==
          PgPartitionFaceSide::Lower);
    CHECK(left_reference.reference.partition_record != right_reference.reference.partition_record);
    CHECK(scalar(left_reference.reference.position_metres[2]) ==
          scalar(right_reference.reference.position_metres[2]));

    CHECK(q(take(ready.vertex(1, 2, 1, limits.chart, plan)).aerodynamic[0]) == 1);
    const auto misplaced = ready.vertex(0, 2, 0, limits.chart, plan);
    REQUIRE_FALSE(misplaced);
    CHECK(misplaced.error().code == core::FailureCode::InvalidInput);
    const auto wrong = [&](core::FailureCode expected) {
        const auto si = checked_case(wire);
        const auto failure = c::CaseMaterialChart::prepare(si, 1, limits, workspace, plan);
        REQUIRE_FALSE(failure);
        CHECK(failure.error().code == expected);
    };
    SECTION("Legacy syntax needs an explicit map at physical preparation") {
        attachment.erase("material_map");
        wrong(core::FailureCode::InvalidInput);
    }
    SECTION("A constant shift on only one side cannot match the seam") {
        attachment.field("material_map")
            .field("pieces")
            .elements.at(1)
            .field("terms")
            .elements.at(0)
            .field("coefficient")
            .elements.at(0)
            .set("numerator", a::string("2"));
        wrong(core::FailureCode::GeometryInvalid);
    }
    SECTION("A point on the wrong layer side cannot qualify") {
        attachment.field("material_map")
            .field("pieces")
            .elements.at(0)
            .set("side", a::string("upper"));
        wrong(core::FailureCode::GeometryInvalid);
    }
    SECTION("An outside corner proves that the full panel does not fit the support domain") {
        for (auto& row : attachment.field("material_map").field("pieces").elements) {
            auto& width = row.field("terms").elements.at(1).field("coefficient").elements.at(0);
            width.set("numerator", a::string("7"));
            width.set("denominator", a::string("8"));
        }
        // The centre is inside the PG square, while xi=1 maps to s=9/8.
        wrong(core::FailureCode::GeometryInvalid);
    }
    SECTION("A hidden zero Jacobian invalidates the entire panel") {
        // s = 1/4 + xi/2 - xi^2/2. The derivative vanishes at the panel centre.
        const auto quadratic = a::parse(R"({"powers":[2,0],"coefficient":[
          {"numerator":"-1","denominator":"2"},{"numerator":"0","denominator":"1"},
          {"numerator":"0","denominator":"1"}]})");
        for (auto& row : attachment.field("material_map").field("pieces").elements)
            row.field("terms").elements.push_back(quadratic);
        wrong(core::FailureCode::GeometryInvalid);
    }
    SECTION("The conditioned original row carries position velocity and transpose work") {
        const auto si = checked_case(wire);
        auto physical = take(c::prepare_case_pg_reference(si, 1, limits.support, workspace, plan));
        INFO("Case PG physical reference prepared");
        const auto original = take(physical.source());
        const PgRankPreparationLimits rank_limits{profile, steps, 128, {}};
        auto rank = take(PgRankPreparation::prepare(original.dictionary, original.partition,
                                                    rank_limits, workspace, plan));
        auto basis = take(PgReferenceBasis::prepare(std::move(physical), std::move(rank),
                                                    rank_limits, workspace, plan));
        auto gram = take(
            PgFlatGram::prepare(std::move(basis), 256, {limits.chart, steps}, workspace, plan));
        auto accepted = take(std::move(gram).generate_and_certify({320, 246}, {profile, steps, {}},
                                                                  workspace, plan));
        const auto prepared =
            take(PgPreparedPatch::prepare(std::move(accepted), limits.chart, workspace, plan));
        const auto patch_id = take(prepared.info()).identity;
        constexpr auto max_error = 1e-9;
        for (const auto precision : {core::Precision::Binary64, core::Precision::TwoWord}) {
            CAPTURE(static_cast<int>(precision));
            const PgPreparedReferenceProjection projection{precision, 1e-13, 1e-13,
                                                           1e-13,     1e-13, 1e-13};
            const auto rows = take(c::CaseMaterialRows::prepare(
                ready, thin, prepared, c::CaseMaterialCollocationSite{0, 1}, projection, 1e-13,
                limits.chart, plan));
            INFO("Collocation rows prepared");
            namespace aero = uvlas::aero::detail;
            const c::CaseMaterialLoadKey original_key{aero::PhysicalFaceId{71}, aero::FacetId{81},
                                                      0, aero::PointId{91}};
            const c::CaseMaterialLoadKey other_key{original_key.face, original_key.facet, 1,
                                                   aero::PointId{92}};
            const auto original_vertex = take(c::CaseMaterialRows::prepare(
                ready, thin, prepared, c::CaseMaterialVertexSite{0, 1, 1}, projection, 1e-13,
                limits.chart, plan, original_key));
            const auto other_vertex = take(c::CaseMaterialRows::prepare(
                ready, thin, prepared, c::CaseMaterialVertexSite{1, 1, 1}, projection, 1e-13,
                limits.chart, plan, other_key));
            INFO("Original vertex row prepared");
            const bool paired = precision == core::Precision::TwoWord;
            const auto capacity =
                static_cast<core::Index>(prepared.reference_coordinates().size() * 3);
            const auto destination =
                take(PgVolumePositionResult::allocate(capacity, precision, plan));
            const auto memory = std::max(destination.workspace_bytes(), core::Index{65536});
            core::Workspace action_workspace{
                take(plan.allocate<std::byte>(memory, core::Operation::EvaluateStep))};
            std::vector<double> q_high, q_low, v_high, v_low;
            const std::array speed{Q{1}, Q{-2}, Q{3}};
            const auto append = [&](const Q& exact, std::vector<double>& high,
                                    std::vector<double>& low) {
                high.push_back(exact.convert_to<double>());
                if (paired) {
                    const auto first =
                        take(ExactRational::from_binary64(high.back(), limits.chart.algebra, plan));
                    const Q residual = exact - q(first);
                    low.push_back(residual.convert_to<double>());
                }
            };
            for (std::size_t j = 0; j < prepared.reference_coordinates().size(); ++j)
                for (std::size_t axis = 0; axis < speed.size(); ++axis) {
                    append(q(prepared.reference_coordinates()[j][axis]), q_high, q_low);
                    append(q(prepared.translation()[j]) * speed[axis], v_high, v_low);
                }
            const PgCartesianWords q_state{precision, q_high, q_low};
            const PgCartesianWords v_state{precision, v_high, v_low};
            const c::CaseMaterialActionLimits apply{steps, max_error, max_error, max_error, {}};
            const auto position =
                take(rows.evaluate(patch_id, q_state, v_state, apply, action_workspace, plan));
            const auto word = [](PgCartesianWords words, std::size_t index) {
                return words.high[index] + (words.low.empty() ? 0.0 : words.low[index]);
            };
            CHECK(std::abs(word(position.position.words(), 0) - 17. / 4.) < max_error);
            CHECK(std::abs(word(position.position.words(), 1) + 5. / 2.) < max_error);
            CHECK(std::abs(word(position.position.words(), 2) - 161. / 32.) < max_error);
            CHECK(std::abs(word(position.velocity.words(), 0) - 1.) < max_error);
            CHECK(std::abs(word(position.velocity.words(), 1) + 2.) < max_error);
            CHECK(std::abs(word(position.velocity.words(), 2) - 3.) < max_error);
            if (!position.tangents) {
                FAIL("Prepared rows must retain the collocation tangent pair");
                return;
            }
            CHECK(std::abs(word(position.tangents->words(), 0) - 1.) < max_error);
            CHECK(std::abs(word(position.tangents->words(), 4) - 2.) < max_error);
            CHECK(std::abs(word(position.tangents->words(), 5) - 1. / 8.) < max_error);
            const std::array force{0., 0., 2.};
            const std::array zero{0., 0., 0.};
            const PgCartesianWords loads{precision, force,
                                         paired ? std::span<const double>{zero}
                                                : std::span<const double>{}};
            const auto back = take(original_vertex.pullback(patch_id, loads, {steps, max_error, {}},
                                                            action_workspace, plan));
            CHECK(std::get<c::CaseMaterialVertexSite>(back.original_site).original_span_cell == 1);
            CHECK(back.surface.shell_id == 1);
            CHECK(back.patch == patch_id);
            const auto dual_words = back.force.words();
            double work = 0;
            for (std::size_t k = 0; k < dual_words.high.size(); ++k)
                work += word(dual_words, k) * word(v_state, k);
            CHECK(std::abs(work - 6.) < max_error);
            namespace interval = aero::p1;
            constexpr auto force_radius = 0x1p-35;
            aero::LocalRootScratch enclosure_scratch;
            const auto zero_force =
                aero::enclose_p1_value(interval::integer(0), precision, enclosure_scratch);
            const interval::Interval z_force{aero::decode_p1_word(2. - force_radius).lower,
                                             aero::decode_p1_word(2. + force_radius).upper};
            const c::DiscreteLoadVector enclosed_force{
                zero_force, zero_force,
                aero::enclose_p1_value(z_force, precision, enclosure_scratch)};
            const auto contribution = [&](const c::CaseMaterialLoadKey& key) {
                return c::DiscreteOriginalLoadContribution{
                    key.face, key.facet, key.local_vertex, key.canonical, {}, enclosed_force};
            };
            const auto first_load = contribution(original_key);
            const auto second_load = contribution(other_key);
            const auto enclosed = take(original_vertex.pullback_enclosed(
                patch_id, first_load, {steps, max_error, {}}, plan));
            const auto second = take(other_vertex.pullback_enclosed(patch_id, second_load,
                                                                    {steps, max_error, {}}, plan));
            CHECK(second.original_site.chord_column == 1);
            const auto before_swap = plan.allocated_bytes();
            check_failure(original_vertex.pullback_enclosed(patch_id, second_load,
                                                            {steps, max_error, {}}, plan),
                          core::FailureCode::StaleGeneration);
            check_failure(
                other_vertex.pullback_enclosed(patch_id, first_load, {steps, max_error, {}}, plan),
                core::FailureCode::StaleGeneration);
            CHECK(plan.allocated_bytes() == before_swap);
            CHECK(enclosed.original_site.original_span_cell == 1);
            double center_work = 0., work_radius = 0.;
            const auto scaled = [](const aero::ScaledEnclosure& value) {
                const auto center = std::visit(
                    [](const auto& x) {
                        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(x)>,
                                                     aero::SingleWord>)
                            return x.value;
                        else
                            return x.high + x.low;
                    },
                    value.center());
                return std::scalbn(center, value.exponent());
            };
            for (std::size_t column = 0; column < enclosed.force.size(); ++column) {
                const auto& row = enclosed.force.values()[column];
                if (!row) {
                    FAIL("A structural column force must hold its enclosure");
                    return;
                }
                for (std::size_t axis = 0; axis < row->size(); ++axis) {
                    const auto& value = (*row)[axis];
                    const auto motion = word(v_state, 3 * column + axis);
                    center_work += motion * scaled(value);
                    work_radius += std::abs(motion) * std::scalbn(value.radius(), value.exponent());
                }
            }
            CHECK(std::abs(center_work - 6.) < max_error);
            CHECK(work_radius >= 3 * force_radius);

            auto stale = patch_id;
            stale.bytes.front() ^= std::byte{1};
            check_failure(original_vertex.pullback(stale, loads, {steps, max_error, {}},
                                                   action_workspace, plan),
                          core::FailureCode::StaleGeneration);
        }
    }

    SECTION("Only the second triangle on the shared edge needs retained storage") {
        const auto check_shared_edge = [&](bool subdivide) {
            auto modified = wire;
            auto& thickness =
                modified.field("model").field("geometry").elements.at(0).field("thickness");
            auto knots = a::parse(R"([
              [{"numerator":"0","denominator":"1"},{"numerator":"0","denominator":"1"}],
              [{"numerator":"5","denominator":"8"},{"numerator":"0","denominator":"1"}],
              [{"numerator":"1","denominator":"1"},{"numerator":"0","denominator":"1"}],
              [{"numerator":"0","denominator":"1"},{"numerator":"1","denominator":"1"}],
              [{"numerator":"5","denominator":"8"},{"numerator":"1","denominator":"1"}],
              [{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"1"}]])");
            if (subdivide)
                knots.elements.push_back(a::parse(R"(
                  [{"numerator":"7","denominator":"8"},{"numerator":"1","denominator":"4"}])"));
            thickness.set("vertices", knots);
            thickness.set("triangles",
                          a::parse(subdivide ? "[[0,1,4],[0,4,3],[1,2,6],[2,5,6],[5,1,6],[1,5,4]]"
                                             : "[[0,1,4],[0,4,3],[1,2,5],[1,5,4]]"));
            thickness.set("values", a::parse(subdivide ? "[0.25,0.25,0.25,0.25,0.25,0.25,0.25]"
                                                       : "[0.25,0.25,0.25,0.25,0.25,0.25]"));
            const auto si = checked_case(modified);
            const auto support =
                take(c::CaseMaterialChart::prepare(si, 1, limits, workspace, plan));
            const auto source = take(take(support.support())->source());
            const auto triangles = take(source.partition.triangles());
            const auto result = take(
                c::evaluate_case_material_collocation_reference(support, 0, 1, limits.chart, plan));
            CHECK(q(result.reference.s) == Q{5} / 8);
            CHECK(q(result.reference.t) == Q{3} / 4);
            REQUIRE(result.alternate_count == 1);
            CHECK(result.alternate.size() == 1);
            CHECK(result.alternate.values().front().partition_record !=
                  result.reference.partition_record);
            return triangles.size();
        };
        const auto original = check_shared_edge(false);
        CHECK(check_shared_edge(true) > original);
    }

    SECTION("Matching point position cannot hide an incorrect material tangent") {
        auto& independent = model.field("geometry").elements.at(1);
        independent.set("coefficients", a::parse(R"([
          {"powers":[0,0],"coefficient":{"numerator":"-3","denominator":"64"}},
          {"powers":[0,1],"coefficient":{"numerator":"1","denominator":"8"}}])"));
        const auto si = checked_case(wire);
        const auto new_ready = take(c::CaseMaterialChart::prepare(si, 1, limits, workspace, plan));
        const auto new_thin = take(c::CaseReferenceSurface::prepare(si, 2, limits.chart, plan));
        const auto at_same_position = take(new_thin.evaluate(
            located.chart.aerodynamic[0], located.chart.aerodynamic[1], limits.chart, plan));
        CHECK(q(at_same_position.position[2]) == Q{161} / 32);
        const auto stale = c::qualify_case_material_collocation_reference(ready, new_thin, 0, 1,
                                                                          limits.chart, plan);
        REQUIRE_FALSE(stale);
        CHECK(stale.error().code == core::FailureCode::StaleGeneration);
        const auto mismatch = c::qualify_case_material_collocation_reference(new_ready, new_thin, 0,
                                                                             1, limits.chart, plan);
        REQUIRE_FALSE(mismatch);
        CHECK(mismatch.error().code == core::FailureCode::GeometryInvalid);
    }
    SECTION("The subdivision allowance distinguishes uncertainty from a proved domain") {
        const auto nonlinear = a::parse(R"([
          {"powers":[2,0],"coefficient":[{"numerator":"-1","denominator":"2"},
             {"numerator":"0","denominator":"1"},{"numerator":"0","denominator":"1"}]},
          {"powers":[3,0],"coefficient":[{"numerator":"1","denominator":"3"},
             {"numerator":"0","denominator":"1"},{"numerator":"0","denominator":"1"}]}])");
        // ds/dxi = (xi-1/2)^2 + 1/64. The degree-six whole-panel Bernstein
        // coefficients straddle zero, while the two half-panels certify positivity.
        for (auto& row : attachment.field("material_map").field("pieces").elements) {
            auto& terms = row.field("terms").elements;
            terms.at(1).field("coefficient").elements.at(0).set("numerator", a::string("17"));
            terms.at(1).field("coefficient").elements.at(0).set("denominator", a::string("64"));
            terms.insert(terms.end(), nonlinear.elements.begin(), nonlinear.elements.end());
        }
        const auto si = checked_case(wire);
        auto no_split = limits;
        no_split.maximum_subdivisions = 0;
        const auto unresolved = c::CaseMaterialChart::prepare(si, 1, no_split, workspace, plan);
        REQUIRE_FALSE(unresolved);
        CHECK(unresolved.error().code == core::FailureCode::PrecisionUnresolved);
        const auto proved = take(c::CaseMaterialChart::prepare(si, 1, limits, workspace, plan));
        const auto right_mid = take(proved.collocation(0, 1, limits.chart, plan));
        CHECK(q(right_mid.chord_tangent[0]) == Q{5} / 64);
    }
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
