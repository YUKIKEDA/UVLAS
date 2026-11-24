#include "../../src/structure/pg_reference_arithmetic.hpp"
#include "../structure/pg_partition_test_support.hpp"
#include "external_pg_regional_thickness.hpp"
#include "source_digest.hpp"

#include <algorithm>
#include <utility>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using namespace uvlas::structure;
using partition_test::take;
using Code = core::FailureCode;

class Reader final : public ReferenceFileReader {
public:
    explicit Reader(std::span<const std::byte> bytes) noexcept : bytes_(bytes) {}

    core::Result<ReferenceFileInfo> info() noexcept override {
        return ReferenceFileInfo{{0, {}}, bytes_.size(), {}, true};
    }

    core::Result<std::size_t> read(std::span<std::byte> destination) noexcept override {
        const auto count = std::min(bytes_.size() - offset_, destination.size());
        std::ranges::copy(bytes_.subspan(offset_, count), destination.begin());
        offset_ += count;
        return count;
    }

private:
    std::span<const std::byte> bytes_;
    std::size_t offset_{};
};

[[nodiscard]] PgReferenceGeometry pg(partition_test::Fixture& input, bool varying,
                                     const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    if (varying) {
        input.heights[3] = input.r("2");
    }
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto dictionary = take(PgCoefficientDictionary::prepare(
        input.chart_graph(), {limits.algebra, 2000, 20'000'000, 100'000, {}}, workspace, plan));
    auto partition = take(PgIntegrationPartition::prepare(input.request(), partition_test::limits(),
                                                          workspace, plan));
    const auto length = input.r("1");
    const std::array position{input.r("0"), input.r("0"), input.r("0")};
    const std::array rotation{input.r("1"), input.r("0"), input.r("0"), input.r("0"), input.r("1"),
                              input.r("0"), input.r("0"), input.r("0"), input.r("1")};
    return take(PgReferenceGeometry::prepare(
        {std::move(input.chart_graph()), std::move(dictionary), std::move(partition)},
        {rotation, position, length}, limits, plan));
}

[[nodiscard]] std::string wire(bool crossing, bool reversed) {
    constexpr std::string_view prefix = R"mesh($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
2 17 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 2 17 41
2 17 2 2
)mesh";
    const auto elements = crossing   ? "17 1 2 4\n41 2 3 4\n"
                          : reversed ? "17 1 3 2\n41 1 4 3\n"
                                     : "17 1 2 3\n41 1 3 4\n";
    return std::string{prefix} + elements + "$EndElements\n";
}

[[nodiscard]] GmshMesh source(std::string_view msh, core::ResourcePlan& plan) {
    const auto bytes = std::as_bytes(std::span{msh});
    Reader reader{bytes};
    auto file =
        take(read_reference_file(reader, source_digest(bytes),
                                 {msh.size(), msh.size() + 1, 3, msh.size()}, plan, {}, nullptr));
    const AcquiredRootSnapshot proof{
        {{{SourceId{90}, 100, {}}, 1, 1}, SourceId{99}, ReferenceRole::GmshMesh, 31},
        std::move(file),
        {}};
    return take(GmshMesh::prepare(proof, {4096, 0, 0, 0, 8, 8, 32'000}, plan, {}));
}

[[nodiscard]] ExternalPlanarPg geometry(const PgReferenceGeometry& reference,
                                        const partition_test::Fixture& coordinates, bool crossing,
                                        bool reversed, const ExternalPgLimits& limits,
                                        core::ResourcePlan& plan) {
    const auto msh = source(wire(crossing, reversed), plan);
    const auto& pts = coordinates.domain;
    constexpr std::array<std::uint64_t, 2> tags{17, 41};
    const auto topology =
        take(GmshSurfaceTopology::prepare(msh, tags, {8, 16, 16, 4, 32'000}, plan, {}));
    const std::array<ExternalPgCorner, 4> points{
        ExternalPgCorner{1, pts[0]}, {2, pts[1]}, {3, pts[2]}, {4, pts[3]}};
    const std::array first{points[0], points[reversed ? 2U : 1U],
                           points[crossing   ? 3U
                                  : reversed ? 1U
                                             : 2U]};
    const std::array second{points[crossing ? 1U : 0U], points[reversed ? 3U : 2U],
                            points[reversed ? 2U : 3U]};
    const auto winding = reversed ? ExternalPgWinding::Reversed : ExternalPgWinding::Forward;
    const std::array faces{ExternalPgFace{17, winding, first}, ExternalPgFace{41, winding, second}};
    return take(ExternalPlanarPg::prepare(reference, msh, topology, SiLengthUnit::Metre, faces,
                                          limits, plan));
}

TEST_CASE("Regional exterior thickness is compared only over positive exact overlap",
          "[FD-03][IO4-04][R-05][N4-02][integration][verification][negative]") {
    bool crossing = false;
    bool reversed = false;
    SECTION("The common diagonal touches but the adjacent thickness field is outside") {}
    SECTION("The reverse source ring does not change the assigned region") {
        reversed = true;
    }
    SECTION("A different exterior diagonal cuts across varying thickness") {
        crossing = true;
    }

    partition_test::Fixture input;
    core::ResourcePlan plan{partition_test::cap};
    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const ExternalPgLimits geometric{arithmetic, 4};
    const ExternalPgThicknessLimits limits{arithmetic, 2, crossing ? 4U : 3U};
    const auto reference = pg(input, true, arithmetic, plan);
    const auto evidence = geometry(reference, input, crossing, reversed, geometric, plan);
    const auto one = input.r("1");
    const std::array claim{ExternalPgThicknessClaim{17, one, core::Index{4}}};
    std::optional<ExternalPgThicknessDiagnostic> diagnostic;
    const auto result = compare_external_pg_regional_thickness(evidence, reference, claim, limits,
                                                               plan, &diagnostic);
    if (crossing) {
        REQUIRE_FALSE(result);
        CHECK(result.error().code == Code::InvalidInput);
        if (const auto& detail = diagnostic) {
            CHECK(detail->reason == ExternalPgThicknessFailure::Difference);
            CHECK(detail->element == 17);
            CHECK(detail->sidecar_node == 4);
            CHECK(detail->case_triangle == 1);
            CHECK(detail->overlap_vertex);
        } else {
            FAIL("Missing positive-overlap failure");
        }
        return;
    }
    REQUIRE(result);
    REQUIRE(result->matches.size() == 1);
    const auto& verified = result->matches.values().front();
    CHECK(verified.element == 17);
    CHECK(verified.sidecar_node == 4);
    CHECK(verified.winding ==
          (reversed ? ExternalPgWinding::Reversed : ExternalPgWinding::Forward));
    CHECK(std::get<ExactRational>(verified.value_metres) == one);
    CHECK(verified.twice_chart_area == one);
    CHECK(verified.pieces > 0);
    CHECK(verified.vertices >= 3 * verified.pieces);
    CHECK_FALSE(diagnostic);
    CHECK(result->pg == take(reference.identity()));

    auto exact = limits;
    exact.arithmetic.maximum_partition_visits = result->work.partition_visits;
    exact.arithmetic.maximum_scalar_operations = result->work.scalar_operations;
    REQUIRE(compare_external_pg_regional_thickness(evidence, reference, claim, exact, plan));
    --exact.arithmetic.maximum_partition_visits;
    const auto shortage =
        compare_external_pg_regional_thickness(evidence, reference, claim, exact, plan);
    REQUIRE_FALSE(shortage);
    CHECK(shortage.error().code == Code::ResourceExceeded);
    exact = limits;
    exact.case_vertices = 2;
    const auto vertices =
        compare_external_pg_regional_thickness(evidence, reference, claim, exact, plan);
    REQUIRE_FALSE(vertices);
    CHECK(vertices.error().code == Code::ResourceExceeded);
    exact = limits;
    exact.maximum_claims = 0;
    const auto layouts =
        compare_external_pg_regional_thickness(evidence, reference, claim, exact, plan);
    REQUIRE_FALSE(layouts);
    CHECK(layouts.error().code == Code::ResourceExceeded);
    exact = limits;
    exact.arithmetic.cancellation = partition_test::requested_stop_token();
    const auto cancelled = compare_external_pg_regional_thickness(evidence, reference, claim, exact,
                                                                  plan, &diagnostic);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Code::Cancelled);
    CHECK_FALSE(diagnostic);

    const std::array duplicate{claim.front(), claim.front()};
    const auto repeated =
        compare_external_pg_regional_thickness(evidence, reference, duplicate, limits, plan);
    REQUIRE_FALSE(repeated);
    CHECK(repeated.error().code == Code::InvalidInput);
    partition_test::Fixture other;
    const auto foreign = pg(other, false, arithmetic, plan);
    const auto stale =
        compare_external_pg_regional_thickness(evidence, foreign, claim, limits, plan);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == Code::StaleGeneration);
}

[[nodiscard]] ExternalPlanarPg quad_geometry(const PgReferenceGeometry& reference,
                                             const partition_test::Fixture& coordinates, bool skew,
                                             bool reversed, const ExternalPgLimits& limits,
                                             core::ResourcePlan& plan) {
    const auto raw = std::string{"$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 4 1 4\n"
                                 "2 17 0 4\n1\n2\n3\n4\n0 0 0\n1 0 0\n"} +
                     (skew ? "1 0.5 0\n" : "1 1 0\n") +
                     "0 1 0\n$EndNodes\n$Elements\n1 1 19 19\n2 17 3 1\n" +
                     (reversed ? "19 1 4 3 2\n" : "19 1 2 3 4\n") + "$EndElements\n";
    const auto msh = source(raw, plan);
    const std::array<std::uint64_t, 1> tags{19};
    const auto topology =
        take(GmshSurfaceTopology::prepare(msh, tags, {8, 16, 16, 4, 32'000}, plan, {}));
    const auto& pts = coordinates.domain;
    const auto first = reversed ? 3U : 1U;
    const auto last = reversed ? 1U : 3U;
    const std::array ring{ExternalPgCorner{1, pts[0]}, ExternalPgCorner{first + 1, pts[first]},
                          ExternalPgCorner{3, pts[2]}, ExternalPgCorner{last + 1, pts[last]}};
    const auto winding = reversed ? ExternalPgWinding::Reversed : ExternalPgWinding::Forward;
    const std::array faces{ExternalPgFace{19, winding, ring}};
    return take(ExternalPlanarPg::prepare(reference, msh, topology, SiLengthUnit::Metre, faces,
                                          limits, plan));
}

TEST_CASE("Original quadrilateral thickness must agree as a bilinear field",
          "[FD-03][IO4-03][R-05][N4-02][integration][verification][negative]") {
    bool split_would_hide = false;
    bool reversed = false;
    bool wrong_order = false;
    SECTION("A non-parallelogram represents a common affine physical thickness") {}
    SECTION("Reverse winding retains the original node thickness order") {
        reversed = true;
    }
    SECTION("Triangulating would hide an interior mismatch with identical corner values") {
        split_would_hide = true;
    }
    SECTION("A reordered corner claim cannot attach source thickness to the wrong node") {
        wrong_order = true;
    }

    partition_test::Fixture input;
    if (!split_would_hide) {
        input.domain[2] = {input.r("1"), input.r("1", "2")};
        input.heights[1] = input.r("2");
        input.heights[2] = input.r("3");
        input.heights[3] = input.r("2");
    }
    core::ResourcePlan plan{partition_test::cap};
    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const auto reference = pg(input, split_would_hide, arithmetic, plan);
    const auto geometry_record =
        quad_geometry(reference, input, !split_would_hide, reversed, {arithmetic, 4}, plan);
    const std::array thickness{input.r("1"), input.r(split_would_hide ? "1" : "2"),
                               split_would_hide ? input.r("1") : input.r("5", "2"), input.r("2")};
    constexpr std::array<std::uint64_t, 4> direct_order{1, 2, 3, 4};
    constexpr std::array<std::uint64_t, 4> reverse_order{1, 4, 3, 2};
    const auto& order = reversed ? reverse_order : direct_order;
    std::array<ExternalPgCornerThicknessClaim, 4> source_fields{};
    for (std::size_t i = 0; i < order.size(); ++i) {
        source_fields[i] = {order[i], &thickness[static_cast<std::size_t>(order[i] - 1)],
                            NastranBulkPosition{100 + i, 2, 4 + i}};
    }
    if (wrong_order) {
        source_fields[0].node = 2;
    }
    const std::array claims{ExternalPgThicknessClaim{19, source_fields, 27}};
    std::optional<ExternalPgThicknessDiagnostic> diagnostic;
    const auto result = compare_external_pg_regional_thickness(
        geometry_record, reference, claims, {arithmetic, 1, 4}, plan, &diagnostic);
    if (split_would_hide || wrong_order) {
        REQUIRE_FALSE(result);
        CHECK(result.error().code == Code::InvalidInput);
        if (const auto& failure = diagnostic) {
            CHECK(failure->element == 19);
            CHECK(failure->sidecar_node == 27);
            if (wrong_order) {
                CHECK(failure->reason == ExternalPgThicknessFailure::Selection);
                CHECK(failure->node == 2);
                CHECK(failure->field == source_fields.front().field);
            } else {
                // On the lower Case triangle at (3/4,1/4), the Case is 1 and
                // the original bilinear field is 1+(1-3/4)*(1/4)=17/16.
                CHECK(failure->reason == ExternalPgThicknessFailure::Difference);
                CHECK(failure->case_triangle);
                CHECK_FALSE(failure->overlap_vertex);
                CHECK(failure->field);
                CHECK((failure->node == 2 || failure->node == 4));
            }
        } else {
            FAIL("Missing original thickness field diagnostic");
        }
        return;
    }
    REQUIRE(result);
    CHECK_FALSE(diagnostic);
    REQUIRE(result->matches.size() == 1);
    const auto& face = result->matches.values().front();
    CHECK(face.twice_chart_area == input.r("3", "2"));
    const auto* nodal = std::get_if<ExternalPgNodalThickness>(&face.value_metres);
    REQUIRE(nodal != nullptr);
    CHECK(nodal->count == 4);
    for (std::size_t i = 0; i < order.size(); ++i) {
        CHECK(nodal->values[i].node == order[i]);
        CHECK(nodal->values[i].value_metres == *source_fields[i].value_metres);
        CHECK(nodal->values[i].field == source_fields[i].field);
    }
}

[[nodiscard]] ExternalPlanarPg two_quads(const PgReferenceGeometry& reference,
                                         partition_test::Fixture& inputs,
                                         const PgReferenceLimits& arithmetic,
                                         core::ResourcePlan& plan) {
    constexpr std::string_view mesh = R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 6 1 6
2 17 0 6
1 2 3 4 5 6
0 0 0
0.5 0 0
1 0 0
0 1 0
0.5 1 0
1 1 0
$EndNodes
$Elements
1 2 19 20
2 17 3 2
19 1 2 5 4
20 2 3 6 5
$EndElements
)msh";
    const auto msh = source(mesh, plan);
    constexpr std::array<std::uint64_t, 2> tags{19, 20};
    const auto topology =
        take(GmshSurfaceTopology::prepare(msh, tags, {8, 16, 16, 4, 32'000}, plan, {}));
    const auto& corners = inputs.domain;
    const PgPartitionPoint bottom{inputs.r("1", "2"), inputs.r("0")};
    const PgPartitionPoint top{inputs.r("1", "2"), inputs.r("1")};
    const std::array<ExternalPgCorner, 4> left{
        {{1, corners[0]}, {2, bottom}, {5, top}, {4, corners[3]}}};
    const std::array<ExternalPgCorner, 4> right{
        {{2, bottom}, {3, corners[1]}, {6, corners[2]}, {5, top}}};
    const std::array faces{ExternalPgFace{19, ExternalPgWinding::Forward, left},
                           ExternalPgFace{20, ExternalPgWinding::Forward, right}};
    return take(ExternalPlanarPg::prepare(reference, msh, topology, SiLengthUnit::Metre, faces,
                                          {arithmetic, 2}, plan));
}

[[nodiscard]] std::uint64_t required_estimate(const ExactLimits& limits, std::uint64_t retained) {
    using Limb = boost::multiprecision::limb_type;
    constexpr auto width = std::numeric_limits<Limb>::digits;
    const auto limbs = limits.maximum_work_bits / width;
    const auto live = structure::detail::reference_live_rationals + retained;
    // A held rational consists of two integer components. The algebra also
    // reserves an integer slot for add/sub and a separate leaf workspace.
    const auto held = (2 * live + 1) * 4 * limbs * sizeof(Limb);
    return take(structure::detail::exact_rational_admission(limits, (limbs - 2) * width, held))
        .estimated_external_payload_bytes;
}

TEST_CASE("Staged thickness memory admits the whole face batch before generation",
          "[FD-03][IO4-03][IO4-04][N4-02][integration][verification][negative]") {
    core::ResourcePlan plan{partition_test::cap};
    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    auto limited = arithmetic;
    // One homogeneous PSHELL region plus the four original corners of one face.
    const auto budget = required_estimate(arithmetic.algebra, 1 + 4);
    limited.algebra.external_memory = ExactEstimateAdmission{budget};
    const auto admitted = admit_external_pg_thickness_inputs(limited, 1, 1, plan);
    REQUIRE(admitted);
    CHECK(std::get<ExactEstimateAdmission>(admitted->external_memory).bytes < budget);
    const auto before = plan.allocated_bytes();
    const auto excess = admit_external_pg_thickness_inputs(limited, 1, 2, plan);
    REQUIRE_FALSE(excess);
    CHECK(excess.error().code == Code::ResourceExceeded);
    CHECK(plan.allocated_bytes() == before);
    limited.algebra.external_memory = ExactEstimateAdmission{budget - 1};
    const auto shortage = admit_external_pg_thickness_inputs(limited, 1, 1, plan);
    REQUIRE_FALSE(shortage);
    CHECK(shortage.error().code == Code::ResourceExceeded);
    CHECK(plan.allocated_bytes() == before);
}

TEST_CASE("Whole-face thickness comparison charges inputs and reported copies together",
          "[FD-03][IO4-03][IO4-04][N4-02][integration][verification][negative]") {
    partition_test::Fixture input;
    core::ResourcePlan plan{partition_test::cap};
    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const auto reference = pg(input, false, arithmetic, plan);
    const auto evidence = two_quads(reference, input, arithmetic, plan);
    const std::array values{input.r("1"), input.r("1"), input.r("1"), input.r("1"),
                            input.r("1"), input.r("1"), input.r("1"), input.r("1")};
    const std::array<ExternalPgCornerThicknessClaim, 4> left{
        {{1, &values[0]}, {2, &values[1]}, {5, &values[2]}, {4, &values[3]}}};
    const std::array<ExternalPgCornerThicknessClaim, 4> right{
        {{2, &values[4]}, {3, &values[5]}, {6, &values[6]}, {5, &values[7]}}};
    const std::array claims{ExternalPgThicknessClaim{19, left},
                            ExternalPgThicknessClaim{20, right}};
    const auto first = std::span{claims}.first(1);
    // Four staged inputs and four independent copies coexist with the chart-area result.
    const auto budget = required_estimate(arithmetic.algebra, 4 + 4 + 1);
    auto limited = arithmetic;
    limited.algebra.external_memory = ExactEstimateAdmission{budget};
    const auto admitted =
        compare_external_pg_regional_thickness(evidence, reference, first, {limited, 2, 4}, plan);
    REQUIRE(admitted);
    REQUIRE(admitted->matches.size() == 1);
    CHECK(
        std::get<ExternalPgNodalThickness>(admitted->matches.values().front().value_metres).count ==
        4);
    const auto before = plan.allocated_bytes();
    const auto excess =
        compare_external_pg_regional_thickness(evidence, reference, claims, {limited, 2, 4}, plan);
    REQUIRE_FALSE(excess);
    CHECK(excess.error().code == Code::ResourceExceeded);
    CHECK(plan.allocated_bytes() == before);
    limited.algebra.external_memory = ExactEstimateAdmission{budget - 1};
    const auto shortage =
        compare_external_pg_regional_thickness(evidence, reference, first, {limited, 2, 4}, plan);
    REQUIRE_FALSE(shortage);
    CHECK(shortage.error().code == Code::ResourceExceeded);
    CHECK(plan.allocated_bytes() == before);
    const auto again = compare_external_pg_regional_thickness(evidence, reference, claims,
                                                              {arithmetic, 2, 4}, plan);
    REQUIRE(again);
    CHECK(again->matches.size() == 2);
}
} // namespace
