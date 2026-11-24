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
    CHECK(verified.value_metres == one);
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

} // namespace
