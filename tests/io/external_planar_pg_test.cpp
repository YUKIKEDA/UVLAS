#include "../structure/pg_partition_test_support.hpp"
#include "external_planar_pg.hpp"
#include "source_digest.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using namespace uvlas::structure;
using partition_test::Fixture;
using partition_test::take;

constexpr std::uint64_t budget = std::uint64_t{1} << 33U;
using Winding = ExternalPgWinding;
using Reason = ExternalPgFailure;
using Code = core::FailureCode;

template <class T>
T take(std::optional<T> value) {
    if (!value) {
        throw std::runtime_error("Missing validated input view");
    }
    return std::move(*value);
}

[[nodiscard]] ExactLimits algebra() noexcept {
    return {
        256, 256, 256, 8, 8, 4096, 8192, ExactEstimateAdmission{budget}, ExactTimingMode::Disabled};
}

[[nodiscard]] PgReferenceLimits reference_limits() noexcept {
    auto exact = algebra();
    exact.maximum_value_bits = 8192;
    exact.maximum_work_bits = 32768;
    return {exact, 20'000'000, 10'000'000, {}};
}

class Reader final : public ReferenceFileReader {
public:
    explicit Reader(std::span<const std::byte> bytes) noexcept : bytes_(bytes) {}

    core::Result<ReferenceFileInfo> info() noexcept override {
        return ReferenceFileInfo{{0, {}}, bytes_.size(), {}, true};
    }

    core::Result<std::size_t> read(std::span<std::byte> out) noexcept override {
        const auto count = std::min(bytes_.size() - offset_, out.size());
        std::ranges::copy(bytes_.subspan(offset_, count), out.begin());
        offset_ += count;
        return count;
    }

private:
    std::span<const std::byte> bytes_;
    std::size_t offset_{};
};

[[nodiscard]] AcquiredRootSnapshot snapshot(std::string_view text, SourceId source,
                                            ReferenceRole role, core::ResourcePlan& plan) {
    const auto raw = std::as_bytes(std::span{text});
    Reader reader{raw};
    auto file = take(read_reference_file(reader, source_digest(raw),
                                         {text.size(), text.size() + 1, 2, text.size()}, plan, {},
                                         nullptr));
    return {{{{SourceId{100}, 20, {}}, 1, 1}, source, role, 10}, std::move(file), {}};
}

[[nodiscard]] std::string mesh_text(std::string_view first_coordinate = "1500") {
    return "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 4 10 95\n2 17 0 4\n10\n40\n90\n95\n" +
           std::string{first_coordinate} +
           " 2500 3000\n1750 6500 3000\n1625 6250 4000\n1250 2000 5000\n"
           "$EndNodes\n$Elements\n1 1 701 701\n2 17 3 1\n701 10 95 90 40\n$EndElements\n";
}

[[nodiscard]] GmshMesh mesh(std::string_view input, core::ResourcePlan& plan) {
    const auto source = snapshot(input, SourceId{17}, ReferenceRole::GmshMesh, plan);
    return take(GmshMesh::prepare(source, {4096, 0, 0, 0, 12, 12, 64000}, plan, {}));
}

[[nodiscard]] NastranShellCards cards(std::string_view connectivity, core::ResourcePlan& plan,
                                      bool nonbasic = false) {
    std::string text = nonbasic ? "CORD2R,5,,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n"
                                  "GRID,10,5,1500.0,2500.0,3000.0\n"
                                : "GRID,10,,1500.0,2500.0,3000.0\n";
    text += "GRID,40,,1750.0,6500.0,3000.0\nGRID,90,,1625.0,6250.0,4000.0\n"
            "GRID,95,,1250.0,2000.0,5000.0\nPSHELL,20\n";
    text += connectivity;
    auto file = snapshot(text, SourceId{20}, ReferenceRole::NastranDeck, plan);
    auto bulk =
        take(NastranBulkDeck::prepare(std::move(file), {4096, 80, 32, 32, 256, 64000}, plan, {}));
    auto coordinate =
        take(NastranCoordinateCards::prepare(std::move(bulk), {12, 8, 8, 20000}, plan, {}));
    return take(NastranShellCards::prepare(std::move(coordinate), {12, 12, 20000}, plan, {}));
}

[[nodiscard]] PgReferenceGeometry pg(Fixture& f, core::ResourcePlan& plan, bool curved = false,
                                     bool triangle = false) {
    INFO("Preparing the target PG chart, curved=" << curved << " triangle=" << triangle);
    f.domain[2] = {f.r("1"), f.r("1", "2")};
    if (triangle) {
        f.domain.resize(3);
    }
    const auto quarter = f.r("1", "4"), width = f.r("2"), skew = f.r("-1", "4");
    const auto fs = curved ? f.r("1249", "10000") : f.r("1", "8");
    const auto ft = f.r("-1", "8"), curvature = f.r("1", "10000");
    const std::array terms{ExactInputTerm{{0, 0}, quarter}, ExactInputTerm{{1, 0}, fs},
                           ExactInputTerm{{0, 1}, ft}, ExactInputTerm{{2, 0}, curvature}};
    const auto height =
        take(ExactChartPolynomial::make(std::span{terms}.first(curved ? 4 : 3), algebra(), plan));
    f.graph = take(PgGraph::prepare({quarter, width, skew, height, 1, PgGraphConvention::Version1},
                                    {algebra(), 64, {}}, plan));
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto dictionary = take(PgCoefficientDictionary::prepare(
        f.chart_graph(), {partition_test::exact_limits(), 2000, 20000000, 100000, {}}, workspace,
        plan));
    auto partition = take(
        PgIntegrationPartition::prepare(f.request(), partition_test::limits(), workspace, plan));
    const auto length = f.r("2");
    const std::array origin{f.r("1"), f.r("2"), f.r("3")};
    std::array<ExactRational, 9> rotation;
    rotation[2] = f.r("1");
    rotation[3] = f.r("1");
    rotation[7] = f.r("1");
    PgReferenceInputs input{std::move(f.chart_graph()), std::move(dictionary),
                            std::move(partition)};
    return take(PgReferenceGeometry::prepare(std::move(input), {rotation, origin, length},
                                             reference_limits(), plan));
}

struct Corners {
    Corners(Fixture& f)
        : a(f.point("0", "0")),
          b(f.point("1", "0")),
          c{f.r("1"), f.r("1", "2")},
          d(f.point("0", "1")) {}
    PgPartitionPoint a, b, c, d;
};

void rejected(const core::Result<ExternalPlanarPg>& result, Code code,
              const std::optional<ExternalPgDiagnostic>& diagnostic, Reason reason) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    CHECK(take(diagnostic).reason == reason);
}

TEST_CASE("Gmsh and Nastran preserve differently oriented source quads on a rotated SI PG plane",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][verification]") {
    const auto result = [] {
        Fixture fixture;
        core::ResourcePlan plan{budget};
        const auto target = pg(fixture, plan);
        const auto identity = take(target.identity());
        const Corners xy{fixture};
        const std::array from_gmsh{ExternalPgCorner{10, xy.a}, ExternalPgCorner{95, xy.d},
                                   ExternalPgCorner{90, xy.c}, ExternalPgCorner{40, xy.b}};
        const std::array from_nastran{ExternalPgCorner{10, xy.a}, ExternalPgCorner{40, xy.b},
                                      ExternalPgCorner{90, xy.c}, ExternalPgCorner{95, xy.d}};
        const std::array gmsh_map{ExternalPgFace{701, Winding::Reversed, from_gmsh}};
        const std::array nastran_map{ExternalPgFace{101, Winding::Forward, from_nastran}};
        const auto msh = mesh(mesh_text(), plan);
        const auto bdf = cards("CQUAD4,101,20,10,40,90,95\n", plan);
        constexpr std::array<std::uint64_t, 1> selected_tags{701};
        const auto topology =
            take(GmshSurfaceTopology::prepare(msh, selected_tags, {4, 16, 16, 4, 4096}, plan, {}));
        const auto source_gmsh = take(msh.view());
        const auto source_bdf = take(bdf.view()).input.bulk;
        const ExternalPgLimits limits{reference_limits(), 2};
        auto gmsh = take(ExternalPlanarPg::prepare(target, msh, topology, SiLengthUnit::Millimetre,
                                                   gmsh_map, limits, plan));
        auto nastran = take(ExternalPlanarPg::prepare(target, bdf, SiLengthUnit::Millimetre,
                                                      nastran_map, limits, plan));
        CHECK(take(gmsh.view()).patch == identity);
        CHECK(take(gmsh.view()).digest == source_gmsh.digest);
        CHECK(take(nastran.view()).digest == source_bdf.digest);
        return std::pair{std::move(gmsh), std::move(nastran)};
    }();

    const auto gmsh = take(result.first.view());
    const auto nastran = take(result.second.view());
    REQUIRE(gmsh.faces.size() == 1);
    REQUIRE(nastran.faces.size() == 1);
    const auto& g = gmsh.faces.front();
    const auto& n = nastran.faces.front();
    REQUIRE(g.ring().size() == 4);
    REQUIRE(n.ring().size() == 4);
    CHECK(g.winding == Winding::Reversed);
    CHECK(n.winding == Winding::Forward);
    CHECK(g.element == 701);
    CHECK(n.element == 101);
    CHECK(gmsh.units == SiLengthUnit::Millimetre);
    CHECK(nastran.patch == gmsh.patch);
    CHECK(g.corners[1].node == 95);
    CHECK(n.corners[1].node == 40);
    CHECK(g.corners[1].position_metres == n.corners[3].position_metres);
    CHECK(n.corners[2].source_bytes[1].has_value());
    Fixture expected;
    const auto& c = n.corners[2];
    // X=1.5+s/4-t/4, Y=2.5+4s-t/2 and Z=3+2t.
    CHECK(c.chart.s == expected.r("1"));
    CHECK(c.chart.t == expected.r("1", "2"));
    CHECK(c.position_metres[0] == expected.r("13", "8"));
    CHECK(c.position_metres[1] == expected.r("25", "4"));
    CHECK(c.position_metres[2] == expected.r("4"));
}

TEST_CASE("Planar certification separates winding, coverage and overlap without a geometry fit",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][negative][verification]") {
    Fixture fixture, alternate, subset;
    core::ResourcePlan plan{budget};
    const auto target = pg(fixture, plan);
    const auto curved = pg(alternate, plan, true);
    const auto clipped = pg(subset, plan, false, true);
    const Corners xy{fixture};
    const auto deck = cards("CQUAD4,101,20,10,40,90,95\n"
                            "CTRIA3,201,20,10,40,90\nCTRIA3,202,20,10,40,95\n"
                            "CTRIA3,203,20,10,90,95\nCQUAD4,301,20,10,90,40,95\n",
                            plan);
    const auto nonbasic = cards("CQUAD4,101,20,10,40,90,95\n", plan, true);
    const auto offset = cards("CQUAD4,101,20,10,40,90,95,,1.0\n", plan);
    const std::array abcd{ExternalPgCorner{10, xy.a}, ExternalPgCorner{40, xy.b},
                          ExternalPgCorner{90, xy.c}, ExternalPgCorner{95, xy.d}};
    const std::array acbd{abcd[0], abcd[2], abcd[1], abcd[3]};
    const std::array abd{abcd[0], abcd[1], abcd[3]};
    const std::array acd{abcd[0], abcd[2], abcd[3]};
    const std::array full{ExternalPgFace{101, Winding::Forward, abcd}};
    const std::array reversed{ExternalPgFace{101, Winding::Reversed, abcd}};
    const std::array folded{ExternalPgFace{301, Winding::Forward, acbd}};
    const auto abc = std::span{abcd}.first(3);
    const std::array gap{ExternalPgFace{201, Winding::Forward, abc}};
    const std::array overlap{gap[0], ExternalPgFace{202, Winding::Forward, abd}};
    const std::array joined{gap[0], ExternalPgFace{203, Winding::Forward, acd}};
    const ExternalPgLimits limits{reference_limits(), 2};
    std::optional<ExternalPgDiagnostic> diagnostic;
    const auto qualify = [&](const PgReferenceGeometry& patch, const NastranShellCards& cards,
                             std::span<const ExternalPgFace> charts) {
        return ExternalPlanarPg::prepare(patch, cards, SiLengthUnit::Millimetre, charts, limits,
                                         plan, &diagnostic);
    };
    // The extra term in the curved graph is (s*s-s)/10000, zero at all four corners.
    rejected(qualify(curved, deck, full), Code::Unsupported, diagnostic, Reason::CurvedPatch);
    rejected(qualify(target, deck, reversed), Code::GeometryInvalid, diagnostic, Reason::Winding);
    rejected(qualify(target, deck, folded), Code::GeometryInvalid, diagnostic, Reason::Winding);
    rejected(qualify(clipped, deck, full), Code::GeometryInvalid, diagnostic, Reason::Outside);
    rejected(qualify(target, deck, gap), Code::GeometryInvalid, diagnostic, Reason::Coverage);
    rejected(qualify(target, deck, overlap), Code::GeometryInvalid, diagnostic, Reason::Overlap);
    CHECK(take(diagnostic).other == 201);
    rejected(qualify(target, nonbasic, full), Code::Unsupported, diagnostic,
             Reason::PositionSystem);
    CHECK(take(diagnostic).node == 10);
    rejected(qualify(target, offset, full), Code::Unsupported, diagnostic, Reason::Offset);
    CHECK(take(diagnostic).byte.has_value());
    REQUIRE(qualify(target, deck, joined));

    constexpr std::array<std::uint64_t, 1> selection{701};
    const auto msh = mesh(mesh_text(), plan);
    const auto changed = mesh(mesh_text("1500.0000000000000000000001"), plan);
    const auto selected =
        take(GmshSurfaceTopology::prepare(msh, selection, {4, 16, 16, 4, 4096}, plan, {}));
    const auto changed_topology =
        take(GmshSurfaceTopology::prepare(changed, selection, {4, 16, 16, 4, 4096}, plan, {}));
    const std::array ga{abcd[0], abcd[3], abcd[2], abcd[1]};
    const std::array mapping{ExternalPgFace{701, Winding::Reversed, ga}};
    const std::array order_error{ExternalPgFace{701, Winding::Reversed, abcd}};
    rejected(ExternalPlanarPg::prepare(target, changed, selected, SiLengthUnit::Millimetre, mapping,
                                       limits, plan, &diagnostic),
             Code::StaleGeneration, diagnostic, Reason::Source);
    rejected(ExternalPlanarPg::prepare(target, changed, changed_topology, SiLengthUnit::Millimetre,
                                       mapping, limits, plan, &diagnostic),
             Code::GeometryInvalid, diagnostic, Reason::Position);
    CHECK(take(diagnostic).byte.has_value());
    rejected(ExternalPlanarPg::prepare(target, msh, selected, SiLengthUnit::Millimetre, order_error,
                                       limits, plan, &diagnostic),
             Code::InvalidInput, diagnostic, Reason::SourceOrder);

    auto short_limits = limits;
    short_limits.maximum_faces = 0;
    const auto limit = ExternalPlanarPg::prepare(target, deck, SiLengthUnit::Millimetre, full,
                                                 short_limits, plan, &diagnostic);
    REQUIRE_FALSE(limit);
    CHECK(limit.error().code == Code::ResourceExceeded);
    short_limits.arithmetic.cancellation = partition_test::requested_stop_token();
    const auto cancelled = ExternalPlanarPg::prepare(target, deck, SiLengthUnit::Millimetre, full,
                                                     short_limits, plan, &diagnostic);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Code::Cancelled);
    CHECK_FALSE(diagnostic);
    REQUIRE(qualify(target, deck, full));
}

} // namespace
