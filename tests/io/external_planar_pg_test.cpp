#include "../structure/pg_partition_test_support.hpp"
#include "external_planar_pg.hpp"
#include "source_digest.hpp"

#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <compare>
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

[[nodiscard]] NastranShellCards read_cards(std::string_view text, core::ResourcePlan& plan) {
    auto file = snapshot(text, SourceId{20}, ReferenceRole::NastranDeck, plan);
    auto bulk =
        take(NastranBulkDeck::prepare(std::move(file), {4096, 80, 32, 32, 256, 64000}, plan, {}));
    auto coordinate =
        take(NastranCoordinateCards::prepare(std::move(bulk), {12, 8, 8, 20000}, plan, {}));
    return take(NastranShellCards::prepare(std::move(coordinate), {12, 12, 20000}, plan, {}));
}

[[nodiscard]] NastranShellCards cards(std::string_view connectivity, core::ResourcePlan& plan,
                                      bool nonbasic = false) {
    std::string text = nonbasic ? "CORD2R,5,,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n"
                                  "GRID,10,5,1500.0,2500.0,3000.0\n"
                                : "GRID,10,,1500.0,2500.0,3000.0\n";
    text += "GRID,40,,1750.0,6500.0,3000.0\nGRID,90,,1625.0,6250.0,4000.0\n"
            "GRID,95,,1250.0,2000.0,5000.0\nPSHELL,20\n";
    text += connectivity;
    return read_cards(text, plan);
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
    Fixture fixture, alternate, subset, notched;
    core::ResourcePlan plan{budget};
    const auto target = pg(fixture, plan);
    const auto curved = pg(alternate, plan, true);
    const auto clipped = pg(subset, plan, false, true);
    notched.domain.insert(notched.domain.begin() + 3,
                          PgPartitionPoint{notched.r("1", "2"), notched.r("1", "3")});
    const auto concave = pg(notched, plan);
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
    // Every original corner belongs to this patch, but the edge CD crosses the notch.
    rejected(qualify(concave, deck, full), Code::GeometryInvalid, diagnostic, Reason::Outside);
    rejected(qualify(target, deck, gap), Code::GeometryInvalid, diagnostic, Reason::Coverage);
    rejected(qualify(target, deck, overlap), Code::GeometryInvalid, diagnostic, Reason::Overlap);
    CHECK(take(diagnostic).other == 201);
    REQUIRE(qualify(target, nonbasic, full));
    rejected(qualify(target, offset, full), Code::PrecisionUnsupported, diagnostic, Reason::Offset);
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

PgReferenceGeometry flat_target(Fixture& f, const PgReferenceFrame& frame,
                                core::ResourcePlan& plan) {
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto dictionary = take(PgCoefficientDictionary::prepare(
        f.chart_graph(), {partition_test::exact_limits(), 2000, 20000000, 100000, {}}, workspace,
        plan));
    auto partition = take(
        PgIntegrationPartition::prepare(f.request(), partition_test::limits(), workspace, plan));
    PgReferenceInputs input{std::move(f.chart_graph()), std::move(dictionary),
                            std::move(partition)};
    return take(PgReferenceGeometry::prepare(std::move(input), frame, reference_limits(), plan));
}

TEST_CASE("Nastran rectangular ancestors and signed shell offsets match one exact PG plane",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][negative][verification]") {
    Fixture f;
    core::ResourcePlan plan{budget};
    const auto length = f.r("1");
    const std::array origin{f.r("61", "50"), f.r("-27", "50"), f.r("1")};
    const std::array rotation{f.r("4", "5"),  f.r("-3", "5"), f.r("0"),
                              f.r("-3", "5"), f.r("-4", "5"), f.r("0"),
                              f.r("0"),       f.r("0"),       f.r("-1")};
    const auto target = flat_target(f, {rotation, origin, length}, plan);
    const std::array points{f.point("0", "0"), f.point("1", "0"), f.point("1", "1"),
                            f.point("0", "1")};
    const std::array abcd{ExternalPgCorner{10, points[0]}, ExternalPgCorner{40, points[1]},
                          ExternalPgCorner{90, points[2]}, ExternalPgCorner{95, points[3]}};
    const std::array dcba{abcd[0], abcd[3], abcd[2], abcd[1]};
    const std::array acd{abcd[0], abcd[2], abcd[3]};
    const std::array forward{ExternalPgFace{201, Winding::Forward, abcd}};
    const std::array reversed{ExternalPgFace{202, Winding::Reversed, dcba}};
    const std::array wrong_sign{ExternalPgFace{203, Winding::Forward, abcd}};
    const std::array triangles{ExternalPgFace{301, Winding::Forward, std::span{abcd}.first(3)},
                               ExternalPgFace{302, Winding::Forward, acd}};
    constexpr std::string_view nested =
        "CORD2R,5,12,100.0,200.0,300.0,100.0,200.0,400.0,100.0,300.0,300.0\n"
        "CORD2R,12,0,1000.0,-500.0,1000.0,1000.0,-500.0,900.0,1000.3,-499.6,1000.0\n";
    constexpr std::string_view nodes = "GRID,10,5,0.0,0.0,0.0\nGRID,40,5,1000.0,0.0,0.0\n"
                                       "GRID,90,5,1000.0,1000.0,0.0\nGRID,95,5,0.0,1000.0,0.0\n"
                                       "PSHELL,20,8,2.54,8,,8\nMAT1,8,1.0\n";
    constexpr std::string_view shells =
        "CQUAD4,201,20,10,40,90,95,,-300.0\nCQUAD4,202,20,10,95,90,40,,300.0\n"
        "CQUAD4,203,20,10,40,90,95,,300.0\n"
        "CTRIA3,301,20,10,40,90,,-300.0\nCTRIA3,302,20,10,90,95,,-300.0\n";
    const auto deck =
        read_cards(std::string{nested} + std::string{nodes} + std::string{shells}, plan);
    const auto view = take(deck.view());
    const ExternalPgLimits limits{reference_limits(), 2};
    std::optional<ExternalPgDiagnostic> diagnostic;
    const auto qualify = [&](const NastranShellCards& input,
                             std::span<const ExternalPgFace> mapping) {
        return ExternalPlanarPg::prepare(target, input, SiLengthUnit::Millimetre, mapping, limits,
                                         plan, &diagnostic);
    };

    const auto forward_match = take(qualify(deck, forward));
    const auto forward_view = take(forward_match.view());
    const auto& face = forward_view.faces.front();
    const auto signed_offset =
        face.nastran_offset
            ? std::reference_wrapper<const NastranReferenceOffset>{*face.nastran_offset}
            : throw std::runtime_error("Missing ZOFFS source");
    CHECK(signed_offset.get().signed_metres == f.r("-3", "10"));
    CHECK(signed_offset.get().byte == take(view.elements[0].offset).field.position.byte);
    const auto& node = face.ring()[1];
    if (const auto& placement = node.nastran) {
        CHECK(placement->position_system == NastranCoordinateId{5});
        CHECK(placement->cp_byte.has_value());
        if (const auto& grid = placement->unshifted_metres) {
            // In world millimetres, (x,y,z) -> (1220+4x/5-3y/5,-540-3x/5-4y/5,700-z).
            CHECK((*grid)[0] == f.r("101", "50"));
            CHECK((*grid)[1] == f.r("-57", "50"));
            CHECK((*grid)[2] == f.r("7", "10"));
        } else {
            FAIL("Original basic GRID coordinates were lost");
        }
    } else {
        FAIL("GRID CP provenance was lost");
    }
    CHECK(node.position_metres[0] == f.r("101", "50"));
    CHECK(node.position_metres[1] == f.r("-57", "50"));
    CHECK(node.position_metres[2] == f.r("1"));

    const auto reverse_match = take(qualify(deck, reversed));
    const auto reverse_view = take(reverse_match.view());
    const auto& reverse_face = reverse_view.faces.front();
    if (const auto& source = reverse_face.nastran_offset) {
        CHECK(source->signed_metres == f.r("3", "10"));
    } else {
        FAIL("Reversed element lost its ZOFFS declaration");
    }
    CHECK(reverse_face.ring()[3].position_metres == node.position_metres);
    REQUIRE(qualify(deck, triangles));
    rejected(qualify(deck, wrong_sign), Code::GeometryInvalid, diagnostic, Reason::Position);
    const auto sign_failure = take(diagnostic);
    CHECK(sign_failure.byte == take(view.elements[2].offset).field.position.byte);
    CHECK(sign_failure.coordinate_system == NastranCoordinateId{5});

    auto warped_nodes = std::string{nodes};
    constexpr std::string_view flat_grid = "GRID,95,5,0.0,1000.0,0.0";
    const auto grid_byte = warped_nodes.find(flat_grid);
    REQUIRE(grid_byte != std::string::npos);
    warped_nodes.replace(grid_byte, flat_grid.size(), "GRID,95,5,0.0,1000.0,100.0");
    const auto warped = read_cards(std::string{nested} + warped_nodes + std::string{shells}, plan);
    rejected(qualify(warped, forward), Code::GeometryInvalid, diagnostic, Reason::Offset);
    CHECK(take(diagnostic).node == 95);

    const auto abnormal = [&](std::string_view frame) {
        return read_cards(std::string{frame} + std::string{nodes} + std::string{shells}, plan);
    };
    const auto irrational = abnormal("CORD2R,5,0,0.0,0.0,0.0,1.0,0.0,1.0,1.0,1.0,0.0\n");
    rejected(qualify(irrational, forward), Code::PrecisionUnsupported, diagnostic,
             Reason::CoordinateFrame);
    CHECK(take(diagnostic).coordinate_system == NastranCoordinateId{5});
    CHECK(take(diagnostic).node == 10);
    CHECK(take(diagnostic).byte.has_value());

    const auto invalid = abnormal("CORD2R,5,0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,2.0\n");
    rejected(qualify(invalid, forward), Code::GeometryInvalid, diagnostic, Reason::CoordinateFrame);
    const auto curvilinear = abnormal("CORD2C,5,0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n");
    rejected(qualify(curvilinear, forward), Code::Unsupported, diagnostic, Reason::PositionSystem);
}

TEST_CASE("Nastran ZOFFS moves both components of an oblique planar normal",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][verification]") {
    Fixture f;
    core::ResourcePlan plan{budget};
    const auto length = f.r("1");
    const std::array origin{f.r("0"), f.r("3", "10"), f.r("2", "5")};
    const std::array rotation{f.r("1"),      f.r("0"), f.r("0"),       f.r("0"),     f.r("4", "5"),
                              f.r("3", "5"), f.r("0"), f.r("-3", "5"), f.r("4", "5")};
    const auto target = flat_target(f, {rotation, origin, length}, plan);
    const std::array points{f.point("0", "0"), f.point("1", "0"), f.point("1", "1"),
                            f.point("0", "1")};
    const std::array ring{ExternalPgCorner{10, points[0]}, ExternalPgCorner{40, points[1]},
                          ExternalPgCorner{90, points[2]}, ExternalPgCorner{95, points[3]}};
    const std::array selection{ExternalPgFace{201, Winding::Forward, ring}};
    const auto deck =
        read_cards("GRID,10,0,0.0,0.0,0.0\nGRID,40,0,1000.0,0.0,0.0\n"
                   "GRID,90,0,1000.0,800.0,-600.0\nGRID,95,0,0.0,800.0,-600.0\n"
                   "PSHELL,20,8,2.54,8,,8\nMAT1,8,1.0\nCQUAD4,201,20,10,40,90,95,,500.0\n",
                   plan);
    const ExternalPgLimits limits{reference_limits(), 1};
    const auto result = take(
        ExternalPlanarPg::prepare(target, deck, SiLengthUnit::Millimetre, selection, limits, plan));
    const auto face = take(result.view());
    const auto& point = face.faces.front().ring()[2];
    // A source point at (1, .8, -.6) plus .5*(0, .6, .8).
    CHECK(point.position_metres[0] == f.r("1"));
    CHECK(point.position_metres[1] == f.r("11", "10"));
    CHECK(point.position_metres[2] == f.r("-1", "5"));
}

TEST_CASE("Irrational Nastran shell placement bounds the whole oriented planar PG map",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][negative][verification]") {
    Fixture f;
    core::ResourcePlan plan{budget};
    const auto length = f.r("1");
    const std::array<ExactRational, 3> origin;
    // The PG map is (s,t,0) with a +Z normal.
    const std::array rotation{f.r("1"), f.r("0"), f.r("0"), f.r("0"), f.r("1"),
                              f.r("0"), f.r("0"), f.r("0"), f.r("1")};
    const auto target = flat_target(f, {rotation, origin, length}, plan);
    const auto a = f.point("0", "0"), b = f.point("1", "0"), c = f.point("1", "1"),
               d = f.point("0", "1");
    const std::array ring{ExternalPgCorner{10, a}, ExternalPgCorner{40, b}, ExternalPgCorner{90, c},
                          ExternalPgCorner{95, d}};
    const std::array reverse{ring[0], ring[3], ring[2], ring[1]};
    const std::array second_triangle{ring[0], ring[2], ring[3]};
    const std::array forward{ExternalPgFace{101, Winding::Forward, ring}};
    const std::array reversed{ExternalPgFace{102, Winding::Reversed, reverse}};
    const std::array opposite{ExternalPgFace{103, Winding::Reversed, reverse}};
    const std::array triangles{ExternalPgFace{201, Winding::Forward, std::span{ring}.first(3)},
                               ExternalPgFace{202, Winding::Forward, second_triangle}};
    constexpr std::string_view source =
        "CORD2R,5,0,0.0,0.0,0.0,1.0,0.0,1.0,1.0,0.0,-1.0\n"
        "GRID,10,5,0.0,0.0,0.0\nGRID,40,5,1000.0,0.0,0.0\n"
        "GRID,90,5,1000.0,1000.0,0.0\nGRID,95,5,0.0,1000.0,0.0\nPSHELL,20\n"
        "CQUAD4,101,20,10,40,90,95,,500.0\nCQUAD4,102,20,10,95,90,40,,-500.0\n"
        "CQUAD4,103,20,10,95,90,40,,500.0\nCTRIA3,201,20,10,40,90,,500.0\n"
        "CTRIA3,202,20,10,90,95,,500.0\n";
    const ExternalPgBoundLimits limits{{reference_limits(), 2}, 200'000};
    std::optional<ExternalPgDiagnostic> diagnostic;
    const auto measure = [&](const NastranShellCards& deck,
                             std::span<const ExternalPgFace> selection) {
        return NastranPlanarPgBounds::prepare(target, deck, SiLengthUnit::Millimetre, selection,
                                              limits, plan, &diagnostic);
    };
    const auto recorded = [&] {
        const auto deck = read_cards(source, plan);
        const auto exact = ExternalPlanarPg::prepare(target, deck, SiLengthUnit::Millimetre,
                                                     forward, limits.region, plan, &diagnostic);
        rejected(exact, Code::PrecisionUnsupported, diagnostic, Reason::CoordinateFrame);
        auto first = take(measure(deck, forward));
        const auto inspected = take(first.view());
        const auto parsed = take(deck.view());
        CHECK(inspected.digest == parsed.input.bulk.digest);
        CHECK(inspected.patch == take(target.identity()));
        CHECK(inspected.faces.front().element_field == parsed.elements[0].id_field.position);
        CHECK(inspected.directed_work > 0);
        const auto result_rev = take(measure(deck, reversed));
        const auto reverse_view = take(result_rev.view());
        REQUIRE(reverse_view.faces.front().nastran_offset);
        if (const auto& reverse_offset = reverse_view.faces.front().nastran_offset) {
            CHECK(reverse_offset->signed_metres == f.r("-1", "2"));
        }
        CHECK(take(reverse_view.faces.front().corners[3].cp_field).line == 3);
        CHECK(reverse_view.faces.front().corners[3].cp == NastranCoordinateId{5});

        const auto other = take(measure(deck, opposite));
        const auto other_view = take(other.view());
        const auto& large = other_view.faces.front().distance_upper_metres;
        CHECK(std::is_gt(take(large.compare(f.r("61", "50"), algebra(), plan))));
        CHECK(std::is_lt(take(large.compare(f.r("5", "4"), algebra(), plan))));
        const auto split = take(measure(deck, triangles));
        const auto split_view = take(split.view());
        CHECK(split_view.faces.size() == 2);

        const auto baseline = plan.remaining_bytes();
        const auto partial = measure(deck, std::span{triangles}.first(1));
        REQUIRE_FALSE(partial);
        CHECK(take(diagnostic).reason == Reason::Coverage);
        CHECK(plan.remaining_bytes() < baseline);
        const auto after_partial = plan.remaining_bytes();
        const ExternalPgBoundLimits no_work{limits.region, 0};
        const auto exhausted = NastranPlanarPgBounds::prepare(
            target, deck, SiLengthUnit::Millimetre, forward, no_work, plan);
        REQUIRE_FALSE(exhausted);
        CHECK(exhausted.error().code == Code::ResourceExceeded);
        CHECK(plan.remaining_bytes() < after_partial);
        CHECK(take(deck.view()).input.bulk.digest == inspected.digest);
        CHECK(take(target.identity()) == inspected.patch);
        return first;
    }();

    const auto recorded_view = take(recorded.view());
    const auto& row = recorded_view.faces.front();
    REQUIRE(row.nastran_offset);
    if (const auto& placed = row.nastran_offset) {
        CHECK(placed->signed_metres == f.r("1", "2"));
    }
    // Independent map before comparison:
    // ((s+1/2)/sqrt(2), t, (1/2-s)/sqrt(2)).
    // The squared residual is convex and achieves 1/4 at s=0.
    CHECK(std::is_gteq(take(row.distance_upper_metres.compare(f.r("1", "2"), algebra(), plan))));
    CHECK(std::is_lt(take(
        row.distance_upper_metres.compare(f.r("500000000001", "1000000000000"), algebra(), plan))));
    // Chord from (1,0,1)/sqrt(2) to +Z has squared length 2-sqrt(2).
    const auto square =
        take(row.normal_chord_upper.multiply(row.normal_chord_upper, algebra(), plan));
    const auto excess = take(f.r("2").subtract(square, algebra(), plan));
    REQUIRE(excess.sign() > 0);
    const auto exact_test = take(excess.multiply(excess, algebra(), plan));
    CHECK(std::is_lteq(take(exact_test.compare(f.r("2"), algebra(), plan))));
    CHECK(std::is_lt(take(row.normal_chord_upper.compare(f.r("766", "1000"), algebra(), plan))));

    std::string warped{source};
    const auto node = warped.find("GRID,95,5,0.0,1000.0,0.0");
    REQUIRE(node != std::string::npos);
    warped.replace(node, std::string_view{"GRID,95,5,0.0,1000.0,0.0"}.size(),
                   "GRID,95,5,0.0,1000.0,1.0");
    const auto nonplanar = read_cards(warped, plan);
    const auto rejected_quad = measure(nonplanar, forward);
    REQUIRE_FALSE(rejected_quad);
    CHECK(rejected_quad.error().code == Code::GeometryInvalid);
    CHECK(take(diagnostic).reason == Reason::Planarity);
    CHECK(take(diagnostic).node == 95);
}

TEST_CASE("Nastran face bounds include the interior gap from a curved PG graph",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][negative][verification]") {
    Fixture fixture;
    core::ResourcePlan plan{budget};
    const auto target = pg(fixture, plan, true);
    const Corners c{fixture};
    const std::array ring{ExternalPgCorner{10, c.a}, ExternalPgCorner{40, c.b},
                          ExternalPgCorner{90, c.c}, ExternalPgCorner{95, c.d}};
    const std::array second{ring[0], ring[2], ring[3]};
    const std::array reversed{ring[0], ring[3], ring[2], ring[1]};
    const std::array quad{ExternalPgFace{101, Winding::Forward, ring}};
    const std::array split{ExternalPgFace{201, Winding::Forward, std::span{ring}.first(3)},
                           ExternalPgFace{202, Winding::Forward, second}};
    const std::array reverse{ExternalPgFace{102, Winding::Reversed, reversed}};
    const std::array wrong{ExternalPgFace{102, Winding::Forward, reversed}};
    const auto deck = cards("CQUAD4,101,20,10,40,90,95\nCQUAD4,102,20,10,95,90,40\n"
                            "CTRIA3,201,20,10,40,90\nCTRIA3,202,20,10,90,95\n",
                            plan);
    const ExternalPgBoundLimits limits{{reference_limits(), 2}, 200'000};
    std::optional<ExternalPgDiagnostic> diagnostic;
    const auto measure = [&](std::span<const ExternalPgFace> selection) {
        return NastranPlanarPgBounds::prepare(target, deck, SiLengthUnit::Millimetre, selection,
                                              limits, plan, &diagnostic);
    };
    const auto compare = [&](const ExactRational& actual, const ExactRational& bound) {
        return take(actual.compare(bound, algebra(), plan));
    };
    // The two positions agree at every corner; the extra height is
    // 2*(s^2-s)/10000 m, so the whole-face maximum is attained at s=1/2.
    for (const auto selection :
         {std::span<const ExternalPgFace>{quad}, std::span<const ExternalPgFace>{split},
          std::span<const ExternalPgFace>{reverse}}) {
        const auto result = take(measure(selection));
        const auto view = take(result.view());
        CHECK(view.patch == take(target.identity()));
        CHECK(view.faces.size() == selection.size());
        for (const auto& face : view.faces) {
            CHECK(std::is_gteq(compare(face.distance_upper_metres, fixture.r("1", "20000"))));
            CHECK(std::is_lt(
                compare(face.distance_upper_metres, fixture.r("500000001", "10000000000000"))));
            // Local unnormalized normals at s=0 differ by (1,1/4,0)/10000.
            // After normalization, the x-component difference exceeds 1/30000.
            CHECK(std::is_gt(compare(face.normal_chord_upper, fixture.r("1", "30000"))));
            CHECK(std::is_lt(compare(face.normal_chord_upper, fixture.r("1", "1000"))));
        }
    }
    const auto rejected_winding = measure(wrong);
    REQUIRE_FALSE(rejected_winding);
    CHECK(rejected_winding.error().code == Code::GeometryInvalid);
    CHECK(take(diagnostic).reason == Reason::SourceWinding);
    CHECK(take(diagnostic).element == 102);
}

[[nodiscard]] PgReferenceGeometry gmsh_bound_target(Fixture& f, core::ResourcePlan& plan,
                                                    bool curved) {
    const auto zero = f.r("0");
    const auto one = f.r("1");
    if (curved) {
        const auto tenth = f.r("1", "10");
        const auto minus_tenth = f.r("-1", "10");
        const std::array terms{ExactInputTerm{{1, 0}, minus_tenth}, ExactInputTerm{{1, 1}, tenth},
                               ExactInputTerm{{2, 0}, tenth}};
        const auto height = take(ExactChartPolynomial::make(terms, algebra(), plan));
        f.graph = take(PgGraph::prepare({zero, one, zero, height, 1, PgGraphConvention::Version1},
                                        {algebra(), 64, {}}, plan));
    }
    const std::array<ExactRational, 3> origin;
    std::array<ExactRational, 9> rotation;
    rotation[0] = f.r("1");
    rotation[4] = f.r("1");
    rotation[8] = f.r("1");
    return flat_target(f, {rotation, origin, one}, plan);
}

[[nodiscard]] std::string bilinear_mesh(std::string_view last_node, std::string_view elements) {
    return "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 4 10 95\n2 17 0 4\n"
           "10\n40\n90\n95\n0 0 0\n1000 0 0\n1000 1000 100\n" +
           std::string{last_node} + "\n$EndNodes\n$Elements\n" + std::string{elements} +
           "$EndElements\n";
}

TEST_CASE("Gmsh quad bounds include curved interior distance and the original bilinear normals",
          "[IO4-04][SH-D02][N4-02][R-05][integration][negative][verification]") {
    const bool curved = GENERATE(false, true);
    INFO("Curved PG: " << curved);
    Fixture f;
    core::ResourcePlan plan{budget};
    const auto a = f.point("0", "0"), b = f.point("1", "0"), c = f.point("1", "1"),
               d = f.point("0", "1");
    const std::array ring{ExternalPgCorner{10, a}, ExternalPgCorner{40, b}, ExternalPgCorner{90, c},
                          ExternalPgCorner{95, d}};
    const std::array reversed{ring[0], ring[3], ring[2], ring[1]};
    const std::array charts{ExternalPgFace{701, Winding::Forward, ring}};
    const std::array reverse_charts{ExternalPgFace{702, Winding::Reversed, reversed}};
    constexpr std::array<std::uint64_t, 1> selected{701};
    constexpr std::array<std::uint64_t, 1> selected_reverse{702};
    const ExternalPgBoundLimits limits{{reference_limits(), 2}, 50'000};
    const auto text = bilinear_mesh("0 1000 0", "1 2 701 702\n2 17 3 2\n701 10 40 90 95\n"
                                                "702 10 95 90 40\n");
    const auto prepared = [&] {
        const auto target = gmsh_bound_target(f, plan, curved);
        const auto msh = mesh(text, plan);
        const auto make_topology = [&](std::span<const std::uint64_t> faces) {
            return take(GmshSurfaceTopology::prepare(msh, faces, {2, 8, 8, 2, 4096}, plan, {}));
        };
        const auto topology = make_topology(selected);
        const auto topology_reverse = make_topology(selected_reverse);
        auto bounds = take(GmshPgBounds::prepare(target, msh, topology, SiLengthUnit::Millimetre,
                                                 charts, limits, plan));
        const auto inspected = take(bounds.view());
        CHECK(inspected.patch == take(target.identity()));
        CHECK(inspected.digest == take(msh.view()).digest);
        CHECK(inspected.source.source == SourceId{17});
        CHECK(inspected.directed_work > 0);
        const auto mirror = take(GmshPgBounds::prepare(
            target, msh, topology_reverse, SiLengthUnit::Millimetre, reverse_charts, limits, plan));
        const auto mirror_view = take(mirror.view());
        const auto& mir = mirror_view.faces.front();
        CHECK(mir.winding == Winding::Reversed);
        CHECK(mir.ring()[1].node == 95);
        CHECK(mir.distance_upper_metres == inspected.faces.front().distance_upper_metres);
        CHECK(mir.normal_chord_upper == inspected.faces.front().normal_chord_upper);

        const ExternalPgBoundLimits exhausted{limits.region, inspected.directed_work - 1};
        const auto denied = GmshPgBounds::prepare(target, msh, topology, SiLengthUnit::Millimetre,
                                                  charts, exhausted, plan);
        REQUIRE_FALSE(denied);
        CHECK(denied.error().code == Code::ResourceExceeded);
        const auto bounded =
            GmshPgBounds::prepare(target, msh, topology, SiLengthUnit::Millimetre, charts,
                                  {limits.region, inspected.directed_work}, plan);
        REQUIRE(bounded);
        return bounds;
    }();

    const auto output = take(prepared.view());
    REQUIRE(output.faces.size() == 1);
    const auto& row = output.faces.front();
    CHECK(row.element == 701);
    CHECK((row.entity == GmshKey{2, 17}));
    REQUIRE(row.ring().size() == 4);
    const auto& lifted = row.ring()[2];
    CHECK(lifted.node == 90);
    CHECK(lifted.chart.s == f.r("1"));
    CHECK(lifted.chart.t == f.r("1"));
    CHECK(lifted.position_metres[2] == f.r("1", "10"));
    CHECK(std::string_view{text}.substr(lifted.source_bytes[2], 3) == "100");
    const auto compare = [&](const ExactRational& actual, const ExactRational& bound) {
        return take(actual.compare(bound, algebra(), plan));
    };
    if (curved) {
        // Actual bilinear source z=st/10. PG z=(st+s^2-s)/10.
        // All corners agree, but the exact maximum gap at s=1/2 is 1/40 m.
        CHECK(std::is_gteq(compare(row.distance_upper_metres, f.r("1", "40"))));
        CHECK(std::is_lt(compare(row.distance_upper_metres, f.r("501", "10000"))));
        // At (0,0), normals are (0,0,1) and (1,0,10)/sqrt(101).
        const auto square =
            take(row.normal_chord_upper.multiply(row.normal_chord_upper, algebra(), plan));
        const auto gap = take(f.r("2").subtract(square, algebra(), plan));
        REQUIRE(gap.sign() > 0);
        const auto check = take(gap.multiply(gap, algebra(), plan));
        CHECK(std::is_lteq(compare(check, f.r("400", "101"))));
    } else {
        // Bilinear z=st/10 versus z=0. The greatest normal chord is at (1,1),
        // where the source direction is (-1,-1,10)/sqrt(102). Both fixed
        // diagonal triangles, or the flat first corner, give smaller chords.
        CHECK(std::is_gteq(compare(row.distance_upper_metres, f.r("1", "10"))));
        CHECK(std::is_lt(compare(row.distance_upper_metres, f.r("1001", "10000"))));
        const auto square =
            take(row.normal_chord_upper.multiply(row.normal_chord_upper, algebra(), plan));
        const auto gap = take(f.r("2").subtract(square, algebra(), plan));
        REQUIRE(gap.sign() > 0);
        const auto check = take(gap.multiply(gap, algebra(), plan));
        CHECK(std::is_lteq(compare(check, f.r("200", "51"))));
        CHECK(std::is_lt(compare(row.normal_chord_upper, f.r("143", "1000"))));
    }
}

TEST_CASE("Gmsh PG bounds reject orientation reversal and unattested selections",
          "[IO4-04][SH-D02][R-05][integration][negative][verification]") {
    Fixture f;
    core::ResourcePlan plan{budget};
    const auto target = gmsh_bound_target(f, plan, false);
    const auto a = f.point("0", "0"), b = f.point("1", "0"), c = f.point("1", "1"),
               d = f.point("0", "1");
    const std::array ring{ExternalPgCorner{10, a}, ExternalPgCorner{40, b}, ExternalPgCorner{90, c},
                          ExternalPgCorner{95, d}};
    const std::array charts{ExternalPgFace{701, Winding::Forward, ring}};
    constexpr std::array<std::uint64_t, 1> selected{701};
    constexpr std::string_view elements = "1 1 701 701\n2 17 3 1\n701 10 40 90 95\n";
    const auto text = bilinear_mesh("2000 1000 0", elements);
    const auto folded = mesh(text, plan);
    const auto regular = mesh(bilinear_mesh("0 1000 0", elements), plan);
    const auto topology =
        take(GmshSurfaceTopology::prepare(folded, selected, {1, 4, 4, 1, 4096}, plan, {}));
    const ExternalPgBoundLimits limits{{reference_limits(), 2}, 50'000};
    std::optional<ExternalPgDiagnostic> diagnostic;
    const auto measure = [&](const GmshMesh& source, std::span<const ExternalPgFace> faces) {
        return GmshPgBounds::prepare(target, source, topology, SiLengthUnit::Millimetre, faces,
                                     limits, plan, &diagnostic);
    };
    const auto stale = measure(regular, charts);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == Code::StaleGeneration);
    CHECK(take(diagnostic).reason == Reason::Source);

    const auto swapped = std::array{ring[0], ring[1], ring[3], ring[2]};
    const std::array wrong{ExternalPgFace{701, Winding::Forward, swapped}};
    const auto source_order = measure(folded, wrong);
    REQUIRE_FALSE(source_order);
    CHECK(source_order.error().code == Code::InvalidInput);
    CHECK(take(diagnostic).reason == Reason::SourceOrder);

    // The chart is a valid square, but the physical Jacobian's z component
    // is 1-2v. Its projection reverses within the source ring.
    const auto bad_geometry = measure(folded, charts);
    REQUIRE_FALSE(bad_geometry);
    CHECK(bad_geometry.error().code == Code::GeometryInvalid);
    const auto reason = take(diagnostic);
    CHECK(reason.reason == Reason::SourceWinding);
    CHECK(reason.node == 90);
    CHECK(reason.byte == text.find("1000 1000 100"));
    CHECK(take(folded.view()).digest == take(topology.view()).source.digest);
}
TEST_CASE("Gmsh source triangles keep their own normals and cover the PG chart",
          "[IO4-04][SH-D02][R-05][integration][verification]") {
    Fixture f;
    core::ResourcePlan plan{budget};
    const auto target = gmsh_bound_target(f, plan, false);
    const auto a = f.point("0", "0"), b = f.point("1", "0"), c = f.point("1", "1"),
               d = f.point("0", "1");
    const std::array first{ExternalPgCorner{10, a}, ExternalPgCorner{40, b},
                           ExternalPgCorner{90, c}};
    const std::array second{first[0], first[2], ExternalPgCorner{95, d}};
    const std::array charts{ExternalPgFace{801, Winding::Forward, first},
                            ExternalPgFace{802, Winding::Forward, second}};
    constexpr std::array<std::uint64_t, 2> selected{801, 802};
    const auto source = mesh(bilinear_mesh("0 1000 0", "1 2 801 802\n2 17 2 2\n"
                                                       "801 10 40 90\n802 10 90 95\n"),
                             plan);
    const auto topology =
        take(GmshSurfaceTopology::prepare(source, selected, {2, 6, 4, 1, 4096}, plan, {}));
    const ExternalPgBoundLimits limits{{reference_limits(), 2}, 50'000};
    std::optional<ExternalPgDiagnostic> diagnostic;
    const auto measure = [&](const GmshSurfaceTopology& faces,
                             std::span<const ExternalPgFace> given) {
        return GmshPgBounds::prepare(target, source, faces, SiLengthUnit::Millimetre, given, limits,
                                     plan, &diagnostic);
    };
    const auto owned = take(measure(topology, charts));
    const auto result = take(owned.view());
    REQUIRE(result.faces.size() == 2);
    for (const auto& row : result.faces) {
        REQUIRE(row.ring().size() == 3);
        // The two original triangles have z=t/10 and z=s/10 respectively.
        // In either case the maximum distance is .1 m and chord is just below .1.
        CHECK(
            std::is_gteq(take(row.distance_upper_metres.compare(f.r("1", "10"), algebra(), plan))));
        CHECK(std::is_gt(take(row.normal_chord_upper.compare(f.r("99", "1000"), algebra(), plan))));
        CHECK(std::is_lt(take(row.normal_chord_upper.compare(f.r("1", "10"), algebra(), plan))));
    }
    const auto partial = take(GmshSurfaceTopology::prepare(source, std::span{selected}.first(1),
                                                           {1, 3, 3, 1, 4096}, plan, {}));
    const auto uncovered = measure(partial, std::span{charts}.first(1));
    REQUIRE_FALSE(uncovered);
    CHECK(uncovered.error().code == Code::GeometryInvalid);
    CHECK(take(diagnostic).reason == Reason::Coverage);
}

TEST_CASE("Nastran triangle distance bound needs the mixed PG height derivative",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][verification]") {
    Fixture f;
    core::ResourcePlan plan{budget};
    const auto zero = f.r("0"), one = f.r("1");
    const std::array terms{ExactInputTerm{{1, 1}, one}};
    const auto height = take(ExactChartPolynomial::make(terms, algebra(), plan));
    f.graph = take(PgGraph::prepare({zero, one, zero, height, 1, PgGraphConvention::Version1},
                                    {algebra(), 64, {}}, plan));
    f.domain[2] = f.point("0", "1");
    f.domain.resize(3);
    for (auto& thickness : f.heights) {
        thickness = f.r("1", "10");
    }
    const std::array<ExactRational, 3> origin;
    std::array<ExactRational, 9> rotation;
    rotation[0] = f.r("1");
    rotation[4] = f.r("1");
    rotation[8] = f.r("1");
    const auto target = flat_target(f, {rotation, origin, one}, plan);

    const auto deck = read_cards("GRID,10,,0.0,0.0,0.0\nGRID,40,,1.0,0.0,0.0\n"
                                 "GRID,90,,0.0,1.0,0.0\nPSHELL,20\nCTRIA3,101,20,10,40,90\n",
                                 plan);
    const std::array ring{ExternalPgCorner{10, f.domain[0]}, ExternalPgCorner{40, f.domain[1]},
                          ExternalPgCorner{90, f.domain[2]}};
    const std::array selection{ExternalPgFace{101, Winding::Forward, ring}};
    const NastranPgBoundLimits limits{{reference_limits(), 1}, 200'000};
    const auto result = take(
        NastranPlanarPgBounds::prepare(target, deck, SiLengthUnit::Metre, selection, limits, plan));
    const auto faces = take(result.view()).faces;
    REQUIRE(faces.size() == 1);
    REQUIRE(faces.front().ring().size() == 3);
    // Source (s,t,0), PG (s,t,st). All vertex errors and both pure second
    // derivatives are zero. The maximum gap is 1/4 at (1/2,1/2).
    const auto quarter = f.r("1", "4");
    CHECK(
        std::is_gteq(take(faces.front().distance_upper_metres.compare(quarter, algebra(), plan))));
}

template <class Measure>
void check_nastran_edge_conformance(Measure measure, Reason slit_reason,
                                    std::uint64_t slit_element) {
    Fixture f;
    core::ResourcePlan plan{budget};
    // This named subset splits the atomic PG bottom edge but not the source edge.
    f.named.push_back(f.point("0", "0"));
    f.named.push_back({f.r("1", "4"), f.r("0")});
    f.names.push_back({"part of bottom", f.named});
    const auto length = f.r("1");
    const std::array<ExactRational, 3> origin;
    const std::array rotation{f.r("1"), f.r("0"), f.r("0"), f.r("0"), f.r("1"),
                              f.r("0"), f.r("0"), f.r("0"), f.r("1")};
    const auto target = flat_target(f, {rotation, origin, length}, plan);
    const PgPartitionPoint a = f.point("0", "0"), b = f.point("1", "0"), c = f.point("1", "1"),
                           d = f.point("0", "1");
    const PgPartitionPoint bottom_mid{f.r("1", "2"), f.r("0")};
    const PgPartitionPoint right_mid{f.r("1"), f.r("1", "2")};
    const PgPartitionPoint top_mid{f.r("1", "2"), f.r("1")};
    const PgPartitionPoint left_mid{f.r("0"), f.r("1", "2")};
    const PgPartitionPoint centre{f.r("1", "2"), f.r("1", "2")};
    const ExternalPgCorner bl{1, a}, br{3, b}, tr{6, c}, tl{8, d}, bm{2, bottom_mid},
        rm{5, right_mid}, tm{7, top_mid}, lm{9, left_mid}, cm{4, centre};
    const auto deck =
        read_cards("GRID,1,,0.0,0.0,0.0\nGRID,2,,0.5,0.0,0.0\nGRID,3,,1.0,0.0,0.0\n"
                   "GRID,4,,0.5,0.5,0.0\nGRID,5,,1.0,0.5,0.0\nGRID,6,,1.0,1.0,0.0\n"
                   "GRID,7,,0.5,1.0,0.0\nGRID,8,,0.0,1.0,0.0\nGRID,9,,0.0,0.5,0.0\n"
                   "GRID,11,,0.0,0.0,0.0\nGRID,16,,1.0,1.0,0.0\nPSHELL,20\n"
                   "CTRIA3,101,20,1,3,6\nCTRIA3,102,20,1,6,8\nCTRIA3,103,20,11,16,8\n"
                   "CTRIA3,104,20,1,8,6\n"
                   "CQUAD4,201,20,1,2,7,8\nCQUAD4,202,20,2,3,5,4\nCQUAD4,203,20,4,5,6,7\n"
                   "CQUAD4,204,20,1,2,4,9\nCQUAD4,205,20,9,4,7,8\n",
                   plan);
    std::optional<ExternalPgDiagnostic> diagnostic;
    const auto check = [&](std::span<const ExternalPgFace> charts) {
        return measure(target, deck, charts, plan, &diagnostic);
    };
    const std::array lower_tri{bl, br, tr}, upper_tri{bl, tr, tl};
    const std::array slit_tri{ExternalPgCorner{11, a}, ExternalPgCorner{16, c}, tl};
    const std::array backwards_tri{bl, tl, tr};
    const std::array attached{ExternalPgFace{101, Winding::Forward, lower_tri},
                              ExternalPgFace{102, Winding::Forward, upper_tri}};
    const std::array slit{attached[0], ExternalPgFace{103, Winding::Forward, slit_tri}};
    const std::array reversed_edge{attached[0],
                                   ExternalPgFace{104, Winding::Reversed, backwards_tri}};

    SECTION("Same coordinates with different GRID tags") {
        REQUIRE(check(attached));
        const auto split = check(slit);
        REQUIRE_FALSE(split);
        CHECK(split.error().code == Code::GeometryInvalid);
        CHECK(take(diagnostic).element == slit_element);
        CHECK(take(diagnostic).reason == slit_reason);
    }
    SECTION("Same raw direction after one declared reversal") {
        const auto wrong_winding = check(reversed_edge);
        REQUIRE_FALSE(wrong_winding);
        CHECK(wrong_winding.error().code == Code::GeometryInvalid);
        CHECK(take(diagnostic).reason == Reason::EdgeIncidence);
        CHECK(take(diagnostic).other == 104);
    }
    SECTION("Interior T-junction") {
        const std::array left{bl, bm, tm, tl}, right_bottom{bm, br, rm, cm},
            right_top{cm, rm, tr, tm};
        const std::array left_bottom{bl, bm, cm, lm}, left_top{lm, cm, tm, tl};
        const std::array t_junction{ExternalPgFace{201, Winding::Forward, left},
                                    ExternalPgFace{202, Winding::Forward, right_bottom},
                                    ExternalPgFace{203, Winding::Forward, right_top}};
        const std::array conforming{ExternalPgFace{204, Winding::Forward, left_bottom},
                                    ExternalPgFace{205, Winding::Forward, left_top}, t_junction[1],
                                    t_junction[2]};
        REQUIRE(check(conforming));
        const auto unmatched = check(t_junction);
        REQUIRE_FALSE(unmatched);
        CHECK(unmatched.error().code == Code::GeometryInvalid);
        CHECK(take(diagnostic).reason == Reason::EdgeIncidence);
        CHECK(take(diagnostic).element == 201);
        CHECK(take(diagnostic).node == 2);
    }
}

TEST_CASE("External planar certification requires source edge incidence across a PG patch",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][negative][verification]") {
    check_nastran_edge_conformance(
        [](const PgReferenceGeometry& target, const NastranShellCards& deck,
           std::span<const ExternalPgFace> charts, core::ResourcePlan& plan,
           std::optional<ExternalPgDiagnostic>* diagnostic) {
            return ExternalPlanarPg::prepare(target, deck, SiLengthUnit::Metre, charts,
                                             {reference_limits(), 4}, plan, diagnostic);
        },
        Reason::EdgeIncidence, 101);
}

TEST_CASE("Nastran bounded placement requires source edge incidence across a PG patch",
          "[IO4-03][IO4-04][SH-D02][N4-02][integration][negative][verification]") {
    check_nastran_edge_conformance(
        [](const PgReferenceGeometry& target, const NastranShellCards& deck,
           std::span<const ExternalPgFace> charts, core::ResourcePlan& plan,
           std::optional<ExternalPgDiagnostic>* diagnostic) {
            const NastranPgBoundLimits limits{{reference_limits(), 4}, 200'000};
            return NastranPlanarPgBounds::prepare(target, deck, SiLengthUnit::Metre, charts, limits,
                                                  plan, diagnostic);
        },
        Reason::SourceOrder, 103);
}

} // namespace
