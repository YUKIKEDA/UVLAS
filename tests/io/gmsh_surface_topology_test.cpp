#include "gmsh_surface_topology.hpp"
#include "source_digest.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using enum core::FailureCode;
using Reason = GmshTopologyFailure;

constexpr std::uint64_t memory = 64000;
constexpr GmshTopologyLimits limits{24, 96, 96, 24, 20000};
constexpr auto sparse = std::uint64_t{9007199254740993};

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw std::logic_error("Expected a value, code=" +
                               std::to_string(static_cast<unsigned>(result.error().code)));
    return std::move(*result);
}

template <class T>
T take(std::optional<T> value) {
    if (!value)
        throw std::logic_error("Expected a view");
    return *value;
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
    std::size_t offset_ = 0;
};

GmshMesh read_mesh(std::string_view bytes) {
    const auto raw = std::as_bytes(std::span{bytes});
    Reader reader{raw};
    core::ResourcePlan plan{memory};
    constexpr RootSnapshotLocation location{
        {{SourceId{19}, 27, {}}, 1, 1}, SourceId{81}, ReferenceRole::GmshMesh, 5};
    AcquiredRootSnapshot source{
        location,
        take(read_reference_file(reader, source_digest(raw),
                                 {bytes.size(), bytes.size() + 1, 3, bytes.size()}, plan, {},
                                 nullptr)),
        {}};
    return take(GmshMesh::prepare(source, {bytes.size(), 0, 0, 0, 128, 64, 32000}, plan, {}));
}

struct Element {
    std::uint64_t tag;
    std::initializer_list<std::uint64_t> nodes;
};

struct NodePosition {
    std::uint64_t tag;
    std::string_view coordinates;
};

std::string mesh_text(std::initializer_list<Element> records,
                      std::initializer_list<NodePosition> coordinates = {}) {
    std::vector<std::uint64_t> nodes;
    for (const auto& record : records)
        nodes.insert(nodes.end(), record.nodes.begin(), record.nodes.end());
    std::ranges::sort(nodes);
    const auto duplicate = std::ranges::unique(nodes);
    nodes.erase(duplicate.begin(), duplicate.end());
    const auto extremal = std::ranges::minmax_element(records, {}, &Element::tag);
    std::string identity;
    std::string positions;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        identity += std::to_string(nodes[i]) + "\n";
        const auto point = std::ranges::find(coordinates, nodes[i], &NodePosition::tag);
        if (point == coordinates.end())
            positions += std::to_string(i) + " " + std::to_string(i * i) + " 0\n";
        else
            positions += std::string{point->coordinates} + "\n";
    }
    std::size_t blocks = 0;
    std::string elements;
    for (const auto width : {2U, 3U, 4U}) {
        const auto count = std::ranges::count_if(
            records, [width](const Element& row) { return row.nodes.size() == width; });
        if (count == 0)
            continue;
        ++blocks;
        elements += (width == 2 ? "1 9 " : "2 17 ") + std::to_string(width - 1) + " " +
                    std::to_string(count) + "\n";
        for (const auto& row : records) {
            if (row.nodes.size() != width)
                continue;
            elements += std::to_string(row.tag);
            for (const auto node : row.nodes)
                elements += " " + std::to_string(node);
            elements += "\n";
        }
    }
    const auto count = std::to_string(nodes.size());
    return "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 " + count + " " +
           std::to_string(nodes.front()) + " " + std::to_string(nodes.back()) + "\n2 17 0 " +
           count + "\n" + identity + positions + "$EndNodes\n$Elements\n" + std::to_string(blocks) +
           " " + std::to_string(records.size()) + " " + std::to_string(extremal.min->tag) + " " +
           std::to_string(extremal.max->tag) + "\n" + elements + "$EndElements\n";
}

std::string mixed_mesh() {
    return mesh_text({{90, {11, 42, sparse}},
                      {7, {sparse, 42, 30, 29}},
                      {501, {101, 204, 1010}},
                      {78, {11, 101}}},
                     {{11, "0 0 0"},
                      {42, "1 0 0"},
                      {sparse, "1 1 0"},
                      {30, "2 0 0"},
                      {29, "2 1 0"},
                      {101, "3 0 0"},
                      {204, "4 0 0"},
                      {1010, "3 1 0"}});
}

TEST_CASE("Gmsh topology owns chosen mixed face rings, reciprocal edges and interleaved components",
          "[IO4-04][IO4-01][verification]") {
    const auto source = mixed_mesh();
    constexpr std::array<std::uint64_t, 3> selection{7, 501, 90};
    core::ResourcePlan plan{memory};
    auto result = [&] {
        const auto mesh = read_mesh(source);
        const auto before = take(mesh.view());
        auto topology = take(GmshSurfaceTopology::prepare(mesh, selection, limits, plan, {}));
        const auto after = take(mesh.view());
        CHECK(after.bytes == source);
        CHECK(after.bytes.data() == before.bytes.data());
        CHECK(after.digest == before.digest);
        return topology;
    }();
    const auto view = take(result.view());
    CHECK(view.source.digest == source_digest(std::as_bytes(std::span{source})));
    CHECK(view.source.location.source == SourceId{81});
    REQUIRE(view.faces.size() == 3);
    REQUIRE(view.edges.size() == 10);
    REQUIRE(view.components.size() == 2);
    CHECK(view.faces[0].element.tag == 7);
    CHECK(view.faces[0].element.entity == GmshKey{2, 17});
    CHECK(view.faces[0].element.kind == GmshElementKind::Quadrangle);
    CHECK(view.faces[1].element.tag == 501);
    CHECK(view.faces[2].element.tag == 90);
    CHECK(view.faces[0].first_edge == 0);
    CHECK(view.faces[1].first_edge == 4);
    CHECK(view.faces[2].first_edge == 7);
    CHECK(view.faces[0].component == 0);
    CHECK(view.faces[1].component == 1);
    CHECK(view.faces[2].component == 0);
    CHECK(view.components[0].first_face == 0);
    CHECK(view.components[0].face_count == 2);
    CHECK(view.components[0].boundary_count == 5);
    CHECK(view.components[1].first_face == 1);
    CHECK(view.components[1].face_count == 1);
    CHECK(view.components[1].boundary_count == 3);
    CHECK(view.edges[0].from == sparse);
    CHECK(view.edges[0].to == 42);
    CHECK(view.edges[0].opposite == 8);
    CHECK(view.edges[8].opposite == 0);
    CHECK(view.edges[3].next == 0);
    CHECK(view.edges[9].next == 7);
    constexpr std::array<std::size_t, 8> boundary{1, 2, 3, 4, 5, 6, 7, 9};
    CHECK(std::ranges::equal(view.boundary, boundary));
    const auto moved = std::move(result);
    // NOLINTNEXTLINE(bugprone-use-after-move): moving clears the previous view.
    CHECK_FALSE(result.view());
    CHECK(take(moved.view()).edges.size() == 10);
}

TEST_CASE("Gmsh topology accepts closed vertex fans and keeps a chosen subcomplex boundary",
          "[IO4-04][verification]") {
    const auto source =
        mesh_text({{11, {1, 2, 3}}, {22, {1, 4, 2}}, {33, {2, 4, 3}}, {44, {3, 4, 1}}},
                  {{1, "0 0 0"}, {2, "1 0 0"}, {3, "0 1 0"}, {4, "0 0 1"}});
    const auto mesh = read_mesh(source);
    constexpr std::array<std::uint64_t, 4> shell{44, 22, 11, 33};
    constexpr std::array<std::uint64_t, 1> face{22};
    core::ResourcePlan plan{memory};
    const auto closed = take(GmshSurfaceTopology::prepare(mesh, shell, limits, plan, {}));
    const auto view = take(closed.view());
    REQUIRE(view.components.size() == 1);
    CHECK(view.components[0].face_count == 4);
    CHECK(view.components[0].boundary_count == 0);
    CHECK(view.boundary.empty());
    for (std::size_t i = 0; i < view.edges.size(); ++i) {
        const auto& edge = view.edges[i];
        const auto opposite = take(edge.opposite);
        CHECK(view.edges[opposite].opposite == i);
        CHECK(view.edges[opposite].from == edge.to);
        CHECK(view.edges[opposite].to == edge.from);
    }
    const auto open = take(GmshSurfaceTopology::prepare(mesh, face, limits, plan, {}));
    CHECK(take(open.view()).boundary.size() == 3);
    core::ResourcePlan no_storage{0};
    const auto empty =
        take(GmshSurfaceTopology::prepare(mesh, {}, {0, 0, 0, 0, 0}, no_storage, {}));
    CHECK(take(empty.view()).faces.empty());
    CHECK(take(empty.view()).components.empty());
}

TEST_CASE("Gmsh topology keeps coincident but differently tagged faces distinct",
          "[IO4-04][verification]") {
    const auto source = mesh_text(
        {{11, {1, 2, 3}}, {22, {4, 5, 6}}},
        {{1, "0 0 0"}, {2, "1 0 0"}, {3, "0 1 0"}, {4, "0 0 0"}, {5, "1 0 0"}, {6, "0 1 0"}});
    const auto mesh = read_mesh(source);
    constexpr std::array<std::uint64_t, 2> selection{11, 22};
    core::ResourcePlan plan{memory};
    const auto topology = take(GmshSurfaceTopology::prepare(mesh, selection, limits, plan, {}));
    const auto view = take(topology.view());
    REQUIRE(view.components.size() == 2);
    CHECK(view.boundary.size() == 6);
    CHECK(view.faces[0].component != view.faces[1].component);
}

TEST_CASE("Gmsh topology diagnoses invalid incidence without changing the source mesh",
          "[IO4-04][negative][verification]") {
    const auto check = [](const std::string& source, std::initializer_list<std::uint64_t> selection,
                          GmshTopologyDiagnostic expected) {
        const auto mesh = read_mesh(source);
        core::ResourcePlan plan{memory};
        GmshTopologyDiagnostic found;
        const auto result = GmshSurfaceTopology::prepare(mesh, selection, limits, plan, {}, &found);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == InvalidInput);
        CHECK(found.reason == expected.reason);
        CHECK(found.face == expected.face);
        CHECK(found.other_face == expected.other_face);
        CHECK(found.vertex == expected.vertex);
        CHECK(found.other_vertex == expected.other_vertex);
        CHECK(take(mesh.view()).bytes == source);
    };
    SECTION("unselected elements do not repair a missing face tag") {
        check(mixed_mesh(), {7, 404}, {Reason::UndefinedFace, 404});
    }
    SECTION("lines are not selected surface faces") {
        check(mixed_mesh(), {78}, {Reason::NonSurfaceElement, 78});
    }
    SECTION("a source face cannot be selected twice") {
        check(mixed_mesh(), {7, 90, 7}, {Reason::RepeatedSelection, 7, 7});
    }
    SECTION("a face cannot repeat a vertex in its ring") {
        check(mesh_text({{1, {10, 20, 30, 20}}}), {1}, {Reason::RepeatedVertex, 1, 0, 20});
    }
    SECTION("a reversed cyclic copy of a face has the same boundary") {
        check(mesh_text({{1, {11, 20, 30, sparse}}, {2, {30, 20, 11, sparse}}}), {1, 2},
              {Reason::DuplicateFace, 2, 1});
    }
    SECTION("three faces cannot share either direction of one edge") {
        check(mesh_text({{1, {10, 20, 30}}, {2, {20, 10, 40}}, {3, {10, 20, 50}}}), {1, 2, 3},
              {Reason::NonManifoldEdge, 1, 3, 10, 20});
    }
    SECTION("a shared edge must have opposite directions") {
        check(mesh_text({{1, {10, 20, 30}}, {2, {10, 20, 40}}}), {1, 2},
              {Reason::InconsistentOrientation, 1, 2, 10, 20});
    }
    SECTION("a connected sheet may still have two fans at one vertex") {
        check(mesh_text({{10, {1, 2, 3}}, {20, {3, 2, 4}}, {30, {4, 2, 5}}, {40, {1, 4, 5}}}),
              {10, 20, 30, 40}, {Reason::DisconnectedVertexFan, 10, 40, 1});
    }
}

TEST_CASE("Gmsh topology accepts exact planned limits and leaves the mesh usable on rejection",
          "[IO4-04][N4-02][negative][verification]") {
    const auto source = mixed_mesh();
    const auto mesh = read_mesh(source);
    constexpr std::array<std::uint64_t, 3> selection{7, 501, 90};
    core::ResourcePlan plan{memory};
    const auto measured = take(GmshSurfaceTopology::prepare(mesh, selection, limits, plan, {}));
    const auto work = take(measured.view()).work;
    const auto bytes = plan.allocated_bytes();
    REQUIRE(work > 0);
    REQUIRE(bytes > 0);
    core::ResourcePlan exact{bytes};
    CHECK(GmshSurfaceTopology::prepare(mesh, selection, {3, 10, 8, 2, work}, exact, {}));
    CHECK(exact.remaining_bytes() == 0);

    constexpr auto work_limit = limits.work;
    for (const auto cap :
         {GmshTopologyLimits{2, 10, 8, 2, work_limit}, GmshTopologyLimits{3, 9, 8, 2, work_limit},
          GmshTopologyLimits{3, 10, 7, 2, work_limit}, GmshTopologyLimits{3, 10, 8, 1, work_limit},
          GmshTopologyLimits{3, 10, 8, 2, work - 1}}) {
        core::ResourcePlan trial{memory};
        const auto limited = GmshSurfaceTopology::prepare(mesh, selection, cap, trial, {});
        REQUIRE_FALSE(limited);
        CHECK(limited.error().code == ResourceExceeded);
    }
    core::ResourcePlan one_less{bytes - 1};
    const auto limited = GmshSurfaceTopology::prepare(mesh, selection, limits, one_less, {});
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == ResourceExceeded);
    core::ResourcePlan no_storage{0};
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() is non-const.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    const auto cancelled =
        GmshSurfaceTopology::prepare(mesh, selection, limits, no_storage, cancel.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK(no_storage.allocated_bytes() == 0);
    CHECK(take(mesh.view()).bytes == source);
}

} // namespace
