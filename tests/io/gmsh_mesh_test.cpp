#include "../support/allocation_audit.hpp"
#include "gmsh_mesh.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using enum core::FailureCode;

constexpr std::string_view mesh = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
3
0 7 "mass point"
1 7 "clamp line"
2 7 "wing skin"
$EndPhysicalNames
$Entities
2 1 1 0
1 0 0 0 1 7
9 1 1 0 0
3 0 0 0 1 1 0 1 7 2 1 -9
17 0 0 0 1 1 0 1 7 1 -3
$EndEntities
$Nodes
3 4 11 9007199254740993
2 17 1 2
9007199254740993 400
1 1 0 0.75 -0.25
0 1 0 0.125 +.250
0 1 0 1
11
-0 0 0
0 9 0 1
42
1. 0 0
$EndNodes
$Elements
4 4 50 990
2 17 3 1
990 11 42 9007199254740993 400
0 1 15 1
50 11
1 3 1 1
75 11 42
2 17 2 1
84 11 42 400
$EndElements
)";

constexpr GmshLimits limits{mesh.size(), 3, 4, 6, 4, 4, 10000};
constexpr auto huge_node = std::uint64_t{9007199254740993};

consteval SourceDigest fixture_digest() {
    constexpr std::string_view hex = "0123456789abcdef";
    constexpr std::string_view hash =
        "1e75f5772cc5680f0111656b39e12937afb277221c7ede759c059f22e2606f13";
    SourceDigest digest{};
    for (std::size_t i = 0; i < digest.bytes.size(); ++i)
        digest.bytes[i] =
            static_cast<std::byte>(hex.find(hash[2 * i]) * 16 + hex.find(hash[2 * i + 1]));
    return digest;
}

class Reader final : public ReferenceFileReader {
public:
    explicit Reader(std::string_view text) noexcept
        : input_(std::as_bytes(std::span{text.data(), text.size()})) {}

    core::Result<ReferenceFileInfo> info() noexcept override {
        return ReferenceFileInfo{{0, {}}, input_.size(), {}, true};
    }

    core::Result<std::size_t> read(std::span<std::byte> bytes) noexcept override {
        const auto size = std::min(input_.size() - offset_, bytes.size());
        std::ranges::copy(input_.subspan(offset_, size), bytes.begin());
        offset_ += size;
        return size;
    }

private:
    std::span<const std::byte> input_;
    std::size_t offset_ = 0;
};

AcquiredRootSnapshot snapshot(std::string_view text, const SourceDigest* expected = nullptr) {
    Reader reader{text};
    core::ResourcePlan plan{text.size()};
    const ReferenceFileLimits read_limits{text.size(), text.size() + 1, 3, text.size()};
    auto result = expected == nullptr
                      ? read_manifest_file(reader, read_limits, plan, {}, nullptr)
                      : read_reference_file(reader, *expected, read_limits, plan, {}, nullptr);
    REQUIRE(result);
    constexpr RootSnapshotLocation location{
        {{SourceId{19}, 27, {}}, 1, 1}, SourceId{81}, ReferenceRole::GmshMesh, 5};
    return {location, std::move(*result), {}};
}

template <class T>
T take_view(std::optional<T> value) {
    if (!value)
        throw std::logic_error("Missing mesh or source view");
    return *value;
}
std::string edited(std::string_view from, std::string_view to) {
    auto copy = std::string{mesh};
    const auto offset = copy.find(from);
    REQUIRE(offset != std::string::npos);
    copy.replace(offset, from.size(), to);
    return copy;
}

} // namespace

TEST_CASE("Gmsh prepares the adopted mesh types and keeps original tags, order and decimal words",
          "[IO4-04][IO4-01][N4-02][verification][property]") {
    using namespace uvlas::test_support;
    const auto value = [] {
        constexpr auto digest = fixture_digest();
        const auto input = snapshot(mesh, &digest);
        const auto before = take_view(input.file.view());
        core::ResourcePlan memory{10000};
        AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
        auto output = GmshMesh::prepare(input, limits, memory, {});
        CHECK(audit.finish().cpp_allocation_attempts == 7);
        REQUIRE(output);
        CHECK(std::ranges::equal(take_view(input.file.view()).bytes, before.bytes));
        CHECK(take_view(output->view()).digest == before.digest);
        return std::move(*output);
    }();
    const auto view = take_view(value.view());
    CHECK(view.bytes == mesh);
    CHECK(view.location.parent.source.source == SourceId{19});
    CHECK(view.location.source == SourceId{81});
    CHECK(view.location.declaration_node == 5);
    CHECK(view.entities_declared);
    CHECK(view.physical_names.size() == 3);
    CHECK(view.physical_names[0].group == GmshKey{0, 7});
    CHECK(view.physical_names[1].group == GmshKey{1, 7});
    CHECK(view.physical_names[2].group == GmshKey{2, 7});
    CHECK(view.text(view.physical_names[2].name) == "wing skin");
    REQUIRE(view.entities.size() == 4);
    const auto& surface = view.entities[3];
    CHECK(surface.key == GmshKey{2, 17});
    REQUIRE(surface.physical_count == 1);
    CHECK(view.entity_links[surface.first_physical].key == GmshKey{2, 7});
    REQUIRE(surface.boundary_count == 1);
    CHECK(view.entity_links[surface.first_boundary].key == GmshKey{1, 3});
    CHECK(view.entity_links[surface.first_boundary].reversed);
    REQUIRE(view.nodes.size() == 4);
    CHECK(view.nodes[0].tag == 11);
    CHECK(view.text(view.nodes[0].coordinates[0]) == "-0");
    CHECK(view.nodes[1].tag == 42);
    CHECK(view.text(view.nodes[1].coordinates[0]) == "1.");
    CHECK(view.nodes[2].tag == 400);
    CHECK(view.text(view.nodes[2].parameters[1]) == "+.250");
    CHECK(view.nodes[3].tag == huge_node);
    CHECK(view.nodes[3].parametric);
    CHECK(view.text(view.nodes[3].parameters[1]) == "-0.25");
    REQUIRE(view.elements.size() == 4);
    CHECK(view.elements[0].kind == GmshElementKind::Point);
    CHECK(view.elements[1].kind == GmshElementKind::Line);
    CHECK(view.elements[2].kind == GmshElementKind::Triangle);
    CHECK(view.elements[3].kind == GmshElementKind::Quadrangle);
    constexpr std::array quad{std::uint64_t{11}, std::uint64_t{42}, huge_node, std::uint64_t{400}};
    CHECK(std::ranges::equal(view.elements[3].connectivity(), quad));
}

TEST_CASE("Gmsh accepts sparse repeated sections and does not invent optional groups",
          "[IO4-04][verification]") {
    constexpr std::string_view sparse = R"($MeshFormat
4.1 0 4
$EndMeshFormat
$Nodes
1 1 42 42
1 17 0 1
42
0 0 0
$EndNodes
$Nodes
1 1 11 11
1 17 0 1
11
1e-500 0 0
$EndNodes
$Elements
1 1 19 19
1 17 1 1
19 42 11
$EndElements
)";
    const auto input = snapshot(sparse);
    core::ResourcePlan plan{10000};
    const GmshLimits minimal{sparse.size(), 0, 0, 0, 2, 1, 10000};
    const auto parsed = GmshMesh::prepare(input, minimal, plan, {});
    REQUIRE(parsed);
    const auto view = take_view(parsed->view());
    CHECK_FALSE(view.entities_declared);
    CHECK(view.physical_names.empty());
    CHECK(view.entity_links.empty());
    REQUIRE(view.nodes.size() == 2);
    CHECK(view.nodes[0].tag == 11);
    CHECK(view.text(view.nodes[0].coordinates[0]) == "1e-500");
    CHECK(view.nodes[1].tag == 42);
}

TEST_CASE("Gmsh counts empty entity blocks separately from stored nodes and elements",
          "[IO4-04][N4-02][verification][property]") {
    constexpr std::string_view sparse = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
1
0 2147483647 "at the tag limit"
$EndPhysicalNames
$Entities
1 1 0 0
2147483647 0 0 0 1 2147483647
2147483647 0 0 0 1 1 1 0 1 -2147483647
$EndEntities
$Nodes
2 1 18446744073709551615 18446744073709551615
1 2147483647 0 0
0 2147483647 0 1
18446744073709551615
0 0 0
$EndNodes
$Elements
2 1 18446744073709551615 18446744073709551615
1 2147483647 1 0
0 2147483647 15 1
18446744073709551615 18446744073709551615
$EndElements
)";
    const auto input = snapshot(sparse);
    const GmshLimits minimal{sparse.size(), 1, 2, 2, 1, 1, 10000};
    core::ResourcePlan plan{10000};
    const auto parsed = GmshMesh::prepare(input, minimal, plan, {});
    REQUIRE(parsed);
    const auto view = take_view(parsed->view());
    REQUIRE(view.nodes.size() == 1);
    REQUIRE(view.elements.size() == 1);
    CHECK(view.nodes[0].tag == std::numeric_limits<std::uint64_t>::max());
    CHECK(view.elements[0].tag == view.nodes[0].tag);
    REQUIRE(view.elements[0].connectivity().size() == 1);
    CHECK(view.elements[0].connectivity()[0] == view.nodes[0].tag);
    REQUIRE(view.physical_names.size() == 1);
    CHECK(view.physical_names[0].group.tag == std::numeric_limits<std::int32_t>::max());
    REQUIRE(view.entity_links.size() == 2);
    CHECK(view.entity_links[1].reversed);
    CHECK(view.entity_links[1].key.tag == std::numeric_limits<std::int32_t>::max());

    auto out_of_range = std::string{sparse};
    const auto offset = out_of_range.find("2147483647");
    REQUIRE(offset != out_of_range.npos);
    out_of_range.replace(offset, 10, "2147483648");
    const auto bad = snapshot(out_of_range);
    core::ResourcePlan rejected_plan{10000};
    std::size_t rejected_offset = 0;
    const auto rejected = GmshMesh::prepare(bad, minimal, rejected_plan, {}, &rejected_offset);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == InvalidInput);
    CHECK(rejected_offset == offset);
    CHECK(rejected_plan.allocated_bytes() == 0);

    // Optional sections can occur after the node and element blocks.
    const auto begin = sparse.find("$Entities\n");
    const auto end = sparse.find("$Nodes\n");
    REQUIRE(begin < end);
    auto late_entities = std::string{sparse.substr(0, begin)};
    late_entities.append(sparse.substr(end));
    late_entities.append(sparse.substr(begin, end - begin));
    const auto late_source = snapshot(late_entities);
    core::ResourcePlan late_plan{10000};
    CHECK(GmshMesh::prepare(late_source, minimal, late_plan, {}));

    for (const auto header : {"1 2147483647 0 0", "1 2147483647 1 0"}) {
        CAPTURE(header);
        auto dangling = late_entities;
        const auto start = dangling.find(header);
        REQUIRE(start != dangling.npos);
        const auto tag_offset = start + 2;
        dangling.replace(tag_offset, 10, "99");
        const auto dangling_source = snapshot(dangling);
        core::ResourcePlan dangling_plan{10000};
        std::size_t diagnostic_offset = 0;
        const auto attempt =
            GmshMesh::prepare(dangling_source, minimal, dangling_plan, {}, &diagnostic_offset);
        CHECK_FALSE(attempt);
        if (!attempt) {
            CHECK(attempt.error().code == InvalidInput);
            CHECK(diagnostic_offset == tag_offset);
        }
    }
}

TEST_CASE("Gmsh rejects overflowing int tags in every definition and reference",
          "[IO4-04][negative][verification]") {
    constexpr std::array replacements{
        std::pair{"2 7 \"wing skin\"", "2 2147483648 \"wing skin\""},
        std::pair{"3 0 0 0 1 1 0", "2147483648 0 0 0 1 1 0"},
        std::pair{"0 1 1 0 1 7", "0 1 1 0 1 2147483648"},
        std::pair{"1 7 2 1 -9", "1 7 2 1 -2147483648"},
        std::pair{"2 17 1 2", "2 2147483648 1 2"},
        std::pair{"2 17 3 1", "2 2147483648 3 1"},
    };
    for (const auto& [before, after] : replacements) {
        CAPTURE(before, after);
        const auto bytes = edited(before, after);
        const auto input = snapshot(bytes);
        auto quota = limits;
        quota.file_bytes = bytes.size();
        core::ResourcePlan plan{10000};
        std::size_t offset = 0;
        const auto result = GmshMesh::prepare(input, quota, plan, {}, &offset);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == InvalidInput);
        CHECK(plan.allocated_bytes() == 0);
        const auto digits = bytes.find("2147483648");
        REQUIRE(digits != std::string::npos);
        CHECK(offset == (bytes[digits - 1] == '-' ? digits - 1 : digits));
    }
}

TEST_CASE("Gmsh distinguishes inconsistent mesh topology from unsupported sections and kinds",
          "[IO4-04][negative][verification]") {
    struct Case {
        std::string_view from, to;
        core::FailureCode code;
    };
    constexpr std::array cases{
        Case{"4.1 0 8", "4.0 0 8", Unsupported},
        Case{"4.1 0 8", "4.1 1 8", Unsupported},
        Case{"$PhysicalNames", "$Periodic", Unsupported},
        Case{"2 17 3 1", "2 17 9 1", Unsupported},
        Case{"990 11 42", "990 12 42", InvalidInput},
        Case{"84 11 42", "75 11 42", InvalidInput},
        Case{"9007199254740993 400\n1", "9007199254740993 42\n1", InvalidInput},
        Case{"2 17 1 2", "2 18 1 2", InvalidInput},
        Case{"1 7 1 -3", "1 7 1 -4", InvalidInput},
        Case{"2 17 3 1", "1 17 3 1", InvalidInput},
        Case{"4 4 50 990", "4 3 50 990", InvalidInput},
        Case{"+.250", "1e+", InvalidInput},
        Case{"$EndNodes", "$EndElements", InvalidInput},
        Case{"11 9007199254740993\n2", "12 9007199254740993\n2", InvalidInput},
    };
    for (const auto& trial : cases) {
        CAPTURE(trial.from, trial.to);
        const auto bytes = edited(trial.from, trial.to);
        const auto input = snapshot(bytes);
        auto quota = limits;
        quota.file_bytes = bytes.size();
        core::ResourcePlan plan{10000};
        std::size_t offset = 0;
        const auto result = GmshMesh::prepare(input, quota, plan, {}, &offset);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == trial.code);
        CHECK(offset <= bytes.size());
        CHECK(take_view(input.file.view()).bytes.size() == bytes.size());
    }
}

TEST_CASE("Gmsh mesh resource bounds include the saved bytes and both syntax and link work",
          "[IO4-04][N4-02][negative][verification]") {
    const auto input = snapshot(mesh);
    core::ResourcePlan memory{10000};
    const auto parsed = GmshMesh::prepare(input, limits, memory, {});
    REQUIRE(parsed);
    const auto view = take_view(parsed->view());
    const auto bytes = memory.allocated_bytes();
    REQUIRE(bytes > 0);
    REQUIRE(view.work > 3 * mesh.size());

    auto quota = limits;
    quota.work = view.work;
    core::ResourcePlan exact{bytes};
    CHECK(GmshMesh::prepare(input, quota, exact, {}));
    CHECK(exact.allocated_bytes() == bytes);

    core::ResourcePlan short_memory{bytes - 1};
    const auto rejected_memory = GmshMesh::prepare(input, quota, short_memory, {});
    REQUIRE_FALSE(rejected_memory);
    CHECK(rejected_memory.error().code == ResourceExceeded);
    CHECK(short_memory.allocated_bytes() == 0);

    --quota.work;
    core::ResourcePlan short_work{bytes};
    const auto rejected_work = GmshMesh::prepare(input, quota, short_work, {});
    REQUIRE_FALSE(rejected_work);
    CHECK(rejected_work.error().code == ResourceExceeded);

    quota = limits;
    --quota.entity_links;
    core::ResourcePlan short_count{bytes};
    const auto rejected_count = GmshMesh::prepare(input, quota, short_count, {});
    REQUIRE_FALSE(rejected_count);
    CHECK(rejected_count.error().code == ResourceExceeded);
    CHECK(short_count.allocated_bytes() == 0);

    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() requires a mutable source.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    core::ResourcePlan cancelled{0};
    const auto aborted = GmshMesh::prepare(input, quota, cancelled, stop.get_token());
    REQUIRE_FALSE(aborted);
    CHECK(aborted.error().code == Cancelled);
    CHECK(cancelled.allocated_bytes() == 0);

    auto wrong_role = snapshot(mesh);
    wrong_role.location.role = ReferenceRole::GmshSidecar;
    core::ResourcePlan role_storage{bytes};
    const auto not_mesh = GmshMesh::prepare(wrong_role, limits, role_storage, {});
    REQUIRE_FALSE(not_mesh);
    CHECK(not_mesh.error().code == InvalidInput);
    CHECK(role_storage.allocated_bytes() == 0);
}
