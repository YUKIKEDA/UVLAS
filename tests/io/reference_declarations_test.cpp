#include "case_schema_detail.hpp"
#include "reference_declarations_fixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <set>

namespace {
using namespace uvlas;
using namespace uvlas::test::declaration_fixture;
constexpr std::array<std::string_view, 8> roles{
    "gmsh-msh4.1-ascii", "gmsh-sidecar",        "reference-geometry",    "nastran-deck",
    "nastran-sidecar",   "preparation-history", "observation-reference", "initial-history"};
constexpr std::array<std::uint64_t, 8> ids{300, 301, 302, 200, 201, 700, 900, 800};
constexpr std::array<std::string_view, 8> pointers{"/model/geometry/0/mesh",
                                                   "/model/geometry/0/sidecar",
                                                   "/model/geometry/0/reference/file",
                                                   "/model/geometry/1/deck",
                                                   "/model/geometry/1/sidecar",
                                                   "/analysis/branch/seed/history",
                                                   "/observations/0/comparison/reference/source",
                                                   "/initial/file"};
constexpr std::array<std::string_view, 8> paths{
    R"(folder\./mesh.msh)", "folder/../sidecar.json", "参照/geometry.json", "model.bdf",
    "model-sidecar.json",   "seed.history",           "observation.dat",    "initial.history"};
void provenance(const rd::DeclaredReferences& data, const JsonToken& fixture) {
    const auto doc = parsed(fixture);
    for (core::Index i = 0; i < 8; ++i) {
        const auto row = take(data.entry(i));
        auto node = take(doc.view()).root;
        auto path = row.declaration_pointer;
        while (!path.empty()) {
            path.remove_prefix(1);
            const auto end = path.find('/');
            const auto component = path.substr(0, end);
            if (auto object = node.object()) {
                const auto child = object->find(component);
                if (!child) {
                    FAIL("Expected owned value");
                    return;
                }
                node = *child;
            } else {
                const auto array = node.array();
                if (!array) {
                    FAIL("Expected owned value");
                    return;
                }
                node =
                    take(array->at(static_cast<core::Index>(std::stoull(std::string{component}))));
            }
            if (end == std::string_view::npos)
                break;
            path.remove_prefix(end);
        }
        CHECK(row.declaration_node == static_cast<std::uint64_t>(node.ordinal().value()));
    }
}
void expected(const rd::DeclaredReferences& data) {
    REQUIRE(take(data.size()) == 8);
    for (std::size_t i = 0; i < 8; ++i) {
        const auto row = take(data.entry(static_cast<core::Index>(i)));
        CHECK(row.source_id.value() == ids[i]);
        CHECK(rd::reference_role_name(row.role) == roles[i]);
        CHECK(row.original_path == paths[i]);
        CHECK(row.declaration_pointer == pointers[i]);
        const auto digit = i < 3 ? i + 3 : i < 5 ? i - 2 : i == 5 ? 6 : i == 6 ? 8 : 7;
        for (auto byte : row.declared_digest.bytes)
            CHECK(std::to_integer<unsigned>(byte) == digit * 17);
        CHECK(row.digest_spelling == std::string(64, static_cast<char>('0' + digit)));
        CHECK(row.declaration_node > 0);
    }
    CHECK(take(data.entry(0)).source_id_spelling == "0300");
}
} // namespace
TEST_CASE("All eight declarations retain dictionary order and owned original provenance",
          "[IO4-01][unit][property]") {
    auto fixture = all_roles();
    const auto source = checked(fixture);
    core::ResourcePlan plan{budget};
    auto output = take(rd::collect_case_file_references(source, limits(), plan, {}, nullptr));
    expected(output);
    provenance(output, fixture);
    CHECK(same_identity(take(output.identity()), take(source.identity())));
    reverse_keys(fixture);
    const auto reversed = checked(fixture);
    core::ResourcePlan other{budget};
    const auto second =
        take(rd::collect_case_file_references(reversed, limits(), other, {}, nullptr));
    expected(second);
    provenance(second, fixture);
    CHECK_FALSE(same_identity(take(second.identity()), take(output.identity())));
    CHECK(take(second.entry(0)).declaration_node != take(output.entry(0)).declaration_node);
    const auto before = snapshot(output);
    auto moved = std::move(output);
    // NOLINTNEXTLINE(bugprone-use-after-move): the documented moved-from state is under test.
    CHECK_FALSE(output.size());
    CHECK(snapshot(moved) == before);
    CHECK_FALSE(moved.entry(-1));
    CHECK_FALSE(moved.entry(8));
    CHECK_FALSE(limits().get(rd::DeclarationQuota::Count));
    CHECK_FALSE(limits().get(static_cast<rd::DeclarationQuota>(255)));
    CHECK_FALSE(rd::reference_role_name(static_cast<rd::ReferenceRole>(255)));
}
TEST_CASE("Declaration extraction keeps every repeated occurrence and array order",
          "[IO4-01][property]") {
    auto fixture = all_roles();
    auto duplicate = fixture.field("observations").elements[0];
    duplicate.set("id", json("\"2\""));
    duplicate.field("comparison")
        .field("reference")
        .field("source")
        .set("source_id", json("\"0900\""));
    fixture.field("observations").elements.push_back(duplicate);
    const auto source = checked(fixture);
    core::ResourcePlan plan{budget};
    const auto data = take(rd::collect_case_file_references(source, limits(), plan, {}, nullptr));
    REQUIRE(take(data.size()) == 9);
    CHECK(take(data.entry(6)).source_id.value() == 900);
    CHECK(take(data.entry(7)).source_id.value() == 900);
    CHECK(take(data.entry(7)).source_id_spelling == "0900");
    CHECK(take(data.entry(7)).declaration_pointer == "/observations/1/comparison/reference/source");
    std::swap(fixture.field("model").field("geometry").elements[0],
              fixture.field("model").field("geometry").elements[1]);
    const auto swapped = checked(fixture);
    core::ResourcePlan second{budget};
    const auto reordered =
        take(rd::collect_case_file_references(swapped, limits(), second, {}, nullptr));
    CHECK(take(reordered.entry(0)).source_id.value() == 200);
    CHECK(take(reordered.entry(2)).source_id.value() == 300);
}
TEST_CASE("Empty declaration batches and detached owners remain valid",
          "[IO4-01][negative][property]") {
    auto detached = [] {
        const auto source = checked(all_roles());
        core::ResourcePlan plan{budget};
        return take(rd::collect_case_file_references(source, limits(), plan, {}, nullptr));
    }();
    expected(detached);
    const auto source = checked(root());
    core::ResourcePlan plan{budget};
    rd::DeclarationDiagnostic slot;
    const auto output = take(rd::collect_case_file_references(source, limits(), plan, {}, &slot));
    CHECK(take(output.size()) == 0);
    CHECK_FALSE(output.entry(0));
    CHECK_FALSE(slot.detail());
    CHECK_FALSE(slot.path());
    const auto report = take(output.memory_report());
    CHECK(report.work <= report.planned_work);
    CHECK(report.owned_bytes == report.cumulative_bytes);
    auto invalid_source = checked(root());
    const auto valid_source = std::move(invalid_source);
    // NOLINTNEXTLINE(bugprone-use-after-move): the documented moved-from state is under test.
    CHECK_FALSE(rd::collect_case_file_references(invalid_source, limits(), plan, {}, &slot));
    const auto detail = slot.detail();
    if (!detail) {
        FAIL("Expected owned value");
        return;
    }
    CHECK(detail->reason == rd::DeclarationReason::InvalidOwner);
    auto moved = std::move(slot);
    // NOLINTNEXTLINE(bugprone-use-after-move): the documented moved-from state is under test.
    CHECK_FALSE(slot.detail());
    CHECK_FALSE(slot.path());
    CHECK(moved.detail().has_value());
}
TEST_CASE("Concrete FileRef dictionary inventory matches all eight implemented roles",
          "[IO4-01][verification]") {
    std::set<std::string> actual, expected_types;
    for (const auto& type : rd::schema_types())
        if (type.spelling.starts_with("FileRef<"))
            actual.emplace(type.spelling);
    for (const auto role : roles)
        expected_types.emplace("FileRef<" + std::string{role} + ">");
    CHECK(actual == expected_types);
    CHECK(actual.size() == 8);
}

TEST_CASE("External source IDs preserve zero and uint64 maximum without indexing",
          "[IO4-01][property]") {
    auto fixture = all_roles();
    auto& geometry = fixture.field("model").field("geometry").elements[0];
    geometry.field("mesh").set("source_id", json("\"0\""));
    geometry.field("sidecar").set("source_id", json("\"18446744073709551615\""));
    const auto source = checked(fixture);
    core::ResourcePlan plan{budget};
    const auto data = take(rd::collect_case_file_references(source, limits(), plan, {}, nullptr));
    CHECK(take(data.entry(0)).source_id.value() == 0);
    CHECK(take(data.entry(1)).source_id.value() == UINT64_MAX);
}
