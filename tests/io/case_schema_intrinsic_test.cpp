#include "case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string_view>
#include <uvlas/io/case_schema.hpp>

namespace {
using uvlas::test::schema_adapter::JsonToken;
using namespace uvlas::test::schema_adapter;

void expect_schema(const JsonToken& value, bool accepted) {
    using namespace uvlas;
    const auto wire = value.wire();
    const io::InputLimits syntax{
        1000000,
        256,
        100000,
        100000,
        1000000,
        1000000,
        1000000,
        {io::ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
    core::ResourcePlan source_plan{64ULL * 1024 * 1024};
    std::optional<io::InputDiagnostic> source_diagnostic;
    auto document =
        io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                io::SourceId{31}, syntax, source_plan, source_diagnostic);
    if (!document) {
        FAIL("Intrinsic controls must pass the syntax layer");
        return;
    }
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> limits{};
    limits.fill(1000000);
    limits[static_cast<std::size_t>(io::SchemaQuota::TypedBytes)] = 64ULL * 1024 * 1024;
    limits[static_cast<std::size_t>(io::SchemaQuota::ScratchBytes)] = 16ULL * 1024 * 1024;
    limits[static_cast<std::size_t>(io::SchemaQuota::Work)] = 1000000000;
    limits[static_cast<std::size_t>(io::SchemaQuota::ExactDigits)] = 4096;
    core::ResourcePlan plan{128ULL * 1024 * 1024};
    io::SchemaDiagnostic diagnostic;
    const auto result =
        io::validate_case_schema(*document, io::SchemaLimits{limits}, plan, &diagnostic);
    if (const auto detail = diagnostic.detail()) {
        INFO("reason " << static_cast<unsigned>(detail->reason));
        if (const auto path = diagnostic.path()) {
            INFO("path " << *path);
            REQUIRE(result.has_value() == accepted);
        } else {
            REQUIRE(result.has_value() == accepted);
        }
    } else {
        REQUIRE(result.has_value() == accepted);
    }
    if (!accepted) {
        CHECK(result.error().code == core::FailureCode::InvalidInput);
    }
    CHECK(document->view().has_value());
}

JsonToken& geometry(JsonToken& value) {
    return value.field("model").field("geometry").elements.at(0);
}
JsonToken& support(JsonToken& value) {
    return value.field("model").field("supports").elements.at(0);
}
void add_second_geometry(JsonToken& value) {
    auto second = geometry(value);
    second.set("id", string("2"));
    second.set("continuity", array());
    value.field("model").field("geometry").elements.push_back(std::move(second));
}
} // namespace

TEST_CASE("Schema graph bounds and polygon endpoints use exact rational values",
          "[verification][negative][IO4-01]") {
    auto control = root();
    expect_schema(control, true);
    SECTION("A term inside the declared bidegree remains valid") {
        geometry(control).set("degree", parse("[1,0]"));
        geometry(control).set(
            "coefficients",
            parse(R"([{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}])"));
        expect_schema(control, true);
        geometry(control).set("degree", parse("[0,0]"));
        expect_schema(control, false);
    }
    SECTION("An explicit closing endpoint is rejected") {
        auto& domain = geometry(control).field("domain").elements;
        domain.push_back(domain.front());
        expect_schema(control, false);
    }
    SECTION("Equal endpoints cannot evade rejection through another rational spelling") {
        auto& domain = geometry(control).field("domain").elements;
        domain.front() =
            parse(R"([{"numerator":"1","denominator":"2"},{"numerator":"0","denominator":"1"}])");
        expect_schema(control, true);
        domain.push_back(
            parse(R"([{"numerator":"2","denominator":"4"},{"numerator":"0","denominator":"3"}])"));
        expect_schema(control, false);
    }
}

TEST_CASE("Schema static analysis rejects negative prescribed speed",
          "[verification][negative][IO4-01]") {
    auto control = root("aerodynamic_steady");
    control.field("analysis").field("speed").set("value", parse("0"));
    expect_schema(control, true);
    control.field("analysis").field("speed").set("value", parse("-1e-400"));
    expect_schema(control, false);
}

TEST_CASE("Schema flutter requires branches intervals and a positive step",
          "[verification][negative][IO4-01]") {
    auto control = root("coupled_flutter");
    expect_schema(control, true);
    SECTION("Zero step") {
        control.field("analysis").set("step", parse("0"));
        expect_schema(control, false);
    }
    SECTION("Negative step") {
        control.field("analysis").set("step", parse("-1"));
        expect_schema(control, false);
    }
    SECTION("Empty speed intervals") {
        control.field("analysis").set("intervals", array());
        expect_schema(control, false);
    }
    SECTION("Empty branches cannot qualify through a dangling coverage reference") {
        control.field("analysis").set("branches", array());
        expect_schema(control, false);
    }
}

TEST_CASE("Schema regions preserve target geometry and valid boundary indices",
          "[verification][negative][IO4-01]") {
    auto control = root();
    expect_schema(control, true);
    SECTION("Point regions require the unused boundary to be exact zero") {
        support(control).field("region").set("kind", string("point"));
        expect_schema(control, true);
        support(control).field("region").set("boundary", parse("1"));
        expect_schema(control, false);
    }
    SECTION("An edge index must exist") {
        support(control).field("region").set("boundary", parse("4"));
        expect_schema(control, false);
    }
    SECTION("An existing but unrelated geometry is not the shell target geometry") {
        add_second_geometry(control);
        expect_schema(control, true);
        support(control).field("region").set("geometry", string("2"));
        expect_schema(control, false);
    }
}

TEST_CASE("Schema continuity requires a distinct graph and valid indices on both sides",
          "[verification][negative][IO4-01]") {
    auto control = root();
    add_second_geometry(control);
    geometry(control).set(
        "continuity",
        parse(
            R"([{"other":"2","kind":"smooth","own_boundary":0,"other_boundary":0,"orientation":"same"}])"));
    expect_schema(control, true);
    SECTION("Self reference") {
        geometry(control).field("continuity").elements[0].set("other", string("1"));
        expect_schema(control, false);
    }
    SECTION("Own boundary") {
        geometry(control).field("continuity").elements[0].set("own_boundary", parse("4"));
        expect_schema(control, false);
    }
    SECTION("Other boundary") {
        geometry(control).field("continuity").elements[0].set("other_boundary", parse("4"));
        expect_schema(control, false);
    }
}

TEST_CASE("Schema shell points cannot acquire independent rotational supports or couples",
          "[verification][negative][IO4-01]") {
    auto control = root();
    SECTION("Rotation support is not hidden by an unrelated scale failure") {
        support(control) = fragment("Support.constraint");
        expect_schema(control, true);
        support(control).set("components", parse("[3]"));
        // Supply the rotational scale so only the target's physical capability changes.
        control.field("numerics")
            .field("scales")
            .elements.back()
            .field("constraint")
            .set("orientation", parse(R"({"kind":"values","values":[1,1,1]})"));
        control.field("numerics")
            .field("scales")
            .elements.back()
            .field("constraint")
            .set("position", parse(R"({"kind":"none"})"));
        expect_schema(control, false);
    }
    SECTION("A shell point force does not imply an independent couple capability") {
        auto load = fragment("Load.couple");
        load.field("target").field("target").set("kind", string("shell"));
        load.set("kind", string("point_force"));
        control.field("analysis").field("loads").elements.push_back(load);
        expect_schema(control, true);
        control.field("analysis").field("loads").elements[0].set("kind", string("couple"));
        expect_schema(control, false);
    }
}

TEST_CASE("Schema aerodynamic sharing checks both boundary owners",
          "[verification][negative][IO4-01]") {
    auto control = root("aerodynamic_steady");
    auto& surfaces = control.field("model").field("aero_surfaces").elements;
    auto second = surfaces.at(0);
    second.set("id", string("2"));
    second.set("sharing", array());
    surfaces.push_back(std::move(second));
    surfaces[0].set(
        "sharing",
        parse(R"([{"other":"2","own_boundary":0,"other_boundary":0,"orientation":"same"}])"));
    expect_schema(control, true);
    SECTION("Own boundary") {
        surfaces[0].field("sharing").elements[0].set("own_boundary", parse("4"));
        expect_schema(control, false);
    }
    SECTION("Other boundary") {
        surfaces[0].field("sharing").elements[0].set("other_boundary", parse("4"));
        expect_schema(control, false);
    }
}

TEST_CASE("Schema transient intervals terminate at explicit motion and flow events",
          "[verification][negative][IO4-01]") {
    auto control = root("coupled_transient");
    control.field("analysis")
        .set(
            "intervals",
            parse(
                R"([{"start":0,"end":0.5,"step":0.5,"restart":"initial"},{"start":0.5,"end":1,"step":0.5,"restart":"event"}])"));
    SECTION("Support motion event") {
        support(control).field("motion").set("events", parse("[0.5]"));
        expect_schema(control, true);
        control.field("analysis")
            .set("intervals", parse(R"([{"start":0,"end":1,"step":1,"restart":"initial"}])"));
        expect_schema(control, false);
    }
    SECTION("Known flow event") {
        control.field("flow").set("events", parse("[0.5]"));
        expect_schema(control, true);
        control.field("analysis")
            .set("intervals", parse(R"([{"start":0,"end":1,"step":1,"restart":"initial"}])"));
        expect_schema(control, false);
    }
}

TEST_CASE("Schema external geometry retains boundary obligations for the importer",
          "[verification][IO4-01][IO4-03]") {
    for (const auto context : {"structural_static", "aerodynamic_steady"}) {
        INFO(context);
        auto control = root(context);
        expect_schema(control, true);
        auto imported = fragment("Geometry.gmsh");
        imported.set("id", string("1"));
        // Different physical files must have distinct source identities.
        imported.field("reference").field("file").set("source_id", string("3"));
        geometry(control) = std::move(imported);
        // No external file is opened here. Its boundary coverage is a prepare/import obligation.
        expect_schema(control, true);
    }
}
