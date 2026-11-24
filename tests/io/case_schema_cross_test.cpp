#include "case_schema_snapshot.hpp"
#include "case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <uvlas/io/case_schema.hpp>
#include <variant>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;

io::SchemaLimits cross_limits() {
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> values{};
    values.fill(1000000);
    values[static_cast<std::size_t>(io::SchemaQuota::TypedBytes)] = 64ULL * 1024 * 1024;
    values[static_cast<std::size_t>(io::SchemaQuota::ScratchBytes)] = 16ULL * 1024 * 1024;
    values[static_cast<std::size_t>(io::SchemaQuota::Work)] = 1000000000;
    values[static_cast<std::size_t>(io::SchemaQuota::ExactDigits)] = 4096;
    return io::SchemaLimits{values};
}
core::Result<io::SchemaCase> validate_cross(const JsonToken& value,
                                            io::SchemaDiagnostic* diagnostic) {
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
    core::ResourcePlan input_plan{64ULL * 1024 * 1024};
    std::optional<io::InputDiagnostic> input_diagnostic;
    const auto document =
        io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                io::SourceId{31}, syntax, input_plan, input_diagnostic);
    if (!document)
        return std::unexpected(document.error());
    core::ResourcePlan plan{64ULL * 1024 * 1024};
    // The returned case and diagnostic must outlive this source, wire, and both plans.
    return io::validate_case_schema(*document, cross_limits(), plan, diagnostic);
}
void require_cross(const JsonToken& value, bool accepted) {
    io::SchemaDiagnostic diagnostic;
    const auto result = validate_cross(value, &diagnostic);
    const auto detail = diagnostic.detail();
    if (detail) {
        INFO("reason " << static_cast<unsigned>(detail->reason));
        const auto path = diagnostic.path();
        if (path)
            INFO("path " << *path);
        CHECK(result.has_value() == accepted);
    } else
        CHECK(result.has_value() == accepted);
    if (!accepted && !result)
        CHECK(result.error().code == core::FailureCode::InvalidInput);
    if (accepted && result)
        CHECK(result->root());
}
JsonToken domain(std::size_t kind) {
    if (kind == 0)
        return parse(
            R"({"kind":"reference","coordinates":"chart","lower":[0,0,-1],"upper":[1,1,1]})");
    if (kind == 1)
        return parse(R"({"kind":"time","lower":0,"upper":1,"events":[]})");
    return parse(
        R"({"kind":"current_position_time","position_lower":[-1,-1,-1],"position_upper":[2,2,2],"time_lower":0,"time_upper":1,"events":[]})");
}
JsonToken expression_case(std::size_t consumer, std::size_t selected_domain) {
    auto value = root(consumer == 1 ? "coupled_transient" : "structural_static");
    auto expression = parse(
        R"({"domain":{},"nodes":[{"op":"constant","value":1,"dimensions":[0,0,0,0]}],"outputs":[0]})");
    expression.set("domain", domain(selected_domain));
    if (consumer == 0) {
        expression.set("outputs", parse("[0,0,0,0,0,0]"));
        auto field = parse(R"({"kind":"expression","expression":{},"error":[0,0,0,0,0,0]})");
        field.set("expression", std::move(expression));
        value.field("model").field("shells").elements[0].field("layers").elements[0].set(
            "initial_green", std::move(field));
    } else if (consumer == 1) {
        expression.field("nodes").elements[0].set("dimensions", parse("[1,0,-1,0]"));
        expression.set("outputs", parse("[0,0,0]"));
        auto field = parse(R"({"kind":"expression","expression":{}})");
        field.set("expression", std::move(expression));
        value.field("flow").set("uniform", std::move(field));
    } else {
        expression.field("nodes").elements[0].set("dimensions", parse("[-1,1,-2,0]"));
        auto load = fragment("Load.current_pressure");
        load.set("expression", std::move(expression));
        value.field("analysis").field("loads").elements.push_back(std::move(load));
    }
    return value;
}
void check_component_values(const io::schema::ComponentScale3Values& component,
                            const io::SchemaDimension& expected) {
    for (core::Index index = 0; index < 3; ++index) {
        const auto field = component.values.at(index);
        REQUIRE(field);
        const auto exact = field->number();
        CHECK(exact.dimensions == expected);
        CHECK(exact.sign == 1);
        CHECK(exact.digits == std::to_string(index + 1));
        CHECK(exact.exponent10 == 0);
    }
}
void check_constraint_dimensions(const io::SchemaCase& value, bool has_orientation) {
    const auto case_root = value.root();
    REQUIRE(case_root);
    const auto scales = case_root->numerics.record().scales;
    bool found = false;
    for (std::uint64_t index = 0; index < scales.size(); ++index) {
        const auto field = scales.at(static_cast<core::Index>(index));
        REQUIRE(field);
        const auto selected = field->record();
        const auto* policy = std::get_if<io::schema::ScalePhysicalLayoutPolicy>(&selected);
        if (policy == nullptr)
            continue;
        const auto constraint = policy->constraint.record();
        const auto* physical = std::get_if<io::schema::ConstraintScalePolicyPhysical>(&constraint);
        REQUIRE(physical != nullptr);
        const auto position = physical->position.record();
        const auto* translation = std::get_if<io::schema::ComponentScale3Values>(&position);
        REQUIRE(translation != nullptr);
        check_component_values(*translation, io::SchemaDimension{{1, 0, 0, 0}});
        const auto orientation = physical->orientation.record();
        if (has_orientation) {
            const auto* rotation = std::get_if<io::schema::ComponentScale3Values>(&orientation);
            REQUIRE(rotation != nullptr);
            check_component_values(*rotation, io::SchemaDimension{{0, 0, 0, 1}});
        } else
            CHECK(std::holds_alternative<io::schema::ComponentScale3None>(orientation));
        found = true;
    }
    CHECK(found);
}
} // namespace

TEST_CASE("All expression consumers require their own domain even for constant outputs",
          "[verification][negative][IO4-01]") {
    // Every domain is independently well formed, and each constant has the consumer's
    // correct dimensions. Only the forbidden dependence changes in a negative case.
    for (std::size_t consumer = 0; consumer < 3; ++consumer) {
        for (std::size_t selected = 0; selected < 3; ++selected) {
            INFO(consumer << "/" << selected);
            require_cross(expression_case(consumer, selected), consumer == selected);
        }
    }
}

TEST_CASE("Static observations reject history windows and retain singleton controls",
          "[verification][negative][IO4-01]") {
    for (const auto singleton : {"[0,0]", "[2,2]"}) {
        auto value = root();
        for (auto& observation : value.field("observations").elements)
            observation.field("comparison").set("window", parse(singleton));
        require_cross(value, true);
    }
    for (const auto statistic : {"instantaneous", "history", "peak", "rms"}) {
        auto value = root();
        auto& comparison = value.field("observations").elements[0].field("comparison");
        comparison.set("window", parse("[0,1]"));
        comparison.set("statistic", string(statistic));
        INFO(statistic);
        require_cross(value, false);
    }
}

TEST_CASE("Mixed physical scales retain length and angle dimensions before SI conversion",
          "[verification][property][IO4-01]") {
    auto shell = root();
    shell.field("units").set("length", string("mm"));
    shell.field("units").set("angle", string("deg"));
    shell.field("numerics")
        .field("scales")
        .elements.back()
        .field("constraint")
        .field("position")
        .set("values", parse("[1,2,3]"));
    io::SchemaDiagnostic diagnostic;
    const auto checked_shell = validate_cross(shell, &diagnostic);
    REQUIRE(checked_shell);
    check_constraint_dimensions(*checked_shell, false);
    bool exercised_orientation = false;
    for (const auto& example : uvlas::test::case_schema_corpus::physical) {
        const auto context = parse(example.context_json);
        if (!context.field("shell_ids").elements.empty())
            continue;
        auto value = physical_case(example.context_json, example.perturbation_json);
        value.field("units").set("length", string("mm"));
        value.field("units").set("angle", string("deg"));
        auto& constraint =
            value.field("numerics").field("scales").elements.back().field("constraint");
        constraint.field("position").set("values", parse("[1,2,3]"));
        constraint.field("orientation").set("values", parse("[1,2,3]"));
        const auto checked = validate_cross(value, &diagnostic);
        REQUIRE(checked);
        check_constraint_dimensions(*checked, true);
        exercised_orientation = true;
    }
    CHECK(exercised_orientation);
}

TEST_CASE("Canonical case traversal observes values beyond counts and ignores object key order",
          "[verification][property][negative][IO4-01]") {
    auto original = root();
    auto changed = original;
    changed.field("model")
        .field("shells")
        .elements[0]
        .field("layers")
        .elements[0]
        .field("initial_green")
        .set("value", parse("[0,0,0,0,0,1]"));
    auto reordered = original;
    std::reverse(reordered.members.begin(), reordered.members.end());
    io::SchemaDiagnostic diagnostic;
    const auto a = validate_cross(original, &diagnostic);
    const auto b = validate_cross(changed, &diagnostic);
    const auto c = validate_cross(reordered, &diagnostic);
    REQUIRE(a);
    REQUIRE(b);
    REQUIRE(c);
    const auto a_root = a->root();
    const auto b_root = b->root();
    const auto c_root = c->root();
    REQUIRE(a_root);
    REQUIRE(b_root);
    REQUIRE(c_root);
    const auto a_bytes = uvlas::test::schema_snapshot::serialize(*a_root);
    const auto b_bytes = uvlas::test::schema_snapshot::serialize(*b_root);
    const auto c_bytes = uvlas::test::schema_snapshot::serialize(*c_root);
    CHECK(a_root->model.record().shells.size() == b_root->model.record().shells.size());
    CHECK(a_bytes != b_bytes);
    CHECK(a_bytes == c_bytes);
    CHECK(a->source_identity()->digest != c->source_identity()->digest);
}

TEST_CASE("Flutter requires a step-map eigenproblem while preserving both valid backends",
          "[verification][negative][IO4-01]") {
    const auto flutter = root("coupled_flutter");
    const auto& requirements = flutter.field("numerics").field("eigen").field("requirements");
    auto full = flutter;
    auto spectrum = fragment("EigenPolicy.full");
    spectrum.set("requirements", requirements);
    full.field("numerics").set("eigen", std::move(spectrum));
    require_cross(full, true);

    auto arnoldi = flutter;
    auto iterative = fragment("EigenPolicy.arnoldi");
    iterative.set("requirements", requirements);
    arnoldi.field("numerics").set("eigen", std::move(iterative));
    require_cross(arnoldi, true);
    auto generalized = flutter;
    auto structural = fragment("EigenPolicy.generalized");
    structural.set("requirements", requirements);
    generalized.field("numerics").set("eigen", std::move(structural));
    require_cross(generalized, false);
    // The generalized policy remains valid for a structural-mode request outside flutter.
    auto static_case = root();
    static_case.field("numerics").set("eigen", fragment("EigenPolicy.generalized"));
    require_cross(static_case, true);
}

namespace {
void check_missing_scale(const JsonToken& value) {
    io::SchemaDiagnostic diagnostic;
    const auto invalid = validate_cross(value, &diagnostic);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
    const auto detail = diagnostic.detail();
    const bool matching_reason = detail && detail->reason == io::SchemaReason::InvalidValue;
    CHECK(matching_reason);
    CHECK(diagnostic.path() == std::string_view{"/numerics/scales"});
    const auto without_diagnostic = validate_cross(value, nullptr);
    REQUIRE_FALSE(without_diagnostic);
    CHECK(without_diagnostic.error().code == invalid.error().code);
}
} // namespace

TEST_CASE("A missing structural scale identifies the scales array",
          "[verification][negative][IO4-01]") {
    auto value = root();
    REQUIRE(validate_cross(value, nullptr));

    auto& scales = value.field("numerics").field("scales").elements;
    const auto removed = std::erase_if(
        scales, [](const auto& scale) { return scale.field("field").text() == "shell_position"; });
    REQUIRE(removed == 1);
    check_missing_scale(value);
}

TEST_CASE("Each missing explicit layout scale identifies the scales array",
          "[verification][negative][IO4-01]") {
    auto control = root();
    REQUIRE_FALSE(control.field("observations").elements.empty());
    auto& scales = control.field("numerics").field("scales").elements;
    const auto removed = std::erase_if(scales, [](const auto& scale) {
        return scale.field("kind").text() == "physical_layout_policy";
    });
    REQUIRE(removed == 1);
    for (const auto field : {"newton_correction", "constraint"}) {
        auto layout = fragment("Scale.explicit_layout");
        layout.set("field", string(field));
        layout.set("targets", array());
        scales.push_back(std::move(layout));
    }
    REQUIRE(validate_cross(control, nullptr));

    for (const auto missing : {"newton_correction", "constraint"}) {
        CAPTURE(missing);
        auto value = control;
        auto& entries = value.field("numerics").field("scales").elements;
        REQUIRE(std::erase_if(entries, [missing](const auto& scale) {
                    return scale.field("field").text() == missing;
                }) == 1);
        check_missing_scale(value);
    }
}

TEST_CASE("Canonical text traversal preserves non-ASCII bytes and embedded zero",
          "[verification][property][IO4-01]") {
    const std::string text{"\xC3\xA9\xE3\x81\x82\0", 6};
    uvlas::test::schema_snapshot::Bytes actual;
    uvlas::test::schema_snapshot::text(actual, text);
    // Eight little-endian length bytes followed by the original UTF-8 bytes and zero.
    const uvlas::test::schema_snapshot::Bytes expected{6, 0,    0,    0,    0,    0,    0,
                                                       0, 0xC3, 0xA9, 0xE3, 0x81, 0x82, 0};
    CHECK(actual == expected);
}

namespace {
JsonToken material_map_case() {
    auto value = root("coupled_flutter");
    auto& surface = value.field("model").field("aero_surfaces").elements.at(0);
    surface.set("chord", parse("[0,0.5,1]"));
    surface.field("attachment").set("material_map", parse(R"({
            "pieces":[
                {"chord_cells":[0,1],"span_cells":[0,1],"layer":"1","side":"interior","terms":[
                    {"powers":[0,0],"coefficient":[
                        {"numerator":"1","denominator":"7"},{"numerator":"0","denominator":"1"},
                        {"numerator":"-1","denominator":"11"}]},
                    {"powers":[1,0],"coefficient":[
                        {"numerator":"1","denominator":"2"},{"numerator":"0","denominator":"1"},
                        {"numerator":"0","denominator":"1"}]}]},
                {"chord_cells":[1,2],"span_cells":[0,1],"layer":"1","side":"lower","terms":[]}
            ]})"));
    return value;
}

JsonToken& material_map_pieces(JsonToken& value) {
    return value.field("model")
        .field("aero_surfaces")
        .elements.at(0)
        .field("attachment")
        .field("material_map")
        .field("pieces");
}

void expect_material_map_rejection(const JsonToken& value, std::string_view path,
                                   io::SchemaReason reason = io::SchemaReason::InvalidValue) {
    io::SchemaDiagnostic diagnostic;
    const auto rejected = validate_cross(value, &diagnostic);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    const auto detail = diagnostic.detail();
    REQUIRE(detail);
    CHECK(detail->reason == reason);
    CHECK(diagnostic.path() == path);
}
} // namespace

TEST_CASE("Structural attachment stores exact material map fields and preserves legacy syntax",
          "[verification][IO4-01][FD-06][AE-02]") {
    const auto loaded = validate_cross(material_map_case(), nullptr);
    REQUIRE(loaded);
    const auto root_case = loaded->root();
    REQUIRE(root_case);
    const auto surface = root_case->model.record().aero_surfaces.at(0);
    REQUIRE(surface);
    const auto attachment = surface->record().attachment.record();
    const auto* structural = std::get_if<io::schema::AeroAttachmentStructural>(&attachment);
    REQUIRE(structural != nullptr);
    const auto map = structural->material_map;
    REQUIRE(map);
    const auto pieces = map->record().pieces;
    REQUIRE(pieces.size() == 2);
    const auto first = pieces.at(0);
    const auto second = pieces.at(1);
    REQUIRE(first);
    REQUIRE(second);
    CHECK(first->record().side.text() == "interior");
    CHECK(second->record().side.text() == "lower");
    CHECK(first->record().layer.unsigned_value() == structural->host.unsigned_value());
    const auto term = first->record().terms.at(0);
    REQUIRE(term);
    const auto first_coefficient = term->record().coefficient.at(0);
    const auto third_coefficient = term->record().coefficient.at(2);
    REQUIRE(first_coefficient);
    REQUIRE(third_coefficient);
    CHECK(first_coefficient->record().numerator.text() == "1");
    CHECK(first_coefficient->record().denominator.text() == "7");
    CHECK(third_coefficient->record().numerator.text() == "-1");
    CHECK(third_coefficient->record().denominator.text() == "11");

    const auto legacy = validate_cross(root("coupled_flutter"), nullptr);
    REQUIRE(legacy);
    const auto legacy_root = legacy->root();
    REQUIRE(legacy_root);
    const auto legacy_surface = legacy_root->model.record().aero_surfaces.at(0);
    REQUIRE(legacy_surface);
    const auto legacy_attachment = legacy_surface->record().attachment.record();
    const auto* legacy_structural =
        std::get_if<io::schema::AeroAttachmentStructural>(&legacy_attachment);
    REQUIRE(legacy_structural != nullptr);
    CHECK_FALSE(legacy_structural->material_map);
}

TEST_CASE("Material map rejects misplaced layers gaps overlaps and repeated powers",
          "[verification][negative][IO4-01][FD-06][AE-02]") {
    auto value = material_map_case();
    auto& pieces = material_map_pieces(value).elements;
    constexpr std::string_view base{"/model/aero_surfaces/0/attachment/material_map/pieces"};
    SECTION("A missing cell is rejected even with valid individual pieces") {
        pieces.pop_back();
        expect_material_map_rejection(value, base);
    }
    SECTION("Equal total area cannot hide overlap and a hole") {
        pieces.at(1).set("chord_cells", parse("[0,1]"));
        expect_material_map_rejection(value, std::string{base} + "/1");
    }
    SECTION("Half-open ranges use cell count rather than knot count") {
        pieces.at(1).set("chord_cells", parse("[1,3]"));
        expect_material_map_rejection(value, std::string{base} + "/1/chord_cells");
    }
    SECTION("Empty map cannot attach a surface") {
        pieces.clear();
        expect_material_map_rejection(value, base, io::SchemaReason::WrongShape);
    }
    SECTION("A resolved layer under another shell cannot be used") {
        auto& shells = value.field("model").field("shells").elements;
        auto other = shells.at(0);
        other.set("id", string("2"));
        other.field("layers").elements.at(0).set("id", string("2"));
        shells.push_back(std::move(other));
        pieces.at(0).set("layer", string("2"));
        expect_material_map_rejection(value, std::string{base} + "/0/layer");
    }
    SECTION("Duplicate exponent pairs are invalid even with different coefficients") {
        pieces.at(0).field("terms").elements.at(1).set("powers", parse("[0,0]"));
        expect_material_map_rejection(value, std::string{base} + "/0/terms");
    }
    SECTION("The total degree is limited even when each exponent separately fits") {
        pieces.at(0).field("terms").elements.at(1).set("powers", parse("[2,3]"));
        expect_material_map_rejection(value, std::string{base} + "/0/terms");
    }
}
