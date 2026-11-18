#include "case_schema_snapshot.hpp"
#include "case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

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
                                            io::SchemaDiagnostic& diagnostic) {
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
    return io::validate_case_schema(*document, cross_limits(), plan, &diagnostic);
}
void require_cross(const JsonToken& value, bool accepted) {
    io::SchemaDiagnostic diagnostic;
    const auto result = validate_cross(value, diagnostic);
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
    const auto checked_shell = validate_cross(shell, diagnostic);
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
        const auto checked = validate_cross(value, diagnostic);
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
    const auto a = validate_cross(original, diagnostic);
    const auto b = validate_cross(changed, diagnostic);
    const auto c = validate_cross(reordered, diagnostic);
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
