#include "case_schema_expected_tags.hpp"
#include "case_schema_payload.hpp"
#include "case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <uvlas/io/input_document.hpp>
#include <variant>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;

struct FragmentInput {
    JsonToken context;
    JsonToken payload;
    std::string type;
    std::uint64_t components;
    io::SchemaDimension dimensions;
};

std::string_view family(std::string_view name) {
    return name.substr(0, name.find('.'));
}
void remove_identifier(JsonToken& collection, std::string_view id) {
    std::erase_if(collection.elements,
                  [&](const auto& item) { return item.field("id").text() == id; });
}
void branch_id(JsonToken& context, std::string_view id) {
    auto& analysis = context.field("analysis");
    if (analysis.contains("branch")) {
        analysis.field("branch").set("id", string(id));
    } else {
        analysis.field("branches").elements.at(0).set("id", string(id));
        analysis.field("search")
            .field("coverage")
            .set("branches", parse(std::string{"[\""} + std::string{id} + "\"]"));
    }
    if (context.contains("initial") && context.field("initial").contains("branch"))
        context.field("initial").set("branch", string(id));
}
FragmentInput input_for(std::string_view name, JsonToken payload) {
    auto context = root();
    const auto kind = family(name);
    if (name == "Analysis.flutter")
        context = root("coupled_flutter");
    else if (name == "Analysis.transient" || name == "Numerics")
        context = root("coupled_transient");
    auto& model = context.field("model");
    if (kind == "Analysis") {
        branch_id(context, "2");
        // The frozen standalone static branch explicitly names this shell target.
        if (kind == "Analysis")
            context.field("numerics")
                .field("scales")
                .elements.at(0)
                .set("targets", parse(R"([{"kind":"shell","id":"1"}])"));
    }
    if (kind == "Geometry" || kind == "Material" || kind == "Shell" || kind == "Beam" ||
        kind == "Support" || kind == "Mass" || kind == "AeroSurface" || kind == "Observation") {
        const auto collection = kind == "Geometry"      ? "geometry"
                                : kind == "Material"    ? "materials"
                                : kind == "Shell"       ? "shells"
                                : kind == "Beam"        ? "beams"
                                : kind == "Support"     ? "supports"
                                : kind == "Mass"        ? "masses"
                                : kind == "AeroSurface" ? "aero_surfaces"
                                                        : "";
        if (kind == "Observation")
            remove_identifier(context.field("observations"), payload.field("id").text());
        else
            remove_identifier(model.field(collection), payload.field("id").text());
    }
    if (kind == "Analysis") {
        auto& analysis = context.field("analysis");
        auto& branch = analysis.contains("branch") ? analysis.field("branch")
                                                   : analysis.field("branches").elements.at(0);
        if (branch.field("continuation").contains("arclength"))
            branch.field("continuation")
                .field("arclength")
                .field("state_scales")
                .elements.at(0)
                .set("targets", parse(R"([{"kind":"shell","id":"1"}])"));
    }
    if (kind == "Support") {
        auto registry_support = fragment("Support.trace_constraint");
        registry_support.set("id", string("2"));
        model.field("supports").elements.push_back(std::move(registry_support));
    }
    if (name == "Joint.adhesive") {
        auto other_geometry = model.field("geometry").elements.at(0);
        other_geometry.set("id", string("2"));
        model.field("geometry").elements.push_back(std::move(other_geometry));
        auto other_shell = model.field("shells").elements.at(0);
        other_shell.set("id", string("2"));
        other_shell.set("geometry", string("2"));
        other_shell.field("layers").elements.at(0).set("id", string("2"));
        model.field("shells").elements.push_back(std::move(other_shell));
    }
    if (kind == "Layer")
        model.field("shells").elements.at(0).field("layers").elements.at(0).set("id", string("2"));
    if (name == "Joint.rigid_attachment" || name == "Joint.weld" || name == "Load.couple") {
        model.field("beams").elements.push_back(fragment("Beam"));
    }
    if (name == "Joint.weld") {
        model.field("masses").elements.push_back(fragment("Mass.rigid"));
    }
    if (name == "Ledger.replace") {
        auto mass = fragment("Mass.particle");
        mass.set("id", string("2"));
        model.field("masses").elements.push_back(std::move(mass));
    }
    if (name == "Observation.aero" || name == "SurfaceWake") {
        model.field("aero_surfaces").elements.push_back(fragment("AeroSurface"));
    }
    if (!model.field("beams").elements.empty())
        for (const auto field : {"beam_translation", "beam_rotation", "moment_balance"})
            context.field("numerics").field("scales").elements.push_back(scale(field));
    if (name == "Joint.weld")
        for (const auto field : {"rigid_translation", "rigid_rotation"})
            context.field("numerics").field("scales").elements.push_back(scale(field));
    if (name == "Observation.aero" || name == "SurfaceWake")
        for (const auto field : {"circulation", "wake_position", "wake_strength", "aero_residual",
                                 "transport_residual"}) {
            auto item = scale(field);
            if (field == std::string_view{"circulation"} ||
                field == std::string_view{"wake_strength"} ||
                field == std::string_view{"aero_residual"})
                item.set("values", parse("[1]"));
            context.field("numerics").field("scales").elements.push_back(std::move(item));
        }
    io::SchemaDimension dimension{{0, 0, 0, 0}};
    std::string type{kind};
    if (kind == "ExprNode") {
        const auto op = payload.field("op").text();
        auto expression = parse(
            R"({"domain":{"kind":"time","lower":0,"upper":1,"events":[]},"nodes":[],"outputs":[]})");
        auto& nodes = expression.field("nodes").elements;
        if (op != "constant" && op != "variable") {
            const auto operand_dimension = (op == "sin" || op == "cos") ? "[0,0,0,1]" : "[0,0,0,0]";
            auto constant = parse(R"({"op":"constant","value":1,"dimensions":[0,0,0,0]})");
            constant.set("dimensions", parse(operand_dimension));
            nodes.push_back(constant);
            if (payload.contains("second"))
                nodes.push_back(constant);
        }
        nodes.push_back(std::move(payload));
        expression.set("outputs", parse("[" + std::to_string(nodes.size() - 1) + "]"));
        if (op == "variable")
            dimension = {{0, 0, 1, 0}};
        payload = std::move(expression);
        type = "Expression";
    }
    return {std::move(context), std::move(payload), std::move(type), 1, dimension};
}

io::SchemaLimits fragment_limits() {
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> values{};
    values.fill(1000000);
    values[static_cast<std::size_t>(io::SchemaQuota::TypedBytes)] = 64ULL * 1024 * 1024;
    values[static_cast<std::size_t>(io::SchemaQuota::ScratchBytes)] = 16ULL * 1024 * 1024;
    values[static_cast<std::size_t>(io::SchemaQuota::Work)] = 1000000000;
    values[static_cast<std::size_t>(io::SchemaQuota::ExactDigits)] = 4096;
    return io::SchemaLimits{values};
}
void expect_fragment(const FragmentInput& input, bool accepted,
                     core::FailureCode failure_code = core::FailureCode::InvalidInput) {
    auto context = input.context;
    if (input.type == "Case")
        context = input.payload;
    else
        context.set("units", input.payload);
    const auto wire = context.wire();
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
    const auto document =
        io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                io::SourceId{37}, syntax, source_plan, source_diagnostic);
    if (!document) {
        if (input.type == "Case" && !accepted) {
            CHECK(document.error().code == failure_code);
            return;
        }
        FAIL("Fragment wrapper must pass syntax without binary number conversion");
        return;
    }
    const auto root_view = document->view();
    if (!root_view) {
        FAIL("Fragment context root must exist");
        return;
    }
    const auto object = root_view->root.object();
    if (!object) {
        FAIL("Fragment context must remain an object");
        return;
    }
    const auto payload =
        input.type == "Case" ? std::optional{root_view->root} : object->find("units");
    if (!payload) {
        FAIL("Explicit fragment payload must exist");
        return;
    }
    const auto types = io::detail::schema_types();
    const auto selected =
        std::ranges::find_if(types, [&](const auto& item) { return item.spelling == input.type; });
    REQUIRE(selected != types.end());
    const auto type = static_cast<io::detail::FieldType>(selected - types.begin());
    core::ResourcePlan plan{128ULL * 1024 * 1024};
    io::SchemaDiagnostic diagnostic;
    const auto checked = io::detail::validate_schema_payload(
        *document, {*payload, type, input.components, input.dimensions, input.type != "Case"},
        fragment_limits(), plan, &diagnostic);
    if (const auto detail = diagnostic.detail()) {
        INFO("reason " << static_cast<unsigned>(detail->reason));
        if (const auto path = diagnostic.path()) {
            INFO("path " << *path);
            CHECK(checked.has_value() == accepted);
        } else {
            CHECK(checked.has_value() == accepted);
        }
    } else {
        CHECK(checked.has_value() == accepted);
    }
    if (!accepted && !checked)
        CHECK(checked.error().code == failure_code);
}
} // namespace

TEST_CASE("All 144 frozen fragments use the common typed lowering with explicit context",
          "[verification][negative][IO4-01]") {
    std::uint64_t visited = 0;
    std::uint64_t accepted = 0;
    std::uint64_t refused = 0;
    for (const auto& example : uvlas::test::case_schema_corpus::fragments) {
        INFO(example.name);
        auto input = input_for(example.name, parse(example.json));
        const bool original_valid = !example.superseded && example.name != "Geometry.gmsh";
        expect_fragment(input, original_valid);
        ++visited;
        if (original_valid)
            ++accepted;
        else
            ++refused;
        if (example.name == "Geometry.gmsh") {
            // The published raw example deliberately remains unchanged and rejected.
            input.payload.field("reference").field("file").set("source_id", string("3"));
            expect_fragment(input, true);
        }
    }
    CHECK(visited == 144);
    CHECK(accepted == 142);
    CHECK(refused == 2);
}

namespace {
using PathStep = std::variant<std::string, std::size_t>;
using Path = std::vector<PathStep>;

JsonToken& locate(JsonToken& root_value, const Path& path) {
    auto* current = &root_value;
    for (const auto& step : path) {
        if (const auto* key = std::get_if<std::string>(&step))
            current = &current->field(*key);
        else
            current = &current->elements.at(std::get<std::size_t>(step));
    }
    return *current;
}
bool alternative(std::string_view choices, std::string_view value) {
    while (!choices.empty()) {
        const auto comma = choices.find(',');
        if (choices.substr(0, comma) == value)
            return true;
        if (comma == std::string_view::npos)
            return false;
        choices.remove_prefix(comma + 1);
    }
    return false;
}
std::optional<std::string> selected_tag(const JsonToken& value, std::string_view selector) {
    const auto* current = &value;
    while (!selector.empty()) {
        const auto dot = selector.find('.');
        const auto key = selector.substr(0, dot);
        if (!current->contains(key))
            return std::nullopt;
        current = &current->field(key);
        if (dot == std::string_view::npos)
            return current->text();
        selector.remove_prefix(dot + 1);
    }
    return std::nullopt;
}
std::string_view declared_choices(std::string_view spelling) {
    const auto body = spelling.substr(4, spelling.size() - 5);
    // FD-09 covers refers to the twelve explicitly declared ErrorGroup kinds.
    if (body == "同じ集合")
        return "geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,"
               "algebra,newton,eigen_range,search_width,load_mapping";
    return body;
}
io::detail::FieldType type_named(std::string_view spelling) {
    const auto types = io::detail::schema_types();
    const auto found =
        std::ranges::find_if(types, [&](const auto& type) { return type.spelling == spelling; });
    if (found == types.end())
        throw std::runtime_error("Explicit fragment type is missing");
    return static_cast<io::detail::FieldType>(found - types.begin());
}
struct SeriesCounts {
    std::set<std::string> records;
    std::set<std::string> fields;
    std::set<std::string> shapes;
};

void dictionary_series(const FragmentInput& control, io::detail::FieldType field_type,
                       const JsonToken& value, const Path& path, std::uint64_t components,
                       SeriesCounts& counts) {
    const auto& type = io::detail::schema_types()[static_cast<std::size_t>(field_type)];
    if (type.components != 0)
        components = type.components;
    if (type.spelling.starts_with("PolyTerm<"))
        components = 3;
    if (type.kind == io::detail::FieldKind::Array) {
        auto extent = type.extent;
        if (extent == 0 && !type.spelling.starts_with("List<"))
            extent = components;
        if (extent != 0 && counts.shapes.insert(std::string{type.spelling}).second) {
            INFO("fixed shape " << type.spelling);
            auto changed = control;
            auto& array_value = locate(changed.payload, path);
            if (array_value.elements.empty())
                array_value.elements.push_back(parse("0"));
            else
                array_value.elements.pop_back();
            expect_fragment(changed, false);
        }
        const auto element = type.element;
        if (!element) {
            FAIL("Array field is missing its declared element type");
            return;
        }
        const auto element_type = element.value();
        for (std::size_t index = 0; index < value.elements.size(); ++index) {
            auto child_path = path;
            child_path.emplace_back(index);
            dictionary_series(control, element_type, value.elements[index], child_path, components,
                              counts);
        }
        return;
    }
    if (type.kind != io::detail::FieldKind::Object)
        return;
    const auto records = io::detail::schema_records();
    const auto selected = std::ranges::find_if(records, [&](const auto& record) {
        if (record.family != type.family)
            return false;
        if (record.selector.empty())
            return true;
        const auto tag = selected_tag(value, record.selector);
        return tag.has_value() && alternative(record.alternatives, tag.value());
    });
    REQUIRE(selected != records.end());
    if (counts.records.insert(std::string{selected->name}).second) {
        INFO("unknown key in " << selected->name);
        auto changed = control;
        locate(changed.payload, path).set("__unknown_schema_field__", parse("0"));
        expect_fragment(changed, false);
    }
    for (const auto& field : selected->fields) {
        if (!value.contains(field.name))
            continue;
        auto child_path = path;
        child_path.emplace_back(std::string{field.name});
        const auto identity = std::string{selected->name} + ":" + std::string{field.name};
        if (counts.fields.insert(identity).second) {
            INFO(identity);
            auto wrong_kind = control;
            locate(wrong_kind.payload, child_path) = parse("null");
            expect_fragment(wrong_kind, false);
            if (io::detail::schema_types()[static_cast<std::size_t>(field.type)].kind ==
                io::detail::FieldKind::Enumeration) {
                auto wrong_tag = control;
                locate(wrong_tag.payload, child_path) = string("__unsupported_schema_tag__");
                expect_fragment(wrong_tag, false, core::FailureCode::Unsupported);
            }
            if (!field.optional) {
                auto missing = control;
                locate(missing.payload, path).erase(field.name);
                expect_fragment(missing, false);
            }
        }
        dictionary_series(control, field.type, value.field(field.name), child_path, components,
                          counts);
    }
}
} // namespace

TEST_CASE("Every closed dictionary enumeration accepts its declared tags and rejects another",
          "[verification][negative][IO4-01]") {
    std::uint64_t enumeration_types = 0;
    std::set<std::string> actual_tag_sets;
    std::uint64_t declared_tags = 0;
    for (const auto& type : io::detail::schema_types()) {
        if (type.kind != io::detail::FieldKind::Enumeration)
            continue;
        INFO(type.spelling);
        ++enumeration_types;
        std::vector<std::string> sorted;
        auto raw_choices = declared_choices(type.spelling);
        while (!raw_choices.empty()) {
            const auto comma = raw_choices.find(',');
            sorted.emplace_back(raw_choices.substr(0, comma));
            if (comma == std::string_view::npos)
                break;
            raw_choices.remove_prefix(comma + 1);
        }
        std::ranges::sort(sorted);
        std::string normalized;
        for (const auto& tag : sorted) {
            if (!normalized.empty())
                normalized += ',';
            normalized += tag;
        }
        actual_tag_sets.insert(std::move(normalized));
        auto choices = declared_choices(type.spelling);
        while (!choices.empty()) {
            const auto comma = choices.find(',');
            const auto tag = choices.substr(0, comma);
            const FragmentInput control{
                root(), string(tag), std::string{type.spelling}, 1, {{0, 0, 0, 0}}};
            expect_fragment(control, true);
            ++declared_tags;
            if (comma == std::string_view::npos)
                break;
            choices.remove_prefix(comma + 1);
        }
        FragmentInput rejected{root(),
                               string("__unsupported_schema_tag__"),
                               std::string{type.spelling},
                               1,
                               {{0, 0, 0, 0}}};
        expect_fragment(rejected, false, core::FailureCode::Unsupported);
        rejected.payload = parse("null");
        expect_fragment(rejected, false);
    }
    std::set<std::string> expected;
    for (const auto tags : uvlas::test::schema_inventory::expected_tag_sets)
        expected.emplace(tags);
    CHECK(actual_tag_sets == expected);
    CHECK(enumeration_types > 100);
    CHECK(declared_tags > enumeration_types);
}

TEST_CASE("Frozen controls exercise recursive dictionary key kind and fixed shape series",
          "[verification][negative][IO4-01]") {
    SeriesCounts counts;
    for (const auto& example : uvlas::test::case_schema_corpus::fragments) {
        INFO(example.name);
        if (example.superseded)
            continue;
        auto control = input_for(example.name, parse(example.json));
        if (example.name == "Geometry.gmsh")
            control.payload.field("reference").field("file").set("source_id", string("3"));
        dictionary_series(control, type_named(control.type), control.payload, {},
                          control.components, counts);
    }
    for (const auto& example : uvlas::test::case_schema_corpus::physical) {
        INFO(example.name);
        FragmentInput control{root(),
                              physical_case(example.context_json, example.perturbation_json),
                              "Case",
                              1,
                              {{0, 0, 0, 0}}};
        const auto physical_context = parse(example.context_json);
        if (!physical_context.field("beam_ids").elements.empty()) {
            auto point = parse(
                R"({"point":{"target":{"kind":"beam","id":"1"},"coordinates":[0,0,0],"side":"interior"},"value":[0,0,0],"tolerance":[0,0,0]})");
            point.field("point").erase("side");
            point.field("point").field("target").set(
                "id", physical_context.field("beam_ids").elements.at(0));
            control.payload.field("initial")
                .field("perturbation")
                .field("displacements")
                .elements.push_back(std::move(point));
        }
        expect_fragment(control, true);
        dictionary_series(control, type_named("Case"), control.payload, {}, 1, counts);
    }
    {
        auto wake = input_for(
            "SurfaceWake",
            parse(
                R"({"surface":"1","rows":1,"points":[[0,0,0],[1,0,0]],"coefficients":[0],"row_ids":["1"],"point_ids":["1","2"],"band_ids":["1"],"birth_times":[0]})"));
        expect_fragment(wake, true);
        dictionary_series(wake, type_named("SurfaceWake"), wake.payload, {}, 1, counts);
        auto cluster =
            input_for("EigenCluster", parse(R"({"members":[0],"expected_dimension":1})"));
        expect_fragment(cluster, true);
        dictionary_series(cluster, type_named("EigenCluster"), cluster.payload, {}, 1, counts);
    }
    {
        INFO("hand-written nonempty polynomial continuity and aerodynamic sharing");
        auto value = root("coupled_flutter");
        value.set("restart", fragment("Restart"));
        value.field("model")
            .field("aero_surfaces")
            .elements.at(0)
            .field("attachment")
            .set(
                "material_map",
                parse(
                    R"({"pieces":[{"chord_cells":[0,1],"span_cells":[0,1],"layer":"1","side":"interior","terms":[{"powers":[1,0],"coefficient":[{"numerator":"1","denominator":"1"},{"numerator":"0","denominator":"1"},{"numerator":"0","denominator":"1"}]}]}]})"));
        value.set(
            "provenance",
            parse(
                R"({"description":"explicit dictionary control","references":[],"assumptions":[]})"));
        auto& graph = value.field("model").field("geometry").elements.at(0);
        graph.set("degree", parse("[1,0]"));
        graph.set("coefficients",
                  parse(R"([{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}])"));
        auto other_graph = graph;
        other_graph.set("id", string("2"));
        value.field("model").field("geometry").elements.push_back(std::move(other_graph));
        value.field("model")
            .field("geometry")
            .elements.at(0)
            .set(
                "continuity",
                parse(
                    R"([{"other":"2","kind":"smooth","own_boundary":0,"other_boundary":0,"orientation":"same"}])"));
        auto& surfaces = value.field("model").field("aero_surfaces").elements;
        auto other_surface = surfaces.at(0);
        other_surface.set("id", string("2"));
        other_surface.set("geometry", string("2"));
        surfaces.push_back(std::move(other_surface));
        surfaces.at(0).set(
            "sharing",
            parse(R"([{"other":"2","own_boundary":0,"other_boundary":0,"orientation":"same"}])"));
        const FragmentInput control{root(), std::move(value), "Case", 1, {{0, 0, 0, 0}}};
        expect_fragment(control, true);
        dictionary_series(control, type_named("Case"), control.payload, {}, 1, counts);
    }
    for (const auto& record : io::detail::schema_records()) {
        INFO("record coverage " << record.name);
        CHECK(counts.records.contains(std::string{record.name}));
    }
    CHECK(counts.records.size() == 194);
    CHECK(counts.fields.size() == 815);
    CHECK(counts.fields.size() > 600);
    CHECK(counts.shapes.size() > 20);
}
