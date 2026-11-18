#pragma once
// Generated syntax/reference contexts, not numerical expected values.
// Source SHA256: 1f4821dfaa3e75cf9e9ec2aa776e6028cf820cd0843147905249a2d86153379e
#include "case_schema_test_adapter.hpp"

#include <algorithm>
#include <string>
#include <string_view>
#include <utility>
namespace uvlas::test::si_fragment_adapter {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
struct FragmentInput {
    JsonToken context;
    JsonToken payload;
    std::string type;
    std::uint64_t components;
    io::SchemaDimension dimensions;
};

inline std::string_view family(std::string_view name) {
    return name.substr(0, name.find('.'));
}
inline void remove_identifier(JsonToken& collection, std::string_view id) {
    std::erase_if(collection.elements,
                  [&](const auto& item) { return item.field("id").text() == id; });
}
inline void branch_id(JsonToken& context, std::string_view id) {
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
inline FragmentInput input_for(std::string_view name, JsonToken payload) {
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

} // namespace uvlas::test::si_fragment_adapter
