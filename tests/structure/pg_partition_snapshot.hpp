#pragma once
#include "exact_polynomial_detail.hpp"
#include "pg_partition_test_support.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <type_traits>
namespace partition_test {
using Json = nlohmann::json;
/// Serialize the actual mathematical result, excluding padding, clocks and resource forecasts.
[[nodiscard]] inline std::string mathematical_payload(const PgIntegrationPartition& result,
                                                      bool include_source = true) {
    Json out = Json::array();
    const auto value = [](const ExactRational& x) {
        const auto& raw = detail::ExactFactory::value(x).backend();
        return Json::array(
            {raw.num().str(0, std::ios_base::dec), raw.denom().str(0, std::ios_base::dec)});
    };
    for (const auto& x : take(result.vertices()))
        out.push_back(Json::array({"vertex", value(x.s), value(x.t)}));
    for (const auto& x : take(result.activity_map()))
        out.push_back(Json::array(
            {"cell", x.cell.first, x.cell.second, value(x.area), x.components, x.active}));
    for (const auto& x : take(result.pieces()))
        out.push_back(
            Json::array({"piece", x.cell.first, x.cell.second, x.component, x.thickness_triangle,
                         x.first_vertex, x.vertex_count, value(x.area)}));
    for (const auto x : take(result.piece_vertex_indices()))
        out.push_back(Json::array({"index", x}));
    for (const auto& x : take(result.edges()))
        out.push_back(
            Json::array({"edge", x.vertices[0], x.vertices[1],
                         x.left.has_value() ? Json{x.left.value_or(core::Index{-1})} : Json{},
                         x.right.has_value() ? Json{x.right.value_or(core::Index{-1})} : Json{},
                         x.first_source, x.source_count, x.external}));
    for (const auto x : take(result.edge_source_indices()))
        out.push_back(Json::array({"source-index", x}));
    for (const auto& x : take(result.triangles()))
        out.push_back(Json::array({"triangle", x.piece, x.vertices[0], x.vertices[1], x.vertices[2],
                                   value(x.twice_area)}));
    for (const auto& x : take(result.volumes()))
        out.push_back(Json::array({"volume", x.triangle, x.layer}));
    for (const auto& x : take(result.horizontal_faces()))
        out.push_back(
            Json::array({"horizontal", x.triangle, x.layer, static_cast<unsigned>(x.side)}));
    for (const auto& x : take(result.lateral_faces()))
        out.push_back(Json::array({"lateral", x.edge, x.layer, x.canonical_direction}));
    for (const auto& x : take(result.edge_traces()))
        out.push_back(Json::array({"trace", x.edge, x.source, x.piece, x.canonical_direction}));
    for (const auto& x : take(result.contacts()))
        out.push_back(Json::array(
            {"contact", x.cell.first, x.cell.second, x.source, x.vertices[0], x.vertices[1]}));
    if (!include_source)
        return out.dump();
    const auto graph = take(result.graph());
    out.push_back(Json::array({"graph", value(graph.x0), value(graph.width), value(graph.skew),
                               graph.n, static_cast<unsigned>(graph.convention)}));
    for (const auto& term : graph.height.terms())
        out.push_back(Json::array(
            {"height", term.powers.first, term.powers.second, value(term.coefficient)}));
    for (const auto cell : graph.cells)
        out.push_back(Json::array({"graph-cell", cell.first, cell.second}));
    for (const auto& x : take(result.source_points()))
        out.push_back(Json::array({"source-point", value(x.s), value(x.t)}));
    for (const auto& x : take(result.sources())) {
        Json desc = Json::array({x.descriptor.index()});
        std::visit(
            [&](const auto& source) {
                using S = std::decay_t<decltype(source)>;
                if constexpr (std::is_same_v<S, PgPartitionDomainSource>)
                    desc.push_back(source.edge);
                else if constexpr (std::is_same_v<S, PgPartitionThicknessSource>) {
                    desc.push_back(source.triangle);
                    desc.push_back(source.edge);
                } else if constexpr (std::is_same_v<S, PgPartitionNamedSource>) {
                    desc.push_back(source.boundary);
                    desc.push_back(source.segment);
                } else {
                    desc.push_back(source.ordinal);
                    desc.push_back(source.id);
                    desc.push_back(source.version);
                }
            },
            x.descriptor);
        out.push_back(Json::array({"source", desc, x.endpoints}));
    }
    const auto thickness = take(result.thickness());
    for (const auto& x : thickness.vertices)
        out.push_back(Json::array({"thickness-point", value(x.s), value(x.t)}));
    for (const auto& x : thickness.triangles)
        out.push_back(Json::array({"thickness-triangle", x}));
    for (const auto& x : thickness.values)
        out.push_back(Json::array({"thickness-value", value(x)}));
    for (const auto& x : take(result.layers()))
        out.push_back(Json::array({"layer", x.id.value(), value(x.lower), value(x.upper)}));
    for (const auto& x : take(result.boundaries()))
        out.push_back(
            Json::array({"boundary", x.text_offset, x.text_count, x.first_point, x.point_count}));
    const auto text = take(result.boundary_text());
    out.push_back(Json::array({"boundary-text", std::string(text.begin(), text.end())}));
    for (const auto b : take(result.identity()).bytes)
        out.push_back(Json::array({"identity", std::to_integer<unsigned>(b)}));
    for (const auto b : graph.identity.bytes)
        out.push_back(Json::array({"graph-identity", std::to_integer<unsigned>(b)}));
    return out.dump();
}
} // namespace partition_test
