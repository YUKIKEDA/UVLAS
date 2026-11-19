#pragma once
#include "aero/global_geometry.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace uvlas::test_global_geometry {
namespace s = aero::detail;
struct PairFixture {
    std::vector<s::PointId> points;
    std::vector<s::PointBinding> bindings;
    std::vector<double> high, low;
    std::vector<s::EdgeInput> edges;
    std::vector<s::EdgeBindingInput> edge_bindings;
    std::array<s::FaceInput, 2> faces{s::FaceInput{s::PhysicalFaceId{50}, s::BandId{100}},
                                      s::FaceInput{s::PhysicalFaceId{51}, s::BandId{101}}};
    std::array<s::BandId, 2> bands{s::BandId{100}, s::BandId{101}};
    std::vector<s::TriangleInput> triangles;
    std::array<s::FaceRoleInput, 2> roles{s::FaceRoleInput{faces[0].id, s::AeroFaceRole::Bound},
                                          s::FaceRoleInput{faces[1].id, s::AeroFaceRole::Bound}};
    PairFixture(const nlohmann::json& record, std::array<std::size_t, 3> pa = {0, 1, 2},
                std::array<std::size_t, 3> pb = {0, 1, 2}, int exponent = 0,
                bool common_origin = false) {
        using Position = std::array<double, 3>;
        const auto a = record.at("a").get<std::array<Position, 3>>(),
                   b = record.at("b").get<std::array<Position, 3>>();
        std::array<std::uint64_t, 3> ids_a{1, 2, 3}, ids_b{4, 5, 6};
        const auto claimed = [&](const Position& value) {
            for (const auto& p : record.at("allowed_points"))
                if (p.get<Position>() == value)
                    return true;
            if (!record.at("allowed_edge").is_null())
                for (const auto& p : record.at("allowed_edge"))
                    if (p.get<Position>() == value)
                        return true;
            return false;
        };
        for (std::size_t j = 0; j < 3; ++j)
            if (claimed(b[j]))
                for (std::size_t i = 0; i < 3; ++i)
                    if (a[i] == b[j])
                        ids_b[j] = ids_a[i];
        const auto insert_point = [&](std::uint64_t id, const Position& value) {
            if (std::find(points.begin(), points.end(), s::PointId{id}) != points.end())
                return;
            points.emplace_back(id);
            bindings.push_back(test_geometry::binding());
            for (const auto x : value) {
                const auto scaled = std::ldexp(x, exponent);
                high.push_back(common_origin ? 0x1p900 : scaled);
                low.push_back(common_origin ? scaled : 0);
            }
        };
        for (std::size_t i = 0; i < 3; ++i) {
            insert_point(ids_a[i], a[i]);
            insert_point(ids_b[i], b[i]);
        }
        std::array<std::uint64_t, 3> ea{10, 11, 12}, eb{20, 21, 22};
        for (std::size_t i = 0; i < 3; ++i)
            edges.emplace_back(s::EdgeId{ea[i]},
                               std::array{s::PointId{ids_a[i]}, s::PointId{ids_a[(i + 1) % 3]}});
        for (std::size_t j = 0; j < 3; ++j) {
            bool shared = false;
            if (!record.at("allowed_edge").is_null())
                for (std::size_t i = 0; i < 3; ++i) {
                    if ((ids_b[j] == ids_a[i] && ids_b[(j + 1) % 3] == ids_a[(i + 1) % 3]) ||
                        (ids_b[j] == ids_a[(i + 1) % 3] && ids_b[(j + 1) % 3] == ids_a[i])) {
                        eb[j] = ea[i];
                        shared = true;
                    }
                }
            if (!shared)
                edges.emplace_back(s::EdgeId{eb[j]}, std::array{s::PointId{ids_b[j]},
                                                                s::PointId{ids_b[(j + 1) % 3]}});
        }
        for (const auto& edge : edges)
            edge_bindings.push_back({edge.id, s::CanonicalEdge{}});
        const auto add_triangle = [&](std::array<std::uint64_t, 3> ids,
                                      std::array<std::uint64_t, 3> edge_ids,
                                      std::array<std::size_t, 3> order, std::size_t which) {
            std::array<s::PointId, 3> const vertex{
                s::PointId{ids[order[0]]}, s::PointId{ids[order[1]]}, s::PointId{ids[order[2]]}};
            std::array<s::EdgeId, 3> sides{s::EdgeId{0}, s::EdgeId{0}, s::EdgeId{0}};
            for (std::size_t k = 0; k < 3; ++k) {
                const auto x = order[k], y = order[(k + 1) % 3];
                for (std::size_t i = 0; i < 3; ++i)
                    if ((x == i && y == (i + 1) % 3) || (y == i && x == (i + 1) % 3))
                        sides[k] = s::EdgeId{edge_ids[i]};
            }
            triangles.emplace_back(s::FacetId{70 + which}, faces[which].id, vertex, sides);
        };
        add_triangle(ids_a, ea, pa, 0);
        add_triangle(ids_b, eb, pb, 1);
    }
    [[nodiscard]] s::GeometryInput input() const {
        return {points,
                bindings,
                s::DoublePositions{high, low},
                s::DoublePositions{high, low},
                edges,
                faces,
                bands,
                triangles,
                {},
                s::SingleWord{0.25},
                test_geometry::context()};
    }
    [[nodiscard]] s::GlobalTopologyInput topology() const {
        return {roles, edge_bindings, {}};
    }
};
} // namespace uvlas::test_global_geometry
