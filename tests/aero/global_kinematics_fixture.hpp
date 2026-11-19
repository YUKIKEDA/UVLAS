#pragma once
#include "aero/global_geometry.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <vector>

namespace uvlas::test_global_geometry {
namespace s = aero::detail;
struct KinematicsFixture {
    test_geometry::TriangleFixture base;
    std::vector<s::PointId> ids{base.points.begin(), base.points.end()};
    std::vector<s::PointBinding> bindings{base.bindings.begin(), base.bindings.end()};
    std::vector<double> positions{base.reference.begin(), base.reference.end()};
    const std::array<s::ScalarWords, 3> axis{s::SingleWord{1}, s::SingleWord{0}, s::SingleWord{0}};
    std::array<s::GlobalDofInput, 3> dofs{
        s::GlobalDofInput{
            core::DofId{1}, {1, 0, 0, 0, 0, 0, 0}, s::DofBasisId{1}, axis, s::SingleWord{1}},
        s::GlobalDofInput{
            core::DofId{2}, {1, 0, 0, 0, 0, 0, 0}, s::DofBasisId{1}, axis, s::SingleWord{0x1p-200}},
        s::GlobalDofInput{core::DofId{3},
                          {1, 0, 0, 0, 0, 0, 0},
                          s::DofBasisId{1},
                          axis,
                          s::SingleWord{0x1p-400}}};
    std::array<s::GlobalJacobianEntry, 3> coefficients{
        s::GlobalJacobianEntry{core::DofId{1}, s::SingleWord{1}},
        s::GlobalJacobianEntry{core::DofId{2}, s::SingleWord{1}},
        s::GlobalJacobianEntry{core::DofId{3}, s::SingleWord{1}}};
    std::array<s::GlobalJacobianEntry, 3> alias_coefficients = coefficients;
    const std::array<s::ScalarWords, 3> zero{s::SingleWord{0}, s::SingleWord{0}, s::SingleWord{0}};
    std::vector<s::PointKinematicsInput> points;
    KinematicsFixture() {
        ids.emplace_back(0);
        bindings.emplace_back(s::ReferenceMapId{1}, s::MappingGeneration{1},
                              s::AliasOf{base.points[1]});
        positions.insert(positions.end(), base.reference.begin() + 3, base.reference.begin() + 6);
        for (const auto id : ids)
            points.push_back({id, {coefficients, {}, {}}, zero, {}});
        points.back().jacobian[0] = alias_coefficients;
    }
    [[nodiscard]] s::GeometryInput geometry() const {
        return {ids,
                bindings,
                s::SinglePositions{positions},
                s::SinglePositions{positions},
                base.edges,
                base.faces,
                base.bands,
                base.triangles,
                {},
                s::SingleWord{.25},
                test_geometry::context()};
    }
    [[nodiscard]] s::AliasKinematicsInput input() const {
        return {s::DofLayoutGeneration{19}, dofs, points};
    }
};
} // namespace uvlas::test_global_geometry
