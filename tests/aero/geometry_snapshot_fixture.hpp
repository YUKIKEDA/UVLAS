#pragma once

#include "aero/geometry_snapshot.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <stdexcept>
#include <utility>

namespace uvlas::test_geometry {
namespace subject = uvlas::aero::detail;
constexpr std::uint64_t budget = 16ULL * 1024 * 1024;
inline subject::GeometryContext context() {
    core::FailureDetails details;
    details.parent = core::AcceptedStateId{2};
    details.trial = core::TrialId{3};
    details.last_accepted = core::AcceptedStateId{2};
    return {core::ModelFingerprint{},
            core::ModelGeneration{1},
            subject::ReferenceGeometryGeneration{2},
            core::GeometryGeneration{3},
            subject::MappingGeneration{4},
            subject::ConnectivityGeneration{5},
            subject::KernelPolicyId{6},
            subject::KernelGeneration{7},
            subject::RadiusGeneration{8},
            core::PrecisionGeneration{9},
            core::NumericalGeneration{10},
            1.25,
            details};
}
inline subject::GeometryLimits limits() {
    return {
        {100, 100, 100, 100, 100, 100}, 10000, {numerics::linear::EstimateAdmission{28672}, 100}};
}
template <class T>
T take(subject::AeroResult<T> result) {
    if (!result) {
        FAIL_CHECK("Unexpected aero failure " << static_cast<unsigned>(result.error().cause.code));
        throw std::runtime_error("Unexpected aero failure");
    }
    return std::move(*result);
}
inline std::unique_ptr<subject::GeometryWorkspace> workspace(core::ResourcePlan& plan) {
    auto result = subject::GeometryWorkspace::create(plan);
    if (!result) {
        FAIL_CHECK("Workspace creation failed");
        throw std::runtime_error("Workspace creation failed");
    }
    return std::move(*result);
}
inline subject::PointBinding binding() {
    return {subject::ReferenceMapId{1}, subject::MappingGeneration{1}, subject::CanonicalPoint{}};
}
struct TriangleFixture {
    std::array<subject::PointId, 3> points{subject::PointId{11}, subject::PointId{3},
                                           subject::PointId{7}};
    std::array<subject::PointBinding, 3> bindings{binding(), binding(), binding()};
    std::array<double, 9> reference{0, 0, 0, 3, 0, 0, 0, 4, 0}, current = reference;
    std::array<double, 9> reference_low{}, current_low{};
    std::array<subject::EdgeInput, 3> edges{
        subject::EdgeInput{subject::EdgeId{10}, {points[0], points[1]}},
        subject::EdgeInput{subject::EdgeId{20}, {points[1], points[2]}},
        subject::EdgeInput{subject::EdgeId{30}, {points[2], points[0]}}};
    std::array<subject::FaceInput, 1> faces{
        subject::FaceInput{subject::PhysicalFaceId{50}, subject::BandId{100}}};
    std::array<subject::BandId, 1> bands{subject::BandId{100}};
    std::array<subject::TriangleInput, 1> triangles{subject::TriangleInput{
        subject::FacetId{70}, faces[0].id, points, {edges[0].id, edges[1].id, edges[2].id}}};
    subject::ScalarWords radius = subject::SingleWord{.25};
    subject::GeometryContext stamp = context();
    bool two_word = false;
    [[nodiscard]] subject::GeometryInput input() const {
        const subject::PositionView r =
            two_word ? subject::PositionView{subject::DoublePositions{reference, reference_low}}
                     : subject::PositionView{subject::SinglePositions{reference}};
        const subject::PositionView c =
            two_word ? subject::PositionView{subject::DoublePositions{current, current_low}}
                     : subject::PositionView{subject::SinglePositions{current}};
        return {points, bindings, r, c, edges, faces, bands, triangles, {}, radius, stamp};
    }
};
template <class T>
void failure(const subject::AeroResult<T>& result, core::FailureCode code) {
    if (result) {
        FAIL_CHECK("Unexpected aero success");
        return;
    }
    CHECK(result.error().cause.code == code);
    CHECK(result.error().cause.operation == core::Operation::EvaluateAero);
}
} // namespace uvlas::test_geometry
