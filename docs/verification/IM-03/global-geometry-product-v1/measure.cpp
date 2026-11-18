// Fixed small-input evidence probe; not an R-03 performance qualification.
#include "aero/global_geometry.hpp"

#include <array>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <uvlas/core/instrumentation.hpp>

namespace {
using namespace uvlas;
namespace s = aero::detail;
constexpr std::uint64_t bytes = 64ULL * 1024 * 1024;
constexpr unsigned samples = 16;
constexpr auto operation = core::Operation::EvaluateAero;
template <class T>
T required(s::GlobalResult<T> result) {
    if (!result)
        throw std::runtime_error("global measure failed");
    return std::move(*result);
}
template <class T>
T required(core::Result<T> result) {
    if (!result)
        throw std::runtime_error("core measure failed");
    return std::move(*result);
}
template <class T>
T required(s::AeroResult<T> result) {
    if (!result)
        throw std::runtime_error("geometry measure failed");
    return std::move(*result);
}
} // namespace
int main(int argc, char** argv) {
    try {
        if (argc != 4)
            throw std::runtime_error(
                "CPU label, logical processors, and installed RAM bytes required");
        const core::MachineDescription machine{"windows-msvc-local", argv[1],
                                               static_cast<std::uint32_t>(std::stoul(argv[2])),
                                               std::stoull(argv[3])};
        core::ResourcePlan initial{bytes};
        auto old = required(s::GeometryWorkspace::create(initial));
        auto work = required(
            s::GlobalWorkspace::create(initial, numerics::linear::EstimateAdmission{4096}));
        const std::array ids{s::PointId{1}, s::PointId{2}, s::PointId{3}, s::PointId{4},
                             s::PointId{5}};
        const s::PointBinding canonical{s::ReferenceMapId{1}, s::MappingGeneration{1},
                                        s::CanonicalPoint{}};
        const std::array binding{
            canonical, canonical, canonical, canonical,
            s::PointBinding{s::ReferenceMapId{1}, s::MappingGeneration{1}, s::AliasOf{ids[1]}}};
        const std::array<double, 15> positions{0, 0, 0, 2, 0, 0, 0, 2, 0, 0, -2, 1, 2, 0, 0};
        const std::array edges{s::EdgeInput{s::EdgeId{10}, {ids[0], ids[1]}},
                               s::EdgeInput{s::EdgeId{11}, {ids[1], ids[2]}},
                               s::EdgeInput{s::EdgeId{12}, {ids[2], ids[0]}},
                               s::EdgeInput{s::EdgeId{13}, {ids[0], ids[3]}},
                               s::EdgeInput{s::EdgeId{14}, {ids[3], ids[1]}}};
        const std::array faces{s::FaceInput{s::PhysicalFaceId{50}, s::BandId{100}},
                               s::FaceInput{s::PhysicalFaceId{51}, s::BandId{101}}};
        const std::array bands{s::BandId{100}, s::BandId{101}};
        const std::array triangles{s::TriangleInput{s::FacetId{70},
                                                    faces[0].id,
                                                    {ids[0], ids[1], ids[2]},
                                                    {edges[0].id, edges[1].id, edges[2].id}},
                                   s::TriangleInput{s::FacetId{71},
                                                    faces[1].id,
                                                    {ids[1], ids[0], ids[3]},
                                                    {edges[0].id, edges[3].id, edges[4].id}}};
        const s::GeometryContext context{core::ModelFingerprint{},
                                         core::ModelGeneration{1},
                                         s::ReferenceGeometryGeneration{1},
                                         core::GeometryGeneration{1},
                                         s::MappingGeneration{1},
                                         s::ConnectivityGeneration{1},
                                         s::KernelPolicyId{1},
                                         s::KernelGeneration{1},
                                         s::RadiusGeneration{1},
                                         core::PrecisionGeneration{1},
                                         core::NumericalGeneration{1},
                                         0};
        auto source = required(s::prepare_owned_geometry(
            {ids,
             binding,
             s::SinglePositions{positions},
             s::SinglePositions{positions},
             edges,
             faces,
             bands,
             triangles,
             {},
             s::SingleWord{0.25},
             context},
            {{10, 10, 10, 10, 10, 10}, 100000, {numerics::linear::EstimateAdmission{28672}, 100}},
            *old, initial));
        const s::GlobalLimits limits{UINT64_MAX / 4, numerics::linear::EstimateAdmission{4096}};
        const std::array roles{s::FaceRoleInput{faces[0].id, s::AeroFaceRole::Bound},
                               s::FaceRoleInput{faces[1].id, s::AeroFaceRole::Bound}};
        const std::array aliases{s::EdgeBindingInput{edges[0].id, s::CanonicalEdge{}},
                                 s::EdgeBindingInput{edges[1].id, s::CanonicalEdge{}},
                                 s::EdgeBindingInput{edges[2].id, s::CanonicalEdge{}},
                                 s::EdgeBindingInput{edges[3].id, s::CanonicalEdge{}},
                                 s::EdgeBindingInput{edges[4].id, s::CanonicalEdge{}}};
        const std::array te{s::TrailingEdgeInput{edges[1].id, {ids[1], ids[2]}},
                            s::TrailingEdgeInput{edges[4].id, {ids[1], ids[3]}}};
        const s::GlobalTopologyInput topology_input{roles, aliases, te};
        const std::array<s::ScalarWords, 3> axis{s::SingleWord{1}, s::SingleWord{0},
                                                 s::SingleWord{0}},
            zero{s::SingleWord{0}, s::SingleWord{0}, s::SingleWord{0}};
        const std::array dofs{s::GlobalDofInput{core::DofId{1},
                                                {1, 0, 0, 0, 0, 0, 0},
                                                s::DofBasisId{1},
                                                axis,
                                                s::DoubleWord{1, 0x1p-200}}};
        const std::array coefficients{s::GlobalJacobianEntry{core::DofId{1}, s::SingleWord{1}}};
        const std::array points{s::PointKinematicsInput{ids[0], {coefficients, {}, {}}, zero},
                                s::PointKinematicsInput{ids[1], {coefficients, {}, {}}, zero},
                                s::PointKinematicsInput{ids[2], {coefficients, {}, {}}, zero},
                                s::PointKinematicsInput{ids[3], {coefficients, {}, {}}, zero},
                                s::PointKinematicsInput{ids[4], {coefficients, {}, {}}, zero}};
        const s::AliasKinematicsInput kinematic_input{s::DofLayoutGeneration{1}, dofs, points};
        auto topology =
            required(s::build_global_topology(source, topology_input, limits, *work, initial));
        core::StageRecorder topology_time{true, operation}, contact_time{true, operation},
            kinematics_time{true, operation};
        core::ResourcePlan topology_memory{bytes}, contact_memory{bytes}, kinematics_memory{bytes};
        std::array<s::GlobalWorkPlan, 3> plans{};
        for (unsigned i = 0; i < samples; ++i) {
            std::optional<s::GlobalTopology> result;
            {
                core::StageScope timer{topology_time, core::MeasurementStage::Preparation};
                result.emplace(required(s::build_global_topology(source, topology_input, limits,
                                                                 *work, topology_memory)));
            }
            plans[0] = required(result->work_plan());
        }
        for (unsigned i = 0; i < samples; ++i) {
            std::optional<s::GlobalContactState> result;
            {
                core::StageScope timer{contact_time, core::MeasurementStage::Aerodynamics};
                result.emplace(required(
                    s::check_global_geometry(topology, source, limits, *work, contact_memory)));
            }
            plans[1] = required(result->work_plan());
        }
        for (unsigned i = 0; i < samples; ++i) {
            std::optional<s::AliasKinematicsState> result;
            {
                core::StageScope timer{kinematics_time, core::MeasurementStage::Preparation};
                result.emplace(required(s::check_alias_kinematics(source, kinematic_input, limits,
                                                                  *work, kinematics_memory)));
            }
            plans[2] = required(result->work_plan());
        }
        const auto memory = required(core::sample_process_memory(operation));
        const std::string_view settings =
            "{\"samples\":16,\"facets\":2,\"points\":5,\"physical_owners\":4,\"TE_edges\":2,"
            "\"DOFs\":1,\"CSR_terms\":5,\"word_planes\":1,\"workers\":1}";
        const auto report = [&](std::string_view name, const auto& recorder,
                                const auto& resources) {
            return required(core::measurement_report_json({machine, name, 1, settings}, recorder,
                                                          resources, memory, operation));
        };
        std::cout << "{\"topology\":"
                  << report("global-topology-small", topology_time, topology_memory)
                  << ",\"contact\":" << report("global-contact-small", contact_time, contact_memory)
                  << ",\"kinematics\":"
                  << report("global-alias-small", kinematics_time, kinematics_memory)
                  << ",\"work_plans\":[";
        for (std::size_t i = 0; i < plans.size(); ++i) {
            if (i != 0)
                std::cout << ',';
            const auto& p = plans[i];
            std::cout << "{\"payload\":" << p.payload_bytes << ",\"reserved\":" << p.reserved_steps
                      << ",\"executed\":" << p.executed_steps << ",\"candidates\":" << p.candidates
                      << ",\"metadata_forecast\":" << p.metadata_forecast
                      << ",\"traversal_reserved\":" << p.traversal_reserved_steps
                      << ",\"traversal_executed\":" << p.traversal_executed_steps << '}';
        }
        std::cout << "]}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
