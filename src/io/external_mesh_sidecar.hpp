#pragma once

#include "gmsh_mesh.hpp"
#include "nastran_isotropic_shell_cards.hpp"
#include "nastran_support_cards.hpp"

#include <type_traits>
#include <uvlas/io/si_case.hpp>
#include <uvlas/structure/pg_integration_partition.hpp>

namespace uvlas::io::detail {

template <SchemaNamespace Space>
using MeshCaseId = core::Identifier<std::integral_constant<SchemaNamespace, Space>>;
using MeshGeometryId = MeshCaseId<SchemaNamespace::Geometry>;
using MeshShellId = MeshCaseId<SchemaNamespace::Shell>;
using MeshMaterialId = MeshCaseId<SchemaNamespace::Material>;
using MeshSupportId = MeshCaseId<SchemaNamespace::Support>;
using MeshJointId = MeshCaseId<SchemaNamespace::Joint>;

enum class MeshSidecarOrientation : std::uint8_t { Same, Reversed };

struct MeshSidecarChart {
    std::array<structure::PgPartitionPoint, 4> vertices;
    std::uint8_t count{};
    MeshSidecarOrientation orientation = MeshSidecarOrientation::Same;
};

struct MeshConversionObservation {
    core::ObservationId id{0};
    core::Index node{};
};

struct MeshConversionBudget {
    structure::ExactRational distance;
    structure::ExactRational normal;
    core::PlannedBuffer<MeshConversionObservation> observations;
};

struct GmshSidecarFace {
    std::uint64_t element{};
    MeshSidecarChart chart;
    core::Index node{};
};

struct GmshSidecarLayer {
    core::LayerId layer{0};
    MeshMaterialId material{0};
    core::Index node{};
};

struct GmshSidecarShell {
    GmshKey physical;
    MeshShellId shell{0};
    core::PlannedBuffer<GmshSidecarFace> faces;
    std::optional<core::PlannedBuffer<GmshSidecarLayer>> material_layers;
    std::optional<structure::ExactRational> full_thickness;
    std::optional<core::Index> full_thickness_node;
    core::Index node{};
};

struct GmshSidecarSupport {
    GmshKey physical;
    MeshSupportId support{0};
    core::Index node{};
};

struct GmshSidecarJoint {
    GmshKey first, second;
    MeshJointId joint{0};
    core::Index node{};
};

struct GmshSidecar {
    SiLengthUnit length = SiLengthUnit::Metre;
    core::PlannedBuffer<GmshSidecarShell> shells;
    core::PlannedBuffer<GmshSidecarSupport> supports;
    core::PlannedBuffer<GmshSidecarJoint> joints;
};

struct NastranSidecarShell {
    NastranElementId element{0};
    MeshShellId shell{0};
    MeshSidecarChart chart;
    core::Index node{};
};

template <SchemaNamespace Space>
struct NastranSidecarElementBinding {
    NastranElementId element{0};
    MeshCaseId<Space> target{0};
    core::Index node{};
};

struct NastranSidecarMaterial {
    NastranMaterialId mid{0};
    MeshMaterialId material{0};
    core::Index node{};
};

struct NastranSidecarSupport {
    NastranSupportSetId set{0};
    NastranGridId grid{0};
    std::uint8_t component_bits{};
    MeshSupportId support{0};
    core::Index node{};
};

struct NastranSidecar {
    SiUnits units{SiLengthUnit::Metre, SiMassUnit::Kilogram, SiTimeUnit::Second,
                  SiAngleUnit::Radian};
    core::PlannedBuffer<char> generator;
    core::PlannedBuffer<NastranSidecarShell> shells;
    core::PlannedBuffer<NastranSidecarMaterial> materials;
    core::PlannedBuffer<NastranSidecarElementBinding<SchemaNamespace::Beam>> beams;
    core::PlannedBuffer<NastranSidecarElementBinding<SchemaNamespace::Mass>> masses;
    core::PlannedBuffer<NastranSidecarSupport> supports;
    core::PlannedBuffer<NastranSidecarElementBinding<SchemaNamespace::Joint>> joints;
};

enum class MeshSidecarFailure : std::uint8_t {
    Snapshot,
    Json,
    Schema,
    Shape,
    Key,
    Scalar,
    Chart,
    Resource,
    Binding,
    CaseReference,
    Selection
};

struct MeshSidecarDiagnostic {
    MeshSidecarFailure reason = MeshSidecarFailure::Snapshot;
    std::optional<RootSnapshotLocation> source;
    std::optional<core::Index> node;
    std::optional<std::uint64_t> byte;
    std::optional<GmshKey> physical;
    std::optional<std::uint64_t> element;
    std::optional<SchemaNamespace> target_namespace;
    std::optional<std::uint64_t> target;
};

struct MeshSidecarLimits {
    MeshSidecarLimits(InputLimits document, structure::ExactLimits algebra, std::uint64_t rows,
                      std::uint64_t faces, std::uint64_t work) noexcept
        : json(document),
          exact(algebra),
          maximum_rows(rows),
          maximum_faces(faces),
          maximum_work(work) {}

    InputLimits json;
    structure::ExactLimits exact;
    std::uint64_t maximum_rows;
    std::uint64_t maximum_faces;
    std::uint64_t maximum_work;
};

using MeshSidecarPayload = std::variant<GmshSidecar, NastranSidecar>;

/// Owns typed declarations and source identity independently of the snapshot.
class ParsedMeshSidecar {
public:
    ParsedMeshSidecar(ParsedMeshSidecar&&) noexcept = default;
    ParsedMeshSidecar& operator=(ParsedMeshSidecar&&) = delete;

    [[nodiscard]] RootSnapshotLocation source() const noexcept {
        return source_;
    }
    [[nodiscard]] SourceDigest digest() const noexcept {
        return digest_;
    }
    [[nodiscard]] MeshGeometryId geometry() const noexcept {
        return geometry_;
    }
    [[nodiscard]] std::uint64_t work() const noexcept {
        return work_;
    }
    [[nodiscard]] const MeshSidecarPayload& payload() const& noexcept {
        return payload_;
    }
    [[nodiscard]] const MeshSidecarPayload& payload() const&& = delete;
    [[nodiscard]] const std::optional<MeshConversionBudget>& conversion() const& noexcept {
        return conversion_;
    }
    [[nodiscard]] const std::optional<MeshConversionBudget>& conversion() const&& = delete;

private:
    friend core::Result<ParsedMeshSidecar>
    read_mesh_sidecar(const AcquiredRootSnapshot&, const MeshSidecarLimits&, core::ResourcePlan&,
                      const std::stop_token&, std::optional<MeshSidecarDiagnostic>*) noexcept;

    ParsedMeshSidecar(RootSnapshotLocation source, SourceDigest digest, MeshGeometryId geometry,
                      MeshSidecarPayload&& payload,
                      std::optional<MeshConversionBudget>&& conversion, std::uint64_t work) noexcept
        : source_(source),
          digest_(digest),
          geometry_(geometry),
          payload_(std::move(payload)),
          conversion_(std::move(conversion)),
          work_(work) {}
    RootSnapshotLocation source_;
    SourceDigest digest_;
    MeshGeometryId geometry_;
    MeshSidecarPayload payload_;
    std::optional<MeshConversionBudget> conversion_;
    std::uint64_t work_;
};

/// Read the adopted closed sidecar schema from the already authenticated bytes.
[[nodiscard]] core::Result<ParsedMeshSidecar>
read_mesh_sidecar(const AcquiredRootSnapshot&, const MeshSidecarLimits&, core::ResourcePlan&,
                  const std::stop_token&, std::optional<MeshSidecarDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
