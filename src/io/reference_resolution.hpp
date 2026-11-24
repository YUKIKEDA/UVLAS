#pragma once

#include "external_mesh_sidecar.hpp"
#include "gmsh_mesh.hpp"
#include "nastran_bulk.hpp"

namespace uvlas::io::detail {

enum class LeafReferenceFailure : std::uint8_t { UnsupportedRole, Format, Graph, Sidecar };

struct LeafFormatLimits {
    GmshLimits gmsh{};
    NastranBulkLimits nastran{};
    std::optional<MeshSidecarLimits> sidecar{};
};

struct LeafReferenceDiagnostic {
    std::optional<LeafReferenceFailure> reason;
    std::optional<RootSnapshotLocation> location;
    std::optional<std::size_t> format_offset;
    std::optional<SourceId> declaring_source;
    std::optional<ReferencePathDiagnostic> path;
    DependencyWalkDiagnostic walk;
    std::optional<MeshSidecarDiagnostic> sidecar;
    core::PlannedBuffer<char> declaration_pointer;

    [[nodiscard]] std::string_view pointer() const& noexcept;
    [[nodiscard]] std::string_view pointer() const&& = delete;
    void reset() noexcept;
};

// A format supplies these from the authenticated source, in declaration order. The byte
// offset is in that source. The root-relative spelling follows any role-specific path rule.
struct ChildReferenceDeclaration {
    SourceId source{0};
    ReferenceRole role{};
    SourceDigest digest{};
    std::uint64_t declaration_byte{};
    std::string_view original_path;
    std::string_view root_relative_path;
};

// Declaration text views point into the authenticated bytes or the returned text buffer.
struct ExtractedReferenceDependencies {
    core::PlannedBuffer<ChildReferenceDeclaration> declarations;
    core::PlannedBuffer<char> text;
};

class ReferenceDependencyExtractor {
public:
    virtual ~ReferenceDependencyExtractor() = default;
    [[nodiscard]] virtual core::Result<ExtractedReferenceDependencies>
    extract(const AcquiredRootSnapshot&) noexcept = 0;
};

struct ReferenceGraphLimits {
    RootInventoryLimits inventory;
    DependencyWalkLimits walk;
    // Root limits are checked in PreparedRootReferences. These bounds also cover new children.
    std::uint64_t normalized_path_bytes;
    ReferencePathLimits child_path;
};

struct ReferenceGraphEdgeView {
    DependencyEdge topology;
    std::string_view original_path;
    std::string_view declaration_pointer; // Present only for a manifest occurrence.
};

class ReferenceGraphBuilder;

class ReferencedCase {
public:
    ReferencedCase(ReferencedCase&&) noexcept = default;
    ReferencedCase& operator=(ReferencedCase&&) = delete;

    [[nodiscard]] core::Result<SchemaIdentity> identity() const noexcept;
    [[nodiscard]] std::span<const core::Index> dependency_first() const& noexcept;
    [[nodiscard]] std::span<const core::Index> dependency_first() const&& = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
    snapshot(core::Index source) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
        snapshot(core::Index) const&& = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
    find(const SchemaIdentity&, SourceId, ReferenceRole, const std::stop_token&) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
    find(const SchemaIdentity&, SourceId, ReferenceRole, const std::stop_token&) const&& = delete;
    [[nodiscard]] core::Result<std::string_view>
    normalized_path(core::Index source) const& noexcept;
    [[nodiscard]] core::Result<std::string_view> normalized_path(core::Index) const&& = delete;
    // Consumer slot zero is the manifest, source ordinal n uses slot n+1.
    [[nodiscard]] core::Result<core::Index>
        dependency_count(std::optional<core::Index>) const noexcept;
    [[nodiscard]] core::Result<ReferenceGraphEdgeView>
    dependency(std::optional<core::Index> consumer, core::Index ordinal) const& noexcept;
    [[nodiscard]] core::Result<ReferenceGraphEdgeView> dependency(std::optional<core::Index>,
                                                                  core::Index) const&& = delete;

private:
    struct Vertex {
        SourceId id{0};
        ReferenceRole role{};
        SourceDigest digest{};
        std::optional<std::uint64_t> declaring_mesh_geometry;
        std::optional<NormalizedReferencePath> child_path;
        core::PlannedBuffer<char> original_child_path;
        std::optional<AcquiredRootSnapshot> snapshot;
    };

    struct Occurrence {
        DependencyEdge edge{{}, 0, SourceId{0}, {}, 0};
        core::Index root_declaration = -1;
        core::PlannedBuffer<char> original_child_path;
    };

    friend class ReferenceGraphBuilder;
    ReferencedCase(PreparedRootReferences&&, core::PlannedBuffer<Vertex>&&, core::Index,
                   core::Index, core::PlannedBuffer<Occurrence>&&,
                   core::PlannedBuffer<core::Index>&&, core::PlannedBuffer<core::Index>&&) noexcept;
    [[nodiscard]] core::Result<std::span<const Occurrence>>
        occurrences(std::optional<core::Index>) const noexcept;

    PreparedRootReferences prepared_;
    core::PlannedBuffer<Vertex> vertices_;
    core::Index vertex_count_{};
    core::Index root_count_{};
    core::PlannedBuffer<Occurrence> occurrences_;
    core::PlannedBuffer<core::Index> offsets_;
    core::PlannedBuffer<core::Index> dependency_first_;
};

// Internal collection engine. Production supplies the role dispatch below; a role must actually
// inspect the authenticated bytes before it can assert that its child list is empty.
[[nodiscard]] core::Result<ReferencedCase>
resolve_reference_graph(const ReferenceInput&, PreparedRootReferences&&, ReferenceGraphLimits,
                        ReferenceDependencyExtractor&, core::ResourcePlan&, const std::stop_token&,
                        RootSnapshotDiagnostic*, LeafReferenceDiagnostic*) noexcept;

[[nodiscard]] core::Result<ReferencedCase>
resolve_leaf_references(const ReferenceInput&, PreparedRootReferences&&, RootInventoryLimits,
                        DependencyWalkLimits, LeafFormatLimits, core::ResourcePlan&,
                        const std::stop_token&, RootSnapshotDiagnostic*,
                        LeafReferenceDiagnostic*) noexcept;

// Resolve the authenticated deck from the same case and retain its immutable bytes for
// typed clients. The input graph can be destroyed after a successful call.
[[nodiscard]] core::Result<NastranBulkDeck>
prepare_referenced_nastran_deck(const ReferencedCase&, const SchemaIdentity&, SourceId,
                                NastranBulkLimits, core::ResourcePlan&, const std::stop_token&,
                                std::optional<NastranBulkDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
