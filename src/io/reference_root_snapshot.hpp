#pragma once

#include "reference_dependency_walk.hpp"
#include "reference_input.hpp"
#include "reference_root_plan.hpp"

#include <functional>

namespace uvlas::io::detail {

struct RootSnapshotLimits {
    ReferenceFileLimits file;
    std::uint64_t maximum_identity_comparisons;
    std::uint64_t diagnostic_bytes;
};

struct RootSnapshotLocation {
    SchemaIdentity parent;
    SourceId source;
    ReferenceRole role;
    std::uint64_t declaration_node;
};

enum class RootSnapshotFailure : std::uint8_t {
    InvalidSelection,
    StaleParent,
    Open,
    Read,
    ProvenanceStorage,
    SourceLimit,
    InventoryStorage,
};

struct RootSnapshotDiagnostic {
    std::optional<RootSnapshotFailure> reason;
    std::optional<RootSnapshotLocation> location;
    std::optional<ManifestOpenDiagnostic> open;
    std::optional<ReferenceFileDiagnostic> read;
    std::optional<ReferenceFileQuota> source_quota;
    core::PlannedBuffer<char> declaration_pointer;

    [[nodiscard]] std::string_view pointer() const& noexcept;
    [[nodiscard]] std::string_view pointer() const&& = delete;
    void reset() noexcept;
};

struct AcquiredRootSnapshot {
    RootSnapshotLocation location;
    ReferenceFileSnapshot file;
    core::PlannedBuffer<char> declaration_pointer;

    [[nodiscard]] std::string_view pointer() const& noexcept;
    [[nodiscard]] std::string_view pointer() const&& = delete;
};

// Acquires one previously unseen source ordinal from the same checked case and bound root.
// Prior IDs must be from already accepted distinct sources in that same origin. This stage
// consumes neither input owner and does not make any claim about a source's dependencies.
// Only the snapshot and its declaration provenance are retained when the handle closes.
[[nodiscard]] core::Result<AcquiredRootSnapshot>
acquire_root_snapshot(const ReferenceInput& input, const PreparedRootReferences& prepared,
                      core::Index source_ordinal, std::span<const ReferenceFileId> earlier_sources,
                      RootSnapshotLimits limits, core::ResourcePlan& plan,
                      const std::stop_token& stop, RootSnapshotDiagnostic* diagnostic) noexcept;

struct RootInventoryLimits {
    RootSnapshotLimits per_source;
    std::uint64_t maximum_sources;
    std::uint64_t total_payload_bytes;
    std::uint64_t total_read_bytes;
    std::uint64_t total_read_calls;
    std::uint64_t total_hash_bytes;
    std::uint64_t total_identity_comparisons;
};

// Owns all authenticated files declared directly by one manifest and their declaration records.
// Views borrow this inventory until it moves or is destroyed. Dependency extraction and
// role-specific qualification are separate operations.
class AcquiredRootReferences {
public:
    AcquiredRootReferences(AcquiredRootReferences&&) noexcept = default;
    AcquiredRootReferences& operator=(AcquiredRootReferences&&) = delete;

    [[nodiscard]] core::Result<SchemaIdentity> identity() const noexcept;
    [[nodiscard]] core::Result<RootReferenceReport> report() const noexcept;
    [[nodiscard]] core::Result<RootReferenceDeclarationView>
        declaration(core::Index) const& noexcept;
    [[nodiscard]] core::Result<RootReferenceDeclarationView>
        declaration(core::Index) const&& = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
    find(const SchemaIdentity& parent, SourceId source, ReferenceRole role,
         const std::stop_token& stop) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
    find(const SchemaIdentity&, SourceId, ReferenceRole, const std::stop_token&) const&& = delete;

private:
    friend core::Result<AcquiredRootReferences>
    acquire_root_references(const ReferenceInput&, PreparedRootReferences&&, RootInventoryLimits,
                            core::ResourcePlan&, const std::stop_token&,
                            RootSnapshotDiagnostic*) noexcept;
    AcquiredRootReferences(
        PreparedRootReferences declarations,
        core::PlannedBuffer<std::optional<AcquiredRootSnapshot>> sources) noexcept;

    PreparedRootReferences declarations_;
    core::PlannedBuffer<std::optional<AcquiredRootSnapshot>> sources_;
};

// Authenticates each distinct directly declared source once, in first declaration order.
// Prepared metadata is moved only on success. All successful charges stay on the shared Plan.
[[nodiscard]] core::Result<AcquiredRootReferences>
acquire_root_references(const ReferenceInput& input, PreparedRootReferences&& prepared,
                        RootInventoryLimits limits, core::ResourcePlan& plan,
                        const std::stop_token& stop, RootSnapshotDiagnostic* diagnostic) noexcept;

enum class ReferenceWalkStartFailure : std::uint8_t {
    Binding,
    Sources,
    WalkStorage,
    Snapshot,
    Edge
};

struct ReferenceWalkStartDiagnostic {
    std::optional<ReferenceWalkStartFailure> reason;
    std::optional<ReferenceFileQuota> sources;
    RootSnapshotDiagnostic snapshot;
    DependencyWalkDiagnostic walk;
    core::PlannedBuffer<char> edge_pointer;

    [[nodiscard]] std::string_view pointer() const& noexcept;
    [[nodiscard]] std::string_view pointer() const&& = delete;
    void reset() noexcept;
};

struct PendingRootDependencies {
    DependencyEdge edge;
    AcquiredRootSnapshot snapshot;
};

struct ReferenceWalkStartView {
    SchemaIdentity parent{{SourceId{0}, 0, {}}, 0, 0};
    DependencyWalkView walk;
    // Present only after the first reference has been authenticated and entered.
    const PendingRootDependencies* pending_root = nullptr;
};

// Owns an intermediate walk before role-specific dependency extraction. Views borrow this
// owner until it moves or is destroyed. No operation can complete an unchecked source.
class ReferenceWalkStart {
public:
    ReferenceWalkStart(ReferenceWalkStart&&) noexcept = default;
    ReferenceWalkStart& operator=(ReferenceWalkStart&&) = delete;

    [[nodiscard]] std::optional<ReferenceWalkStartView> view() const& noexcept;
    [[nodiscard]] std::optional<ReferenceWalkStartView> view() const&& = delete;

private:
    friend core::Result<ReferenceWalkStart>
    start_reference_walk(const ReferenceInput&, const PreparedRootReferences&, DependencyWalkLimits,
                         RootSnapshotLimits, core::ResourcePlan&, const std::stop_token&,
                         ReferenceWalkStartDiagnostic*) noexcept;

    ReferenceWalkStart(SchemaIdentity parent, DependencyWalk walk,
                       std::optional<PendingRootDependencies> pending) noexcept;
    SchemaIdentity parent_;
    DependencyWalk walk_;
    std::optional<PendingRootDependencies> pending_;
};

// Preflights storage for known root ordinals; discovering new ordinals is a later stage. Only the
// first declaration is opened, verified, and entered. Later declarations are left for the
// role-specific continuation.
[[nodiscard]] core::Result<ReferenceWalkStart>
start_reference_walk(const ReferenceInput& input, const PreparedRootReferences& prepared,
                     DependencyWalkLimits walk_limits, RootSnapshotLimits snapshot_limits,
                     core::ResourcePlan& plan, const std::stop_token& stop,
                     ReferenceWalkStartDiagnostic* diagnostic) noexcept;

} // namespace uvlas::io::detail
