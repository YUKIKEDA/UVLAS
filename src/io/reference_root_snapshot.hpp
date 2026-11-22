#pragma once

#include "reference_input.hpp"
#include "reference_root_plan.hpp"

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
};

struct RootSnapshotDiagnostic {
    std::optional<RootSnapshotFailure> reason;
    std::optional<RootSnapshotLocation> location;
    std::optional<ManifestOpenDiagnostic> open;
    std::optional<ReferenceFileDiagnostic> read;
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

} // namespace uvlas::io::detail
