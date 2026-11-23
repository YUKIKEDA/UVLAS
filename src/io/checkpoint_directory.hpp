#pragma once

#include "checkpoint_manifest.hpp"
#include "checkpoint_publication.hpp"
#include "reference_file_snapshot.hpp"

#include <filesystem>

namespace uvlas::io::detail {

struct LoadedCheckpointDirectory {
    std::uint64_t number;
    CheckpointManifest manifest;
    CheckpointPayloadSet payloads;
};

struct CheckpointDiscoveryLimits {
    std::uint64_t maximum_entries;
    std::uint64_t maximum_candidates;
};

struct CheckpointDirectoryLimits {
    ReferenceFileLimits file;
    std::uint64_t maximum_entries;
};

/// Build the manifest directly from the typed chunks, then publish the native
/// transaction under a new number. The same checked bytes feed both operations.
[[nodiscard]] core::Result<void>
save_checkpoint_directory(const std::filesystem::path& parent, std::uint64_t number,
                          std::span<const NamedCheckpointChunk>, const CheckpointFingerprints&,
                          const CheckpointManifestLimits&, CheckpointPublicationLimits,
                          core::ResourcePlan&, const std::stop_token&,
                          std::optional<CheckpointPublicationDiagnostic>* = nullptr) noexcept;

/// Read a selected completed checkpoint on Windows or Linux. The last path component
/// must be its canonical checkpoint number. Every child is opened under one root handle.
[[nodiscard]] core::Result<LoadedCheckpointDirectory>
read_checkpoint_directory(const std::filesystem::path&, SourceId, const CheckpointFingerprints&,
                          const CheckpointManifestLimits&, CheckpointDirectoryLimits,
                          core::ResourcePlan&, const std::stop_token&) noexcept;

/// Enumerate the caller-selected parent under explicit quotas. Consider only canonical
/// numbered names, greatest first. Ignore corrupt/incompatible candidates and stop on
/// budget exhaustion or cancellation. The caller binds the parent pathname for the call.
[[nodiscard]] core::Result<LoadedCheckpointDirectory>
read_latest_checkpoint(const std::filesystem::path& parent, SourceId, const CheckpointFingerprints&,
                       const CheckpointManifestLimits&, CheckpointDirectoryLimits,
                       CheckpointDiscoveryLimits, core::ResourcePlan&,
                       const std::stop_token&) noexcept;

} // namespace uvlas::io::detail
