#pragma once

#include "checkpoint_manifest.hpp"
#include "checkpoint_publication.hpp"
#include "reference_input.hpp"

#include <filesystem>

namespace uvlas::io::detail {

struct LoadedCheckpointDirectory {
    std::uint64_t number;
    CheckpointManifest manifest;
    CheckpointPayloadSet payloads;
    ReferenceFileSnapshot manifest_file;
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
/// A supplied manifest digest pins the exact original bytes before parsing or reading chunks.
[[nodiscard]] core::Result<LoadedCheckpointDirectory>
read_checkpoint_directory(const std::filesystem::path&, SourceId, const CheckpointFingerprints&,
                          const CheckpointManifestLimits&, CheckpointDirectoryLimits,
                          core::ResourcePlan&, const std::stop_token&,
                          std::optional<SourceDigest> manifest_digest = {}) noexcept;

struct CaseCheckpointDirectoryLimits {
    ReferencePathLimits path;
    CheckpointDirectoryLimits container;
};

enum class CaseCheckpointDirectoryStage : std::uint8_t {
    Restart,
    Path,
    Open,
    ManifestDigest,
    Envelope,
};

struct CaseCheckpointDirectoryDiagnostic {
    std::optional<CaseCheckpointDirectoryStage> stage;
    std::optional<ReferencePathDiagnostic> path;
    std::optional<ManifestOpenDiagnostic> open;

    [[nodiscard]] std::string_view pointer() const noexcept;
};

/// Open only the canonical checkpoint named by the Case under its acquired root.
/// Authenticate the original manifest bytes before parsing or opening payloads.
/// The caller supplies expected fingerprints from the receiving context. The
/// resulting envelope is independent of the root handle and of the source files.
[[nodiscard]] core::Result<LoadedCheckpointDirectory> read_case_checkpoint_directory(
    const ReferenceInput&, SourceId manifest_source, const CheckpointFingerprints&,
    const CheckpointManifestLimits&, CaseCheckpointDirectoryLimits, core::ResourcePlan&,
    const std::stop_token&, CaseCheckpointDirectoryDiagnostic* = nullptr) noexcept;

/// Enumerate the caller-selected parent under explicit quotas. Consider only canonical
/// numbered names, greatest first. Ignore corrupt/incompatible candidates and stop on
/// budget exhaustion or cancellation. The caller binds the parent pathname for the call.
[[nodiscard]] core::Result<LoadedCheckpointDirectory>
read_latest_checkpoint(const std::filesystem::path& parent, SourceId, const CheckpointFingerprints&,
                       const CheckpointManifestLimits&, CheckpointDirectoryLimits,
                       CheckpointDiscoveryLimits, core::ResourcePlan&,
                       const std::stop_token&) noexcept;

} // namespace uvlas::io::detail
