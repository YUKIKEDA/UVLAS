#pragma once

#include "checkpoint_chunk.hpp"

#include <filesystem>

namespace uvlas::io::detail {

struct CheckpointPublicationChunk {
    std::string_view logical_name;
    CheckpointChunkView chunk;
};

struct CheckpointPublicationLimits {
    std::uint64_t maximum_chunks;
    std::uint64_t maximum_file_bytes;
    std::uint64_t maximum_total_bytes;
};

enum class CheckpointPublicationStage : std::uint8_t {
    Preflight,
    OpenParent,
    StageDirectory,
    Create,
    Write,
    Flush,
    Close,
    ReadBack,
    Verify,
    Publish,
    DestinationExists
};

struct CheckpointPublicationDiagnostic {
    CheckpointPublicationStage stage;
    std::optional<std::size_t> chunk;
    std::optional<std::uint32_t> native_error;
};

// Publish an already validated container in a caller-managed output directory. The caller
// keeps the input stable until return and supplies the committed sequence number for the directory
// name. The native publisher does not interpret the manifest or advance a solver state.
// Interruption before the rename can leave a temporary name but never publishes a
// checkpoint. The operation does not claim retention across hardware or volume failures.
[[nodiscard]] core::Result<void> publish_checkpoint_directory(
    const std::filesystem::path& parent, std::uint64_t checkpoint_number,
    std::span<const CheckpointPublicationChunk> chunks, std::span<const std::byte> manifest,
    CheckpointPublicationLimits limits, core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<CheckpointPublicationDiagnostic>* diagnostic = nullptr) noexcept;

} // namespace uvlas::io::detail
