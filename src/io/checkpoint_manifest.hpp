#pragma once

#include "checkpoint_chunk.hpp"

#include <functional>

namespace uvlas::io::detail {

/// Opaque SHA-256 identifiers derived from the independently prepared receiving context.
struct CheckpointFingerprints {
    SourceDigest model;
    SourceDigest basis;
    SourceDigest material;
    SourceDigest constraint;
    SourceDigest algorithm;
    SourceDigest precision;
    SourceDigest version;
    bool operator==(const CheckpointFingerprints&) const = default;
};

struct CheckpointManifestLimits {
    InputLimits json;
    CheckpointChunkLimits chunk;
    std::uint64_t maximum_chunks;
    std::uint64_t maximum_total_payload_bytes;
};

struct NamedCheckpointChunk {
    std::string_view name;
    std::reference_wrapper<const CheckpointChunk> chunk;
};

struct CheckpointPayload {
    std::string_view name;
    std::span<const std::byte> bytes;
};

/// Entry borrows remain valid while the manifest has ownership. Payload files
/// are named chunk-<name>.bin and may be acquired separately by native adapters.
struct CheckpointManifestEntry {
    std::string_view name;
    CheckpointChunkHeader header;
};

struct CheckpointManifestStorage;

class CheckpointManifest {
public:
    CheckpointManifest(CheckpointManifest&&) noexcept;
    CheckpointManifest& operator=(CheckpointManifest&&) noexcept;
    ~CheckpointManifest();

    [[nodiscard]] static core::Result<CheckpointManifest>
    read(std::span<const std::byte>, SourceId, const CheckpointFingerprints& expected,
         const CheckpointManifestLimits&, core::ResourcePlan&, const std::stop_token&) noexcept;

    [[nodiscard]] core::Index size() const noexcept;
    [[nodiscard]] core::Result<CheckpointManifestEntry> at(core::Index) const& noexcept;
    [[nodiscard]] core::Result<CheckpointManifestEntry> at(core::Index) const&& = delete;

private:
    explicit CheckpointManifest(core::PlannedBuffer<CheckpointManifestStorage>) noexcept;
    core::PlannedBuffer<CheckpointManifestStorage> storage_;
};

class CheckpointPayloadSet {
public:
    CheckpointPayloadSet(CheckpointPayloadSet&&) noexcept = default;
    CheckpointPayloadSet& operator=(CheckpointPayloadSet&&) noexcept = default;

    [[nodiscard]] core::Index size() const noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CheckpointChunk>>
        at(core::Index) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CheckpointChunk>>
        at(core::Index) const&& = delete;

private:
    friend struct CheckpointDirectoryReader;
    friend core::Result<CheckpointPayloadSet>
    authenticate_checkpoint_payloads(const CheckpointManifest&, std::span<const CheckpointPayload>,
                                     CheckpointChunkLimits, core::ResourcePlan&,
                                     const std::stop_token&) noexcept;
    explicit CheckpointPayloadSet(core::PlannedBuffer<std::optional<CheckpointChunk>>) noexcept;
    core::PlannedBuffer<std::optional<CheckpointChunk>> chunks_;
};

/// Serialize verified chunks in lexical name order. The source byte cap also bounds
/// output; the remaining JSON decoder quotas apply when the receiver reads it.
/// Borrowed labels and owners must stay stable for the duration of the call.
[[nodiscard]] core::Result<core::PlannedBuffer<char>>
write_checkpoint_manifest(std::span<const NamedCheckpointChunk>, const CheckpointFingerprints&,
                          const CheckpointManifestLimits&, core::ResourcePlan&,
                          const std::stop_token&) noexcept;

/// Verify the whole declared set and publish independent chunks in manifest order.
/// No partial result escapes. Borrowed bytes and labels must stay stable during the call.
[[nodiscard]] core::Result<CheckpointPayloadSet>
authenticate_checkpoint_payloads(const CheckpointManifest&, std::span<const CheckpointPayload>,
                                 CheckpointChunkLimits, core::ResourcePlan&,
                                 const std::stop_token&) noexcept;

} // namespace uvlas::io::detail
