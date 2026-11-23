#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stop_token>
#include <uvlas/io/input_document.hpp>

namespace uvlas::io::detail {

// File IDs are comparable only within the same native adapter and filesystem scope.
// The adapter populates the full object identifier, including the filesystem namespace.
struct ReferenceFileId {
    std::uint64_t filesystem;
    std::array<std::byte, 16> object;
    bool operator==(const ReferenceFileId&) const = default;
};

struct ReferenceFileInfo {
    ReferenceFileId identity;
    std::uint64_t byte_count;
    // Opaque adapter representation of all qualified modification indicators. No zero fallback.
    std::array<std::uint64_t, 4> change_indicators;
    bool regular_file;
    bool operator==(const ReferenceFileInfo&) const = default;
};

// Private native adapter boundary. Methods operate on the same already-open file, with its
// cursor initially at byte zero. The caller retains and closes the handle. read() writes
// only inside the supplied span and returns its actual byte count, with zero denoting EOF.
class ReferenceFileReader {
public:
    virtual ~ReferenceFileReader() = default;
    [[nodiscard]] virtual core::Result<ReferenceFileInfo> info() noexcept = 0;
    [[nodiscard]] virtual core::Result<std::size_t> read(std::span<std::byte>) noexcept = 0;
};

struct ReferenceFileLimits {
    std::uint64_t maximum_file_bytes;
    std::uint64_t maximum_read_bytes; // Includes one byte of capacity for the EOF probe.
    std::uint64_t maximum_read_calls; // Counts adapter reads, including the probe. The adapter
                                      // budgets any native retries separately.
    std::uint64_t maximum_hash_bytes;
};

struct ReferenceFileIdentityChecks {
    // IDs are from the same native origin and from already authenticated sources.
    std::optional<ReferenceFileId> manifest{};
    std::span<const ReferenceFileId> earlier_distinct_sources{};
    std::uint64_t maximum_comparisons = 0;
};

enum class ReferenceFileFailure : std::uint8_t {
    Metadata,
    NotRegular,
    FileBytes,
    ReadBytes,
    ReadCalls,
    HashBytes,
    PlanBytes,
    Read,
    InvalidReadCount,
    PrematureEof,
    TrailingByte,
    Modified,
    DigestMismatch,
    Cycle,
    DuplicateFile,
    IdentityComparisons,
};

struct ReferenceFileQuota {
    std::uint64_t required;
    std::uint64_t limit;
};

struct ReferenceFileDiagnostic {
    ReferenceFileDiagnostic(ReferenceFileFailure failure, std::uint64_t offset,
                            std::optional<ReferenceFileInfo> initial = {},
                            std::optional<ReferenceFileInfo> final = {},
                            std::optional<ReferenceFileQuota> limit = {}) noexcept
        : reason(failure), bytes_read(offset), before(initial), after(final), quota(limit) {}

    ReferenceFileFailure reason;
    std::uint64_t bytes_read; // Payload bytes read before failure, excluding the one-byte probe.
    std::optional<ReferenceFileInfo> before{};
    std::optional<ReferenceFileInfo> after{};
    std::optional<ReferenceFileQuota> quota{}; // Units are given by reason: bytes or calls.
};

struct ReferenceFileView {
    ReferenceFileView(std::span<const std::byte> payload, SourceDigest verified_digest,
                      ReferenceFileInfo verified_info) noexcept
        : bytes(payload), digest(verified_digest), file_info(verified_info) {}

    std::span<const std::byte> bytes;
    SourceDigest digest;
    ReferenceFileInfo file_info;
};

// Owns the exact payload read from a stable handle and the SHA-256 of those bytes.
// Views borrow the owner until it is moved or destroyed.
class ReferenceFileSnapshot {
public:
    ReferenceFileSnapshot(ReferenceFileSnapshot&& other) noexcept;
    ReferenceFileSnapshot& operator=(ReferenceFileSnapshot&& other) noexcept;
    ReferenceFileSnapshot(const ReferenceFileSnapshot&) = delete;
    ReferenceFileSnapshot& operator=(const ReferenceFileSnapshot&) = delete;
    ~ReferenceFileSnapshot() = default;
    [[nodiscard]] std::optional<ReferenceFileView> view() const& noexcept;
    [[nodiscard]] std::optional<ReferenceFileView> view() const&& = delete;

private:
    friend class CheckpointChunk;
    friend core::Result<ReferenceFileSnapshot>
    read_reference_file(ReferenceFileReader&, const SourceDigest&, ReferenceFileLimits,
                        core::ResourcePlan&, const std::stop_token&,
                        std::optional<ReferenceFileDiagnostic>*,
                        const ReferenceFileIdentityChecks&) noexcept;
    friend core::Result<ReferenceFileSnapshot>
    read_manifest_file(ReferenceFileReader&, ReferenceFileLimits, core::ResourcePlan&,
                       const std::stop_token&, std::optional<ReferenceFileDiagnostic>*) noexcept;

    [[nodiscard]] static core::Result<ReferenceFileSnapshot>
    read_from_handle(ReferenceFileReader&, const SourceDigest* declared_digest, ReferenceFileLimits,
                     core::ResourcePlan&, const std::stop_token&,
                     std::optional<ReferenceFileDiagnostic>*,
                     const ReferenceFileIdentityChecks&) noexcept;

    ReferenceFileSnapshot(core::PlannedBuffer<std::byte> payload, SourceDigest verified_digest,
                          ReferenceFileInfo verified_info) noexcept;
    core::PlannedBuffer<std::byte> bytes_;
    SourceDigest digest_;
    std::optional<ReferenceFileInfo> info_;
};

// Reads and verifies one handle. On failure, successful allocations stay charged but no
// partial payload escapes. Diagnostic reports observed progress and clears on success or cancel.
// Manifest and earlier source identities are checked against the very first metadata, before
// payload admission, allocation or reading.
// Adapter and ResourcePlan must have exclusive use during this synchronous call.
[[nodiscard]] core::Result<ReferenceFileSnapshot>
read_reference_file(ReferenceFileReader& reader, const SourceDigest& declared_digest,
                    ReferenceFileLimits limits, core::ResourcePlan& plan,
                    const std::stop_token& stop, std::optional<ReferenceFileDiagnostic>* diagnostic,
                    const ReferenceFileIdentityChecks& identities = {}) noexcept;

// Reads the manifest itself, which has no enclosing declaration to supply an expected digest.
// The same length, identity, modification, cancellation and resource checks apply. The result
// owns the original bytes, their computed digest and the observed file identity. The caller
// keeps the reader's handle alive separately when the manifest identity must stay pinned.
[[nodiscard]] core::Result<ReferenceFileSnapshot>
read_manifest_file(ReferenceFileReader& reader, ReferenceFileLimits limits,
                   core::ResourcePlan& plan, const std::stop_token& stop,
                   std::optional<ReferenceFileDiagnostic>* diagnostic) noexcept;

} // namespace uvlas::io::detail
