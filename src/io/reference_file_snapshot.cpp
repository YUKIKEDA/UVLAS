#include "reference_file_snapshot.hpp"

#include "source_digest.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
constexpr std::size_t read_chunk = 65'536;

[[nodiscard]] std::unexpected<core::Failure> cancelled() noexcept {
    return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
}

} // namespace

ReferenceFileSnapshot::ReferenceFileSnapshot(core::PlannedBuffer<std::byte> payload,
                                             SourceDigest verified_digest,
                                             ReferenceFileInfo verified_info) noexcept
    : bytes_(std::move(payload)), digest_(verified_digest), info_(verified_info) {}

ReferenceFileSnapshot::ReferenceFileSnapshot(ReferenceFileSnapshot&& other) noexcept
    : bytes_(std::move(other.bytes_)),
      digest_(other.digest_),
      info_(std::exchange(other.info_, std::nullopt)) {}

ReferenceFileSnapshot& ReferenceFileSnapshot::operator=(ReferenceFileSnapshot&& other) noexcept {
    if (this != &other) {
        bytes_ = std::move(other.bytes_);
        digest_ = other.digest_;
        info_ = std::exchange(other.info_, std::nullopt);
    }
    return *this;
}

std::optional<ReferenceFileView> ReferenceFileSnapshot::view() const& noexcept {
    const auto stable_info = info_;
    if (!stable_info) {
        return std::nullopt;
    }
    return ReferenceFileView{bytes_.values(), digest_, *stable_info};
}

core::Result<ReferenceFileSnapshot> ReferenceFileSnapshot::read_from_handle(
    ReferenceFileReader& reader, const SourceDigest* declared_digest, ReferenceFileLimits limits,
    core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<ReferenceFileDiagnostic>* diagnostic,
    const ReferenceFileIdentityChecks& identities) noexcept {
    if (diagnostic != nullptr) {
        diagnostic->reset();
    }
    std::size_t offset = 0;
    std::optional<ReferenceFileInfo> observed_before;
    std::optional<ReferenceFileInfo> observed_after;
    const auto describe = [&](ReferenceFileFailure reason,
                              std::optional<ReferenceFileQuota> quota = std::nullopt) noexcept {
        if (diagnostic != nullptr) {
            *diagnostic =
                ReferenceFileDiagnostic{reason, offset, observed_before, observed_after, quota};
        }
    };
    const auto reject = [&](ReferenceFileFailure reason) noexcept {
        describe(reason);
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    };
    const auto budget = [&](ReferenceFileFailure reason, std::uint64_t required,
                            std::uint64_t limit) noexcept {
        describe(reason, ReferenceFileQuota{required, limit});
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, operation});
    };
    const auto adapter_error = [&](ReferenceFileFailure reason, core::Failure failure) noexcept {
        if (failure.code != core::FailureCode::Cancelled) {
            describe(reason);
        }
        failure.operation = operation;
        return std::unexpected(failure);
    };

    if (stop.stop_requested()) {
        return cancelled();
    }
    const auto before = reader.info();
    if (stop.stop_requested()) {
        return cancelled();
    }
    if (!before) {
        return adapter_error(ReferenceFileFailure::Metadata, before.error());
    }
    observed_before = *before;
    if (!before->regular_file) {
        return reject(ReferenceFileFailure::NotRegular);
    }

    // Bound all native-ID comparisons without hiding a mismatch already found. The manifest
    // precedes the prior sources, and both precede the payload budgets below.
    const std::uint64_t manifest_comparisons = identities.manifest.has_value() ? 1 : 0;
    if (manifest_comparisons > identities.maximum_comparisons) {
        return budget(ReferenceFileFailure::IdentityComparisons, manifest_comparisons,
                      identities.maximum_comparisons);
    }
    if (identities.manifest && before->identity == *identities.manifest) {
        return reject(ReferenceFileFailure::Cycle);
    }
    const auto prior = identities.earlier_distinct_sources;
    const auto available = identities.maximum_comparisons - manifest_comparisons;
    const auto admitted = std::min<std::uint64_t>(prior.size(), available);
    for (const auto& seen : prior.first(static_cast<std::size_t>(admitted))) {
        if (stop.stop_requested()) {
            return cancelled();
        }
        if (before->identity == seen) {
            return reject(ReferenceFileFailure::DuplicateFile);
        }
    }
    if (admitted < prior.size()) {
        // The number of addressable ID records leaves room for the manifest and this probe.
        return budget(ReferenceFileFailure::IdentityComparisons,
                      manifest_comparisons + admitted + 1, identities.maximum_comparisons);
    }
    if (stop.stop_requested()) {
        return cancelled();
    }

    const auto length = before->byte_count;
    constexpr auto index_max = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    constexpr auto host_max = static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max());
    const auto file_limit = std::min({limits.maximum_file_bytes, index_max, host_max});
    if (length > file_limit) {
        return budget(ReferenceFileFailure::FileBytes, length, file_limit);
    }
    if (length + 1 > limits.maximum_read_bytes) {
        return budget(ReferenceFileFailure::ReadBytes, length + 1, limits.maximum_read_bytes);
    }
    // Each successful chunk advances by at most read_chunk. Budget at least the inevitable
    // calls before allocating; short reads may consume more and are checked as they occur.
    const auto minimum_calls = length / read_chunk + (length % read_chunk != 0 ? 1U : 0U) + 1;
    if (minimum_calls > limits.maximum_read_calls) {
        return budget(ReferenceFileFailure::ReadCalls, minimum_calls, limits.maximum_read_calls);
    }
    // SHA-256 encodes the original message length as an unsigned 64-bit bit count.
    constexpr auto sha256_max = std::numeric_limits<std::uint64_t>::max() / 8;
    const auto hash_limit = std::min(limits.maximum_hash_bytes, sha256_max);
    if (length > hash_limit) {
        return budget(ReferenceFileFailure::HashBytes, length, hash_limit);
    }
    auto storage = plan.allocate<std::byte>(static_cast<core::Index>(length), operation);
    if (!storage) {
        const auto usage = storage.error().details.resources;
        const auto quota =
            usage ? std::optional{ReferenceFileQuota{usage->required, plan.byte_limit()}}
                  : std::nullopt;
        describe(ReferenceFileFailure::PlanBytes, quota);
        return std::unexpected(storage.error());
    }
    const auto bytes = storage->values();
    auto remaining_calls = limits.maximum_read_calls;
    const auto read = [&](std::span<std::byte> destination) noexcept -> core::Result<std::size_t> {
        if (stop.stop_requested()) {
            return cancelled();
        }
        if (remaining_calls == 0) {
            // At most one read per payload byte plus the probe can advance this far.
            return budget(ReferenceFileFailure::ReadCalls, limits.maximum_read_calls + 1,
                          limits.maximum_read_calls);
        }
        --remaining_calls;
        const auto count = reader.read(destination);
        if (stop.stop_requested()) {
            return cancelled();
        }
        if (!count) {
            return adapter_error(ReferenceFileFailure::Read, count.error());
        }
        if (*count > destination.size()) {
            return reject(ReferenceFileFailure::InvalidReadCount);
        }
        return *count;
    };
    while (offset < bytes.size()) {
        const auto requested = std::min(read_chunk, bytes.size() - offset);
        const auto count = read(bytes.subspan(offset, requested));
        if (!count) {
            return std::unexpected(count.error());
        }
        if (*count == 0) {
            return reject(ReferenceFileFailure::PrematureEof);
        }
        offset += *count;
    }
    std::array<std::byte, 1> excess{};
    const auto tail = read(excess);
    if (!tail) {
        return std::unexpected(tail.error());
    }
    if (*tail != 0) {
        return reject(ReferenceFileFailure::TrailingByte);
    }

    if (stop.stop_requested()) {
        return cancelled();
    }
    const auto after = reader.info();
    if (stop.stop_requested()) {
        return cancelled();
    }
    if (!after) {
        return adapter_error(ReferenceFileFailure::Metadata, after.error());
    }
    observed_after = *after;
    if (*after != *before) {
        return reject(ReferenceFileFailure::Modified);
    }
    const auto digest = source_digest(bytes);
    if (stop.stop_requested()) {
        return cancelled();
    }
    if (declared_digest != nullptr && digest != *declared_digest) {
        return reject(ReferenceFileFailure::DigestMismatch);
    }
    return ReferenceFileSnapshot{std::move(*storage), digest, *before};
}

core::Result<ReferenceFileSnapshot>
read_reference_file(ReferenceFileReader& reader, const SourceDigest& declared_digest,
                    ReferenceFileLimits limits, core::ResourcePlan& plan,
                    const std::stop_token& stop, std::optional<ReferenceFileDiagnostic>* diagnostic,
                    const ReferenceFileIdentityChecks& identities) noexcept {
    return ReferenceFileSnapshot::read_from_handle(reader, &declared_digest, limits, plan, stop,
                                                   diagnostic, identities);
}

core::Result<ReferenceFileSnapshot>
read_manifest_file(ReferenceFileReader& reader, ReferenceFileLimits limits,
                   core::ResourcePlan& plan, const std::stop_token& stop,
                   std::optional<ReferenceFileDiagnostic>* diagnostic) noexcept {
    return ReferenceFileSnapshot::read_from_handle(reader, nullptr, limits, plan, stop, diagnostic,
                                                   {});
}

} // namespace uvlas::io::detail
