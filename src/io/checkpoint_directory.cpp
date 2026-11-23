#include "checkpoint_directory.hpp"

#include "checkpoint_directory_name.hpp"

#if defined(_WIN32)
#include "windows_reference_root.hpp"
#elif defined(__linux__)
#include "linux_reference_root.hpp"
#endif

#include <algorithm>
#include <array>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto recover = core::Operation::Recover;
constexpr auto corrupt = core::FailureCode::CheckpointCorrupt;

#if defined(_WIN32)
using NativeRoot = WindowsReferenceRoot;
#else
using NativeRoot = LinuxReferenceRoot;
#endif

[[nodiscard]] core::Failure checkpoint_error(core::Failure failure) noexcept {
    failure.operation = recover;
    if (failure.code == core::FailureCode::InvalidInput) {
        failure.code = corrupt;
    }
    return failure;
}

template <class T>
[[nodiscard]] T checked(core::Result<T> input) {
    if (!input)
        throw checkpoint_error(input.error());
    return std::move(*input);
}

template <class T>
[[nodiscard]] T checked(std::optional<T> input) {
    if (!input)
        throw core::Failure{corrupt, recover, {}};
    return std::move(*input);
}

[[nodiscard]] auto checked_file(std::string_view name, const NativeRoot& root,
                                core::ResourcePlan& plan, const std::stop_token& stop) {
    const auto relative =
        checked(normalize_reference_path(name, {name.size(), name.size(), 1}, plan, stop, nullptr));
    return checked(root.open_file(relative, plan, stop));
}

[[nodiscard]] ReferenceFileLimits bound_size(ReferenceFileLimits limits,
                                             std::uint64_t maximum) noexcept {
    limits.maximum_file_bytes = std::min(limits.maximum_file_bytes, maximum);
    return limits;
}

[[nodiscard]] std::optional<std::uint64_t> checkpoint_number(const std::filesystem::path& path) {
    const auto name = path.filename();
    const auto& native = name.native();
    std::array<char, 32> ascii{};
    if (native.size() >= ascii.size())
        return std::nullopt;
    for (std::size_t i = 0; i < native.size(); ++i) {
        if (static_cast<std::uint32_t>(native[i]) > 127)
            return std::nullopt;
        ascii[i] = static_cast<char>(native[i]);
    }
    return parse_checkpoint_directory_name({ascii.data(), native.size()});
}

[[nodiscard]] ReferenceFileSnapshot
read_payload(const NativeRoot& root, const CheckpointManifestEntry& entry, ReferenceFileId manifest,
             ReferenceFileLimits limits, core::ResourcePlan& plan, const std::stop_token& stop) {
    constexpr auto prefix = std::string_view{"chunk-"};
    constexpr auto suffix = std::string_view{".bin"};
    std::array<char, 74> filename{};
    auto end = std::ranges::copy(prefix, filename.begin()).out;
    end = std::ranges::copy(entry.name, end).out;
    end = std::ranges::copy(suffix, end).out;
    const std::string_view name{filename.data(), static_cast<std::size_t>(end - filename.begin())};
    auto file = checked_file(name, root, plan, stop);
    const auto info = checked(file.info());
    if (info.byte_count != entry.header.byte_count)
        throw core::Failure{corrupt, recover, {}};
    return checked(read_reference_file(file, entry.header.digest,
                                       bound_size(limits, entry.header.byte_count), plan, stop,
                                       nullptr, {manifest, {}, 1}));
}

} // namespace

struct CheckpointDirectoryReader {
    [[nodiscard]] static LoadedCheckpointDirectory
    read(const std::filesystem::path& directory, SourceId source,
         const CheckpointFingerprints& expected, const CheckpointManifestLimits& limits,
         ReferenceFileLimits files, core::ResourcePlan& plan, const std::stop_token& stop) {
        if (stop.stop_requested())
            throw core::Failure{core::FailureCode::Cancelled, recover, {}};
        const auto number = checked(checkpoint_number(directory));
        const auto root = checked(NativeRoot::open(directory, stop));
        auto manifest_handle = checked_file("manifest.json", root, plan, stop);
        const auto snapshot = checked(
            read_manifest_file(manifest_handle, bound_size(files, limits.json.maximum_source_bytes),
                               plan, stop, nullptr));
        const auto original = checked(snapshot.view());
        auto manifest =
            checked(CheckpointManifest::read(original.bytes, source, expected, limits, plan, stop));
        auto chunks =
            checked(plan.allocate<std::optional<CheckpointChunk>>(manifest.size(), recover));
        for (core::Index i = 0; i < manifest.size(); ++i) {
            const auto entry = checked(manifest.at(i));
            auto file = read_payload(root, entry, original.file_info.identity, files, plan, stop);
            chunks.values()[static_cast<std::size_t>(i)].emplace(
                checked(CheckpointChunk::from_verified_file(entry.header, std::move(file),
                                                            limits.chunk, plan, stop)));
        }
        if (stop.stop_requested())
            throw core::Failure{core::FailureCode::Cancelled, recover, {}};
        return {number, std::move(manifest), CheckpointPayloadSet{std::move(chunks)}};
    }
};

core::Result<void> save_checkpoint_directory(
    const std::filesystem::path& parent, std::uint64_t number,
    std::span<const NamedCheckpointChunk> chunks, const CheckpointFingerprints& identity,
    const CheckpointManifestLimits& limits, CheckpointPublicationLimits publication,
    core::ResourcePlan& plan, const std::stop_token& stop) noexcept {
    constexpr auto save = core::Operation::SaveCheckpoint;
    constexpr auto invalid = core::FailureCode::InvalidInput;
    auto manifest = write_checkpoint_manifest(chunks, identity, limits, plan, stop);
    if (!manifest)
        return std::unexpected(manifest.error());
    auto entries =
        plan.allocate<CheckpointPublicationChunk>(static_cast<core::Index>(chunks.size()), save);
    if (!entries)
        return std::unexpected(entries.error());
    for (std::size_t i = 0; i < chunks.size(); ++i) {
        const auto& chunk = chunks[i];
        const auto view = chunk.chunk.get().view();
        if (!view)
            return std::unexpected(core::Failure{invalid, save, {}});
        entries->values()[i] = {chunk.name, *view};
    }
    return publish_checkpoint_directory(parent, number, entries->values(),
                                        std::as_bytes(manifest->values()), publication, plan, stop);
}

core::Result<LoadedCheckpointDirectory>
read_checkpoint_directory(const std::filesystem::path& directory, SourceId source,
                          const CheckpointFingerprints& expected,
                          const CheckpointManifestLimits& limits, ReferenceFileLimits files,
                          core::ResourcePlan& plan, const std::stop_token& stop) noexcept {
    return core::guard_entry(recover, corrupt, {},
                             [&]() -> core::Result<LoadedCheckpointDirectory> {
                                 try {
                                     return CheckpointDirectoryReader::read(
                                         directory, source, expected, limits, files, plan, stop);
                                 } catch (const core::Failure& failure) {
                                     return std::unexpected(failure);
                                 }
                             });
}

core::Result<LoadedCheckpointDirectory>
read_latest_checkpoint(const std::filesystem::path& parent, SourceId source,
                       const CheckpointFingerprints& expected,
                       const CheckpointManifestLimits& limits, ReferenceFileLimits files,
                       CheckpointDiscoveryLimits discovery, core::ResourcePlan& plan,
                       const std::stop_token& stop) noexcept {
    using core::FailureCode;
    return core::guard_entry(
        recover, corrupt, {}, [&]() -> core::Result<LoadedCheckpointDirectory> {
            if (stop.stop_requested())
                return std::unexpected(core::Failure{FailureCode::Cancelled, recover, {}});
            const auto capacity = std::min(discovery.maximum_entries, discovery.maximum_candidates);
            if (!std::in_range<core::Index>(capacity))
                return std::unexpected(core::Failure{FailureCode::ResourceExceeded, recover, {}});
            auto candidates =
                plan.allocate<std::uint64_t>(static_cast<core::Index>(capacity), recover);
            if (!candidates)
                return std::unexpected(candidates.error());
            std::uint64_t examined = 0;
            std::size_t used = 0;
            for (const auto& entry : std::filesystem::directory_iterator(parent)) {
                if (stop.stop_requested())
                    return std::unexpected(core::Failure{FailureCode::Cancelled, recover, {}});
                if (examined == discovery.maximum_entries)
                    return std::unexpected(
                        core::Failure{FailureCode::ResourceExceeded, recover, {}});
                ++examined;
                const auto number = checkpoint_number(entry.path());
                if (number) {
                    if (used == candidates->size())
                        return std::unexpected(
                            core::Failure{FailureCode::ResourceExceeded, recover, {}});
                    candidates->values()[used++] = *number;
                }
            }
            auto numbered = candidates->values().first(used);
            std::ranges::sort(numbered, std::greater<>{});
            std::optional<core::Failure> newest_failure;
            for (const auto number : numbered) {
                const auto name = make_checkpoint_directory_name(number);
                auto candidate = read_checkpoint_directory(parent / name.data(), source, expected,
                                                           limits, files, plan, stop);
                if (candidate)
                    return std::move(*candidate);
                const auto failure = candidate.error();
                if (failure.code != corrupt && failure.code != FailureCode::CheckpointIncompatible)
                    return std::unexpected(failure);
                if (!newest_failure)
                    newest_failure = failure;
            }
            if (stop.stop_requested())
                return std::unexpected(core::Failure{FailureCode::Cancelled, recover, {}});
            return std::unexpected(newest_failure.value_or(core::Failure{corrupt, recover, {}}));
        });
}

} // namespace uvlas::io::detail
