#include "checkpoint_directory.hpp"

#include "checkpoint_directory_name.hpp"
#include "source_digest.hpp"

#if defined(_WIN32)
#include "windows_checkpoint_query.hpp"
#include "windows_reference_root.hpp"

#include <cstring>
#include <windows.h>
#elif defined(__linux__)
#include "linux_reference_root.hpp"

#include <cerrno>
#include <dirent.h>
#include <fcntl.h>
#include <memory>
#include <unistd.h>
#endif

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto recover = core::Operation::Recover;
constexpr auto corrupt = core::FailureCode::CheckpointCorrupt;
constexpr auto exceeded = core::FailureCode::ResourceExceeded;
constexpr auto cancelled = core::FailureCode::Cancelled;
constexpr std::string_view chunk_prefix = "chunk-";
constexpr std::string_view chunk_suffix = ".bin";

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

template <class Char>
[[nodiscard]] unsigned unsigned_code_unit(Char value) noexcept {
    return static_cast<std::make_unsigned_t<Char>>(value);
}

[[nodiscard]] std::optional<std::uint64_t> checkpoint_number(const std::filesystem::path& path) {
    const auto name = path.filename();
    const auto& native = name.native();
    std::array<char, 32> ascii{};
    if (native.size() >= ascii.size())
        return std::nullopt;
    for (std::size_t i = 0; i < native.size(); ++i) {
        if (unsigned_code_unit(native[i]) > 127)
            return std::nullopt;
        ascii[i] = static_cast<char>(native[i]);
    }
    return parse_checkpoint_directory_name({ascii.data(), native.size()});
}

[[noreturn]] void reject(core::FailureCode code) {
    throw core::Failure{code, recover, {}};
}

template <class Char>
[[nodiscard]] bool resembles_chunk(std::basic_string_view<Char> name) noexcept {
    const auto same_letter = [](Char actual, char expected) {
        const auto value = unsigned_code_unit(actual);
        const auto lower = value >= 'A' && value <= 'Z' ? value + ('a' - 'A') : value;
        return lower == unsigned_code_unit(expected);
    };
    return name.size() >= chunk_prefix.size() + chunk_suffix.size() &&
           std::ranges::equal(name.substr(0, chunk_prefix.size()), chunk_prefix, same_letter) &&
           std::ranges::equal(name.substr(name.size() - chunk_suffix.size()), chunk_suffix,
                              same_letter);
}

class DirectoryEntries {
public:
    DirectoryEntries(const CheckpointManifest& manifest, std::uint64_t limit,
                     const std::stop_token& stop) noexcept
        : manifest_(manifest), limit_(limit), stop_(stop) {}

    template <class Char>
    void inspect(std::basic_string_view<Char> name) {
        poll();
        if ((name.size() == 1 && name[0] == '.') ||
            (name.size() == 2 && name[0] == '.' && name[1] == '.')) {
            return;
        }
        if (examined_ == limit_)
            reject(exceeded);
        ++examined_;
        if (!resembles_chunk(name))
            return;
        std::array<char, 74> ascii{};
        if (name.size() > ascii.size())
            reject(corrupt);
        for (std::size_t i = 0; i < name.size(); ++i) {
            if (unsigned_code_unit(name[i]) > 127)
                reject(corrupt);
            ascii[i] = static_cast<char>(name[i]);
        }
        const std::string_view filename{ascii.data(), name.size()};
        if (!filename.starts_with(chunk_prefix) || !filename.ends_with(chunk_suffix))
            reject(corrupt);
        const auto label = filename.substr(
            chunk_prefix.size(), filename.size() - chunk_prefix.size() - chunk_suffix.size());
        core::Index first = 0;
        auto last = manifest_.size();
        while (first < last) {
            const auto middle = first + (last - first) / 2;
            if (checked(manifest_.at(middle)).name < label)
                first = middle + 1;
            else
                last = middle;
        }
        if (first == manifest_.size() || checked(manifest_.at(first)).name != label ||
            counted_ == manifest_.size())
            reject(corrupt);
        ++counted_;
    }

    void poll() const {
        if (stop_.stop_requested())
            reject(cancelled);
    }

    void finish() const {
        poll();
        if (counted_ != manifest_.size())
            reject(corrupt);
    }

private:
    const CheckpointManifest& manifest_;
    std::uint64_t limit_;
    const std::stop_token& stop_;
    std::uint64_t examined_ = 0;
    core::Index counted_ = 0;
};

#ifdef _WIN32
void inspect_directory(HANDLE directory, DirectoryEntries& entries) {
    constexpr auto header_bytes = offsetof(FILE_ID_BOTH_DIR_INFO, FileName);
    alignas(FILE_ID_BOTH_DIR_INFO) std::array<std::byte, 4096> buffer{};
    auto information = FileIdBothDirectoryRestartInfo;
    for (;;) {
        entries.poll();
        if (!GetFileInformationByHandleEx(directory, information, buffer.data(),
                                          static_cast<DWORD>(buffer.size()))) {
            const auto error = GetLastError();
            const auto failure = windows_checkpoint_query_failure(error);
            if (!failure)
                return;
            reject(*failure);
        }
        information = FileIdBothDirectoryInfo;
        std::span<const std::byte> record{buffer};
        for (;;) {
            if (record.size() < header_bytes)
                reject(corrupt);
            FILE_ID_BOTH_DIR_INFO info{};
            std::memcpy(&info, record.data(), header_bytes);
            const auto bytes = info.FileNameLength;
            if (bytes == 0 || bytes % sizeof(wchar_t) != 0 || bytes > record.size() - header_bytes)
                reject(corrupt);
            std::array<wchar_t, buffer.size() / sizeof(wchar_t)> name{};
            std::memcpy(name.data(), record.subspan(header_bytes).data(), bytes);
            entries.inspect(std::wstring_view{name.data(), bytes / sizeof(wchar_t)});
            const auto offset = info.NextEntryOffset;
            if (offset == 0)
                break;
            if (offset < header_bytes + bytes || offset >= record.size())
                reject(corrupt);
            record = record.subspan(offset);
        }
    }
}
#else
void inspect_directory(int directory, DirectoryEntries& entries) {
    const auto fd = ::openat(directory, ".", O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (fd < 0)
        reject(errno == ENOMEM || errno == EMFILE || errno == ENFILE ? exceeded : corrupt);
    auto* stream = ::fdopendir(fd);
    if (stream == nullptr) {
        const auto error = errno;
        const LinuxReferenceHandle owned{fd};
        reject(error == ENOMEM ? exceeded : corrupt);
    }
    const auto close = [](DIR* handle) noexcept { static_cast<void>(::closedir(handle)); };
    const std::unique_ptr<DIR, decltype(close)> owned{stream, close};
    for (;;) {
        entries.poll();
        errno = 0;
        // glibc readdir(3) is thread-safe for distinct streams. This scan owns its
        // stream exclusively and consumes the returned entry before the next call.
        // NOLINTNEXTLINE(concurrency-mt-unsafe)
        const auto* entry = ::readdir(owned.get());
        if (entry == nullptr) {
            const auto error = errno;
            if (error != 0)
                reject(error == ENOMEM ? exceeded : corrupt);
            return;
        }
        entries.inspect(std::string_view{entry->d_name});
    }
}
#endif

[[nodiscard]] ReferenceFileSnapshot
read_payload(const NativeRoot& root, const CheckpointManifestEntry& entry, ReferenceFileId manifest,
             ReferenceFileLimits limits, core::ResourcePlan& plan, const std::stop_token& stop) {
    std::array<char, 74> filename{};
    auto end = std::ranges::copy(chunk_prefix, filename.begin()).out;
    end = std::ranges::copy(entry.name, end).out;
    end = std::ranges::copy(chunk_suffix, end).out;
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
         CheckpointDirectoryLimits files, core::ResourcePlan& plan, const std::stop_token& stop,
         const std::optional<SourceDigest>& manifest_digest) {
        if (stop.stop_requested())
            throw core::Failure{core::FailureCode::Cancelled, recover, {}};
        const auto number = checked(checkpoint_number(directory));
        const auto root = checked(NativeRoot::open(directory, stop));
        auto manifest_handle = checked_file("manifest.json", root, plan, stop);
        const auto snapshot = checked(read_manifest_file(
            manifest_handle, bound_size(files.file, limits.json.maximum_source_bytes), plan, stop,
            nullptr));
        const auto original = checked(snapshot.view());
        if (manifest_digest && original.digest != *manifest_digest)
            throw core::Failure{corrupt, recover, {}};
        auto manifest =
            checked(CheckpointManifest::read(original.bytes, source, expected, limits, plan, stop));
        DirectoryEntries entries{manifest, files.maximum_entries, stop};
        inspect_directory(root.directory_.get(), entries);
        entries.finish();
        auto chunks =
            checked(plan.allocate<std::optional<CheckpointChunk>>(manifest.size(), recover));
        for (core::Index i = 0; i < manifest.size(); ++i) {
            const auto entry = checked(manifest.at(i));
            auto file =
                read_payload(root, entry, original.file_info.identity, files.file, plan, stop);
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
    core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<CheckpointPublicationDiagnostic>* diagnostic) noexcept {
    constexpr auto save = core::Operation::SaveCheckpoint;
    constexpr auto invalid = core::FailureCode::InvalidInput;
    if (diagnostic != nullptr)
        diagnostic->reset();
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
                                        std::as_bytes(manifest->values()), publication, plan, stop,
                                        diagnostic);
}

core::Result<LoadedCheckpointDirectory>
read_checkpoint_directory(const std::filesystem::path& directory, SourceId source,
                          const CheckpointFingerprints& expected,
                          const CheckpointManifestLimits& limits, CheckpointDirectoryLimits files,
                          core::ResourcePlan& plan, const std::stop_token& stop,
                          std::optional<SourceDigest> manifest_digest) noexcept {
    return core::guard_entry(
        recover, corrupt, {}, [&]() -> core::Result<LoadedCheckpointDirectory> {
            try {
                return CheckpointDirectoryReader::read(directory, source, expected, limits, files,
                                                       plan, stop, manifest_digest);
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            }
        });
}

core::Result<LoadedCheckpointDirectory> read_case_checkpoint_directory(
    const SchemaCase& input, const std::filesystem::path& resolved_directory,
    SourceId manifest_source, const CheckpointFingerprints& expected,
    const CheckpointManifestLimits& limits, CheckpointDirectoryLimits files,
    core::ResourcePlan& plan, const std::stop_token& stop) noexcept {
    if (stop.stop_requested())
        return std::unexpected(core::Failure{cancelled, recover, {}});
    const auto root = input.root();
    if (!root)
        return std::unexpected(root.error());
    if (!root->restart)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, recover, {}});
    const auto digest = decode_source_digest(root->restart->record().manifest_sha256.text());
    if (!digest)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, recover, {}});
    return read_checkpoint_directory(resolved_directory, manifest_source, expected, limits, files,
                                     plan, stop, digest);
}

core::Result<LoadedCheckpointDirectory>
read_latest_checkpoint(const std::filesystem::path& parent, SourceId source,
                       const CheckpointFingerprints& expected,
                       const CheckpointManifestLimits& limits, CheckpointDirectoryLimits files,
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
