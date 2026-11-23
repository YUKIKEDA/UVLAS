#include "checkpoint_publication.hpp"

#include "checkpoint_directory_name.hpp"
#include "source_digest.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <limits>
#include <utility>
#ifdef _WIN32
#include <windows.h>
#include <winternl.h>

// The user-mode export is in ntdll; the SDK's winternl.h omits its declaration.
extern "C" NTSYSAPI NTSTATUS NTAPI NtSetInformationFile(HANDLE, PIO_STATUS_BLOCK, PVOID, ULONG,
                                                        FILE_INFORMATION_CLASS);
#else
#include <cerrno>
#include <fcntl.h>
#include <linux/fs.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

namespace uvlas::io::detail {
namespace {

constexpr auto save = core::Operation::SaveCheckpoint;
using Code = core::FailureCode;
using Stage = CheckpointPublicationStage;
using Name = std::array<char, 80>;
constexpr std::size_t block_size = 64 * 1024;

[[nodiscard]] Name name_with(std::string_view prefix, std::string_view middle,
                             std::string_view suffix) noexcept {
    Name name{};
    auto end = std::ranges::copy(prefix, name.begin()).out;
    end = std::ranges::copy(middle, end).out;
    std::ranges::copy(suffix, end);
    return name;
}

[[nodiscard]] Name chunk_name(std::string_view logical_name) noexcept {
    return name_with("chunk-", logical_name, ".bin");
}

[[nodiscard]] bool safe_chunk_name(std::string_view name) noexcept {
    const auto letter = [](char c) { return c >= 'a' && c <= 'z'; };
    return !name.empty() && name.size() <= 64 && letter(name.front()) &&
           std::ranges::all_of(name, [&](char c) {
               return letter(c) || (c >= '0' && c <= '9') || c == '_' || c == '-';
           });
}

template <class T>
using NativeResult = std::expected<T, std::uint32_t>;
using NativeStatus = NativeResult<void>;

[[nodiscard]] bool interrupted(std::uint32_t error) noexcept {
#ifdef _WIN32
    return error == ERROR_OPERATION_ABORTED;
#else
    return error == EINTR;
#endif
}

#ifdef _WIN32
using RawHandle = HANDLE;
const RawHandle invalid_handle = nullptr;

[[nodiscard]] std::uint32_t native_error() noexcept {
    return GetLastError();
}

[[nodiscard]] NativeStatus close_native(RawHandle handle) noexcept {
    if (!CloseHandle(handle)) {
        return std::unexpected(native_error());
    }
    return {};
}

[[nodiscard]] bool exists_error(std::uint32_t error) noexcept {
    return error == ERROR_ALREADY_EXISTS || error == ERROR_FILE_EXISTS;
}

[[nodiscard]] Code native_code(std::uint32_t error) noexcept {
    if (error == ERROR_DISK_FULL || error == ERROR_HANDLE_DISK_FULL ||
        error == ERROR_NOT_ENOUGH_MEMORY || error == ERROR_OUTOFMEMORY) {
        return Code::ResourceExceeded;
    }
    if (error == ERROR_NOT_SUPPORTED || error == ERROR_INVALID_FUNCTION ||
        error == ERROR_CALL_NOT_IMPLEMENTED) {
        return Code::Unsupported;
    }
    return Code::InvalidInput;
}

[[nodiscard]] std::array<wchar_t, Name{}.size()> wide_name(const Name& name) noexcept {
    std::array<wchar_t, Name{}.size()> wide{};
    std::ranges::transform(name, wide.begin(), [](char c) { return static_cast<wchar_t>(c); });
    return wide;
}

[[nodiscard]] bool valid_parent(RawHandle root) noexcept {
    FILE_ATTRIBUTE_TAG_INFO info{};
    return GetFileType(root) == FILE_TYPE_DISK &&
           GetFileInformationByHandleEx(root, FileAttributeTagInfo, &info, sizeof(info)) &&
           (info.FileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0 &&
           (info.FileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) == 0;
}

[[nodiscard]] NativeResult<RawHandle> open_parent(const std::filesystem::path& parent) noexcept {
    constexpr DWORD access = FILE_LIST_DIRECTORY | FILE_ADD_SUBDIRECTORY | FILE_TRAVERSE |
                             FILE_READ_ATTRIBUTES | SYNCHRONIZE;
    const auto handle = CreateFileW(
        parent.c_str(), access, FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, OPEN_EXISTING,
        FILE_FLAG_OPEN_REPARSE_POINT | FILE_FLAG_BACKUP_SEMANTICS, nullptr);
    if (handle == INVALID_HANDLE_VALUE) {
        return std::unexpected(native_error());
    }
    if (!valid_parent(handle)) {
        static_cast<void>(close_native(handle));
        return std::unexpected(ERROR_DIRECTORY);
    }
    return handle;
}

[[nodiscard]] NativeResult<RawHandle> open_relative(RawHandle parent, const Name& name,
                                                    ACCESS_MASK access, ULONG disposition,
                                                    ULONG type) noexcept {
    auto text = wide_name(name);
    const auto length =
        static_cast<USHORT>(std::char_traits<char>::length(name.data()) * sizeof(wchar_t));
    UNICODE_STRING native{length, length, text.data()};
    OBJECT_ATTRIBUTES attributes{sizeof(OBJECT_ATTRIBUTES),
                                 parent,
                                 &native,
                                 OBJ_CASE_INSENSITIVE | OBJ_DONT_REPARSE,
                                 nullptr,
                                 nullptr};
    IO_STATUS_BLOCK info{};
    HANDLE handle = nullptr;
    const auto status =
        NtCreateFile(&handle, access | SYNCHRONIZE, &attributes, &info, nullptr,
                     FILE_ATTRIBUTE_NORMAL, FILE_SHARE_READ | FILE_SHARE_WRITE, disposition,
                     type | FILE_OPEN_REPARSE_POINT | FILE_SYNCHRONOUS_IO_NONALERT, nullptr, 0);
    if (status < 0) {
        return std::unexpected(RtlNtStatusToDosError(status));
    }
    return handle;
}

[[nodiscard]] auto create_stage(RawHandle parent, const Name& name) noexcept {
    return open_relative(parent, name,
                         FILE_LIST_DIRECTORY | FILE_ADD_FILE | FILE_TRAVERSE |
                             FILE_READ_ATTRIBUTES | DELETE,
                         FILE_CREATE, FILE_DIRECTORY_FILE);
}

[[nodiscard]] auto open_output(RawHandle parent, const Name& name) noexcept {
    return open_relative(parent, name, GENERIC_WRITE | FILE_READ_ATTRIBUTES, FILE_CREATE,
                         FILE_NON_DIRECTORY_FILE);
}

[[nodiscard]] auto open_readback(RawHandle parent, const Name& name) noexcept {
    return open_relative(parent, name, FILE_READ_DATA | FILE_READ_ATTRIBUTES, FILE_OPEN,
                         FILE_NON_DIRECTORY_FILE);
}

struct FileInfo {
    FILE_ID_INFO id;
    std::uint64_t length;
    [[nodiscard]] bool matches(const FileInfo& other) const noexcept {
        return id.VolumeSerialNumber == other.id.VolumeSerialNumber &&
               std::ranges::equal(id.FileId.Identifier, other.id.FileId.Identifier);
    }
};

[[nodiscard]] NativeResult<FileInfo> file_info(RawHandle file) noexcept {
    FILE_STANDARD_INFO standard{};
    FILE_ATTRIBUTE_TAG_INFO attributes{};
    FileInfo info{};
    if (!GetFileInformationByHandleEx(file, FileIdInfo, &info.id, sizeof(info.id)) ||
        !GetFileInformationByHandleEx(file, FileStandardInfo, &standard, sizeof(standard)) ||
        !GetFileInformationByHandleEx(file, FileAttributeTagInfo, &attributes,
                                      sizeof(attributes))) {
        return std::unexpected(native_error());
    }
    if (standard.Directory || standard.EndOfFile.QuadPart < 0 ||
        (attributes.FileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0 ||
        GetFileType(file) != FILE_TYPE_DISK) {
        return std::unexpected(ERROR_INVALID_DATA);
    }
    info.length = static_cast<std::uint64_t>(standard.EndOfFile.QuadPart);
    return info;
}

[[nodiscard]] NativeResult<std::size_t> write_block(RawHandle file,
                                                    std::span<const std::byte> bytes) noexcept {
    DWORD count = 0;
    if (!WriteFile(file, bytes.data(), static_cast<DWORD>(bytes.size()), &count, nullptr)) {
        return std::unexpected(native_error());
    }
    return count;
}

[[nodiscard]] NativeResult<std::size_t> read_block(RawHandle file,
                                                   std::span<std::byte> bytes) noexcept {
    DWORD count = 0;
    if (!ReadFile(file, bytes.data(), static_cast<DWORD>(bytes.size()), &count, nullptr)) {
        return std::unexpected(native_error());
    }
    return count;
}

[[nodiscard]] NativeStatus flush(RawHandle file) noexcept {
    if (!FlushFileBuffers(file)) {
        return std::unexpected(native_error());
    }
    return {};
}

void mark_for_removal(RawHandle file) noexcept {
    FILE_DISPOSITION_INFO disposition{TRUE};
    static_cast<void>(
        SetFileInformationByHandle(file, FileDispositionInfo, &disposition, sizeof(disposition)));
}

void remove_stage_file(RawHandle stage, const Name& name) noexcept {
    auto file = open_relative(stage, name, DELETE, FILE_OPEN, FILE_NON_DIRECTORY_FILE);
    if (file) {
        mark_for_removal(*file);
        static_cast<void>(close_native(*file));
    }
}

void remove_stage(RawHandle /*parent*/, RawHandle stage, const Name& /*name*/) noexcept {
    mark_for_removal(stage);
}

[[nodiscard]] NativeStatus publish(RawHandle parent, RawHandle stage, const Name& /*temporary*/,
                                   const Name& final_name) noexcept {
    struct RenameBuffer {
        DWORD flags;
        HANDLE root;
        DWORD length;
        std::array<wchar_t, Name{}.size()> text;
    };
    static_assert(offsetof(RenameBuffer, text) == offsetof(FILE_RENAME_INFO, FileName));
    RenameBuffer rename_info{
        0, parent,
        static_cast<DWORD>(std::char_traits<char>::length(final_name.data()) * sizeof(wchar_t)),
        wide_name(final_name)};
    // FileRenameInformation supports a directory-handle-relative destination.
    constexpr auto rename_class = static_cast<FILE_INFORMATION_CLASS>(10);
    IO_STATUS_BLOCK info{};
    const auto length = static_cast<ULONG>(offsetof(RenameBuffer, text)) + rename_info.length;
    const auto status = NtSetInformationFile(stage, &info, &rename_info, length, rename_class);
    if (status < 0) {
        return std::unexpected(RtlNtStatusToDosError(status));
    }
    return {};
}

[[nodiscard]] auto process_id() noexcept {
    return GetCurrentProcessId();
}

#else
using RawHandle = int;
constexpr RawHandle invalid_handle = -1;

[[nodiscard]] std::uint32_t native_error() noexcept {
    return static_cast<std::uint32_t>(errno);
}

[[nodiscard]] NativeStatus close_native(RawHandle handle) noexcept {
    // Linux releases the descriptor even when close reports EINTR.
    if (close(handle) != 0) {
        return std::unexpected(native_error());
    }
    return {};
}

[[nodiscard]] bool exists_error(std::uint32_t error) noexcept {
    return error == EEXIST;
}

[[nodiscard]] Code native_code(std::uint32_t error) noexcept {
    if (error == ENOSPC || error == EDQUOT || error == ENOMEM || error == EMFILE ||
        error == ENFILE) {
        return Code::ResourceExceeded;
    }
    if (error == ENOSYS || error == EOPNOTSUPP || error == EXDEV) {
        return Code::Unsupported;
    }
    return Code::InvalidInput;
}

[[nodiscard]] NativeResult<RawHandle> own_fd(int handle) noexcept {
    if (handle < 0) {
        return std::unexpected(native_error());
    }
    return handle;
}

[[nodiscard]] auto open_parent(const std::filesystem::path& parent) noexcept {
    return own_fd(open(parent.c_str(), O_RDONLY | O_DIRECTORY | O_NOFOLLOW | O_CLOEXEC));
}

[[nodiscard]] NativeResult<RawHandle> create_stage(RawHandle parent, const Name& name) noexcept {
    if (mkdirat(parent, name.data(), 0700) != 0) {
        return std::unexpected(native_error());
    }
    auto stage =
        own_fd(openat(parent, name.data(), O_RDONLY | O_DIRECTORY | O_NOFOLLOW | O_CLOEXEC));
    if (!stage) {
        static_cast<void>(unlinkat(parent, name.data(), AT_REMOVEDIR));
    }
    return stage;
}

[[nodiscard]] auto open_output(RawHandle parent, const Name& name) noexcept {
    return own_fd(
        openat(parent, name.data(), O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW | O_CLOEXEC, 0600));
}

[[nodiscard]] auto open_readback(RawHandle parent, const Name& name) noexcept {
    return own_fd(openat(parent, name.data(), O_RDONLY | O_NONBLOCK | O_NOFOLLOW | O_CLOEXEC));
}

struct FileInfo {
    dev_t device;
    ino_t inode;
    std::uint64_t length;
    [[nodiscard]] bool matches(const FileInfo& other) const noexcept {
        return device == other.device && inode == other.inode;
    }
};

[[nodiscard]] NativeResult<FileInfo> file_info(RawHandle file) noexcept {
    struct stat info{};
    if (fstat(file, &info) != 0) {
        return std::unexpected(native_error());
    }
    if (!S_ISREG(info.st_mode) || info.st_size < 0) {
        return std::unexpected(EINVAL);
    }
    return FileInfo{info.st_dev, info.st_ino, static_cast<std::uint64_t>(info.st_size)};
}

[[nodiscard]] NativeResult<std::size_t> write_block(RawHandle file,
                                                    std::span<const std::byte> bytes) noexcept {
    const auto count = write(file, bytes.data(), bytes.size());
    if (count < 0) {
        return std::unexpected(native_error());
    }
    return static_cast<std::size_t>(count);
}

[[nodiscard]] NativeResult<std::size_t> read_block(RawHandle file,
                                                   std::span<std::byte> bytes) noexcept {
    const auto count = read(file, bytes.data(), bytes.size());
    if (count < 0) {
        return std::unexpected(native_error());
    }
    return static_cast<std::size_t>(count);
}

[[nodiscard]] NativeStatus flush(RawHandle file) noexcept {
    if (fsync(file) != 0) {
        return std::unexpected(native_error());
    }
    return {};
}

void remove_stage_file(RawHandle stage, const Name& name) noexcept {
    static_cast<void>(unlinkat(stage, name.data(), 0));
}

[[nodiscard]] bool still_named(RawHandle parent, RawHandle stage, const Name& name) noexcept {
    struct stat pinned{}, named{};
    return fstat(stage, &pinned) == 0 &&
           fstatat(parent, name.data(), &named, AT_SYMLINK_NOFOLLOW) == 0 &&
           S_ISDIR(named.st_mode) && pinned.st_dev == named.st_dev && pinned.st_ino == named.st_ino;
}

void remove_stage(RawHandle parent, RawHandle stage, const Name& name) noexcept {
    if (still_named(parent, stage, name)) {
        static_cast<void>(unlinkat(parent, name.data(), AT_REMOVEDIR));
    }
}

[[nodiscard]] NativeStatus publish(RawHandle parent, RawHandle stage, const Name& temporary,
                                   const Name& final_name) noexcept {
    if (!still_named(parent, stage, temporary)) {
        return std::unexpected(ENOENT);
    }
    if (syscall(SYS_renameat2, parent, temporary.data(), parent, final_name.data(),
                RENAME_NOREPLACE) != 0) {
        return std::unexpected(errno == EINVAL ? static_cast<std::uint32_t>(EOPNOTSUPP)
                                               : native_error());
    }
    return {};
}

[[nodiscard]] auto process_id() noexcept {
    return getpid();
}
#endif

class Handle {
public:
    explicit Handle(RawHandle raw) noexcept : raw_(raw) {}
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    ~Handle() {
        if (raw_ != invalid_handle) {
            static_cast<void>(close_native(raw_));
        }
    }
    [[nodiscard]] RawHandle get() const noexcept {
        return raw_;
    }
    [[nodiscard]] NativeStatus close() noexcept {
        return close_native(std::exchange(raw_, invalid_handle));
    }

private:
    RawHandle raw_;
};

struct Publication {
    const std::stop_token& stop;
    std::optional<CheckpointPublicationDiagnostic>* diagnostic;
    std::optional<std::size_t> chunk;

    [[nodiscard]] auto fail(Stage stage, Code code,
                            std::optional<std::uint32_t> error = {}) const noexcept {
        if (diagnostic != nullptr) {
            diagnostic->emplace(CheckpointPublicationDiagnostic{stage, chunk, error});
        }
        return std::unexpected(core::Failure{code, save});
    }

    [[nodiscard]] auto fail_io(Stage stage, std::uint32_t error) const noexcept {
        return fail(stage, native_code(error), error);
    }
};

struct Staging {
    RawHandle parent;
    Handle directory;
    Name name;
    std::span<const CheckpointPublicationChunk> chunks;
    std::size_t created_chunks = 0;
    bool created_manifest = false;
    bool published = false;

    ~Staging() {
        if (!published) {
            if (created_manifest) {
                remove_stage_file(directory.get(), name_with("", "manifest.json", ""));
            }
            for (const auto& chunk : chunks.first(created_chunks)) {
                remove_stage_file(directory.get(), chunk_name(chunk.logical_name));
            }
            remove_stage(parent, directory.get(), name);
        }
    }
};

[[nodiscard]] NativeResult<RawHandle> fresh_stage(RawHandle parent, Name& name) noexcept {
    static std::atomic<std::uint64_t> counter{0};
    for (int attempt = 0; attempt != 128; ++attempt) {
        name = name_with("", ".uvlas-checkpoint-", "");
        const auto end = name.data() + name.size() - 1;
        auto split = std::to_chars(name.data() + std::char_traits<char>::length(name.data()), end,
                                   process_id())
                         .ptr;
        *split++ = '-';
        static_cast<void>(
            std::to_chars(split, end, counter.fetch_add(1, std::memory_order_relaxed)));
        auto stage = create_stage(parent, name);
        if (stage || !exists_error(stage.error())) {
            return stage;
        }
    }
#ifdef _WIN32
    return std::unexpected(ERROR_ALREADY_EXISTS);
#else
    return std::unexpected(EEXIST);
#endif
}

[[nodiscard]] core::Result<void> readback(RawHandle directory, const Name& name,
                                          const FileInfo& written, SourceDigest expected,
                                          std::span<std::byte> scratch,
                                          const Publication& state) noexcept {
    auto reopened = open_readback(directory, name);
    if (!reopened) {
        return state.fail_io(Stage::ReadBack, reopened.error());
    }
    Handle reader{*reopened};
    const auto before = file_info(reader.get());
    if (!before) {
        return state.fail_io(Stage::ReadBack, before.error());
    }
    if (!before->matches(written) || before->length != written.length) {
        return state.fail(Stage::Verify, Code::CheckpointCorrupt);
    }
    boost::hash2::sha2_256 hash;
    auto remaining = written.length;
    for (;;) {
        if (state.stop.stop_requested()) {
            return state.fail(Stage::ReadBack, Code::Cancelled);
        }
        const auto count = read_block(reader.get(), scratch);
        if (!count) {
            if (interrupted(count.error())) {
                continue;
            }
            return state.fail_io(Stage::ReadBack, count.error());
        }
        if (*count == 0) {
            break;
        }
        if (*count > remaining) {
            return state.fail(Stage::Verify, Code::CheckpointCorrupt);
        }
        hash.update(scratch.data(), *count);
        remaining -= *count;
    }
    const auto after = file_info(reader.get());
    if (!after) {
        return state.fail_io(Stage::ReadBack, after.error());
    }
    if (const auto closed = reader.close(); !closed) {
        return state.fail_io(Stage::Close, closed.error());
    }
    if (remaining != 0 || after->length != written.length ||
        finish_source_digest(hash) != expected) {
        return state.fail(Stage::Verify, Code::CheckpointCorrupt);
    }
    return {};
}

[[nodiscard]] core::Result<void> store(RawHandle directory, const Name& name,
                                       std::span<const std::byte> bytes,
                                       std::optional<SourceDigest> expected,
                                       std::span<std::byte> scratch, const Publication& state,
                                       bool& created) noexcept {
    if (state.stop.stop_requested()) {
        return state.fail(Stage::Create, Code::Cancelled);
    }
    auto output = open_output(directory, name);
    if (!output) {
        return state.fail_io(Stage::Create, output.error());
    }
    created = true;
    Handle file{*output};
    const auto length = bytes.size();
    boost::hash2::sha2_256 source;
    while (!bytes.empty()) {
        if (state.stop.stop_requested()) {
            return state.fail(Stage::Write, Code::Cancelled);
        }
        const auto count = write_block(file.get(), bytes.first(std::min(block_size, bytes.size())));
        if (!count) {
            if (interrupted(count.error())) {
                continue;
            }
            return state.fail_io(Stage::Write, count.error());
        }
        if (*count == 0) {
            return state.fail(Stage::Write, Code::InvalidInput);
        }
        if (!expected) {
            source.update(bytes.data(), *count);
        }
        bytes = bytes.subspan(*count);
    }
    if (const auto flushed = flush(file.get()); !flushed) {
        return state.fail_io(Stage::Flush, flushed.error());
    }
    const auto info = file_info(file.get());
    if (!info) {
        return state.fail_io(Stage::Verify, info.error());
    }
    if (const auto closed = file.close(); !closed) {
        return state.fail_io(Stage::Close, closed.error());
    }
    if (info->length != length) {
        return state.fail(Stage::Verify, Code::CheckpointCorrupt);
    }
    if (!expected) {
        expected = finish_source_digest(source);
    }
    return readback(directory, name, *info, *expected, scratch, state);
}

} // namespace

core::Result<void> publish_checkpoint_directory(
    const std::filesystem::path& parent, std::uint64_t checkpoint_number,
    std::span<const CheckpointPublicationChunk> chunks, std::span<const std::byte> manifest,
    CheckpointPublicationLimits limits, core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<CheckpointPublicationDiagnostic>* diagnostic) noexcept {
    if (diagnostic != nullptr) {
        diagnostic->reset();
    }
    Publication state{stop, diagnostic, {}};
    if (stop.stop_requested()) {
        return state.fail(Stage::Preflight, Code::Cancelled);
    }
    if (parent.empty() || manifest.empty() ||
        parent.native().find(std::filesystem::path::value_type{}) !=
            std::filesystem::path::string_type::npos) {
        return state.fail(Stage::Preflight, Code::InvalidInput);
    }
    if (chunks.size() > limits.maximum_chunks || !std::in_range<core::Index>(chunks.size()) ||
        manifest.size() > limits.maximum_file_bytes ||
        manifest.size() > limits.maximum_total_bytes) {
        return state.fail(Stage::Preflight, Code::ResourceExceeded);
    }
    std::uint64_t total = manifest.size();
    std::size_t largest = manifest.size();
    for (std::size_t i = 0; i < chunks.size(); ++i) {
        state.chunk = i;
        if (stop.stop_requested()) {
            return state.fail(Stage::Preflight, Code::Cancelled);
        }
        const auto& [name, input] = chunks[i];
        if (!safe_chunk_name(name) || input.bytes.size() != input.header.byte_count) {
            return state.fail(Stage::Preflight, Code::InvalidInput);
        }
        if (input.bytes.size() > limits.maximum_file_bytes ||
            input.bytes.size() > limits.maximum_total_bytes - total) {
            return state.fail(Stage::Preflight, Code::ResourceExceeded);
        }
        total += input.bytes.size();
        largest = std::max(largest, input.bytes.size());
    }
    state.chunk.reset();
    auto names = plan.allocate<std::string_view>(static_cast<core::Index>(chunks.size()), save);
    if (!names) {
        return std::unexpected(names.error());
    }
    std::ranges::transform(chunks, names->values().begin(),
                           [](const auto& chunk) { return chunk.logical_name; });
    std::ranges::sort(names->values());
    if (std::ranges::adjacent_find(names->values()) != names->values().end()) {
        return state.fail(Stage::Preflight, Code::InvalidInput);
    }
    auto scratch =
        plan.allocate<std::byte>(static_cast<core::Index>(std::min(largest, block_size)), save);
    if (!scratch) {
        return std::unexpected(scratch.error());
    }
    if (stop.stop_requested()) {
        return state.fail(Stage::OpenParent, Code::Cancelled);
    }
    const auto raw_parent = open_parent(parent);
    if (!raw_parent) {
        return state.fail_io(Stage::OpenParent, raw_parent.error());
    }
    const Handle root{*raw_parent};
    Name temporary{};
    const auto raw_stage = fresh_stage(root.get(), temporary);
    if (!raw_stage) {
        return state.fail_io(Stage::StageDirectory, raw_stage.error());
    }
    Staging stage{root.get(), Handle{*raw_stage}, temporary, chunks};
    for (std::size_t i = 0; i < chunks.size(); ++i) {
        state.chunk = i;
        const auto& [name, input] = chunks[i];
        bool created = false;
        const auto stored = store(stage.directory.get(), chunk_name(name), input.bytes,
                                  input.header.digest, scratch->values(), state, created);
        if (created) {
            ++stage.created_chunks;
        }
        if (!stored) {
            return stored;
        }
    }
    state.chunk.reset();
    const auto stored = store(stage.directory.get(), name_with("", "manifest.json", ""), manifest,
                              {}, scratch->values(), state, stage.created_manifest);
    if (!stored) {
        return stored;
    }
#ifndef _WIN32
    if (const auto flushed = flush(stage.directory.get()); !flushed) {
        return state.fail_io(Stage::Flush, flushed.error());
    }
#endif
    if (stop.stop_requested()) {
        return state.fail(Stage::Publish, Code::Cancelled);
    }
    const auto final_name = make_checkpoint_directory_name(checkpoint_number);
    const auto published =
        publish(root.get(), stage.directory.get(), temporary, name_with("", final_name.data(), ""));
    if (!published) {
        return state.fail_io(exists_error(published.error()) ? Stage::DestinationExists
                                                             : Stage::Publish,
                             published.error());
    }
    stage.published = true;
    return {};
}

} // namespace uvlas::io::detail
