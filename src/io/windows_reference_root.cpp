#ifdef _WIN32

#include "windows_reference_root.hpp"

#include "windows_checkpoint_query.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <ranges>
#include <utility>
#include <windows.h>
#include <winternl.h>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
constexpr DWORD read_share = FILE_SHARE_READ | FILE_SHARE_WRITE;
constexpr ACCESS_MASK directory_access =
    FILE_LIST_DIRECTORY | FILE_TRAVERSE | FILE_READ_ATTRIBUTES | SYNCHRONIZE;
constexpr ACCESS_MASK file_access = FILE_READ_DATA | FILE_READ_ATTRIBUTES | SYNCHRONIZE;
constexpr auto max_name_units = std::numeric_limits<USHORT>::max() / sizeof(wchar_t);
using Reason = WindowsReferenceOpenReason;

[[nodiscard]] auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}

[[nodiscard]] std::optional<Reason> wrong_kind(HANDLE handle, bool directory) noexcept {
    FILE_ATTRIBUTE_TAG_INFO tags{};
    if (!GetFileInformationByHandleEx(handle, FileAttributeTagInfo, &tags, sizeof(tags))) {
        return Reason::Attributes;
    }
    if ((tags.FileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0) {
        return Reason::Reparse;
    }
    if (((tags.FileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0) != directory ||
        GetFileType(handle) != FILE_TYPE_DISK) {
        return Reason::WrongType;
    }
    return std::nullopt;
}

} // namespace

void WindowsHandleCloser::operator()(void* handle) const noexcept {
    static_cast<void>(CloseHandle(handle));
}

WindowsReferenceRoot::WindowsReferenceRoot(WindowsHandle directory,
                                           core::PlannedBuffer<WindowsHandle> parents) noexcept
    : parents_(std::move(parents)), directory_(std::move(directory)) {}

core::Result<WindowsReferenceRoot>
WindowsReferenceRoot::open(const std::filesystem::path& directory,
                           const std::stop_token& stop) noexcept {
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    if (directory.empty() || directory.native().find(L'\0') != std::wstring::npos) {
        return fail(core::FailureCode::InvalidInput);
    }
    const auto raw =
        CreateFileW(directory.c_str(), directory_access, read_share, nullptr, OPEN_EXISTING,
                    FILE_FLAG_BACKUP_SEMANTICS | FILE_FLAG_OPEN_REPARSE_POINT, nullptr);
    WindowsHandle held{raw == INVALID_HANDLE_VALUE ? nullptr : raw};
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    if (!held) {
        return fail(core::FailureCode::InvalidInput);
    }
    const auto reason = wrong_kind(held.get(), true);
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    if (reason) {
        return fail(core::FailureCode::InvalidInput);
    }
    DWORD volume_flags = 0;
    const bool queried = GetVolumeInformationByHandleW(held.get(), nullptr, 0, nullptr, nullptr,
                                                       &volume_flags, nullptr, 0) != 0;
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    if (!queried || (volume_flags & FILE_SUPPORTS_OPEN_BY_FILE_ID) == 0) {
        return fail(core::FailureCode::Unsupported);
    }
    return WindowsReferenceRoot{std::move(held)};
}

core::Result<WindowsReferenceRoot::RelativeEntry> WindowsReferenceRoot::open_relative(
    const NormalizedReferencePath& path, bool directory, core::ResourcePlan& plan,
    const std::stop_token& stop,
    std::optional<WindowsReferenceOpenDiagnostic>* diagnostic) const noexcept {
    if (diagnostic != nullptr) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    const auto relative = path.text();
    if (!directory_ || relative.empty()) {
        return fail(core::FailureCode::InvalidInput);
    }
    auto components = relative | std::views::split('/');
    core::Index parent_count = -1;
    std::size_t longest = 0;
    for (const auto component : components) {
        ++parent_count;
        longest = std::max(longest, static_cast<std::size_t>(std::ranges::distance(component)));
    }
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    const auto parent_bytes =
        core::checked_byte_count(parent_count, sizeof(WindowsHandle), operation);
    if (!parent_bytes) {
        return std::unexpected(parent_bytes.error());
    }
    const auto name_capacity = std::min(longest, max_name_units);
    const auto name_bytes = static_cast<std::uint64_t>(name_capacity * sizeof(wchar_t));
    if (*parent_bytes > std::numeric_limits<std::uint64_t>::max() - name_bytes) {
        return fail(core::FailureCode::ResourceExceeded);
    }
    const auto required = *parent_bytes + name_bytes;
    const auto remaining = plan.remaining_bytes();
    if (required > remaining) {
        core::FailureDetails details{};
        const auto charged = plan.byte_limit() - remaining;
        if (required <= std::numeric_limits<std::uint64_t>::max() - charged) {
            details.resources =
                core::ResourceUsage{core::ResourceKind::Bytes, charged + required, charged};
        }
        return std::unexpected(
            core::Failure{core::FailureCode::ResourceExceeded, operation, details});
    }
    auto parents = plan.allocate<WindowsHandle>(parent_count, operation);
    if (!parents) {
        return std::unexpected(parents.error());
    }
    auto name = plan.allocate<wchar_t>(static_cast<core::Index>(name_capacity), operation);
    if (!name) {
        return std::unexpected(name.error());
    }

    const auto fail_component = [&](Reason reason, std::size_t offset,
                                    core::FailureCode code =
                                        core::FailureCode::InvalidInput) noexcept {
        if (diagnostic != nullptr) {
            diagnostic->emplace(WindowsReferenceOpenDiagnostic{reason, offset});
        }
        return fail(code);
    };
    HANDLE parent = directory_.get();
    std::size_t index = 0;
    for (const auto part : components) {
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled);
        }
        const std::string_view component{part.begin(), part.end()};
        const auto offset = static_cast<std::size_t>(part.begin() - relative.begin());
        if (component.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            return fail_component(Reason::Encoding, offset);
        }
        const auto units = MultiByteToWideChar(
            CP_UTF8, MB_ERR_INVALID_CHARS, component.data(), static_cast<int>(component.size()),
            name->values().data(), static_cast<int>(name->size()));
        if (units == 0) {
            return fail_component(Reason::Encoding, offset);
        }

        // One component and one held directory per lookup. The exclusion of FILE_SHARE_DELETE
        // applies to each directory for the lifetime of the returned entry.
        const bool is_parent = index < parents->size();
        const bool expects_directory = is_parent || directory;
        const auto bytes = static_cast<USHORT>(static_cast<std::size_t>(units) * sizeof(wchar_t));
        UNICODE_STRING native{bytes, bytes, name->values().data()};
        OBJECT_ATTRIBUTES attributes{sizeof(OBJECT_ATTRIBUTES),
                                     parent,
                                     &native,
                                     OBJ_CASE_INSENSITIVE | OBJ_DONT_REPARSE,
                                     nullptr,
                                     nullptr};
        IO_STATUS_BLOCK status_block{};
        HANDLE raw = nullptr;
        const auto status =
            NtCreateFile(&raw, expects_directory ? directory_access : file_access, &attributes,
                         &status_block, nullptr, 0, read_share, FILE_OPEN,
                         FILE_OPEN_REPARSE_POINT | FILE_SYNCHRONOUS_IO_NONALERT, nullptr, 0);
        WindowsHandle entry{status >= 0 ? raw : nullptr};
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled);
        }
        if (!entry) {
            const auto error = RtlNtStatusToDosError(status);
            if (error == ERROR_REPARSE_POINT_ENCOUNTERED) {
                return fail_component(Reason::Reparse, offset);
            }
            if (error == ERROR_INVALID_FUNCTION || error == ERROR_NOT_SUPPORTED ||
                error == ERROR_CALL_NOT_IMPLEMENTED) {
                return fail_component(Reason::Open, offset, core::FailureCode::Unsupported);
            }
            if (error == ERROR_TOO_MANY_OPEN_FILES ||
                windows_checkpoint_query_failure(error) == core::FailureCode::ResourceExceeded) {
                return fail_component(Reason::Open, offset, core::FailureCode::ResourceExceeded);
            }
            return fail_component(Reason::Open, offset);
        }
        const auto reason = wrong_kind(entry.get(), expects_directory);
        if (stop.stop_requested()) {
            return fail(core::FailureCode::Cancelled);
        }
        if (reason) {
            return fail_component(*reason, offset);
        }
        if (!is_parent) {
            return RelativeEntry{std::move(*parents), std::move(entry)};
        }
        parent = entry.get();
        parents->values()[index] = std::move(entry);
        ++index;
    }
    return fail(core::FailureCode::InvalidInput);
}

core::Result<WindowsReferenceFile> WindowsReferenceRoot::open_file(
    const NormalizedReferencePath& path, core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<WindowsReferenceOpenDiagnostic>* diagnostic) const noexcept {
    auto entry = open_relative(path, false, plan, stop, diagnostic);
    if (!entry) {
        return std::unexpected(entry.error());
    }
    return WindowsReferenceFile{std::move(entry->parents), std::move(entry->handle)};
}

core::Result<WindowsReferenceRoot> WindowsReferenceRoot::open_directory(
    const NormalizedReferencePath& path, core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<WindowsReferenceOpenDiagnostic>* diagnostic) const noexcept {
    auto entry = open_relative(path, true, plan, stop, diagnostic);
    if (!entry) {
        return std::unexpected(entry.error());
    }
    return WindowsReferenceRoot{std::move(entry->handle), std::move(entry->parents)};
}

WindowsReferenceFile::WindowsReferenceFile(core::PlannedBuffer<WindowsHandle> parents,
                                           WindowsHandle file) noexcept
    : parents_(std::move(parents)), file_(std::move(file)) {}

core::Result<ReferenceFileInfo> WindowsReferenceFile::info() noexcept {
    if (!file_) {
        return fail(core::FailureCode::InvalidInput);
    }
    FILE_ID_INFO id{};
    if (!GetFileInformationByHandleEx(file_.get(), FileIdInfo, &id, sizeof(id)) ||
        std::ranges::all_of(id.FileId.Identifier, [](BYTE value) { return value == 0; })) {
        return fail(core::FailureCode::Unsupported);
    }
    FILE_STANDARD_INFO standard{};
    FILE_BASIC_INFO basic{};
    if (!GetFileInformationByHandleEx(file_.get(), FileStandardInfo, &standard, sizeof(standard)) ||
        !GetFileInformationByHandleEx(file_.get(), FileBasicInfo, &basic, sizeof(basic)) ||
        standard.EndOfFile.QuadPart < 0) {
        return fail(core::FailureCode::InvalidInput);
    }
    const ReferenceFileId identity{id.VolumeSerialNumber,
                                   std::bit_cast<std::array<std::byte, 16>>(id.FileId.Identifier)};
    // Access time can change because of our own reads. Only content and metadata indicators
    // are used for before/after comparison.
    return ReferenceFileInfo{identity,
                             static_cast<std::uint64_t>(standard.EndOfFile.QuadPart),
                             {static_cast<std::uint64_t>(basic.LastWriteTime.QuadPart),
                              static_cast<std::uint64_t>(basic.ChangeTime.QuadPart),
                              static_cast<std::uint64_t>(basic.CreationTime.QuadPart),
                              basic.FileAttributes},
                             !standard.Directory && !standard.DeletePending &&
                                 (basic.FileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) == 0};
}

core::Result<std::size_t> WindowsReferenceFile::read(std::span<std::byte> out) noexcept {
    if (!file_ || out.size() > std::numeric_limits<DWORD>::max()) {
        return fail(core::FailureCode::InvalidInput);
    }
    DWORD received = 0;
    if (!ReadFile(file_.get(), out.data(), static_cast<DWORD>(out.size()), &received, nullptr)) {
        if (GetLastError() == ERROR_HANDLE_EOF) {
            return std::size_t{0};
        }
        return fail(core::FailureCode::InvalidInput);
    }
    return static_cast<std::size_t>(received);
}

} // namespace uvlas::io::detail

#endif
