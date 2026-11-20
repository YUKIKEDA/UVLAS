#ifdef _WIN32

#include "windows_reference_root.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <windows.h>
#include <winioctl.h>

namespace {

namespace fs = std::filesystem;
namespace core = uvlas::core;
using uvlas::io::detail::normalize_reference_path;
using uvlas::io::detail::read_reference_file;
using uvlas::io::detail::WindowsHandle;
using uvlas::io::detail::WindowsReferenceOpenDiagnostic;
using uvlas::io::detail::WindowsReferenceOpenReason;
using uvlas::io::detail::WindowsReferenceRoot;

struct Tree {
    const fs::path directory;
    explicit Tree(std::wstring_view name)
        : directory(fs::temp_directory_path() /
                    (L"uvlas-reference-" + std::to_wstring(GetCurrentProcessId()) + L"-" +
                     std::wstring{name})) {
        REQUIRE(fs::create_directory(directory));
    }
    ~Tree() {
        // Unlink the junction itself before removing the exclusively created fixture tree.
        static_cast<void>(RemoveDirectoryW((directory / L"junction").c_str()));
        std::error_code error;
        fs::remove_all(directory, error);
    }
    Tree(const Tree&) = delete;
    Tree& operator=(const Tree&) = delete;
};

void write_file(const fs::path& path, std::string_view bytes) {
    std::ofstream stream(path, std::ios::binary);
    stream << bytes;
    stream.close();
    REQUIRE(stream.good());
}

void check_delete_sharing(const fs::path& path) {
    const auto raw =
        CreateFileW(path.c_str(), DELETE, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                    nullptr, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, nullptr);
    const auto error = raw == INVALID_HANDLE_VALUE ? GetLastError() : ERROR_SUCCESS;
    const WindowsHandle handle{raw == INVALID_HANDLE_VALUE ? nullptr : raw};
    CHECK(error == ERROR_SHARING_VIOLATION);
}

void create_junction(const fs::path& junction, const fs::path& target) {
    REQUIRE(fs::create_directory(junction));
    const auto raw =
        CreateFileW(junction.c_str(), GENERIC_WRITE, 0, nullptr, OPEN_EXISTING,
                    FILE_FLAG_BACKUP_SEMANTICS | FILE_FLAG_OPEN_REPARSE_POINT, nullptr);
    REQUIRE(raw != INVALID_HANDLE_VALUE);
    const WindowsHandle handle{raw};
    // FSCTL_SET_REPARSE_POINT mount-point wire format. It also works on CI hosts where
    // creating symbolic links requires a privilege the test account does not have.
    struct MountPoint {
        DWORD tag = IO_REPARSE_TAG_MOUNT_POINT;
        WORD data_bytes = 0;
        WORD reserved = 0;
        WORD substitute_offset = 0;
        WORD substitute_bytes = 0;
        WORD print_offset = 0;
        WORD print_bytes = 0;
        std::array<wchar_t, 1024> text{};
    } buffer;
    static_assert(offsetof(MountPoint, text) == 16);
    const auto native = L"\\??\\" + target.native();
    REQUIRE(native.size() + 2 <= buffer.text.size());
    std::ranges::copy(native, buffer.text.begin());
    buffer.substitute_bytes = static_cast<WORD>(native.size() * sizeof(wchar_t));
    buffer.print_offset = static_cast<WORD>((native.size() + 1) * sizeof(wchar_t));
    buffer.data_bytes = static_cast<WORD>(8 + (native.size() + 2) * sizeof(wchar_t));
    DWORD returned = 0;
    REQUIRE(DeviceIoControl(handle.get(), FSCTL_SET_REPARSE_POINT, &buffer, 8U + buffer.data_bytes,
                            nullptr, 0, &returned, nullptr));
}

[[nodiscard]] auto normalized(std::string_view path, core::ResourcePlan& plan) {
    return normalize_reference_path(path, {256, 256, 8}, plan, {}, nullptr);
}

} // namespace

TEST_CASE("Windows reference file holds the verified bytes and the complete file identity",
          "[IO4-01][N4-02][verification][property]") {
    const Tree tree(L"identity");
    REQUIRE(fs::create_directory(tree.directory / L"sub"));
    const auto original = tree.directory / L"sub" / L"\u00e9.bin";
    write_file(original, "abc");
    fs::create_hard_link(original, tree.directory / L"alias.bin");
    write_file(tree.directory / L"copy.bin", "abc");
    // SHA-256 of the three ASCII bytes abc (independent hashlib oracle).
    constexpr uvlas::io::SourceDigest digest{
        {std::byte{0xba}, std::byte{0x78}, std::byte{0x16}, std::byte{0xbf}, std::byte{0x8f},
         std::byte{0x01}, std::byte{0xcf}, std::byte{0xea}, std::byte{0x41}, std::byte{0x41},
         std::byte{0x40}, std::byte{0xde}, std::byte{0x5d}, std::byte{0xae}, std::byte{0x22},
         std::byte{0x23}, std::byte{0xb0}, std::byte{0x03}, std::byte{0x61}, std::byte{0xa3},
         std::byte{0x96}, std::byte{0x17}, std::byte{0x7a}, std::byte{0x9c}, std::byte{0xb4},
         std::byte{0x10}, std::byte{0xff}, std::byte{0x61}, std::byte{0xf2}, std::byte{0x00},
         std::byte{0x15}, std::byte{0xad}}};
    const auto snapshot = [&] {
        auto root = WindowsReferenceRoot::open(tree.directory, {});
        REQUIRE(root.has_value());
        core::ResourcePlan plan{4096};
        const auto path = normalized("sub/./\xc3\xa9.bin", plan);
        const auto alias = normalized("alias.bin", plan);
        const auto copy = normalized("copy.bin", plan);
        REQUIRE(path.has_value());
        REQUIRE(alias.has_value());
        REQUIRE(copy.has_value());
        auto data = root->open_file(*path, plan, {});
        REQUIRE(data.has_value());
        auto link = root->open_file(*alias, plan, {});
        REQUIRE(link.has_value());
        auto duplicate = root->open_file(*copy, plan, {});
        REQUIRE(duplicate.has_value());
        const auto file_info = data->info();
        const auto link_info = link->info();
        const auto duplicate_info = duplicate->info();
        REQUIRE(file_info.has_value());
        REQUIRE(link_info.has_value());
        REQUIRE(duplicate_info.has_value());
        CHECK(file_info->identity == link_info->identity);
        CHECK_FALSE(file_info->identity == duplicate_info->identity);
        CHECK(file_info->byte_count == 3);
        const auto manifest =
            uvlas::io::detail::read_manifest_file(*link, {3, 4, 2, 3}, plan, {}, nullptr);
        REQUIRE(manifest.has_value());
        if (const auto view = manifest->view()) {
            CHECK(view->file_info == *link_info);
            CHECK(view->digest == digest);
        } else {
            FAIL("Manifest has no owned view");
        }
        return read_reference_file(*data, digest, {3, 4, 2, 3}, plan, {}, nullptr);
    }();
    REQUIRE(snapshot.has_value());
    REQUIRE(fs::remove(original));
    constexpr std::string_view contents = "abc";
    const auto expected = std::as_bytes(std::span{contents});
    if (const auto view = snapshot->view()) {
        CHECK(std::ranges::equal(view->bytes, expected));
        CHECK(view->digest == digest);
    } else {
        FAIL("The owning snapshot must expose a view after file removal");
    }
}

TEST_CASE("Windows reference lookup refuses reparse points and pins intermediate directories",
          "[IO4-01][verification][negative]") {
    const Tree tree(L"lookup");
    const auto outer = tree.directory / L"outer";
    const auto inner = outer / L"inner";
    REQUIRE(fs::create_directories(inner));
    write_file(inner / L"data.bin", "abc");
    create_junction(tree.directory / L"junction", inner);
    const auto root = WindowsReferenceRoot::open(tree.directory, {});
    REQUIRE(root.has_value());
    check_delete_sharing(tree.directory);
    const auto rejected_root = WindowsReferenceRoot::open(tree.directory / L"junction", {});
    REQUIRE_FALSE(rejected_root.has_value());
    CHECK(rejected_root.error().code == core::FailureCode::InvalidInput);
    core::ResourcePlan plan{4096};
    for (const std::string_view name : {"junction", "junction/data.bin"}) {
        CAPTURE(name);
        const auto path = normalized(name, plan);
        REQUIRE(path.has_value());
        std::optional<WindowsReferenceOpenDiagnostic> diagnostic;
        const auto rejected = root->open_file(*path, plan, {}, &diagnostic);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error().code == core::FailureCode::InvalidInput);
        if (diagnostic) {
            CHECK(diagnostic->reason == WindowsReferenceOpenReason::Reparse);
            CHECK(diagnostic->normalized_offset == 0);
        } else {
            FAIL("Expected an owned diagnostic for the reparse point");
        }
    }
    const auto path = normalized("outer/inner/data.bin", plan);
    REQUIRE(path.has_value());
    const auto moved = tree.directory / L"moved";
    {
        const auto file = root->open_file(*path, plan, {});
        REQUIRE(file.has_value());
        check_delete_sharing(outer);
        const auto renamed = MoveFileExW(outer.c_str(), moved.c_str(), 0);
        const auto error = renamed ? ERROR_SUCCESS : GetLastError();
        CHECK_FALSE(renamed);
        CHECK(error == ERROR_SHARING_VIOLATION);
    }
    REQUIRE(MoveFileExW(outer.c_str(), moved.c_str(), 0));
}

TEST_CASE("Windows reference lookup handles cancellation, path faults and handle metadata changes",
          "[IO4-01][N4-02][verification][negative]") {
    const Tree tree(L"failures");
    REQUIRE(fs::create_directory(tree.directory / L"sub"));
    const auto filename = tree.directory / L"sub" / L"one.bin";
    write_file(filename, "abc");
    const auto root = WindowsReferenceRoot::open(tree.directory, {});
    REQUIRE(root.has_value());
    core::ResourcePlan plan{4096};
    const auto data = normalized("sub/one.bin", plan);
    REQUIRE(data.has_value());
    std::stop_source source; // NOLINT(misc-const-correctness): MSVC request_stop is non-const.
    REQUIRE(source.request_stop());
    core::ResourcePlan empty{0};
    std::optional<WindowsReferenceOpenDiagnostic> diagnostic{
        WindowsReferenceOpenDiagnostic{WindowsReferenceOpenReason::Open, 0}};
    const auto cancelled = root->open_file(*data, empty, source.get_token(), &diagnostic);
    REQUIRE_FALSE(cancelled.has_value());
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK_FALSE(diagnostic.has_value());
    CHECK(empty.allocated_bytes() == 0);
    const auto exhausted = root->open_file(*data, empty, {});
    REQUIRE_FALSE(exhausted.has_value());
    CHECK(exhausted.error().code == core::FailureCode::ResourceExceeded);
    CHECK(empty.allocated_bytes() == 0);
    const auto workspace_bytes = sizeof(WindowsHandle) + 7 * sizeof(wchar_t); // "one.bin"
    core::ResourcePlan short_plan{workspace_bytes - 1};
    const auto short_workspace = root->open_file(*data, short_plan, {});
    REQUIRE_FALSE(short_workspace.has_value());
    CHECK(short_workspace.error().code == core::FailureCode::ResourceExceeded);
    CHECK(short_plan.allocated_bytes() == 0);
    {
        core::ResourcePlan exact_plan{workspace_bytes};
        const auto exact = root->open_file(*data, exact_plan, {});
        REQUIRE(exact.has_value());
        CHECK(exact_plan.allocated_bytes() == workspace_bytes);
    }
    const auto missing = normalized("sub/missing", plan);
    REQUIRE(missing.has_value());
    const auto rejected = root->open_file(*missing, plan, {}, &diagnostic);
    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    if (diagnostic) {
        CHECK(diagnostic->reason == WindowsReferenceOpenReason::Open);
        CHECK(diagnostic->normalized_offset == 4);
    } else {
        FAIL("Expected an owned diagnostic for the missing path");
    }
    const auto moved = tree.directory / L"unheld";
    REQUIRE(MoveFileExW((tree.directory / L"sub").c_str(), moved.c_str(), 0));
    REQUIRE(MoveFileExW(moved.c_str(), (tree.directory / L"sub").c_str(), 0));

    auto file = root->open_file(*data, plan, {});
    REQUIRE(file.has_value());
    const auto before = file->info();
    REQUIRE(before.has_value());
    const auto initial_time = fs::last_write_time(filename);
    fs::last_write_time(filename, initial_time + std::chrono::seconds{20});
    const auto after = file->info();
    REQUIRE(after.has_value());
    CHECK(after->identity == before->identity);
    CHECK(after->byte_count == before->byte_count);
    CHECK_FALSE(after->change_indicators == before->change_indicators);
}

#endif
