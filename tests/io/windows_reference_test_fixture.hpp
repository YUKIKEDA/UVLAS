#pragma once

#include "windows_reference_root.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <windows.h>
#include <winioctl.h>

namespace uvlas::test::windows_reference_fixture {
namespace fs = std::filesystem;
using io::detail::WindowsHandle;

inline void create_junction(const fs::path& junction, const fs::path& target) {
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

class Junction {
public:
    Junction(const fs::path& link, const fs::path& target) : link_(link) {
        create_junction(link_, target);
    }
    Junction(const Junction&) = delete;
    Junction& operator=(const Junction&) = delete;
    ~Junction() {
        static_cast<void>(RemoveDirectoryW(link_.c_str()));
    }

private:
    fs::path link_;
};

} // namespace uvlas::test::windows_reference_fixture
