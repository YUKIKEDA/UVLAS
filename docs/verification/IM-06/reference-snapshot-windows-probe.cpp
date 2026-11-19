#define WIN32_LEAN_AND_MEAN
#include <cstdio>
#include <cstring>
#include <cwchar>
#include <windows.h>
#include <winternl.h>

HANDLE open_relative(HANDLE parent, const wchar_t* name, ULONG type) {
    UNICODE_STRING relative{};
    relative.Length = static_cast<USHORT>(std::wcslen(name) * sizeof(wchar_t));
    relative.MaximumLength = relative.Length;
    relative.Buffer = const_cast<wchar_t*>(name);
    OBJECT_ATTRIBUTES attributes{};
    attributes.Length = sizeof(attributes);
    attributes.RootDirectory = parent;
    attributes.ObjectName = &relative;
    attributes.Attributes = OBJ_CASE_INSENSITIVE;
    IO_STATUS_BLOCK status_block{};
    HANDLE child{};
    ULONG access = FILE_READ_ATTRIBUTES | SYNCHRONIZE;
#ifdef DO_LIST
    if (type == FILE_DIRECTORY_FILE)
        access |= FILE_LIST_DIRECTORY;
#endif
    const auto status = NtCreateFile(
        &child, access, &attributes, &status_block, nullptr, 0, FILE_SHARE_READ | FILE_SHARE_WRITE,
        FILE_OPEN, type | FILE_OPEN_REPARSE_POINT | FILE_SYNCHRONOUS_IO_NONALERT, nullptr, 0);
    if (status < 0) {
        std::fprintf(stderr, "relative open failed: 0x%08lX\n", static_cast<unsigned long>(status));
        return INVALID_HANDLE_VALUE;
    }
    return child;
}

bool marked_reparse(HANDLE file, bool expected) {
    FILE_ATTRIBUTE_TAG_INFO info{};
    return GetFileInformationByHandleEx(file, FileAttributeTagInfo, &info, sizeof(info)) &&
           ((info.FileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0) == expected;
}

int wmain(int argc, wchar_t** argv) {
    if (argc != 4)
        return 1;
    const auto root =
        CreateFileW(argv[1], FILE_READ_ATTRIBUTES | FILE_TRAVERSE | SYNCHRONIZE,
                    FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, OPEN_EXISTING,
                    FILE_FLAG_BACKUP_SEMANTICS | FILE_FLAG_OPEN_REPARSE_POINT, nullptr);
    if (root == INVALID_HANDLE_VALUE || !marked_reparse(root, false))
        return 2;
    const auto directory = open_relative(root, L"parent", FILE_DIRECTORY_FILE);
    const auto junction = open_relative(root, L"junction", FILE_DIRECTORY_FILE);
    if (directory == INVALID_HANDLE_VALUE || junction == INVALID_HANDLE_VALUE ||
        !marked_reparse(directory, false) || !marked_reparse(junction, true))
        return 3;
    CloseHandle(junction);
    const auto file = open_relative(directory, L"original.txt", FILE_NON_DIRECTORY_FILE);
    const auto duplicate = open_relative(directory, L"alias.txt", FILE_NON_DIRECTORY_FILE);
    if (file == INVALID_HANDLE_VALUE || duplicate == INVALID_HANDLE_VALUE ||
        !marked_reparse(file, false))
        return 4;
    FILE_ID_INFO first{}, second{};
    if (!GetFileInformationByHandleEx(file, FileIdInfo, &first, sizeof(first)) ||
        !GetFileInformationByHandleEx(duplicate, FileIdInfo, &second, sizeof(second)) ||
        first.VolumeSerialNumber != second.VolumeSerialNumber ||
        std::memcmp(first.FileId.Identifier, second.FileId.Identifier,
                    sizeof(first.FileId.Identifier)) != 0)
        return 5;
    CloseHandle(file);
    CloseHandle(duplicate);
    const auto first_move = MoveFileExW(argv[2], argv[3], 0);
    const auto move_error = first_move ? ERROR_SUCCESS : GetLastError();
    std::printf("first rename: success=%d error=%lu\n", first_move, move_error);
    if (first_move || move_error != ERROR_SHARING_VIOLATION)
        return 6;
    CloseHandle(directory);
    if (!MoveFileExW(argv[2], argv[3], 0))
        return 7;
    CloseHandle(root);
    std::puts("Normal directory opened; junction identified; hard-link ID matched; parent rename "
              "blocked until handle release.");
    return 0;
}
