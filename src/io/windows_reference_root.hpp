#pragma once

#include "reference_file_snapshot.hpp"
#include "reference_path.hpp"

#include <filesystem>
#include <memory>

namespace uvlas::io::detail {

struct WindowsHandleCloser {
    void operator()(void* handle) const noexcept;
};

using WindowsHandle = std::unique_ptr<void, WindowsHandleCloser>;

enum class WindowsReferenceOpenReason : std::uint8_t {
    Open,
    Encoding,
    Attributes,
    Reparse,
    WrongType,
};

struct WindowsReferenceOpenDiagnostic {
    WindowsReferenceOpenReason reason;
    std::size_t normalized_offset;
};

// Owns the synchronous file cursor and pins the intermediate directories until destruction.
// Keep the originating root alive during the operation. A moved-from reader fails with
// InvalidInput.
class WindowsReferenceFile final : public ReferenceFileReader {
public:
    WindowsReferenceFile(WindowsReferenceFile&&) noexcept = default;
    WindowsReferenceFile& operator=(WindowsReferenceFile&&) = delete;
    [[nodiscard]] core::Result<ReferenceFileInfo> info() noexcept override;
    [[nodiscard]] core::Result<std::size_t> read(std::span<std::byte> out) noexcept override;

private:
    friend class WindowsReferenceRoot;
    WindowsReferenceFile(core::PlannedBuffer<WindowsHandle> parents, WindowsHandle file) noexcept;

    core::PlannedBuffer<WindowsHandle> parents_;
    WindowsHandle file_;
};

// Opens only a caller-selected directory and anchors subsequent lookups to its handle.
class WindowsReferenceRoot {
public:
    [[nodiscard]] static core::Result<WindowsReferenceRoot>
    open(const std::filesystem::path& directory, const std::stop_token& stop) noexcept;
    [[nodiscard]] bool is_open() const noexcept {
        return directory_ != nullptr;
    }

    // Path components and their count come from normalize_reference_path. Charges the UTF-16
    // scratch and directory owner array to plan before opening descendants. A diagnostic, if
    // available, gives the failing component's byte offset in the normalized UTF-8 spelling.
    [[nodiscard]] core::Result<WindowsReferenceFile>
    open_file(const NormalizedReferencePath& path, core::ResourcePlan& plan,
              const std::stop_token& stop,
              std::optional<WindowsReferenceOpenDiagnostic>* diagnostic = nullptr) const noexcept;

private:
    friend struct CheckpointDirectoryReader;
    explicit WindowsReferenceRoot(WindowsHandle directory) noexcept;
    WindowsHandle directory_;
};

} // namespace uvlas::io::detail
