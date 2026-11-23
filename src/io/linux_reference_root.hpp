#pragma once

#include "reference_file_snapshot.hpp"
#include "reference_path.hpp"

#include <filesystem>

namespace uvlas::io::detail {

// Own exactly one native Linux descriptor. A move leaves the source closed.
class LinuxReferenceHandle {
public:
    explicit LinuxReferenceHandle(int descriptor) noexcept;
    LinuxReferenceHandle(LinuxReferenceHandle&& other) noexcept;
    LinuxReferenceHandle& operator=(LinuxReferenceHandle&&) = delete;
    LinuxReferenceHandle(const LinuxReferenceHandle&) = delete;
    LinuxReferenceHandle& operator=(const LinuxReferenceHandle&) = delete;
    ~LinuxReferenceHandle();
    [[nodiscard]] int get() const noexcept;

private:
    int descriptor_;
};

enum class LinuxReferenceOpenReason : std::uint8_t {
    Open,
    Link,
    Boundary,
    Attributes,
    WrongType,
};

// Owns the cursor and metadata source for one already-open regular file.
// Both operations reject a moved-from reader.
class LinuxReferenceFile final : public ReferenceFileReader {
public:
    LinuxReferenceFile(LinuxReferenceFile&&) noexcept = default;
    LinuxReferenceFile& operator=(LinuxReferenceFile&&) = delete;
    [[nodiscard]] core::Result<ReferenceFileInfo> info() noexcept override;
    [[nodiscard]] core::Result<std::size_t> read(std::span<std::byte> out) noexcept override;

private:
    friend class LinuxReferenceRoot;
    explicit LinuxReferenceFile(LinuxReferenceHandle file) noexcept;
    LinuxReferenceHandle file_;
};

// Anchors lookups to a caller-selected directory handle. No path-string fallback is used.
class LinuxReferenceRoot {
public:
    [[nodiscard]] static core::Result<LinuxReferenceRoot>
    open(const std::filesystem::path& directory, const std::stop_token& stop) noexcept;
    [[nodiscard]] bool is_open() const noexcept {
        return directory_.get() >= 0;
    }

    // The plan pays for the terminated native path. The kernel resolves the whole path at once.
    // A diagnostic gives the native rejection category when one is observed; Linux does not
    // report the offending component. Success and cancellation clear the diagnostic.
    [[nodiscard]] core::Result<LinuxReferenceFile>
    open_file(const NormalizedReferencePath& path, core::ResourcePlan& plan,
              const std::stop_token& stop,
              std::optional<LinuxReferenceOpenReason>* diagnostic = nullptr) const noexcept;

private:
    friend struct CheckpointDirectoryReader;
    explicit LinuxReferenceRoot(LinuxReferenceHandle directory) noexcept;
    LinuxReferenceHandle directory_;
};

} // namespace uvlas::io::detail
