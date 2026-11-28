#ifdef __linux__

#include "linux_reference_root.hpp"

#include <algorithm>
#include <bit>
#include <cerrno>
#include <fcntl.h>
#include <limits>
#include <linux/openat2.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
using Reason = LinuxReferenceOpenReason;

[[nodiscard]] auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}

[[nodiscard]] core::Result<LinuxReferenceHandle>
open_path(int directory, const char* path, const open_how& how, const std::stop_token& stop,
          std::optional<Reason>* diagnostic) noexcept {
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    const auto opened = ::syscall(SYS_openat2, directory, path, &how, sizeof(how));
    const int error = opened < 0 ? errno : 0;
    LinuxReferenceHandle held{static_cast<int>(opened)};
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    if (held.get() >= 0) {
        return held;
    }
    if (diagnostic != nullptr) {
        if (error == ELOOP) {
            *diagnostic = Reason::Link;
        } else if (error == EXDEV || error == EAGAIN) {
            *diagnostic = Reason::Boundary;
        } else {
            *diagnostic = Reason::Open;
        }
    }
    if (error == ENOSYS || error == EINVAL || error == E2BIG || error == EOPNOTSUPP) {
        return fail(core::FailureCode::Unsupported);
    }
    if (error == EMFILE || error == ENFILE || error == ENOMEM) {
        return fail(core::FailureCode::ResourceExceeded);
    }
    return fail(core::FailureCode::InvalidInput);
}

} // namespace

LinuxReferenceHandle::LinuxReferenceHandle(int descriptor) noexcept : descriptor_(descriptor) {}

LinuxReferenceHandle::LinuxReferenceHandle(LinuxReferenceHandle&& other) noexcept
    : descriptor_(std::exchange(other.descriptor_, -1)) {}

LinuxReferenceHandle::~LinuxReferenceHandle() {
    if (descriptor_ >= 0) {
        // Linux releases the descriptor even on EINTR. A retry could close a reused fd.
        static_cast<void>(::close(descriptor_));
    }
}

int LinuxReferenceHandle::get() const noexcept {
    return descriptor_;
}

LinuxReferenceRoot::LinuxReferenceRoot(LinuxReferenceHandle directory) noexcept
    : directory_(std::move(directory)) {}

core::Result<LinuxReferenceRoot> LinuxReferenceRoot::open(const std::filesystem::path& directory,
                                                          const std::stop_token& stop) noexcept {
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    if (directory.empty() || directory.native().find('\0') != std::string::npos) {
        return fail(core::FailureCode::InvalidInput);
    }
    constexpr open_how how{O_PATH | O_DIRECTORY | O_CLOEXEC, 0, RESOLVE_NO_SYMLINKS};
    auto root = open_path(AT_FDCWD, directory.c_str(), how, stop, nullptr);
    if (!root) {
        return std::unexpected(root.error());
    }
    return LinuxReferenceRoot{std::move(*root)};
}

core::Result<LinuxReferenceHandle> LinuxReferenceRoot::open_relative(
    const NormalizedReferencePath& path, bool directory, core::ResourcePlan& plan,
    const std::stop_token& stop,
    std::optional<LinuxReferenceOpenReason>* diagnostic) const noexcept {
    if (diagnostic != nullptr) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    const auto relative = path.text();
    if (directory_.get() < 0 || relative.empty()) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (std::cmp_greater_equal(relative.size(), std::numeric_limits<core::Index>::max())) {
        return fail(core::FailureCode::ResourceExceeded);
    }
    auto path_buffer =
        plan.allocate<char>(static_cast<core::Index>(relative.size()) + 1, operation);
    if (!path_buffer) {
        return std::unexpected(path_buffer.error());
    }
    const auto native_path = path_buffer->values();
    std::ranges::copy(relative, native_path.begin());
    native_path[relative.size()] = '\0';

    // Nonblocking open keeps an unexpected FIFO from waiting for a writer before fstat.
    const open_how how{
        static_cast<std::uint64_t>(directory ? O_PATH | O_DIRECTORY | O_CLOEXEC
                                             : O_RDONLY | O_CLOEXEC | O_NONBLOCK | O_NOCTTY),
        0, RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS};
    return open_path(directory_.get(), native_path.data(), how, stop, diagnostic);
}

core::Result<LinuxReferenceFile>
LinuxReferenceRoot::open_file(const NormalizedReferencePath& path, core::ResourcePlan& plan,
                              const std::stop_token& stop,
                              std::optional<LinuxReferenceOpenReason>* diagnostic) const noexcept {
    auto file = open_relative(path, false, plan, stop, diagnostic);
    if (!file) {
        return std::unexpected(file.error());
    }
    struct stat metadata{};
    const auto observed = ::fstat(file->get(), &metadata);
    if (stop.stop_requested()) {
        return fail(core::FailureCode::Cancelled);
    }
    if (observed != 0 || !S_ISREG(metadata.st_mode)) {
        if (diagnostic != nullptr) {
            *diagnostic = observed == 0 ? Reason::WrongType : Reason::Attributes;
        }
        return fail(core::FailureCode::InvalidInput);
    }
    return LinuxReferenceFile{std::move(*file)};
}

core::Result<LinuxReferenceRoot> LinuxReferenceRoot::open_directory(
    const NormalizedReferencePath& path, core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<LinuxReferenceOpenReason>* diagnostic) const noexcept {
    auto directory = open_relative(path, true, plan, stop, diagnostic);
    if (!directory) {
        return std::unexpected(directory.error());
    }
    return LinuxReferenceRoot{std::move(*directory)};
}

LinuxReferenceFile::LinuxReferenceFile(LinuxReferenceHandle file) noexcept
    : file_(std::move(file)) {}

core::Result<ReferenceFileInfo> LinuxReferenceFile::info() noexcept {
    struct stat metadata{};
    if (file_.get() < 0 || ::fstat(file_.get(), &metadata) != 0 ||
        !std::in_range<std::uint64_t>(metadata.st_size)) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (!std::in_range<std::uint64_t>(metadata.st_dev) ||
        !std::in_range<std::uint64_t>(metadata.st_ino)) {
        return fail(core::FailureCode::Unsupported);
    }
    ReferenceFileId identity{static_cast<std::uint64_t>(metadata.st_dev), {}};
    const auto inode = std::bit_cast<std::array<std::byte, sizeof(std::uint64_t)>>(
        static_cast<std::uint64_t>(metadata.st_ino));
    std::ranges::copy(inode, identity.object.begin());
    // Reading may update atime. Keep both components of the content and status change times.
    return ReferenceFileInfo{identity,
                             static_cast<std::uint64_t>(metadata.st_size),
                             {static_cast<std::uint64_t>(metadata.st_mtim.tv_sec),
                              static_cast<std::uint64_t>(metadata.st_mtim.tv_nsec),
                              static_cast<std::uint64_t>(metadata.st_ctim.tv_sec),
                              static_cast<std::uint64_t>(metadata.st_ctim.tv_nsec)},
                             S_ISREG(metadata.st_mode)};
}

core::Result<std::size_t> LinuxReferenceFile::read(std::span<std::byte> out) noexcept {
    if (file_.get() < 0 || std::cmp_greater(out.size(), std::numeric_limits<ssize_t>::max())) {
        return fail(core::FailureCode::InvalidInput);
    }
    // One syscall per metered adapter call. An interrupted read fails without a hidden retry.
    const auto count = ::read(file_.get(), out.data(), out.size());
    if (count < 0) {
        return fail(core::FailureCode::InvalidInput);
    }
    return static_cast<std::size_t>(count);
}

} // namespace uvlas::io::detail

#endif
