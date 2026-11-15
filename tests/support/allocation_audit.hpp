#pragma once

#include <Eigen/Core>

#include <cstdint>

namespace uvlas::test_support {

enum class AllocationBackend : std::uint8_t { Cpp = 1, Eigen = 2, CHeap = 4, External = 8 };
enum class AuditThreadDomain : std::uint8_t { CallingThread, AllThreads };
enum class AllocationAuditStatus : std::uint8_t { Clean, AllocationObserved, Unsupported };
struct AllocationAuditResult {
    AllocationAuditStatus status;
    std::uint64_t cpp_allocation_attempts;
};

/// Audit declared allocation paths in a dedicated test process, never in a shipped binary.
/// C++ allocation attempts are process-wide; Eigen permission is per calling thread.
/// A caller must declare every backend used by the tested code. Unsupported backends
/// or an all-thread Eigen guarantee produce Unsupported, never a false zero-allocation pass.
/// Construct, finish and destroy a guard on the same thread. Each Eigen worker needs
/// its own guard; the C++ counter also observes allocations made by other workers.
class AllocationAudit {
public:
    /// Query coverage before executing a kernel; unsupported guards have no side effects.
    [[nodiscard]] static bool supports(unsigned required_backends,
                                       AuditThreadDomain domain) noexcept;
    explicit AllocationAudit(unsigned required_backends, AuditThreadDomain domain) noexcept;
    ~AllocationAudit() noexcept;
    AllocationAudit(const AllocationAudit&) = delete;
    AllocationAudit& operator=(const AllocationAudit&) = delete;
    [[nodiscard]] AllocationAuditResult finish() noexcept;

private:
    std::uint64_t start_;
    bool supported_;
    bool eigen_guard_;
    bool active_ = true;
    AllocationAuditResult result_{AllocationAuditStatus::Unsupported, 0};
};

constexpr unsigned cpp_backend = static_cast<unsigned>(AllocationBackend::Cpp);
constexpr unsigned eigen_backend = static_cast<unsigned>(AllocationBackend::Eigen);
constexpr unsigned strict_kernel_backends = cpp_backend | eigen_backend;

} // namespace uvlas::test_support
