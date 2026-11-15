#include "allocation_audit.hpp"

#include <atomic>
#include <cstdlib>
#include <limits>
#include <new>
#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
constinit std::atomic<std::uint64_t> attempts{0};
constinit std::atomic<unsigned> active_scopes{0};
constinit std::atomic<bool> count_overflow{false};
thread_local unsigned eigen_scopes = 0;
thread_local bool original_eigen_allowed = true;

void count_attempt() noexcept {
    if (active_scopes.load(std::memory_order_relaxed) != 0) {
        if (attempts.fetch_add(1, std::memory_order_relaxed) ==
            std::numeric_limits<std::uint64_t>::max()) {
            count_overflow.store(true, std::memory_order_relaxed);
        }
    }
}

void* allocate(std::size_t bytes) {
    count_attempt();
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    if (void* pointer = std::malloc(bytes == 0 ? 1 : bytes)) {
        return pointer;
    }
    throw std::bad_alloc{};
}

void* allocate_aligned(std::size_t bytes, std::size_t alignment) {
    count_attempt();
    bytes = bytes == 0 ? 1 : bytes;
#if defined(_WIN32)
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* pointer = _aligned_malloc(bytes, alignment);
#else
    if (bytes > std::numeric_limits<std::size_t>::max() - (alignment - 1)) {
        throw std::bad_alloc{};
    }
    const auto rounded = ((bytes + alignment - 1) / alignment) * alignment;
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* pointer = std::aligned_alloc(alignment, rounded);
#endif
    if (pointer != nullptr) {
        return pointer;
    }
    throw std::bad_alloc{};
}

void free_aligned(void* pointer) noexcept {
#if defined(_WIN32)
    _aligned_free(pointer);
#else
    std::free(pointer);
#endif
}
} // namespace

// Test-process replacement only: cover scalar/array, aligned, sized and nothrow forms.
void* operator new(std::size_t bytes) {
    return allocate(bytes);
}
void* operator new[](std::size_t bytes) {
    return allocate(bytes);
}
void operator delete(void* pointer) noexcept {
    std::free(pointer);
}
void operator delete[](void* pointer) noexcept {
    std::free(pointer);
}
void operator delete(void* pointer, std::size_t) noexcept {
    std::free(pointer);
}
void operator delete[](void* pointer, std::size_t) noexcept {
    std::free(pointer);
}
void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept {
    try {
        return allocate(bytes);
    } catch (...) {
        return nullptr;
    }
}
void* operator new[](std::size_t bytes, const std::nothrow_t&) noexcept {
    try {
        return allocate(bytes);
    } catch (...) {
        return nullptr;
    }
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
    std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
    std::free(pointer);
}
void* operator new(std::size_t bytes, std::align_val_t alignment) {
    return allocate_aligned(bytes, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t bytes, std::align_val_t alignment) {
    return allocate_aligned(bytes, static_cast<std::size_t>(alignment));
}
void operator delete(void* pointer, std::align_val_t) noexcept {
    free_aligned(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
    free_aligned(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
    free_aligned(pointer);
}
void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept {
    free_aligned(pointer);
}
void* operator new(std::size_t bytes, std::align_val_t alignment, const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(bytes, static_cast<std::size_t>(alignment));
    } catch (...) {
        return nullptr;
    }
}
void* operator new[](std::size_t bytes, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(bytes, static_cast<std::size_t>(alignment));
    } catch (...) {
        return nullptr;
    }
}
void operator delete(void* pointer, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(pointer);
}
void operator delete[](void* pointer, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(pointer);
}

namespace uvlas::test_support {
bool AllocationAudit::supports(unsigned required, AuditThreadDomain domain) noexcept {
    if (required == 0 || (required & ~strict_kernel_backends) != 0) {
        return false;
    }
    if ((required & eigen_backend) == 0) {
        return true;
    }
#if defined(EIGEN_RUNTIME_NO_MALLOC) && !defined(EIGEN_NO_DEBUG)
    return domain == AuditThreadDomain::CallingThread;
#else
    (void)domain;
    return false;
#endif
}

AllocationAudit::AllocationAudit(unsigned required, AuditThreadDomain domain) noexcept
    : start_(0),
      supported_(supports(required, domain)),
      eigen_guard_(supported_ && (required & eigen_backend) != 0) {
    if (!supported_) {
        active_ = false;
        return;
    }
#if defined(EIGEN_RUNTIME_NO_MALLOC) && !defined(EIGEN_NO_DEBUG)
    if (eigen_guard_) {
        if (eigen_scopes == 0) {
            original_eigen_allowed = Eigen::internal::is_malloc_allowed();
        }
        ++eigen_scopes;
        Eigen::internal::set_is_malloc_allowed(false);
    }
#endif
    active_scopes.fetch_add(1, std::memory_order_relaxed);
    start_ = attempts.load(std::memory_order_relaxed);
}

AllocationAudit::~AllocationAudit() noexcept {
    (void)finish();
}

AllocationAuditResult AllocationAudit::finish() noexcept {
    if (!active_) {
        return result_;
    }
    const auto count = attempts.load(std::memory_order_relaxed) - start_;
    active_scopes.fetch_sub(1, std::memory_order_relaxed);
#if defined(EIGEN_RUNTIME_NO_MALLOC) && !defined(EIGEN_NO_DEBUG)
    if (eigen_guard_ && --eigen_scopes == 0) {
        Eigen::internal::set_is_malloc_allowed(original_eigen_allowed);
    }
#endif
    active_ = false;
    result_ = {supported_ && !count_overflow.load(std::memory_order_relaxed)
                   ? (count == 0 ? AllocationAuditStatus::Clean
                                 : AllocationAuditStatus::AllocationObserved)
                   : AllocationAuditStatus::Unsupported,
               count};
    return result_;
}
} // namespace uvlas::test_support
