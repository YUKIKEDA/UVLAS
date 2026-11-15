#include "../support/allocation_audit.hpp"
#include "../support/thread_determinism.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <memory>
#include <new>
#include <thread>
#include <uvlas/core/instrumentation.hpp>

namespace {
using namespace uvlas::test_support;
}

TEST_CASE("Audit catches all C++ allocation forms", "[negative][N4-09]") {
    for (const int kind : {0, 1, 2, 3, 4, 5, 6, 7}) {
        AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
        // Direct operator calls cannot be elided like a new-expression by the optimizer.
        void* pointer = nullptr;
        switch (kind) {
        case 0:
            pointer = ::operator new(16);
            ::operator delete(pointer);
            break;
        case 1:
            pointer = ::operator new[](16);
            ::operator delete[](pointer);
            break;
        case 2:
            pointer = ::operator new(16, std::align_val_t{64});
            ::operator delete(pointer, std::align_val_t{64});
            break;
        case 3:
            pointer = ::operator new[](16, std::align_val_t{64});
            ::operator delete[](pointer, std::align_val_t{64});
            break;
        case 4:
            pointer = ::operator new(16, std::nothrow);
            ::operator delete(pointer);
            break;
        case 5:
            pointer = ::operator new[](16, std::nothrow);
            ::operator delete[](pointer);
            break;
        case 6:
            pointer = ::operator new(16, std::align_val_t{64}, std::nothrow);
            ::operator delete(pointer, std::align_val_t{64});
            break;
        default:
            pointer = ::operator new[](16, std::align_val_t{64}, std::nothrow);
            ::operator delete[](pointer, std::align_val_t{64});
            break;
        }
        const auto result = audit.finish();
        CHECK(result.status == AllocationAuditStatus::AllocationObserved);
        CHECK(result.cpp_allocation_attempts > 0);
    }
}

TEST_CASE("Audit rejects unsupported coverage and sees workers", "[negative][N4-09]") {
    for (const unsigned requested :
         {0U, static_cast<unsigned>(AllocationBackend::CHeap),
          static_cast<unsigned>(AllocationBackend::External),
          strict_kernel_backends | static_cast<unsigned>(AllocationBackend::CHeap)}) {
        AllocationAudit audit{requested, AuditThreadDomain::CallingThread};
        CHECK(audit.finish().status == AllocationAuditStatus::Unsupported);
    }
    AllocationAudit unsupported_threads{strict_kernel_backends, AuditThreadDomain::AllThreads};
    CHECK_FALSE(AllocationAudit::supports(strict_kernel_backends, AuditThreadDomain::AllThreads));
    CHECK(Eigen::internal::is_malloc_allowed());
    Eigen::VectorXd unguarded = Eigen::VectorXd::Ones(4);
    CHECK(unguarded[0] == 1);
    CHECK(unsupported_threads.finish().status == AllocationAuditStatus::Unsupported);
    std::atomic<bool> start{false};
    std::atomic<bool> done{false};
    std::jthread worker{[&] {
        while (!start.load()) {
            std::this_thread::yield();
        }
        void* data = ::operator new(8);
        ::operator delete(data);
        done.store(true);
    }};
    AllocationAudit workers{cpp_backend, AuditThreadDomain::AllThreads};
    start.store(true);
    while (!done.load()) {
        std::this_thread::yield();
    }
    const auto observed = workers.finish();
    worker.join();
    CHECK(observed.status == AllocationAuditStatus::AllocationObserved);
}

TEST_CASE("Nested audit guards retain Eigen prohibition", "[negative][N4-09]") {
    AllocationAudit outer{strict_kernel_backends, AuditThreadDomain::CallingThread};
    AllocationAudit inner{strict_kernel_backends, AuditThreadDomain::CallingThread};
    const auto outer_result = outer.finish();
    const bool still_forbidden = !Eigen::internal::is_malloc_allowed();
    const auto inner_result = inner.finish();
    CHECK(outer_result.status == AllocationAuditStatus::Clean);
    CHECK(inner_result.status == AllocationAuditStatus::Clean);
    CHECK(still_forbidden);
    CHECK(Eigen::internal::is_malloc_allowed());
}

TEST_CASE("Prepared numerical scopes allocate no C++ or Eigen storage", "[property][N4-09]") {
    const Eigen::VectorXd input = Eigen::VectorXd::Ones(8);
    Eigen::VectorXd output = Eigen::VectorXd::Zero(8);
    uvlas::core::StageRecorder recorder{true, uvlas::core::Operation::PrepareModel};
    AllocationAudit audit{strict_kernel_backends, AuditThreadDomain::CallingThread};
    {
        const uvlas::core::StageScope scope{recorder, uvlas::core::MeasurementStage::Structure};
        output.noalias() = input * 2;
    }
    const auto result = audit.finish();
    CHECK(result.status == AllocationAuditStatus::Clean);
    CHECK(output[0] == 2);
    CHECK(Eigen::internal::is_malloc_allowed());
}

TEST_CASE("Thread determinism detects changed result bytes", "[property][N4-03]") {
    const std::array<std::size_t, 3> counts{1, 2, 4};
    const auto run = [](std::size_t threads) {
        std::array<double, 32> output{};
        std::vector<std::jthread> workers;
        workers.reserve(threads);
        for (std::size_t worker = 0; worker < threads; ++worker) {
            workers.emplace_back([&, worker] {
                for (std::size_t i = worker; i < output.size(); i += threads) {
                    output[i] = static_cast<double>(i) * 0.5;
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
        const auto bytes = std::as_bytes(std::span{output});
        return std::vector<std::byte>{bytes.begin(), bytes.end()};
    };
    CHECK(identical_across_threads(counts, run));
    CHECK_FALSE(identical_across_threads(counts, [](std::size_t threads) {
        return std::vector<std::byte>{static_cast<std::byte>(threads)};
    }));
    CHECK_FALSE(identical_across_threads(std::span<const std::size_t>{}, run));
    CHECK_FALSE(identical_across_threads(std::array<std::size_t, 1>{1}, run));
    CHECK_FALSE(identical_across_threads(std::array<std::size_t, 2>{1, 1}, run));
    CHECK_FALSE(identical_across_threads(std::array<std::size_t, 2>{1, 0}, run));
}
