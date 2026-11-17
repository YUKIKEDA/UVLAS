#include "input_document_detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string_view>
#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {
// Process-local injection belongs only to this dedicated test executable.
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
// Replacement new injects only its standard allocation failure. Throwing
// arbitrary exceptions here can violate standard-library noexcept contracts;
// core/failure_test.cpp tests guard_entry's arbitrary-exception guarantee.
void allocation_attempt() {
    const auto next = remaining.load(std::memory_order_relaxed);
    if (next != -1) {
        attempts.fetch_add(1, std::memory_order_relaxed);
        if (next == 0) {
            // Disable before exception construction and stack unwinding.
            remaining.store(-1, std::memory_order_relaxed);
            throw std::bad_alloc{};
        }
        if (next > 0)
            remaining.fetch_sub(1, std::memory_order_relaxed);
    }
}
void* allocate(std::size_t count) {
    allocation_attempt();
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable pointee storage.
    void* const value = std::malloc(count == 0 ? 1 : count);
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
void* allocate_aligned(std::size_t count, std::size_t alignment) {
    allocation_attempt();
    count = count == 0 ? 1 : count;
#if defined(_WIN32)
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* const value = _aligned_malloc(count, alignment);
#else
    if (count > std::numeric_limits<std::size_t>::max() - (alignment - 1))
        throw std::bad_alloc{};
    const auto rounded = ((count + alignment - 1) / alignment) * alignment;
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable storage.
    void* const value = std::aligned_alloc(alignment, rounded);
#endif
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
void free_aligned(void* value) noexcept {
#if defined(_WIN32)
    _aligned_free(value);
#else
    std::free(value);
#endif
}
struct FaultScope {
    explicit FaultScope(std::int64_t count) noexcept {
        attempts.store(0, std::memory_order_relaxed);
        remaining.store(count, std::memory_order_relaxed);
    }
    ~FaultScope() {
        remaining.store(-1, std::memory_order_relaxed);
    }
    FaultScope(const FaultScope&) = delete;
    FaultScope& operator=(const FaultScope&) = delete;
};
} // namespace
void* operator new(std::size_t count) {
    return allocate(count);
}
void* operator new[](std::size_t count) {
    return allocate(count);
}
void operator delete(void* value) noexcept {
    std::free(value);
}
void operator delete[](void* value) noexcept {
    std::free(value);
}
void operator delete(void* value, std::size_t) noexcept {
    std::free(value);
}
void operator delete[](void* value, std::size_t) noexcept {
    std::free(value);
}
// Cover the forms used by library temporary buffers as well as user new
// expressions. Otherwise ASan's default nothrow new can pair with our free.
void* operator new(std::size_t count, const std::nothrow_t&) noexcept {
    try {
        return allocate(count);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void* operator new[](std::size_t count, const std::nothrow_t&) noexcept {
    try {
        return allocate(count);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void operator delete(void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void operator delete[](void* value, const std::nothrow_t&) noexcept {
    std::free(value);
}
void* operator new(std::size_t count, std::align_val_t alignment) {
    return allocate_aligned(count, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t count, std::align_val_t alignment) {
    return allocate_aligned(count, static_cast<std::size_t>(alignment));
}
void operator delete(void* value, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete(void* value, std::size_t, std::align_val_t) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::size_t, std::align_val_t) noexcept {
    free_aligned(value);
}
void* operator new(std::size_t count, std::align_val_t alignment, const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(count, static_cast<std::size_t>(alignment));
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void* operator new[](std::size_t count, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
    try {
        return allocate_aligned(count, static_cast<std::size_t>(alignment));
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}
void operator delete(void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}
void operator delete[](void* value, std::align_val_t, const std::nothrow_t&) noexcept {
    free_aligned(value);
}

namespace {
using namespace uvlas;
using namespace uvlas::io;
constexpr std::uint64_t budget = 64ULL * 1024 * 1024;
constexpr std::string_view source =
    R"({"schema":"uvlas.case/1","units":null,"model":["日本語",1.0,{"x":null}],"flow":null,"analysis":null,"observations":null,"numerics":null,"resources":null,"output":null})";
InputLimits limits() {
    return {100000,
            64,
            10000,
            10000,
            10000,
            10000,
            100000,
            {ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
}
core::Result<InputDocument> read(const InputLimits& request, core::ResourcePlan& resources,
                                 std::optional<InputDiagnostic>& diagnostic,
                                 core::FailureDetails details = {}) {
    return read_input_document(std::as_bytes(std::span{source.data(), source.size()}), SourceId{7},
                               request, resources, diagnostic, details);
}
} // namespace

TEST_CASE("Document fault allocator covers matching allocation and deallocation forms",
          "[unit][negative][N4-02][N4-08]") {
    constexpr std::size_t bytes = 32;
    constexpr std::align_val_t alignment{64};
    struct Route {
        const char* name;
        void* (*allocate)();
        void (*release)(void*);
        bool aligned;
        bool nothrow;
    };
    constexpr auto scalar = [] { return ::operator new(bytes); };
    constexpr auto array = [] { return ::operator new[](bytes); };
    constexpr auto scalar_nothrow = [] { return ::operator new(bytes, std::nothrow); };
    constexpr auto array_nothrow = [] { return ::operator new[](bytes, std::nothrow); };
    constexpr auto scalar_aligned = [] { return ::operator new(bytes, alignment); };
    constexpr auto array_aligned = [] { return ::operator new[](bytes, alignment); };
    constexpr auto scalar_aligned_nothrow = [] {
        return ::operator new(bytes, alignment, std::nothrow);
    };
    constexpr auto array_aligned_nothrow = [] {
        return ::operator new[](bytes, alignment, std::nothrow);
    };
    constexpr auto delete_scalar = [](void* value) { ::operator delete(value); };
    constexpr auto delete_array = [](void* value) { ::operator delete[](value); };
    constexpr auto delete_scalar_sized = [](void* value) { ::operator delete(value, bytes); };
    constexpr auto delete_array_sized = [](void* value) { ::operator delete[](value, bytes); };
    constexpr auto delete_scalar_nothrow = [](void* value) {
        ::operator delete(value, std::nothrow);
    };
    constexpr auto delete_array_nothrow = [](void* value) {
        ::operator delete[](value, std::nothrow);
    };
    constexpr auto delete_scalar_aligned = [](void* value) { ::operator delete(value, alignment); };
    constexpr auto delete_array_aligned = [](void* value) {
        ::operator delete[](value, alignment);
    };
    constexpr auto delete_scalar_aligned_sized = [](void* value) {
        ::operator delete(value, bytes, alignment);
    };
    constexpr auto delete_array_aligned_sized = [](void* value) {
        ::operator delete[](value, bytes, alignment);
    };
    constexpr auto delete_scalar_aligned_nothrow = [](void* value) {
        ::operator delete(value, alignment, std::nothrow);
    };
    constexpr auto delete_array_aligned_nothrow = [](void* value) {
        ::operator delete[](value, alignment, std::nothrow);
    };
    const std::array routes{
        Route{"scalar", scalar, delete_scalar, false, false},
        Route{"array", array, delete_array, false, false},
        Route{"scalar sized", scalar, delete_scalar_sized, false, false},
        Route{"array sized", array, delete_array_sized, false, false},
        Route{"scalar nothrow", scalar_nothrow, delete_scalar_nothrow, false, true},
        Route{"array nothrow", array_nothrow, delete_array_nothrow, false, true},
        Route{"scalar aligned", scalar_aligned, delete_scalar_aligned, true, false},
        Route{"array aligned", array_aligned, delete_array_aligned, true, false},
        Route{"scalar aligned sized", scalar_aligned, delete_scalar_aligned_sized, true, false},
        Route{"array aligned sized", array_aligned, delete_array_aligned_sized, true, false},
        Route{"scalar aligned nothrow", scalar_aligned_nothrow, delete_scalar_aligned_nothrow, true,
              true},
        Route{"array aligned nothrow", array_aligned_nothrow, delete_array_aligned_nothrow, true,
              true},
        // libstdc++ temporary buffers allocate with nothrow new, then use a
        // normal sized delete. This must remain in the same allocation family.
        Route{"scalar nothrow sized", scalar_nothrow, delete_scalar_sized, false, true},
        Route{"array nothrow sized", array_nothrow, delete_array_sized, false, true},
        Route{"scalar aligned nothrow sized", scalar_aligned_nothrow, delete_scalar_aligned_sized,
              true, true},
        Route{"array aligned nothrow sized", array_aligned_nothrow, delete_array_aligned_sized,
              true, true}};
    for (const auto& route : routes) {
        CAPTURE(route.name);
        void* const memory = [&] {
            const FaultScope count{-2};
            return route.allocate();
        }();
        REQUIRE(memory != nullptr);
        CHECK(attempts.load(std::memory_order_relaxed) == 1);
        if (route.aligned) {
            // NOLINTNEXTLINE(misc-const-correctness): std::align takes void*&.
            void* candidate = memory;
            auto available = bytes;
            CHECK(std::align(static_cast<std::size_t>(alignment), bytes, candidate, available) ==
                  memory);
        }
        auto* const writable = static_cast<std::byte*>(memory);
        writable[0] = std::byte{0x5a};
        writable[bytes - 1] = std::byte{0xa5};
        CHECK(writable[0] == std::byte{0x5a});
        CHECK(writable[bytes - 1] == std::byte{0xa5});
        route.release(memory);
        bool threw = false;
        void* const failed = [&] {
            const FaultScope fault{0};
            try {
                return route.allocate();
            } catch (const std::bad_alloc&) {
                threw = true;
                return static_cast<void*>(nullptr);
            }
        }();
        route.release(failed);
        CHECK(failed == nullptr);
        CHECK(threw == !route.nothrow);
        CHECK(attempts.load(std::memory_order_relaxed) == 1);
    }
    void* const empty = [&] {
        const FaultScope count{-2};
        return ::operator new(0, std::nothrow);
    }();
    CHECK(empty != nullptr);
    CHECK(attempts.load(std::memory_order_relaxed) == 1);
    ::operator delete(empty, std::nothrow);
}

TEST_CASE("Document allocation failures preserve input old owners and cumulative plan charges",
          "[negative][IO4-01][N4-02][N4-03][N4-08]") {
    core::ResourcePlan original_plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    const auto old = read(limits(), original_plan, diagnostic);
    REQUIRE(old);
    const auto old_view = old->view();
    REQUIRE(old_view);
    const auto identity = old_view->identity;
    const auto original_bytes = std::as_bytes(std::span{source.data(), source.size()});
    core::ResourcePlan control_plan{budget};
    const auto control = [&] {
        const FaultScope count{-2};
        return read(limits(), control_plan, diagnostic);
    }();
    REQUIRE(control);
    const auto count = attempts.load(std::memory_order_relaxed);
    REQUIRE(count > 8); // Product arrays plus the actual external SAX path.
    const core::FailureDetails context{core::ModelGeneration{1},
                                       core::AcceptedStateId{2},
                                       core::TrialId{3},
                                       core::NumericalGeneration{4},
                                       core::ComponentId{5},
                                       core::ProblemId{core::PanelId{6}},
                                       std::nullopt,
                                       core::AcceptedStateId{9}};
    bool retained_successful_charge = false;
    for (std::uint64_t site = 0; site < count; ++site) {
        core::ResourcePlan plan{budget};
        const auto prior = plan.allocate<std::byte>(3, core::Operation::PrepareModel);
        REQUIRE(prior);
        const auto failed = [&] {
            const FaultScope fault{static_cast<std::int64_t>(site)};
            return read(limits(), plan, diagnostic, context);
        }();
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
        CHECK(failed.error().operation == core::Operation::PrepareModel);
        CHECK(failed.error().details.model == context.model);
        CHECK(failed.error().details.parent == context.parent);
        CHECK(failed.error().details.trial == context.trial);
        CHECK(failed.error().details.numerical_generation == context.numerical_generation);
        CHECK(failed.error().details.component == context.component);
        CHECK(failed.error().details.problem == context.problem);
        CHECK(failed.error().details.last_accepted == context.last_accepted);
        if (!diagnostic) {
            FAIL_CHECK("Missing allocation diagnostic");
            return;
        }
        CHECK(diagnostic->reason == InputDiagnosticReason::Allocation);
        CHECK(diagnostic->source.value() == 7);
        CHECK(plan.allocated_bytes() >= 3);
        retained_successful_charge = retained_successful_charge || plan.allocated_bytes() > 3;
        CHECK(old->view()->identity == identity);
        CHECK(std::equal(old_view->source_bytes.begin(), old_view->source_bytes.end(),
                         original_bytes.begin()));
    }
    CHECK(retained_successful_charge);
    core::ResourcePlan recovery{budget};
    CHECK(read(limits(), recovery, diagnostic));
    CHECK_FALSE(diagnostic);
    std::printf("input-document actual C++new fault sites=%llu\n",
                static_cast<unsigned long long>(count));
}

TEST_CASE("Document entry refuses known policies without allocation",
          "[negative][IO4-01][N4-02][N4-08]") {
    core::ResourcePlan plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    auto strict = limits();
    strict.external_memory.mode = ParserAllocationMode::RequireProvenBound;
    const auto rejected = [&] {
        const FaultScope count{-2};
        return read(strict, plan, diagnostic);
    }();
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::Unsupported);
    CHECK(attempts.load(std::memory_order_relaxed) == 0);
    CHECK(plan.allocated_bytes() == 0);
    auto shortage = limits();
    shortage.external_memory.maximum_bytes = 0;
    const auto too_small = [&] {
        const FaultScope count{-2};
        return read(shortage, plan, diagnostic);
    }();
    REQUIRE_FALSE(too_small);
    CHECK(too_small.error().code == core::FailureCode::ResourceExceeded);
    CHECK(attempts.load(std::memory_order_relaxed) == 0);
}

TEST_CASE(
    "Parser message fallback classification owns no third party message and allocates nothing",
    "[unit][negative][IO4-01][N4-08]") {
    using Kind = io::detail::ParserFailureClass;
    using Code = core::FailureCode;
    // This reaches the real private classifier. It does NOT inject CRT malloc.
    for (const auto kind : {Kind::Syntax, Kind::NumberRange}) {
        const auto result = [&] {
            const FaultScope count{-2};
            return io::detail::classify_parser_failure(kind, kind == Kind::Syntax ? 101 : 406,
                                                       "Unknown exception");
        }();
        CHECK(result.code == Code::ResourceExceeded);
        CHECK(result.reason == InputDiagnosticReason::Allocation);
        CHECK(attempts.load(std::memory_order_relaxed) == 0);
    }
    const auto foreign = io::detail::classify_parser_failure(Kind::Other, 101, "Unknown exception");
    CHECK(foreign.code == Code::InvalidInput);
    CHECK(foreign.reason == InputDiagnosticReason::UnexpectedException);
    const auto range = io::detail::classify_parser_failure(
        Kind::NumberRange, 406, "[json.exception.out_of_range.406] fixed");
    CHECK(range.code == Code::Nonfinite);
    CHECK(range.reason == InputDiagnosticReason::NumberRange);
    const auto malformed =
        io::detail::classify_parser_failure(Kind::Syntax, 101, "unexpected message");
    CHECK(malformed.reason == InputDiagnosticReason::UnexpectedException);
}
