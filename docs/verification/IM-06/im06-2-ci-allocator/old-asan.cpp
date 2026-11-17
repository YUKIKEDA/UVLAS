#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <memory>
#include <new>
namespace {
// Process-local injection belongs only to this dedicated test executable.
enum class InjectedFailure : std::uint8_t { Allocation, Unexpected, Nonstandard };
struct NonstandardException {};
constinit std::atomic<std::int64_t> remaining{-1};
constinit std::atomic<std::uint64_t> attempts{0};
constinit InjectedFailure injected = InjectedFailure::Allocation;
class ForeignException final : public std::exception {
public:
    [[nodiscard]] const char* what() const noexcept override {
        return "fixture exception";
    }
};
void* allocate(std::size_t count) {
    const auto next = remaining.load(std::memory_order_relaxed);
    if (next != -1) {
        attempts.fetch_add(1, std::memory_order_relaxed);
        if (next == 0) {
            // Disable before exception construction and stack unwinding.
            remaining.store(-1, std::memory_order_relaxed);
            if (injected == InjectedFailure::Unexpected)
                throw ForeignException{};
            if (injected == InjectedFailure::Nonstandard)
                throw NonstandardException{};
            throw std::bad_alloc{};
        }
        if (next > 0)
            remaining.fetch_sub(1, std::memory_order_relaxed);
    }
    // NOLINTNEXTLINE(misc-const-correctness): operator new must return writable pointee storage.
    void* const value = std::malloc(count == 0 ? 1 : count);
    if (value == nullptr)
        throw std::bad_alloc{};
    return value;
}
struct FaultScope {
    explicit FaultScope(std::int64_t count,
                        InjectedFailure kind = InjectedFailure::Allocation) noexcept {
        attempts.store(0, std::memory_order_relaxed);
        injected = kind;
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


unsigned check_count = 0;
void check(bool value, unsigned line) {
    ++check_count;
    if (!value) {
        std::fprintf(stderr, "failed control at line %u\n", line);
        std::exit(42);
    }
}
#define TEST_CASE(...) void check_routes()
#define CHECK(...) check((__VA_ARGS__), __LINE__)
#define REQUIRE(...) check((__VA_ARGS__), __LINE__)
#define CAPTURE(...) static_cast<void>(__VA_ARGS__)
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


int main() {
    // Match the library's nothrow temporary-buffer allocation + sized delete.
    void* const memory = ::operator new(32, std::nothrow);
    ::operator delete(memory, std::size_t{32});
    check_routes();
    std::printf("allocation control checks=%u\n", check_count);
}
