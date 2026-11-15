#pragma once
#include <uvlas/numerics/backend_memory.hpp>

namespace uvlas::numerics::linear::detail {
// Only backend factories include this private header. Public callers cannot manufacture
// a proof by supplying an arbitrary number to the admission function.
class BackendPlanFactory {
public:
    [[nodiscard]] static BackendMemoryPlan
    make(std::uint64_t estimate, std::optional<std::uint64_t> proven = std::nullopt) noexcept {
        return BackendMemoryPlan{estimate, proven};
    }
};
} // namespace uvlas::numerics::linear::detail
