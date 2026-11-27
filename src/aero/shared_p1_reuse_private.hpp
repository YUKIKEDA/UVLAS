#pragma once

#include "exact_triangle_predicates.hpp"

#include <optional>

namespace uvlas::aero::detail {
class LocalPacketEngine;
class P1Engine;
class SharedP1Engine;

/// A slot private to a mixed-action assembly. Its previous numerical contents
/// can be overwritten on failure because no accessor publishes this slot.
template <class Result>
class SharedP1ReuseSlot {
    friend class LocalPacketEngine;
    friend class P1Engine;
    friend class SharedP1Engine;

    explicit SharedP1ReuseSlot(std::optional<Result>& value,
                               TriangleOutcome* predicate = nullptr) noexcept
        : value_(value), predicate_(predicate) {}

    std::optional<Result>& value_;
    TriangleOutcome* predicate_;
};
} // namespace uvlas::aero::detail
