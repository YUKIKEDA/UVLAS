#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace uvlas::test_support {

/// Execute the same computation at each positive worker count and compare owned result bytes.
/// The runner must actually honor its worker count and return only initialized numerical
/// bytes, not pointer metadata or structure padding. At least two distinct positive counts
/// are required; otherwise no comparison across worker counts has been established.
template <class Runner>
[[nodiscard]] bool identical_across_threads(std::span<const std::size_t> counts, Runner&& runner) {
    if (counts.empty() || std::ranges::find(counts, 0) != counts.end()) {
        return false;
    }
    if (std::ranges::all_of(
            counts, [first = counts.front()](std::size_t count) { return count == first; })) {
        return false;
    }
    const std::vector<std::byte> reference = std::invoke(runner, counts.front());
    for (const auto count : counts.subspan(1)) {
        const std::vector<std::byte> result = std::invoke(runner, count);
        if (result != reference) {
            return false;
        }
    }
    return true;
}
} // namespace uvlas::test_support
