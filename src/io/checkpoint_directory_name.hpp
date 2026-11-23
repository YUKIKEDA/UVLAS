#pragma once

#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <limits>
#include <optional>
#include <string_view>

namespace uvlas::io::detail {

inline constexpr std::string_view checkpoint_directory_prefix = "checkpoint-";
using CheckpointDirectoryName =
    std::array<char, checkpoint_directory_prefix.size() +
                         std::numeric_limits<std::uint64_t>::digits10 + 2>;

[[nodiscard]] inline CheckpointDirectoryName
make_checkpoint_directory_name(std::uint64_t number) noexcept {
    CheckpointDirectoryName result{};
    const auto end = std::ranges::copy(checkpoint_directory_prefix, result.data()).out;
    static_cast<void>(std::to_chars(end, result.data() + result.size() - 1, number));
    return result;
}

[[nodiscard]] inline std::optional<std::uint64_t>
parse_checkpoint_directory_name(std::string_view name) noexcept {
    if (!name.starts_with(checkpoint_directory_prefix)) {
        return std::nullopt;
    }
    name.remove_prefix(checkpoint_directory_prefix.size());
    if (name.empty() || (name.front() == '0' && name.size() > 1)) {
        return std::nullopt;
    }
    std::uint64_t number = 0;
    const auto end = name.data() + name.size();
    const auto converted = std::from_chars(name.data(), end, number);
    if (converted.ec != std::errc{} || converted.ptr != end) {
        return std::nullopt;
    }
    return number;
}

} // namespace uvlas::io::detail
