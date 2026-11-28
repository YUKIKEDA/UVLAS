#pragma once

#include <boost/hash2/sha2.hpp>

#include <optional>
#include <span>
#include <string_view>
#include <uvlas/io/input_document.hpp>

namespace uvlas::io::detail {

[[nodiscard]] inline std::optional<SourceDigest>
decode_source_digest(std::string_view text) noexcept {
    constexpr std::string_view digits = "0123456789abcdef";
    SourceDigest result{};
    if (text.size() != result.bytes.size() * 2)
        return std::nullopt;
    for (std::size_t i = 0; i < result.bytes.size(); ++i) {
        const auto high = digits.find(text[2 * i]);
        const auto low = digits.find(text[2 * i + 1]);
        if (high == digits.npos || low == digits.npos)
            return std::nullopt;
        result.bytes[i] = static_cast<std::byte>((high << 4) | low);
    }
    return result;
}

[[nodiscard]] inline SourceDigest finish_source_digest(boost::hash2::sha2_256& hash) noexcept {
    const auto result = hash.result();
    SourceDigest digest{};
    for (std::size_t i = 0; i < digest.bytes.size(); ++i) {
        digest.bytes[i] = static_cast<std::byte>(result[i]);
    }
    return digest;
}

[[nodiscard]] inline SourceDigest source_digest(std::span<const std::byte> bytes) noexcept {
    boost::hash2::sha2_256 hash;
    if (!bytes.empty()) {
        hash.update(bytes.data(), bytes.size());
    }
    return finish_source_digest(hash);
}

} // namespace uvlas::io::detail
