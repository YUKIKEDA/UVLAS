#pragma once

#include <boost/hash2/sha2.hpp>

#include <span>
#include <uvlas/io/input_document.hpp>

namespace uvlas::io::detail {

[[nodiscard]] inline SourceDigest source_digest(std::span<const std::byte> bytes) noexcept {
    boost::hash2::sha2_256 hash;
    if (!bytes.empty()) {
        hash.update(bytes.data(), bytes.size());
    }
    const auto result = hash.result();
    SourceDigest digest{};
    for (std::size_t i = 0; i < digest.bytes.size(); ++i) {
        digest.bytes[i] = static_cast<std::byte>(result[i]);
    }
    return digest;
}

} // namespace uvlas::io::detail
