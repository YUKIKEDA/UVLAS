#pragma once

#include <string_view>

namespace uvlas::io::detail {

[[nodiscard]] inline bool valid_utf8(std::string_view text) noexcept {
    while (!text.empty()) {
        const auto first = static_cast<unsigned char>(text.front());
        if (first < 0x80) {
            text.remove_prefix(1);
            continue;
        }

        std::size_t width = 0;
        unsigned char low = 0x80, high = 0xbf;
        if (first >= 0xc2 && first <= 0xdf) {
            width = 2;
        } else if (first >= 0xe0 && first <= 0xef) {
            width = 3;
            low = first == 0xe0 ? 0xa0 : low;
            high = first == 0xed ? 0x9f : high;
        } else if (first >= 0xf0 && first <= 0xf4) {
            width = 4;
            low = first == 0xf0 ? 0x90 : low;
            high = first == 0xf4 ? 0x8f : high;
        } else {
            return false;
        }

        if (text.size() < width) {
            return false;
        }
        const auto second = static_cast<unsigned char>(text[1]);
        if (second < low || second > high) {
            return false;
        }
        for (std::size_t i = 2; i < width; ++i) {
            const auto byte = static_cast<unsigned char>(text[i]);
            if (byte < 0x80 || byte > 0xbf) {
                return false;
            }
        }
        text.remove_prefix(width);
    }
    return true;
}

} // namespace uvlas::io::detail
