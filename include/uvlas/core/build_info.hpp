#pragma once

#include <string_view>

namespace uvlas {

struct BuildInfo {
    std::string_view version;
    std::string_view compiler;
    std::string_view platform;
    std::string_view configuration;
};

// String views refer to static storage. This is an internal, unstable API.
[[nodiscard]] BuildInfo build_info() noexcept;

} // namespace uvlas
