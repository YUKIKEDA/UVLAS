#pragma once

#include <string>
#include <string_view>

namespace uvlas {

// Where a binary came from (QA-06). String views refer to static storage.
// This is an internal, unstable API.
struct BuildInfo {
    std::string_view version;
    std::string_view commit;
    bool uncommitted_changes;
    std::string_view compiler;
    std::string_view platform;
    std::string_view configuration;
    std::string_view floating_point;
    std::string_view instruction_set;
    // Operating system, architecture, compiler with major.minor, instruction
    // set: the unit under which numerical results are qualified (QA-07).
    std::string_view platform_key;
};

[[nodiscard]] BuildInfo build_info() noexcept;

// The build information and the versions of the numerical dependencies, as JSON.
[[nodiscard]] std::string build_info_json();

} // namespace uvlas
