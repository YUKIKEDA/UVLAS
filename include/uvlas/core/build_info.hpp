#pragma once

#include <string>
#include <string_view>

namespace uvlas {

/// Version, source revision and build settings compiled into the running binary.
/// All string views borrow static storage and outlive the returned BuildInfo value.
struct BuildInfo {
    /// Project version compiled into this binary.
    std::string_view version;
    /// Source commit recorded at build time.
    std::string_view commit;
    /// Whether the build detected uncommitted source changes.
    bool uncommitted_changes;
    /// Compiler identity and version.
    std::string_view compiler;
    /// Target operating system and architecture.
    std::string_view platform;
    /// Build configuration, such as Debug or Release.
    std::string_view configuration;
    /// Floating-point compilation policy.
    std::string_view floating_point;
    /// Target instruction-set level.
    std::string_view instruction_set;
    /// Key combining OS, architecture, compiler major.minor and instruction set.
    /// Use it to distinguish environments when comparing numerical results.
    std::string_view platform_key;
};

/// Return binary provenance by value without allocation; views borrow static storage.
[[nodiscard]] BuildInfo build_info() noexcept;

/// Return owned JSON containing provenance, dependency versions and source-repair revisions.
/// @throws std::bad_alloc If serialization cannot allocate memory.
/// Allocation exceptions propagate to the caller; catch them if invoking this
/// diagnostic function from an exception-free application boundary.
[[nodiscard]] std::string build_info_json();

} // namespace uvlas
