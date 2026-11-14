#pragma once

#include <string>
#include <string_view>

namespace uvlas {

/// Binary provenance (QA-06); an internal API whose compatibility is not yet fixed.
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
    /// OS, architecture, compiler major.minor and instruction-set qualification key (QA-07).
    std::string_view platform_key;
};

/// Return binary provenance by value without allocation; views borrow static storage.
[[nodiscard]] BuildInfo build_info() noexcept;

/// Return owned JSON containing provenance and numerical dependency versions.
/// @throws std::bad_alloc If serialization cannot allocate memory.
/// This diagnostic helper is not a solver entry; it does not return Result<T>.
[[nodiscard]] std::string build_info_json();

} // namespace uvlas
