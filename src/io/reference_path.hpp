#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stop_token>
#include <string_view>
#include <uvlas/core/resource_plan.hpp>

namespace uvlas::io::detail {

struct ReferencePathLimits {
    std::size_t input_bytes;
    std::size_t normalized_bytes;
    std::size_t normalized_components;
};

enum class ReferencePathReason : std::uint8_t {
    EscapeRoot,
    RootedPath,
    StreamPath,
    MissingFile,
    InvalidPath,
    InputBytes,
    NormalizedBytes,
    NormalizedDepth,
    PlanBytes,
};

struct ReferencePathDiagnostic {
    ReferencePathDiagnostic(ReferencePathReason failure, std::optional<std::size_t> offset) noexcept
        : reason(failure), input_offset(offset) {}

    ReferencePathReason reason;
    std::optional<std::size_t> input_offset;
};

// The owner contains one nonempty slash-separated path and no OS-specific state.
// A moved-from owner is empty. The caller keeps the original declaration separately.
class NormalizedReferencePath {
public:
    [[nodiscard]] std::string_view text() const& noexcept;
    [[nodiscard]] std::string_view text() const&& = delete;

private:
    friend core::Result<NormalizedReferencePath>
    normalize_reference_path(std::string_view, ReferencePathLimits, core::ResourcePlan&,
                             const std::stop_token&,
                             std::optional<ReferencePathDiagnostic>*) noexcept;

    explicit NormalizedReferencePath(core::PlannedBuffer<char> bytes) noexcept;
    core::PlannedBuffer<char> bytes_;
};

// Input is UTF-8 already validated by the manifest or role-specific decoder. Only the owned
// normalized bytes are charged to plan. No OS resource is opened or inspected. On a lexical or
// resource failure, diagnostic receives an owned reason and optional zero-based input byte offset.
// Cancellation uses the common failure code. Success and cancellation clear the diagnostic.
[[nodiscard]] core::Result<NormalizedReferencePath>
normalize_reference_path(std::string_view input, ReferencePathLimits limits,
                         core::ResourcePlan& plan, const std::stop_token& stop,
                         std::optional<ReferencePathDiagnostic>* diagnostic) noexcept;

} // namespace uvlas::io::detail
