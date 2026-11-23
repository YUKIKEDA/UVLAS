#pragma once

#include <uvlas/io/input_document.hpp>

namespace uvlas::io::detail {

/// Read a strict JSON object with the existing source and parser quotas. The caller
/// validates the document-specific envelope before passing it to any consumer.
[[nodiscard]] core::Result<InputDocument>
read_json_object(std::span<const std::byte>, SourceId, const InputLimits&, core::ResourcePlan&,
                 std::optional<InputDiagnostic>&) noexcept;

} // namespace uvlas::io::detail
