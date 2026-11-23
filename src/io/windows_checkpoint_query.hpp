#pragma once

#include <optional>
#include <uvlas/core/failure.hpp>
#include <windows.h>

namespace uvlas::io::detail {

/// Classify a failed directory query using the error captured at the call site.
/// No failure means enumeration reached its end.
[[nodiscard]] constexpr std::optional<core::FailureCode>
windows_checkpoint_query_failure(DWORD error) noexcept {
    switch (error) {
    case ERROR_NO_MORE_FILES:
        return std::nullopt;
    case ERROR_NOT_ENOUGH_MEMORY:
    case ERROR_OUTOFMEMORY:
    case ERROR_NO_SYSTEM_RESOURCES:
    case ERROR_NOT_ENOUGH_QUOTA:
    case ERROR_COMMITMENT_LIMIT:
        return core::FailureCode::ResourceExceeded;
    default:
        return core::FailureCode::CheckpointCorrupt;
    }
}

} // namespace uvlas::io::detail
