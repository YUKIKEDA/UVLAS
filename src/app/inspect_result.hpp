#pragma once

#include "reference_file_snapshot.hpp"

#include <filesystem>
#include <stop_token>
#include <string_view>

namespace uvlas::app {

struct ResultInspectionFailure {
    core::Failure failure;
    std::string_view subject;
};

/// Return the original status bytes of a result whose metadata and optional resolved case agree.
[[nodiscard]] std::expected<io::detail::ReferenceFileSnapshot, ResultInspectionFailure>
inspect_result_file(const std::filesystem::path&, const std::stop_token& = {}) noexcept;

int inspect_result(const std::filesystem::path&);

} // namespace uvlas::app
