#pragma once

#include "reference_resolution.hpp"

#include <filesystem>
#include <optional>
#include <stop_token>
#include <string>
#include <string_view>
#include <uvlas/core/failure.hpp>

namespace uvlas::app {

struct CaseValidationFailure {
    core::Failure failure;
    std::string_view stage;
    std::string location;
    std::optional<std::size_t> byte;
};

struct ValidatedCase {
    io::detail::ReferenceInput input;
    io::detail::ReferencedCase references;
    io::SchemaIdentity identity;
};

/// Check the case schema and the external reference graph without writing files.
[[nodiscard]] std::expected<ValidatedCase, CaseValidationFailure>
validate_case_file(const std::filesystem::path&, const std::stop_token& = {}) noexcept;

int validate_case(const std::filesystem::path&);

} // namespace uvlas::app
