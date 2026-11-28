#pragma once

#include "reference_resolution.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <span>
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
    std::optional<io::InputNodeOrdinal> reference_node{};
};

struct ValidatedCase {
    io::detail::ReferenceInput input;
    io::detail::ReferencedCase references;
    io::SchemaIdentity identity;
    io::detail::ObservationLeafResults observations;
};

inline constexpr std::size_t case_manifest_max_bytes = std::size_t{8} * 1024 * 1024;
inline constexpr std::size_t case_manifest_max_depth = 128;
inline constexpr std::size_t case_manifest_max_nodes = std::size_t{256} * 1024;
inline constexpr std::size_t case_manifest_max_string_bytes = std::size_t{64} * 1024;
inline constexpr std::size_t case_manifest_max_number_bytes = 1024;

/// Check the case schema and the external reference graph without writing files.
[[nodiscard]] std::expected<ValidatedCase, CaseValidationFailure>
validate_case_file(const std::filesystem::path&, const std::stop_token& = {}) noexcept;

/// Check a copied manifest with an optional directory for its external references.
[[nodiscard]] std::expected<ValidatedCase, CaseValidationFailure>
validate_case_memory(std::span<const std::byte>, const std::optional<std::filesystem::path>&,
                     const std::stop_token& = {}) noexcept;

int validate_case(const std::filesystem::path&);

} // namespace uvlas::app
