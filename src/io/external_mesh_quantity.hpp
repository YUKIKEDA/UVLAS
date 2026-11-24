#pragma once

#include "nastran_laminate_cards.hpp"

#include <stop_token>
#include <string_view>
#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail {

// A thickness comparison is admitted only after correspondence of the physical regions.
// The external decimal is a single full-thickness value declared over the entire region.
struct ImportedUniformThickness {
    std::string_view decimal;
    SiLengthUnit length;
    std::optional<NastranBulkPosition> position{};
};

enum class ImportedQuantityField : std::uint8_t {
    Young1,
    Young2,
    Shear23,
    Shear13,
    Shear12,
    Poisson12,
    Density,
    FullThickness,
};

struct ImportedQuantityMismatch {
    ImportedQuantityField quantity;
    std::optional<NastranBulkPosition> card_position;
    // For a case thickness field, the index of the first different vertex.
    std::optional<core::Index> vertex;
};

using ImportedCaseMaterial =
    std::variant<schema::MaterialIsotropic, schema::MaterialOrthotropic, schema::MaterialEffective>;

[[nodiscard]] core::Result<void> check_imported_uniform_thickness(
    ImportedUniformThickness, const schema::Thickness&, SiLengthUnit case_length,
    std::uint64_t maximum_vertices, const structure::ExactLimits&, core::ResourcePlan&,
    const std::stop_token&, std::optional<ImportedQuantityMismatch>* = nullptr) noexcept;

/// Compare the uniform PSHELL thickness for an element without explicit corner overrides.
/// Element and property belong to the same authenticated deck. The admitting caller
/// has established the physical correspondence with the case thickness region.
[[nodiscard]] core::Result<void> check_imported_pshell_thickness(
    const NastranShellElement&, const NastranPshell&, SiLengthUnit deck_length,
    const schema::Thickness&, SiLengthUnit case_length, std::uint64_t maximum_vertices,
    const structure::ExactLimits&, core::ResourcePlan&, const std::stop_token&,
    std::optional<ImportedQuantityMismatch>* = nullptr) noexcept;

[[nodiscard]] core::Result<void>
check_imported_material(const NastranMat1&, const ImportedCaseMaterial&, const SiUnits& deck_units,
                        const SiUnits& case_units, const structure::ExactLimits&,
                        core::ResourcePlan&, const std::stop_token&,
                        std::optional<ImportedQuantityMismatch>* = nullptr) noexcept;

[[nodiscard]] core::Result<void>
check_imported_material(const NastranMat8&, const ImportedCaseMaterial&, const SiUnits& deck_units,
                        const SiUnits& case_units, const structure::ExactLimits&,
                        core::ResourcePlan&, const std::stop_token&,
                        std::optional<ImportedQuantityMismatch>* = nullptr) noexcept;

} // namespace uvlas::io::detail
