#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <stop_token>
#include <uvlas/structure/pg_coefficients.hpp>
#include <uvlas/structure/pg_integration_partition.hpp>

namespace uvlas::structure {

/// Allocations, addressable capacity and work of a single local rank preparation.
struct PgRankPreparationLimits {
    std::uint64_t maximum_payload_bytes;
    std::uint64_t maximum_steps;
    std::uint32_t maximum_digits;
    std::stop_token cancellation;
};

/// One activity-row coefficient, in the dictionary's unscaled rational coordinates.
struct PgRankRow {
    PgCoefficientField field = PgCoefficientField::R;
    PgCellIndex cell;
    ExactPowers powers;
};

struct PgRankBasisColumn {
    core::Index value;
};
struct PgRankOriginalColumn {
    core::Index value;
};
struct PgRankRequiredColumn {
    core::Index value;
};
struct PgRankAnchorColumn {
    core::Index value;
};

/// Success records actual work, capacities and the finite-field witness separately.
struct PgRankPreparationInfo {
    std::uint64_t active_rows = 0;
    std::uint64_t anchor_rows = 0;
    std::uint64_t candidates = 0;
    std::uint64_t normalized_terms = 0;
    std::uint64_t rank = 0;
    std::uint64_t rank_c = 0;
    std::uint64_t independent_required = 0;
    std::uint64_t proof_rows = 0;
    std::uint64_t verified_anchor_rows = 0;
    std::uint64_t steps = 0;
    std::uint64_t preflight_bytes = 0;
    std::uint64_t allocated_payload_bytes = 0;
    std::uint32_t prime = 0;
};

/// Cumulative validation work and replay payload; the original rank remains owned.
struct PgRankValidationInfo {
    std::uint64_t steps = 0;
    std::uint64_t preflight_bytes = 0;
    std::uint64_t allocated_payload_bytes = 0;
};

/// The chart, source polygon and per-cell activity retained with the basis.
struct PgRankWitnessInfo {
    core::Index n = 0;
    PgGraphConvention convention = PgGraphConvention::Version1;
    std::uint32_t recipe_version = 0;
    std::uint32_t reduction_version = 0;
    std::uint64_t height_terms = 0;
    std::uint64_t domain_vertices = 0;
    std::uint64_t activity_cells = 0;
    PgGraphContentIdentity graph_identity{{}};
};

struct PgRankChartParameters {
    ExactRational x0, width, skew;
};
struct PgRankHeightTerm {
    ExactPowers powers;
    ExactRational coefficient;
};
struct PgRankDomainPoint {
    ExactRational s, t;
};
struct PgRankActivityCell {
    PgCellIndex cell;
    ExactRational area;
    bool active;
};

/// An ordered term of G or J in raw dictionary coordinates, with its algebraic sign.
struct PgRankRecipeTerm {
    PgRankOriginalColumn source;
    ExactRational coefficient;
};

/// Distinct Cartesian directions of required translation columns 30, 31 and 32.
enum class PgRankTranslationAxis : std::uint8_t { X, Y, Z };

namespace detail {
struct PgRankState;
}

/// Owns one bounded, unconditioned local basis and its exact reconstruction maps.
/// Source owners remain unchanged. Independent certificates are test inputs only.
class PgRankPreparation {
public:
    PgRankPreparation(PgRankPreparation&&) noexcept;
    PgRankPreparation& operator=(PgRankPreparation&&) noexcept;
    PgRankPreparation(const PgRankPreparation&) = delete;
    PgRankPreparation& operator=(const PgRankPreparation&) = delete;
    ~PgRankPreparation();

    [[nodiscard]] static core::Result<PgRankPreparation>
    prepare(const PgCoefficientDictionary&, const PgIntegrationPartition&,
            const PgRankPreparationLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    /// Recheck the source witness and every reconstruction against live factory inputs.
    /// The caller's plan admits new work while the original result remains resident.
    [[nodiscard]] core::Result<PgRankValidationInfo>
    validate_for(const PgCoefficientDictionary&, const PgIntegrationPartition&,
                 const PgRankPreparationLimits&, core::Workspace&,
                 core::ResourcePlan&) const noexcept;

    [[nodiscard]] core::Result<std::span<const PgRankRow>> rows() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgRankRow>> rows() const&& = delete;
    /// The selected ordinal in [G,I,J], in greedy order. J never extends the basis.
    [[nodiscard]] core::Result<std::span<const std::uint32_t>> selected() const& noexcept;
    [[nodiscard]] core::Result<std::span<const std::uint32_t>> selected() const&& = delete;
    [[nodiscard]] core::Result<PgRankPreparationInfo> preparation() const noexcept;

    /// Source values remain inspectable after the dictionary and partition are released.
    [[nodiscard]] core::Result<PgRankWitnessInfo> witness() const noexcept;
    [[nodiscard]] core::Result<PgRankChartParameters>
    source_chart(const ExactLimits&, core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<PgRankHeightTerm> source_height(core::Index, const ExactLimits&,
                                                               core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<PgRankDomainPoint> source_domain(core::Index, const ExactLimits&,
                                                                core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<PgRankActivityCell>
    source_activity(core::Index, const ExactLimits&, core::ResourcePlan&) const noexcept;

    /// Ordered provenance for every raw and appended dictionary column.
    [[nodiscard]] core::Result<std::span<const PgCoefficientSource>> provenance() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgCoefficientSource>> provenance() const&& = delete;

    /// Raw G and J recipes. Repeated source terms retain their order and sign.
    [[nodiscard]] core::Result<std::uint64_t> recipe_terms(PgRankRequiredColumn) const noexcept;
    [[nodiscard]] core::Result<std::uint64_t> recipe_terms(PgRankAnchorColumn) const noexcept;
    [[nodiscard]] core::Result<PgRankRecipeTerm> recipe(PgRankRequiredColumn, core::Index term,
                                                        const ExactLimits&,
                                                        core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<PgRankRecipeTerm> recipe(PgRankAnchorColumn, core::Index term,
                                                        const ExactLimits&,
                                                        core::ResourcePlan&) const noexcept;

    /// Distinguish translations that share one scalar recipe. Other targets have no axis.
    [[nodiscard]] core::Result<std::optional<PgRankTranslationAxis>>
        translation_axis(PgRankRequiredColumn) const noexcept;

    /// B in the original field, cell and polynomial units. Zero is also exact.
    [[nodiscard]] core::Result<ExactRational> basis(PgRankBasisColumn, core::Index row,
                                                    const ExactLimits&,
                                                    core::ResourcePlan&) const noexcept;
    /// T pulls one basis column from the original source dictionary.
    [[nodiscard]] core::Result<ExactRational> source_map(PgRankBasisColumn, PgRankOriginalColumn,
                                                         const ExactLimits&,
                                                         core::ResourcePlan&) const noexcept;
    /// S, W and Z expand raw, required and full-background-anchor-kernel candidates.
    [[nodiscard]] core::Result<ExactRational> original(PgRankBasisColumn, PgRankOriginalColumn,
                                                       const ExactLimits&,
                                                       core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<ExactRational> required(PgRankBasisColumn, PgRankRequiredColumn,
                                                       const ExactLimits&,
                                                       core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<ExactRational> anchor(PgRankBasisColumn, PgRankAnchorColumn,
                                                     const ExactLimits&,
                                                     core::ResourcePlan&) const noexcept;

private:
    explicit PgRankPreparation(core::PlannedBuffer<detail::PgRankState>) noexcept;
    core::PlannedBuffer<detail::PgRankState> state_;
};

} // namespace uvlas::structure
