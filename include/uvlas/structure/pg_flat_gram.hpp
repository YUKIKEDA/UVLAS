#pragma once

#include <uvlas/structure/pg_conditioning.hpp>

namespace uvlas::structure {

struct PgFlatGramLimits {
    PgReferenceLimits arithmetic;
    std::uint64_t maximum_integer_steps = 0;
};

/// Density-one geometric Gram enclosure for a full, flat chart of constant thickness.
/// Preparation consumes the matching physical basis only on success.
class PgFlatGram {
public:
    PgFlatGram(PgFlatGram&&) noexcept = default;
    PgFlatGram& operator=(PgFlatGram&&) noexcept = default;

    [[nodiscard]] static core::Result<PgFlatGram> prepare(PgReferenceBasis&&,
                                                          std::uint32_t fractional_bits,
                                                          const PgFlatGramLimits&, core::Workspace&,
                                                          core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgDiscreteGramEnclosure> enclosure() const& noexcept;
    [[nodiscard]] core::Result<PgDiscreteGramEnclosure> enclosure() const&& = delete;

    /// Certify the candidate against this enclosure and its actual owned reference.
    [[nodiscard]] core::Result<PgConditioningCertificate> certify(PgConditioningCandidate,
                                                                  const PgConditioningLimits&,
                                                                  core::Workspace&,
                                                                  core::ResourcePlan&) && noexcept;

private:
    PgFlatGram(PgReferenceBasis&&, core::PlannedBuffer<PgGramEntry>, std::uint32_t) noexcept;
    PgReferenceBasis basis_;
    core::PlannedBuffer<PgGramEntry> entries_;
    std::uint32_t fractional_bits_;
};

} // namespace uvlas::structure
