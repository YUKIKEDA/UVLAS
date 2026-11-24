#pragma once

#include <uvlas/structure/pg_conditioning.hpp>

namespace uvlas::structure {

/// Weight in the unit triangle and unit depth coordinates of the selected volume.
/// The exact physical Jacobian already includes the partition and layer map.
struct PgDiscreteVolumePoint {
    PgReferenceVolumePoint point;
    const ExactRational& weight;
};

/// Exact immutable snapshot supporting the finite-rule qualification.
struct PgOwnedDiscreteVolumePoint {
    core::Index volume{};
    ExactRational u, v, depth, weight;
};

/// Conditioning qualified against exactly the finite rule retained here.
class PgFixedDiscreteGramCertificate {
public:
    PgFixedDiscreteGramCertificate(PgFixedDiscreteGramCertificate&&) noexcept = default;
    PgFixedDiscreteGramCertificate& operator=(PgFixedDiscreteGramCertificate&&) noexcept = default;
    PgFixedDiscreteGramCertificate(const PgFixedDiscreteGramCertificate&) = delete;
    PgFixedDiscreteGramCertificate& operator=(const PgFixedDiscreteGramCertificate&) = delete;

    [[nodiscard]] const PgConditioningCertificate& conditioning() const& noexcept {
        return conditioning_;
    }
    [[nodiscard]] const PgConditioningCertificate& conditioning() const&& = delete;
    [[nodiscard]] std::span<const PgOwnedDiscreteVolumePoint> rule() const& noexcept {
        return rule_.values();
    }
    [[nodiscard]] std::span<const PgOwnedDiscreteVolumePoint> rule() const&& = delete;

private:
    friend class PgDiscreteGram;
    PgFixedDiscreteGramCertificate(PgConditioningCertificate&&,
                                   core::PlannedBuffer<PgOwnedDiscreteVolumePoint>) noexcept;

    PgConditioningCertificate conditioning_;
    core::PlannedBuffer<PgOwnedDiscreteVolumePoint> rule_;
};

struct PgDiscreteGramGrid {
    std::uint32_t working_fractional_bits;
    std::uint32_t output_fractional_bits;
};

struct PgDiscreteGramLimits {
    PgReferenceLimits arithmetic;
    std::uint64_t maximum_integer_steps = 0;
};

struct PgDiscreteGramInfo {
    PgReferenceCallInfo arithmetic;
    std::uint64_t integer_steps;
    std::size_t points;
};

/// Encloses the density-one inner product for an explicit finite quadrature rule.
class PgDiscreteGram {
public:
    PgDiscreteGram(PgDiscreteGram&&) noexcept = default;
    PgDiscreteGram& operator=(PgDiscreteGram&&) noexcept = default;
    PgDiscreteGram(const PgDiscreteGram&) = delete;
    PgDiscreteGram& operator=(const PgDiscreteGram&) = delete;

    /// Point records are grouped in increasing volume order, cover every
    /// owned volume, and have strictly positive exact weights.
    /// Transfer the basis only after successful preparation.
    [[nodiscard]] static core::Result<PgDiscreteGram>
    prepare(PgReferenceBasis&&, std::span<const PgDiscreteVolumePoint>, PgDiscreteGramGrid,
            const PgDiscreteGramLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgDiscreteGramEnclosure> enclosure() const& noexcept;
    [[nodiscard]] core::Result<PgDiscreteGramEnclosure> enclosure() const&& = delete;
    [[nodiscard]] std::span<const PgOwnedDiscreteVolumePoint> rule() const& noexcept {
        return rule_.values();
    }
    [[nodiscard]] std::span<const PgOwnedDiscreteVolumePoint> rule() const&& = delete;
    [[nodiscard]] PgDiscreteGramInfo info() const noexcept {
        return info_;
    }

    /// Generate the transform and check its congruence and required fields against
    /// this enclosure. Leave the basis and enclosure available on failure.
    [[nodiscard]] core::Result<PgFixedDiscreteGramCertificate>
    generate_and_certify(PgConditioningGeneration, const PgConditioningLimits&, core::Workspace&,
                         core::ResourcePlan&) && noexcept;

private:
    PgDiscreteGram(PgReferenceBasis&&, core::PlannedBuffer<PgGramEntry>, std::uint32_t,
                   core::PlannedBuffer<PgOwnedDiscreteVolumePoint>, PgDiscreteGramInfo) noexcept;

    PgReferenceBasis basis_;
    core::PlannedBuffer<PgGramEntry> entries_;
    std::uint32_t fractional_bits_;
    core::PlannedBuffer<PgOwnedDiscreteVolumePoint> rule_;
    PgDiscreteGramInfo info_;
};

} // namespace uvlas::structure
