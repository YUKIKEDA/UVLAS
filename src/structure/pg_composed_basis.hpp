#pragma once

#include "pg_rank_storage.hpp"
#include "pg_reference_arithmetic.hpp"

#include <uvlas/structure/pg_conditioning.hpp>

namespace uvlas::structure::detail {

/// A per-call cache from one immutable conditioning certificate. It borrows that
/// certificate, whose lifetime encloses every evaluation and the cache itself.
class PgComposedBasis {
public:
    [[nodiscard]] static std::optional<PgComposedBasis> prepare(const PgConditioningCertificate&,
                                                                ReferenceArithmetic&);

    /// Nullopt leaves the point on the existing exact reference route.
    [[nodiscard]] std::optional<PgConditionedReferenceBatch>
    evaluate(const PgConditioningCertificate&, const PgReferenceVolumePoint&,
             ReferenceArithmetic&) const;

private:
    struct SourceColumn {
        std::span<const std::uint32_t> rows;
        std::span<const RankFraction<RankCell>> values;
    };

    [[nodiscard]] static SourceColumn source_column(const PgRankPreparation&, PgRankBasisColumn);

    using Units = RankInteger<512>;
    static constexpr std::size_t groups = 12;
    static constexpr std::uint64_t coefficient_bits = 512;
    static constexpr std::uint64_t coordinate_bits = 32;
    static constexpr std::size_t degree = 21;

    const PgConditioningCertificate* owner_ = nullptr;
    const PgReferenceBasis* basis_ = nullptr;
    std::span<const PgRankRow> rows_;
    core::PlannedBuffer<Units> entries_;
    std::array<Units, groups> denominators_{};
    std::size_t rank_ = 0;
};

} // namespace uvlas::structure::detail
