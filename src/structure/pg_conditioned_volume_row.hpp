#pragma once

#include <uvlas/structure/pg_conditioning.hpp>

namespace uvlas::structure::detail {

/// Internal reference volume row for mass, which consumes physical measure and N.
struct PgConditionedVolumeRow {
    PgReferenceValue reference;
    core::PlannedBuffer<PgReferenceScalar> values;
    PgReferenceCallInfo work;

    [[nodiscard]] static core::Result<PgConditionedVolumeRow>
    evaluate(const PgConditioningCertificate&, const PgReferenceVolumePoint&,
             const PgReferenceLimits&, core::ResourcePlan&) noexcept;
};

} // namespace uvlas::structure::detail
