#pragma once

#include <cstdint>
#include <stop_token>
#include <uvlas/io/si_case.hpp>
#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::coupling::detail {

struct CasePgReferenceLimits {
    structure::PgGraphLimits graph;
    structure::PgCoefficientLimits coefficients;
    structure::PgPartitionLimits partition;
    structure::PgReferenceLimits geometry;
    std::uint64_t maximum_input_work;
    std::stop_token cancellation;
};

/// Prepare the local reference of a Case shell whose geometry is an explicit PG.
/// Cross-patch continuity, basis qualification and constraints are separate stages.
/// The result owns its graph, dictionary, partition and physical frame.
[[nodiscard]] core::Result<structure::PgReferenceGeometry>
prepare_case_pg_reference(const io::SiCase&, std::uint64_t shell_id, const CasePgReferenceLimits&,
                          core::Workspace&, core::ResourcePlan&) noexcept;
} // namespace uvlas::coupling::detail
