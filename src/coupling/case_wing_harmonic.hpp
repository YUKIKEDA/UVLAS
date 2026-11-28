#pragma once

#include "case_wing_geometry.hpp"
#include "discrete_nonpenetration.hpp"

namespace uvlas::coupling::detail {

struct CaseWingHarmonicGroup {
    std::uint64_t force_observation, moment_observation;
};

struct CaseWingHarmonicRequest {
    std::array<CaseWingHarmonicGroup, 2> observations;
    std::array<std::span<const HarmonicCirculationRhs>, 4> drives;
    double omega_radians_per_second;
    std::span<const double> age_panel_starts_seconds;
    std::uint64_t maximum_age_samples;
};

struct CaseWingHarmonicProjection {
    core::Precision representation;
    CaseFlowProjectionLimits limits;
    CaseFlowComponentBudget density;
    std::array<CaseFlowComponentBudget, 3> velocity, origin;
};

struct CaseWingHarmonicResult {
    io::SchemaIdentity source;
    std::array<CaseWingHarmonicGroup, 2> observations;
    DiscreteHarmonicCirculation solved;
};

/// Fixed-base finite harmonic verification. Uses the Case wake retention, uniform ambient
/// and convection, density, and the paired disjoint world force/moment observations.
[[nodiscard]] core::Result<CaseWingHarmonicResult>
solve_case_wing_harmonic(const io::SiCase&, const CaseWingGeometry&,
                         const aero::detail::GeometryContext&, const CaseWingHarmonicRequest&,
                         const CaseWingHarmonicProjection&, const structure::PgReferenceLimits&,
                         const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest&,
                         core::ResourcePlan&) noexcept;

} // namespace uvlas::coupling::detail
