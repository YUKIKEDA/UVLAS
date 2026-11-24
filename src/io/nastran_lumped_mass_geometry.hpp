#pragma once

#include "nastran_lumped_mass_cards.hpp"
#include "nastran_position_enclosures.hpp"

#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::io::detail {

enum class NastranBasicMassKind : std::uint8_t { Particle, Rigid };
enum class NastranMassGeometryFailure : std::uint8_t {
    Source,
    Selection,
    Frame,
    Mass,
    Inertia,
    Position,
    Grid,
    Arithmetic,
};

struct NastranMassGeometryDiagnostic {
    NastranMassGeometryFailure reason = NastranMassGeometryFailure::Source;
    std::optional<NastranElementId> element;
    std::optional<NastranBulkPosition> position;
    std::optional<NastranCoordinateId> coordinate;
};

struct NastranMassGeometryLimits {
    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_directed_work = 0;
};

struct NastranBasicMassGeometry {
    NastranLumpedMassSource source;
    NastranLumpedMass card;
    SiUnits input_units;
    NastranBasicMassKind kind = NastranBasicMassKind::Particle;
    structure::ExactRational mass_kg;
    std::array<NastranCartesianInterval, 3> grid_metres, centre_metres, grid_to_centre_metres;
    // Symmetric tensors in basic axes, lower triangle in row order (11,21,22,31,32,33).
    std::array<structure::ExactRational, 6> centre_inertia_kg_m2;
    std::array<NastranCartesianInterval, 6> grid_inertia_kg_m2;
    structure::PgReferenceCallInfo exact_work;
    std::uint64_t directed_work = 0;
};

// Converts CID zero/blank and CID -1 from an authenticated deck. The listed
// diagonal/inertia products and the tensor about the centre remain available
// independently of the single parallel-axis shift back to the GRID.
[[nodiscard]] core::Result<NastranBasicMassGeometry>
prepare_nastran_basic_mass(const NastranShellCards&, const NastranLumpedMassCards&,
                           NastranElementId, const SiUnits&, const NastranMassGeometryLimits&,
                           core::ResourcePlan&,
                           std::optional<NastranMassGeometryDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
