#pragma once

#include "nastran_pbeam_cards.hpp"
#include "nastran_position_enclosures.hpp"

#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::io::detail {

using NastranSiBeamVector = std::array<NastranCartesianInterval, 3>;

enum class NastranBeamSiFailure : std::uint8_t {
    Source,
    Selection,
    Release,
    Warping,
    Recovery,
    Orientation,
    EndOffset,
    Axis,
    Material,
    Section,
    Modifier,
};

struct NastranBeamSiDiagnostic {
    NastranBeamSiFailure reason = NastranBeamSiFailure::Source;
    std::optional<NastranElementId> beam;
    std::optional<NastranBulkPosition> position;
    std::optional<NastranGridId> grid;
    std::optional<NastranCoordinateId> coordinate;
};

struct NastranBeamSiLimits {
    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_directed_work = 0;
    std::uint64_t maximum_stations = 0; // Includes end A.
};

struct NastranSiBeamEnd {
    NastranSiBeamVector grid_metres, offset_metres, shear_centre_metres;
};

struct NastranSiBeamSection {
    structure::ExactRational area_m2, plane1_m4, plane2_m4, product_m4, torsion_m4,
        nonstructural_kg_per_m;
};

struct NastranSiBeamStiffnessScalars {
    // The plane numbers are native PBEAM plane numbers.
    structure::ExactRational axial_n, shear_plane1_n, shear_plane2_n, torsion_nm2, bend_plane1_nm2,
        bend_plane2_nm2, bend_product_nm2;
};

struct NastranSiBeamStation {
    structure::ExactRational fraction;
    NastranPbeamSection source;
    std::optional<NastranBulkPosition> fraction_field;
    NastranSiBeamSection section;
    NastranSiBeamStiffnessScalars stiffness;
    std::optional<structure::ExactRational> line_mass_kg_per_m;
};

struct NastranSiBeamMaterial {
    NastranMat1 source;
    structure::ExactRational young_pa, shear_pa, poisson;
    std::optional<structure::ExactRational> density_kg_per_m3;
};

struct NastranSiBeamSectionOffsets {
    // Pairs of (y,z) at end A, then end B, from the shear centre.
    std::array<structure::ExactRational, 4> nonstructural_mass_centre_metres, neutral_axis_metres;
    std::array<structure::ExactRational, 2> nonstructural_rotary_kg_m;
};

struct NastranSiBeamInput {
    NastranBeamOrigin source;
    NastranCbeam beam;
    NastranPbeam property;
    SiUnits units;
    NastranSiBeamMaterial material;
    NastranSiBeamSectionOffsets section_offsets;
    NastranSiBeamEnd a, b;
    // Columns of the right-handed shear-centre frame in basic axes.
    std::array<NastranSiBeamVector, 3> axes;
    NastranCartesianInterval length_metres;
    // End A first, followed by the native declaration order, including a declared end B.
    core::PlannedBuffer<NastranSiBeamStation> stations;
    structure::PgReferenceCallInfo exact_work;
    std::uint64_t directed_work = 0;
};

// Owns point/section inputs of one physical straight beam. The products per
// native section plane are not a condensed 6x6 law or a recovered mass tensor.
[[nodiscard]] core::Result<NastranSiBeamInput>
prepare_nastran_si_beam(const NastranIsotropicShellCards&, const NastranBeamCards&,
                        const NastranPbeamCards&, NastranElementId, const SiUnits&,
                        const NastranBeamSiLimits&, core::ResourcePlan&,
                        std::optional<NastranBeamSiDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
