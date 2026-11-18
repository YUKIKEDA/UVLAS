#include "case_schema_storage.hpp"

#include <array>
namespace uvlas::io::detail {
namespace {
constexpr std::array<FieldDefinition, 12> fields_0{{
    {"schema", FieldType::T0, false},
    {"units", FieldType::T1, false},
    {"model", FieldType::T2, false},
    {"flow", FieldType::T3, false},
    {"analysis", FieldType::T4, false},
    {"observations", FieldType::T5, false},
    {"numerics", FieldType::T6, false},
    {"resources", FieldType::T7, false},
    {"output", FieldType::T8, false},
    {"initial", FieldType::T9, true},
    {"restart", FieldType::T10, true},
    {"provenance", FieldType::T11, true},
}};
constexpr std::array<FieldDefinition, 4> fields_1{{
    {"length", FieldType::T12, false},
    {"mass", FieldType::T13, false},
    {"time", FieldType::T14, false},
    {"angle", FieldType::T15, false},
}};
constexpr std::array<FieldDefinition, 8> fields_2{{
    {"geometry", FieldType::T16, false},
    {"materials", FieldType::T17, false},
    {"shells", FieldType::T18, false},
    {"beams", FieldType::T19, false},
    {"joints", FieldType::T20, false},
    {"supports", FieldType::T21, false},
    {"masses", FieldType::T22, false},
    {"aero_surfaces", FieldType::T23, false},
}};
constexpr std::array<FieldDefinition, 3> fields_3{{
    {"description", FieldType::T24, false},
    {"references", FieldType::T25, false},
    {"assumptions", FieldType::T25, false},
}};
constexpr std::array<FieldDefinition, 2> fields_4{{
    {"origin", FieldType::T26, false},
    {"rotation", FieldType::T27, false},
}};
constexpr std::array<FieldDefinition, 2> fields_5{{
    {"kind", FieldType::T28, false},
    {"id", FieldType::T29, false},
}};
constexpr std::array<FieldDefinition, 3> fields_6{{
    {"target", FieldType::T30, false},
    {"coordinates", FieldType::T31, false},
    {"side", FieldType::T32, false},
}};
constexpr std::array<FieldDefinition, 5> fields_7{{
    {"target", FieldType::T30, false},
    {"kind", FieldType::T33, false},
    {"geometry", FieldType::T34, false},
    {"boundary", FieldType::T35, false},
    {"side", FieldType::T32, false},
}};
constexpr std::array<FieldDefinition, 2> fields_8{{
    {"powers", FieldType::T36, false},
    {"coefficient", FieldType::T37, false},
}};
constexpr std::array<FieldDefinition, 2> fields_9{{
    {"kind", FieldType::T38, false},
    {"value", FieldType::T39, false},
}};
constexpr std::array<FieldDefinition, 6> fields_10{{
    {"kind", FieldType::T40, false},
    {"coordinates", FieldType::T41, false},
    {"origin", FieldType::T31, false},
    {"lengths", FieldType::T31, false},
    {"components", FieldType::T42, false},
    {"error", FieldType::T39, false},
}};
constexpr std::array<FieldDefinition, 2> fields_11{{
    {"kind", FieldType::T38, false},
    {"value", FieldType::T43, false},
}};
constexpr std::array<FieldDefinition, 5> fields_12{{
    {"kind", FieldType::T40, false},
    {"origin", FieldType::T44, false},
    {"scale", FieldType::T44, false},
    {"coefficients", FieldType::T45, false},
    {"domain", FieldType::T46, false},
}};
constexpr std::array<FieldDefinition, 4> fields_13{{
    {"kind", FieldType::T47, false},
    {"times", FieldType::T48, false},
    {"values", FieldType::T45, false},
    {"outside", FieldType::T49, false},
}};
constexpr std::array<FieldDefinition, 1> fields_14{{
    {"kind", FieldType::T50, false},
}};
constexpr std::array<FieldDefinition, 2> fields_15{{
    {"kind", FieldType::T51, false},
    {"period", FieldType::T44, false},
}};
constexpr std::array<FieldDefinition, 4> fields_16{{
    {"origin", FieldType::T52, false},
    {"translation", FieldType::T53, false},
    {"rotation_vector", FieldType::T53, false},
    {"events", FieldType::T48, false},
}};
constexpr std::array<FieldDefinition, 3> fields_17{{
    {"kind", FieldType::T54, false},
    {"direction", FieldType::T55, false},
    {"angle", FieldType::T56, false},
}};
constexpr std::array<FieldDefinition, 2> fields_18{{
    {"kind", FieldType::T57, false},
    {"columns", FieldType::T58, false},
}};
constexpr std::array<FieldDefinition, 15> fields_19{{
    {"kind", FieldType::T59, false},
    {"id", FieldType::T60, false},
    {"frame", FieldType::T52, false},
    {"reference_length", FieldType::T61, false},
    {"x0", FieldType::T62, false},
    {"width", FieldType::T62, false},
    {"skew", FieldType::T62, false},
    {"degree", FieldType::T63, false},
    {"coefficients", FieldType::T64, false},
    {"domain", FieldType::T65, false},
    {"thickness", FieldType::T66, false},
    {"n", FieldType::T67, false},
    {"convention", FieldType::T68, false},
    {"boundaries", FieldType::T69, false},
    {"continuity", FieldType::T70, false},
}};
constexpr std::array<FieldDefinition, 2> fields_20{{
    {"powers", FieldType::T63, false},
    {"coefficient", FieldType::T62, false},
}};
constexpr std::array<FieldDefinition, 3> fields_21{{
    {"vertices", FieldType::T65, false},
    {"triangles", FieldType::T71, false},
    {"values", FieldType::T72, false},
}};
constexpr std::array<FieldDefinition, 2> fields_22{{
    {"name", FieldType::T24, false},
    {"vertices", FieldType::T65, false},
}};
constexpr std::array<FieldDefinition, 5> fields_23{{
    {"other", FieldType::T34, false},
    {"kind", FieldType::T73, false},
    {"own_boundary", FieldType::T35, false},
    {"other_boundary", FieldType::T35, false},
    {"orientation", FieldType::T74, false},
}};
constexpr std::array<FieldDefinition, 12> fields_24{{
    {"kind", FieldType::T75, false},
    {"id", FieldType::T60, false},
    {"frame", FieldType::T52, false},
    {"reference_length", FieldType::T61, false},
    {"x0", FieldType::T62, false},
    {"width", FieldType::T62, false},
    {"skew", FieldType::T62, false},
    {"degree", FieldType::T63, false},
    {"coefficients", FieldType::T64, false},
    {"domain", FieldType::T65, false},
    {"boundaries", FieldType::T69, false},
    {"continuity", FieldType::T70, false},
}};
constexpr std::array<FieldDefinition, 7> fields_25{{
    {"kind", FieldType::T76, false},
    {"id", FieldType::T60, false},
    {"deck", FieldType::T77, false},
    {"sidecar", FieldType::T78, false},
    {"dialect", FieldType::T79, false},
    {"mode", FieldType::T80, false},
    {"excluded_cards", FieldType::T25, false},
}};
constexpr std::array<FieldDefinition, 5> fields_26{{
    {"kind", FieldType::T81, false},
    {"id", FieldType::T60, false},
    {"mesh", FieldType::T82, false},
    {"sidecar", FieldType::T83, false},
    {"reference", FieldType::T84, false},
}};
constexpr std::array<FieldDefinition, 2> fields_27{{
    {"kind", FieldType::T85, false},
    {"file", FieldType::T86, false},
}};
constexpr std::array<FieldDefinition, 5> fields_28{{
    {"kind", FieldType::T87, false},
    {"degree", FieldType::T67, false},
    {"distance_budget", FieldType::T61, false},
    {"normal_budget", FieldType::T88, false},
    {"conversion_observations", FieldType::T89, false},
}};
constexpr std::array<FieldDefinition, 7> fields_29{{
    {"kind", FieldType::T90, false},
    {"id", FieldType::T91, false},
    {"young", FieldType::T92, false},
    {"poisson", FieldType::T88, false},
    {"density", FieldType::T93, false},
    {"provenance", FieldType::T11, false},
    {"applicability", FieldType::T94, false},
}};
constexpr std::array<FieldDefinition, 8> fields_30{{
    {"kind", FieldType::T95, false},
    {"id", FieldType::T91, false},
    {"young", FieldType::T96, false},
    {"shear", FieldType::T96, false},
    {"poisson", FieldType::T97, false},
    {"density", FieldType::T93, false},
    {"provenance", FieldType::T11, false},
    {"applicability", FieldType::T94, false},
}};
constexpr std::array<FieldDefinition, 7> fields_31{{
    {"kind", FieldType::T98, false},
    {"id", FieldType::T91, false},
    {"stiffness", FieldType::T99, false},
    {"components", FieldType::T100, false},
    {"density", FieldType::T93, false},
    {"provenance", FieldType::T11, false},
    {"applicability", FieldType::T94, false},
}};
constexpr std::array<FieldDefinition, 1> fields_32{{
    {"kind", FieldType::T101, false},
}};
constexpr std::array<FieldDefinition, 4> fields_33{{
    {"kind", FieldType::T102, false},
    {"strain_lower", FieldType::T103, false},
    {"strain_upper", FieldType::T103, false},
    {"reference", FieldType::T24, false},
}};
constexpr std::array<FieldDefinition, 6> fields_34{{
    {"id", FieldType::T104, false},
    {"interval", FieldType::T105, false},
    {"material", FieldType::T106, false},
    {"axis", FieldType::T107, false},
    {"initial_green", FieldType::T108, false},
    {"viscosity", FieldType::T109, false},
}};
constexpr std::array<FieldDefinition, 5> fields_35{{
    {"id", FieldType::T110, false},
    {"geometry", FieldType::T34, false},
    {"layers", FieldType::T111, false},
    {"damping_components", FieldType::T100, false},
    {"loads", FieldType::T112, false},
}};
constexpr std::array<FieldDefinition, 8> fields_36{{
    {"id", FieldType::T113, false},
    {"nodes", FieldType::T114, false},
    {"coordinate", FieldType::T115, false},
    {"section", FieldType::T116, false},
    {"initial_strain", FieldType::T103, false},
    {"residual_energy", FieldType::T117, false},
    {"viscosity", FieldType::T118, false},
    {"loads", FieldType::T112, false},
}};
constexpr std::array<FieldDefinition, 2> fields_37{{
    {"kind", FieldType::T119, false},
    {"length", FieldType::T61, false},
}};
constexpr std::array<FieldDefinition, 3> fields_38{{
    {"kind", FieldType::T120, false},
    {"length", FieldType::T61, false},
    {"arc_scale", FieldType::T88, false},
}};
constexpr std::array<FieldDefinition, 9> fields_39{{
    {"kind", FieldType::T121, false},
    {"stiffness", FieldType::T118, false},
    {"line_density", FieldType::T122, false},
    {"center", FieldType::T26, false},
    {"inertia", FieldType::T123, false},
    {"envelope", FieldType::T124, false},
    {"warping", FieldType::T125, false},
    {"recovery", FieldType::T126, false},
    {"source", FieldType::T24, false},
}};
constexpr std::array<FieldDefinition, 7> fields_40{{
    {"kind", FieldType::T127, false},
    {"radius", FieldType::T61, false},
    {"young", FieldType::T92, false},
    {"poisson", FieldType::T88, false},
    {"density", FieldType::T93, false},
    {"effective_shear", FieldType::T128, false},
    {"warping", FieldType::T125, false},
}};
constexpr std::array<FieldDefinition, 8> fields_41{{
    {"kind", FieldType::T129, false},
    {"sides", FieldType::T124, false},
    {"young", FieldType::T92, false},
    {"poisson", FieldType::T88, false},
    {"density", FieldType::T93, false},
    {"effective_shear", FieldType::T128, false},
    {"torsion", FieldType::T130, false},
    {"warping", FieldType::T125, false},
}};
constexpr std::array<FieldDefinition, 6> fields_42{{
    {"kind", FieldType::T131, false},
    {"id", FieldType::T132, false},
    {"first", FieldType::T133, false},
    {"second", FieldType::T133, false},
    {"mapping", FieldType::T134, false},
    {"split_at_attachment", FieldType::T135, false},
}};
constexpr std::array<FieldDefinition, 6> fields_43{{
    {"first_cell", FieldType::T136, false},
    {"second_cell", FieldType::T136, false},
    {"domain", FieldType::T65, false},
    {"first", FieldType::T137, false},
    {"second", FieldType::T137, false},
    {"orientation", FieldType::T74, false},
}};
constexpr std::array<FieldDefinition, 2> fields_44{{
    {"powers", FieldType::T63, false},
    {"coefficient", FieldType::T62, false},
}};
constexpr std::array<FieldDefinition, 5> fields_45{{
    {"kind", FieldType::T138, false},
    {"id", FieldType::T132, false},
    {"attached", FieldType::T133, false},
    {"frame", FieldType::T30, false},
    {"offset", FieldType::T52, false},
}};
constexpr std::array<FieldDefinition, 5> fields_46{{
    {"kind", FieldType::T139, false},
    {"id", FieldType::T132, false},
    {"first", FieldType::T30, false},
    {"second", FieldType::T30, false},
    {"relative", FieldType::T52, false},
}};
constexpr std::array<FieldDefinition, 5> fields_47{{
    {"kind", FieldType::T140, false},
    {"id", FieldType::T141, false},
    {"target", FieldType::T142, false},
    {"components", FieldType::T143, false},
    {"motion", FieldType::T144, false},
}};
constexpr std::array<FieldDefinition, 5> fields_48{{
    {"kind", FieldType::T145, false},
    {"id", FieldType::T141, false},
    {"region", FieldType::T133, false},
    {"components", FieldType::T143, false},
    {"motion", FieldType::T144, false},
}};
constexpr std::array<FieldDefinition, 7> fields_49{{
    {"kind", FieldType::T146, false},
    {"id", FieldType::T141, false},
    {"target", FieldType::T142, false},
    {"reference", FieldType::T144, false},
    {"stiffness", FieldType::T118, false},
    {"damping", FieldType::T118, false},
    {"frame", FieldType::T52, false},
}};
constexpr std::array<FieldDefinition, 5> fields_50{{
    {"kind", FieldType::T147, false},
    {"id", FieldType::T148, false},
    {"target", FieldType::T142, false},
    {"mass", FieldType::T149, false},
    {"ledger", FieldType::T150, false},
}};
constexpr std::array<FieldDefinition, 8> fields_51{{
    {"kind", FieldType::T151, false},
    {"id", FieldType::T148, false},
    {"target", FieldType::T142, false},
    {"mass", FieldType::T149, false},
    {"center", FieldType::T26, false},
    {"inertia", FieldType::T152, false},
    {"frame", FieldType::T52, false},
    {"ledger", FieldType::T150, false},
}};
constexpr std::array<FieldDefinition, 3> fields_52{{
    {"kind", FieldType::T153, false},
    {"source", FieldType::T24, false},
    {"owner", FieldType::T30, false},
}};
constexpr std::array<FieldDefinition, 4> fields_53{{
    {"kind", FieldType::T154, false},
    {"source", FieldType::T24, false},
    {"owner", FieldType::T30, false},
    {"replaced", FieldType::T155, false},
}};
constexpr std::array<FieldDefinition, 5> fields_54{{
    {"kind", FieldType::T156, false},
    {"id", FieldType::T157, false},
    {"target", FieldType::T142, false},
    {"frame", FieldType::T158, false},
    {"series", FieldType::T53, false},
}};
constexpr std::array<FieldDefinition, 5> fields_55{{
    {"kind", FieldType::T159, false},
    {"id", FieldType::T157, false},
    {"target", FieldType::T142, false},
    {"frame", FieldType::T158, false},
    {"series", FieldType::T53, false},
}};
constexpr std::array<FieldDefinition, 3> fields_56{{
    {"kind", FieldType::T160, false},
    {"id", FieldType::T157, false},
    {"acceleration", FieldType::T53, false},
}};
constexpr std::array<FieldDefinition, 7> fields_57{{
    {"kind", FieldType::T161, false},
    {"id", FieldType::T157, false},
    {"region", FieldType::T133, false},
    {"area", FieldType::T162, false},
    {"frame", FieldType::T158, false},
    {"spatial", FieldType::T163, false},
    {"multiplier", FieldType::T164, false},
}};
constexpr std::array<FieldDefinition, 5> fields_58{{
    {"kind", FieldType::T165, false},
    {"id", FieldType::T157, false},
    {"region", FieldType::T133, false},
    {"spatial", FieldType::T166, false},
    {"multiplier", FieldType::T164, false},
}};
constexpr std::array<FieldDefinition, 5> fields_59{{
    {"kind", FieldType::T167, false},
    {"id", FieldType::T157, false},
    {"region", FieldType::T133, false},
    {"spatial", FieldType::T168, false},
    {"multiplier", FieldType::T164, false},
}};
constexpr std::array<FieldDefinition, 9> fields_60{{
    {"id", FieldType::T169, false},
    {"geometry", FieldType::T34, false},
    {"attachment", FieldType::T170, false},
    {"chord", FieldType::T171, false},
    {"span", FieldType::T171, false},
    {"leading", FieldType::T35, false},
    {"trailing", FieldType::T35, false},
    {"orientation", FieldType::T172, false},
    {"sharing", FieldType::T173, false},
}};
constexpr std::array<FieldDefinition, 1> fields_61{{
    {"kind", FieldType::T174, false},
}};
constexpr std::array<FieldDefinition, 2> fields_62{{
    {"kind", FieldType::T175, false},
    {"motion", FieldType::T144, false},
}};
constexpr std::array<FieldDefinition, 2> fields_63{{
    {"kind", FieldType::T176, false},
    {"host", FieldType::T177, false},
}};
constexpr std::array<FieldDefinition, 4> fields_64{{
    {"other", FieldType::T178, false},
    {"own_boundary", FieldType::T35, false},
    {"other_boundary", FieldType::T35, false},
    {"orientation", FieldType::T74, false},
}};
constexpr std::array<FieldDefinition, 1> fields_65{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 6> fields_66{{
    {"kind", FieldType::T180, false},
    {"density", FieldType::T93, false},
    {"uniform", FieldType::T53, false},
    {"gust", FieldType::T181, false},
    {"reference", FieldType::T182, false},
    {"events", FieldType::T48, false},
}};
constexpr std::array<FieldDefinition, 1> fields_67{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 8> fields_68{{
    {"kind", FieldType::T183, false},
    {"amplitude", FieldType::T164, false},
    {"direction", FieldType::T31, false},
    {"wavevector", FieldType::T184, false},
    {"angular_frequency", FieldType::T185, false},
    {"phase", FieldType::T186, false},
    {"origin", FieldType::T26, false},
    {"time_origin", FieldType::T44, false},
}};
constexpr std::array<FieldDefinition, 6> fields_69{{
    {"density", FieldType::T93, false},
    {"speed", FieldType::T187, false},
    {"area", FieldType::T188, false},
    {"chord", FieldType::T61, false},
    {"span", FieldType::T61, false},
    {"moment_origin", FieldType::T26, false},
}};
constexpr std::array<FieldDefinition, 1> fields_70{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 6> fields_71{{
    {"kind", FieldType::T189, false},
    {"convection", FieldType::T190, false},
    {"retention", FieldType::T191, false},
    {"core_radius", FieldType::T61, false},
    {"model", FieldType::T192, false},
    {"mapping", FieldType::T193, false},
}};
constexpr std::array<FieldDefinition, 1> fields_72{{
    {"kind", FieldType::T125, false},
}};
constexpr std::array<FieldDefinition, 2> fields_73{{
    {"kind", FieldType::T194, false},
    {"velocity", FieldType::T53, false},
}};
constexpr std::array<FieldDefinition, 1> fields_74{{
    {"kind", FieldType::T195, false},
}};
constexpr std::array<FieldDefinition, 4> fields_75{{
    {"kind", FieldType::T196, false},
    {"duration", FieldType::T44, false},
    {"distance", FieldType::T61, false},
    {"error_budget", FieldType::T88, false},
}};
constexpr std::array<FieldDefinition, 6> fields_76{{
    {"kind", FieldType::T197, false},
    {"physics", FieldType::T198, false},
    {"branch", FieldType::T199, false},
    {"speed", FieldType::T200, false},
    {"wake", FieldType::T201, false},
    {"loads", FieldType::T112, false},
}};
constexpr std::array<FieldDefinition, 6> fields_77{{
    {"kind", FieldType::T202, false},
    {"physics", FieldType::T198, false},
    {"branch", FieldType::T199, false},
    {"intervals", FieldType::T203, false},
    {"wake", FieldType::T201, false},
    {"loads", FieldType::T112, false},
}};
constexpr std::array<FieldDefinition, 8> fields_78{{
    {"kind", FieldType::T204, false},
    {"physics", FieldType::T205, false},
    {"branches", FieldType::T206, false},
    {"intervals", FieldType::T207, false},
    {"step", FieldType::T44, false},
    {"wake", FieldType::T201, false},
    {"search", FieldType::T208, false},
    {"loads", FieldType::T112, false},
}};
constexpr std::array<FieldDefinition, 1> fields_79{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 2> fields_80{{
    {"kind", FieldType::T209, false},
    {"value", FieldType::T187, false},
}};
constexpr std::array<FieldDefinition, 3> fields_81{{
    {"id", FieldType::T210, false},
    {"seed", FieldType::T211, false},
    {"continuation", FieldType::T212, false},
}};
constexpr std::array<FieldDefinition, 1> fields_82{{
    {"kind", FieldType::T213, false},
}};
constexpr std::array<FieldDefinition, 2> fields_83{{
    {"kind", FieldType::T214, false},
    {"history", FieldType::T215, false},
}};
constexpr std::array<FieldDefinition, 1> fields_84{{
    {"kind", FieldType::T216, false},
}};
constexpr std::array<FieldDefinition, 7> fields_85{{
    {"kind", FieldType::T217, false},
    {"initial_increment", FieldType::T187, false},
    {"minimum_increment", FieldType::T187, false},
    {"maximum_increment", FieldType::T187, false},
    {"arclength", FieldType::T218, false},
    {"maximum_points", FieldType::T219, false},
    {"maximum_retries", FieldType::T67, false},
}};
constexpr std::array<FieldDefinition, 4> fields_86{{
    {"state_scales", FieldType::T220, false},
    {"speed_scale", FieldType::T187, false},
    {"orientation", FieldType::T221, false},
    {"enable_at_limit", FieldType::T135, false},
}};
constexpr std::array<FieldDefinition, 4> fields_87{{
    {"start", FieldType::T44, false},
    {"end", FieldType::T44, false},
    {"step", FieldType::T44, false},
    {"restart", FieldType::T222, false},
}};
constexpr std::array<FieldDefinition, 2> fields_88{{
    {"lower", FieldType::T187, false},
    {"upper", FieldType::T187, false},
}};
constexpr std::array<FieldDefinition, 8> fields_89{{
    {"queue", FieldType::T223, false},
    {"speed_absolute", FieldType::T187, false},
    {"speed_relative", FieldType::T88, false},
    {"growth_budget", FieldType::T224, false},
    {"mode_association", FieldType::T225, false},
    {"maximum_points", FieldType::T219, false},
    {"maximum_retries", FieldType::T67, false},
    {"coverage", FieldType::T226, false},
}};
constexpr std::array<FieldDefinition, 5> fields_90{{
    {"physical_observations", FieldType::T89, false},
    {"mac_minimum", FieldType::T88, false},
    {"subspace_angle_maximum", FieldType::T186, false},
    {"eigenvalue_absolute", FieldType::T224, false},
    {"eigenvalue_relative", FieldType::T88, false},
}};
constexpr std::array<FieldDefinition, 3> fields_91{{
    {"branches", FieldType::T227, false},
    {"mode_observations", FieldType::T89, false},
    {"require_physical", FieldType::T135, false},
}};
constexpr std::array<FieldDefinition, 3> fields_92{{
    {"kind", FieldType::T228, false},
    {"branch", FieldType::T229, false},
    {"perturbation", FieldType::T230, false},
}};
constexpr std::array<FieldDefinition, 4> fields_93{{
    {"kind", FieldType::T231, false},
    {"state", FieldType::T232, false},
    {"circulation", FieldType::T233, false},
    {"wake", FieldType::T234, false},
}};
constexpr std::array<FieldDefinition, 2> fields_94{{
    {"kind", FieldType::T235, false},
    {"file", FieldType::T236, false},
}};
constexpr std::array<FieldDefinition, 1> fields_95{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 3> fields_96{{
    {"shells", FieldType::T237, false},
    {"beams", FieldType::T238, false},
    {"rigids", FieldType::T239, false},
}};
constexpr std::array<FieldDefinition, 5> fields_97{{
    {"shell", FieldType::T177, false},
    {"reference", FieldType::T240, false},
    {"representation", FieldType::T241, false},
    {"position", FieldType::T242, false},
    {"velocity", FieldType::T242, false},
}};
constexpr std::array<FieldDefinition, 2> fields_98{{
    {"kind", FieldType::T243, false},
    {"high", FieldType::T171, false},
}};
constexpr std::array<FieldDefinition, 3> fields_99{{
    {"kind", FieldType::T244, false},
    {"high", FieldType::T171, false},
    {"low", FieldType::T171, false},
}};
constexpr std::array<FieldDefinition, 2> fields_100{{
    {"beam", FieldType::T245, false},
    {"nodes", FieldType::T246, false},
}};
constexpr std::array<FieldDefinition, 2> fields_101{{
    {"rigid", FieldType::T247, false},
    {"state", FieldType::T248, false},
}};
constexpr std::array<FieldDefinition, 4> fields_102{{
    {"position", FieldType::T26, false},
    {"rotation", FieldType::T27, false},
    {"linear_velocity", FieldType::T249, false},
    {"angular_velocity", FieldType::T250, false},
}};
constexpr std::array<FieldDefinition, 1> fields_103{{
    {"surfaces", FieldType::T251, false},
}};
constexpr std::array<FieldDefinition, 8> fields_104{{
    {"surface", FieldType::T178, false},
    {"rows", FieldType::T67, false},
    {"points", FieldType::T252, false},
    {"coefficients", FieldType::T233, false},
    {"row_ids", FieldType::T253, false},
    {"point_ids", FieldType::T254, false},
    {"band_ids", FieldType::T255, false},
    {"birth_times", FieldType::T48, false},
}};
constexpr std::array<FieldDefinition, 4> fields_105{{
    {"directory", FieldType::T24, false},
    {"manifest_sha256", FieldType::T240, false},
    {"compatibility", FieldType::T256, false},
    {"budget", FieldType::T257, false},
}};
constexpr std::array<FieldDefinition, 7> fields_106{{
    {"kind", FieldType::T258, false},
    {"id", FieldType::T259, false},
    {"point", FieldType::T142, false},
    {"quantity", FieldType::T260, false},
    {"basis", FieldType::T261, false},
    {"comparison", FieldType::T262, false},
    {"budget", FieldType::T263, false},
}};
constexpr std::array<FieldDefinition, 10> fields_107{{
    {"kind", FieldType::T264, false},
    {"id", FieldType::T259, false},
    {"point", FieldType::T142, false},
    {"layer", FieldType::T265, false},
    {"interface", FieldType::T266, false},
    {"measure", FieldType::T267, false},
    {"contribution", FieldType::T268, false},
    {"basis", FieldType::T261, false},
    {"comparison", FieldType::T262, false},
    {"budget", FieldType::T263, false},
}};
constexpr std::array<FieldDefinition, 8> fields_108{{
    {"kind", FieldType::T269, false},
    {"id", FieldType::T259, false},
    {"region", FieldType::T133, false},
    {"origin", FieldType::T26, false},
    {"basis", FieldType::T261, false},
    {"quantity", FieldType::T270, false},
    {"comparison", FieldType::T262, false},
    {"budget", FieldType::T263, false},
}};
constexpr std::array<FieldDefinition, 7> fields_109{{
    {"kind", FieldType::T271, false},
    {"id", FieldType::T259, false},
    {"support", FieldType::T272, false},
    {"origin", FieldType::T26, false},
    {"basis", FieldType::T261, false},
    {"comparison", FieldType::T262, false},
    {"budget", FieldType::T263, false},
}};
constexpr std::array<FieldDefinition, 8> fields_110{{
    {"kind", FieldType::T273, false},
    {"id", FieldType::T259, false},
    {"surfaces", FieldType::T274, false},
    {"origin", FieldType::T26, false},
    {"basis", FieldType::T261, false},
    {"quantity", FieldType::T275, false},
    {"comparison", FieldType::T262, false},
    {"budget", FieldType::T263, false},
}};
constexpr std::array<FieldDefinition, 8> fields_111{{
    {"kind", FieldType::T276, false},
    {"id", FieldType::T259, false},
    {"targets", FieldType::T277, false},
    {"points", FieldType::T278, false},
    {"basis", FieldType::T261, false},
    {"quantity", FieldType::T279, false},
    {"comparison", FieldType::T262, false},
    {"budget", FieldType::T263, false},
}};
constexpr std::array<FieldDefinition, 1> fields_112{{
    {"kind", FieldType::T280, false},
}};
constexpr std::array<FieldDefinition, 2> fields_113{{
    {"kind", FieldType::T57, false},
    {"rotation", FieldType::T27, false},
}};
constexpr std::array<FieldDefinition, 4> fields_114{{
    {"window", FieldType::T46, false},
    {"statistic", FieldType::T281, false},
    {"side", FieldType::T282, false},
    {"reference", FieldType::T283, false},
}};
constexpr std::array<FieldDefinition, 1> fields_115{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 4> fields_116{{
    {"kind", FieldType::T284, false},
    {"source", FieldType::T285, false},
    {"columns", FieldType::T25, false},
    {"uncertainty", FieldType::T171, false},
}};
constexpr std::array<FieldDefinition, 4> fields_117{{
    {"absolute", FieldType::T171, false},
    {"relative", FieldType::T171, false},
    {"scale", FieldType::T171, false},
    {"groups", FieldType::T286, false},
}};
constexpr std::array<FieldDefinition, 3> fields_118{{
    {"status", FieldType::T287, false},
    {"kind", FieldType::T288, false},
    {"reason", FieldType::T24, false},
}};
constexpr std::array<FieldDefinition, 3> fields_119{{
    {"status", FieldType::T289, false},
    {"kind", FieldType::T288, false},
    {"representative", FieldType::T288, false},
}};
constexpr std::array<FieldDefinition, 5> fields_120{{
    {"status", FieldType::T290, false},
    {"kind", FieldType::T288, false},
    {"covers", FieldType::T291, false},
    {"absolute", FieldType::T171, false},
    {"relative", FieldType::T171, false},
}};
constexpr std::array<FieldDefinition, 9> fields_121{{
    {"patch", FieldType::T292, false},
    {"quadrature", FieldType::T293, false},
    {"linear", FieldType::T294, false},
    {"newton", FieldType::T295, false},
    {"precision", FieldType::T296, false},
    {"eigen", FieldType::T297, false},
    {"wake_rows", FieldType::T298, false},
    {"scales", FieldType::T220, false},
    {"error_groups", FieldType::T286, false},
}};
constexpr std::array<FieldDefinition, 1> fields_122{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 2> fields_123{{
    {"kind", FieldType::T299, false},
    {"maximum_rows", FieldType::T219, false},
}};
constexpr std::array<FieldDefinition, 1> fields_124{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 6> fields_125{{
    {"kind", FieldType::T189, false},
    {"maximum_degree", FieldType::T67, false},
    {"maximum_partition", FieldType::T219, false},
    {"geometry_error", FieldType::T61, false},
    {"trace_error", FieldType::T61, false},
    {"rank", FieldType::T300, false},
}};
constexpr std::array<FieldDefinition, 2> fields_126{{
    {"structure", FieldType::T301, false},
    {"aero", FieldType::T302, false},
}};
constexpr std::array<FieldDefinition, 1> fields_127{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 8> fields_128{{
    {"kind", FieldType::T189, false},
    {"surface_degree", FieldType::T219, false},
    {"thickness_points", FieldType::T219, false},
    {"independent_surface_degree", FieldType::T219, false},
    {"independent_thickness_points", FieldType::T219, false},
    {"absolute", FieldType::T88, false},
    {"relative", FieldType::T88, false},
    {"maximum_subdivisions", FieldType::T67, false},
}};
constexpr std::array<FieldDefinition, 1> fields_129{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 3> fields_130{{
    {"kind", FieldType::T303, false},
    {"absolute", FieldType::T88, false},
    {"relative", FieldType::T88, false},
}};
constexpr std::array<FieldDefinition, 4> fields_131{{
    {"kind", FieldType::T304, false},
    {"field", FieldType::T305, false},
    {"targets", FieldType::T277, false},
    {"values", FieldType::T171, false},
}};
constexpr std::array<FieldDefinition, 3> fields_132{{
    {"kind", FieldType::T306, false},
    {"field", FieldType::T307, false},
    {"constraint", FieldType::T308, false},
}};
constexpr std::array<FieldDefinition, 1> fields_133{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 6> fields_134{{
    {"kind", FieldType::T309, false},
    {"position", FieldType::T310, false},
    {"orientation", FieldType::T310, false},
    {"multiplier_work_scale", FieldType::T311, false},
    {"reaction_force", FieldType::T312, false},
    {"reaction_moment", FieldType::T313, false},
}};
constexpr std::array<FieldDefinition, 1> fields_135{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 2> fields_136{{
    {"kind", FieldType::T314, false},
    {"values", FieldType::T31, false},
}};
constexpr std::array<FieldDefinition, 6> fields_137{{
    {"kind", FieldType::T315, false},
    {"field", FieldType::T305, false},
    {"targets", FieldType::T277, false},
    {"layout", FieldType::T240, false},
    {"dimensions", FieldType::T316, false},
    {"values", FieldType::T171, false},
}};
constexpr std::array<FieldDefinition, 4> fields_138{{
    {"backend", FieldType::T317, false},
    {"accuracy", FieldType::T318, false},
    {"maximum_refinements", FieldType::T67, false},
    {"backend_memory", FieldType::T319, false},
}};
constexpr std::array<FieldDefinition, 9> fields_139{{
    {"backend", FieldType::T320, false},
    {"accuracy", FieldType::T318, false},
    {"maximum_iterations", FieldType::T219, false},
    {"restart", FieldType::T219, false},
    {"maximum_applications", FieldType::T67, false},
    {"maximum_adjoint", FieldType::T67, false},
    {"maximum_preconditioner", FieldType::T67, false},
    {"maximum_true_residual_checks", FieldType::T67, false},
    {"backend_memory", FieldType::T319, false},
}};
constexpr std::array<FieldDefinition, 4> fields_140{{
    {"absolute", FieldType::T88, false},
    {"relative", FieldType::T88, false},
    {"backward", FieldType::T88, false},
    {"norm", FieldType::T321, false},
}};
constexpr std::array<FieldDefinition, 6> fields_141{{
    {"maximum_iterations", FieldType::T219, false},
    {"maximum_line_searches", FieldType::T219, false},
    {"absolute", FieldType::T88, false},
    {"relative", FieldType::T88, false},
    {"correction", FieldType::T88, false},
    {"minimum_step", FieldType::T88, false},
}};
constexpr std::array<FieldDefinition, 5> fields_142{{
    {"input", FieldType::T322, false},
    {"retained", FieldType::T323, false},
    {"returned", FieldType::T324, false},
    {"projection_absolute", FieldType::T88, false},
    {"projection_relative", FieldType::T88, false},
}};
constexpr std::array<FieldDefinition, 2> fields_143{{
    {"admission", FieldType::T325, false},
    {"maximum_bytes", FieldType::T67, false},
}};
constexpr std::array<FieldDefinition, 1> fields_144{{
    {"kind", FieldType::T179, false},
}};
constexpr std::array<FieldDefinition, 6> fields_145{{
    {"kind", FieldType::T326, false},
    {"target", FieldType::T327, false},
    {"eigenvalue_scale", FieldType::T328, false},
    {"packet", FieldType::T329, false},
    {"requirements", FieldType::T330, false},
    {"limits", FieldType::T331, false},
}};
constexpr std::array<FieldDefinition, 5> fields_146{{
    {"kind", FieldType::T332, false},
    {"target", FieldType::T333, false},
    {"eigenvalue_scale", FieldType::T334, false},
    {"requirements", FieldType::T330, false},
    {"limits", FieldType::T331, false},
}};
constexpr std::array<FieldDefinition, 11> fields_147{{
    {"kind", FieldType::T335, false},
    {"target", FieldType::T327, false},
    {"eigenvalue_scale", FieldType::T328, false},
    {"packet", FieldType::T336, false},
    {"selector", FieldType::T337, false},
    {"nev", FieldType::T219, false},
    {"ncv", FieldType::T219, false},
    {"seed", FieldType::T338, false},
    {"backend_tolerance", FieldType::T88, false},
    {"requirements", FieldType::T330, false},
    {"limits", FieldType::T331, false},
}};
constexpr std::array<FieldDefinition, 2> fields_148{{
    {"real", FieldType::T328, false},
    {"imaginary", FieldType::T328, false},
}};
constexpr std::array<FieldDefinition, 2> fields_149{{
    {"kind", FieldType::T339, false},
    {"value", FieldType::T340, false},
}};
constexpr std::array<FieldDefinition, 3> fields_150{{
    {"kind", FieldType::T315, false},
    {"layout", FieldType::T240, false},
    {"values", FieldType::T341, false},
}};
constexpr std::array<FieldDefinition, 1> fields_151{{
    {"kind", FieldType::T342, false},
}};
constexpr std::array<FieldDefinition, 2> fields_152{{
    {"kind", FieldType::T343, false},
    {"shift", FieldType::T340, false},
}};
constexpr std::array<FieldDefinition, 12> fields_153{{
    {"residual", FieldType::T88, false},
    {"residual_observation", FieldType::T344, false},
    {"normalization", FieldType::T88, false},
    {"normalization_observation", FieldType::T344, false},
    {"orthogonality", FieldType::T88, false},
    {"orthogonality_observation", FieldType::T344, false},
    {"backward", FieldType::T345, false},
    {"association_absolute", FieldType::T88, false},
    {"association_relative", FieldType::T88, false},
    {"left", FieldType::T135, false},
    {"clusters", FieldType::T346, false},
    {"physical_coverage", FieldType::T135, false},
}};
constexpr std::array<FieldDefinition, 1> fields_154{{
    {"kind", FieldType::T347, false},
}};
constexpr std::array<FieldDefinition, 3> fields_155{{
    {"kind", FieldType::T290, false},
    {"budget", FieldType::T88, false},
    {"observation", FieldType::T344, false},
}};
constexpr std::array<FieldDefinition, 2> fields_156{{
    {"members", FieldType::T143, false},
    {"expected_dimension", FieldType::T219, false},
}};
constexpr std::array<FieldDefinition, 9> fields_157{{
    {"maximum_order", FieldType::T67, false},
    {"maximum_results", FieldType::T67, false},
    {"maximum_ncv", FieldType::T67, false},
    {"maximum_backend_iterations", FieldType::T67, false},
    {"maximum_attempts", FieldType::T67, false},
    {"maximum_forward", FieldType::T67, false},
    {"maximum_adjoint", FieldType::T67, false},
    {"maximum_inverse", FieldType::T67, false},
    {"backend_memory", FieldType::T319, false},
}};
constexpr std::array<FieldDefinition, 8> fields_158{{
    {"profile", FieldType::T348, false},
    {"threads", FieldType::T219, false},
    {"maximum_memory_bytes", FieldType::T67, false},
    {"maximum_wall_seconds", FieldType::T44, false},
    {"maximum_steps", FieldType::T67, false},
    {"maximum_total_iterations", FieldType::T67, false},
    {"checkpoint_temporary_bytes", FieldType::T67, false},
    {"input", FieldType::T349, false},
}};
constexpr std::array<FieldDefinition, 19> fields_159{{
    {"maximum_source_bytes", FieldType::T67, false},
    {"maximum_depth", FieldType::T67, false},
    {"maximum_nodes", FieldType::T67, false},
    {"maximum_members", FieldType::T67, false},
    {"maximum_decoded_string_bytes", FieldType::T67, false},
    {"maximum_decimal_token_bytes", FieldType::T67, false},
    {"maximum_total_text_bytes", FieldType::T67, false},
    {"maximum_typed_bytes", FieldType::T67, false},
    {"maximum_entities", FieldType::T67, false},
    {"maximum_references", FieldType::T67, false},
    {"maximum_scalars", FieldType::T67, false},
    {"maximum_polynomial_terms", FieldType::T67, false},
    {"maximum_table_samples", FieldType::T67, false},
    {"maximum_expression_nodes", FieldType::T67, false},
    {"maximum_expression_edges", FieldType::T67, false},
    {"maximum_expression_depth", FieldType::T67, false},
    {"maximum_expression_power", FieldType::T67, false},
    {"maximum_expression_points", FieldType::T67, false},
    {"maximum_diagnostic_path_bytes", FieldType::T67, false},
}};
constexpr std::array<FieldDefinition, 6> fields_160{{
    {"directory", FieldType::T24, false},
    {"observations", FieldType::T350, false},
    {"visualization", FieldType::T351, false},
    {"checkpoint", FieldType::T352, false},
    {"overwrite", FieldType::T353, false},
    {"resume", FieldType::T354, false},
}};
constexpr std::array<FieldDefinition, 2> fields_161{{
    {"every_steps", FieldType::T219, false},
    {"observations", FieldType::T89, false},
}};
constexpr std::array<FieldDefinition, 6> fields_162{{
    {"enabled", FieldType::T135, false},
    {"every_steps", FieldType::T219, false},
    {"layers", FieldType::T355, false},
    {"interfaces", FieldType::T356, false},
    {"reference_geometry", FieldType::T135, false},
    {"current_geometry", FieldType::T135, false},
}};
constexpr std::array<FieldDefinition, 1> fields_163{{
    {"kind", FieldType::T357, false},
}};
constexpr std::array<FieldDefinition, 4> fields_164{{
    {"kind", FieldType::T358, false},
    {"every_steps", FieldType::T219, false},
    {"directory", FieldType::T24, false},
    {"retain", FieldType::T219, false},
}};
constexpr std::array<FieldDefinition, 3> fields_165{{
    {"kind", FieldType::T359, false},
    {"expression", FieldType::T360, false},
    {"error", FieldType::T39, false},
}};
constexpr std::array<FieldDefinition, 2> fields_166{{
    {"kind", FieldType::T359, false},
    {"expression", FieldType::T360, false},
}};
constexpr std::array<FieldDefinition, 4> fields_167{{
    {"kind", FieldType::T361, false},
    {"id", FieldType::T157, false},
    {"region", FieldType::T133, false},
    {"expression", FieldType::T362, false},
}};
constexpr std::array<FieldDefinition, 3> fields_168{{
    {"domain", FieldType::T363, false},
    {"nodes", FieldType::T364, false},
    {"outputs", FieldType::T365, false},
}};
constexpr std::array<FieldDefinition, 4> fields_169{{
    {"kind", FieldType::T213, false},
    {"coordinates", FieldType::T41, false},
    {"lower", FieldType::T31, false},
    {"upper", FieldType::T31, false},
}};
constexpr std::array<FieldDefinition, 4> fields_170{{
    {"kind", FieldType::T366, false},
    {"lower", FieldType::T44, false},
    {"upper", FieldType::T44, false},
    {"events", FieldType::T48, false},
}};
constexpr std::array<FieldDefinition, 6> fields_171{{
    {"kind", FieldType::T367, false},
    {"position_lower", FieldType::T26, false},
    {"position_upper", FieldType::T26, false},
    {"time_lower", FieldType::T44, false},
    {"time_upper", FieldType::T44, false},
    {"events", FieldType::T48, false},
}};
constexpr std::array<FieldDefinition, 3> fields_172{{
    {"op", FieldType::T38, false},
    {"value", FieldType::T328, false},
    {"dimensions", FieldType::T368, false},
}};
constexpr std::array<FieldDefinition, 2> fields_173{{
    {"op", FieldType::T369, false},
    {"axis", FieldType::T35, false},
}};
constexpr std::array<FieldDefinition, 3> fields_174{{
    {"op", FieldType::T370, false},
    {"first", FieldType::T35, false},
    {"second", FieldType::T35, false},
}};
constexpr std::array<FieldDefinition, 2> fields_175{{
    {"op", FieldType::T371, false},
    {"operand", FieldType::T35, false},
}};
constexpr std::array<FieldDefinition, 3> fields_176{{
    {"op", FieldType::T372, false},
    {"operand", FieldType::T35, false},
    {"exponent", FieldType::T373, false},
}};
constexpr std::array<FieldDefinition, 2> fields_177{{
    {"op", FieldType::T374, false},
    {"operand", FieldType::T35, false},
}};
constexpr std::array<FieldDefinition, 6> fields_178{{
    {"kind", FieldType::T309, false},
    {"shells", FieldType::T375, false},
    {"beams", FieldType::T238, false},
    {"rigids", FieldType::T239, false},
    {"displacements", FieldType::T376, false},
    {"velocities", FieldType::T376, false},
}};
constexpr std::array<FieldDefinition, 3> fields_179{{
    {"shell", FieldType::T177, false},
    {"displacement", FieldType::T377, false},
    {"velocity", FieldType::T377, false},
}};
constexpr std::array<FieldDefinition, 4> fields_180{{
    {"kind", FieldType::T378, false},
    {"r", FieldType::T379, false},
    {"p", FieldType::T379, false},
    {"c", FieldType::T379, false},
}};
constexpr std::array<FieldDefinition, 2> fields_181{{
    {"kind", FieldType::T380, false},
    {"cells", FieldType::T381, false},
}};
constexpr std::array<FieldDefinition, 4> fields_182{{
    {"cell", FieldType::T136, false},
    {"r", FieldType::T379, false},
    {"p", FieldType::T379, false},
    {"c", FieldType::T379, false},
}};
constexpr std::array<FieldDefinition, 1> fields_183{{
    {"terms", FieldType::T382, false},
}};
constexpr std::array<FieldDefinition, 2> fields_184{{
    {"powers", FieldType::T63, false},
    {"value", FieldType::T383, false},
}};
constexpr std::array<FieldDefinition, 3> fields_185{{
    {"point", FieldType::T384, false},
    {"value", FieldType::T31, false},
    {"tolerance", FieldType::T31, false},
}};
constexpr std::array<FieldDefinition, 4> fields_186{{
    {"target", FieldType::T385, false},
    {"coordinates", FieldType::T31, false},
    {"side_s", FieldType::T32, false},
    {"side_t", FieldType::T32, false},
}};
constexpr std::array<FieldDefinition, 2> fields_187{{
    {"target", FieldType::T386, false},
    {"coordinates", FieldType::T31, false},
}};
constexpr std::array<FieldDefinition, 2> fields_188{{
    {"target", FieldType::T387, false},
    {"coordinates", FieldType::T31, false},
}};
constexpr std::array<FieldDefinition, 2> fields_189{{
    {"numerator", FieldType::T24, false},
    {"denominator", FieldType::T24, false},
}};
constexpr std::array<FieldDefinition, 4> fields_190{{
    {"source_id", FieldType::T388, false},
    {"path", FieldType::T24, false},
    {"sha256", FieldType::T240, false},
    {"role", FieldType::T446, false},
}};
const std::array<RecordDefinition, 191> records{{
    {"Case", "Case", fields_0, "", ""},
    {"Units", "Units", fields_1, "", ""},
    {"Model", "Model", fields_2, "", ""},
    {"Provenance", "Provenance", fields_3, "", ""},
    {"Frame", "Frame", fields_4, "", ""},
    {"Target", "Target", fields_5, "kind", "shell,beam,rigid"},
    {"MaterialPoint", "MaterialPoint", fields_6, "", ""},
    {"Region", "Region", fields_7, "kind", "volume,face,edge,point"},
    {"PolyTerm", "PolyTerm", fields_8, "", ""},
    {"ReferenceFieldConstant", "ReferenceField", fields_9, "kind", "constant"},
    {"ReferenceFieldPolynomial", "ReferenceField", fields_10, "kind", "polynomial"},
    {"TimeSeriesConstant", "TimeSeries", fields_11, "kind", "constant"},
    {"TimeSeriesPolynomial", "TimeSeries", fields_12, "kind", "polynomial"},
    {"TimeSeriesPchip", "TimeSeries", fields_13, "kind", "pchip"},
    {"OutsideRejectHold", "Outside", fields_14, "kind", "reject,hold"},
    {"OutsidePeriodic", "Outside", fields_15, "kind", "periodic"},
    {"Motion", "Motion", fields_16, "", ""},
    {"AxisProjected", "Axis", fields_17, "kind", "projected"},
    {"AxisExplicit", "Axis", fields_18, "kind", "explicit"},
    {"GeometryPg", "Geometry", fields_19, "kind", "pg"},
    {"PgTerm", "PgTerm", fields_20, "", ""},
    {"Thickness", "Thickness", fields_21, "", ""},
    {"Boundary", "Boundary", fields_22, "", ""},
    {"Continuity", "Continuity", fields_23, "kind", "smooth,fold"},
    {"GeometrySurfaceGraph", "Geometry", fields_24, "kind", "surface_graph"},
    {"GeometryNastran", "Geometry", fields_25, "kind", "nastran"},
    {"GeometryGmsh", "Geometry", fields_26, "kind", "gmsh"},
    {"ReferenceGeometrySupplied", "ReferenceGeometry", fields_27, "kind", "supplied"},
    {"ReferenceGeometryFit", "ReferenceGeometry", fields_28, "kind", "fit"},
    {"MaterialIsotropic", "Material", fields_29, "kind", "isotropic"},
    {"MaterialOrthotropic", "Material", fields_30, "kind", "orthotropic3d"},
    {"MaterialEffective", "Material", fields_31, "kind", "effective3d"},
    {"ApplicabilityUnevaluated", "Applicability", fields_32, "kind", "unevaluated"},
    {"ApplicabilityDeclared", "Applicability", fields_33, "kind", "declared"},
    {"Layer", "Layer", fields_34, "", ""},
    {"Shell", "Shell", fields_35, "", ""},
    {"Beam", "Beam", fields_36, "", ""},
    {"BeamCoordinateArclength", "BeamCoordinate", fields_37, "kind", "arclength"},
    {"BeamCoordinateScaled", "BeamCoordinate", fields_38, "kind", "scaled"},
    {"SectionMatrix", "Section", fields_39, "kind", "matrix"},
    {"SectionCircle", "Section", fields_40, "kind", "circle"},
    {"SectionRectangle", "Section", fields_41, "kind", "rectangle"},
    {"JointAdhesive", "Joint", fields_42, "kind", "adhesive"},
    {"TraceMap", "TraceMap", fields_43, "", ""},
    {"MapTerm", "MapTerm", fields_44, "", ""},
    {"JointRigidAttachment", "Joint", fields_45, "kind", "rigid_attachment"},
    {"JointWeld", "Joint", fields_46, "kind", "weld"},
    {"SupportConstraint", "Support", fields_47, "kind", "constraint"},
    {"SupportTraceConstraint", "Support", fields_48, "kind", "trace_constraint"},
    {"SupportElastic", "Support", fields_49, "kind", "elastic"},
    {"MassParticle", "Mass", fields_50, "kind", "particle"},
    {"MassRigid", "Mass", fields_51, "kind", "rigid"},
    {"LedgerAdd", "Ledger", fields_52, "kind", "add"},
    {"LedgerReplace", "Ledger", fields_53, "kind", "replace"},
    {"LoadPointForce", "Load", fields_54, "kind", "point_force"},
    {"LoadCouple", "Load", fields_55, "kind", "couple"},
    {"LoadGravity", "Load", fields_56, "kind", "gravity"},
    {"LoadTraction", "Load", fields_57, "kind", "traction"},
    {"LoadPressure", "Load", fields_58, "kind", "pressure"},
    {"LoadBody", "Load", fields_59, "kind", "body_force"},
    {"AeroSurface", "AeroSurface", fields_60, "", ""},
    {"AeroAttachmentFixed", "AeroAttachment", fields_61, "kind", "fixed_reference"},
    {"AeroAttachmentPrescribed", "AeroAttachment", fields_62, "kind", "prescribed_motion"},
    {"AeroAttachmentStructural", "AeroAttachment", fields_63, "kind", "structural"},
    {"AeroSharing", "AeroSharing", fields_64, "", ""},
    {"FlowNone", "Flow", fields_65, "kind", "none"},
    {"FlowKnown", "Flow", fields_66, "kind", "known"},
    {"GustNone", "Gust", fields_67, "kind", "none"},
    {"GustSine", "Gust", fields_68, "kind", "sine"},
    {"FlowReference", "FlowReference", fields_69, "", ""},
    {"WakeNone", "Wake", fields_70, "kind", "none"},
    {"WakeActive", "Wake", fields_71, "kind", "active"},
    {"ConvectionFree", "Convection", fields_72, "kind", "free"},
    {"ConvectionPrescribed", "Convection", fields_73, "kind", "prescribed"},
    {"RetentionAll", "Retention", fields_74, "kind", "all"},
    {"RetentionWindow", "Retention", fields_75, "kind", "closed_window"},
    {"AnalysisStatic", "Analysis", fields_76, "kind", "static"},
    {"AnalysisTransient", "Analysis", fields_77, "kind", "transient"},
    {"AnalysisFlutter", "Analysis", fields_78, "kind", "flutter"},
    {"AnalysisSpeedNone", "AnalysisSpeed", fields_79, "kind", "none"},
    {"AnalysisSpeedValue", "AnalysisSpeed", fields_80, "kind", "value"},
    {"Branch", "Branch", fields_81, "", ""},
    {"BranchSeedReference", "BranchSeed", fields_82, "kind", "reference"},
    {"BranchSeedPrepared", "BranchSeed", fields_83, "kind", "prepared"},
    {"ContinuationFixed", "Continuation", fields_84, "kind", "fixed"},
    {"ContinuationSweep", "Continuation", fields_85, "kind", "sweep"},
    {"ArcLength", "ArcLength", fields_86, "", ""},
    {"TimeInterval", "TimeInterval", fields_87, "", ""},
    {"SpeedInterval", "SpeedInterval", fields_88, "", ""},
    {"Search", "Search", fields_89, "", ""},
    {"Association", "Association", fields_90, "", ""},
    {"Coverage", "Coverage", fields_91, "", ""},
    {"InitialEquilibrium", "Initial", fields_92, "kind", "equilibrium"},
    {"InitialCold", "Initial", fields_93, "kind", "cold_start"},
    {"InitialPrepared", "Initial", fields_94, "kind", "prepared_history"},
    {"PerturbationNone", "Perturbation", fields_95, "kind", "none"},
    {"InitialState", "InitialState", fields_96, "", ""},
    {"CoefficientState", "CoefficientState", fields_97, "", ""},
    {"PlanesBinary64", "Planes", fields_98, "kind", "binary64"},
    {"PlanesTwoWord", "Planes", fields_99, "kind", "two_word"},
    {"BeamState", "BeamState", fields_100, "", ""},
    {"RigidState", "RigidState", fields_101, "", ""},
    {"FrameStatePayload", "FrameStatePayload", fields_102, "", ""},
    {"InitialWake", "InitialWake", fields_103, "", ""},
    {"SurfaceWake", "SurfaceWake", fields_104, "", ""},
    {"Restart", "Restart", fields_105, "", ""},
    {"ObservationPoint", "Observation", fields_106, "kind", "point"},
    {"ObservationStress", "Observation", fields_107, "kind", "stress"},
    {"ObservationSection", "Observation", fields_108, "kind", "section"},
    {"ObservationReaction", "Observation", fields_109, "kind", "reaction"},
    {"ObservationAero", "Observation", fields_110, "kind", "aero"},
    {"ObservationMode", "Observation", fields_111, "kind", "mode"},
    {"BasisWorldMaterial", "Basis", fields_112, "kind",
     "world,reference_material,current_material"},
    {"BasisExplicit", "Basis", fields_113, "kind", "explicit"},
    {"Comparison", "Comparison", fields_114, "", ""},
    {"ComparisonReferenceNone", "ComparisonReference", fields_115, "kind", "none"},
    {"ComparisonReferenceFile", "ComparisonReference", fields_116, "kind", "file"},
    {"ObservationBudget", "ObservationBudget", fields_117, "", ""},
    {"ErrorGroupNotApplicable", "ErrorGroup", fields_118, "status", "not_applicable"},
    {"ErrorGroupCoveredBy", "ErrorGroup", fields_119, "status", "covered_by"},
    {"ErrorGroupRequired", "ErrorGroup", fields_120, "status", "required"},
    {"Numerics", "Numerics", fields_121, "", ""},
    {"WakeCapacityNone", "WakeCapacity", fields_122, "kind", "none"},
    {"WakeCapacityLimited", "WakeCapacity", fields_123, "kind", "limited"},
    {"PatchPolicyNone", "PatchPolicy", fields_124, "kind", "none"},
    {"PatchPolicyActive", "PatchPolicy", fields_125, "kind", "active"},
    {"Quadrature", "Quadrature", fields_126, "", ""},
    {"StructuralQuadratureNone", "StructuralQuadrature", fields_127, "kind", "none"},
    {"StructuralQuadratureActive", "StructuralQuadrature", fields_128, "kind", "active"},
    {"AeroQuadratureNone", "AeroQuadrature", fields_129, "kind", "none"},
    {"AeroQuadratureActive", "AeroQuadrature", fields_130, "kind", "wp1_induction_v1"},
    {"ScaleUniformComponents", "Scale", fields_131, "kind", "uniform_components"},
    {"ScalePhysicalLayoutPolicy", "Scale", fields_132, "kind", "physical_layout_policy"},
    {"ConstraintScalePolicyNone", "ConstraintScalePolicy", fields_133, "kind", "none"},
    {"ConstraintScalePolicyPhysical", "ConstraintScalePolicy", fields_134, "kind", "physical"},
    {"ComponentScale3None", "ComponentScale3", fields_135, "kind", "none"},
    {"ComponentScale3Values", "ComponentScale3", fields_136, "kind", "values"},
    {"ScaleExplicitLayout", "Scale", fields_137, "kind", "explicit_layout"},
    {"LinearPolicyDirect", "LinearPolicy", fields_138, "backend", "dense_lu,sparse_lu,sparse_llt"},
    {"LinearPolicyGmres", "LinearPolicy", fields_139, "backend", "gmres"},
    {"LinearAccuracy", "LinearAccuracy", fields_140, "", ""},
    {"NewtonPolicy", "NewtonPolicy", fields_141, "", ""},
    {"PrecisionPolicy", "PrecisionPolicy", fields_142, "", ""},
    {"BackendMemory", "BackendMemory", fields_143, "", ""},
    {"EigenPolicyNone", "EigenPolicy", fields_144, "kind", "none"},
    {"EigenPolicyFull", "EigenPolicy", fields_145, "kind", "full"},
    {"EigenPolicyGeneralized", "EigenPolicy", fields_146, "kind", "generalized"},
    {"EigenPolicyArnoldi", "EigenPolicy", fields_147, "kind", "arnoldi"},
    {"Complex", "Complex", fields_148, "", ""},
    {"ArnoldiSeedUniform", "ArnoldiSeed", fields_149, "kind", "uniform_nonzero"},
    {"ArnoldiSeedExplicit", "ArnoldiSeed", fields_150, "kind", "explicit_layout"},
    {"EigenSelectorLargest", "EigenSelector", fields_151, "kind", "largest_magnitude"},
    {"EigenSelectorShift", "EigenSelector", fields_152, "kind", "shift_invert"},
    {"EigenRequirements", "EigenRequirements", fields_153, "", ""},
    {"BackwardNone", "Backward", fields_154, "kind", "not_requested"},
    {"BackwardRequired", "Backward", fields_155, "kind", "required"},
    {"EigenCluster", "EigenCluster", fields_156, "", ""},
    {"EigenLimits", "EigenLimits", fields_157, "", ""},
    {"Resources", "Resources", fields_158, "", ""},
    {"InputQuotas", "InputQuotas", fields_159, "", ""},
    {"Output", "Output", fields_160, "", ""},
    {"OutputCadence", "OutputCadence", fields_161, "", ""},
    {"Visualization", "Visualization", fields_162, "", ""},
    {"CheckpointOutputOff", "CheckpointOutput", fields_163, "kind", "off"},
    {"CheckpointOutputOn", "CheckpointOutput", fields_164, "kind", "on"},
    {"ReferenceFieldExpression", "ReferenceField", fields_165, "kind", "expression"},
    {"TimeSeriesExpression", "TimeSeries", fields_166, "kind", "expression"},
    {"LoadCurrentExpression", "Load", fields_167, "kind", "current_pressure"},
    {"Expression", "Expression", fields_168, "", ""},
    {"ExpressionDomainReference", "ExpressionDomain", fields_169, "kind", "reference"},
    {"ExpressionDomainTime", "ExpressionDomain", fields_170, "kind", "time"},
    {"ExpressionDomainCurrentPositionTime", "ExpressionDomain", fields_171, "kind",
     "current_position_time"},
    {"ExprNodeConstant", "ExprNode", fields_172, "op", "constant"},
    {"ExprNodeVariable", "ExprNode", fields_173, "op", "variable"},
    {"ExprNodeBinary", "ExprNode", fields_174, "op", "add,subtract,multiply,divide"},
    {"ExprNodeNegate", "ExprNode", fields_175, "op", "negate"},
    {"ExprNodePower", "ExprNode", fields_176, "op", "integer_power"},
    {"ExprNodeAnalytic", "ExprNode", fields_177, "op", "sin,cos,exp,log,sqrt"},
    {"PerturbationPhysical", "Perturbation", fields_178, "kind", "physical"},
    {"ShellPhysicalState", "ShellPhysicalState", fields_179, "", ""},
    {"FullShellFieldGlobalPolynomial", "FullShellField", fields_180, "kind", "global_polynomial"},
    {"FullShellFieldHalfCells", "FullShellField", fields_181, "kind", "half_cells"},
    {"CellField", "CellField", fields_182, "", ""},
    {"VectorPolynomial", "VectorPolynomial", fields_183, "", ""},
    {"VectorTerm", "VectorTerm", fields_184, "", ""},
    {"PointCheck", "PointCheck", fields_185, "", ""},
    {"InitialPointShell", "InitialPoint", fields_186, "target.kind", "shell"},
    {"InitialPointBeam", "InitialPoint", fields_187, "target.kind", "beam"},
    {"InitialPointRigid", "InitialPoint", fields_188, "target.kind", "rigid"},
    {"Rational", "Rational", fields_189, "", ""},
    {"FileReference", "FileReference", fields_190, "", ""},
}};
const std::array<TypeDefinition, 447> type_table{{
    {"Tag{uvlas.case/1}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Units", FieldKind::Object, std::nullopt, 0, "Units", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Model", FieldKind::Object, std::nullopt, 0, "Model", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Flow", FieldKind::Object, std::nullopt, 0, "Flow", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Analysis", FieldKind::Object, std::nullopt, 0, "Analysis", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Observation>", FieldKind::Array, FieldType::T415, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Numerics", FieldKind::Object, std::nullopt, 0, "Numerics", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Resources", FieldKind::Object, std::nullopt, 0, "Resources", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Output", FieldKind::Object, std::nullopt, 0, "Output", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Initial", FieldKind::Object, std::nullopt, 0, "Initial", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Restart", FieldKind::Object, std::nullopt, 0, "Restart", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Provenance", FieldKind::Object, std::nullopt, 0, "Provenance", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{m,mm,cm,km,in,ft}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{kg,g,t,lbm}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{s,ms,min}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{rad,deg}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Geometry>", FieldKind::Array, FieldType::T393, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Material>", FieldKind::Array, FieldType::T397, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Shell>", FieldKind::Array, FieldType::T399, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Beam>", FieldKind::Array, FieldType::T400, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Joint>", FieldKind::Array, FieldType::T401, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Support>", FieldKind::Array, FieldType::T404, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Mass>", FieldKind::Array, FieldType::T405, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<AeroSurface>", FieldKind::Array, FieldType::T407, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Text", FieldKind::Text, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<Text>", FieldKind::Array, FieldType::T24, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<Q<L>,3>", FieldKind::Array, FieldType::T61, 3, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Mat<Number,3,3>", FieldKind::Array, FieldType::T31, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{shell,beam,rigid}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"対応Ref<Shell/Beam/Mass>", FieldKind::Reference, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Target", FieldKind::Object, std::nullopt, 0, "Target", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<Number,3>", FieldKind::Array, FieldType::T328, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{interior,lower,upper}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{volume,face,edge,point}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Ref<Geometry>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Index", FieldKind::Integer, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Vec<Count,k>", FieldKind::Array, FieldType::T67, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Q<d>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, true, 0},
    {"Tag{constant}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Q<d>,n>", FieldKind::Array, FieldType::T37, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{polynomial}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{chart,material,world_reference}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Vec<List<PolyTerm<d,3>>,n>", FieldKind::Array, FieldType::T426, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Vec<Number,n>", FieldKind::Array, FieldType::T328, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<T>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, 1, 0}}, false, 0},
    {"List<Vec<Number,n>>", FieldKind::Array, FieldType::T43, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Q<T>,2>", FieldKind::Array, FieldType::T44, 2, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{pchip}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Q<T>>", FieldKind::Array, FieldType::T44, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Outside", FieldKind::Object, std::nullopt, 0, "Outside", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{reject,hold}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{periodic}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Frame", FieldKind::Object, std::nullopt, 0, "Frame", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"TimeSeries<3>", FieldKind::Object, std::nullopt, 0, "TimeSeries",
     SchemaDimension{{0, 0, 0, 0}}, true, 3},
    {"Tag{projected}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ReferenceField<1,3>", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{0, 0, 0, 0}}, false, 3},
    {"ReferenceField<angle,1>", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{0, 0, 0, 1}}, false, 1},
    {"Tag{explicit}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ReferenceField<1,9>", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{0, 0, 0, 0}}, false, 9},
    {"Tag{pg}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Id<Geometry>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<L>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{1, 0, 0, 0}}, false, 0},
    {"Rational", FieldKind::Object, std::nullopt, 0, "Rational", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Count,2>", FieldKind::Array, FieldType::T67, 2, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<PgTerm>", FieldKind::Array, FieldType::T394, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Vec<Rational,2>>", FieldKind::Array, FieldType::T105, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Thickness", FieldKind::Object, std::nullopt, 0, "Thickness", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Count", FieldKind::Integer, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{pg-graph-v1}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Boundary>", FieldKind::Array, FieldType::T395, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Continuity>", FieldKind::Array, FieldType::T396, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Vec<Index,3>>", FieldKind::Array, FieldType::T427, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Q<L>>", FieldKind::Array, FieldType::T61, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{smooth,fold}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{same,reversed}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{surface_graph}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{nastran}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"FileRef<nastran-deck>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"FileRef<nastran-sidecar>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{msc-2022.1}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{structural-extraction}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{gmsh}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"FileRef<gmsh-msh4.1-ascii>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"FileRef<gmsh-sidecar>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"ReferenceGeometry", FieldKind::Object, std::nullopt, 0, "ReferenceGeometry",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{supplied}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"FileRef<reference-geometry>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{fit}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Q<1>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<Ref<Observation>>", FieldKind::Array, FieldType::T344, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{isotropic}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Material>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<pressure>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{-1, 1, -2, 0}}, false,
     0},
    {"Q<M/L³>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{-3, 1, 0, 0}}, false, 0},
    {"Applicability", FieldKind::Object, std::nullopt, 0, "Applicability",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{orthotropic3d}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Vec<Q<pressure>,3>", FieldKind::Array, FieldType::T92, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Q<1>,3>", FieldKind::Array, FieldType::T88, 3, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{effective3d}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Mat<Q<pressure>,6,6>", FieldKind::Array, FieldType::T428, 6, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{mandel,engineering}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{unevaluated}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{declared}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Number,6>", FieldKind::Array, FieldType::T328, 6, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Layer>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<Rational,2>", FieldKind::Array, FieldType::T62, 2, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Ref<Material>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Axis", FieldKind::Object, std::nullopt, 0, "Axis", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"ReferenceField<1,6>", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{0, 0, 0, 0}}, false, 6},
    {"Mat<Q<pressure*T>,6,6>", FieldKind::Array, FieldType::T430, 6, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Id<Shell>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Layer>", FieldKind::Array, FieldType::T398, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Load>", FieldKind::Array, FieldType::T406, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Id<Beam>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<Frame,2>", FieldKind::Array, FieldType::T52, 2, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"BeamCoordinate", FieldKind::Object, std::nullopt, 0, "BeamCoordinate",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Section", FieldKind::Object, std::nullopt, 0, "Section", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Q<energy/L>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{1, 1, -2, 0}}, false,
     0},
    {"Mat<Number,6,6>", FieldKind::Array, FieldType::T103, 6, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{arclength}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{scaled}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{matrix}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<M/L>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{-1, 1, 0, 0}}, false, 0},
    {"Mat<Q<M*L>,3,3>", FieldKind::Array, FieldType::T432, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Q<L>,2>", FieldKind::Array, FieldType::T61, 2, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{free}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{unavailable}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{circle}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Q<L²>,2>", FieldKind::Array, FieldType::T188, 2, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{rectangle}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<L⁴>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{4, 0, 0, 0}}, false, 0},
    {"Tag{adhesive}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Joint>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Region", FieldKind::Object, std::nullopt, 0, "Region", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<TraceMap>", FieldKind::Array, FieldType::T402, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Bool", FieldKind::Boolean, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Vec<Index,2>", FieldKind::Array, FieldType::T35, 2, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<List<MapTerm>,3>", FieldKind::Array, FieldType::T433, 3, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{rigid_attachment}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{weld}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{constraint}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Support>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"MaterialPoint", FieldKind::Object, std::nullopt, 0, "MaterialPoint",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<Index>", FieldKind::Array, FieldType::T35, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Motion", FieldKind::Object, std::nullopt, 0, "Motion", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{trace_constraint}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{elastic}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{particle}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Mass>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Q<M>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 1, 0, 0}}, false, 0},
    {"Ledger", FieldKind::Object, std::nullopt, 0, "Ledger", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{rigid}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Mat<Q<M*L²>,3,3>", FieldKind::Array, FieldType::T435, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{add}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{replace}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Ref<Mass>>", FieldKind::Array, FieldType::T247, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{point_force}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Load>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{world,body}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{couple}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{gravity}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{traction}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{reference,current}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"ReferenceField<pressure,3>", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{-1, 1, -2, 0}}, false, 3},
    {"TimeSeries<1>", FieldKind::Object, std::nullopt, 0, "TimeSeries",
     SchemaDimension{{0, 0, 0, 0}}, true, 1},
    {"Tag{pressure}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ReferenceField<pressure,1>", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{-1, 1, -2, 0}}, false, 1},
    {"Tag{body_force}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ReferenceField<force/L³,3>", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{-2, 1, -2, 0}}, false, 3},
    {"Id<AeroSurface>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"AeroAttachment", FieldKind::Object, std::nullopt, 0, "AeroAttachment",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<Number>", FieldKind::Array, FieldType::T328, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{positive,reversed}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<AeroSharing>", FieldKind::Array, FieldType::T408, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{fixed_reference}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{prescribed_motion}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{structural}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Ref<Shell>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Ref<AeroSurface>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{none}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{known}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Gust", FieldKind::Object, std::nullopt, 0, "Gust", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"FlowReference", FieldKind::Object, std::nullopt, 0, "FlowReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{sine}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<Q<angle/L>,3>", FieldKind::Array, FieldType::T436, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<angle/T>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, -1, 1}}, false,
     0},
    {"Q<angle>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 1}}, false, 0},
    {"Q<L/T>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{1, 0, -1, 0}}, false, 0},
    {"Q<L²>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{2, 0, 0, 0}}, false, 0},
    {"Tag{active}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Convection", FieldKind::Object, std::nullopt, 0, "Convection", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Retention", FieldKind::Object, std::nullopt, 0, "Retention", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{wp1-stream-p0-span}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{dual-quarter-v1}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{prescribed}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{all}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{closed_window}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{static}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{aerodynamic,structural,coupled}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Branch", FieldKind::Object, std::nullopt, 0, "Branch", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"AnalysisSpeed", FieldKind::Object, std::nullopt, 0, "AnalysisSpeed",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Wake", FieldKind::Object, std::nullopt, 0, "Wake", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{transient}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<TimeInterval>", FieldKind::Array, FieldType::T409, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{flutter}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{coupled}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Branch>", FieldKind::Array, FieldType::T199, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<SpeedInterval>", FieldKind::Array, FieldType::T410, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Search", FieldKind::Object, std::nullopt, 0, "Search", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{value}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Branch>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"BranchSeed", FieldKind::Object, std::nullopt, 0, "BranchSeed", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Continuation", FieldKind::Object, std::nullopt, 0, "Continuation",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{reference}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{prepared}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"FileRef<preparation-history>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{fixed}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{sweep}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ArcLength", FieldKind::Object, std::nullopt, 0, "ArcLength", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"PositiveCount", FieldKind::Integer, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Scale>", FieldKind::Array, FieldType::T417, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{previous_tangent}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{initial,event,step_change,chart_change}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{lowest_interval_first}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Q<1/T>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, -1, 0}}, false, 0},
    {"Association", FieldKind::Object, std::nullopt, 0, "Association",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Coverage", FieldKind::Object, std::nullopt, 0, "Coverage", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Ref<Branch>>", FieldKind::Array, FieldType::T229, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{equilibrium}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Ref<Branch>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Perturbation", FieldKind::Object, std::nullopt, 0, "Perturbation",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{cold_start}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"InitialState", FieldKind::Object, std::nullopt, 0, "InitialState",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<Q<L²/T>>", FieldKind::Array, FieldType::T437, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"InitialWake", FieldKind::Object, std::nullopt, 0, "InitialWake",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{prepared_history}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"FileRef<initial-history>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<CoefficientState>", FieldKind::Array, FieldType::T411, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<BeamState>", FieldKind::Array, FieldType::T412, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<RigidState>", FieldKind::Array, FieldType::T413, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Digest", FieldKind::Text, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{absolute,reference_increment}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Planes", FieldKind::Object, std::nullopt, 0, "Planes", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{binary64}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{two_word}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Ref<Beam>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<FrameStatePayload,2>", FieldKind::Array, FieldType::T248, 2, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Ref<Mass>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"FrameStatePayload", FieldKind::Object, std::nullopt, 0, "FrameStatePayload",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Vec<Q<L/T>,3>", FieldKind::Array, FieldType::T187, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Q<angle/T>,3>", FieldKind::Array, FieldType::T185, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<SurfaceWake>", FieldKind::Array, FieldType::T414, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Vec<Q<L>,3>>", FieldKind::Array, FieldType::T26, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Id<WakeRow>>", FieldKind::Array, FieldType::T438, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Id<WakePoint>>", FieldKind::Array, FieldType::T439, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Id<WakeBand>>", FieldKind::Array, FieldType::T440, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{same_environment_bitwise,qualified_cross_environment}", FieldKind::Enumeration,
     std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<ObservationBudget>", FieldKind::Array, FieldType::T263, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{point}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<Observation>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{position,displacement,velocity,acceleration,force}", FieldKind::Enumeration, std::nullopt,
     0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Basis", FieldKind::Object, std::nullopt, 0, "Basis", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Comparison", FieldKind::Object, std::nullopt, 0, "Comparison", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ObservationBudget", FieldKind::Object, std::nullopt, 0, "ObservationBudget",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{stress}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Ref<Layer>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{interior,lower,upper,both}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{green,second_pk,first_pk,cauchy}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{elastic,viscous,total}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{section}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{force,moment}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{reaction}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Ref<Support>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{aero}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Ref<AeroSurface>>", FieldKind::Array, FieldType::T178, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{force,moment,coefficient}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{mode}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Target>", FieldKind::Array, FieldType::T30, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<MaterialPoint>", FieldKind::Array, FieldType::T142, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{frequency,growth,subspace}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{world,reference_material,current_material}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{instantaneous,history,peak,rms}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{left,right,both}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"ComparisonReference", FieldKind::Object, std::nullopt, 0, "ComparisonReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{file}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"FileRef<observation-reference>", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<ErrorGroup>", FieldKind::Array, FieldType::T416, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{not_applicable}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,"
     "newton,eigen_range,search_width,load_mapping}",
     FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{covered_by}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{required}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<Tag{同じ集合}>", FieldKind::Array, FieldType::T441, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"PatchPolicy", FieldKind::Object, std::nullopt, 0, "PatchPolicy",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Quadrature", FieldKind::Object, std::nullopt, 0, "Quadrature", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"LinearPolicy", FieldKind::Object, std::nullopt, 0, "LinearPolicy",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"NewtonPolicy", FieldKind::Object, std::nullopt, 0, "NewtonPolicy",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"PrecisionPolicy", FieldKind::Object, std::nullopt, 0, "PrecisionPolicy",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"EigenPolicy", FieldKind::Object, std::nullopt, 0, "EigenPolicy",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"WakeCapacity", FieldKind::Object, std::nullopt, 0, "WakeCapacity",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{limited}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{exact_no_drop}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"StructuralQuadrature", FieldKind::Object, std::nullopt, 0, "StructuralQuadrature",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"AeroQuadrature", FieldKind::Object, std::nullopt, 0, "AeroQuadrature",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{wp1_induction_v1}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{uniform_components}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{shell_position,beam_translation,beam_rotation,rigid_translation,rigid_rotation,"
     "structural_velocity,structural_angular_velocity,circulation,wake_position,wake_strength,"
     "constraint,force_balance,moment_balance,aero_residual,transport_residual,newton_correction,"
     "arclength_state}",
     FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{physical_layout_policy}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{mixed_blocks}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"ConstraintScalePolicy", FieldKind::Object, std::nullopt, 0, "ConstraintScalePolicy",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{physical}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ComponentScale3", FieldKind::Object, std::nullopt, 0, "ComponentScale3",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Q<energy>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{2, 1, -2, 0}}, false, 0},
    {"Vec<Q<force>,3>", FieldKind::Array, FieldType::T442, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Q<force*L>,3>", FieldKind::Array, FieldType::T443, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{values}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{explicit_layout}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<Vec<signed-native-integer,4>>", FieldKind::Array, FieldType::T368, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{dense_lu,sparse_lu,sparse_llt}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"LinearAccuracy", FieldKind::Object, std::nullopt, 0, "LinearAccuracy",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"BackendMemory", FieldKind::Object, std::nullopt, 0, "BackendMemory",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{gmres}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{infinity,l2}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{preserve}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{adaptive_binary64_two_word}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{retained,binary64_projection}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{estimate,require_proven_bound}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{full}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{step_map}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Number", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{spectrum,pairs,subspace}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"EigenRequirements", FieldKind::Object, std::nullopt, 0, "EigenRequirements",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"EigenLimits", FieldKind::Object, std::nullopt, 0, "EigenLimits",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{generalized}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{structural_modes}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Q<1/T²>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{0, 0, -2, 0}}, false, 0},
    {"Tag{arnoldi}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{pairs,subspace}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"EigenSelector", FieldKind::Object, std::nullopt, 0, "EigenSelector",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"ArnoldiSeed", FieldKind::Object, std::nullopt, 0, "ArnoldiSeed",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{uniform_nonzero}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Complex", FieldKind::Object, std::nullopt, 0, "Complex", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"List<Complex>", FieldKind::Array, FieldType::T340, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{largest_magnitude}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{shift_invert}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Ref<Observation>", FieldKind::Reference, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Backward", FieldKind::Object, std::nullopt, 0, "Backward", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<EigenCluster>", FieldKind::Array, FieldType::T418, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{not_requested}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{24gib,48gib,custom}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"InputQuotas", FieldKind::Object, std::nullopt, 0, "InputQuotas",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"OutputCadence", FieldKind::Object, std::nullopt, 0, "OutputCadence",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Visualization", FieldKind::Object, std::nullopt, 0, "Visualization",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"CheckpointOutput", FieldKind::Object, std::nullopt, 0, "CheckpointOutput",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{reject_existing}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{explicit_only}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<Ref<Layer>>", FieldKind::Array, FieldType::T265, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{both}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{off}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{on}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{expression}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Expression<n>", FieldKind::Object, std::nullopt, 0, "Expression",
     SchemaDimension{{0, 0, 0, 0}}, true, 0},
    {"Tag{current_pressure}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Expression<1>", FieldKind::Object, std::nullopt, 0, "Expression",
     SchemaDimension{{0, 0, 0, 0}}, true, 1},
    {"ExpressionDomain", FieldKind::Object, std::nullopt, 0, "ExpressionDomain",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<ExprNode>", FieldKind::Array, FieldType::T420, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Index,n>", FieldKind::Array, FieldType::T35, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{time}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Tag{current_position_time}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Vec<signed-native-integer,4>", FieldKind::Array, FieldType::T373, 4, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{variable}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{add,subtract,multiply,divide}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{negate}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{integer_power}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"signed-native-integer", FieldKind::Integer, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{sin,cos,exp,log,sqrt}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<ShellPhysicalState>", FieldKind::Array, FieldType::T421, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<PointCheck>", FieldKind::Array, FieldType::T424, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"FullShellField", FieldKind::Object, std::nullopt, 0, "FullShellField",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{global_polynomial}", FieldKind::Enumeration, std::nullopt, 0, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"VectorPolynomial", FieldKind::Object, std::nullopt, 0, "VectorPolynomial",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Tag{half_cells}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<CellField>", FieldKind::Array, FieldType::T422, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<VectorTerm>", FieldKind::Array, FieldType::T423, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<ExactQuantity,3>", FieldKind::Array, FieldType::T444, 3, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"InitialPoint", FieldKind::Object, std::nullopt, 0, "InitialPoint",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Target{kind=shell}", FieldKind::Object, std::nullopt, 0, "Target",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Target{kind=beam}", FieldKind::Object, std::nullopt, 0, "Target",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Target{kind=rigid}", FieldKind::Object, std::nullopt, 0, "Target",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Id<Source>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Case", FieldKind::Object, std::nullopt, 0, "Case", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"PolyTerm", FieldKind::Object, std::nullopt, 0, "PolyTerm", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ReferenceField", FieldKind::Object, std::nullopt, 0, "ReferenceField",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"TimeSeries", FieldKind::Object, std::nullopt, 0, "TimeSeries", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Geometry", FieldKind::Object, std::nullopt, 0, "Geometry", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"PgTerm", FieldKind::Object, std::nullopt, 0, "PgTerm", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Boundary", FieldKind::Object, std::nullopt, 0, "Boundary", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Continuity", FieldKind::Object, std::nullopt, 0, "Continuity", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Material", FieldKind::Object, std::nullopt, 0, "Material", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Layer", FieldKind::Object, std::nullopt, 0, "Layer", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Shell", FieldKind::Object, std::nullopt, 0, "Shell", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Beam", FieldKind::Object, std::nullopt, 0, "Beam", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Joint", FieldKind::Object, std::nullopt, 0, "Joint", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"TraceMap", FieldKind::Object, std::nullopt, 0, "TraceMap", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"MapTerm", FieldKind::Object, std::nullopt, 0, "MapTerm", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Support", FieldKind::Object, std::nullopt, 0, "Support", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Mass", FieldKind::Object, std::nullopt, 0, "Mass", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Load", FieldKind::Object, std::nullopt, 0, "Load", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"AeroSurface", FieldKind::Object, std::nullopt, 0, "AeroSurface",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"AeroSharing", FieldKind::Object, std::nullopt, 0, "AeroSharing",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"TimeInterval", FieldKind::Object, std::nullopt, 0, "TimeInterval",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"SpeedInterval", FieldKind::Object, std::nullopt, 0, "SpeedInterval",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"CoefficientState", FieldKind::Object, std::nullopt, 0, "CoefficientState",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"BeamState", FieldKind::Object, std::nullopt, 0, "BeamState", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"RigidState", FieldKind::Object, std::nullopt, 0, "RigidState", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"SurfaceWake", FieldKind::Object, std::nullopt, 0, "SurfaceWake",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Observation", FieldKind::Object, std::nullopt, 0, "Observation",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"ErrorGroup", FieldKind::Object, std::nullopt, 0, "ErrorGroup", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Scale", FieldKind::Object, std::nullopt, 0, "Scale", SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"EigenCluster", FieldKind::Object, std::nullopt, 0, "EigenCluster",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Expression", FieldKind::Object, std::nullopt, 0, "Expression", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ExprNode", FieldKind::Object, std::nullopt, 0, "ExprNode", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"ShellPhysicalState", FieldKind::Object, std::nullopt, 0, "ShellPhysicalState",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"CellField", FieldKind::Object, std::nullopt, 0, "CellField", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"VectorTerm", FieldKind::Object, std::nullopt, 0, "VectorTerm", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"PointCheck", FieldKind::Object, std::nullopt, 0, "PointCheck", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"FileReference", FieldKind::Object, std::nullopt, 0, "FileReference",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"List<PolyTerm<d,3>>", FieldKind::Array, FieldType::T445, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Vec<Index,3>", FieldKind::Array, FieldType::T35, 3, "", SchemaDimension{{0, 0, 0, 0}}, false,
     0},
    {"Vec<Q<pressure>,6>", FieldKind::Array, FieldType::T92, 6, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<pressure*T>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{-1, 1, -1, 0}},
     false, 0},
    {"Vec<Q<pressure*T>,6>", FieldKind::Array, FieldType::T429, 6, "",
     SchemaDimension{{0, 0, 0, 0}}, false, 0},
    {"Q<M*L>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{1, 1, 0, 0}}, false, 0},
    {"Vec<Q<M*L>,3>", FieldKind::Array, FieldType::T431, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"List<MapTerm>", FieldKind::Array, FieldType::T403, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<M*L²>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{2, 1, 0, 0}}, false, 0},
    {"Vec<Q<M*L²>,3>", FieldKind::Array, FieldType::T434, 3, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<angle/L>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{-1, 0, 0, 1}}, false,
     0},
    {"Q<L²/T>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{2, 0, -1, 0}}, false, 0},
    {"Id<WakeRow>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<WakePoint>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Id<WakeBand>", FieldKind::Identifier, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{同じ集合}", FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Q<force>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{1, 1, -2, 0}}, false, 0},
    {"Q<force*L>", FieldKind::Number, std::nullopt, 0, "", SchemaDimension{{2, 1, -2, 0}}, false,
     0},
    {"ExactQuantity", FieldKind::ExactQuantity, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"PolyTerm<d,3>", FieldKind::Object, std::nullopt, 0, "PolyTerm", SchemaDimension{{0, 0, 0, 0}},
     false, 0},
    {"Tag{nastran-deck,nastran-sidecar,gmsh-msh4.1-ascii,gmsh-sidecar,reference-geometry,"
     "preparation-history,initial-history,observation-reference}",
     FieldKind::Enumeration, std::nullopt, 0, "", SchemaDimension{{0, 0, 0, 0}}, false, 0},
}};
} // namespace
std::span<const TypeDefinition> schema_types() noexcept {
    return type_table;
}
std::span<const RecordDefinition> schema_records() noexcept {
    return records;
}
FieldType schema_root_type() noexcept {
    return FieldType::T389;
}
std::uint64_t schema_record_bytes(SchemaContext& context) {
    std::uint64_t total = 0;
    total = context.add(total, context.multiply(context.output->record_counts[0],
                                                sizeof(std::optional<schema::Case>)));
    total = context.add(total, context.multiply(context.output->record_counts[1],
                                                sizeof(std::optional<schema::Units>)));
    total = context.add(total, context.multiply(context.output->record_counts[2],
                                                sizeof(std::optional<schema::Model>)));
    total = context.add(total, context.multiply(context.output->record_counts[3],
                                                sizeof(std::optional<schema::Provenance>)));
    total = context.add(total, context.multiply(context.output->record_counts[4],
                                                sizeof(std::optional<schema::Frame>)));
    total = context.add(total, context.multiply(context.output->record_counts[5],
                                                sizeof(std::optional<schema::Target>)));
    total = context.add(total, context.multiply(context.output->record_counts[6],
                                                sizeof(std::optional<schema::MaterialPoint>)));
    total = context.add(total, context.multiply(context.output->record_counts[7],
                                                sizeof(std::optional<schema::Region>)));
    total = context.add(total, context.multiply(context.output->record_counts[8],
                                                sizeof(std::optional<schema::PolyTerm>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[9],
                                            sizeof(std::optional<schema::ReferenceFieldConstant>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[10],
                                         sizeof(std::optional<schema::ReferenceFieldPolynomial>)));
    total = context.add(total, context.multiply(context.output->record_counts[11],
                                                sizeof(std::optional<schema::TimeSeriesConstant>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[12],
                                            sizeof(std::optional<schema::TimeSeriesPolynomial>)));
    total = context.add(total, context.multiply(context.output->record_counts[13],
                                                sizeof(std::optional<schema::TimeSeriesPchip>)));
    total = context.add(total, context.multiply(context.output->record_counts[14],
                                                sizeof(std::optional<schema::OutsideRejectHold>)));
    total = context.add(total, context.multiply(context.output->record_counts[15],
                                                sizeof(std::optional<schema::OutsidePeriodic>)));
    total = context.add(total, context.multiply(context.output->record_counts[16],
                                                sizeof(std::optional<schema::Motion>)));
    total = context.add(total, context.multiply(context.output->record_counts[17],
                                                sizeof(std::optional<schema::AxisProjected>)));
    total = context.add(total, context.multiply(context.output->record_counts[18],
                                                sizeof(std::optional<schema::AxisExplicit>)));
    total = context.add(total, context.multiply(context.output->record_counts[19],
                                                sizeof(std::optional<schema::GeometryPg>)));
    total = context.add(total, context.multiply(context.output->record_counts[20],
                                                sizeof(std::optional<schema::PgTerm>)));
    total = context.add(total, context.multiply(context.output->record_counts[21],
                                                sizeof(std::optional<schema::Thickness>)));
    total = context.add(total, context.multiply(context.output->record_counts[22],
                                                sizeof(std::optional<schema::Boundary>)));
    total = context.add(total, context.multiply(context.output->record_counts[23],
                                                sizeof(std::optional<schema::Continuity>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[24],
                                            sizeof(std::optional<schema::GeometrySurfaceGraph>)));
    total = context.add(total, context.multiply(context.output->record_counts[25],
                                                sizeof(std::optional<schema::GeometryNastran>)));
    total = context.add(total, context.multiply(context.output->record_counts[26],
                                                sizeof(std::optional<schema::GeometryGmsh>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[27],
                                         sizeof(std::optional<schema::ReferenceGeometrySupplied>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[28],
                                            sizeof(std::optional<schema::ReferenceGeometryFit>)));
    total = context.add(total, context.multiply(context.output->record_counts[29],
                                                sizeof(std::optional<schema::MaterialIsotropic>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[30],
                                            sizeof(std::optional<schema::MaterialOrthotropic>)));
    total = context.add(total, context.multiply(context.output->record_counts[31],
                                                sizeof(std::optional<schema::MaterialEffective>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[32],
                                         sizeof(std::optional<schema::ApplicabilityUnevaluated>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[33],
                                            sizeof(std::optional<schema::ApplicabilityDeclared>)));
    total = context.add(total, context.multiply(context.output->record_counts[34],
                                                sizeof(std::optional<schema::Layer>)));
    total = context.add(total, context.multiply(context.output->record_counts[35],
                                                sizeof(std::optional<schema::Shell>)));
    total = context.add(total, context.multiply(context.output->record_counts[36],
                                                sizeof(std::optional<schema::Beam>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[37],
                                         sizeof(std::optional<schema::BeamCoordinateArclength>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[38],
                                            sizeof(std::optional<schema::BeamCoordinateScaled>)));
    total = context.add(total, context.multiply(context.output->record_counts[39],
                                                sizeof(std::optional<schema::SectionMatrix>)));
    total = context.add(total, context.multiply(context.output->record_counts[40],
                                                sizeof(std::optional<schema::SectionCircle>)));
    total = context.add(total, context.multiply(context.output->record_counts[41],
                                                sizeof(std::optional<schema::SectionRectangle>)));
    total = context.add(total, context.multiply(context.output->record_counts[42],
                                                sizeof(std::optional<schema::JointAdhesive>)));
    total = context.add(total, context.multiply(context.output->record_counts[43],
                                                sizeof(std::optional<schema::TraceMap>)));
    total = context.add(total, context.multiply(context.output->record_counts[44],
                                                sizeof(std::optional<schema::MapTerm>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[45],
                                            sizeof(std::optional<schema::JointRigidAttachment>)));
    total = context.add(total, context.multiply(context.output->record_counts[46],
                                                sizeof(std::optional<schema::JointWeld>)));
    total = context.add(total, context.multiply(context.output->record_counts[47],
                                                sizeof(std::optional<schema::SupportConstraint>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[48],
                                            sizeof(std::optional<schema::SupportTraceConstraint>)));
    total = context.add(total, context.multiply(context.output->record_counts[49],
                                                sizeof(std::optional<schema::SupportElastic>)));
    total = context.add(total, context.multiply(context.output->record_counts[50],
                                                sizeof(std::optional<schema::MassParticle>)));
    total = context.add(total, context.multiply(context.output->record_counts[51],
                                                sizeof(std::optional<schema::MassRigid>)));
    total = context.add(total, context.multiply(context.output->record_counts[52],
                                                sizeof(std::optional<schema::LedgerAdd>)));
    total = context.add(total, context.multiply(context.output->record_counts[53],
                                                sizeof(std::optional<schema::LedgerReplace>)));
    total = context.add(total, context.multiply(context.output->record_counts[54],
                                                sizeof(std::optional<schema::LoadPointForce>)));
    total = context.add(total, context.multiply(context.output->record_counts[55],
                                                sizeof(std::optional<schema::LoadCouple>)));
    total = context.add(total, context.multiply(context.output->record_counts[56],
                                                sizeof(std::optional<schema::LoadGravity>)));
    total = context.add(total, context.multiply(context.output->record_counts[57],
                                                sizeof(std::optional<schema::LoadTraction>)));
    total = context.add(total, context.multiply(context.output->record_counts[58],
                                                sizeof(std::optional<schema::LoadPressure>)));
    total = context.add(total, context.multiply(context.output->record_counts[59],
                                                sizeof(std::optional<schema::LoadBody>)));
    total = context.add(total, context.multiply(context.output->record_counts[60],
                                                sizeof(std::optional<schema::AeroSurface>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[61],
                                            sizeof(std::optional<schema::AeroAttachmentFixed>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[62],
                                         sizeof(std::optional<schema::AeroAttachmentPrescribed>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[63],
                                         sizeof(std::optional<schema::AeroAttachmentStructural>)));
    total = context.add(total, context.multiply(context.output->record_counts[64],
                                                sizeof(std::optional<schema::AeroSharing>)));
    total = context.add(total, context.multiply(context.output->record_counts[65],
                                                sizeof(std::optional<schema::FlowNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[66],
                                                sizeof(std::optional<schema::FlowKnown>)));
    total = context.add(total, context.multiply(context.output->record_counts[67],
                                                sizeof(std::optional<schema::GustNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[68],
                                                sizeof(std::optional<schema::GustSine>)));
    total = context.add(total, context.multiply(context.output->record_counts[69],
                                                sizeof(std::optional<schema::FlowReference>)));
    total = context.add(total, context.multiply(context.output->record_counts[70],
                                                sizeof(std::optional<schema::WakeNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[71],
                                                sizeof(std::optional<schema::WakeActive>)));
    total = context.add(total, context.multiply(context.output->record_counts[72],
                                                sizeof(std::optional<schema::ConvectionFree>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[73],
                                            sizeof(std::optional<schema::ConvectionPrescribed>)));
    total = context.add(total, context.multiply(context.output->record_counts[74],
                                                sizeof(std::optional<schema::RetentionAll>)));
    total = context.add(total, context.multiply(context.output->record_counts[75],
                                                sizeof(std::optional<schema::RetentionWindow>)));
    total = context.add(total, context.multiply(context.output->record_counts[76],
                                                sizeof(std::optional<schema::AnalysisStatic>)));
    total = context.add(total, context.multiply(context.output->record_counts[77],
                                                sizeof(std::optional<schema::AnalysisTransient>)));
    total = context.add(total, context.multiply(context.output->record_counts[78],
                                                sizeof(std::optional<schema::AnalysisFlutter>)));
    total = context.add(total, context.multiply(context.output->record_counts[79],
                                                sizeof(std::optional<schema::AnalysisSpeedNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[80],
                                                sizeof(std::optional<schema::AnalysisSpeedValue>)));
    total = context.add(total, context.multiply(context.output->record_counts[81],
                                                sizeof(std::optional<schema::Branch>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[82],
                                            sizeof(std::optional<schema::BranchSeedReference>)));
    total = context.add(total, context.multiply(context.output->record_counts[83],
                                                sizeof(std::optional<schema::BranchSeedPrepared>)));
    total = context.add(total, context.multiply(context.output->record_counts[84],
                                                sizeof(std::optional<schema::ContinuationFixed>)));
    total = context.add(total, context.multiply(context.output->record_counts[85],
                                                sizeof(std::optional<schema::ContinuationSweep>)));
    total = context.add(total, context.multiply(context.output->record_counts[86],
                                                sizeof(std::optional<schema::ArcLength>)));
    total = context.add(total, context.multiply(context.output->record_counts[87],
                                                sizeof(std::optional<schema::TimeInterval>)));
    total = context.add(total, context.multiply(context.output->record_counts[88],
                                                sizeof(std::optional<schema::SpeedInterval>)));
    total = context.add(total, context.multiply(context.output->record_counts[89],
                                                sizeof(std::optional<schema::Search>)));
    total = context.add(total, context.multiply(context.output->record_counts[90],
                                                sizeof(std::optional<schema::Association>)));
    total = context.add(total, context.multiply(context.output->record_counts[91],
                                                sizeof(std::optional<schema::Coverage>)));
    total = context.add(total, context.multiply(context.output->record_counts[92],
                                                sizeof(std::optional<schema::InitialEquilibrium>)));
    total = context.add(total, context.multiply(context.output->record_counts[93],
                                                sizeof(std::optional<schema::InitialCold>)));
    total = context.add(total, context.multiply(context.output->record_counts[94],
                                                sizeof(std::optional<schema::InitialPrepared>)));
    total = context.add(total, context.multiply(context.output->record_counts[95],
                                                sizeof(std::optional<schema::PerturbationNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[96],
                                                sizeof(std::optional<schema::InitialState>)));
    total = context.add(total, context.multiply(context.output->record_counts[97],
                                                sizeof(std::optional<schema::CoefficientState>)));
    total = context.add(total, context.multiply(context.output->record_counts[98],
                                                sizeof(std::optional<schema::PlanesBinary64>)));
    total = context.add(total, context.multiply(context.output->record_counts[99],
                                                sizeof(std::optional<schema::PlanesTwoWord>)));
    total = context.add(total, context.multiply(context.output->record_counts[100],
                                                sizeof(std::optional<schema::BeamState>)));
    total = context.add(total, context.multiply(context.output->record_counts[101],
                                                sizeof(std::optional<schema::RigidState>)));
    total = context.add(total, context.multiply(context.output->record_counts[102],
                                                sizeof(std::optional<schema::FrameStatePayload>)));
    total = context.add(total, context.multiply(context.output->record_counts[103],
                                                sizeof(std::optional<schema::InitialWake>)));
    total = context.add(total, context.multiply(context.output->record_counts[104],
                                                sizeof(std::optional<schema::SurfaceWake>)));
    total = context.add(total, context.multiply(context.output->record_counts[105],
                                                sizeof(std::optional<schema::Restart>)));
    total = context.add(total, context.multiply(context.output->record_counts[106],
                                                sizeof(std::optional<schema::ObservationPoint>)));
    total = context.add(total, context.multiply(context.output->record_counts[107],
                                                sizeof(std::optional<schema::ObservationStress>)));
    total = context.add(total, context.multiply(context.output->record_counts[108],
                                                sizeof(std::optional<schema::ObservationSection>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[109],
                                            sizeof(std::optional<schema::ObservationReaction>)));
    total = context.add(total, context.multiply(context.output->record_counts[110],
                                                sizeof(std::optional<schema::ObservationAero>)));
    total = context.add(total, context.multiply(context.output->record_counts[111],
                                                sizeof(std::optional<schema::ObservationMode>)));
    total = context.add(total, context.multiply(context.output->record_counts[112],
                                                sizeof(std::optional<schema::BasisWorldMaterial>)));
    total = context.add(total, context.multiply(context.output->record_counts[113],
                                                sizeof(std::optional<schema::BasisExplicit>)));
    total = context.add(total, context.multiply(context.output->record_counts[114],
                                                sizeof(std::optional<schema::Comparison>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[115],
                                         sizeof(std::optional<schema::ComparisonReferenceNone>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[116],
                                         sizeof(std::optional<schema::ComparisonReferenceFile>)));
    total = context.add(total, context.multiply(context.output->record_counts[117],
                                                sizeof(std::optional<schema::ObservationBudget>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[118],
                                         sizeof(std::optional<schema::ErrorGroupNotApplicable>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[119],
                                            sizeof(std::optional<schema::ErrorGroupCoveredBy>)));
    total = context.add(total, context.multiply(context.output->record_counts[120],
                                                sizeof(std::optional<schema::ErrorGroupRequired>)));
    total = context.add(total, context.multiply(context.output->record_counts[121],
                                                sizeof(std::optional<schema::Numerics>)));
    total = context.add(total, context.multiply(context.output->record_counts[122],
                                                sizeof(std::optional<schema::WakeCapacityNone>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[123],
                                            sizeof(std::optional<schema::WakeCapacityLimited>)));
    total = context.add(total, context.multiply(context.output->record_counts[124],
                                                sizeof(std::optional<schema::PatchPolicyNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[125],
                                                sizeof(std::optional<schema::PatchPolicyActive>)));
    total = context.add(total, context.multiply(context.output->record_counts[126],
                                                sizeof(std::optional<schema::Quadrature>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[127],
                                         sizeof(std::optional<schema::StructuralQuadratureNone>)));
    total = context.add(
        total, context.multiply(context.output->record_counts[128],
                                sizeof(std::optional<schema::StructuralQuadratureActive>)));
    total = context.add(total, context.multiply(context.output->record_counts[129],
                                                sizeof(std::optional<schema::AeroQuadratureNone>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[130],
                                            sizeof(std::optional<schema::AeroQuadratureActive>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[131],
                                            sizeof(std::optional<schema::ScaleUniformComponents>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[132],
                                         sizeof(std::optional<schema::ScalePhysicalLayoutPolicy>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[133],
                                         sizeof(std::optional<schema::ConstraintScalePolicyNone>)));
    total = context.add(
        total, context.multiply(context.output->record_counts[134],
                                sizeof(std::optional<schema::ConstraintScalePolicyPhysical>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[135],
                                            sizeof(std::optional<schema::ComponentScale3None>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[136],
                                            sizeof(std::optional<schema::ComponentScale3Values>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[137],
                                            sizeof(std::optional<schema::ScaleExplicitLayout>)));
    total = context.add(total, context.multiply(context.output->record_counts[138],
                                                sizeof(std::optional<schema::LinearPolicyDirect>)));
    total = context.add(total, context.multiply(context.output->record_counts[139],
                                                sizeof(std::optional<schema::LinearPolicyGmres>)));
    total = context.add(total, context.multiply(context.output->record_counts[140],
                                                sizeof(std::optional<schema::LinearAccuracy>)));
    total = context.add(total, context.multiply(context.output->record_counts[141],
                                                sizeof(std::optional<schema::NewtonPolicy>)));
    total = context.add(total, context.multiply(context.output->record_counts[142],
                                                sizeof(std::optional<schema::PrecisionPolicy>)));
    total = context.add(total, context.multiply(context.output->record_counts[143],
                                                sizeof(std::optional<schema::BackendMemory>)));
    total = context.add(total, context.multiply(context.output->record_counts[144],
                                                sizeof(std::optional<schema::EigenPolicyNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[145],
                                                sizeof(std::optional<schema::EigenPolicyFull>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[146],
                                            sizeof(std::optional<schema::EigenPolicyGeneralized>)));
    total = context.add(total, context.multiply(context.output->record_counts[147],
                                                sizeof(std::optional<schema::EigenPolicyArnoldi>)));
    total = context.add(total, context.multiply(context.output->record_counts[148],
                                                sizeof(std::optional<schema::Complex>)));
    total = context.add(total, context.multiply(context.output->record_counts[149],
                                                sizeof(std::optional<schema::ArnoldiSeedUniform>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[150],
                                            sizeof(std::optional<schema::ArnoldiSeedExplicit>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[151],
                                            sizeof(std::optional<schema::EigenSelectorLargest>)));
    total = context.add(total, context.multiply(context.output->record_counts[152],
                                                sizeof(std::optional<schema::EigenSelectorShift>)));
    total = context.add(total, context.multiply(context.output->record_counts[153],
                                                sizeof(std::optional<schema::EigenRequirements>)));
    total = context.add(total, context.multiply(context.output->record_counts[154],
                                                sizeof(std::optional<schema::BackwardNone>)));
    total = context.add(total, context.multiply(context.output->record_counts[155],
                                                sizeof(std::optional<schema::BackwardRequired>)));
    total = context.add(total, context.multiply(context.output->record_counts[156],
                                                sizeof(std::optional<schema::EigenCluster>)));
    total = context.add(total, context.multiply(context.output->record_counts[157],
                                                sizeof(std::optional<schema::EigenLimits>)));
    total = context.add(total, context.multiply(context.output->record_counts[158],
                                                sizeof(std::optional<schema::Resources>)));
    total = context.add(total, context.multiply(context.output->record_counts[159],
                                                sizeof(std::optional<schema::InputQuotas>)));
    total = context.add(total, context.multiply(context.output->record_counts[160],
                                                sizeof(std::optional<schema::Output>)));
    total = context.add(total, context.multiply(context.output->record_counts[161],
                                                sizeof(std::optional<schema::OutputCadence>)));
    total = context.add(total, context.multiply(context.output->record_counts[162],
                                                sizeof(std::optional<schema::Visualization>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[163],
                                            sizeof(std::optional<schema::CheckpointOutputOff>)));
    total = context.add(total, context.multiply(context.output->record_counts[164],
                                                sizeof(std::optional<schema::CheckpointOutputOn>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[165],
                                         sizeof(std::optional<schema::ReferenceFieldExpression>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[166],
                                            sizeof(std::optional<schema::TimeSeriesExpression>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[167],
                                            sizeof(std::optional<schema::LoadCurrentExpression>)));
    total = context.add(total, context.multiply(context.output->record_counts[168],
                                                sizeof(std::optional<schema::Expression>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[169],
                                         sizeof(std::optional<schema::ExpressionDomainReference>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[170],
                                            sizeof(std::optional<schema::ExpressionDomainTime>)));
    total = context.add(
        total,
        context.multiply(context.output->record_counts[171],
                         sizeof(std::optional<schema::ExpressionDomainCurrentPositionTime>)));
    total = context.add(total, context.multiply(context.output->record_counts[172],
                                                sizeof(std::optional<schema::ExprNodeConstant>)));
    total = context.add(total, context.multiply(context.output->record_counts[173],
                                                sizeof(std::optional<schema::ExprNodeVariable>)));
    total = context.add(total, context.multiply(context.output->record_counts[174],
                                                sizeof(std::optional<schema::ExprNodeBinary>)));
    total = context.add(total, context.multiply(context.output->record_counts[175],
                                                sizeof(std::optional<schema::ExprNodeNegate>)));
    total = context.add(total, context.multiply(context.output->record_counts[176],
                                                sizeof(std::optional<schema::ExprNodePower>)));
    total = context.add(total, context.multiply(context.output->record_counts[177],
                                                sizeof(std::optional<schema::ExprNodeAnalytic>)));
    total =
        context.add(total, context.multiply(context.output->record_counts[178],
                                            sizeof(std::optional<schema::PerturbationPhysical>)));
    total = context.add(total, context.multiply(context.output->record_counts[179],
                                                sizeof(std::optional<schema::ShellPhysicalState>)));
    total = context.add(
        total, context.multiply(context.output->record_counts[180],
                                sizeof(std::optional<schema::FullShellFieldGlobalPolynomial>)));
    total = context.add(total,
                        context.multiply(context.output->record_counts[181],
                                         sizeof(std::optional<schema::FullShellFieldHalfCells>)));
    total = context.add(total, context.multiply(context.output->record_counts[182],
                                                sizeof(std::optional<schema::CellField>)));
    total = context.add(total, context.multiply(context.output->record_counts[183],
                                                sizeof(std::optional<schema::VectorPolynomial>)));
    total = context.add(total, context.multiply(context.output->record_counts[184],
                                                sizeof(std::optional<schema::VectorTerm>)));
    total = context.add(total, context.multiply(context.output->record_counts[185],
                                                sizeof(std::optional<schema::PointCheck>)));
    total = context.add(total, context.multiply(context.output->record_counts[186],
                                                sizeof(std::optional<schema::InitialPointShell>)));
    total = context.add(total, context.multiply(context.output->record_counts[187],
                                                sizeof(std::optional<schema::InitialPointBeam>)));
    total = context.add(total, context.multiply(context.output->record_counts[188],
                                                sizeof(std::optional<schema::InitialPointRigid>)));
    total = context.add(total, context.multiply(context.output->record_counts[189],
                                                sizeof(std::optional<schema::Rational>)));
    total = context.add(total, context.multiply(context.output->record_counts[190],
                                                sizeof(std::optional<schema::FileReference>)));
    return total;
}
void schema_prepare_records(SchemaContext& context) {
    auto& owner = *context.output;
    owner.records_Case =
        context.allocate<std::optional<schema::Case>>(owner.record_counts[0], false);
    owner.records_Units =
        context.allocate<std::optional<schema::Units>>(owner.record_counts[1], false);
    owner.records_Model =
        context.allocate<std::optional<schema::Model>>(owner.record_counts[2], false);
    owner.records_Provenance =
        context.allocate<std::optional<schema::Provenance>>(owner.record_counts[3], false);
    owner.records_Frame =
        context.allocate<std::optional<schema::Frame>>(owner.record_counts[4], false);
    owner.records_Target =
        context.allocate<std::optional<schema::Target>>(owner.record_counts[5], false);
    owner.records_MaterialPoint =
        context.allocate<std::optional<schema::MaterialPoint>>(owner.record_counts[6], false);
    owner.records_Region =
        context.allocate<std::optional<schema::Region>>(owner.record_counts[7], false);
    owner.records_PolyTerm =
        context.allocate<std::optional<schema::PolyTerm>>(owner.record_counts[8], false);
    owner.records_ReferenceFieldConstant =
        context.allocate<std::optional<schema::ReferenceFieldConstant>>(owner.record_counts[9],
                                                                        false);
    owner.records_ReferenceFieldPolynomial =
        context.allocate<std::optional<schema::ReferenceFieldPolynomial>>(owner.record_counts[10],
                                                                          false);
    owner.records_TimeSeriesConstant =
        context.allocate<std::optional<schema::TimeSeriesConstant>>(owner.record_counts[11], false);
    owner.records_TimeSeriesPolynomial =
        context.allocate<std::optional<schema::TimeSeriesPolynomial>>(owner.record_counts[12],
                                                                      false);
    owner.records_TimeSeriesPchip =
        context.allocate<std::optional<schema::TimeSeriesPchip>>(owner.record_counts[13], false);
    owner.records_OutsideRejectHold =
        context.allocate<std::optional<schema::OutsideRejectHold>>(owner.record_counts[14], false);
    owner.records_OutsidePeriodic =
        context.allocate<std::optional<schema::OutsidePeriodic>>(owner.record_counts[15], false);
    owner.records_Motion =
        context.allocate<std::optional<schema::Motion>>(owner.record_counts[16], false);
    owner.records_AxisProjected =
        context.allocate<std::optional<schema::AxisProjected>>(owner.record_counts[17], false);
    owner.records_AxisExplicit =
        context.allocate<std::optional<schema::AxisExplicit>>(owner.record_counts[18], false);
    owner.records_GeometryPg =
        context.allocate<std::optional<schema::GeometryPg>>(owner.record_counts[19], false);
    owner.records_PgTerm =
        context.allocate<std::optional<schema::PgTerm>>(owner.record_counts[20], false);
    owner.records_Thickness =
        context.allocate<std::optional<schema::Thickness>>(owner.record_counts[21], false);
    owner.records_Boundary =
        context.allocate<std::optional<schema::Boundary>>(owner.record_counts[22], false);
    owner.records_Continuity =
        context.allocate<std::optional<schema::Continuity>>(owner.record_counts[23], false);
    owner.records_GeometrySurfaceGraph =
        context.allocate<std::optional<schema::GeometrySurfaceGraph>>(owner.record_counts[24],
                                                                      false);
    owner.records_GeometryNastran =
        context.allocate<std::optional<schema::GeometryNastran>>(owner.record_counts[25], false);
    owner.records_GeometryGmsh =
        context.allocate<std::optional<schema::GeometryGmsh>>(owner.record_counts[26], false);
    owner.records_ReferenceGeometrySupplied =
        context.allocate<std::optional<schema::ReferenceGeometrySupplied>>(owner.record_counts[27],
                                                                           false);
    owner.records_ReferenceGeometryFit =
        context.allocate<std::optional<schema::ReferenceGeometryFit>>(owner.record_counts[28],
                                                                      false);
    owner.records_MaterialIsotropic =
        context.allocate<std::optional<schema::MaterialIsotropic>>(owner.record_counts[29], false);
    owner.records_MaterialOrthotropic =
        context.allocate<std::optional<schema::MaterialOrthotropic>>(owner.record_counts[30],
                                                                     false);
    owner.records_MaterialEffective =
        context.allocate<std::optional<schema::MaterialEffective>>(owner.record_counts[31], false);
    owner.records_ApplicabilityUnevaluated =
        context.allocate<std::optional<schema::ApplicabilityUnevaluated>>(owner.record_counts[32],
                                                                          false);
    owner.records_ApplicabilityDeclared =
        context.allocate<std::optional<schema::ApplicabilityDeclared>>(owner.record_counts[33],
                                                                       false);
    owner.records_Layer =
        context.allocate<std::optional<schema::Layer>>(owner.record_counts[34], false);
    owner.records_Shell =
        context.allocate<std::optional<schema::Shell>>(owner.record_counts[35], false);
    owner.records_Beam =
        context.allocate<std::optional<schema::Beam>>(owner.record_counts[36], false);
    owner.records_BeamCoordinateArclength =
        context.allocate<std::optional<schema::BeamCoordinateArclength>>(owner.record_counts[37],
                                                                         false);
    owner.records_BeamCoordinateScaled =
        context.allocate<std::optional<schema::BeamCoordinateScaled>>(owner.record_counts[38],
                                                                      false);
    owner.records_SectionMatrix =
        context.allocate<std::optional<schema::SectionMatrix>>(owner.record_counts[39], false);
    owner.records_SectionCircle =
        context.allocate<std::optional<schema::SectionCircle>>(owner.record_counts[40], false);
    owner.records_SectionRectangle =
        context.allocate<std::optional<schema::SectionRectangle>>(owner.record_counts[41], false);
    owner.records_JointAdhesive =
        context.allocate<std::optional<schema::JointAdhesive>>(owner.record_counts[42], false);
    owner.records_TraceMap =
        context.allocate<std::optional<schema::TraceMap>>(owner.record_counts[43], false);
    owner.records_MapTerm =
        context.allocate<std::optional<schema::MapTerm>>(owner.record_counts[44], false);
    owner.records_JointRigidAttachment =
        context.allocate<std::optional<schema::JointRigidAttachment>>(owner.record_counts[45],
                                                                      false);
    owner.records_JointWeld =
        context.allocate<std::optional<schema::JointWeld>>(owner.record_counts[46], false);
    owner.records_SupportConstraint =
        context.allocate<std::optional<schema::SupportConstraint>>(owner.record_counts[47], false);
    owner.records_SupportTraceConstraint =
        context.allocate<std::optional<schema::SupportTraceConstraint>>(owner.record_counts[48],
                                                                        false);
    owner.records_SupportElastic =
        context.allocate<std::optional<schema::SupportElastic>>(owner.record_counts[49], false);
    owner.records_MassParticle =
        context.allocate<std::optional<schema::MassParticle>>(owner.record_counts[50], false);
    owner.records_MassRigid =
        context.allocate<std::optional<schema::MassRigid>>(owner.record_counts[51], false);
    owner.records_LedgerAdd =
        context.allocate<std::optional<schema::LedgerAdd>>(owner.record_counts[52], false);
    owner.records_LedgerReplace =
        context.allocate<std::optional<schema::LedgerReplace>>(owner.record_counts[53], false);
    owner.records_LoadPointForce =
        context.allocate<std::optional<schema::LoadPointForce>>(owner.record_counts[54], false);
    owner.records_LoadCouple =
        context.allocate<std::optional<schema::LoadCouple>>(owner.record_counts[55], false);
    owner.records_LoadGravity =
        context.allocate<std::optional<schema::LoadGravity>>(owner.record_counts[56], false);
    owner.records_LoadTraction =
        context.allocate<std::optional<schema::LoadTraction>>(owner.record_counts[57], false);
    owner.records_LoadPressure =
        context.allocate<std::optional<schema::LoadPressure>>(owner.record_counts[58], false);
    owner.records_LoadBody =
        context.allocate<std::optional<schema::LoadBody>>(owner.record_counts[59], false);
    owner.records_AeroSurface =
        context.allocate<std::optional<schema::AeroSurface>>(owner.record_counts[60], false);
    owner.records_AeroAttachmentFixed =
        context.allocate<std::optional<schema::AeroAttachmentFixed>>(owner.record_counts[61],
                                                                     false);
    owner.records_AeroAttachmentPrescribed =
        context.allocate<std::optional<schema::AeroAttachmentPrescribed>>(owner.record_counts[62],
                                                                          false);
    owner.records_AeroAttachmentStructural =
        context.allocate<std::optional<schema::AeroAttachmentStructural>>(owner.record_counts[63],
                                                                          false);
    owner.records_AeroSharing =
        context.allocate<std::optional<schema::AeroSharing>>(owner.record_counts[64], false);
    owner.records_FlowNone =
        context.allocate<std::optional<schema::FlowNone>>(owner.record_counts[65], false);
    owner.records_FlowKnown =
        context.allocate<std::optional<schema::FlowKnown>>(owner.record_counts[66], false);
    owner.records_GustNone =
        context.allocate<std::optional<schema::GustNone>>(owner.record_counts[67], false);
    owner.records_GustSine =
        context.allocate<std::optional<schema::GustSine>>(owner.record_counts[68], false);
    owner.records_FlowReference =
        context.allocate<std::optional<schema::FlowReference>>(owner.record_counts[69], false);
    owner.records_WakeNone =
        context.allocate<std::optional<schema::WakeNone>>(owner.record_counts[70], false);
    owner.records_WakeActive =
        context.allocate<std::optional<schema::WakeActive>>(owner.record_counts[71], false);
    owner.records_ConvectionFree =
        context.allocate<std::optional<schema::ConvectionFree>>(owner.record_counts[72], false);
    owner.records_ConvectionPrescribed =
        context.allocate<std::optional<schema::ConvectionPrescribed>>(owner.record_counts[73],
                                                                      false);
    owner.records_RetentionAll =
        context.allocate<std::optional<schema::RetentionAll>>(owner.record_counts[74], false);
    owner.records_RetentionWindow =
        context.allocate<std::optional<schema::RetentionWindow>>(owner.record_counts[75], false);
    owner.records_AnalysisStatic =
        context.allocate<std::optional<schema::AnalysisStatic>>(owner.record_counts[76], false);
    owner.records_AnalysisTransient =
        context.allocate<std::optional<schema::AnalysisTransient>>(owner.record_counts[77], false);
    owner.records_AnalysisFlutter =
        context.allocate<std::optional<schema::AnalysisFlutter>>(owner.record_counts[78], false);
    owner.records_AnalysisSpeedNone =
        context.allocate<std::optional<schema::AnalysisSpeedNone>>(owner.record_counts[79], false);
    owner.records_AnalysisSpeedValue =
        context.allocate<std::optional<schema::AnalysisSpeedValue>>(owner.record_counts[80], false);
    owner.records_Branch =
        context.allocate<std::optional<schema::Branch>>(owner.record_counts[81], false);
    owner.records_BranchSeedReference =
        context.allocate<std::optional<schema::BranchSeedReference>>(owner.record_counts[82],
                                                                     false);
    owner.records_BranchSeedPrepared =
        context.allocate<std::optional<schema::BranchSeedPrepared>>(owner.record_counts[83], false);
    owner.records_ContinuationFixed =
        context.allocate<std::optional<schema::ContinuationFixed>>(owner.record_counts[84], false);
    owner.records_ContinuationSweep =
        context.allocate<std::optional<schema::ContinuationSweep>>(owner.record_counts[85], false);
    owner.records_ArcLength =
        context.allocate<std::optional<schema::ArcLength>>(owner.record_counts[86], false);
    owner.records_TimeInterval =
        context.allocate<std::optional<schema::TimeInterval>>(owner.record_counts[87], false);
    owner.records_SpeedInterval =
        context.allocate<std::optional<schema::SpeedInterval>>(owner.record_counts[88], false);
    owner.records_Search =
        context.allocate<std::optional<schema::Search>>(owner.record_counts[89], false);
    owner.records_Association =
        context.allocate<std::optional<schema::Association>>(owner.record_counts[90], false);
    owner.records_Coverage =
        context.allocate<std::optional<schema::Coverage>>(owner.record_counts[91], false);
    owner.records_InitialEquilibrium =
        context.allocate<std::optional<schema::InitialEquilibrium>>(owner.record_counts[92], false);
    owner.records_InitialCold =
        context.allocate<std::optional<schema::InitialCold>>(owner.record_counts[93], false);
    owner.records_InitialPrepared =
        context.allocate<std::optional<schema::InitialPrepared>>(owner.record_counts[94], false);
    owner.records_PerturbationNone =
        context.allocate<std::optional<schema::PerturbationNone>>(owner.record_counts[95], false);
    owner.records_InitialState =
        context.allocate<std::optional<schema::InitialState>>(owner.record_counts[96], false);
    owner.records_CoefficientState =
        context.allocate<std::optional<schema::CoefficientState>>(owner.record_counts[97], false);
    owner.records_PlanesBinary64 =
        context.allocate<std::optional<schema::PlanesBinary64>>(owner.record_counts[98], false);
    owner.records_PlanesTwoWord =
        context.allocate<std::optional<schema::PlanesTwoWord>>(owner.record_counts[99], false);
    owner.records_BeamState =
        context.allocate<std::optional<schema::BeamState>>(owner.record_counts[100], false);
    owner.records_RigidState =
        context.allocate<std::optional<schema::RigidState>>(owner.record_counts[101], false);
    owner.records_FrameStatePayload =
        context.allocate<std::optional<schema::FrameStatePayload>>(owner.record_counts[102], false);
    owner.records_InitialWake =
        context.allocate<std::optional<schema::InitialWake>>(owner.record_counts[103], false);
    owner.records_SurfaceWake =
        context.allocate<std::optional<schema::SurfaceWake>>(owner.record_counts[104], false);
    owner.records_Restart =
        context.allocate<std::optional<schema::Restart>>(owner.record_counts[105], false);
    owner.records_ObservationPoint =
        context.allocate<std::optional<schema::ObservationPoint>>(owner.record_counts[106], false);
    owner.records_ObservationStress =
        context.allocate<std::optional<schema::ObservationStress>>(owner.record_counts[107], false);
    owner.records_ObservationSection = context.allocate<std::optional<schema::ObservationSection>>(
        owner.record_counts[108], false);
    owner.records_ObservationReaction =
        context.allocate<std::optional<schema::ObservationReaction>>(owner.record_counts[109],
                                                                     false);
    owner.records_ObservationAero =
        context.allocate<std::optional<schema::ObservationAero>>(owner.record_counts[110], false);
    owner.records_ObservationMode =
        context.allocate<std::optional<schema::ObservationMode>>(owner.record_counts[111], false);
    owner.records_BasisWorldMaterial = context.allocate<std::optional<schema::BasisWorldMaterial>>(
        owner.record_counts[112], false);
    owner.records_BasisExplicit =
        context.allocate<std::optional<schema::BasisExplicit>>(owner.record_counts[113], false);
    owner.records_Comparison =
        context.allocate<std::optional<schema::Comparison>>(owner.record_counts[114], false);
    owner.records_ComparisonReferenceNone =
        context.allocate<std::optional<schema::ComparisonReferenceNone>>(owner.record_counts[115],
                                                                         false);
    owner.records_ComparisonReferenceFile =
        context.allocate<std::optional<schema::ComparisonReferenceFile>>(owner.record_counts[116],
                                                                         false);
    owner.records_ObservationBudget =
        context.allocate<std::optional<schema::ObservationBudget>>(owner.record_counts[117], false);
    owner.records_ErrorGroupNotApplicable =
        context.allocate<std::optional<schema::ErrorGroupNotApplicable>>(owner.record_counts[118],
                                                                         false);
    owner.records_ErrorGroupCoveredBy =
        context.allocate<std::optional<schema::ErrorGroupCoveredBy>>(owner.record_counts[119],
                                                                     false);
    owner.records_ErrorGroupRequired = context.allocate<std::optional<schema::ErrorGroupRequired>>(
        owner.record_counts[120], false);
    owner.records_Numerics =
        context.allocate<std::optional<schema::Numerics>>(owner.record_counts[121], false);
    owner.records_WakeCapacityNone =
        context.allocate<std::optional<schema::WakeCapacityNone>>(owner.record_counts[122], false);
    owner.records_WakeCapacityLimited =
        context.allocate<std::optional<schema::WakeCapacityLimited>>(owner.record_counts[123],
                                                                     false);
    owner.records_PatchPolicyNone =
        context.allocate<std::optional<schema::PatchPolicyNone>>(owner.record_counts[124], false);
    owner.records_PatchPolicyActive =
        context.allocate<std::optional<schema::PatchPolicyActive>>(owner.record_counts[125], false);
    owner.records_Quadrature =
        context.allocate<std::optional<schema::Quadrature>>(owner.record_counts[126], false);
    owner.records_StructuralQuadratureNone =
        context.allocate<std::optional<schema::StructuralQuadratureNone>>(owner.record_counts[127],
                                                                          false);
    owner.records_StructuralQuadratureActive =
        context.allocate<std::optional<schema::StructuralQuadratureActive>>(
            owner.record_counts[128], false);
    owner.records_AeroQuadratureNone = context.allocate<std::optional<schema::AeroQuadratureNone>>(
        owner.record_counts[129], false);
    owner.records_AeroQuadratureActive =
        context.allocate<std::optional<schema::AeroQuadratureActive>>(owner.record_counts[130],
                                                                      false);
    owner.records_ScaleUniformComponents =
        context.allocate<std::optional<schema::ScaleUniformComponents>>(owner.record_counts[131],
                                                                        false);
    owner.records_ScalePhysicalLayoutPolicy =
        context.allocate<std::optional<schema::ScalePhysicalLayoutPolicy>>(owner.record_counts[132],
                                                                           false);
    owner.records_ConstraintScalePolicyNone =
        context.allocate<std::optional<schema::ConstraintScalePolicyNone>>(owner.record_counts[133],
                                                                           false);
    owner.records_ConstraintScalePolicyPhysical =
        context.allocate<std::optional<schema::ConstraintScalePolicyPhysical>>(
            owner.record_counts[134], false);
    owner.records_ComponentScale3None =
        context.allocate<std::optional<schema::ComponentScale3None>>(owner.record_counts[135],
                                                                     false);
    owner.records_ComponentScale3Values =
        context.allocate<std::optional<schema::ComponentScale3Values>>(owner.record_counts[136],
                                                                       false);
    owner.records_ScaleExplicitLayout =
        context.allocate<std::optional<schema::ScaleExplicitLayout>>(owner.record_counts[137],
                                                                     false);
    owner.records_LinearPolicyDirect = context.allocate<std::optional<schema::LinearPolicyDirect>>(
        owner.record_counts[138], false);
    owner.records_LinearPolicyGmres =
        context.allocate<std::optional<schema::LinearPolicyGmres>>(owner.record_counts[139], false);
    owner.records_LinearAccuracy =
        context.allocate<std::optional<schema::LinearAccuracy>>(owner.record_counts[140], false);
    owner.records_NewtonPolicy =
        context.allocate<std::optional<schema::NewtonPolicy>>(owner.record_counts[141], false);
    owner.records_PrecisionPolicy =
        context.allocate<std::optional<schema::PrecisionPolicy>>(owner.record_counts[142], false);
    owner.records_BackendMemory =
        context.allocate<std::optional<schema::BackendMemory>>(owner.record_counts[143], false);
    owner.records_EigenPolicyNone =
        context.allocate<std::optional<schema::EigenPolicyNone>>(owner.record_counts[144], false);
    owner.records_EigenPolicyFull =
        context.allocate<std::optional<schema::EigenPolicyFull>>(owner.record_counts[145], false);
    owner.records_EigenPolicyGeneralized =
        context.allocate<std::optional<schema::EigenPolicyGeneralized>>(owner.record_counts[146],
                                                                        false);
    owner.records_EigenPolicyArnoldi = context.allocate<std::optional<schema::EigenPolicyArnoldi>>(
        owner.record_counts[147], false);
    owner.records_Complex =
        context.allocate<std::optional<schema::Complex>>(owner.record_counts[148], false);
    owner.records_ArnoldiSeedUniform = context.allocate<std::optional<schema::ArnoldiSeedUniform>>(
        owner.record_counts[149], false);
    owner.records_ArnoldiSeedExplicit =
        context.allocate<std::optional<schema::ArnoldiSeedExplicit>>(owner.record_counts[150],
                                                                     false);
    owner.records_EigenSelectorLargest =
        context.allocate<std::optional<schema::EigenSelectorLargest>>(owner.record_counts[151],
                                                                      false);
    owner.records_EigenSelectorShift = context.allocate<std::optional<schema::EigenSelectorShift>>(
        owner.record_counts[152], false);
    owner.records_EigenRequirements =
        context.allocate<std::optional<schema::EigenRequirements>>(owner.record_counts[153], false);
    owner.records_BackwardNone =
        context.allocate<std::optional<schema::BackwardNone>>(owner.record_counts[154], false);
    owner.records_BackwardRequired =
        context.allocate<std::optional<schema::BackwardRequired>>(owner.record_counts[155], false);
    owner.records_EigenCluster =
        context.allocate<std::optional<schema::EigenCluster>>(owner.record_counts[156], false);
    owner.records_EigenLimits =
        context.allocate<std::optional<schema::EigenLimits>>(owner.record_counts[157], false);
    owner.records_Resources =
        context.allocate<std::optional<schema::Resources>>(owner.record_counts[158], false);
    owner.records_InputQuotas =
        context.allocate<std::optional<schema::InputQuotas>>(owner.record_counts[159], false);
    owner.records_Output =
        context.allocate<std::optional<schema::Output>>(owner.record_counts[160], false);
    owner.records_OutputCadence =
        context.allocate<std::optional<schema::OutputCadence>>(owner.record_counts[161], false);
    owner.records_Visualization =
        context.allocate<std::optional<schema::Visualization>>(owner.record_counts[162], false);
    owner.records_CheckpointOutputOff =
        context.allocate<std::optional<schema::CheckpointOutputOff>>(owner.record_counts[163],
                                                                     false);
    owner.records_CheckpointOutputOn = context.allocate<std::optional<schema::CheckpointOutputOn>>(
        owner.record_counts[164], false);
    owner.records_ReferenceFieldExpression =
        context.allocate<std::optional<schema::ReferenceFieldExpression>>(owner.record_counts[165],
                                                                          false);
    owner.records_TimeSeriesExpression =
        context.allocate<std::optional<schema::TimeSeriesExpression>>(owner.record_counts[166],
                                                                      false);
    owner.records_LoadCurrentExpression =
        context.allocate<std::optional<schema::LoadCurrentExpression>>(owner.record_counts[167],
                                                                       false);
    owner.records_Expression =
        context.allocate<std::optional<schema::Expression>>(owner.record_counts[168], false);
    owner.records_ExpressionDomainReference =
        context.allocate<std::optional<schema::ExpressionDomainReference>>(owner.record_counts[169],
                                                                           false);
    owner.records_ExpressionDomainTime =
        context.allocate<std::optional<schema::ExpressionDomainTime>>(owner.record_counts[170],
                                                                      false);
    owner.records_ExpressionDomainCurrentPositionTime =
        context.allocate<std::optional<schema::ExpressionDomainCurrentPositionTime>>(
            owner.record_counts[171], false);
    owner.records_ExprNodeConstant =
        context.allocate<std::optional<schema::ExprNodeConstant>>(owner.record_counts[172], false);
    owner.records_ExprNodeVariable =
        context.allocate<std::optional<schema::ExprNodeVariable>>(owner.record_counts[173], false);
    owner.records_ExprNodeBinary =
        context.allocate<std::optional<schema::ExprNodeBinary>>(owner.record_counts[174], false);
    owner.records_ExprNodeNegate =
        context.allocate<std::optional<schema::ExprNodeNegate>>(owner.record_counts[175], false);
    owner.records_ExprNodePower =
        context.allocate<std::optional<schema::ExprNodePower>>(owner.record_counts[176], false);
    owner.records_ExprNodeAnalytic =
        context.allocate<std::optional<schema::ExprNodeAnalytic>>(owner.record_counts[177], false);
    owner.records_PerturbationPhysical =
        context.allocate<std::optional<schema::PerturbationPhysical>>(owner.record_counts[178],
                                                                      false);
    owner.records_ShellPhysicalState = context.allocate<std::optional<schema::ShellPhysicalState>>(
        owner.record_counts[179], false);
    owner.records_FullShellFieldGlobalPolynomial =
        context.allocate<std::optional<schema::FullShellFieldGlobalPolynomial>>(
            owner.record_counts[180], false);
    owner.records_FullShellFieldHalfCells =
        context.allocate<std::optional<schema::FullShellFieldHalfCells>>(owner.record_counts[181],
                                                                         false);
    owner.records_CellField =
        context.allocate<std::optional<schema::CellField>>(owner.record_counts[182], false);
    owner.records_VectorPolynomial =
        context.allocate<std::optional<schema::VectorPolynomial>>(owner.record_counts[183], false);
    owner.records_VectorTerm =
        context.allocate<std::optional<schema::VectorTerm>>(owner.record_counts[184], false);
    owner.records_PointCheck =
        context.allocate<std::optional<schema::PointCheck>>(owner.record_counts[185], false);
    owner.records_InitialPointShell =
        context.allocate<std::optional<schema::InitialPointShell>>(owner.record_counts[186], false);
    owner.records_InitialPointBeam =
        context.allocate<std::optional<schema::InitialPointBeam>>(owner.record_counts[187], false);
    owner.records_InitialPointRigid =
        context.allocate<std::optional<schema::InitialPointRigid>>(owner.record_counts[188], false);
    owner.records_Rational =
        context.allocate<std::optional<schema::Rational>>(owner.record_counts[189], false);
    owner.records_FileReference =
        context.allocate<std::optional<schema::FileReference>>(owner.record_counts[190], false);
}
void schema_build_records(SchemaContext& context) {
    auto& owner = *context.output;
    for (std::uint64_t node = 0; node < owner.nodes.size(); ++node) {
        const auto& typed = owner.nodes.values()[static_cast<std::size_t>(node)];
        if (!typed.assigned || node_kind(owner.snapshot.values()[0], node) != InputNodeKind::Object)
            continue;
        context.work();
        switch (typed.record_kind) {
        case 0:
            owner.records_Case.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Case{
                    SchemaFactory::required_field<FieldType::T0>(context, node, "schema"),
                    SchemaFactory::required_field<FieldType::T1>(context, node, "units"),
                    SchemaFactory::required_field<FieldType::T2>(context, node, "model"),
                    SchemaFactory::required_field<FieldType::T3>(context, node, "flow"),
                    SchemaFactory::required_field<FieldType::T4>(context, node, "analysis"),
                    SchemaFactory::required_field<FieldType::T5>(context, node, "observations"),
                    SchemaFactory::required_field<FieldType::T6>(context, node, "numerics"),
                    SchemaFactory::required_field<FieldType::T7>(context, node, "resources"),
                    SchemaFactory::required_field<FieldType::T8>(context, node, "output"),
                    SchemaFactory::optional_field<FieldType::T9>(context, node, "initial"),
                    SchemaFactory::optional_field<FieldType::T10>(context, node, "restart"),
                    SchemaFactory::optional_field<FieldType::T11>(context, node, "provenance"),
                });
            break;
        case 1:
            owner.records_Units.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Units{
                    SchemaFactory::required_field<FieldType::T12>(context, node, "length"),
                    SchemaFactory::required_field<FieldType::T13>(context, node, "mass"),
                    SchemaFactory::required_field<FieldType::T14>(context, node, "time"),
                    SchemaFactory::required_field<FieldType::T15>(context, node, "angle"),
                });
            break;
        case 2:
            owner.records_Model.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Model{
                    SchemaFactory::required_field<FieldType::T16>(context, node, "geometry"),
                    SchemaFactory::required_field<FieldType::T17>(context, node, "materials"),
                    SchemaFactory::required_field<FieldType::T18>(context, node, "shells"),
                    SchemaFactory::required_field<FieldType::T19>(context, node, "beams"),
                    SchemaFactory::required_field<FieldType::T20>(context, node, "joints"),
                    SchemaFactory::required_field<FieldType::T21>(context, node, "supports"),
                    SchemaFactory::required_field<FieldType::T22>(context, node, "masses"),
                    SchemaFactory::required_field<FieldType::T23>(context, node, "aero_surfaces"),
                });
            break;
        case 3:
            owner.records_Provenance.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Provenance{
                    SchemaFactory::required_field<FieldType::T24>(context, node, "description"),
                    SchemaFactory::required_field<FieldType::T25>(context, node, "references"),
                    SchemaFactory::required_field<FieldType::T25>(context, node, "assumptions"),
                });
            break;
        case 4:
            owner.records_Frame.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Frame{
                    SchemaFactory::required_field<FieldType::T26>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T27>(context, node, "rotation"),
                });
            break;
        case 5:
            owner.records_Target.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Target{
                    SchemaFactory::required_field<FieldType::T28>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T29>(context, node, "id"),
                });
            break;
        case 6:
            owner.records_MaterialPoint.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::MaterialPoint{
                    SchemaFactory::required_field<FieldType::T30>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "coordinates"),
                    SchemaFactory::required_field<FieldType::T32>(context, node, "side"),
                });
            break;
        case 7:
            owner.records_Region.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Region{
                    SchemaFactory::required_field<FieldType::T30>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T33>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T34>(context, node, "geometry"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "boundary"),
                    SchemaFactory::required_field<FieldType::T32>(context, node, "side"),
                });
            break;
        case 8:
            owner.records_PolyTerm.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::PolyTerm{
                    SchemaFactory::required_field<FieldType::T36>(context, node, "powers"),
                    SchemaFactory::required_field<FieldType::T37>(context, node, "coefficient"),
                });
            break;
        case 9:
            owner.records_ReferenceFieldConstant
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ReferenceFieldConstant{
                    SchemaFactory::required_field<FieldType::T38>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T39>(context, node, "value"),
                });
            break;
        case 10:
            owner.records_ReferenceFieldPolynomial
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ReferenceFieldPolynomial{
                    SchemaFactory::required_field<FieldType::T40>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T41>(context, node, "coordinates"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "lengths"),
                    SchemaFactory::required_field<FieldType::T42>(context, node, "components"),
                    SchemaFactory::required_field<FieldType::T39>(context, node, "error"),
                });
            break;
        case 11:
            owner.records_TimeSeriesConstant.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::TimeSeriesConstant{
                    SchemaFactory::required_field<FieldType::T38>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T43>(context, node, "value"),
                });
            break;
        case 12:
            owner.records_TimeSeriesPolynomial.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::TimeSeriesPolynomial{
                    SchemaFactory::required_field<FieldType::T40>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "scale"),
                    SchemaFactory::required_field<FieldType::T45>(context, node, "coefficients"),
                    SchemaFactory::required_field<FieldType::T46>(context, node, "domain"),
                });
            break;
        case 13:
            owner.records_TimeSeriesPchip.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::TimeSeriesPchip{
                    SchemaFactory::required_field<FieldType::T47>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T48>(context, node, "times"),
                    SchemaFactory::required_field<FieldType::T45>(context, node, "values"),
                    SchemaFactory::required_field<FieldType::T49>(context, node, "outside"),
                });
            break;
        case 14:
            owner.records_OutsideRejectHold.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::OutsideRejectHold{
                    SchemaFactory::required_field<FieldType::T50>(context, node, "kind"),
                });
            break;
        case 15:
            owner.records_OutsidePeriodic.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::OutsidePeriodic{
                    SchemaFactory::required_field<FieldType::T51>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "period"),
                });
            break;
        case 16:
            owner.records_Motion.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Motion{
                    SchemaFactory::required_field<FieldType::T52>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T53>(context, node, "translation"),
                    SchemaFactory::required_field<FieldType::T53>(context, node, "rotation_vector"),
                    SchemaFactory::required_field<FieldType::T48>(context, node, "events"),
                });
            break;
        case 17:
            owner.records_AxisProjected.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AxisProjected{
                    SchemaFactory::required_field<FieldType::T54>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T55>(context, node, "direction"),
                    SchemaFactory::required_field<FieldType::T56>(context, node, "angle"),
                });
            break;
        case 18:
            owner.records_AxisExplicit.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::AxisExplicit{
                    SchemaFactory::required_field<FieldType::T57>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T58>(context, node, "columns"),
                });
            break;
        case 19:
            owner.records_GeometryPg.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::GeometryPg{
                    SchemaFactory::required_field<FieldType::T59>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T60>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T52>(context, node, "frame"),
                    SchemaFactory::required_field<FieldType::T61>(context, node,
                                                                  "reference_length"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "x0"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "width"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "skew"),
                    SchemaFactory::required_field<FieldType::T63>(context, node, "degree"),
                    SchemaFactory::required_field<FieldType::T64>(context, node, "coefficients"),
                    SchemaFactory::required_field<FieldType::T65>(context, node, "domain"),
                    SchemaFactory::required_field<FieldType::T66>(context, node, "thickness"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "n"),
                    SchemaFactory::required_field<FieldType::T68>(context, node, "convention"),
                    SchemaFactory::required_field<FieldType::T69>(context, node, "boundaries"),
                    SchemaFactory::required_field<FieldType::T70>(context, node, "continuity"),
                });
            break;
        case 20:
            owner.records_PgTerm.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::PgTerm{
                    SchemaFactory::required_field<FieldType::T63>(context, node, "powers"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "coefficient"),
                });
            break;
        case 21:
            owner.records_Thickness.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Thickness{
                    SchemaFactory::required_field<FieldType::T65>(context, node, "vertices"),
                    SchemaFactory::required_field<FieldType::T71>(context, node, "triangles"),
                    SchemaFactory::required_field<FieldType::T72>(context, node, "values"),
                });
            break;
        case 22:
            owner.records_Boundary.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Boundary{
                    SchemaFactory::required_field<FieldType::T24>(context, node, "name"),
                    SchemaFactory::required_field<FieldType::T65>(context, node, "vertices"),
                });
            break;
        case 23:
            owner.records_Continuity.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Continuity{
                    SchemaFactory::required_field<FieldType::T34>(context, node, "other"),
                    SchemaFactory::required_field<FieldType::T73>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "own_boundary"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "other_boundary"),
                    SchemaFactory::required_field<FieldType::T74>(context, node, "orientation"),
                });
            break;
        case 24:
            owner.records_GeometrySurfaceGraph.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::GeometrySurfaceGraph{
                    SchemaFactory::required_field<FieldType::T75>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T60>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T52>(context, node, "frame"),
                    SchemaFactory::required_field<FieldType::T61>(context, node,
                                                                  "reference_length"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "x0"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "width"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "skew"),
                    SchemaFactory::required_field<FieldType::T63>(context, node, "degree"),
                    SchemaFactory::required_field<FieldType::T64>(context, node, "coefficients"),
                    SchemaFactory::required_field<FieldType::T65>(context, node, "domain"),
                    SchemaFactory::required_field<FieldType::T69>(context, node, "boundaries"),
                    SchemaFactory::required_field<FieldType::T70>(context, node, "continuity"),
                });
            break;
        case 25:
            owner.records_GeometryNastran.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::GeometryNastran{
                    SchemaFactory::required_field<FieldType::T76>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T60>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T77>(context, node, "deck"),
                    SchemaFactory::required_field<FieldType::T78>(context, node, "sidecar"),
                    SchemaFactory::required_field<FieldType::T79>(context, node, "dialect"),
                    SchemaFactory::required_field<FieldType::T80>(context, node, "mode"),
                    SchemaFactory::required_field<FieldType::T25>(context, node, "excluded_cards"),
                });
            break;
        case 26:
            owner.records_GeometryGmsh.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::GeometryGmsh{
                    SchemaFactory::required_field<FieldType::T81>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T60>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T82>(context, node, "mesh"),
                    SchemaFactory::required_field<FieldType::T83>(context, node, "sidecar"),
                    SchemaFactory::required_field<FieldType::T84>(context, node, "reference"),
                });
            break;
        case 27:
            owner.records_ReferenceGeometrySupplied
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ReferenceGeometrySupplied{
                    SchemaFactory::required_field<FieldType::T85>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T86>(context, node, "file"),
                });
            break;
        case 28:
            owner.records_ReferenceGeometryFit.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ReferenceGeometryFit{
                    SchemaFactory::required_field<FieldType::T87>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "degree"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "distance_budget"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "normal_budget"),
                    SchemaFactory::required_field<FieldType::T89>(context, node,
                                                                  "conversion_observations"),
                });
            break;
        case 29:
            owner.records_MaterialIsotropic.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::MaterialIsotropic{
                    SchemaFactory::required_field<FieldType::T90>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T91>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T92>(context, node, "young"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "poisson"),
                    SchemaFactory::required_field<FieldType::T93>(context, node, "density"),
                    SchemaFactory::required_field<FieldType::T11>(context, node, "provenance"),
                    SchemaFactory::required_field<FieldType::T94>(context, node, "applicability"),
                });
            break;
        case 30:
            owner.records_MaterialOrthotropic.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::MaterialOrthotropic{
                    SchemaFactory::required_field<FieldType::T95>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T91>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T96>(context, node, "young"),
                    SchemaFactory::required_field<FieldType::T96>(context, node, "shear"),
                    SchemaFactory::required_field<FieldType::T97>(context, node, "poisson"),
                    SchemaFactory::required_field<FieldType::T93>(context, node, "density"),
                    SchemaFactory::required_field<FieldType::T11>(context, node, "provenance"),
                    SchemaFactory::required_field<FieldType::T94>(context, node, "applicability"),
                });
            break;
        case 31:
            owner.records_MaterialEffective.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::MaterialEffective{
                    SchemaFactory::required_field<FieldType::T98>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T91>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T99>(context, node, "stiffness"),
                    SchemaFactory::required_field<FieldType::T100>(context, node, "components"),
                    SchemaFactory::required_field<FieldType::T93>(context, node, "density"),
                    SchemaFactory::required_field<FieldType::T11>(context, node, "provenance"),
                    SchemaFactory::required_field<FieldType::T94>(context, node, "applicability"),
                });
            break;
        case 32:
            owner.records_ApplicabilityUnevaluated
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ApplicabilityUnevaluated{
                    SchemaFactory::required_field<FieldType::T101>(context, node, "kind"),
                });
            break;
        case 33:
            owner.records_ApplicabilityDeclared.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ApplicabilityDeclared{
                    SchemaFactory::required_field<FieldType::T102>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T103>(context, node, "strain_lower"),
                    SchemaFactory::required_field<FieldType::T103>(context, node, "strain_upper"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "reference"),
                });
            break;
        case 34:
            owner.records_Layer.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Layer{
                    SchemaFactory::required_field<FieldType::T104>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T105>(context, node, "interval"),
                    SchemaFactory::required_field<FieldType::T106>(context, node, "material"),
                    SchemaFactory::required_field<FieldType::T107>(context, node, "axis"),
                    SchemaFactory::required_field<FieldType::T108>(context, node, "initial_green"),
                    SchemaFactory::required_field<FieldType::T109>(context, node, "viscosity"),
                });
            break;
        case 35:
            owner.records_Shell.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Shell{
                    SchemaFactory::required_field<FieldType::T110>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T34>(context, node, "geometry"),
                    SchemaFactory::required_field<FieldType::T111>(context, node, "layers"),
                    SchemaFactory::required_field<FieldType::T100>(context, node,
                                                                   "damping_components"),
                    SchemaFactory::required_field<FieldType::T112>(context, node, "loads"),
                });
            break;
        case 36:
            owner.records_Beam.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Beam{
                    SchemaFactory::required_field<FieldType::T113>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T114>(context, node, "nodes"),
                    SchemaFactory::required_field<FieldType::T115>(context, node, "coordinate"),
                    SchemaFactory::required_field<FieldType::T116>(context, node, "section"),
                    SchemaFactory::required_field<FieldType::T103>(context, node, "initial_strain"),
                    SchemaFactory::required_field<FieldType::T117>(context, node,
                                                                   "residual_energy"),
                    SchemaFactory::required_field<FieldType::T118>(context, node, "viscosity"),
                    SchemaFactory::required_field<FieldType::T112>(context, node, "loads"),
                });
            break;
        case 37:
            owner.records_BeamCoordinateArclength
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BeamCoordinateArclength{
                    SchemaFactory::required_field<FieldType::T119>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "length"),
                });
            break;
        case 38:
            owner.records_BeamCoordinateScaled.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BeamCoordinateScaled{
                    SchemaFactory::required_field<FieldType::T120>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "length"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "arc_scale"),
                });
            break;
        case 39:
            owner.records_SectionMatrix.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::SectionMatrix{
                    SchemaFactory::required_field<FieldType::T121>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T118>(context, node, "stiffness"),
                    SchemaFactory::required_field<FieldType::T122>(context, node, "line_density"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "center"),
                    SchemaFactory::required_field<FieldType::T123>(context, node, "inertia"),
                    SchemaFactory::required_field<FieldType::T124>(context, node, "envelope"),
                    SchemaFactory::required_field<FieldType::T125>(context, node, "warping"),
                    SchemaFactory::required_field<FieldType::T126>(context, node, "recovery"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "source"),
                });
            break;
        case 40:
            owner.records_SectionCircle.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::SectionCircle{
                    SchemaFactory::required_field<FieldType::T127>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "radius"),
                    SchemaFactory::required_field<FieldType::T92>(context, node, "young"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "poisson"),
                    SchemaFactory::required_field<FieldType::T93>(context, node, "density"),
                    SchemaFactory::required_field<FieldType::T128>(context, node,
                                                                   "effective_shear"),
                    SchemaFactory::required_field<FieldType::T125>(context, node, "warping"),
                });
            break;
        case 41:
            owner.records_SectionRectangle.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::SectionRectangle{
                    SchemaFactory::required_field<FieldType::T129>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T124>(context, node, "sides"),
                    SchemaFactory::required_field<FieldType::T92>(context, node, "young"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "poisson"),
                    SchemaFactory::required_field<FieldType::T93>(context, node, "density"),
                    SchemaFactory::required_field<FieldType::T128>(context, node,
                                                                   "effective_shear"),
                    SchemaFactory::required_field<FieldType::T130>(context, node, "torsion"),
                    SchemaFactory::required_field<FieldType::T125>(context, node, "warping"),
                });
            break;
        case 42:
            owner.records_JointAdhesive.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::JointAdhesive{
                    SchemaFactory::required_field<FieldType::T131>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T132>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "first"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "second"),
                    SchemaFactory::required_field<FieldType::T134>(context, node, "mapping"),
                    SchemaFactory::required_field<FieldType::T135>(context, node,
                                                                   "split_at_attachment"),
                });
            break;
        case 43:
            owner.records_TraceMap.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::TraceMap{
                    SchemaFactory::required_field<FieldType::T136>(context, node, "first_cell"),
                    SchemaFactory::required_field<FieldType::T136>(context, node, "second_cell"),
                    SchemaFactory::required_field<FieldType::T65>(context, node, "domain"),
                    SchemaFactory::required_field<FieldType::T137>(context, node, "first"),
                    SchemaFactory::required_field<FieldType::T137>(context, node, "second"),
                    SchemaFactory::required_field<FieldType::T74>(context, node, "orientation"),
                });
            break;
        case 44:
            owner.records_MapTerm.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::MapTerm{
                    SchemaFactory::required_field<FieldType::T63>(context, node, "powers"),
                    SchemaFactory::required_field<FieldType::T62>(context, node, "coefficient"),
                });
            break;
        case 45:
            owner.records_JointRigidAttachment.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::JointRigidAttachment{
                    SchemaFactory::required_field<FieldType::T138>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T132>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "attached"),
                    SchemaFactory::required_field<FieldType::T30>(context, node, "frame"),
                    SchemaFactory::required_field<FieldType::T52>(context, node, "offset"),
                });
            break;
        case 46:
            owner.records_JointWeld.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::JointWeld{
                    SchemaFactory::required_field<FieldType::T139>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T132>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T30>(context, node, "first"),
                    SchemaFactory::required_field<FieldType::T30>(context, node, "second"),
                    SchemaFactory::required_field<FieldType::T52>(context, node, "relative"),
                });
            break;
        case 47:
            owner.records_SupportConstraint.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::SupportConstraint{
                    SchemaFactory::required_field<FieldType::T140>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T141>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T143>(context, node, "components"),
                    SchemaFactory::required_field<FieldType::T144>(context, node, "motion"),
                });
            break;
        case 48:
            owner.records_SupportTraceConstraint
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::SupportTraceConstraint{
                    SchemaFactory::required_field<FieldType::T145>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T141>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "region"),
                    SchemaFactory::required_field<FieldType::T143>(context, node, "components"),
                    SchemaFactory::required_field<FieldType::T144>(context, node, "motion"),
                });
            break;
        case 49:
            owner.records_SupportElastic.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::SupportElastic{
                    SchemaFactory::required_field<FieldType::T146>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T141>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T144>(context, node, "reference"),
                    SchemaFactory::required_field<FieldType::T118>(context, node, "stiffness"),
                    SchemaFactory::required_field<FieldType::T118>(context, node, "damping"),
                    SchemaFactory::required_field<FieldType::T52>(context, node, "frame"),
                });
            break;
        case 50:
            owner.records_MassParticle.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::MassParticle{
                    SchemaFactory::required_field<FieldType::T147>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T148>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T149>(context, node, "mass"),
                    SchemaFactory::required_field<FieldType::T150>(context, node, "ledger"),
                });
            break;
        case 51:
            owner.records_MassRigid.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::MassRigid{
                    SchemaFactory::required_field<FieldType::T151>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T148>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T149>(context, node, "mass"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "center"),
                    SchemaFactory::required_field<FieldType::T152>(context, node, "inertia"),
                    SchemaFactory::required_field<FieldType::T52>(context, node, "frame"),
                    SchemaFactory::required_field<FieldType::T150>(context, node, "ledger"),
                });
            break;
        case 52:
            owner.records_LedgerAdd.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::LedgerAdd{
                    SchemaFactory::required_field<FieldType::T153>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "source"),
                    SchemaFactory::required_field<FieldType::T30>(context, node, "owner"),
                });
            break;
        case 53:
            owner.records_LedgerReplace.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::LedgerReplace{
                    SchemaFactory::required_field<FieldType::T154>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "source"),
                    SchemaFactory::required_field<FieldType::T30>(context, node, "owner"),
                    SchemaFactory::required_field<FieldType::T155>(context, node, "replaced"),
                });
            break;
        case 54:
            owner.records_LoadPointForce.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::LoadPointForce{
                    SchemaFactory::required_field<FieldType::T156>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T157>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T158>(context, node, "frame"),
                    SchemaFactory::required_field<FieldType::T53>(context, node, "series"),
                });
            break;
        case 55:
            owner.records_LoadCouple.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::LoadCouple{
                    SchemaFactory::required_field<FieldType::T159>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T157>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T158>(context, node, "frame"),
                    SchemaFactory::required_field<FieldType::T53>(context, node, "series"),
                });
            break;
        case 56:
            owner.records_LoadGravity.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::LoadGravity{
                    SchemaFactory::required_field<FieldType::T160>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T157>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T53>(context, node, "acceleration"),
                });
            break;
        case 57:
            owner.records_LoadTraction.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::LoadTraction{
                    SchemaFactory::required_field<FieldType::T161>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T157>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "region"),
                    SchemaFactory::required_field<FieldType::T162>(context, node, "area"),
                    SchemaFactory::required_field<FieldType::T158>(context, node, "frame"),
                    SchemaFactory::required_field<FieldType::T163>(context, node, "spatial"),
                    SchemaFactory::required_field<FieldType::T164>(context, node, "multiplier"),
                });
            break;
        case 58:
            owner.records_LoadPressure.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::LoadPressure{
                    SchemaFactory::required_field<FieldType::T165>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T157>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "region"),
                    SchemaFactory::required_field<FieldType::T166>(context, node, "spatial"),
                    SchemaFactory::required_field<FieldType::T164>(context, node, "multiplier"),
                });
            break;
        case 59:
            owner.records_LoadBody.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::LoadBody{
                    SchemaFactory::required_field<FieldType::T167>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T157>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "region"),
                    SchemaFactory::required_field<FieldType::T168>(context, node, "spatial"),
                    SchemaFactory::required_field<FieldType::T164>(context, node, "multiplier"),
                });
            break;
        case 60:
            owner.records_AeroSurface.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::AeroSurface{
                    SchemaFactory::required_field<FieldType::T169>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T34>(context, node, "geometry"),
                    SchemaFactory::required_field<FieldType::T170>(context, node, "attachment"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "chord"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "span"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "leading"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "trailing"),
                    SchemaFactory::required_field<FieldType::T172>(context, node, "orientation"),
                    SchemaFactory::required_field<FieldType::T173>(context, node, "sharing"),
                });
            break;
        case 61:
            owner.records_AeroAttachmentFixed.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AeroAttachmentFixed{
                    SchemaFactory::required_field<FieldType::T174>(context, node, "kind"),
                });
            break;
        case 62:
            owner.records_AeroAttachmentPrescribed
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AeroAttachmentPrescribed{
                    SchemaFactory::required_field<FieldType::T175>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T144>(context, node, "motion"),
                });
            break;
        case 63:
            owner.records_AeroAttachmentStructural
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AeroAttachmentStructural{
                    SchemaFactory::required_field<FieldType::T176>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T177>(context, node, "host"),
                });
            break;
        case 64:
            owner.records_AeroSharing.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::AeroSharing{
                    SchemaFactory::required_field<FieldType::T178>(context, node, "other"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "own_boundary"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "other_boundary"),
                    SchemaFactory::required_field<FieldType::T74>(context, node, "orientation"),
                });
            break;
        case 65:
            owner.records_FlowNone.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::FlowNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 66:
            owner.records_FlowKnown.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::FlowKnown{
                    SchemaFactory::required_field<FieldType::T180>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T93>(context, node, "density"),
                    SchemaFactory::required_field<FieldType::T53>(context, node, "uniform"),
                    SchemaFactory::required_field<FieldType::T181>(context, node, "gust"),
                    SchemaFactory::required_field<FieldType::T182>(context, node, "reference"),
                    SchemaFactory::required_field<FieldType::T48>(context, node, "events"),
                });
            break;
        case 67:
            owner.records_GustNone.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::GustNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 68:
            owner.records_GustSine.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::GustSine{
                    SchemaFactory::required_field<FieldType::T183>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T164>(context, node, "amplitude"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "direction"),
                    SchemaFactory::required_field<FieldType::T184>(context, node, "wavevector"),
                    SchemaFactory::required_field<FieldType::T185>(context, node,
                                                                   "angular_frequency"),
                    SchemaFactory::required_field<FieldType::T186>(context, node, "phase"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "time_origin"),
                });
            break;
        case 69:
            owner.records_FlowReference.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::FlowReference{
                    SchemaFactory::required_field<FieldType::T93>(context, node, "density"),
                    SchemaFactory::required_field<FieldType::T187>(context, node, "speed"),
                    SchemaFactory::required_field<FieldType::T188>(context, node, "area"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "chord"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "span"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "moment_origin"),
                });
            break;
        case 70:
            owner.records_WakeNone.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::WakeNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 71:
            owner.records_WakeActive.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::WakeActive{
                    SchemaFactory::required_field<FieldType::T189>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T190>(context, node, "convection"),
                    SchemaFactory::required_field<FieldType::T191>(context, node, "retention"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "core_radius"),
                    SchemaFactory::required_field<FieldType::T192>(context, node, "model"),
                    SchemaFactory::required_field<FieldType::T193>(context, node, "mapping"),
                });
            break;
        case 72:
            owner.records_ConvectionFree.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ConvectionFree{
                    SchemaFactory::required_field<FieldType::T125>(context, node, "kind"),
                });
            break;
        case 73:
            owner.records_ConvectionPrescribed.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ConvectionPrescribed{
                    SchemaFactory::required_field<FieldType::T194>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T53>(context, node, "velocity"),
                });
            break;
        case 74:
            owner.records_RetentionAll.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::RetentionAll{
                    SchemaFactory::required_field<FieldType::T195>(context, node, "kind"),
                });
            break;
        case 75:
            owner.records_RetentionWindow.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::RetentionWindow{
                    SchemaFactory::required_field<FieldType::T196>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "duration"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "distance"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "error_budget"),
                });
            break;
        case 76:
            owner.records_AnalysisStatic.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AnalysisStatic{
                    SchemaFactory::required_field<FieldType::T197>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T198>(context, node, "physics"),
                    SchemaFactory::required_field<FieldType::T199>(context, node, "branch"),
                    SchemaFactory::required_field<FieldType::T200>(context, node, "speed"),
                    SchemaFactory::required_field<FieldType::T201>(context, node, "wake"),
                    SchemaFactory::required_field<FieldType::T112>(context, node, "loads"),
                });
            break;
        case 77:
            owner.records_AnalysisTransient.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AnalysisTransient{
                    SchemaFactory::required_field<FieldType::T202>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T198>(context, node, "physics"),
                    SchemaFactory::required_field<FieldType::T199>(context, node, "branch"),
                    SchemaFactory::required_field<FieldType::T203>(context, node, "intervals"),
                    SchemaFactory::required_field<FieldType::T201>(context, node, "wake"),
                    SchemaFactory::required_field<FieldType::T112>(context, node, "loads"),
                });
            break;
        case 78:
            owner.records_AnalysisFlutter.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AnalysisFlutter{
                    SchemaFactory::required_field<FieldType::T204>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T205>(context, node, "physics"),
                    SchemaFactory::required_field<FieldType::T206>(context, node, "branches"),
                    SchemaFactory::required_field<FieldType::T207>(context, node, "intervals"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "step"),
                    SchemaFactory::required_field<FieldType::T201>(context, node, "wake"),
                    SchemaFactory::required_field<FieldType::T208>(context, node, "search"),
                    SchemaFactory::required_field<FieldType::T112>(context, node, "loads"),
                });
            break;
        case 79:
            owner.records_AnalysisSpeedNone.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AnalysisSpeedNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 80:
            owner.records_AnalysisSpeedValue.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AnalysisSpeedValue{
                    SchemaFactory::required_field<FieldType::T209>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T187>(context, node, "value"),
                });
            break;
        case 81:
            owner.records_Branch.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Branch{
                    SchemaFactory::required_field<FieldType::T210>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T211>(context, node, "seed"),
                    SchemaFactory::required_field<FieldType::T212>(context, node, "continuation"),
                });
            break;
        case 82:
            owner.records_BranchSeedReference.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BranchSeedReference{
                    SchemaFactory::required_field<FieldType::T213>(context, node, "kind"),
                });
            break;
        case 83:
            owner.records_BranchSeedPrepared.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BranchSeedPrepared{
                    SchemaFactory::required_field<FieldType::T214>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T215>(context, node, "history"),
                });
            break;
        case 84:
            owner.records_ContinuationFixed.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ContinuationFixed{
                    SchemaFactory::required_field<FieldType::T216>(context, node, "kind"),
                });
            break;
        case 85:
            owner.records_ContinuationSweep.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ContinuationSweep{
                    SchemaFactory::required_field<FieldType::T217>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T187>(context, node,
                                                                   "initial_increment"),
                    SchemaFactory::required_field<FieldType::T187>(context, node,
                                                                   "minimum_increment"),
                    SchemaFactory::required_field<FieldType::T187>(context, node,
                                                                   "maximum_increment"),
                    SchemaFactory::required_field<FieldType::T218>(context, node, "arclength"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "maximum_points"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_retries"),
                });
            break;
        case 86:
            owner.records_ArcLength.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::ArcLength{
                    SchemaFactory::required_field<FieldType::T220>(context, node, "state_scales"),
                    SchemaFactory::required_field<FieldType::T187>(context, node, "speed_scale"),
                    SchemaFactory::required_field<FieldType::T221>(context, node, "orientation"),
                    SchemaFactory::required_field<FieldType::T135>(context, node,
                                                                   "enable_at_limit"),
                });
            break;
        case 87:
            owner.records_TimeInterval.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::TimeInterval{
                    SchemaFactory::required_field<FieldType::T44>(context, node, "start"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "end"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "step"),
                    SchemaFactory::required_field<FieldType::T222>(context, node, "restart"),
                });
            break;
        case 88:
            owner.records_SpeedInterval.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::SpeedInterval{
                    SchemaFactory::required_field<FieldType::T187>(context, node, "lower"),
                    SchemaFactory::required_field<FieldType::T187>(context, node, "upper"),
                });
            break;
        case 89:
            owner.records_Search.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Search{
                    SchemaFactory::required_field<FieldType::T223>(context, node, "queue"),
                    SchemaFactory::required_field<FieldType::T187>(context, node, "speed_absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "speed_relative"),
                    SchemaFactory::required_field<FieldType::T224>(context, node, "growth_budget"),
                    SchemaFactory::required_field<FieldType::T225>(context, node,
                                                                   "mode_association"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "maximum_points"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_retries"),
                    SchemaFactory::required_field<FieldType::T226>(context, node, "coverage"),
                });
            break;
        case 90:
            owner.records_Association.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Association{
                    SchemaFactory::required_field<FieldType::T89>(context, node,
                                                                  "physical_observations"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "mac_minimum"),
                    SchemaFactory::required_field<FieldType::T186>(context, node,
                                                                   "subspace_angle_maximum"),
                    SchemaFactory::required_field<FieldType::T224>(context, node,
                                                                   "eigenvalue_absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node,
                                                                  "eigenvalue_relative"),
                });
            break;
        case 91:
            owner.records_Coverage.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Coverage{
                    SchemaFactory::required_field<FieldType::T227>(context, node, "branches"),
                    SchemaFactory::required_field<FieldType::T89>(context, node,
                                                                  "mode_observations"),
                    SchemaFactory::required_field<FieldType::T135>(context, node,
                                                                   "require_physical"),
                });
            break;
        case 92:
            owner.records_InitialEquilibrium.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::InitialEquilibrium{
                    SchemaFactory::required_field<FieldType::T228>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T229>(context, node, "branch"),
                    SchemaFactory::required_field<FieldType::T230>(context, node, "perturbation"),
                });
            break;
        case 93:
            owner.records_InitialCold.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::InitialCold{
                    SchemaFactory::required_field<FieldType::T231>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T232>(context, node, "state"),
                    SchemaFactory::required_field<FieldType::T233>(context, node, "circulation"),
                    SchemaFactory::required_field<FieldType::T234>(context, node, "wake"),
                });
            break;
        case 94:
            owner.records_InitialPrepared.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::InitialPrepared{
                    SchemaFactory::required_field<FieldType::T235>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T236>(context, node, "file"),
                });
            break;
        case 95:
            owner.records_PerturbationNone.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::PerturbationNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 96:
            owner.records_InitialState.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::InitialState{
                    SchemaFactory::required_field<FieldType::T237>(context, node, "shells"),
                    SchemaFactory::required_field<FieldType::T238>(context, node, "beams"),
                    SchemaFactory::required_field<FieldType::T239>(context, node, "rigids"),
                });
            break;
        case 97:
            owner.records_CoefficientState.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::CoefficientState{
                    SchemaFactory::required_field<FieldType::T177>(context, node, "shell"),
                    SchemaFactory::required_field<FieldType::T240>(context, node, "reference"),
                    SchemaFactory::required_field<FieldType::T241>(context, node, "representation"),
                    SchemaFactory::required_field<FieldType::T242>(context, node, "position"),
                    SchemaFactory::required_field<FieldType::T242>(context, node, "velocity"),
                });
            break;
        case 98:
            owner.records_PlanesBinary64.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::PlanesBinary64{
                    SchemaFactory::required_field<FieldType::T243>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "high"),
                });
            break;
        case 99:
            owner.records_PlanesTwoWord.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::PlanesTwoWord{
                    SchemaFactory::required_field<FieldType::T244>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "high"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "low"),
                });
            break;
        case 100:
            owner.records_BeamState.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::BeamState{
                    SchemaFactory::required_field<FieldType::T245>(context, node, "beam"),
                    SchemaFactory::required_field<FieldType::T246>(context, node, "nodes"),
                });
            break;
        case 101:
            owner.records_RigidState.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::RigidState{
                    SchemaFactory::required_field<FieldType::T247>(context, node, "rigid"),
                    SchemaFactory::required_field<FieldType::T248>(context, node, "state"),
                });
            break;
        case 102:
            owner.records_FrameStatePayload.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::FrameStatePayload{
                    SchemaFactory::required_field<FieldType::T26>(context, node, "position"),
                    SchemaFactory::required_field<FieldType::T27>(context, node, "rotation"),
                    SchemaFactory::required_field<FieldType::T249>(context, node,
                                                                   "linear_velocity"),
                    SchemaFactory::required_field<FieldType::T250>(context, node,
                                                                   "angular_velocity"),
                });
            break;
        case 103:
            owner.records_InitialWake.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::InitialWake{
                    SchemaFactory::required_field<FieldType::T251>(context, node, "surfaces"),
                });
            break;
        case 104:
            owner.records_SurfaceWake.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::SurfaceWake{
                    SchemaFactory::required_field<FieldType::T178>(context, node, "surface"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "rows"),
                    SchemaFactory::required_field<FieldType::T252>(context, node, "points"),
                    SchemaFactory::required_field<FieldType::T233>(context, node, "coefficients"),
                    SchemaFactory::required_field<FieldType::T253>(context, node, "row_ids"),
                    SchemaFactory::required_field<FieldType::T254>(context, node, "point_ids"),
                    SchemaFactory::required_field<FieldType::T255>(context, node, "band_ids"),
                    SchemaFactory::required_field<FieldType::T48>(context, node, "birth_times"),
                });
            break;
        case 105:
            owner.records_Restart.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Restart{
                    SchemaFactory::required_field<FieldType::T24>(context, node, "directory"),
                    SchemaFactory::required_field<FieldType::T240>(context, node,
                                                                   "manifest_sha256"),
                    SchemaFactory::required_field<FieldType::T256>(context, node, "compatibility"),
                    SchemaFactory::required_field<FieldType::T257>(context, node, "budget"),
                });
            break;
        case 106:
            owner.records_ObservationPoint.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ObservationPoint{
                    SchemaFactory::required_field<FieldType::T258>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T259>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "point"),
                    SchemaFactory::required_field<FieldType::T260>(context, node, "quantity"),
                    SchemaFactory::required_field<FieldType::T261>(context, node, "basis"),
                    SchemaFactory::required_field<FieldType::T262>(context, node, "comparison"),
                    SchemaFactory::required_field<FieldType::T263>(context, node, "budget"),
                });
            break;
        case 107:
            owner.records_ObservationStress.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ObservationStress{
                    SchemaFactory::required_field<FieldType::T264>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T259>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T142>(context, node, "point"),
                    SchemaFactory::required_field<FieldType::T265>(context, node, "layer"),
                    SchemaFactory::required_field<FieldType::T266>(context, node, "interface"),
                    SchemaFactory::required_field<FieldType::T267>(context, node, "measure"),
                    SchemaFactory::required_field<FieldType::T268>(context, node, "contribution"),
                    SchemaFactory::required_field<FieldType::T261>(context, node, "basis"),
                    SchemaFactory::required_field<FieldType::T262>(context, node, "comparison"),
                    SchemaFactory::required_field<FieldType::T263>(context, node, "budget"),
                });
            break;
        case 108:
            owner.records_ObservationSection.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ObservationSection{
                    SchemaFactory::required_field<FieldType::T269>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T259>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "region"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T261>(context, node, "basis"),
                    SchemaFactory::required_field<FieldType::T270>(context, node, "quantity"),
                    SchemaFactory::required_field<FieldType::T262>(context, node, "comparison"),
                    SchemaFactory::required_field<FieldType::T263>(context, node, "budget"),
                });
            break;
        case 109:
            owner.records_ObservationReaction.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ObservationReaction{
                    SchemaFactory::required_field<FieldType::T271>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T259>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T272>(context, node, "support"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T261>(context, node, "basis"),
                    SchemaFactory::required_field<FieldType::T262>(context, node, "comparison"),
                    SchemaFactory::required_field<FieldType::T263>(context, node, "budget"),
                });
            break;
        case 110:
            owner.records_ObservationAero.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ObservationAero{
                    SchemaFactory::required_field<FieldType::T273>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T259>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T274>(context, node, "surfaces"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "origin"),
                    SchemaFactory::required_field<FieldType::T261>(context, node, "basis"),
                    SchemaFactory::required_field<FieldType::T275>(context, node, "quantity"),
                    SchemaFactory::required_field<FieldType::T262>(context, node, "comparison"),
                    SchemaFactory::required_field<FieldType::T263>(context, node, "budget"),
                });
            break;
        case 111:
            owner.records_ObservationMode.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ObservationMode{
                    SchemaFactory::required_field<FieldType::T276>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T259>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T277>(context, node, "targets"),
                    SchemaFactory::required_field<FieldType::T278>(context, node, "points"),
                    SchemaFactory::required_field<FieldType::T261>(context, node, "basis"),
                    SchemaFactory::required_field<FieldType::T279>(context, node, "quantity"),
                    SchemaFactory::required_field<FieldType::T262>(context, node, "comparison"),
                    SchemaFactory::required_field<FieldType::T263>(context, node, "budget"),
                });
            break;
        case 112:
            owner.records_BasisWorldMaterial.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BasisWorldMaterial{
                    SchemaFactory::required_field<FieldType::T280>(context, node, "kind"),
                });
            break;
        case 113:
            owner.records_BasisExplicit.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BasisExplicit{
                    SchemaFactory::required_field<FieldType::T57>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T27>(context, node, "rotation"),
                });
            break;
        case 114:
            owner.records_Comparison.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Comparison{
                    SchemaFactory::required_field<FieldType::T46>(context, node, "window"),
                    SchemaFactory::required_field<FieldType::T281>(context, node, "statistic"),
                    SchemaFactory::required_field<FieldType::T282>(context, node, "side"),
                    SchemaFactory::required_field<FieldType::T283>(context, node, "reference"),
                });
            break;
        case 115:
            owner.records_ComparisonReferenceNone
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ComparisonReferenceNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 116:
            owner.records_ComparisonReferenceFile
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ComparisonReferenceFile{
                    SchemaFactory::required_field<FieldType::T284>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T285>(context, node, "source"),
                    SchemaFactory::required_field<FieldType::T25>(context, node, "columns"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "uncertainty"),
                });
            break;
        case 117:
            owner.records_ObservationBudget.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ObservationBudget{
                    SchemaFactory::required_field<FieldType::T171>(context, node, "absolute"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "relative"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "scale"),
                    SchemaFactory::required_field<FieldType::T286>(context, node, "groups"),
                });
            break;
        case 118:
            owner.records_ErrorGroupNotApplicable
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ErrorGroupNotApplicable{
                    SchemaFactory::required_field<FieldType::T287>(context, node, "status"),
                    SchemaFactory::required_field<FieldType::T288>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "reason"),
                });
            break;
        case 119:
            owner.records_ErrorGroupCoveredBy.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ErrorGroupCoveredBy{
                    SchemaFactory::required_field<FieldType::T289>(context, node, "status"),
                    SchemaFactory::required_field<FieldType::T288>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T288>(context, node, "representative"),
                });
            break;
        case 120:
            owner.records_ErrorGroupRequired.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ErrorGroupRequired{
                    SchemaFactory::required_field<FieldType::T290>(context, node, "status"),
                    SchemaFactory::required_field<FieldType::T288>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T291>(context, node, "covers"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "absolute"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "relative"),
                });
            break;
        case 121:
            owner.records_Numerics.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Numerics{
                    SchemaFactory::required_field<FieldType::T292>(context, node, "patch"),
                    SchemaFactory::required_field<FieldType::T293>(context, node, "quadrature"),
                    SchemaFactory::required_field<FieldType::T294>(context, node, "linear"),
                    SchemaFactory::required_field<FieldType::T295>(context, node, "newton"),
                    SchemaFactory::required_field<FieldType::T296>(context, node, "precision"),
                    SchemaFactory::required_field<FieldType::T297>(context, node, "eigen"),
                    SchemaFactory::required_field<FieldType::T298>(context, node, "wake_rows"),
                    SchemaFactory::required_field<FieldType::T220>(context, node, "scales"),
                    SchemaFactory::required_field<FieldType::T286>(context, node, "error_groups"),
                });
            break;
        case 122:
            owner.records_WakeCapacityNone.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::WakeCapacityNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 123:
            owner.records_WakeCapacityLimited.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::WakeCapacityLimited{
                    SchemaFactory::required_field<FieldType::T299>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "maximum_rows"),
                });
            break;
        case 124:
            owner.records_PatchPolicyNone.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::PatchPolicyNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 125:
            owner.records_PatchPolicyActive.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::PatchPolicyActive{
                    SchemaFactory::required_field<FieldType::T189>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_degree"),
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "maximum_partition"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "geometry_error"),
                    SchemaFactory::required_field<FieldType::T61>(context, node, "trace_error"),
                    SchemaFactory::required_field<FieldType::T300>(context, node, "rank"),
                });
            break;
        case 126:
            owner.records_Quadrature.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Quadrature{
                    SchemaFactory::required_field<FieldType::T301>(context, node, "structure"),
                    SchemaFactory::required_field<FieldType::T302>(context, node, "aero"),
                });
            break;
        case 127:
            owner.records_StructuralQuadratureNone
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::StructuralQuadratureNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 128:
            owner.records_StructuralQuadratureActive
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::StructuralQuadratureActive{
                    SchemaFactory::required_field<FieldType::T189>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "surface_degree"),
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "thickness_points"),
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "independent_surface_degree"),
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "independent_thickness_points"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "relative"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_subdivisions"),
                });
            break;
        case 129:
            owner.records_AeroQuadratureNone.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AeroQuadratureNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 130:
            owner.records_AeroQuadratureActive.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::AeroQuadratureActive{
                    SchemaFactory::required_field<FieldType::T303>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "relative"),
                });
            break;
        case 131:
            owner.records_ScaleUniformComponents
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ScaleUniformComponents{
                    SchemaFactory::required_field<FieldType::T304>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T305>(context, node, "field"),
                    SchemaFactory::required_field<FieldType::T277>(context, node, "targets"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "values"),
                });
            break;
        case 132:
            owner.records_ScalePhysicalLayoutPolicy
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ScalePhysicalLayoutPolicy{
                    SchemaFactory::required_field<FieldType::T306>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T307>(context, node, "field"),
                    SchemaFactory::required_field<FieldType::T308>(context, node, "constraint"),
                });
            break;
        case 133:
            owner.records_ConstraintScalePolicyNone
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ConstraintScalePolicyNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 134:
            owner.records_ConstraintScalePolicyPhysical
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ConstraintScalePolicyPhysical{
                    SchemaFactory::required_field<FieldType::T309>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T310>(context, node, "position"),
                    SchemaFactory::required_field<FieldType::T310>(context, node, "orientation"),
                    SchemaFactory::required_field<FieldType::T311>(context, node,
                                                                   "multiplier_work_scale"),
                    SchemaFactory::required_field<FieldType::T312>(context, node, "reaction_force"),
                    SchemaFactory::required_field<FieldType::T313>(context, node,
                                                                   "reaction_moment"),
                });
            break;
        case 135:
            owner.records_ComponentScale3None.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ComponentScale3None{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 136:
            owner.records_ComponentScale3Values.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ComponentScale3Values{
                    SchemaFactory::required_field<FieldType::T314>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "values"),
                });
            break;
        case 137:
            owner.records_ScaleExplicitLayout.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ScaleExplicitLayout{
                    SchemaFactory::required_field<FieldType::T315>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T305>(context, node, "field"),
                    SchemaFactory::required_field<FieldType::T277>(context, node, "targets"),
                    SchemaFactory::required_field<FieldType::T240>(context, node, "layout"),
                    SchemaFactory::required_field<FieldType::T316>(context, node, "dimensions"),
                    SchemaFactory::required_field<FieldType::T171>(context, node, "values"),
                });
            break;
        case 138:
            owner.records_LinearPolicyDirect.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::LinearPolicyDirect{
                    SchemaFactory::required_field<FieldType::T317>(context, node, "backend"),
                    SchemaFactory::required_field<FieldType::T318>(context, node, "accuracy"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_refinements"),
                    SchemaFactory::required_field<FieldType::T319>(context, node, "backend_memory"),
                });
            break;
        case 139:
            owner.records_LinearPolicyGmres.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::LinearPolicyGmres{
                    SchemaFactory::required_field<FieldType::T320>(context, node, "backend"),
                    SchemaFactory::required_field<FieldType::T318>(context, node, "accuracy"),
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "maximum_iterations"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "restart"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_applications"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_adjoint"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_preconditioner"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_true_residual_checks"),
                    SchemaFactory::required_field<FieldType::T319>(context, node, "backend_memory"),
                });
            break;
        case 140:
            owner.records_LinearAccuracy.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::LinearAccuracy{
                    SchemaFactory::required_field<FieldType::T88>(context, node, "absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "relative"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "backward"),
                    SchemaFactory::required_field<FieldType::T321>(context, node, "norm"),
                });
            break;
        case 141:
            owner.records_NewtonPolicy.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::NewtonPolicy{
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "maximum_iterations"),
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "maximum_line_searches"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "relative"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "correction"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "minimum_step"),
                });
            break;
        case 142:
            owner.records_PrecisionPolicy.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::PrecisionPolicy{
                    SchemaFactory::required_field<FieldType::T322>(context, node, "input"),
                    SchemaFactory::required_field<FieldType::T323>(context, node, "retained"),
                    SchemaFactory::required_field<FieldType::T324>(context, node, "returned"),
                    SchemaFactory::required_field<FieldType::T88>(context, node,
                                                                  "projection_absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node,
                                                                  "projection_relative"),
                });
            break;
        case 143:
            owner.records_BackendMemory.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BackendMemory{
                    SchemaFactory::required_field<FieldType::T325>(context, node, "admission"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_bytes"),
                });
            break;
        case 144:
            owner.records_EigenPolicyNone.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::EigenPolicyNone{
                    SchemaFactory::required_field<FieldType::T179>(context, node, "kind"),
                });
            break;
        case 145:
            owner.records_EigenPolicyFull.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::EigenPolicyFull{
                    SchemaFactory::required_field<FieldType::T326>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T327>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T328>(context, node,
                                                                   "eigenvalue_scale"),
                    SchemaFactory::required_field<FieldType::T329>(context, node, "packet"),
                    SchemaFactory::required_field<FieldType::T330>(context, node, "requirements"),
                    SchemaFactory::required_field<FieldType::T331>(context, node, "limits"),
                });
            break;
        case 146:
            owner.records_EigenPolicyGeneralized
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::EigenPolicyGeneralized{
                    SchemaFactory::required_field<FieldType::T332>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T333>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T334>(context, node,
                                                                   "eigenvalue_scale"),
                    SchemaFactory::required_field<FieldType::T330>(context, node, "requirements"),
                    SchemaFactory::required_field<FieldType::T331>(context, node, "limits"),
                });
            break;
        case 147:
            owner.records_EigenPolicyArnoldi.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::EigenPolicyArnoldi{
                    SchemaFactory::required_field<FieldType::T335>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T327>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T328>(context, node,
                                                                   "eigenvalue_scale"),
                    SchemaFactory::required_field<FieldType::T336>(context, node, "packet"),
                    SchemaFactory::required_field<FieldType::T337>(context, node, "selector"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "nev"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "ncv"),
                    SchemaFactory::required_field<FieldType::T338>(context, node, "seed"),
                    SchemaFactory::required_field<FieldType::T88>(context, node,
                                                                  "backend_tolerance"),
                    SchemaFactory::required_field<FieldType::T330>(context, node, "requirements"),
                    SchemaFactory::required_field<FieldType::T331>(context, node, "limits"),
                });
            break;
        case 148:
            owner.records_Complex.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Complex{
                    SchemaFactory::required_field<FieldType::T328>(context, node, "real"),
                    SchemaFactory::required_field<FieldType::T328>(context, node, "imaginary"),
                });
            break;
        case 149:
            owner.records_ArnoldiSeedUniform.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ArnoldiSeedUniform{
                    SchemaFactory::required_field<FieldType::T339>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T340>(context, node, "value"),
                });
            break;
        case 150:
            owner.records_ArnoldiSeedExplicit.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ArnoldiSeedExplicit{
                    SchemaFactory::required_field<FieldType::T315>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T240>(context, node, "layout"),
                    SchemaFactory::required_field<FieldType::T341>(context, node, "values"),
                });
            break;
        case 151:
            owner.records_EigenSelectorLargest.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::EigenSelectorLargest{
                    SchemaFactory::required_field<FieldType::T342>(context, node, "kind"),
                });
            break;
        case 152:
            owner.records_EigenSelectorShift.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::EigenSelectorShift{
                    SchemaFactory::required_field<FieldType::T343>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T340>(context, node, "shift"),
                });
            break;
        case 153:
            owner.records_EigenRequirements.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::EigenRequirements{
                    SchemaFactory::required_field<FieldType::T88>(context, node, "residual"),
                    SchemaFactory::required_field<FieldType::T344>(context, node,
                                                                   "residual_observation"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "normalization"),
                    SchemaFactory::required_field<FieldType::T344>(context, node,
                                                                   "normalization_observation"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "orthogonality"),
                    SchemaFactory::required_field<FieldType::T344>(context, node,
                                                                   "orthogonality_observation"),
                    SchemaFactory::required_field<FieldType::T345>(context, node, "backward"),
                    SchemaFactory::required_field<FieldType::T88>(context, node,
                                                                  "association_absolute"),
                    SchemaFactory::required_field<FieldType::T88>(context, node,
                                                                  "association_relative"),
                    SchemaFactory::required_field<FieldType::T135>(context, node, "left"),
                    SchemaFactory::required_field<FieldType::T346>(context, node, "clusters"),
                    SchemaFactory::required_field<FieldType::T135>(context, node,
                                                                   "physical_coverage"),
                });
            break;
        case 154:
            owner.records_BackwardNone.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::BackwardNone{
                    SchemaFactory::required_field<FieldType::T347>(context, node, "kind"),
                });
            break;
        case 155:
            owner.records_BackwardRequired.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::BackwardRequired{
                    SchemaFactory::required_field<FieldType::T290>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T88>(context, node, "budget"),
                    SchemaFactory::required_field<FieldType::T344>(context, node, "observation"),
                });
            break;
        case 156:
            owner.records_EigenCluster.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::EigenCluster{
                    SchemaFactory::required_field<FieldType::T143>(context, node, "members"),
                    SchemaFactory::required_field<FieldType::T219>(context, node,
                                                                   "expected_dimension"),
                });
            break;
        case 157:
            owner.records_EigenLimits.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::EigenLimits{
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_order"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_results"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_ncv"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_backend_iterations"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_attempts"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_forward"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_adjoint"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_inverse"),
                    SchemaFactory::required_field<FieldType::T319>(context, node, "backend_memory"),
                });
            break;
        case 158:
            owner.records_Resources.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Resources{
                    SchemaFactory::required_field<FieldType::T348>(context, node, "profile"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "threads"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_memory_bytes"),
                    SchemaFactory::required_field<FieldType::T44>(context, node,
                                                                  "maximum_wall_seconds"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_steps"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_total_iterations"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "checkpoint_temporary_bytes"),
                    SchemaFactory::required_field<FieldType::T349>(context, node, "input"),
                });
            break;
        case 159:
            owner.records_InputQuotas.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::InputQuotas{
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_source_bytes"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_depth"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_nodes"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_members"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_decoded_string_bytes"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_decimal_token_bytes"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_total_text_bytes"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_typed_bytes"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_entities"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_references"),
                    SchemaFactory::required_field<FieldType::T67>(context, node, "maximum_scalars"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_polynomial_terms"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_table_samples"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_expression_nodes"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_expression_edges"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_expression_depth"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_expression_power"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_expression_points"),
                    SchemaFactory::required_field<FieldType::T67>(context, node,
                                                                  "maximum_diagnostic_path_bytes"),
                });
            break;
        case 160:
            owner.records_Output.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Output{
                    SchemaFactory::required_field<FieldType::T24>(context, node, "directory"),
                    SchemaFactory::required_field<FieldType::T350>(context, node, "observations"),
                    SchemaFactory::required_field<FieldType::T351>(context, node, "visualization"),
                    SchemaFactory::required_field<FieldType::T352>(context, node, "checkpoint"),
                    SchemaFactory::required_field<FieldType::T353>(context, node, "overwrite"),
                    SchemaFactory::required_field<FieldType::T354>(context, node, "resume"),
                });
            break;
        case 161:
            owner.records_OutputCadence.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::OutputCadence{
                    SchemaFactory::required_field<FieldType::T219>(context, node, "every_steps"),
                    SchemaFactory::required_field<FieldType::T89>(context, node, "observations"),
                });
            break;
        case 162:
            owner.records_Visualization.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::Visualization{
                    SchemaFactory::required_field<FieldType::T135>(context, node, "enabled"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "every_steps"),
                    SchemaFactory::required_field<FieldType::T355>(context, node, "layers"),
                    SchemaFactory::required_field<FieldType::T356>(context, node, "interfaces"),
                    SchemaFactory::required_field<FieldType::T135>(context, node,
                                                                   "reference_geometry"),
                    SchemaFactory::required_field<FieldType::T135>(context, node,
                                                                   "current_geometry"),
                });
            break;
        case 163:
            owner.records_CheckpointOutputOff.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::CheckpointOutputOff{
                    SchemaFactory::required_field<FieldType::T357>(context, node, "kind"),
                });
            break;
        case 164:
            owner.records_CheckpointOutputOn.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::CheckpointOutputOn{
                    SchemaFactory::required_field<FieldType::T358>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "every_steps"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "directory"),
                    SchemaFactory::required_field<FieldType::T219>(context, node, "retain"),
                });
            break;
        case 165:
            owner.records_ReferenceFieldExpression
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ReferenceFieldExpression{
                    SchemaFactory::required_field<FieldType::T359>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T360>(context, node, "expression"),
                    SchemaFactory::required_field<FieldType::T39>(context, node, "error"),
                });
            break;
        case 166:
            owner.records_TimeSeriesExpression.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::TimeSeriesExpression{
                    SchemaFactory::required_field<FieldType::T359>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T360>(context, node, "expression"),
                });
            break;
        case 167:
            owner.records_LoadCurrentExpression.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::LoadCurrentExpression{
                    SchemaFactory::required_field<FieldType::T361>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T157>(context, node, "id"),
                    SchemaFactory::required_field<FieldType::T133>(context, node, "region"),
                    SchemaFactory::required_field<FieldType::T362>(context, node, "expression"),
                });
            break;
        case 168:
            owner.records_Expression.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Expression{
                    SchemaFactory::required_field<FieldType::T363>(context, node, "domain"),
                    SchemaFactory::required_field<FieldType::T364>(context, node, "nodes"),
                    SchemaFactory::required_field<FieldType::T365>(context, node, "outputs"),
                });
            break;
        case 169:
            owner.records_ExpressionDomainReference
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExpressionDomainReference{
                    SchemaFactory::required_field<FieldType::T213>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T41>(context, node, "coordinates"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "lower"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "upper"),
                });
            break;
        case 170:
            owner.records_ExpressionDomainTime.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExpressionDomainTime{
                    SchemaFactory::required_field<FieldType::T366>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "lower"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "upper"),
                    SchemaFactory::required_field<FieldType::T48>(context, node, "events"),
                });
            break;
        case 171:
            owner.records_ExpressionDomainCurrentPositionTime
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExpressionDomainCurrentPositionTime{
                    SchemaFactory::required_field<FieldType::T367>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "position_lower"),
                    SchemaFactory::required_field<FieldType::T26>(context, node, "position_upper"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "time_lower"),
                    SchemaFactory::required_field<FieldType::T44>(context, node, "time_upper"),
                    SchemaFactory::required_field<FieldType::T48>(context, node, "events"),
                });
            break;
        case 172:
            owner.records_ExprNodeConstant.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExprNodeConstant{
                    SchemaFactory::required_field<FieldType::T38>(context, node, "op"),
                    SchemaFactory::required_field<FieldType::T328>(context, node, "value"),
                    SchemaFactory::required_field<FieldType::T368>(context, node, "dimensions"),
                });
            break;
        case 173:
            owner.records_ExprNodeVariable.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExprNodeVariable{
                    SchemaFactory::required_field<FieldType::T369>(context, node, "op"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "axis"),
                });
            break;
        case 174:
            owner.records_ExprNodeBinary.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExprNodeBinary{
                    SchemaFactory::required_field<FieldType::T370>(context, node, "op"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "first"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "second"),
                });
            break;
        case 175:
            owner.records_ExprNodeNegate.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExprNodeNegate{
                    SchemaFactory::required_field<FieldType::T371>(context, node, "op"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "operand"),
                });
            break;
        case 176:
            owner.records_ExprNodePower.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExprNodePower{
                    SchemaFactory::required_field<FieldType::T372>(context, node, "op"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "operand"),
                    SchemaFactory::required_field<FieldType::T373>(context, node, "exponent"),
                });
            break;
        case 177:
            owner.records_ExprNodeAnalytic.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ExprNodeAnalytic{
                    SchemaFactory::required_field<FieldType::T374>(context, node, "op"),
                    SchemaFactory::required_field<FieldType::T35>(context, node, "operand"),
                });
            break;
        case 178:
            owner.records_PerturbationPhysical.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::PerturbationPhysical{
                    SchemaFactory::required_field<FieldType::T309>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T375>(context, node, "shells"),
                    SchemaFactory::required_field<FieldType::T238>(context, node, "beams"),
                    SchemaFactory::required_field<FieldType::T239>(context, node, "rigids"),
                    SchemaFactory::required_field<FieldType::T376>(context, node, "displacements"),
                    SchemaFactory::required_field<FieldType::T376>(context, node, "velocities"),
                });
            break;
        case 179:
            owner.records_ShellPhysicalState.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::ShellPhysicalState{
                    SchemaFactory::required_field<FieldType::T177>(context, node, "shell"),
                    SchemaFactory::required_field<FieldType::T377>(context, node, "displacement"),
                    SchemaFactory::required_field<FieldType::T377>(context, node, "velocity"),
                });
            break;
        case 180:
            owner.records_FullShellFieldGlobalPolynomial
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::FullShellFieldGlobalPolynomial{
                    SchemaFactory::required_field<FieldType::T378>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T379>(context, node, "r"),
                    SchemaFactory::required_field<FieldType::T379>(context, node, "p"),
                    SchemaFactory::required_field<FieldType::T379>(context, node, "c"),
                });
            break;
        case 181:
            owner.records_FullShellFieldHalfCells
                .values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::FullShellFieldHalfCells{
                    SchemaFactory::required_field<FieldType::T380>(context, node, "kind"),
                    SchemaFactory::required_field<FieldType::T381>(context, node, "cells"),
                });
            break;
        case 182:
            owner.records_CellField.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::CellField{
                    SchemaFactory::required_field<FieldType::T136>(context, node, "cell"),
                    SchemaFactory::required_field<FieldType::T379>(context, node, "r"),
                    SchemaFactory::required_field<FieldType::T379>(context, node, "p"),
                    SchemaFactory::required_field<FieldType::T379>(context, node, "c"),
                });
            break;
        case 183:
            owner.records_VectorPolynomial.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::VectorPolynomial{
                    SchemaFactory::required_field<FieldType::T382>(context, node, "terms"),
                });
            break;
        case 184:
            owner.records_VectorTerm.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::VectorTerm{
                    SchemaFactory::required_field<FieldType::T63>(context, node, "powers"),
                    SchemaFactory::required_field<FieldType::T383>(context, node, "value"),
                });
            break;
        case 185:
            owner.records_PointCheck.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::PointCheck{
                    SchemaFactory::required_field<FieldType::T384>(context, node, "point"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "value"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "tolerance"),
                });
            break;
        case 186:
            owner.records_InitialPointShell.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::InitialPointShell{
                    SchemaFactory::required_field<FieldType::T385>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "coordinates"),
                    SchemaFactory::required_field<FieldType::T32>(context, node, "side_s"),
                    SchemaFactory::required_field<FieldType::T32>(context, node, "side_t"),
                });
            break;
        case 187:
            owner.records_InitialPointBeam.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::InitialPointBeam{
                    SchemaFactory::required_field<FieldType::T386>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "coordinates"),
                });
            break;
        case 188:
            owner.records_InitialPointRigid.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::InitialPointRigid{
                    SchemaFactory::required_field<FieldType::T387>(context, node, "target"),
                    SchemaFactory::required_field<FieldType::T31>(context, node, "coordinates"),
                });
            break;
        case 189:
            owner.records_Rational.values()[static_cast<std::size_t>(typed.record_row)].emplace(
                schema::Rational{
                    SchemaFactory::required_field<FieldType::T24>(context, node, "numerator"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "denominator"),
                });
            break;
        case 190:
            owner.records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]
                .emplace(schema::FileReference{
                    SchemaFactory::required_field<FieldType::T388>(context, node, "source_id"),
                    SchemaFactory::required_field<FieldType::T24>(context, node, "path"),
                    SchemaFactory::required_field<FieldType::T240>(context, node, "sha256"),
                    SchemaFactory::required_field<FieldType::T446>(context, node, "role"),
                });
            break;
        default:
            context.reject(SchemaReason::WrongKind);
        }
    }
}
RecordAccess<FieldType::T1>::Value RecordAccess<FieldType::T1>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Units.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T2>::Value RecordAccess<FieldType::T2>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Model.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T3>::Value RecordAccess<FieldType::T3>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 65:
        return ready_record(
            owner->records_FlowNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 66:
        return ready_record(
            owner->records_FlowKnown.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_FlowNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T4>::Value RecordAccess<FieldType::T4>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 76:
        return ready_record(
            owner->records_AnalysisStatic.values()[static_cast<std::size_t>(typed.record_row)]);
    case 77:
        return ready_record(
            owner->records_AnalysisTransient.values()[static_cast<std::size_t>(typed.record_row)]);
    case 78:
        return ready_record(
            owner->records_AnalysisFlutter.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_AnalysisStatic.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T6>::Value RecordAccess<FieldType::T6>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Numerics.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T7>::Value RecordAccess<FieldType::T7>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Resources.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T8>::Value RecordAccess<FieldType::T8>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Output.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T9>::Value RecordAccess<FieldType::T9>::read(const SchemaData* owner,
                                                                     std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 92:
        return ready_record(
            owner->records_InitialEquilibrium.values()[static_cast<std::size_t>(typed.record_row)]);
    case 93:
        return ready_record(
            owner->records_InitialCold.values()[static_cast<std::size_t>(typed.record_row)]);
    case 94:
        return ready_record(
            owner->records_InitialPrepared.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_InitialEquilibrium.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T10>::Value
RecordAccess<FieldType::T10>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Restart.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T11>::Value
RecordAccess<FieldType::T11>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Provenance.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T30>::Value
RecordAccess<FieldType::T30>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Target.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T49>::Value
RecordAccess<FieldType::T49>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 14:
        return ready_record(
            owner->records_OutsideRejectHold.values()[static_cast<std::size_t>(typed.record_row)]);
    case 15:
        return ready_record(
            owner->records_OutsidePeriodic.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_OutsideRejectHold.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T52>::Value
RecordAccess<FieldType::T52>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Frame.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T53>::Value
RecordAccess<FieldType::T53>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 11:
        return ready_record(
            owner->records_TimeSeriesConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    case 12:
        return ready_record(owner->records_TimeSeriesPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 13:
        return ready_record(
            owner->records_TimeSeriesPchip.values()[static_cast<std::size_t>(typed.record_row)]);
    case 166:
        return ready_record(owner->records_TimeSeriesExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_TimeSeriesConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T55>::Value
RecordAccess<FieldType::T55>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T56>::Value
RecordAccess<FieldType::T56>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T58>::Value
RecordAccess<FieldType::T58>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T62>::Value
RecordAccess<FieldType::T62>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Rational.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T66>::Value
RecordAccess<FieldType::T66>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Thickness.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T77>::Value
RecordAccess<FieldType::T77>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T78>::Value
RecordAccess<FieldType::T78>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T82>::Value
RecordAccess<FieldType::T82>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T83>::Value
RecordAccess<FieldType::T83>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T84>::Value
RecordAccess<FieldType::T84>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 27:
        return ready_record(owner->records_ReferenceGeometrySupplied
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 28:
        return ready_record(owner->records_ReferenceGeometryFit
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceGeometrySupplied
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T86>::Value
RecordAccess<FieldType::T86>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T94>::Value
RecordAccess<FieldType::T94>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 32:
        return ready_record(owner->records_ApplicabilityUnevaluated
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 33:
        return ready_record(owner->records_ApplicabilityDeclared
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ApplicabilityUnevaluated
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T107>::Value
RecordAccess<FieldType::T107>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 17:
        return ready_record(
            owner->records_AxisProjected.values()[static_cast<std::size_t>(typed.record_row)]);
    case 18:
        return ready_record(
            owner->records_AxisExplicit.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_AxisProjected.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T108>::Value
RecordAccess<FieldType::T108>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T115>::Value
RecordAccess<FieldType::T115>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 37:
        return ready_record(owner->records_BeamCoordinateArclength
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 38:
        return ready_record(owner->records_BeamCoordinateScaled
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_BeamCoordinateArclength
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T116>::Value
RecordAccess<FieldType::T116>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 39:
        return ready_record(
            owner->records_SectionMatrix.values()[static_cast<std::size_t>(typed.record_row)]);
    case 40:
        return ready_record(
            owner->records_SectionCircle.values()[static_cast<std::size_t>(typed.record_row)]);
    case 41:
        return ready_record(
            owner->records_SectionRectangle.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_SectionMatrix.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T133>::Value
RecordAccess<FieldType::T133>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Region.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T142>::Value
RecordAccess<FieldType::T142>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_MaterialPoint.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T144>::Value
RecordAccess<FieldType::T144>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Motion.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T150>::Value
RecordAccess<FieldType::T150>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 52:
        return ready_record(
            owner->records_LedgerAdd.values()[static_cast<std::size_t>(typed.record_row)]);
    case 53:
        return ready_record(
            owner->records_LedgerReplace.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_LedgerAdd.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T163>::Value
RecordAccess<FieldType::T163>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T164>::Value
RecordAccess<FieldType::T164>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 11:
        return ready_record(
            owner->records_TimeSeriesConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    case 12:
        return ready_record(owner->records_TimeSeriesPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 13:
        return ready_record(
            owner->records_TimeSeriesPchip.values()[static_cast<std::size_t>(typed.record_row)]);
    case 166:
        return ready_record(owner->records_TimeSeriesExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_TimeSeriesConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T166>::Value
RecordAccess<FieldType::T166>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T168>::Value
RecordAccess<FieldType::T168>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T170>::Value
RecordAccess<FieldType::T170>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 61:
        return ready_record(owner->records_AeroAttachmentFixed
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 62:
        return ready_record(owner->records_AeroAttachmentPrescribed
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 63:
        return ready_record(owner->records_AeroAttachmentStructural
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_AeroAttachmentFixed
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T181>::Value
RecordAccess<FieldType::T181>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 67:
        return ready_record(
            owner->records_GustNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 68:
        return ready_record(
            owner->records_GustSine.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_GustNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T182>::Value
RecordAccess<FieldType::T182>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FlowReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T190>::Value
RecordAccess<FieldType::T190>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 72:
        return ready_record(
            owner->records_ConvectionFree.values()[static_cast<std::size_t>(typed.record_row)]);
    case 73:
        return ready_record(owner->records_ConvectionPrescribed
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_ConvectionFree.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T191>::Value
RecordAccess<FieldType::T191>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 74:
        return ready_record(
            owner->records_RetentionAll.values()[static_cast<std::size_t>(typed.record_row)]);
    case 75:
        return ready_record(
            owner->records_RetentionWindow.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_RetentionAll.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T199>::Value
RecordAccess<FieldType::T199>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Branch.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T200>::Value
RecordAccess<FieldType::T200>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 79:
        return ready_record(
            owner->records_AnalysisSpeedNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 80:
        return ready_record(
            owner->records_AnalysisSpeedValue.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_AnalysisSpeedNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T201>::Value
RecordAccess<FieldType::T201>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 70:
        return ready_record(
            owner->records_WakeNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 71:
        return ready_record(
            owner->records_WakeActive.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_WakeNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T208>::Value
RecordAccess<FieldType::T208>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Search.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T211>::Value
RecordAccess<FieldType::T211>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 82:
        return ready_record(owner->records_BranchSeedReference
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 83:
        return ready_record(
            owner->records_BranchSeedPrepared.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_BranchSeedReference
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T212>::Value
RecordAccess<FieldType::T212>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 84:
        return ready_record(
            owner->records_ContinuationFixed.values()[static_cast<std::size_t>(typed.record_row)]);
    case 85:
        return ready_record(
            owner->records_ContinuationSweep.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_ContinuationFixed.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T215>::Value
RecordAccess<FieldType::T215>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T218>::Value
RecordAccess<FieldType::T218>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_ArcLength.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T225>::Value
RecordAccess<FieldType::T225>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Association.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T226>::Value
RecordAccess<FieldType::T226>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Coverage.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T230>::Value
RecordAccess<FieldType::T230>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 95:
        return ready_record(
            owner->records_PerturbationNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 178:
        return ready_record(owner->records_PerturbationPhysical
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_PerturbationNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T232>::Value
RecordAccess<FieldType::T232>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_InitialState.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T234>::Value
RecordAccess<FieldType::T234>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_InitialWake.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T236>::Value
RecordAccess<FieldType::T236>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T242>::Value
RecordAccess<FieldType::T242>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 98:
        return ready_record(
            owner->records_PlanesBinary64.values()[static_cast<std::size_t>(typed.record_row)]);
    case 99:
        return ready_record(
            owner->records_PlanesTwoWord.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_PlanesBinary64.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T248>::Value
RecordAccess<FieldType::T248>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FrameStatePayload.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T261>::Value
RecordAccess<FieldType::T261>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 112:
        return ready_record(
            owner->records_BasisWorldMaterial.values()[static_cast<std::size_t>(typed.record_row)]);
    case 113:
        return ready_record(
            owner->records_BasisExplicit.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_BasisWorldMaterial.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T262>::Value
RecordAccess<FieldType::T262>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Comparison.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T263>::Value
RecordAccess<FieldType::T263>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_ObservationBudget.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T283>::Value
RecordAccess<FieldType::T283>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 115:
        return ready_record(owner->records_ComparisonReferenceNone
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 116:
        return ready_record(owner->records_ComparisonReferenceFile
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ComparisonReferenceNone
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T285>::Value
RecordAccess<FieldType::T285>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T292>::Value
RecordAccess<FieldType::T292>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 124:
        return ready_record(
            owner->records_PatchPolicyNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 125:
        return ready_record(
            owner->records_PatchPolicyActive.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_PatchPolicyNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T293>::Value
RecordAccess<FieldType::T293>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Quadrature.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T294>::Value
RecordAccess<FieldType::T294>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 138:
        return ready_record(
            owner->records_LinearPolicyDirect.values()[static_cast<std::size_t>(typed.record_row)]);
    case 139:
        return ready_record(
            owner->records_LinearPolicyGmres.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_LinearPolicyDirect.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T295>::Value
RecordAccess<FieldType::T295>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_NewtonPolicy.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T296>::Value
RecordAccess<FieldType::T296>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_PrecisionPolicy.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T297>::Value
RecordAccess<FieldType::T297>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 144:
        return ready_record(
            owner->records_EigenPolicyNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 145:
        return ready_record(
            owner->records_EigenPolicyFull.values()[static_cast<std::size_t>(typed.record_row)]);
    case 146:
        return ready_record(owner->records_EigenPolicyGeneralized
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 147:
        return ready_record(
            owner->records_EigenPolicyArnoldi.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_EigenPolicyNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T298>::Value
RecordAccess<FieldType::T298>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 122:
        return ready_record(
            owner->records_WakeCapacityNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 123:
        return ready_record(owner->records_WakeCapacityLimited
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_WakeCapacityNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T301>::Value
RecordAccess<FieldType::T301>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 127:
        return ready_record(owner->records_StructuralQuadratureNone
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 128:
        return ready_record(owner->records_StructuralQuadratureActive
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_StructuralQuadratureNone
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T302>::Value
RecordAccess<FieldType::T302>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 129:
        return ready_record(
            owner->records_AeroQuadratureNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 130:
        return ready_record(owner->records_AeroQuadratureActive
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_AeroQuadratureNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T308>::Value
RecordAccess<FieldType::T308>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 133:
        return ready_record(owner->records_ConstraintScalePolicyNone
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 134:
        return ready_record(owner->records_ConstraintScalePolicyPhysical
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ConstraintScalePolicyNone
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T310>::Value
RecordAccess<FieldType::T310>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 135:
        return ready_record(owner->records_ComponentScale3None
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 136:
        return ready_record(owner->records_ComponentScale3Values
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ComponentScale3None
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T318>::Value
RecordAccess<FieldType::T318>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_LinearAccuracy.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T319>::Value
RecordAccess<FieldType::T319>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_BackendMemory.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T330>::Value
RecordAccess<FieldType::T330>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_EigenRequirements.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T331>::Value
RecordAccess<FieldType::T331>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_EigenLimits.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T337>::Value
RecordAccess<FieldType::T337>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 151:
        return ready_record(owner->records_EigenSelectorLargest
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 152:
        return ready_record(
            owner->records_EigenSelectorShift.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_EigenSelectorLargest
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T338>::Value
RecordAccess<FieldType::T338>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 149:
        return ready_record(
            owner->records_ArnoldiSeedUniform.values()[static_cast<std::size_t>(typed.record_row)]);
    case 150:
        return ready_record(owner->records_ArnoldiSeedExplicit
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_ArnoldiSeedUniform.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T340>::Value
RecordAccess<FieldType::T340>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Complex.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T345>::Value
RecordAccess<FieldType::T345>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 154:
        return ready_record(
            owner->records_BackwardNone.values()[static_cast<std::size_t>(typed.record_row)]);
    case 155:
        return ready_record(
            owner->records_BackwardRequired.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_BackwardNone.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T349>::Value
RecordAccess<FieldType::T349>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_InputQuotas.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T350>::Value
RecordAccess<FieldType::T350>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_OutputCadence.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T351>::Value
RecordAccess<FieldType::T351>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Visualization.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T352>::Value
RecordAccess<FieldType::T352>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 163:
        return ready_record(owner->records_CheckpointOutputOff
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 164:
        return ready_record(
            owner->records_CheckpointOutputOn.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_CheckpointOutputOff
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T360>::Value
RecordAccess<FieldType::T360>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Expression.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T362>::Value
RecordAccess<FieldType::T362>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Expression.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T363>::Value
RecordAccess<FieldType::T363>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 169:
        return ready_record(owner->records_ExpressionDomainReference
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 170:
        return ready_record(owner->records_ExpressionDomainTime
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 171:
        return ready_record(owner->records_ExpressionDomainCurrentPositionTime
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ExpressionDomainReference
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T377>::Value
RecordAccess<FieldType::T377>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 180:
        return ready_record(owner->records_FullShellFieldGlobalPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 181:
        return ready_record(owner->records_FullShellFieldHalfCells
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_FullShellFieldGlobalPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T379>::Value
RecordAccess<FieldType::T379>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_VectorPolynomial.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T384>::Value
RecordAccess<FieldType::T384>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 186:
        return ready_record(
            owner->records_InitialPointShell.values()[static_cast<std::size_t>(typed.record_row)]);
    case 187:
        return ready_record(
            owner->records_InitialPointBeam.values()[static_cast<std::size_t>(typed.record_row)]);
    case 188:
        return ready_record(
            owner->records_InitialPointRigid.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_InitialPointShell.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T385>::Value
RecordAccess<FieldType::T385>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Target.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T386>::Value
RecordAccess<FieldType::T386>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Target.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T387>::Value
RecordAccess<FieldType::T387>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Target.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T389>::Value
RecordAccess<FieldType::T389>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Case.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T390>::Value
RecordAccess<FieldType::T390>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_PolyTerm.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T391>::Value
RecordAccess<FieldType::T391>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 9:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 10:
        return ready_record(owner->records_ReferenceFieldPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 165:
        return ready_record(owner->records_ReferenceFieldExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ReferenceFieldConstant
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T392>::Value
RecordAccess<FieldType::T392>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 11:
        return ready_record(
            owner->records_TimeSeriesConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    case 12:
        return ready_record(owner->records_TimeSeriesPolynomial
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 13:
        return ready_record(
            owner->records_TimeSeriesPchip.values()[static_cast<std::size_t>(typed.record_row)]);
    case 166:
        return ready_record(owner->records_TimeSeriesExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_TimeSeriesConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T393>::Value
RecordAccess<FieldType::T393>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 19:
        return ready_record(
            owner->records_GeometryPg.values()[static_cast<std::size_t>(typed.record_row)]);
    case 24:
        return ready_record(owner->records_GeometrySurfaceGraph
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 25:
        return ready_record(
            owner->records_GeometryNastran.values()[static_cast<std::size_t>(typed.record_row)]);
    case 26:
        return ready_record(
            owner->records_GeometryGmsh.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_GeometryPg.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T394>::Value
RecordAccess<FieldType::T394>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_PgTerm.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T395>::Value
RecordAccess<FieldType::T395>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Boundary.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T396>::Value
RecordAccess<FieldType::T396>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Continuity.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T397>::Value
RecordAccess<FieldType::T397>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 29:
        return ready_record(
            owner->records_MaterialIsotropic.values()[static_cast<std::size_t>(typed.record_row)]);
    case 30:
        return ready_record(owner->records_MaterialOrthotropic
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 31:
        return ready_record(
            owner->records_MaterialEffective.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_MaterialIsotropic.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T398>::Value
RecordAccess<FieldType::T398>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Layer.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T399>::Value
RecordAccess<FieldType::T399>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Shell.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T400>::Value
RecordAccess<FieldType::T400>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(owner->records_Beam.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T401>::Value
RecordAccess<FieldType::T401>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 42:
        return ready_record(
            owner->records_JointAdhesive.values()[static_cast<std::size_t>(typed.record_row)]);
    case 45:
        return ready_record(owner->records_JointRigidAttachment
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 46:
        return ready_record(
            owner->records_JointWeld.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_JointAdhesive.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T402>::Value
RecordAccess<FieldType::T402>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_TraceMap.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T403>::Value
RecordAccess<FieldType::T403>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_MapTerm.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T404>::Value
RecordAccess<FieldType::T404>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 47:
        return ready_record(
            owner->records_SupportConstraint.values()[static_cast<std::size_t>(typed.record_row)]);
    case 48:
        return ready_record(owner->records_SupportTraceConstraint
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 49:
        return ready_record(
            owner->records_SupportElastic.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_SupportConstraint.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T405>::Value
RecordAccess<FieldType::T405>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 50:
        return ready_record(
            owner->records_MassParticle.values()[static_cast<std::size_t>(typed.record_row)]);
    case 51:
        return ready_record(
            owner->records_MassRigid.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_MassParticle.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T406>::Value
RecordAccess<FieldType::T406>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 54:
        return ready_record(
            owner->records_LoadPointForce.values()[static_cast<std::size_t>(typed.record_row)]);
    case 55:
        return ready_record(
            owner->records_LoadCouple.values()[static_cast<std::size_t>(typed.record_row)]);
    case 56:
        return ready_record(
            owner->records_LoadGravity.values()[static_cast<std::size_t>(typed.record_row)]);
    case 57:
        return ready_record(
            owner->records_LoadTraction.values()[static_cast<std::size_t>(typed.record_row)]);
    case 58:
        return ready_record(
            owner->records_LoadPressure.values()[static_cast<std::size_t>(typed.record_row)]);
    case 59:
        return ready_record(
            owner->records_LoadBody.values()[static_cast<std::size_t>(typed.record_row)]);
    case 167:
        return ready_record(owner->records_LoadCurrentExpression
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_LoadPointForce.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T407>::Value
RecordAccess<FieldType::T407>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_AeroSurface.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T408>::Value
RecordAccess<FieldType::T408>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_AeroSharing.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T409>::Value
RecordAccess<FieldType::T409>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_TimeInterval.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T410>::Value
RecordAccess<FieldType::T410>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_SpeedInterval.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T411>::Value
RecordAccess<FieldType::T411>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_CoefficientState.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T412>::Value
RecordAccess<FieldType::T412>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_BeamState.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T413>::Value
RecordAccess<FieldType::T413>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_RigidState.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T414>::Value
RecordAccess<FieldType::T414>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_SurfaceWake.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T415>::Value
RecordAccess<FieldType::T415>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 106:
        return ready_record(
            owner->records_ObservationPoint.values()[static_cast<std::size_t>(typed.record_row)]);
    case 107:
        return ready_record(
            owner->records_ObservationStress.values()[static_cast<std::size_t>(typed.record_row)]);
    case 108:
        return ready_record(
            owner->records_ObservationSection.values()[static_cast<std::size_t>(typed.record_row)]);
    case 109:
        return ready_record(owner->records_ObservationReaction
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 110:
        return ready_record(
            owner->records_ObservationAero.values()[static_cast<std::size_t>(typed.record_row)]);
    case 111:
        return ready_record(
            owner->records_ObservationMode.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_ObservationPoint.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T416>::Value
RecordAccess<FieldType::T416>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 118:
        return ready_record(owner->records_ErrorGroupNotApplicable
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 119:
        return ready_record(owner->records_ErrorGroupCoveredBy
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 120:
        return ready_record(
            owner->records_ErrorGroupRequired.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ErrorGroupNotApplicable
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T417>::Value
RecordAccess<FieldType::T417>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 131:
        return ready_record(owner->records_ScaleUniformComponents
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 132:
        return ready_record(owner->records_ScalePhysicalLayoutPolicy
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    case 137:
        return ready_record(owner->records_ScaleExplicitLayout
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(owner->records_ScaleUniformComponents
                                .values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T418>::Value
RecordAccess<FieldType::T418>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_EigenCluster.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T419>::Value
RecordAccess<FieldType::T419>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_Expression.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T420>::Value
RecordAccess<FieldType::T420>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    switch (typed.record_kind) {
    case 172:
        return ready_record(
            owner->records_ExprNodeConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    case 173:
        return ready_record(
            owner->records_ExprNodeVariable.values()[static_cast<std::size_t>(typed.record_row)]);
    case 174:
        return ready_record(
            owner->records_ExprNodeBinary.values()[static_cast<std::size_t>(typed.record_row)]);
    case 175:
        return ready_record(
            owner->records_ExprNodeNegate.values()[static_cast<std::size_t>(typed.record_row)]);
    case 176:
        return ready_record(
            owner->records_ExprNodePower.values()[static_cast<std::size_t>(typed.record_row)]);
    case 177:
        return ready_record(
            owner->records_ExprNodeAnalytic.values()[static_cast<std::size_t>(typed.record_row)]);
    default:
        return ready_record(
            owner->records_ExprNodeConstant.values()[static_cast<std::size_t>(typed.record_row)]);
    }
}
RecordAccess<FieldType::T421>::Value
RecordAccess<FieldType::T421>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_ShellPhysicalState.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T422>::Value
RecordAccess<FieldType::T422>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_CellField.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T423>::Value
RecordAccess<FieldType::T423>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_VectorTerm.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T424>::Value
RecordAccess<FieldType::T424>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_PointCheck.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T425>::Value
RecordAccess<FieldType::T425>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_FileReference.values()[static_cast<std::size_t>(typed.record_row)]);
}
RecordAccess<FieldType::T444>::Value
RecordAccess<FieldType::T444>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    if (node_kind(owner->snapshot.values()[0], node) != InputNodeKind::Object)
        return schema_number(owner, node);
    return ready_record(owner->records_Rational.values()[static_cast<std::size_t>(
        owner->nodes.values()[static_cast<std::size_t>(node)].record_row)]);
}
RecordAccess<FieldType::T445>::Value
RecordAccess<FieldType::T445>::read(const SchemaData* owner, std::uint64_t node) noexcept {
    const auto& typed = owner->nodes.values()[static_cast<std::size_t>(node)];
    return ready_record(
        owner->records_PolyTerm.values()[static_cast<std::size_t>(typed.record_row)]);
}
} // namespace uvlas::io::detail
