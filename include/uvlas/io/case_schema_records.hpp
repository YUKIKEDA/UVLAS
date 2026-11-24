#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string_view>
#include <uvlas/core/failure.hpp>
#include <uvlas/io/input_document.hpp>
#include <variant>

namespace uvlas::io {
class SiCase;
namespace detail {
struct SchemaData;
struct SchemaFactory;
enum class FieldKind : std::uint8_t {
    Number,
    Integer,
    Identifier,
    Reference,
    Enumeration,
    Text,
    Boolean,
    Array,
    ExactQuantity,
    Object
};
enum class FieldType : std::uint16_t {
    T0,   // Tag{uvlas.case/1}
    T1,   // Units
    T2,   // Model
    T3,   // Flow
    T4,   // Analysis
    T5,   // List<Observation>
    T6,   // Numerics
    T7,   // Resources
    T8,   // Output
    T9,   // Initial
    T10,  // Restart
    T11,  // Provenance
    T12,  // Tag{m,mm,cm,km,in,ft}
    T13,  // Tag{kg,g,t,lbm}
    T14,  // Tag{s,ms,min}
    T15,  // Tag{rad,deg}
    T16,  // List<Geometry>
    T17,  // List<Material>
    T18,  // List<Shell>
    T19,  // List<Beam>
    T20,  // List<Joint>
    T21,  // List<Support>
    T22,  // List<Mass>
    T23,  // List<AeroSurface>
    T24,  // Text
    T25,  // List<Text>
    T26,  // Vec<Q<L>,3>
    T27,  // Mat<Number,3,3>
    T28,  // Tag{shell,beam,rigid}
    T29,  // target-dependent reference to shell, beam or rigid mass
    T30,  // Target
    T31,  // Vec<Number,3>
    T32,  // Tag{interior,lower,upper}
    T33,  // Tag{volume,face,edge,point}
    T34,  // Ref<Geometry>
    T35,  // Index
    T36,  // Vec<Count,k>
    T37,  // Q<d>
    T38,  // Tag{constant}
    T39,  // Vec<Q<d>,n>
    T40,  // Tag{polynomial}
    T41,  // Tag{chart,material,world_reference}
    T42,  // Vec<List<PolyTerm<d,3>>,n>
    T43,  // Vec<Number,n>
    T44,  // Q<T>
    T45,  // List<Vec<Number,n>>
    T46,  // Vec<Q<T>,2>
    T47,  // Tag{pchip}
    T48,  // List<Q<T>>
    T49,  // Outside
    T50,  // Tag{reject,hold}
    T51,  // Tag{periodic}
    T52,  // Frame
    T53,  // TimeSeries<3>
    T54,  // Tag{projected}
    T55,  // ReferenceField<1,3>
    T56,  // ReferenceField<angle,1>
    T57,  // Tag{explicit}
    T58,  // ReferenceField<1,9>
    T59,  // Tag{pg}
    T60,  // Id<Geometry>
    T61,  // Q<L>
    T62,  // Rational
    T63,  // Vec<Count,2>
    T64,  // List<PgTerm>
    T65,  // List<Vec<Rational,2>>
    T66,  // Thickness
    T67,  // Count
    T68,  // Tag{pg-graph-v1}
    T69,  // List<Boundary>
    T70,  // List<Continuity>
    T71,  // List<Vec<Index,3>>
    T72,  // List<Q<L>>
    T73,  // Tag{smooth,fold}
    T74,  // Tag{same,reversed}
    T75,  // Tag{surface_graph}
    T76,  // Tag{nastran}
    T77,  // FileRef<nastran-deck>
    T78,  // FileRef<nastran-sidecar>
    T79,  // Tag{msc-2022.1}
    T80,  // Tag{structural-extraction}
    T81,  // Tag{gmsh}
    T82,  // FileRef<gmsh-msh4.1-ascii>
    T83,  // FileRef<gmsh-sidecar>
    T84,  // ReferenceGeometry
    T85,  // Tag{supplied}
    T86,  // FileRef<reference-geometry>
    T87,  // Tag{fit}
    T88,  // Q<1>
    T89,  // List<Ref<Observation>>
    T90,  // Tag{isotropic}
    T91,  // Id<Material>
    T92,  // Q<pressure>
    T93,  // Q<M/L³>
    T94,  // Applicability
    T95,  // Tag{orthotropic3d}
    T96,  // Vec<Q<pressure>,3>
    T97,  // Vec<Q<1>,3>
    T98,  // Tag{effective3d}
    T99,  // Mat<Q<pressure>,6,6>
    T100, // Tag{mandel,engineering}
    T101, // Tag{unevaluated}
    T102, // Tag{declared}
    T103, // Vec<Number,6>
    T104, // Id<Layer>
    T105, // Vec<Rational,2>
    T106, // Ref<Material>
    T107, // Axis
    T108, // ReferenceField<1,6>
    T109, // Mat<Q<pressure*T>,6,6>
    T110, // Id<Shell>
    T111, // List<Layer>
    T112, // List<Load>
    T113, // Id<Beam>
    T114, // Vec<Frame,2>
    T115, // BeamCoordinate
    T116, // Section
    T117, // Q<energy/L>
    T118, // Mat<Number,6,6>
    T119, // Tag{arclength}
    T120, // Tag{scaled}
    T121, // Tag{matrix}
    T122, // Q<M/L>
    T123, // Mat<Q<M*L>,3,3>
    T124, // Vec<Q<L>,2>
    T125, // Tag{free}
    T126, // Tag{unavailable}
    T127, // Tag{circle}
    T128, // Vec<Q<L²>,2>
    T129, // Tag{rectangle}
    T130, // Q<L⁴>
    T131, // Tag{adhesive}
    T132, // Id<Joint>
    T133, // Region
    T134, // List<TraceMap>
    T135, // Bool
    T136, // Vec<Index,2>
    T137, // Vec<List<MapTerm>,3>
    T138, // Tag{rigid_attachment}
    T139, // Tag{weld}
    T140, // Tag{constraint}
    T141, // Id<Support>
    T142, // MaterialPoint
    T143, // List<Index>
    T144, // Motion
    T145, // Tag{trace_constraint}
    T146, // Tag{elastic}
    T147, // Tag{particle}
    T148, // Id<Mass>
    T149, // Q<M>
    T150, // Ledger
    T151, // Tag{rigid}
    T152, // Mat<Q<M*L²>,3,3>
    T153, // Tag{add}
    T154, // Tag{replace}
    T155, // List<Ref<Mass>>
    T156, // Tag{point_force}
    T157, // Id<Load>
    T158, // Tag{world,body}
    T159, // Tag{couple}
    T160, // Tag{gravity}
    T161, // Tag{traction}
    T162, // Tag{reference,current}
    T163, // ReferenceField<pressure,3>
    T164, // TimeSeries<1>
    T165, // Tag{pressure}
    T166, // ReferenceField<pressure,1>
    T167, // Tag{body_force}
    T168, // ReferenceField<force/L³,3>
    T169, // Id<AeroSurface>
    T170, // AeroAttachment
    T171, // List<Number>
    T172, // Tag{positive,reversed}
    T173, // List<AeroSharing>
    T174, // Tag{fixed_reference}
    T175, // Tag{prescribed_motion}
    T176, // Tag{structural}
    T177, // Ref<Shell>
    T178, // Ref<AeroSurface>
    T179, // Tag{none}
    T180, // Tag{known}
    T181, // Gust
    T182, // FlowReference
    T183, // Tag{sine}
    T184, // Vec<Q<angle/L>,3>
    T185, // Q<angle/T>
    T186, // Q<angle>
    T187, // Q<L/T>
    T188, // Q<L²>
    T189, // Tag{active}
    T190, // Convection
    T191, // Retention
    T192, // Tag{wp1-stream-p0-span}
    T193, // Tag{dual-quarter-v1}
    T194, // Tag{prescribed}
    T195, // Tag{all}
    T196, // Tag{closed_window}
    T197, // Tag{static}
    T198, // Tag{aerodynamic,structural,coupled}
    T199, // Branch
    T200, // AnalysisSpeed
    T201, // Wake
    T202, // Tag{transient}
    T203, // List<TimeInterval>
    T204, // Tag{flutter}
    T205, // Tag{coupled}
    T206, // List<Branch>
    T207, // List<SpeedInterval>
    T208, // Search
    T209, // Tag{value}
    T210, // Id<Branch>
    T211, // BranchSeed
    T212, // Continuation
    T213, // Tag{reference}
    T214, // Tag{prepared}
    T215, // FileRef<preparation-history>
    T216, // Tag{fixed}
    T217, // Tag{sweep}
    T218, // ArcLength
    T219, // PositiveCount
    T220, // List<Scale>
    T221, // Tag{previous_tangent}
    T222, // Tag{initial,event,step_change,chart_change}
    T223, // Tag{lowest_interval_first}
    T224, // Q<1/T>
    T225, // Association
    T226, // Coverage
    T227, // List<Ref<Branch>>
    T228, // Tag{equilibrium}
    T229, // Ref<Branch>
    T230, // Perturbation
    T231, // Tag{cold_start}
    T232, // InitialState
    T233, // List<Q<L²/T>>
    T234, // InitialWake
    T235, // Tag{prepared_history}
    T236, // FileRef<initial-history>
    T237, // List<CoefficientState>
    T238, // List<BeamState>
    T239, // List<RigidState>
    T240, // Digest
    T241, // Tag{absolute,reference_increment}
    T242, // Planes
    T243, // Tag{binary64}
    T244, // Tag{two_word}
    T245, // Ref<Beam>
    T246, // Vec<FrameStatePayload,2>
    T247, // Ref<Mass>
    T248, // FrameStatePayload
    T249, // Vec<Q<L/T>,3>
    T250, // Vec<Q<angle/T>,3>
    T251, // List<SurfaceWake>
    T252, // List<Vec<Q<L>,3>>
    T253, // List<Id<WakeRow>>
    T254, // List<Id<WakePoint>>
    T255, // List<Id<WakeBand>>
    T256, // Tag{same_environment_bitwise,qualified_cross_environment}
    T257, // List<ObservationBudget>
    T258, // Tag{point}
    T259, // Id<Observation>
    T260, // Tag{position,displacement,velocity,acceleration,force}
    T261, // Basis
    T262, // Comparison
    T263, // ObservationBudget
    T264, // Tag{stress}
    T265, // Ref<Layer>
    T266, // Tag{interior,lower,upper,both}
    T267, // Tag{green,second_pk,first_pk,cauchy}
    T268, // Tag{elastic,viscous,total}
    T269, // Tag{section}
    T270, // Tag{force,moment}
    T271, // Tag{reaction}
    T272, // Ref<Support>
    T273, // Tag{aero}
    T274, // List<Ref<AeroSurface>>
    T275, // Tag{force,moment,coefficient}
    T276, // Tag{mode}
    T277, // List<Target>
    T278, // List<MaterialPoint>
    T279, // Tag{frequency,growth,subspace}
    T280, // Tag{world,reference_material,current_material}
    T281, // Tag{instantaneous,history,peak,rms}
    T282, // Tag{left,right,both}
    T283, // ComparisonReference
    T284, // Tag{file}
    T285, // FileRef<observation-reference>
    T286, // List<ErrorGroup>
    T287, // Tag{not_applicable}
    T288, // Tag{geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping}
    T289, // Tag{covered_by}
    T290, // Tag{required}
    T291, // List<Tag{the same twelve error groups}>
    T292, // PatchPolicy
    T293, // Quadrature
    T294, // LinearPolicy
    T295, // NewtonPolicy
    T296, // PrecisionPolicy
    T297, // EigenPolicy
    T298, // WakeCapacity
    T299, // Tag{limited}
    T300, // Tag{exact_no_drop}
    T301, // StructuralQuadrature
    T302, // AeroQuadrature
    T303, // Tag{wp1_induction_v1}
    T304, // Tag{uniform_components}
    T305, // Tag{shell_position,beam_translation,beam_rotation,rigid_translation,rigid_rotation,structural_velocity,structural_angular_velocity,circulation,wake_position,wake_strength,constraint,force_balance,moment_balance,aero_residual,transport_residual,newton_correction,arclength_state}
    T306, // Tag{physical_layout_policy}
    T307, // Tag{mixed_blocks}
    T308, // ConstraintScalePolicy
    T309, // Tag{physical}
    T310, // ComponentScale3
    T311, // Q<energy>
    T312, // Vec<Q<force>,3>
    T313, // Vec<Q<force*L>,3>
    T314, // Tag{values}
    T315, // Tag{explicit_layout}
    T316, // List<Vec<signed-native-integer,4>>
    T317, // Tag{dense_lu,sparse_lu,sparse_llt}
    T318, // LinearAccuracy
    T319, // BackendMemory
    T320, // Tag{gmres}
    T321, // Tag{infinity,l2}
    T322, // Tag{preserve}
    T323, // Tag{adaptive_binary64_two_word}
    T324, // Tag{retained,binary64_projection}
    T325, // Tag{estimate,require_proven_bound}
    T326, // Tag{full}
    T327, // Tag{step_map}
    T328, // Number
    T329, // Tag{spectrum,pairs,subspace}
    T330, // EigenRequirements
    T331, // EigenLimits
    T332, // Tag{generalized}
    T333, // Tag{structural_modes}
    T334, // Q<1/T²>
    T335, // Tag{arnoldi}
    T336, // Tag{pairs,subspace}
    T337, // EigenSelector
    T338, // ArnoldiSeed
    T339, // Tag{uniform_nonzero}
    T340, // Complex
    T341, // List<Complex>
    T342, // Tag{largest_magnitude}
    T343, // Tag{shift_invert}
    T344, // Ref<Observation>
    T345, // Backward
    T346, // List<EigenCluster>
    T347, // Tag{not_requested}
    T348, // Tag{24gib,48gib,custom}
    T349, // InputQuotas
    T350, // OutputCadence
    T351, // Visualization
    T352, // CheckpointOutput
    T353, // Tag{reject_existing}
    T354, // Tag{explicit_only}
    T355, // List<Ref<Layer>>
    T356, // Tag{both}
    T357, // Tag{off}
    T358, // Tag{on}
    T359, // Tag{expression}
    T360, // Expression<n>
    T361, // Tag{current_pressure}
    T362, // Expression<1>
    T363, // ExpressionDomain
    T364, // List<ExprNode>
    T365, // Vec<Index,n>
    T366, // Tag{time}
    T367, // Tag{current_position_time}
    T368, // Vec<signed-native-integer,4>
    T369, // Tag{variable}
    T370, // Tag{add,subtract,multiply,divide}
    T371, // Tag{negate}
    T372, // Tag{integer_power}
    T373, // signed-native-integer
    T374, // Tag{sin,cos,exp,log,sqrt}
    T375, // List<ShellPhysicalState>
    T376, // List<PointCheck>
    T377, // FullShellField
    T378, // Tag{global_polynomial}
    T379, // VectorPolynomial
    T380, // Tag{half_cells}
    T381, // List<CellField>
    T382, // List<VectorTerm>
    T383, // Vec<ExactQuantity,3>
    T384, // InitialPoint
    T385, // Target{kind=shell}
    T386, // Target{kind=beam}
    T387, // Target{kind=rigid}
    T388, // Id<Source>
    T389, // Case
    T390, // PolyTerm
    T391, // ReferenceField
    T392, // TimeSeries
    T393, // Geometry
    T394, // PgTerm
    T395, // Boundary
    T396, // Continuity
    T397, // Material
    T398, // Layer
    T399, // Shell
    T400, // Beam
    T401, // Joint
    T402, // TraceMap
    T403, // MapTerm
    T404, // Support
    T405, // Mass
    T406, // Load
    T407, // AeroSurface
    T408, // AeroSharing
    T409, // TimeInterval
    T410, // SpeedInterval
    T411, // CoefficientState
    T412, // BeamState
    T413, // RigidState
    T414, // SurfaceWake
    T415, // Observation
    T416, // ErrorGroup
    T417, // Scale
    T418, // EigenCluster
    T419, // Expression
    T420, // ExprNode
    T421, // ShellPhysicalState
    T422, // CellField
    T423, // VectorTerm
    T424, // PointCheck
    T425, // FileReference
    T426, // List<PolyTerm<d,3>>
    T427, // Vec<Index,3>
    T428, // Vec<Q<pressure>,6>
    T429, // Q<pressure*T>
    T430, // Vec<Q<pressure*T>,6>
    T431, // Q<M*L>
    T432, // Vec<Q<M*L>,3>
    T433, // List<MapTerm>
    T434, // Q<M*L²>
    T435, // Vec<Q<M*L²>,3>
    T436, // Q<angle/L>
    T437, // Q<L²/T>
    T438, // Id<WakeRow>
    T439, // Id<WakePoint>
    T440, // Id<WakeBand>
    T441, // Tag{the same twelve error groups}
    T442, // Q<force>
    T443, // Q<force*L>
    T444, // ExactQuantity
    T445, // PolyTerm<d,3>
    T446, // Closed external-file role
    T447, // MaterialMap
    T448, // List<MaterialMapPiece>
    T449, // MaterialMapPiece
    T450, // List<MaterialMapTerm>
    T451, // MaterialMapTerm
    T452, // Vec<Rational,3>
};
inline constexpr std::array field_kinds{
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Array,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Array,       FieldKind::Array,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Text,          FieldKind::Array,       FieldKind::Array,
    FieldKind::Array,         FieldKind::Enumeration, FieldKind::Reference,
    FieldKind::Object,        FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Reference,   FieldKind::Integer,
    FieldKind::Array,         FieldKind::Number,      FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Number,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Identifier,    FieldKind::Number,      FieldKind::Object,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Object,        FieldKind::Integer,     FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Array,         FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Object,        FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Number,      FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Identifier,  FieldKind::Number,
    FieldKind::Number,        FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Array,       FieldKind::Identifier,
    FieldKind::Array,         FieldKind::Reference,   FieldKind::Object,
    FieldKind::Object,        FieldKind::Array,       FieldKind::Identifier,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Identifier,
    FieldKind::Array,         FieldKind::Object,      FieldKind::Object,
    FieldKind::Number,        FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Number,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Number,      FieldKind::Enumeration,
    FieldKind::Identifier,    FieldKind::Object,      FieldKind::Array,
    FieldKind::Boolean,       FieldKind::Array,       FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Identifier,    FieldKind::Object,      FieldKind::Array,
    FieldKind::Object,        FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Identifier,  FieldKind::Number,
    FieldKind::Object,        FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Identifier,  FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Object,        FieldKind::Identifier,  FieldKind::Object,
    FieldKind::Array,         FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Reference,     FieldKind::Reference,   FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Array,       FieldKind::Number,
    FieldKind::Number,        FieldKind::Number,      FieldKind::Number,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Array,         FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Identifier,    FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Integer,       FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Number,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Reference,   FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Array,
    FieldKind::Object,        FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Text,          FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Reference,
    FieldKind::Array,         FieldKind::Reference,   FieldKind::Object,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Array,         FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Identifier,  FieldKind::Enumeration,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Reference,   FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Reference,
    FieldKind::Enumeration,   FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Array,       FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Object,        FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Number,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Number,      FieldKind::Enumeration,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Number,      FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Reference,
    FieldKind::Object,        FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Array,       FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Object,        FieldKind::Enumeration, FieldKind::Object,
    FieldKind::Object,        FieldKind::Array,       FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Array,
    FieldKind::Enumeration,   FieldKind::Enumeration, FieldKind::Enumeration,
    FieldKind::Enumeration,   FieldKind::Integer,     FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Object,
    FieldKind::Enumeration,   FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Identifier,  FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Object,        FieldKind::Object,      FieldKind::Object,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Array,
    FieldKind::Number,        FieldKind::Array,       FieldKind::Number,
    FieldKind::Array,         FieldKind::Array,       FieldKind::Number,
    FieldKind::Array,         FieldKind::Number,      FieldKind::Number,
    FieldKind::Identifier,    FieldKind::Identifier,  FieldKind::Identifier,
    FieldKind::Enumeration,   FieldKind::Number,      FieldKind::Number,
    FieldKind::ExactQuantity, FieldKind::Object,      FieldKind::Enumeration,
    FieldKind::Object,        FieldKind::Array,       FieldKind::Object,
    FieldKind::Array,         FieldKind::Object,      FieldKind::Array,
};
inline constexpr std::array field_is_signed{
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, true,  false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
};
template <FieldType>
struct Enumeration {
    enum class Value : std::uint8_t {};
};
template <>
struct Enumeration<FieldType::T446> {
    enum class Value : std::uint8_t {
        ItemNastranDeck,
        ItemNastranSidecar,
        ItemGmshMsh41Ascii,
        ItemGmshSidecar,
        ItemReferenceGeometry,
        ItemPreparationHistory,
        ItemInitialHistory,
        ItemObservationReference
    };
};
template <>
struct Enumeration<FieldType::T0> {
    enum class Value : std::uint8_t { ItemUvlasCase1 };
};
template <>
struct Enumeration<FieldType::T12> {
    enum class Value : std::uint8_t { ItemM, ItemMm, ItemCm, ItemKm, ItemIn, ItemFt };
};
template <>
struct Enumeration<FieldType::T13> {
    enum class Value : std::uint8_t { ItemKg, ItemG, ItemT, ItemLbm };
};
template <>
struct Enumeration<FieldType::T14> {
    enum class Value : std::uint8_t { ItemS, ItemMs, ItemMin };
};
template <>
struct Enumeration<FieldType::T15> {
    enum class Value : std::uint8_t { ItemRad, ItemDeg };
};
template <>
struct Enumeration<FieldType::T28> {
    enum class Value : std::uint8_t { ItemShell, ItemBeam, ItemRigid };
};
template <>
struct Enumeration<FieldType::T32> {
    enum class Value : std::uint8_t { ItemInterior, ItemLower, ItemUpper };
};
template <>
struct Enumeration<FieldType::T33> {
    enum class Value : std::uint8_t { ItemVolume, ItemFace, ItemEdge, ItemPoint };
};
template <>
struct Enumeration<FieldType::T38> {
    enum class Value : std::uint8_t { ItemConstant };
};
template <>
struct Enumeration<FieldType::T40> {
    enum class Value : std::uint8_t { ItemPolynomial };
};
template <>
struct Enumeration<FieldType::T41> {
    enum class Value : std::uint8_t { ItemChart, ItemMaterial, ItemWorldReference };
};
template <>
struct Enumeration<FieldType::T47> {
    enum class Value : std::uint8_t { ItemPchip };
};
template <>
struct Enumeration<FieldType::T50> {
    enum class Value : std::uint8_t { ItemReject, ItemHold };
};
template <>
struct Enumeration<FieldType::T51> {
    enum class Value : std::uint8_t { ItemPeriodic };
};
template <>
struct Enumeration<FieldType::T54> {
    enum class Value : std::uint8_t { ItemProjected };
};
template <>
struct Enumeration<FieldType::T57> {
    enum class Value : std::uint8_t { ItemExplicit };
};
template <>
struct Enumeration<FieldType::T59> {
    enum class Value : std::uint8_t { ItemPg };
};
template <>
struct Enumeration<FieldType::T68> {
    enum class Value : std::uint8_t { ItemPgGraphV1 };
};
template <>
struct Enumeration<FieldType::T73> {
    enum class Value : std::uint8_t { ItemSmooth, ItemFold };
};
template <>
struct Enumeration<FieldType::T74> {
    enum class Value : std::uint8_t { ItemSame, ItemReversed };
};
template <>
struct Enumeration<FieldType::T75> {
    enum class Value : std::uint8_t { ItemSurfaceGraph };
};
template <>
struct Enumeration<FieldType::T76> {
    enum class Value : std::uint8_t { ItemNastran };
};
template <>
struct Enumeration<FieldType::T79> {
    enum class Value : std::uint8_t { ItemMsc20221 };
};
template <>
struct Enumeration<FieldType::T80> {
    enum class Value : std::uint8_t { ItemStructuralExtraction };
};
template <>
struct Enumeration<FieldType::T81> {
    enum class Value : std::uint8_t { ItemGmsh };
};
template <>
struct Enumeration<FieldType::T85> {
    enum class Value : std::uint8_t { ItemSupplied };
};
template <>
struct Enumeration<FieldType::T87> {
    enum class Value : std::uint8_t { ItemFit };
};
template <>
struct Enumeration<FieldType::T90> {
    enum class Value : std::uint8_t { ItemIsotropic };
};
template <>
struct Enumeration<FieldType::T95> {
    enum class Value : std::uint8_t { ItemOrthotropic3D };
};
template <>
struct Enumeration<FieldType::T98> {
    enum class Value : std::uint8_t { ItemEffective3D };
};
template <>
struct Enumeration<FieldType::T100> {
    enum class Value : std::uint8_t { ItemMandel, ItemEngineering };
};
template <>
struct Enumeration<FieldType::T101> {
    enum class Value : std::uint8_t { ItemUnevaluated };
};
template <>
struct Enumeration<FieldType::T102> {
    enum class Value : std::uint8_t { ItemDeclared };
};
template <>
struct Enumeration<FieldType::T119> {
    enum class Value : std::uint8_t { ItemArclength };
};
template <>
struct Enumeration<FieldType::T120> {
    enum class Value : std::uint8_t { ItemScaled };
};
template <>
struct Enumeration<FieldType::T121> {
    enum class Value : std::uint8_t { ItemMatrix };
};
template <>
struct Enumeration<FieldType::T125> {
    enum class Value : std::uint8_t { ItemFree };
};
template <>
struct Enumeration<FieldType::T126> {
    enum class Value : std::uint8_t { ItemUnavailable };
};
template <>
struct Enumeration<FieldType::T127> {
    enum class Value : std::uint8_t { ItemCircle };
};
template <>
struct Enumeration<FieldType::T129> {
    enum class Value : std::uint8_t { ItemRectangle };
};
template <>
struct Enumeration<FieldType::T131> {
    enum class Value : std::uint8_t { ItemAdhesive };
};
template <>
struct Enumeration<FieldType::T138> {
    enum class Value : std::uint8_t { ItemRigidAttachment };
};
template <>
struct Enumeration<FieldType::T139> {
    enum class Value : std::uint8_t { ItemWeld };
};
template <>
struct Enumeration<FieldType::T140> {
    enum class Value : std::uint8_t { ItemConstraint };
};
template <>
struct Enumeration<FieldType::T145> {
    enum class Value : std::uint8_t { ItemTraceConstraint };
};
template <>
struct Enumeration<FieldType::T146> {
    enum class Value : std::uint8_t { ItemElastic };
};
template <>
struct Enumeration<FieldType::T147> {
    enum class Value : std::uint8_t { ItemParticle };
};
template <>
struct Enumeration<FieldType::T151> {
    enum class Value : std::uint8_t { ItemRigid };
};
template <>
struct Enumeration<FieldType::T153> {
    enum class Value : std::uint8_t { ItemAdd };
};
template <>
struct Enumeration<FieldType::T154> {
    enum class Value : std::uint8_t { ItemReplace };
};
template <>
struct Enumeration<FieldType::T156> {
    enum class Value : std::uint8_t { ItemPointForce };
};
template <>
struct Enumeration<FieldType::T158> {
    enum class Value : std::uint8_t { ItemWorld, ItemBody };
};
template <>
struct Enumeration<FieldType::T159> {
    enum class Value : std::uint8_t { ItemCouple };
};
template <>
struct Enumeration<FieldType::T160> {
    enum class Value : std::uint8_t { ItemGravity };
};
template <>
struct Enumeration<FieldType::T161> {
    enum class Value : std::uint8_t { ItemTraction };
};
template <>
struct Enumeration<FieldType::T162> {
    enum class Value : std::uint8_t { ItemReference, ItemCurrent };
};
template <>
struct Enumeration<FieldType::T165> {
    enum class Value : std::uint8_t { ItemPressure };
};
template <>
struct Enumeration<FieldType::T167> {
    enum class Value : std::uint8_t { ItemBodyForce };
};
template <>
struct Enumeration<FieldType::T172> {
    enum class Value : std::uint8_t { ItemPositive, ItemReversed };
};
template <>
struct Enumeration<FieldType::T174> {
    enum class Value : std::uint8_t { ItemFixedReference };
};
template <>
struct Enumeration<FieldType::T175> {
    enum class Value : std::uint8_t { ItemPrescribedMotion };
};
template <>
struct Enumeration<FieldType::T176> {
    enum class Value : std::uint8_t { ItemStructural };
};
template <>
struct Enumeration<FieldType::T179> {
    enum class Value : std::uint8_t { ItemNone };
};
template <>
struct Enumeration<FieldType::T180> {
    enum class Value : std::uint8_t { ItemKnown };
};
template <>
struct Enumeration<FieldType::T183> {
    enum class Value : std::uint8_t { ItemSine };
};
template <>
struct Enumeration<FieldType::T189> {
    enum class Value : std::uint8_t { ItemActive };
};
template <>
struct Enumeration<FieldType::T192> {
    enum class Value : std::uint8_t { ItemWp1StreamP0Span };
};
template <>
struct Enumeration<FieldType::T193> {
    enum class Value : std::uint8_t { ItemDualQuarterV1 };
};
template <>
struct Enumeration<FieldType::T194> {
    enum class Value : std::uint8_t { ItemPrescribed };
};
template <>
struct Enumeration<FieldType::T195> {
    enum class Value : std::uint8_t { ItemAll };
};
template <>
struct Enumeration<FieldType::T196> {
    enum class Value : std::uint8_t { ItemClosedWindow };
};
template <>
struct Enumeration<FieldType::T197> {
    enum class Value : std::uint8_t { ItemStatic };
};
template <>
struct Enumeration<FieldType::T198> {
    enum class Value : std::uint8_t { ItemAerodynamic, ItemStructural, ItemCoupled };
};
template <>
struct Enumeration<FieldType::T202> {
    enum class Value : std::uint8_t { ItemTransient };
};
template <>
struct Enumeration<FieldType::T204> {
    enum class Value : std::uint8_t { ItemFlutter };
};
template <>
struct Enumeration<FieldType::T205> {
    enum class Value : std::uint8_t { ItemCoupled };
};
template <>
struct Enumeration<FieldType::T209> {
    enum class Value : std::uint8_t { ItemValue };
};
template <>
struct Enumeration<FieldType::T213> {
    enum class Value : std::uint8_t { ItemReference };
};
template <>
struct Enumeration<FieldType::T214> {
    enum class Value : std::uint8_t { ItemPrepared };
};
template <>
struct Enumeration<FieldType::T216> {
    enum class Value : std::uint8_t { ItemFixed };
};
template <>
struct Enumeration<FieldType::T217> {
    enum class Value : std::uint8_t { ItemSweep };
};
template <>
struct Enumeration<FieldType::T221> {
    enum class Value : std::uint8_t { ItemPreviousTangent };
};
template <>
struct Enumeration<FieldType::T222> {
    enum class Value : std::uint8_t { ItemInitial, ItemEvent, ItemStepChange, ItemChartChange };
};
template <>
struct Enumeration<FieldType::T223> {
    enum class Value : std::uint8_t { ItemLowestIntervalFirst };
};
template <>
struct Enumeration<FieldType::T228> {
    enum class Value : std::uint8_t { ItemEquilibrium };
};
template <>
struct Enumeration<FieldType::T231> {
    enum class Value : std::uint8_t { ItemColdStart };
};
template <>
struct Enumeration<FieldType::T235> {
    enum class Value : std::uint8_t { ItemPreparedHistory };
};
template <>
struct Enumeration<FieldType::T241> {
    enum class Value : std::uint8_t { ItemAbsolute, ItemReferenceIncrement };
};
template <>
struct Enumeration<FieldType::T243> {
    enum class Value : std::uint8_t { ItemBinary64 };
};
template <>
struct Enumeration<FieldType::T244> {
    enum class Value : std::uint8_t { ItemTwoWord };
};
template <>
struct Enumeration<FieldType::T256> {
    enum class Value : std::uint8_t { ItemSameEnvironmentBitwise, ItemQualifiedCrossEnvironment };
};
template <>
struct Enumeration<FieldType::T258> {
    enum class Value : std::uint8_t { ItemPoint };
};
template <>
struct Enumeration<FieldType::T260> {
    enum class Value : std::uint8_t {
        ItemPosition,
        ItemDisplacement,
        ItemVelocity,
        ItemAcceleration,
        ItemForce
    };
};
template <>
struct Enumeration<FieldType::T264> {
    enum class Value : std::uint8_t { ItemStress };
};
template <>
struct Enumeration<FieldType::T266> {
    enum class Value : std::uint8_t { ItemInterior, ItemLower, ItemUpper, ItemBoth };
};
template <>
struct Enumeration<FieldType::T267> {
    enum class Value : std::uint8_t { ItemGreen, ItemSecondPk, ItemFirstPk, ItemCauchy };
};
template <>
struct Enumeration<FieldType::T268> {
    enum class Value : std::uint8_t { ItemElastic, ItemViscous, ItemTotal };
};
template <>
struct Enumeration<FieldType::T269> {
    enum class Value : std::uint8_t { ItemSection };
};
template <>
struct Enumeration<FieldType::T270> {
    enum class Value : std::uint8_t { ItemForce, ItemMoment };
};
template <>
struct Enumeration<FieldType::T271> {
    enum class Value : std::uint8_t { ItemReaction };
};
template <>
struct Enumeration<FieldType::T273> {
    enum class Value : std::uint8_t { ItemAero };
};
template <>
struct Enumeration<FieldType::T275> {
    enum class Value : std::uint8_t { ItemForce, ItemMoment, ItemCoefficient };
};
template <>
struct Enumeration<FieldType::T276> {
    enum class Value : std::uint8_t { ItemMode };
};
template <>
struct Enumeration<FieldType::T279> {
    enum class Value : std::uint8_t { ItemFrequency, ItemGrowth, ItemSubspace };
};
template <>
struct Enumeration<FieldType::T280> {
    enum class Value : std::uint8_t { ItemWorld, ItemReferenceMaterial, ItemCurrentMaterial };
};
template <>
struct Enumeration<FieldType::T281> {
    enum class Value : std::uint8_t { ItemInstantaneous, ItemHistory, ItemPeak, ItemRms };
};
template <>
struct Enumeration<FieldType::T282> {
    enum class Value : std::uint8_t { ItemLeft, ItemRight, ItemBoth };
};
template <>
struct Enumeration<FieldType::T284> {
    enum class Value : std::uint8_t { ItemFile };
};
template <>
struct Enumeration<FieldType::T287> {
    enum class Value : std::uint8_t { ItemNotApplicable };
};
template <>
struct Enumeration<FieldType::T288> {
    enum class Value : std::uint8_t {
        ItemGeometry,
        ItemStructureGrid,
        ItemAeroGrid,
        ItemTimeWakeAge,
        ItemWakeLength,
        ItemCore,
        ItemQuadrature,
        ItemAlgebra,
        ItemNewton,
        ItemEigenRange,
        ItemSearchWidth,
        ItemLoadMapping
    };
};
template <>
struct Enumeration<FieldType::T289> {
    enum class Value : std::uint8_t { ItemCoveredBy };
};
template <>
struct Enumeration<FieldType::T290> {
    enum class Value : std::uint8_t { ItemRequired };
};
template <>
struct Enumeration<FieldType::T299> {
    enum class Value : std::uint8_t { ItemLimited };
};
template <>
struct Enumeration<FieldType::T300> {
    enum class Value : std::uint8_t { ItemExactNoDrop };
};
template <>
struct Enumeration<FieldType::T303> {
    enum class Value : std::uint8_t { ItemWp1InductionV1 };
};
template <>
struct Enumeration<FieldType::T304> {
    enum class Value : std::uint8_t { ItemUniformComponents };
};
template <>
struct Enumeration<FieldType::T305> {
    enum class Value : std::uint8_t {
        ItemShellPosition,
        ItemBeamTranslation,
        ItemBeamRotation,
        ItemRigidTranslation,
        ItemRigidRotation,
        ItemStructuralVelocity,
        ItemStructuralAngularVelocity,
        ItemCirculation,
        ItemWakePosition,
        ItemWakeStrength,
        ItemConstraint,
        ItemForceBalance,
        ItemMomentBalance,
        ItemAeroResidual,
        ItemTransportResidual,
        ItemNewtonCorrection,
        ItemArclengthState
    };
};
template <>
struct Enumeration<FieldType::T306> {
    enum class Value : std::uint8_t { ItemPhysicalLayoutPolicy };
};
template <>
struct Enumeration<FieldType::T307> {
    enum class Value : std::uint8_t { ItemMixedBlocks };
};
template <>
struct Enumeration<FieldType::T309> {
    enum class Value : std::uint8_t { ItemPhysical };
};
template <>
struct Enumeration<FieldType::T314> {
    enum class Value : std::uint8_t { ItemValues };
};
template <>
struct Enumeration<FieldType::T315> {
    enum class Value : std::uint8_t { ItemExplicitLayout };
};
template <>
struct Enumeration<FieldType::T317> {
    enum class Value : std::uint8_t { ItemDenseLu, ItemSparseLu, ItemSparseLlt };
};
template <>
struct Enumeration<FieldType::T320> {
    enum class Value : std::uint8_t { ItemGmres };
};
template <>
struct Enumeration<FieldType::T321> {
    enum class Value : std::uint8_t { ItemInfinity, ItemL2 };
};
template <>
struct Enumeration<FieldType::T322> {
    enum class Value : std::uint8_t { ItemPreserve };
};
template <>
struct Enumeration<FieldType::T323> {
    enum class Value : std::uint8_t { ItemAdaptiveBinary64TwoWord };
};
template <>
struct Enumeration<FieldType::T324> {
    enum class Value : std::uint8_t { ItemRetained, ItemBinary64Projection };
};
template <>
struct Enumeration<FieldType::T325> {
    enum class Value : std::uint8_t { ItemEstimate, ItemRequireProvenBound };
};
template <>
struct Enumeration<FieldType::T326> {
    enum class Value : std::uint8_t { ItemFull };
};
template <>
struct Enumeration<FieldType::T327> {
    enum class Value : std::uint8_t { ItemStepMap };
};
template <>
struct Enumeration<FieldType::T329> {
    enum class Value : std::uint8_t { ItemSpectrum, ItemPairs, ItemSubspace };
};
template <>
struct Enumeration<FieldType::T332> {
    enum class Value : std::uint8_t { ItemGeneralized };
};
template <>
struct Enumeration<FieldType::T333> {
    enum class Value : std::uint8_t { ItemStructuralModes };
};
template <>
struct Enumeration<FieldType::T335> {
    enum class Value : std::uint8_t { ItemArnoldi };
};
template <>
struct Enumeration<FieldType::T336> {
    enum class Value : std::uint8_t { ItemPairs, ItemSubspace };
};
template <>
struct Enumeration<FieldType::T339> {
    enum class Value : std::uint8_t { ItemUniformNonzero };
};
template <>
struct Enumeration<FieldType::T342> {
    enum class Value : std::uint8_t { ItemLargestMagnitude };
};
template <>
struct Enumeration<FieldType::T343> {
    enum class Value : std::uint8_t { ItemShiftInvert };
};
template <>
struct Enumeration<FieldType::T347> {
    enum class Value : std::uint8_t { ItemNotRequested };
};
template <>
struct Enumeration<FieldType::T348> {
    enum class Value : std::uint8_t { Item24Gib, Item48Gib, ItemCustom };
};
template <>
struct Enumeration<FieldType::T353> {
    enum class Value : std::uint8_t { ItemRejectExisting };
};
template <>
struct Enumeration<FieldType::T354> {
    enum class Value : std::uint8_t { ItemExplicitOnly };
};
template <>
struct Enumeration<FieldType::T356> {
    enum class Value : std::uint8_t { ItemBoth };
};
template <>
struct Enumeration<FieldType::T357> {
    enum class Value : std::uint8_t { ItemOff };
};
template <>
struct Enumeration<FieldType::T358> {
    enum class Value : std::uint8_t { ItemOn };
};
template <>
struct Enumeration<FieldType::T359> {
    enum class Value : std::uint8_t { ItemExpression };
};
template <>
struct Enumeration<FieldType::T361> {
    enum class Value : std::uint8_t { ItemCurrentPressure };
};
template <>
struct Enumeration<FieldType::T366> {
    enum class Value : std::uint8_t { ItemTime };
};
template <>
struct Enumeration<FieldType::T367> {
    enum class Value : std::uint8_t { ItemCurrentPositionTime };
};
template <>
struct Enumeration<FieldType::T369> {
    enum class Value : std::uint8_t { ItemVariable };
};
template <>
struct Enumeration<FieldType::T370> {
    enum class Value : std::uint8_t { ItemAdd, ItemSubtract, ItemMultiply, ItemDivide };
};
template <>
struct Enumeration<FieldType::T371> {
    enum class Value : std::uint8_t { ItemNegate };
};
template <>
struct Enumeration<FieldType::T372> {
    enum class Value : std::uint8_t { ItemIntegerPower };
};
template <>
struct Enumeration<FieldType::T374> {
    enum class Value : std::uint8_t { ItemSin, ItemCos, ItemExp, ItemLog, ItemSqrt };
};
template <>
struct Enumeration<FieldType::T378> {
    enum class Value : std::uint8_t { ItemGlobalPolynomial };
};
template <>
struct Enumeration<FieldType::T380> {
    enum class Value : std::uint8_t { ItemHalfCells };
};
template <>
struct Enumeration<FieldType::T441> {
    enum class Value : std::uint8_t {
        ItemGeometry,
        ItemStructureGrid,
        ItemAeroGrid,
        ItemTimeWakeAge,
        ItemWakeLength,
        ItemCore,
        ItemQuadrature,
        ItemAlgebra,
        ItemNewton,
        ItemEigenRange,
        ItemSearchWidth,
        ItemLoadMapping
    };
};
inline constexpr std::array field_elements{
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T415, FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T393, FieldType::T397, FieldType::T399, FieldType::T400,
    FieldType::T401, FieldType::T404, FieldType::T405, FieldType::T407, FieldType::T0,
    FieldType::T24,  FieldType::T61,  FieldType::T31,  FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T328, FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T67,  FieldType::T0,   FieldType::T0,   FieldType::T37,
    FieldType::T0,   FieldType::T0,   FieldType::T426, FieldType::T328, FieldType::T0,
    FieldType::T43,  FieldType::T44,  FieldType::T0,   FieldType::T44,  FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T67,  FieldType::T394,
    FieldType::T105, FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T395,
    FieldType::T396, FieldType::T427, FieldType::T61,  FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T344,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T92,  FieldType::T88,  FieldType::T0,   FieldType::T428,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T328, FieldType::T0,
    FieldType::T62,  FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T430,
    FieldType::T0,   FieldType::T398, FieldType::T406, FieldType::T0,   FieldType::T52,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T103, FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T432, FieldType::T61,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T188, FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T402,
    FieldType::T0,   FieldType::T35,  FieldType::T433, FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T35,  FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T435, FieldType::T0,   FieldType::T0,
    FieldType::T247, FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T328, FieldType::T0,   FieldType::T408, FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T436,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T409, FieldType::T0,
    FieldType::T0,   FieldType::T199, FieldType::T410, FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T417, FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T229, FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T437, FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T411, FieldType::T412, FieldType::T413,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T248, FieldType::T0,   FieldType::T0,   FieldType::T187,
    FieldType::T185, FieldType::T414, FieldType::T26,  FieldType::T438, FieldType::T439,
    FieldType::T440, FieldType::T0,   FieldType::T263, FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T178,
    FieldType::T0,   FieldType::T0,   FieldType::T30,  FieldType::T142, FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T416, FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T441, FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T442, FieldType::T443, FieldType::T0,
    FieldType::T0,   FieldType::T368, FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T340, FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T418, FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T265, FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T420,
    FieldType::T35,  FieldType::T0,   FieldType::T0,   FieldType::T373, FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T421, FieldType::T424, FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T422, FieldType::T423, FieldType::T444, FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T445, FieldType::T35,  FieldType::T92,  FieldType::T0,
    FieldType::T429, FieldType::T0,   FieldType::T431, FieldType::T403, FieldType::T0,
    FieldType::T434, FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T0,
    FieldType::T0,   FieldType::T0,   FieldType::T0,   FieldType::T449, FieldType::T0,
    FieldType::T451, FieldType::T0,   FieldType::T62,
};
template <FieldType>
struct RecordAccess;
std::string_view schema_text(const SchemaData*, std::uint64_t) noexcept;
std::uint64_t schema_count(const SchemaData*, std::uint64_t) noexcept;
std::int64_t schema_integer(const SchemaData*, std::uint64_t) noexcept;
bool schema_boolean(const SchemaData*, std::uint64_t) noexcept;
core::Result<std::uint64_t> schema_child(const SchemaData*, std::uint64_t, core::Index) noexcept;
} // namespace detail

/// Four exact unit exponents, ordered length, mass, time and angle.
struct SchemaDimension {
    std::array<std::int16_t, 4>
        exponents; ///< Checked exponents; angles remain distinct from dimensionless values.
    bool operator==(const SchemaDimension&) const = default;
};
/// Read-only exact decimal value. Views expire when the owning case is destroyed or replaced.
struct SchemaNumberView {
    int sign;                    ///< Minus one, zero or plus one; signed zero has value zero.
    std::string_view digits;     ///< Normalized ASCII significand, preserving every input digit.
    std::int64_t exponent10;     ///< Value is sign times digits times ten to this exponent.
    InputNodeKind original_kind; ///< Native integer and decimal input remain distinguishable.
    SchemaDimension dimensions;  ///< Consumer-assigned dimensions before SI conversion.
};
namespace detail {
SchemaNumberView schema_number(const SchemaData*, std::uint64_t) noexcept;
}
namespace schema {
/// Immutable owner-bound reference to one fully checked field. It cannot be built from an ordinal.
/// Keep the owning SchemaCase alive; destruction or assignment invalidates all its borrowed fields.
/// A field is never looked up against a separately supplied owner.
template <detail::FieldType Type>
class Field {
public:
    /// Return the exact numeric value without converting through binary floating point.
    [[nodiscard]] SchemaNumberView number() const noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Number)
    {
        return detail::schema_number(owner_, node_);
    }
    /// Borrow decoded text or the selected enumeration spelling until the owner changes.
    [[nodiscard]] std::string_view text() const noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Text ||
                 detail::field_kinds[static_cast<std::size_t>(Type)] ==
                     detail::FieldKind::Enumeration)
    {
        return detail::schema_text(owner_, node_);
    }
    /// Read a checked nonnegative count, identifier or resolved reference's user ID.
    [[nodiscard]] std::uint64_t unsigned_value() const noexcept
        requires(
            (detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Integer &&
             !detail::field_is_signed[static_cast<std::size_t>(Type)]) ||
            detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Identifier ||
            detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Reference)
    {
        return detail::schema_count(owner_, node_);
    }
    /// Read a checked signed integer without narrowing or changing its JSON kind.
    [[nodiscard]] std::int64_t signed_value() const noexcept
        requires(detail::field_is_signed[static_cast<std::size_t>(Type)])
    {
        return detail::schema_integer(owner_, node_);
    }
    /// Return the selected member of this field's closed enumeration type.
    /// Enumerations from different fields cannot be implicitly interchanged.
    [[nodiscard]] typename detail::Enumeration<Type>::Value selected() const noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] ==
                 detail::FieldKind::Enumeration)
    {
        return static_cast<typename detail::Enumeration<Type>::Value>(
            detail::schema_count(owner_, node_));
    }
    /// Read the JSON boolean without numeric coercion.
    [[nodiscard]] bool boolean() const noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Boolean)
    {
        return detail::schema_boolean(owner_, node_);
    }
    /// Number of checked elements; empty arrays remain distinct from absent optional fields.
    [[nodiscard]] std::uint64_t size() const noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Array)
    {
        return detail::schema_count(owner_, node_);
    }
    /// Borrow a checked element of its statically declared type; reject negative/out-of-range
    /// indices.
    [[nodiscard]] auto at(core::Index index) const noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Array)
    {
        using Element = Field<detail::field_elements[static_cast<std::size_t>(Type)]>;
        const auto child = detail::schema_child(owner_, node_, index);
        if (!child)
            return core::Result<Element>{std::unexpected(child.error())};
        return core::Result<Element>{Element{owner_, *child}};
    }
    /// Return the named record or closed variant for this checked object, bound to the same owner.
    [[nodiscard]] auto record() const noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Object ||
                 detail::field_kinds[static_cast<std::size_t>(Type)] ==
                     detail::FieldKind::ExactQuantity)
    {
        return detail::RecordAccess<Type>::read(owner_, node_);
    }

private:
    friend struct detail::SchemaFactory;
    friend class uvlas::io::SiCase;
    template <detail::FieldType>
    friend class Field;
    Field(const detail::SchemaData* owner, std::uint64_t node) noexcept
        : owner_(owner), node_(node) {}
    const detail::SchemaData* owner_;
    std::uint64_t node_;
};
/// Checked Case payload; all fields retain one immutable case owner.
struct Case {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Case(Field<detail::FieldType::T0> schema_value, Field<detail::FieldType::T1> units_value,
         Field<detail::FieldType::T2> model_value, Field<detail::FieldType::T3> flow_value,
         Field<detail::FieldType::T4> analysis_value,
         Field<detail::FieldType::T5> observations_value,
         Field<detail::FieldType::T6> numerics_value, Field<detail::FieldType::T7> resources_value,
         Field<detail::FieldType::T8> output_value,
         std::optional<Field<detail::FieldType::T9>> initial_value,
         std::optional<Field<detail::FieldType::T10>> restart_value,
         std::optional<Field<detail::FieldType::T11>> provenance_value) noexcept
        : schema(schema_value),
          units(units_value),
          model(model_value),
          flow(flow_value),
          analysis(analysis_value),
          observations(observations_value),
          numerics(numerics_value),
          resources(resources_value),
          output(output_value),
          initial(initial_value),
          restart(restart_value),
          provenance(provenance_value) {}
    Field<detail::FieldType::T0> schema; ///< Schema; choose one of uvlas.case/1.
    Field<detail::FieldType::T1>
        units; ///< Units payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T2>
        model; ///< Model payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T3>
        flow; ///< Flow payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T4>
        analysis; ///< Analysis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T5> observations; ///< Ordered observations entries; size() and at()
                                               ///< expose the checked sequence without copying.
    Field<detail::FieldType::T6>
        numerics; ///< Numerics payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T7>
        resources; ///< Resources payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T8>
        output; ///< Output payload; record() selects its checked concrete structure.
    std::optional<Field<detail::FieldType::T9>>
        initial; ///< Initial payload; record() selects its checked concrete structure.
    std::optional<Field<detail::FieldType::T10>>
        restart; ///< Restart payload; record() selects its checked concrete structure.
    std::optional<Field<detail::FieldType::T11>>
        provenance; ///< Provenance payload; record() selects its checked concrete structure.
};
/// Checked Units payload; all fields retain one immutable case owner.
struct Units {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Units(Field<detail::FieldType::T12> length_value, Field<detail::FieldType::T13> mass_value,
          Field<detail::FieldType::T14> time_value,
          Field<detail::FieldType::T15> angle_value) noexcept
        : length(length_value), mass(mass_value), time(time_value), angle(angle_value) {}
    Field<detail::FieldType::T12> length; ///< Length; choose one of m,mm,cm,km,in,ft.
    Field<detail::FieldType::T13> mass;   ///< Mass; choose one of kg,g,t,lbm.
    Field<detail::FieldType::T14> time;   ///< Time; choose one of s,ms,min.
    Field<detail::FieldType::T15> angle;  ///< Angle; choose one of rad,deg.
};
/// Checked Model payload; all fields retain one immutable case owner.
struct Model {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Model(Field<detail::FieldType::T16> geometry_value,
          Field<detail::FieldType::T17> materials_value, Field<detail::FieldType::T18> shells_value,
          Field<detail::FieldType::T19> beams_value, Field<detail::FieldType::T20> joints_value,
          Field<detail::FieldType::T21> supports_value, Field<detail::FieldType::T22> masses_value,
          Field<detail::FieldType::T23> aero_surfaces_value) noexcept
        : geometry(geometry_value),
          materials(materials_value),
          shells(shells_value),
          beams(beams_value),
          joints(joints_value),
          supports(supports_value),
          masses(masses_value),
          aero_surfaces(aero_surfaces_value) {}
    Field<detail::FieldType::T16> geometry;  ///< Ordered geometry entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T17> materials; ///< Ordered materials entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T18> shells;    ///< Ordered shells entries; size() and at() expose the
                                             ///< checked sequence without copying.
    Field<detail::FieldType::T19> beams;     ///< Ordered beams entries; size() and at() expose the
                                             ///< checked sequence without copying.
    Field<detail::FieldType::T20> joints;    ///< Ordered joints entries; size() and at() expose the
                                             ///< checked sequence without copying.
    Field<detail::FieldType::T21> supports;  ///< Ordered supports entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T22> masses;    ///< Ordered masses entries; size() and at() expose the
                                             ///< checked sequence without copying.
    Field<detail::FieldType::T23> aero_surfaces; ///< Ordered aero surfaces entries; size() and at()
                                                 ///< expose the checked sequence without copying.
};
/// Checked Provenance payload; all fields retain one immutable case owner.
struct Provenance {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Provenance(Field<detail::FieldType::T24> description_value,
               Field<detail::FieldType::T25> references_value,
               Field<detail::FieldType::T25> assumptions_value) noexcept
        : description(description_value),
          references(references_value),
          assumptions(assumptions_value) {}
    Field<detail::FieldType::T24>
        description; ///< Description as decoded UTF-8 text; no normalization or implicit
                     ///< replacement is performed.
    Field<detail::FieldType::T25> references;  ///< Ordered references entries; size() and at()
                                               ///< expose the checked sequence without copying.
    Field<detail::FieldType::T25> assumptions; ///< Ordered assumptions entries; size() and at()
                                               ///< expose the checked sequence without copying.
};
/// Checked Frame payload; all fields retain one immutable case owner.
struct Frame {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Frame(Field<detail::FieldType::T26> origin_value,
          Field<detail::FieldType::T27> rotation_value) noexcept
        : origin(origin_value), rotation(rotation_value) {}
    Field<detail::FieldType::T26>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T27> rotation; ///< Rotation matrix with ordered rows; frame regularity
                                            ///< is checked during preparation.
};
/// Checked Target payload; all fields retain one immutable case owner.
struct Target {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Target(Field<detail::FieldType::T28> kind_value,
           Field<detail::FieldType::T29> id_value) noexcept
        : kind(kind_value), id(id_value) {}
    Field<detail::FieldType::T28> kind; ///< Kind; choose one of shell,beam,rigid.
    Field<detail::FieldType::T29> id;   ///< Identifier of the shell, beam or rigid mass selected by
                                        ///< the accompanying target kind.
};
/// Checked MaterialPoint payload; all fields retain one immutable case owner.
struct MaterialPoint {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    MaterialPoint(Field<detail::FieldType::T30> target_value,
                  Field<detail::FieldType::T31> coordinates_value,
                  Field<detail::FieldType::T32> side_value) noexcept
        : target(target_value), coordinates(coordinates_value), side(side_value) {}
    Field<detail::FieldType::T30>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T31>
        coordinates; ///< Reference material coordinates; component dimensions depend on the
                     ///< selected shell, beam or rigid target.
    Field<detail::FieldType::T32> side; ///< Side; choose one of interior,lower,upper.
};
/// Checked Region payload; all fields retain one immutable case owner.
struct Region {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Region(Field<detail::FieldType::T30> target_value, Field<detail::FieldType::T33> kind_value,
           Field<detail::FieldType::T34> geometry_value,
           Field<detail::FieldType::T35> boundary_value,
           Field<detail::FieldType::T32> side_value) noexcept
        : target(target_value),
          kind(kind_value),
          geometry(geometry_value),
          boundary(boundary_value),
          side(side_value) {}
    Field<detail::FieldType::T30>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T33> kind; ///< Kind; choose one of volume,face,edge,point.
    Field<detail::FieldType::T34>
        geometry; ///< Geometry referring to a checked Geometry definition in this case.
    Field<detail::FieldType::T35> boundary; ///< Zero-based index of boundary.
    Field<detail::FieldType::T32> side;     ///< Side; choose one of interior,lower,upper.
};
/// Checked PolyTerm<d,k> payload; all fields retain one immutable case owner.
struct PolyTerm {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PolyTerm(Field<detail::FieldType::T36> powers_value,
             Field<detail::FieldType::T37> coefficient_value) noexcept
        : powers(powers_value), coefficient(coefficient_value) {}
    Field<detail::FieldType::T36> powers; ///< Nonnegative monomial exponents, in coordinate order.
    Field<detail::FieldType::T37> coefficient; ///< Exact d quantity for coefficient, in declared
                                               ///< case units; number() retains its dimensions.
};
/// Checked ReferenceField<d,n> constant payload; all fields retain one immutable case owner.
struct ReferenceFieldConstant {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ReferenceFieldConstant(Field<detail::FieldType::T38> kind_value,
                           Field<detail::FieldType::T39> value_value) noexcept
        : kind(kind_value), value(value_value) {}
    Field<detail::FieldType::T38> kind;  ///< Kind; choose one of constant.
    Field<detail::FieldType::T39> value; ///< Ordered value components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
};
/// Checked ReferenceField<d,n> polynomial payload; all fields retain one immutable case owner.
struct ReferenceFieldPolynomial {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ReferenceFieldPolynomial(Field<detail::FieldType::T40> kind_value,
                             Field<detail::FieldType::T41> coordinates_value,
                             Field<detail::FieldType::T31> origin_value,
                             Field<detail::FieldType::T31> lengths_value,
                             Field<detail::FieldType::T42> components_value,
                             Field<detail::FieldType::T39> error_value) noexcept
        : kind(kind_value),
          coordinates(coordinates_value),
          origin(origin_value),
          lengths(lengths_value),
          components(components_value),
          error(error_value) {}
    Field<detail::FieldType::T40> kind; ///< Kind; choose one of polynomial.
    Field<detail::FieldType::T41>
        coordinates; ///< Coordinate system used by this field; its choices are closed.
    Field<detail::FieldType::T31>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T31>
        lengths; ///< Positive coordinate normalization lengths, ordered x, y, z.
    Field<detail::FieldType::T42>
        components; ///< Ordered components components; each numeric element exposes its physical
                    ///< dimensions and exact input value.
    Field<detail::FieldType::T39> error; ///< Ordered error components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
};
/// Checked TimeSeries<n> constant payload; all fields retain one immutable case owner.
struct TimeSeriesConstant {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    TimeSeriesConstant(Field<detail::FieldType::T38> kind_value,
                       Field<detail::FieldType::T43> value_value) noexcept
        : kind(kind_value), value(value_value) {}
    Field<detail::FieldType::T38> kind;  ///< Kind; choose one of constant.
    Field<detail::FieldType::T43> value; ///< Ordered value components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
};
/// Checked TimeSeries<n> polynomial payload; all fields retain one immutable case owner.
struct TimeSeriesPolynomial {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    TimeSeriesPolynomial(Field<detail::FieldType::T40> kind_value,
                         Field<detail::FieldType::T44> origin_value,
                         Field<detail::FieldType::T44> scale_value,
                         Field<detail::FieldType::T45> coefficients_value,
                         Field<detail::FieldType::T46> domain_value) noexcept
        : kind(kind_value),
          origin(origin_value),
          scale(scale_value),
          coefficients(coefficients_value),
          domain(domain_value) {}
    Field<detail::FieldType::T40> kind; ///< Kind; choose one of polynomial.
    Field<detail::FieldType::T44>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T44> scale; ///< Exact time quantity for scale, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T45> coefficients; ///< Ordered coefficients entries; size() and at()
                                                ///< expose the checked sequence without copying.
    Field<detail::FieldType::T46>
        domain; ///< Ordered domain components; each numeric element exposes its physical dimensions
                ///< and exact input value.
};
/// Checked TimeSeries<n> pchip payload; all fields retain one immutable case owner.
struct TimeSeriesPchip {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    TimeSeriesPchip(Field<detail::FieldType::T47> kind_value,
                    Field<detail::FieldType::T48> times_value,
                    Field<detail::FieldType::T45> values_value,
                    Field<detail::FieldType::T49> outside_value) noexcept
        : kind(kind_value), times(times_value), values(values_value), outside(outside_value) {}
    Field<detail::FieldType::T47> kind;   ///< Kind; choose one of pchip.
    Field<detail::FieldType::T48> times;  ///< Ordered times entries; size() and at() expose the
                                          ///< checked sequence without copying.
    Field<detail::FieldType::T45> values; ///< Ordered values entries; size() and at() expose the
                                          ///< checked sequence without copying.
    Field<detail::FieldType::T49>
        outside; ///< Outside payload; record() selects its checked concrete structure.
};
/// Checked Outside reject/hold payload; all fields retain one immutable case owner.
struct OutsideRejectHold {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    OutsideRejectHold(Field<detail::FieldType::T50> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T50> kind; ///< Kind; choose one of reject,hold.
};
/// Checked Outside periodic payload; all fields retain one immutable case owner.
struct OutsidePeriodic {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    OutsidePeriodic(Field<detail::FieldType::T51> kind_value,
                    Field<detail::FieldType::T44> period_value) noexcept
        : kind(kind_value), period(period_value) {}
    Field<detail::FieldType::T51> kind;   ///< Kind; choose one of periodic.
    Field<detail::FieldType::T44> period; ///< Exact time quantity for period, in declared case
                                          ///< units; number() retains its dimensions.
};
/// Checked Motion payload; all fields retain one immutable case owner.
struct Motion {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Motion(Field<detail::FieldType::T52> origin_value,
           Field<detail::FieldType::T53> translation_value,
           Field<detail::FieldType::T53> rotation_vector_value,
           Field<detail::FieldType::T48> events_value) noexcept
        : origin(origin_value),
          translation(translation_value),
          rotation_vector(rotation_vector_value),
          events(events_value) {}
    Field<detail::FieldType::T52>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T53> translation; ///< Translation as a function of time in declared
                                               ///< case units; evaluation belongs to preparation.
    Field<detail::FieldType::T53>
        rotation_vector;                  ///< Body rotation vector field, in the case angle unit.
    Field<detail::FieldType::T48> events; ///< Ordered events entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked Axis projected payload; all fields retain one immutable case owner.
struct AxisProjected {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AxisProjected(Field<detail::FieldType::T54> kind_value,
                  Field<detail::FieldType::T55> direction_value,
                  Field<detail::FieldType::T56> angle_value) noexcept
        : kind(kind_value), direction(direction_value), angle(angle_value) {}
    Field<detail::FieldType::T54> kind; ///< Kind; choose one of projected.
    Field<detail::FieldType::T55>
        direction; ///< Direction over reference coordinates; evaluation and physical-domain
                   ///< qualification belong to preparation.
    Field<detail::FieldType::T56> angle; ///< Angle over reference coordinates; evaluation and
                                         ///< physical-domain qualification belong to preparation.
};
/// Checked Axis explicit payload; all fields retain one immutable case owner.
struct AxisExplicit {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AxisExplicit(Field<detail::FieldType::T57> kind_value,
                 Field<detail::FieldType::T58> columns_value) noexcept
        : kind(kind_value), columns(columns_value) {}
    Field<detail::FieldType::T57> kind;    ///< Kind; choose one of explicit.
    Field<detail::FieldType::T58> columns; ///< Columns over reference coordinates; evaluation and
                                           ///< physical-domain qualification belong to preparation.
};
/// Checked Geometry pg payload; all fields retain one immutable case owner.
struct GeometryPg {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    GeometryPg(Field<detail::FieldType::T59> kind_value, Field<detail::FieldType::T60> id_value,
               Field<detail::FieldType::T52> frame_value,
               Field<detail::FieldType::T61> reference_length_value,
               Field<detail::FieldType::T62> x0_value, Field<detail::FieldType::T62> width_value,
               Field<detail::FieldType::T62> skew_value, Field<detail::FieldType::T63> degree_value,
               Field<detail::FieldType::T64> coefficients_value,
               Field<detail::FieldType::T65> domain_value,
               Field<detail::FieldType::T66> thickness_value, Field<detail::FieldType::T67> n_value,
               Field<detail::FieldType::T68> convention_value,
               Field<detail::FieldType::T69> boundaries_value,
               Field<detail::FieldType::T70> continuity_value) noexcept
        : kind(kind_value),
          id(id_value),
          frame(frame_value),
          reference_length(reference_length_value),
          x0(x0_value),
          width(width_value),
          skew(skew_value),
          degree(degree_value),
          coefficients(coefficients_value),
          domain(domain_value),
          thickness(thickness_value),
          n(n_value),
          convention(convention_value),
          boundaries(boundaries_value),
          continuity(continuity_value) {}
    Field<detail::FieldType::T59> kind; ///< Kind; choose one of pg.
    Field<detail::FieldType::T60> id;   ///< Id in the Geometry identifier namespace; equal IDs in
                                        ///< other namespaces remain distinct.
    Field<detail::FieldType::T52>
        frame; ///< Frame payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T61> reference_length; ///< Positive physical length used to
                                                    ///< dimensionalize the reference geometry.
    Field<detail::FieldType::T62> x0;               ///< Exact dimensionless chart origin.
    Field<detail::FieldType::T62> width;            ///< Exact dimensionless chart width.
    Field<detail::FieldType::T62> skew;             ///< Exact dimensionless chart skew.
    Field<detail::FieldType::T63>
        degree; ///< Ordered degree components; each numeric element exposes its physical dimensions
                ///< and exact input value.
    Field<detail::FieldType::T64> coefficients; ///< Ordered coefficients entries; size() and at()
                                                ///< expose the checked sequence without copying.
    Field<detail::FieldType::T65> domain; ///< Ordered domain entries; size() and at() expose the
                                          ///< checked sequence without copying.
    Field<detail::FieldType::T66>
        thickness; ///< Thickness payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T67> n;          ///< Nonnegative count of n.
    Field<detail::FieldType::T68> convention; ///< Convention; choose one of pg-graph-v1.
    Field<detail::FieldType::T69> boundaries; ///< Ordered boundaries entries; size() and at()
                                              ///< expose the checked sequence without copying.
    Field<detail::FieldType::T70> continuity; ///< Ordered continuity entries; size() and at()
                                              ///< expose the checked sequence without copying.
};
/// Checked PgTerm payload; all fields retain one immutable case owner.
struct PgTerm {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PgTerm(Field<detail::FieldType::T63> powers_value,
           Field<detail::FieldType::T62> coefficient_value) noexcept
        : powers(powers_value), coefficient(coefficient_value) {}
    Field<detail::FieldType::T63> powers; ///< Nonnegative monomial exponents, in coordinate order.
    Field<detail::FieldType::T62>
        coefficient; ///< Exact rational coefficient; numerator and positive denominator remain
                     ///< available separately.
};
/// Checked Thickness payload; all fields retain one immutable case owner.
struct Thickness {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Thickness(Field<detail::FieldType::T65> vertices_value,
              Field<detail::FieldType::T71> triangles_value,
              Field<detail::FieldType::T72> values_value) noexcept
        : vertices(vertices_value), triangles(triangles_value), values(values_value) {}
    Field<detail::FieldType::T65> vertices;  ///< Ordered vertices entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T71> triangles; ///< Ordered triangles entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T72> values;    ///< Ordered values entries; size() and at() expose the
                                             ///< checked sequence without copying.
};
/// Checked Boundary payload; all fields retain one immutable case owner.
struct Boundary {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Boundary(Field<detail::FieldType::T24> name_value,
             Field<detail::FieldType::T65> vertices_value) noexcept
        : name(name_value), vertices(vertices_value) {}
    Field<detail::FieldType::T24> name;     ///< Name as decoded UTF-8 text; no normalization or
                                            ///< implicit replacement is performed.
    Field<detail::FieldType::T65> vertices; ///< Ordered vertices entries; size() and at() expose
                                            ///< the checked sequence without copying.
};
/// Checked Continuity payload; all fields retain one immutable case owner.
struct Continuity {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Continuity(Field<detail::FieldType::T34> other_value, Field<detail::FieldType::T73> kind_value,
               Field<detail::FieldType::T35> own_boundary_value,
               Field<detail::FieldType::T35> other_boundary_value,
               Field<detail::FieldType::T74> orientation_value) noexcept
        : other(other_value),
          kind(kind_value),
          own_boundary(own_boundary_value),
          other_boundary(other_boundary_value),
          orientation(orientation_value) {}
    Field<detail::FieldType::T34>
        other; ///< Other referring to a checked Geometry definition in this case.
    Field<detail::FieldType::T73> kind;           ///< Kind; choose one of smooth,fold.
    Field<detail::FieldType::T35> own_boundary;   ///< Zero-based index of own boundary.
    Field<detail::FieldType::T35> other_boundary; ///< Zero-based index of other boundary.
    Field<detail::FieldType::T74> orientation;    ///< Orientation; choose one of same,reversed.
};
/// Checked Geometry surface_graph payload; all fields retain one immutable case owner.
struct GeometrySurfaceGraph {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    GeometrySurfaceGraph(
        Field<detail::FieldType::T75> kind_value, Field<detail::FieldType::T60> id_value,
        Field<detail::FieldType::T52> frame_value,
        Field<detail::FieldType::T61> reference_length_value,
        Field<detail::FieldType::T62> x0_value, Field<detail::FieldType::T62> width_value,
        Field<detail::FieldType::T62> skew_value, Field<detail::FieldType::T63> degree_value,
        Field<detail::FieldType::T64> coefficients_value,
        Field<detail::FieldType::T65> domain_value, Field<detail::FieldType::T69> boundaries_value,
        Field<detail::FieldType::T70> continuity_value) noexcept
        : kind(kind_value),
          id(id_value),
          frame(frame_value),
          reference_length(reference_length_value),
          x0(x0_value),
          width(width_value),
          skew(skew_value),
          degree(degree_value),
          coefficients(coefficients_value),
          domain(domain_value),
          boundaries(boundaries_value),
          continuity(continuity_value) {}
    Field<detail::FieldType::T75> kind; ///< Kind; choose one of surface_graph.
    Field<detail::FieldType::T60> id;   ///< Id in the Geometry identifier namespace; equal IDs in
                                        ///< other namespaces remain distinct.
    Field<detail::FieldType::T52>
        frame; ///< Frame payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T61> reference_length; ///< Positive physical length used to
                                                    ///< dimensionalize the reference geometry.
    Field<detail::FieldType::T62> x0;               ///< Exact dimensionless chart origin.
    Field<detail::FieldType::T62> width;            ///< Exact dimensionless chart width.
    Field<detail::FieldType::T62> skew;             ///< Exact dimensionless chart skew.
    Field<detail::FieldType::T63>
        degree; ///< Ordered degree components; each numeric element exposes its physical dimensions
                ///< and exact input value.
    Field<detail::FieldType::T64> coefficients; ///< Ordered coefficients entries; size() and at()
                                                ///< expose the checked sequence without copying.
    Field<detail::FieldType::T65> domain; ///< Ordered domain entries; size() and at() expose the
                                          ///< checked sequence without copying.
    Field<detail::FieldType::T69> boundaries; ///< Ordered boundaries entries; size() and at()
                                              ///< expose the checked sequence without copying.
    Field<detail::FieldType::T70> continuity; ///< Ordered continuity entries; size() and at()
                                              ///< expose the checked sequence without copying.
};
/// Checked Geometry nastran payload; all fields retain one immutable case owner.
struct GeometryNastran {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    GeometryNastran(Field<detail::FieldType::T76> kind_value,
                    Field<detail::FieldType::T60> id_value,
                    Field<detail::FieldType::T77> deck_value,
                    Field<detail::FieldType::T78> sidecar_value,
                    Field<detail::FieldType::T79> dialect_value,
                    Field<detail::FieldType::T80> mode_value,
                    Field<detail::FieldType::T25> excluded_cards_value) noexcept
        : kind(kind_value),
          id(id_value),
          deck(deck_value),
          sidecar(sidecar_value),
          dialect(dialect_value),
          mode(mode_value),
          excluded_cards(excluded_cards_value) {}
    Field<detail::FieldType::T76> kind; ///< Kind; choose one of nastran.
    Field<detail::FieldType::T60> id;   ///< Id in the Geometry identifier namespace; equal IDs in
                                        ///< other namespaces remain distinct.
    Field<detail::FieldType::T77>
        deck; ///< Deck source reference, including owned source ID, path, digest and expected role.
    Field<detail::FieldType::T78> sidecar; ///< Sidecar source reference, including owned source ID,
                                           ///< path, digest and expected role.
    Field<detail::FieldType::T79> dialect; ///< Dialect; choose one of msc-2022.1.
    Field<detail::FieldType::T80> mode;    ///< Mode; choose one of structural-extraction.
    Field<detail::FieldType::T25>
        excluded_cards; ///< Ordered excluded cards entries; size() and at() expose the checked
                        ///< sequence without copying.
};
/// Checked Geometry gmsh payload; all fields retain one immutable case owner.
struct GeometryGmsh {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    GeometryGmsh(Field<detail::FieldType::T81> kind_value, Field<detail::FieldType::T60> id_value,
                 Field<detail::FieldType::T82> mesh_value,
                 Field<detail::FieldType::T83> sidecar_value,
                 Field<detail::FieldType::T84> reference_value) noexcept
        : kind(kind_value),
          id(id_value),
          mesh(mesh_value),
          sidecar(sidecar_value),
          reference(reference_value) {}
    Field<detail::FieldType::T81> kind; ///< Kind; choose one of gmsh.
    Field<detail::FieldType::T60> id;   ///< Id in the Geometry identifier namespace; equal IDs in
                                        ///< other namespaces remain distinct.
    Field<detail::FieldType::T82>
        mesh; ///< Mesh source reference, including owned source ID, path, digest and expected role.
    Field<detail::FieldType::T83> sidecar; ///< Sidecar source reference, including owned source ID,
                                           ///< path, digest and expected role.
    Field<detail::FieldType::T84>
        reference; ///< Reference payload; record() selects its checked concrete structure.
};
/// Checked ReferenceGeometry supplied payload; all fields retain one immutable case owner.
struct ReferenceGeometrySupplied {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ReferenceGeometrySupplied(Field<detail::FieldType::T85> kind_value,
                              Field<detail::FieldType::T86> file_value) noexcept
        : kind(kind_value), file(file_value) {}
    Field<detail::FieldType::T85> kind; ///< Kind; choose one of supplied.
    Field<detail::FieldType::T86>
        file; ///< File source reference, including owned source ID, path, digest and expected role.
};
/// Checked ReferenceGeometry fit payload; all fields retain one immutable case owner.
struct ReferenceGeometryFit {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ReferenceGeometryFit(Field<detail::FieldType::T87> kind_value,
                         Field<detail::FieldType::T67> degree_value,
                         Field<detail::FieldType::T61> distance_budget_value,
                         Field<detail::FieldType::T88> normal_budget_value,
                         Field<detail::FieldType::T89> conversion_observations_value) noexcept
        : kind(kind_value),
          degree(degree_value),
          distance_budget(distance_budget_value),
          normal_budget(normal_budget_value),
          conversion_observations(conversion_observations_value) {}
    Field<detail::FieldType::T87> kind;   ///< Kind; choose one of fit.
    Field<detail::FieldType::T67> degree; ///< Nonnegative count of degree.
    Field<detail::FieldType::T61>
        distance_budget; ///< Exact length quantity for distance budget, in declared case units;
                         ///< number() retains its dimensions.
    Field<detail::FieldType::T88>
        normal_budget; ///< Exact dimensionless quantity for normal budget, in declared case units;
                       ///< number() retains its dimensions.
    Field<detail::FieldType::T89>
        conversion_observations; ///< Ordered conversion observations entries; size() and at()
                                 ///< expose the checked sequence without copying.
};
/// Checked Material isotropic payload; all fields retain one immutable case owner.
struct MaterialIsotropic {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    MaterialIsotropic(Field<detail::FieldType::T90> kind_value,
                      Field<detail::FieldType::T91> id_value,
                      Field<detail::FieldType::T92> young_value,
                      Field<detail::FieldType::T88> poisson_value,
                      Field<detail::FieldType::T93> density_value,
                      Field<detail::FieldType::T11> provenance_value,
                      Field<detail::FieldType::T94> applicability_value) noexcept
        : kind(kind_value),
          id(id_value),
          young(young_value),
          poisson(poisson_value),
          density(density_value),
          provenance(provenance_value),
          applicability(applicability_value) {}
    Field<detail::FieldType::T90> kind;  ///< Kind; choose one of isotropic.
    Field<detail::FieldType::T91> id;    ///< Id in the Material identifier namespace; equal IDs in
                                         ///< other namespaces remain distinct.
    Field<detail::FieldType::T92> young; ///< Exact pressure quantity for young, in declared case
                                         ///< units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        poisson; ///< Exact dimensionless quantity for poisson, in declared case units; number()
                 ///< retains its dimensions.
    Field<detail::FieldType::T93> density; ///< Exact M/L³ quantity for density, in declared case
                                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T11>
        provenance; ///< Provenance payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T94>
        applicability; ///< Applicability payload; record() selects its checked concrete structure.
};
/// Checked Material orthotropic payload; all fields retain one immutable case owner.
struct MaterialOrthotropic {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    MaterialOrthotropic(Field<detail::FieldType::T95> kind_value,
                        Field<detail::FieldType::T91> id_value,
                        Field<detail::FieldType::T96> young_value,
                        Field<detail::FieldType::T96> shear_value,
                        Field<detail::FieldType::T97> poisson_value,
                        Field<detail::FieldType::T93> density_value,
                        Field<detail::FieldType::T11> provenance_value,
                        Field<detail::FieldType::T94> applicability_value) noexcept
        : kind(kind_value),
          id(id_value),
          young(young_value),
          shear(shear_value),
          poisson(poisson_value),
          density(density_value),
          provenance(provenance_value),
          applicability(applicability_value) {}
    Field<detail::FieldType::T95> kind;  ///< Kind; choose one of orthotropic3d.
    Field<detail::FieldType::T91> id;    ///< Id in the Material identifier namespace; equal IDs in
                                         ///< other namespaces remain distinct.
    Field<detail::FieldType::T96> young; ///< Ordered young components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
    Field<detail::FieldType::T96> shear; ///< Ordered shear components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
    Field<detail::FieldType::T97>
        poisson; ///< Ordered poisson components; each numeric element exposes its physical
                 ///< dimensions and exact input value.
    Field<detail::FieldType::T93> density; ///< Exact M/L³ quantity for density, in declared case
                                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T11>
        provenance; ///< Provenance payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T94>
        applicability; ///< Applicability payload; record() selects its checked concrete structure.
};
/// Checked Material effective payload; all fields retain one immutable case owner.
struct MaterialEffective {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    MaterialEffective(Field<detail::FieldType::T98> kind_value,
                      Field<detail::FieldType::T91> id_value,
                      Field<detail::FieldType::T99> stiffness_value,
                      Field<detail::FieldType::T100> components_value,
                      Field<detail::FieldType::T93> density_value,
                      Field<detail::FieldType::T11> provenance_value,
                      Field<detail::FieldType::T94> applicability_value) noexcept
        : kind(kind_value),
          id(id_value),
          stiffness(stiffness_value),
          components(components_value),
          density(density_value),
          provenance(provenance_value),
          applicability(applicability_value) {}
    Field<detail::FieldType::T98> kind; ///< Kind; choose one of effective3d.
    Field<detail::FieldType::T91> id;   ///< Id in the Material identifier namespace; equal IDs in
                                        ///< other namespaces remain distinct.
    Field<detail::FieldType::T99>
        stiffness; ///< Stiffness as ordered matrix rows; component dimensions follow the
                   ///< corresponding input and output quantities.
    Field<detail::FieldType::T100> components; ///< Components; choose one of mandel,engineering.
    Field<detail::FieldType::T93> density; ///< Exact M/L³ quantity for density, in declared case
                                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T11>
        provenance; ///< Provenance payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T94>
        applicability; ///< Applicability payload; record() selects its checked concrete structure.
};
/// Checked Applicability unevaluated payload; all fields retain one immutable case owner.
struct ApplicabilityUnevaluated {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ApplicabilityUnevaluated(Field<detail::FieldType::T101> kind_value) noexcept
        : kind(kind_value) {}
    Field<detail::FieldType::T101> kind; ///< Kind; choose one of unevaluated.
};
/// Checked Applicability declared payload; all fields retain one immutable case owner.
struct ApplicabilityDeclared {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ApplicabilityDeclared(Field<detail::FieldType::T102> kind_value,
                          Field<detail::FieldType::T103> strain_lower_value,
                          Field<detail::FieldType::T103> strain_upper_value,
                          Field<detail::FieldType::T24> reference_value) noexcept
        : kind(kind_value),
          strain_lower(strain_lower_value),
          strain_upper(strain_upper_value),
          reference(reference_value) {}
    Field<detail::FieldType::T102> kind; ///< Kind; choose one of declared.
    Field<detail::FieldType::T103>
        strain_lower; ///< Ordered strain lower components; each numeric element exposes its
                      ///< physical dimensions and exact input value.
    Field<detail::FieldType::T103>
        strain_upper; ///< Ordered strain upper components; each numeric element exposes its
                      ///< physical dimensions and exact input value.
    Field<detail::FieldType::T24> reference; ///< Reference as decoded UTF-8 text; no normalization
                                             ///< or implicit replacement is performed.
};
/// Checked Layer payload; all fields retain one immutable case owner.
struct Layer {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Layer(Field<detail::FieldType::T104> id_value, Field<detail::FieldType::T105> interval_value,
          Field<detail::FieldType::T106> material_value, Field<detail::FieldType::T107> axis_value,
          Field<detail::FieldType::T108> initial_green_value,
          Field<detail::FieldType::T109> viscosity_value) noexcept
        : id(id_value),
          interval(interval_value),
          material(material_value),
          axis(axis_value),
          initial_green(initial_green_value),
          viscosity(viscosity_value) {}
    Field<detail::FieldType::T104> id; ///< Id in the Layer identifier namespace; equal IDs in other
                                       ///< namespaces remain distinct.
    Field<detail::FieldType::T105>
        interval; ///< Ordered exact layer endpoints through the normalized shell thickness.
    Field<detail::FieldType::T106>
        material; ///< Material referring to a checked Material definition in this case.
    Field<detail::FieldType::T107>
        axis; ///< Axis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T108>
        initial_green; ///< Initial green over reference coordinates; evaluation and physical-domain
                       ///< qualification belong to preparation.
    Field<detail::FieldType::T109>
        viscosity; ///< Viscosity as ordered matrix rows; component dimensions follow the
                   ///< corresponding input and output quantities.
};
/// Checked Shell payload; all fields retain one immutable case owner.
struct Shell {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Shell(Field<detail::FieldType::T110> id_value, Field<detail::FieldType::T34> geometry_value,
          Field<detail::FieldType::T111> layers_value,
          Field<detail::FieldType::T100> damping_components_value,
          Field<detail::FieldType::T112> loads_value) noexcept
        : id(id_value),
          geometry(geometry_value),
          layers(layers_value),
          damping_components(damping_components_value),
          loads(loads_value) {}
    Field<detail::FieldType::T110> id; ///< Id in the Shell identifier namespace; equal IDs in other
                                       ///< namespaces remain distinct.
    Field<detail::FieldType::T34>
        geometry; ///< Geometry referring to a checked Geometry definition in this case.
    Field<detail::FieldType::T111> layers; ///< Ordered layers entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T100>
        damping_components;               ///< Damping components; choose one of mandel,engineering.
    Field<detail::FieldType::T112> loads; ///< Ordered loads entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked Beam payload; all fields retain one immutable case owner.
struct Beam {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Beam(Field<detail::FieldType::T113> id_value, Field<detail::FieldType::T114> nodes_value,
         Field<detail::FieldType::T115> coordinate_value,
         Field<detail::FieldType::T116> section_value,
         Field<detail::FieldType::T103> initial_strain_value,
         Field<detail::FieldType::T117> residual_energy_value,
         Field<detail::FieldType::T118> viscosity_value,
         Field<detail::FieldType::T112> loads_value) noexcept
        : id(id_value),
          nodes(nodes_value),
          coordinate(coordinate_value),
          section(section_value),
          initial_strain(initial_strain_value),
          residual_energy(residual_energy_value),
          viscosity(viscosity_value),
          loads(loads_value) {}
    Field<detail::FieldType::T113>
        id; ///< Id in the Beam identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T114>
        nodes; ///< Ordered nodes components; each numeric element exposes its physical dimensions
               ///< and exact input value.
    Field<detail::FieldType::T115>
        coordinate; ///< Coordinate payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T116>
        section; ///< Section payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T103>
        initial_strain; ///< Ordered initial strain components; each numeric element exposes its
                        ///< physical dimensions and exact input value.
    Field<detail::FieldType::T117>
        residual_energy; ///< Exact energy/L quantity for residual energy, in declared case units;
                         ///< number() retains its dimensions.
    Field<detail::FieldType::T118>
        viscosity; ///< Viscosity as ordered matrix rows; component dimensions follow the
                   ///< corresponding input and output quantities.
    Field<detail::FieldType::T112> loads; ///< Ordered loads entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked BeamCoordinate arclength payload; all fields retain one immutable case owner.
struct BeamCoordinateArclength {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BeamCoordinateArclength(Field<detail::FieldType::T119> kind_value,
                            Field<detail::FieldType::T61> length_value) noexcept
        : kind(kind_value), length(length_value) {}
    Field<detail::FieldType::T119> kind;  ///< Kind; choose one of arclength.
    Field<detail::FieldType::T61> length; ///< Exact length quantity for length, in declared case
                                          ///< units; number() retains its dimensions.
};
/// Checked BeamCoordinate scaled payload; all fields retain one immutable case owner.
struct BeamCoordinateScaled {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BeamCoordinateScaled(Field<detail::FieldType::T120> kind_value,
                         Field<detail::FieldType::T61> length_value,
                         Field<detail::FieldType::T88> arc_scale_value) noexcept
        : kind(kind_value), length(length_value), arc_scale(arc_scale_value) {}
    Field<detail::FieldType::T120> kind;  ///< Kind; choose one of scaled.
    Field<detail::FieldType::T61> length; ///< Exact length quantity for length, in declared case
                                          ///< units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        arc_scale; ///< Exact dimensionless quantity for arc scale, in declared case units; number()
                   ///< retains its dimensions.
};
/// Checked Section matrix payload; all fields retain one immutable case owner.
struct SectionMatrix {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SectionMatrix(Field<detail::FieldType::T121> kind_value,
                  Field<detail::FieldType::T118> stiffness_value,
                  Field<detail::FieldType::T122> line_density_value,
                  Field<detail::FieldType::T26> center_value,
                  Field<detail::FieldType::T123> inertia_value,
                  Field<detail::FieldType::T124> envelope_value,
                  Field<detail::FieldType::T125> warping_value,
                  Field<detail::FieldType::T126> recovery_value,
                  Field<detail::FieldType::T24> source_value) noexcept
        : kind(kind_value),
          stiffness(stiffness_value),
          line_density(line_density_value),
          center(center_value),
          inertia(inertia_value),
          envelope(envelope_value),
          warping(warping_value),
          recovery(recovery_value),
          source(source_value) {}
    Field<detail::FieldType::T121> kind; ///< Kind; choose one of matrix.
    Field<detail::FieldType::T118>
        stiffness; ///< Stiffness as ordered matrix rows; component dimensions follow the
                   ///< corresponding input and output quantities.
    Field<detail::FieldType::T122>
        line_density; ///< Exact M/L quantity for line density, in declared case units; number()
                      ///< retains its dimensions.
    Field<detail::FieldType::T26>
        center; ///< Ordered center components; each numeric element exposes its physical dimensions
                ///< and exact input value.
    Field<detail::FieldType::T123>
        inertia; ///< Inertia as ordered matrix rows; component dimensions follow the corresponding
                 ///< input and output quantities.
    Field<detail::FieldType::T124>
        envelope; ///< Ordered envelope components; each numeric element exposes its physical
                  ///< dimensions and exact input value.
    Field<detail::FieldType::T125> warping;  ///< Warping; choose one of free.
    Field<detail::FieldType::T126> recovery; ///< Recovery; choose one of unavailable.
    Field<detail::FieldType::T24> source;    ///< Source as decoded UTF-8 text; no normalization or
                                             ///< implicit replacement is performed.
};
/// Checked Section circle payload; all fields retain one immutable case owner.
struct SectionCircle {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SectionCircle(Field<detail::FieldType::T127> kind_value,
                  Field<detail::FieldType::T61> radius_value,
                  Field<detail::FieldType::T92> young_value,
                  Field<detail::FieldType::T88> poisson_value,
                  Field<detail::FieldType::T93> density_value,
                  Field<detail::FieldType::T128> effective_shear_value,
                  Field<detail::FieldType::T125> warping_value) noexcept
        : kind(kind_value),
          radius(radius_value),
          young(young_value),
          poisson(poisson_value),
          density(density_value),
          effective_shear(effective_shear_value),
          warping(warping_value) {}
    Field<detail::FieldType::T127> kind;  ///< Kind; choose one of circle.
    Field<detail::FieldType::T61> radius; ///< Exact length quantity for radius, in declared case
                                          ///< units; number() retains its dimensions.
    Field<detail::FieldType::T92> young;  ///< Exact pressure quantity for young, in declared case
                                          ///< units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        poisson; ///< Exact dimensionless quantity for poisson, in declared case units; number()
                 ///< retains its dimensions.
    Field<detail::FieldType::T93> density; ///< Exact M/L³ quantity for density, in declared case
                                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T128>
        effective_shear; ///< Ordered effective shear components; each numeric element exposes its
                         ///< physical dimensions and exact input value.
    Field<detail::FieldType::T125> warping; ///< Warping; choose one of free.
};
/// Checked Section rectangle payload; all fields retain one immutable case owner.
struct SectionRectangle {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SectionRectangle(Field<detail::FieldType::T129> kind_value,
                     Field<detail::FieldType::T124> sides_value,
                     Field<detail::FieldType::T92> young_value,
                     Field<detail::FieldType::T88> poisson_value,
                     Field<detail::FieldType::T93> density_value,
                     Field<detail::FieldType::T128> effective_shear_value,
                     Field<detail::FieldType::T130> torsion_value,
                     Field<detail::FieldType::T125> warping_value) noexcept
        : kind(kind_value),
          sides(sides_value),
          young(young_value),
          poisson(poisson_value),
          density(density_value),
          effective_shear(effective_shear_value),
          torsion(torsion_value),
          warping(warping_value) {}
    Field<detail::FieldType::T129> kind; ///< Kind; choose one of rectangle.
    Field<detail::FieldType::T124>
        sides; ///< Ordered sides components; each numeric element exposes its physical dimensions
               ///< and exact input value.
    Field<detail::FieldType::T92> young; ///< Exact pressure quantity for young, in declared case
                                         ///< units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        poisson; ///< Exact dimensionless quantity for poisson, in declared case units; number()
                 ///< retains its dimensions.
    Field<detail::FieldType::T93> density; ///< Exact M/L³ quantity for density, in declared case
                                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T128>
        effective_shear; ///< Ordered effective shear components; each numeric element exposes its
                         ///< physical dimensions and exact input value.
    Field<detail::FieldType::T130> torsion; ///< Exact L⁴ quantity for torsion, in declared case
                                            ///< units; number() retains its dimensions.
    Field<detail::FieldType::T125> warping; ///< Warping; choose one of free.
};
/// Checked Joint adhesive payload; all fields retain one immutable case owner.
struct JointAdhesive {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    JointAdhesive(Field<detail::FieldType::T131> kind_value,
                  Field<detail::FieldType::T132> id_value,
                  Field<detail::FieldType::T133> first_value,
                  Field<detail::FieldType::T133> second_value,
                  Field<detail::FieldType::T134> mapping_value,
                  Field<detail::FieldType::T135> split_at_attachment_value) noexcept
        : kind(kind_value),
          id(id_value),
          first(first_value),
          second(second_value),
          mapping(mapping_value),
          split_at_attachment(split_at_attachment_value) {}
    Field<detail::FieldType::T131> kind; ///< Kind; choose one of adhesive.
    Field<detail::FieldType::T132> id; ///< Id in the Joint identifier namespace; equal IDs in other
                                       ///< namespaces remain distinct.
    Field<detail::FieldType::T133> first;   ///< First earlier operand node index.
    Field<detail::FieldType::T133> second;  ///< Second earlier operand node index.
    Field<detail::FieldType::T134> mapping; ///< Ordered mapping entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T135> split_at_attachment; ///< Whether split at attachment is enabled.
};
/// Checked TraceMap payload; all fields retain one immutable case owner.
struct TraceMap {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    TraceMap(Field<detail::FieldType::T136> first_cell_value,
             Field<detail::FieldType::T136> second_cell_value,
             Field<detail::FieldType::T65> domain_value, Field<detail::FieldType::T137> first_value,
             Field<detail::FieldType::T137> second_value,
             Field<detail::FieldType::T74> orientation_value) noexcept
        : first_cell(first_cell_value),
          second_cell(second_cell_value),
          domain(domain_value),
          first(first_value),
          second(second_value),
          orientation(orientation_value) {}
    Field<detail::FieldType::T136>
        first_cell; ///< Ordered first cell components; each numeric element exposes its physical
                    ///< dimensions and exact input value.
    Field<detail::FieldType::T136>
        second_cell; ///< Ordered second cell components; each numeric element exposes its physical
                     ///< dimensions and exact input value.
    Field<detail::FieldType::T65> domain;  ///< Ordered domain entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T137> first;  ///< First earlier operand node index.
    Field<detail::FieldType::T137> second; ///< Second earlier operand node index.
    Field<detail::FieldType::T74> orientation; ///< Orientation; choose one of same,reversed.
};
/// Checked MapTerm payload; all fields retain one immutable case owner.
struct MapTerm {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    MapTerm(Field<detail::FieldType::T63> powers_value,
            Field<detail::FieldType::T62> coefficient_value) noexcept
        : powers(powers_value), coefficient(coefficient_value) {}
    Field<detail::FieldType::T63> powers; ///< Nonnegative monomial exponents, in coordinate order.
    Field<detail::FieldType::T62>
        coefficient; ///< Exact rational coefficient; numerator and positive denominator remain
                     ///< available separately.
};
/// Checked Joint rigid_attachment payload; all fields retain one immutable case owner.
struct JointRigidAttachment {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    JointRigidAttachment(Field<detail::FieldType::T138> kind_value,
                         Field<detail::FieldType::T132> id_value,
                         Field<detail::FieldType::T133> attached_value,
                         Field<detail::FieldType::T30> frame_value,
                         Field<detail::FieldType::T52> offset_value) noexcept
        : kind(kind_value),
          id(id_value),
          attached(attached_value),
          frame(frame_value),
          offset(offset_value) {}
    Field<detail::FieldType::T138> kind; ///< Kind; choose one of rigid_attachment.
    Field<detail::FieldType::T132> id; ///< Id in the Joint identifier namespace; equal IDs in other
                                       ///< namespaces remain distinct.
    Field<detail::FieldType::T133>
        attached; ///< Attached payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T30>
        frame; ///< Frame payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T52>
        offset; ///< Offset payload; record() selects its checked concrete structure.
};
/// Checked Joint weld payload; all fields retain one immutable case owner.
struct JointWeld {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    JointWeld(Field<detail::FieldType::T139> kind_value, Field<detail::FieldType::T132> id_value,
              Field<detail::FieldType::T30> first_value, Field<detail::FieldType::T30> second_value,
              Field<detail::FieldType::T52> relative_value) noexcept
        : kind(kind_value),
          id(id_value),
          first(first_value),
          second(second_value),
          relative(relative_value) {}
    Field<detail::FieldType::T139> kind; ///< Kind; choose one of weld.
    Field<detail::FieldType::T132> id; ///< Id in the Joint identifier namespace; equal IDs in other
                                       ///< namespaces remain distinct.
    Field<detail::FieldType::T30> first;  ///< First earlier operand node index.
    Field<detail::FieldType::T30> second; ///< Second earlier operand node index.
    Field<detail::FieldType::T52>
        relative; ///< Relative payload; record() selects its checked concrete structure.
};
/// Checked Support constraint payload; all fields retain one immutable case owner.
struct SupportConstraint {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SupportConstraint(Field<detail::FieldType::T140> kind_value,
                      Field<detail::FieldType::T141> id_value,
                      Field<detail::FieldType::T142> target_value,
                      Field<detail::FieldType::T143> components_value,
                      Field<detail::FieldType::T144> motion_value) noexcept
        : kind(kind_value),
          id(id_value),
          target(target_value),
          components(components_value),
          motion(motion_value) {}
    Field<detail::FieldType::T140> kind; ///< Kind; choose one of constraint.
    Field<detail::FieldType::T141> id;   ///< Id in the Support identifier namespace; equal IDs in
                                         ///< other namespaces remain distinct.
    Field<detail::FieldType::T142>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T143> components; ///< Ordered components entries; size() and at()
                                               ///< expose the checked sequence without copying.
    Field<detail::FieldType::T144>
        motion; ///< Motion payload; record() selects its checked concrete structure.
};
/// Checked Support trace_constraint payload; all fields retain one immutable case owner.
struct SupportTraceConstraint {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SupportTraceConstraint(Field<detail::FieldType::T145> kind_value,
                           Field<detail::FieldType::T141> id_value,
                           Field<detail::FieldType::T133> region_value,
                           Field<detail::FieldType::T143> components_value,
                           Field<detail::FieldType::T144> motion_value) noexcept
        : kind(kind_value),
          id(id_value),
          region(region_value),
          components(components_value),
          motion(motion_value) {}
    Field<detail::FieldType::T145> kind; ///< Kind; choose one of trace_constraint.
    Field<detail::FieldType::T141> id;   ///< Id in the Support identifier namespace; equal IDs in
                                         ///< other namespaces remain distinct.
    Field<detail::FieldType::T133>
        region; ///< Region payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T143> components; ///< Ordered components entries; size() and at()
                                               ///< expose the checked sequence without copying.
    Field<detail::FieldType::T144>
        motion; ///< Motion payload; record() selects its checked concrete structure.
};
/// Checked Support elastic payload; all fields retain one immutable case owner.
struct SupportElastic {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SupportElastic(Field<detail::FieldType::T146> kind_value,
                   Field<detail::FieldType::T141> id_value,
                   Field<detail::FieldType::T142> target_value,
                   Field<detail::FieldType::T144> reference_value,
                   Field<detail::FieldType::T118> stiffness_value,
                   Field<detail::FieldType::T118> damping_value,
                   Field<detail::FieldType::T52> frame_value) noexcept
        : kind(kind_value),
          id(id_value),
          target(target_value),
          reference(reference_value),
          stiffness(stiffness_value),
          damping(damping_value),
          frame(frame_value) {}
    Field<detail::FieldType::T146> kind; ///< Kind; choose one of elastic.
    Field<detail::FieldType::T141> id;   ///< Id in the Support identifier namespace; equal IDs in
                                         ///< other namespaces remain distinct.
    Field<detail::FieldType::T142>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T144>
        reference; ///< Reference payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T118>
        stiffness; ///< Stiffness as ordered matrix rows; component dimensions follow the
                   ///< corresponding input and output quantities.
    Field<detail::FieldType::T118>
        damping; ///< Damping as ordered matrix rows; component dimensions follow the corresponding
                 ///< input and output quantities.
    Field<detail::FieldType::T52>
        frame; ///< Frame payload; record() selects its checked concrete structure.
};
/// Checked Mass particle payload; all fields retain one immutable case owner.
struct MassParticle {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    MassParticle(Field<detail::FieldType::T147> kind_value, Field<detail::FieldType::T148> id_value,
                 Field<detail::FieldType::T142> target_value,
                 Field<detail::FieldType::T149> mass_value,
                 Field<detail::FieldType::T150> ledger_value) noexcept
        : kind(kind_value),
          id(id_value),
          target(target_value),
          mass(mass_value),
          ledger(ledger_value) {}
    Field<detail::FieldType::T147> kind; ///< Kind; choose one of particle.
    Field<detail::FieldType::T148>
        id; ///< Id in the Mass identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T142>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T149> mass; ///< Exact mass quantity for mass, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T150>
        ledger; ///< Ledger payload; record() selects its checked concrete structure.
};
/// Checked Mass rigid payload; all fields retain one immutable case owner.
struct MassRigid {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    MassRigid(Field<detail::FieldType::T151> kind_value, Field<detail::FieldType::T148> id_value,
              Field<detail::FieldType::T142> target_value,
              Field<detail::FieldType::T149> mass_value, Field<detail::FieldType::T26> center_value,
              Field<detail::FieldType::T152> inertia_value,
              Field<detail::FieldType::T52> frame_value,
              Field<detail::FieldType::T150> ledger_value) noexcept
        : kind(kind_value),
          id(id_value),
          target(target_value),
          mass(mass_value),
          center(center_value),
          inertia(inertia_value),
          frame(frame_value),
          ledger(ledger_value) {}
    Field<detail::FieldType::T151> kind; ///< Kind; choose one of rigid.
    Field<detail::FieldType::T148>
        id; ///< Id in the Mass identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T142>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T149> mass; ///< Exact mass quantity for mass, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T26>
        center; ///< Ordered center components; each numeric element exposes its physical dimensions
                ///< and exact input value.
    Field<detail::FieldType::T152>
        inertia; ///< Inertia as ordered matrix rows; component dimensions follow the corresponding
                 ///< input and output quantities.
    Field<detail::FieldType::T52>
        frame; ///< Frame payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T150>
        ledger; ///< Ledger payload; record() selects its checked concrete structure.
};
/// Checked Ledger add payload; all fields retain one immutable case owner.
struct LedgerAdd {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LedgerAdd(Field<detail::FieldType::T153> kind_value, Field<detail::FieldType::T24> source_value,
              Field<detail::FieldType::T30> owner_value) noexcept
        : kind(kind_value), source(source_value), owner(owner_value) {}
    Field<detail::FieldType::T153> kind;  ///< Kind; choose one of add.
    Field<detail::FieldType::T24> source; ///< Source as decoded UTF-8 text; no normalization or
                                          ///< implicit replacement is performed.
    Field<detail::FieldType::T30>
        owner; ///< Owner payload; record() selects its checked concrete structure.
};
/// Checked Ledger replace payload; all fields retain one immutable case owner.
struct LedgerReplace {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LedgerReplace(Field<detail::FieldType::T154> kind_value,
                  Field<detail::FieldType::T24> source_value,
                  Field<detail::FieldType::T30> owner_value,
                  Field<detail::FieldType::T155> replaced_value) noexcept
        : kind(kind_value), source(source_value), owner(owner_value), replaced(replaced_value) {}
    Field<detail::FieldType::T154> kind;  ///< Kind; choose one of replace.
    Field<detail::FieldType::T24> source; ///< Source as decoded UTF-8 text; no normalization or
                                          ///< implicit replacement is performed.
    Field<detail::FieldType::T30>
        owner; ///< Owner payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T155> replaced; ///< Ordered replaced entries; size() and at() expose
                                             ///< the checked sequence without copying.
};
/// Checked Load point_force payload; all fields retain one immutable case owner.
struct LoadPointForce {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LoadPointForce(Field<detail::FieldType::T156> kind_value,
                   Field<detail::FieldType::T157> id_value,
                   Field<detail::FieldType::T142> target_value,
                   Field<detail::FieldType::T158> frame_value,
                   Field<detail::FieldType::T53> series_value) noexcept
        : kind(kind_value),
          id(id_value),
          target(target_value),
          frame(frame_value),
          series(series_value) {}
    Field<detail::FieldType::T156> kind; ///< Kind; choose one of point_force.
    Field<detail::FieldType::T157>
        id; ///< Id in the Load identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T142>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T158> frame; ///< Frame; choose one of world,body.
    Field<detail::FieldType::T53> series; ///< Series as a function of time in declared case units;
                                          ///< evaluation belongs to preparation.
};
/// Checked Load couple payload; all fields retain one immutable case owner.
struct LoadCouple {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LoadCouple(Field<detail::FieldType::T159> kind_value, Field<detail::FieldType::T157> id_value,
               Field<detail::FieldType::T142> target_value,
               Field<detail::FieldType::T158> frame_value,
               Field<detail::FieldType::T53> series_value) noexcept
        : kind(kind_value),
          id(id_value),
          target(target_value),
          frame(frame_value),
          series(series_value) {}
    Field<detail::FieldType::T159> kind; ///< Kind; choose one of couple.
    Field<detail::FieldType::T157>
        id; ///< Id in the Load identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T142>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T158> frame; ///< Frame; choose one of world,body.
    Field<detail::FieldType::T53> series; ///< Series as a function of time in declared case units;
                                          ///< evaluation belongs to preparation.
};
/// Checked Load gravity payload; all fields retain one immutable case owner.
struct LoadGravity {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LoadGravity(Field<detail::FieldType::T160> kind_value, Field<detail::FieldType::T157> id_value,
                Field<detail::FieldType::T53> acceleration_value) noexcept
        : kind(kind_value), id(id_value), acceleration(acceleration_value) {}
    Field<detail::FieldType::T160> kind; ///< Kind; choose one of gravity.
    Field<detail::FieldType::T157>
        id; ///< Id in the Load identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T53> acceleration; ///< Acceleration as a function of time in declared
                                                ///< case units; evaluation belongs to preparation.
};
/// Checked Load traction payload; all fields retain one immutable case owner.
struct LoadTraction {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LoadTraction(Field<detail::FieldType::T161> kind_value, Field<detail::FieldType::T157> id_value,
                 Field<detail::FieldType::T133> region_value,
                 Field<detail::FieldType::T162> area_value,
                 Field<detail::FieldType::T158> frame_value,
                 Field<detail::FieldType::T163> spatial_value,
                 Field<detail::FieldType::T164> multiplier_value) noexcept
        : kind(kind_value),
          id(id_value),
          region(region_value),
          area(area_value),
          frame(frame_value),
          spatial(spatial_value),
          multiplier(multiplier_value) {}
    Field<detail::FieldType::T161> kind; ///< Kind; choose one of traction.
    Field<detail::FieldType::T157>
        id; ///< Id in the Load identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T133>
        region; ///< Region payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T162> area;  ///< Area; choose one of reference,current.
    Field<detail::FieldType::T158> frame; ///< Frame; choose one of world,body.
    Field<detail::FieldType::T163>
        spatial; ///< Spatial over reference coordinates; evaluation and physical-domain
                 ///< qualification belong to preparation.
    Field<detail::FieldType::T164> multiplier; ///< Multiplier as a function of time in declared
                                               ///< case units; evaluation belongs to preparation.
};
/// Checked Load pressure payload; all fields retain one immutable case owner.
struct LoadPressure {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LoadPressure(Field<detail::FieldType::T165> kind_value, Field<detail::FieldType::T157> id_value,
                 Field<detail::FieldType::T133> region_value,
                 Field<detail::FieldType::T166> spatial_value,
                 Field<detail::FieldType::T164> multiplier_value) noexcept
        : kind(kind_value),
          id(id_value),
          region(region_value),
          spatial(spatial_value),
          multiplier(multiplier_value) {}
    Field<detail::FieldType::T165> kind; ///< Kind; choose one of pressure.
    Field<detail::FieldType::T157>
        id; ///< Id in the Load identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T133>
        region; ///< Region payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T166>
        spatial; ///< Spatial over reference coordinates; evaluation and physical-domain
                 ///< qualification belong to preparation.
    Field<detail::FieldType::T164> multiplier; ///< Multiplier as a function of time in declared
                                               ///< case units; evaluation belongs to preparation.
};
/// Checked Load body payload; all fields retain one immutable case owner.
struct LoadBody {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LoadBody(Field<detail::FieldType::T167> kind_value, Field<detail::FieldType::T157> id_value,
             Field<detail::FieldType::T133> region_value,
             Field<detail::FieldType::T168> spatial_value,
             Field<detail::FieldType::T164> multiplier_value) noexcept
        : kind(kind_value),
          id(id_value),
          region(region_value),
          spatial(spatial_value),
          multiplier(multiplier_value) {}
    Field<detail::FieldType::T167> kind; ///< Kind; choose one of body_force.
    Field<detail::FieldType::T157>
        id; ///< Id in the Load identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T133>
        region; ///< Region payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T168>
        spatial; ///< Spatial over reference coordinates; evaluation and physical-domain
                 ///< qualification belong to preparation.
    Field<detail::FieldType::T164> multiplier; ///< Multiplier as a function of time in declared
                                               ///< case units; evaluation belongs to preparation.
};
/// Checked AeroSurface payload; all fields retain one immutable case owner.
struct AeroSurface {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AeroSurface(Field<detail::FieldType::T169> id_value,
                Field<detail::FieldType::T34> geometry_value,
                Field<detail::FieldType::T170> attachment_value,
                Field<detail::FieldType::T171> chord_value,
                Field<detail::FieldType::T171> span_value,
                Field<detail::FieldType::T35> leading_value,
                Field<detail::FieldType::T35> trailing_value,
                Field<detail::FieldType::T172> orientation_value,
                Field<detail::FieldType::T173> sharing_value) noexcept
        : id(id_value),
          geometry(geometry_value),
          attachment(attachment_value),
          chord(chord_value),
          span(span_value),
          leading(leading_value),
          trailing(trailing_value),
          orientation(orientation_value),
          sharing(sharing_value) {}
    Field<detail::FieldType::T169> id; ///< Id in the AeroSurface identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T34>
        geometry; ///< Geometry referring to a checked Geometry definition in this case.
    Field<detail::FieldType::T170>
        attachment; ///< Attachment payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T171> chord;   ///< Ordered chord entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T171> span;    ///< Ordered span entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T35> leading;  ///< Zero-based index of leading.
    Field<detail::FieldType::T35> trailing; ///< Zero-based index of trailing.
    Field<detail::FieldType::T172> orientation; ///< Orientation; choose one of positive,reversed.
    Field<detail::FieldType::T173> sharing; ///< Ordered sharing entries; size() and at() expose the
                                            ///< checked sequence without copying.
};
/// Checked AeroAttachment fixed payload; all fields retain one immutable case owner.
struct AeroAttachmentFixed {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AeroAttachmentFixed(Field<detail::FieldType::T174> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T174> kind; ///< Kind; choose one of fixed_reference.
};
/// Checked AeroAttachment prescribed payload; all fields retain one immutable case owner.
struct AeroAttachmentPrescribed {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AeroAttachmentPrescribed(Field<detail::FieldType::T175> kind_value,
                             Field<detail::FieldType::T144> motion_value) noexcept
        : kind(kind_value), motion(motion_value) {}
    Field<detail::FieldType::T175> kind; ///< Kind; choose one of prescribed_motion.
    Field<detail::FieldType::T144>
        motion; ///< Motion payload; record() selects its checked concrete structure.
};
/// Checked AeroAttachment structural payload; all fields retain one immutable case owner.
struct AeroAttachmentStructural {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AeroAttachmentStructural(
        Field<detail::FieldType::T176> kind_value, Field<detail::FieldType::T177> host_value,
        std::optional<Field<detail::FieldType::T447>> material_map_value) noexcept
        : kind(kind_value), host(host_value), material_map(material_map_value) {}
    Field<detail::FieldType::T176> kind; ///< Kind; choose one of structural.
    Field<detail::FieldType::T177>
        host; ///< Host referring to a checked Shell definition in this case.
    std::optional<Field<detail::FieldType::T447>>
        material_map; ///< Fixed material chart; absent in legacy host-only inputs.
};
/// Checked AeroSharing payload; all fields retain one immutable case owner.
struct AeroSharing {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AeroSharing(Field<detail::FieldType::T178> other_value,
                Field<detail::FieldType::T35> own_boundary_value,
                Field<detail::FieldType::T35> other_boundary_value,
                Field<detail::FieldType::T74> orientation_value) noexcept
        : other(other_value),
          own_boundary(own_boundary_value),
          other_boundary(other_boundary_value),
          orientation(orientation_value) {}
    Field<detail::FieldType::T178>
        other; ///< Other referring to a checked AeroSurface definition in this case.
    Field<detail::FieldType::T35> own_boundary;   ///< Zero-based index of own boundary.
    Field<detail::FieldType::T35> other_boundary; ///< Zero-based index of other boundary.
    Field<detail::FieldType::T74> orientation;    ///< Orientation; choose one of same,reversed.
};
/// Checked Flow none payload; all fields retain one immutable case owner.
struct FlowNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    FlowNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked Flow known payload; all fields retain one immutable case owner.
struct FlowKnown {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    FlowKnown(Field<detail::FieldType::T180> kind_value,
              Field<detail::FieldType::T93> density_value,
              Field<detail::FieldType::T53> uniform_value,
              Field<detail::FieldType::T181> gust_value,
              Field<detail::FieldType::T182> reference_value,
              Field<detail::FieldType::T48> events_value) noexcept
        : kind(kind_value),
          density(density_value),
          uniform(uniform_value),
          gust(gust_value),
          reference(reference_value),
          events(events_value) {}
    Field<detail::FieldType::T180> kind;   ///< Kind; choose one of known.
    Field<detail::FieldType::T93> density; ///< Exact M/L³ quantity for density, in declared case
                                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T53> uniform; ///< Uniform as a function of time in declared case
                                           ///< units; evaluation belongs to preparation.
    Field<detail::FieldType::T181>
        gust; ///< Gust payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T182>
        reference; ///< Reference payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T48> events; ///< Ordered events entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked Gust none payload; all fields retain one immutable case owner.
struct GustNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    GustNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked Gust sine payload; all fields retain one immutable case owner.
struct GustSine {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    GustSine(Field<detail::FieldType::T183> kind_value,
             Field<detail::FieldType::T164> amplitude_value,
             Field<detail::FieldType::T31> direction_value,
             Field<detail::FieldType::T184> wavevector_value,
             Field<detail::FieldType::T185> angular_frequency_value,
             Field<detail::FieldType::T186> phase_value, Field<detail::FieldType::T26> origin_value,
             Field<detail::FieldType::T44> time_origin_value) noexcept
        : kind(kind_value),
          amplitude(amplitude_value),
          direction(direction_value),
          wavevector(wavevector_value),
          angular_frequency(angular_frequency_value),
          phase(phase_value),
          origin(origin_value),
          time_origin(time_origin_value) {}
    Field<detail::FieldType::T183> kind;      ///< Kind; choose one of sine.
    Field<detail::FieldType::T164> amplitude; ///< Amplitude as a function of time in declared case
                                              ///< units; evaluation belongs to preparation.
    Field<detail::FieldType::T31>
        direction; ///< Ordered direction components; each numeric element exposes its physical
                   ///< dimensions and exact input value.
    Field<detail::FieldType::T184>
        wavevector; ///< Ordered wavevector components; each numeric element exposes its physical
                    ///< dimensions and exact input value.
    Field<detail::FieldType::T185>
        angular_frequency; ///< Exact angle/T quantity for angular frequency, in declared case
                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T186> phase; ///< Exact angle quantity for phase, in declared case
                                          ///< units; number() retains its dimensions.
    Field<detail::FieldType::T26>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T44> time_origin; ///< Exact time quantity for time origin, in declared
                                               ///< case units; number() retains its dimensions.
};
/// Checked FlowReference payload; all fields retain one immutable case owner.
struct FlowReference {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    FlowReference(Field<detail::FieldType::T93> density_value,
                  Field<detail::FieldType::T187> speed_value,
                  Field<detail::FieldType::T188> area_value,
                  Field<detail::FieldType::T61> chord_value,
                  Field<detail::FieldType::T61> span_value,
                  Field<detail::FieldType::T26> moment_origin_value) noexcept
        : density(density_value),
          speed(speed_value),
          area(area_value),
          chord(chord_value),
          span(span_value),
          moment_origin(moment_origin_value) {}
    Field<detail::FieldType::T93> density; ///< Exact M/L³ quantity for density, in declared case
                                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T187> speed; ///< Exact L/T quantity for speed, in declared case units;
                                          ///< number() retains its dimensions.
    Field<detail::FieldType::T188> area;  ///< Exact L² quantity for area, in declared case units;
                                          ///< number() retains its dimensions.
    Field<detail::FieldType::T61> chord;  ///< Exact length quantity for chord, in declared case
                                          ///< units; number() retains its dimensions.
    Field<detail::FieldType::T61> span; ///< Exact length quantity for span, in declared case units;
                                        ///< number() retains its dimensions.
    Field<detail::FieldType::T26>
        moment_origin; ///< Ordered moment origin components; each numeric element exposes its
                       ///< physical dimensions and exact input value.
};
/// Checked Wake none payload; all fields retain one immutable case owner.
struct WakeNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    WakeNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked Wake active payload; all fields retain one immutable case owner.
struct WakeActive {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    WakeActive(Field<detail::FieldType::T189> kind_value,
               Field<detail::FieldType::T190> convection_value,
               Field<detail::FieldType::T191> retention_value,
               Field<detail::FieldType::T61> core_radius_value,
               Field<detail::FieldType::T192> model_value,
               Field<detail::FieldType::T193> mapping_value) noexcept
        : kind(kind_value),
          convection(convection_value),
          retention(retention_value),
          core_radius(core_radius_value),
          model(model_value),
          mapping(mapping_value) {}
    Field<detail::FieldType::T189> kind; ///< Kind; choose one of active.
    Field<detail::FieldType::T190>
        convection; ///< Convection payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T191>
        retention; ///< Retention payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T61>
        core_radius; ///< Exact length quantity for core radius, in declared case units; number()
                     ///< retains its dimensions.
    Field<detail::FieldType::T192> model;   ///< Model; choose one of wp1-stream-p0-span.
    Field<detail::FieldType::T193> mapping; ///< Mapping; choose one of dual-quarter-v1.
};
/// Checked Convection free payload; all fields retain one immutable case owner.
struct ConvectionFree {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ConvectionFree(Field<detail::FieldType::T125> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T125> kind; ///< Kind; choose one of free.
};
/// Checked Convection prescribed payload; all fields retain one immutable case owner.
struct ConvectionPrescribed {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ConvectionPrescribed(Field<detail::FieldType::T194> kind_value,
                         Field<detail::FieldType::T53> velocity_value) noexcept
        : kind(kind_value), velocity(velocity_value) {}
    Field<detail::FieldType::T194> kind;    ///< Kind; choose one of prescribed.
    Field<detail::FieldType::T53> velocity; ///< Velocity as a function of time in declared case
                                            ///< units; evaluation belongs to preparation.
};
/// Checked Retention all payload; all fields retain one immutable case owner.
struct RetentionAll {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    RetentionAll(Field<detail::FieldType::T195> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T195> kind; ///< Kind; choose one of all.
};
/// Checked Retention window payload; all fields retain one immutable case owner.
struct RetentionWindow {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    RetentionWindow(Field<detail::FieldType::T196> kind_value,
                    Field<detail::FieldType::T44> duration_value,
                    Field<detail::FieldType::T61> distance_value,
                    Field<detail::FieldType::T88> error_budget_value) noexcept
        : kind(kind_value),
          duration(duration_value),
          distance(distance_value),
          error_budget(error_budget_value) {}
    Field<detail::FieldType::T196> kind;    ///< Kind; choose one of closed_window.
    Field<detail::FieldType::T44> duration; ///< Exact time quantity for duration, in declared case
                                            ///< units; number() retains its dimensions.
    Field<detail::FieldType::T61> distance; ///< Exact length quantity for distance, in declared
                                            ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        error_budget; ///< Exact dimensionless quantity for error budget, in declared case units;
                      ///< number() retains its dimensions.
};
/// Checked Analysis static payload; all fields retain one immutable case owner.
struct AnalysisStatic {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AnalysisStatic(Field<detail::FieldType::T197> kind_value,
                   Field<detail::FieldType::T198> physics_value,
                   Field<detail::FieldType::T199> branch_value,
                   Field<detail::FieldType::T200> speed_value,
                   Field<detail::FieldType::T201> wake_value,
                   Field<detail::FieldType::T112> loads_value) noexcept
        : kind(kind_value),
          physics(physics_value),
          branch(branch_value),
          speed(speed_value),
          wake(wake_value),
          loads(loads_value) {}
    Field<detail::FieldType::T197> kind; ///< Kind; choose one of static.
    Field<detail::FieldType::T198>
        physics; ///< Physics; choose one of aerodynamic,structural,coupled.
    Field<detail::FieldType::T199>
        branch; ///< Branch payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T200>
        speed; ///< Speed payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T201>
        wake; ///< Wake payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T112> loads; ///< Ordered loads entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked Analysis transient payload; all fields retain one immutable case owner.
struct AnalysisTransient {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AnalysisTransient(Field<detail::FieldType::T202> kind_value,
                      Field<detail::FieldType::T198> physics_value,
                      Field<detail::FieldType::T199> branch_value,
                      Field<detail::FieldType::T203> intervals_value,
                      Field<detail::FieldType::T201> wake_value,
                      Field<detail::FieldType::T112> loads_value) noexcept
        : kind(kind_value),
          physics(physics_value),
          branch(branch_value),
          intervals(intervals_value),
          wake(wake_value),
          loads(loads_value) {}
    Field<detail::FieldType::T202> kind; ///< Kind; choose one of transient.
    Field<detail::FieldType::T198>
        physics; ///< Physics; choose one of aerodynamic,structural,coupled.
    Field<detail::FieldType::T199>
        branch; ///< Branch payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T203> intervals; ///< Ordered intervals entries; size() and at() expose
                                              ///< the checked sequence without copying.
    Field<detail::FieldType::T201>
        wake; ///< Wake payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T112> loads; ///< Ordered loads entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked Analysis flutter payload; all fields retain one immutable case owner.
struct AnalysisFlutter {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AnalysisFlutter(Field<detail::FieldType::T204> kind_value,
                    Field<detail::FieldType::T205> physics_value,
                    Field<detail::FieldType::T206> branches_value,
                    Field<detail::FieldType::T207> intervals_value,
                    Field<detail::FieldType::T44> step_value,
                    Field<detail::FieldType::T201> wake_value,
                    Field<detail::FieldType::T208> search_value,
                    Field<detail::FieldType::T112> loads_value) noexcept
        : kind(kind_value),
          physics(physics_value),
          branches(branches_value),
          intervals(intervals_value),
          step(step_value),
          wake(wake_value),
          search(search_value),
          loads(loads_value) {}
    Field<detail::FieldType::T204> kind;      ///< Kind; choose one of flutter.
    Field<detail::FieldType::T205> physics;   ///< Physics; choose one of coupled.
    Field<detail::FieldType::T206> branches;  ///< Ordered branches entries; size() and at() expose
                                              ///< the checked sequence without copying.
    Field<detail::FieldType::T207> intervals; ///< Ordered intervals entries; size() and at() expose
                                              ///< the checked sequence without copying.
    Field<detail::FieldType::T44> step; ///< Exact time quantity for step, in declared case units;
                                        ///< number() retains its dimensions.
    Field<detail::FieldType::T201>
        wake; ///< Wake payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T208>
        search; ///< Search payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T112> loads; ///< Ordered loads entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked AnalysisSpeed none payload; all fields retain one immutable case owner.
struct AnalysisSpeedNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AnalysisSpeedNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked AnalysisSpeed value payload; all fields retain one immutable case owner.
struct AnalysisSpeedValue {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AnalysisSpeedValue(Field<detail::FieldType::T209> kind_value,
                       Field<detail::FieldType::T187> value_value) noexcept
        : kind(kind_value), value(value_value) {}
    Field<detail::FieldType::T209> kind;  ///< Kind; choose one of value.
    Field<detail::FieldType::T187> value; ///< Exact L/T quantity for value, in declared case units;
                                          ///< number() retains its dimensions.
};
/// Checked Branch payload; all fields retain one immutable case owner.
struct Branch {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Branch(Field<detail::FieldType::T210> id_value, Field<detail::FieldType::T211> seed_value,
           Field<detail::FieldType::T212> continuation_value) noexcept
        : id(id_value), seed(seed_value), continuation(continuation_value) {}
    Field<detail::FieldType::T210> id; ///< Id in the Branch identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T211>
        seed; ///< Seed payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T212>
        continuation; ///< Continuation payload; record() selects its checked concrete structure.
};
/// Checked BranchSeed reference payload; all fields retain one immutable case owner.
struct BranchSeedReference {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BranchSeedReference(Field<detail::FieldType::T213> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T213> kind; ///< Kind; choose one of reference.
};
/// Checked BranchSeed prepared payload; all fields retain one immutable case owner.
struct BranchSeedPrepared {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BranchSeedPrepared(Field<detail::FieldType::T214> kind_value,
                       Field<detail::FieldType::T215> history_value) noexcept
        : kind(kind_value), history(history_value) {}
    Field<detail::FieldType::T214> kind;    ///< Kind; choose one of prepared.
    Field<detail::FieldType::T215> history; ///< History source reference, including owned source
                                            ///< ID, path, digest and expected role.
};
/// Checked Continuation fixed payload; all fields retain one immutable case owner.
struct ContinuationFixed {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ContinuationFixed(Field<detail::FieldType::T216> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T216> kind; ///< Kind; choose one of fixed.
};
/// Checked Continuation sweep payload; all fields retain one immutable case owner.
struct ContinuationSweep {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ContinuationSweep(Field<detail::FieldType::T217> kind_value,
                      Field<detail::FieldType::T187> initial_increment_value,
                      Field<detail::FieldType::T187> minimum_increment_value,
                      Field<detail::FieldType::T187> maximum_increment_value,
                      Field<detail::FieldType::T218> arclength_value,
                      Field<detail::FieldType::T219> maximum_points_value,
                      Field<detail::FieldType::T67> maximum_retries_value) noexcept
        : kind(kind_value),
          initial_increment(initial_increment_value),
          minimum_increment(minimum_increment_value),
          maximum_increment(maximum_increment_value),
          arclength(arclength_value),
          maximum_points(maximum_points_value),
          maximum_retries(maximum_retries_value) {}
    Field<detail::FieldType::T217> kind; ///< Kind; choose one of sweep.
    Field<detail::FieldType::T187>
        initial_increment; ///< Exact L/T quantity for initial increment, in declared case units;
                           ///< number() retains its dimensions.
    Field<detail::FieldType::T187>
        minimum_increment; ///< Exact L/T quantity for minimum increment, in declared case units;
                           ///< number() retains its dimensions.
    Field<detail::FieldType::T187>
        maximum_increment; ///< Exact L/T quantity for maximum increment, in declared case units;
                           ///< number() retains its dimensions.
    Field<detail::FieldType::T218>
        arclength; ///< Arclength payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T219> maximum_points; ///< Maximum points payload; record() selects its
                                                   ///< checked concrete structure.
    Field<detail::FieldType::T67> maximum_retries; ///< Nonnegative count of maximum retries.
};
/// Checked ArcLength payload; all fields retain one immutable case owner.
struct ArcLength {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ArcLength(Field<detail::FieldType::T220> state_scales_value,
              Field<detail::FieldType::T187> speed_scale_value,
              Field<detail::FieldType::T221> orientation_value,
              Field<detail::FieldType::T135> enable_at_limit_value) noexcept
        : state_scales(state_scales_value),
          speed_scale(speed_scale_value),
          orientation(orientation_value),
          enable_at_limit(enable_at_limit_value) {}
    Field<detail::FieldType::T220> state_scales; ///< Ordered state scales entries; size() and at()
                                                 ///< expose the checked sequence without copying.
    Field<detail::FieldType::T187> speed_scale; ///< Exact L/T quantity for speed scale, in declared
                                                ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T221> orientation; ///< Orientation; choose one of previous_tangent.
    Field<detail::FieldType::T135> enable_at_limit; ///< Whether enable at limit is enabled.
};
/// Checked TimeInterval payload; all fields retain one immutable case owner.
struct TimeInterval {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    TimeInterval(Field<detail::FieldType::T44> start_value, Field<detail::FieldType::T44> end_value,
                 Field<detail::FieldType::T44> step_value,
                 Field<detail::FieldType::T222> restart_value) noexcept
        : start(start_value), end(end_value), step(step_value), restart(restart_value) {}
    Field<detail::FieldType::T44> start; ///< Exact time quantity for start, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T44> end;   ///< Exact time quantity for end, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T44> step;  ///< Exact time quantity for step, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T222>
        restart; ///< Restart; choose one of initial,event,step_change,chart_change.
};
/// Checked SpeedInterval payload; all fields retain one immutable case owner.
struct SpeedInterval {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SpeedInterval(Field<detail::FieldType::T187> lower_value,
                  Field<detail::FieldType::T187> upper_value) noexcept
        : lower(lower_value), upper(upper_value) {}
    Field<detail::FieldType::T187> lower; ///< Exact L/T quantity for lower, in declared case units;
                                          ///< number() retains its dimensions.
    Field<detail::FieldType::T187> upper; ///< Exact L/T quantity for upper, in declared case units;
                                          ///< number() retains its dimensions.
};
/// Checked Search payload; all fields retain one immutable case owner.
struct Search {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Search(Field<detail::FieldType::T223> queue_value,
           Field<detail::FieldType::T187> speed_absolute_value,
           Field<detail::FieldType::T88> speed_relative_value,
           Field<detail::FieldType::T224> growth_budget_value,
           Field<detail::FieldType::T225> mode_association_value,
           Field<detail::FieldType::T219> maximum_points_value,
           Field<detail::FieldType::T67> maximum_retries_value,
           Field<detail::FieldType::T226> coverage_value) noexcept
        : queue(queue_value),
          speed_absolute(speed_absolute_value),
          speed_relative(speed_relative_value),
          growth_budget(growth_budget_value),
          mode_association(mode_association_value),
          maximum_points(maximum_points_value),
          maximum_retries(maximum_retries_value),
          coverage(coverage_value) {}
    Field<detail::FieldType::T223> queue; ///< Queue; choose one of lowest_interval_first.
    Field<detail::FieldType::T187>
        speed_absolute; ///< Exact L/T quantity for speed absolute, in declared case units; number()
                        ///< retains its dimensions.
    Field<detail::FieldType::T88>
        speed_relative; ///< Exact dimensionless quantity for speed relative, in declared case
                        ///< units; number() retains its dimensions.
    Field<detail::FieldType::T224>
        growth_budget; ///< Exact 1/T quantity for growth budget, in declared case units; number()
                       ///< retains its dimensions.
    Field<detail::FieldType::T225> mode_association; ///< Mode association payload; record() selects
                                                     ///< its checked concrete structure.
    Field<detail::FieldType::T219> maximum_points; ///< Maximum points payload; record() selects its
                                                   ///< checked concrete structure.
    Field<detail::FieldType::T67> maximum_retries; ///< Nonnegative count of maximum retries.
    Field<detail::FieldType::T226>
        coverage; ///< Coverage payload; record() selects its checked concrete structure.
};
/// Checked Association payload; all fields retain one immutable case owner.
struct Association {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Association(Field<detail::FieldType::T89> physical_observations_value,
                Field<detail::FieldType::T88> mac_minimum_value,
                Field<detail::FieldType::T186> subspace_angle_maximum_value,
                Field<detail::FieldType::T224> eigenvalue_absolute_value,
                Field<detail::FieldType::T88> eigenvalue_relative_value) noexcept
        : physical_observations(physical_observations_value),
          mac_minimum(mac_minimum_value),
          subspace_angle_maximum(subspace_angle_maximum_value),
          eigenvalue_absolute(eigenvalue_absolute_value),
          eigenvalue_relative(eigenvalue_relative_value) {}
    Field<detail::FieldType::T89>
        physical_observations; ///< Ordered physical observations entries; size() and at() expose
                               ///< the checked sequence without copying.
    Field<detail::FieldType::T88>
        mac_minimum; ///< Exact dimensionless quantity for mac minimum, in declared case units;
                     ///< number() retains its dimensions.
    Field<detail::FieldType::T186>
        subspace_angle_maximum; ///< Exact angle quantity for subspace angle maximum, in declared
                                ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T224>
        eigenvalue_absolute; ///< Exact 1/T quantity for eigenvalue absolute, in declared case
                             ///< units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        eigenvalue_relative; ///< Exact dimensionless quantity for eigenvalue relative, in declared
                             ///< case units; number() retains its dimensions.
};
/// Checked Coverage payload; all fields retain one immutable case owner.
struct Coverage {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Coverage(Field<detail::FieldType::T227> branches_value,
             Field<detail::FieldType::T89> mode_observations_value,
             Field<detail::FieldType::T135> require_physical_value) noexcept
        : branches(branches_value),
          mode_observations(mode_observations_value),
          require_physical(require_physical_value) {}
    Field<detail::FieldType::T227> branches; ///< Ordered branches entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T89>
        mode_observations; ///< Ordered mode observations entries; size() and at() expose the
                           ///< checked sequence without copying.
    Field<detail::FieldType::T135> require_physical; ///< Whether require physical is enabled.
};
/// Checked Initial equilibrium payload; all fields retain one immutable case owner.
struct InitialEquilibrium {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialEquilibrium(Field<detail::FieldType::T228> kind_value,
                       Field<detail::FieldType::T229> branch_value,
                       Field<detail::FieldType::T230> perturbation_value) noexcept
        : kind(kind_value), branch(branch_value), perturbation(perturbation_value) {}
    Field<detail::FieldType::T228> kind; ///< Kind; choose one of equilibrium.
    Field<detail::FieldType::T229>
        branch; ///< Branch referring to a checked Branch definition in this case.
    Field<detail::FieldType::T230>
        perturbation; ///< Perturbation payload; record() selects its checked concrete structure.
};
/// Checked Initial cold payload; all fields retain one immutable case owner.
struct InitialCold {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialCold(Field<detail::FieldType::T231> kind_value,
                Field<detail::FieldType::T232> state_value,
                Field<detail::FieldType::T233> circulation_value,
                Field<detail::FieldType::T234> wake_value) noexcept
        : kind(kind_value), state(state_value), circulation(circulation_value), wake(wake_value) {}
    Field<detail::FieldType::T231> kind; ///< Kind; choose one of cold_start.
    Field<detail::FieldType::T232>
        state; ///< State payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T233> circulation; ///< Ordered circulation entries; size() and at()
                                                ///< expose the checked sequence without copying.
    Field<detail::FieldType::T234>
        wake; ///< Wake payload; record() selects its checked concrete structure.
};
/// Checked Initial prepared payload; all fields retain one immutable case owner.
struct InitialPrepared {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialPrepared(Field<detail::FieldType::T235> kind_value,
                    Field<detail::FieldType::T236> file_value) noexcept
        : kind(kind_value), file(file_value) {}
    Field<detail::FieldType::T235> kind; ///< Kind; choose one of prepared_history.
    Field<detail::FieldType::T236>
        file; ///< File source reference, including owned source ID, path, digest and expected role.
};
/// Checked Perturbation none payload; all fields retain one immutable case owner.
struct PerturbationNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PerturbationNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked InitialState payload; all fields retain one immutable case owner.
struct InitialState {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialState(Field<detail::FieldType::T237> shells_value,
                 Field<detail::FieldType::T238> beams_value,
                 Field<detail::FieldType::T239> rigids_value) noexcept
        : shells(shells_value), beams(beams_value), rigids(rigids_value) {}
    Field<detail::FieldType::T237> shells; ///< Ordered shells entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T238> beams;  ///< Ordered beams entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T239> rigids; ///< Ordered rigids entries; size() and at() expose the
                                           ///< checked sequence without copying.
};
/// Checked CoefficientState payload; all fields retain one immutable case owner.
struct CoefficientState {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    CoefficientState(Field<detail::FieldType::T177> shell_value,
                     Field<detail::FieldType::T240> reference_value,
                     Field<detail::FieldType::T241> representation_value,
                     Field<detail::FieldType::T242> position_value,
                     Field<detail::FieldType::T242> velocity_value) noexcept
        : shell(shell_value),
          reference(reference_value),
          representation(representation_value),
          position(position_value),
          velocity(velocity_value) {}
    Field<detail::FieldType::T177>
        shell; ///< Shell referring to a checked Shell definition in this case.
    Field<detail::FieldType::T240> reference; ///< Reference as decoded UTF-8 text; no normalization
                                              ///< or implicit replacement is performed.
    Field<detail::FieldType::T241>
        representation; ///< Representation; choose one of absolute,reference_increment.
    Field<detail::FieldType::T242>
        position; ///< Position coefficient planes in the declared absolute or reference-increment
                  ///< representation; layout is qualified during preparation.
    Field<detail::FieldType::T242> velocity; ///< Velocity coefficient planes using the same basis
                                             ///< and independent ordering as position.
};
/// Checked Planes binary64 payload; all fields retain one immutable case owner.
struct PlanesBinary64 {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PlanesBinary64(Field<detail::FieldType::T243> kind_value,
                   Field<detail::FieldType::T171> high_value) noexcept
        : kind(kind_value), high(high_value) {}
    Field<detail::FieldType::T243> kind; ///< Kind; choose one of binary64.
    Field<detail::FieldType::T171> high; ///< Ordered high entries; size() and at() expose the
                                         ///< checked sequence without copying.
};
/// Checked Planes two_word payload; all fields retain one immutable case owner.
struct PlanesTwoWord {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PlanesTwoWord(Field<detail::FieldType::T244> kind_value,
                  Field<detail::FieldType::T171> high_value,
                  Field<detail::FieldType::T171> low_value) noexcept
        : kind(kind_value), high(high_value), low(low_value) {}
    Field<detail::FieldType::T244> kind; ///< Kind; choose one of two_word.
    Field<detail::FieldType::T171> high; ///< Ordered high entries; size() and at() expose the
                                         ///< checked sequence without copying.
    Field<detail::FieldType::T171>
        low; ///< Ordered low entries; size() and at() expose the checked sequence without copying.
};
/// Checked BeamState payload; all fields retain one immutable case owner.
struct BeamState {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BeamState(Field<detail::FieldType::T245> beam_value,
              Field<detail::FieldType::T246> nodes_value) noexcept
        : beam(beam_value), nodes(nodes_value) {}
    Field<detail::FieldType::T245>
        beam; ///< Beam referring to a checked Beam definition in this case.
    Field<detail::FieldType::T246>
        nodes; ///< Ordered nodes components; each numeric element exposes its physical dimensions
               ///< and exact input value.
};
/// Checked RigidState payload; all fields retain one immutable case owner.
struct RigidState {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    RigidState(Field<detail::FieldType::T247> rigid_value,
               Field<detail::FieldType::T248> state_value) noexcept
        : rigid(rigid_value), state(state_value) {}
    Field<detail::FieldType::T247>
        rigid; ///< Rigid referring to a checked Mass definition in this case.
    Field<detail::FieldType::T248>
        state; ///< State payload; record() selects its checked concrete structure.
};
/// Checked FrameStatePayload payload; all fields retain one immutable case owner.
struct FrameStatePayload {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    FrameStatePayload(Field<detail::FieldType::T26> position_value,
                      Field<detail::FieldType::T27> rotation_value,
                      Field<detail::FieldType::T249> linear_velocity_value,
                      Field<detail::FieldType::T250> angular_velocity_value) noexcept
        : position(position_value),
          rotation(rotation_value),
          linear_velocity(linear_velocity_value),
          angular_velocity(angular_velocity_value) {}
    Field<detail::FieldType::T26> position; ///< Reference-frame position, in the case length unit.
    Field<detail::FieldType::T27> rotation; ///< Rotation matrix with ordered rows; frame regularity
                                            ///< is checked during preparation.
    Field<detail::FieldType::T249>
        linear_velocity; ///< Ordered linear velocity components; each numeric element exposes its
                         ///< physical dimensions and exact input value.
    Field<detail::FieldType::T250>
        angular_velocity; ///< Body-frame angular velocity components, in angle per case time unit.
};
/// Checked InitialWake payload; all fields retain one immutable case owner.
struct InitialWake {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialWake(Field<detail::FieldType::T251> surfaces_value) noexcept
        : surfaces(surfaces_value) {}
    Field<detail::FieldType::T251> surfaces; ///< Ordered surfaces entries; size() and at() expose
                                             ///< the checked sequence without copying.
};
/// Checked SurfaceWake payload; all fields retain one immutable case owner.
struct SurfaceWake {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    SurfaceWake(Field<detail::FieldType::T178> surface_value,
                Field<detail::FieldType::T67> rows_value,
                Field<detail::FieldType::T252> points_value,
                Field<detail::FieldType::T233> coefficients_value,
                Field<detail::FieldType::T253> row_ids_value,
                Field<detail::FieldType::T254> point_ids_value,
                Field<detail::FieldType::T255> band_ids_value,
                Field<detail::FieldType::T48> birth_times_value) noexcept
        : surface(surface_value),
          rows(rows_value),
          points(points_value),
          coefficients(coefficients_value),
          row_ids(row_ids_value),
          point_ids(point_ids_value),
          band_ids(band_ids_value),
          birth_times(birth_times_value) {}
    Field<detail::FieldType::T178>
        surface; ///< Surface referring to a checked AeroSurface definition in this case.
    Field<detail::FieldType::T67> rows;    ///< Nonnegative count of rows.
    Field<detail::FieldType::T252> points; ///< Ordered points entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T233> coefficients; ///< Ordered coefficients entries; size() and at()
                                                 ///< expose the checked sequence without copying.
    Field<detail::FieldType::T253> row_ids; ///< Ordered row ids entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T254> point_ids; ///< Ordered point ids entries; size() and at() expose
                                              ///< the checked sequence without copying.
    Field<detail::FieldType::T255> band_ids;  ///< Ordered band ids entries; size() and at() expose
                                              ///< the checked sequence without copying.
    Field<detail::FieldType::T48> birth_times; ///< Ordered birth times entries; size() and at()
                                               ///< expose the checked sequence without copying.
};
/// Checked Restart payload; all fields retain one immutable case owner.
struct Restart {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Restart(Field<detail::FieldType::T24> directory_value,
            Field<detail::FieldType::T240> manifest_sha256_value,
            Field<detail::FieldType::T256> compatibility_value,
            Field<detail::FieldType::T257> budget_value) noexcept
        : directory(directory_value),
          manifest_sha256(manifest_sha256_value),
          compatibility(compatibility_value),
          budget(budget_value) {}
    Field<detail::FieldType::T24> directory; ///< Directory as decoded UTF-8 text; no normalization
                                             ///< or implicit replacement is performed.
    Field<detail::FieldType::T240>
        manifest_sha256; ///< Manifest sha256 as decoded UTF-8 text; no normalization or implicit
                         ///< replacement is performed.
    Field<detail::FieldType::T256>
        compatibility; ///< Compatibility; choose one of
                       ///< same_environment_bitwise,qualified_cross_environment.
    Field<detail::FieldType::T257> budget; ///< Ordered budget entries; size() and at() expose the
                                           ///< checked sequence without copying.
};
/// Checked Observation point payload; all fields retain one immutable case owner.
struct ObservationPoint {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ObservationPoint(Field<detail::FieldType::T258> kind_value,
                     Field<detail::FieldType::T259> id_value,
                     Field<detail::FieldType::T142> point_value,
                     Field<detail::FieldType::T260> quantity_value,
                     Field<detail::FieldType::T261> basis_value,
                     Field<detail::FieldType::T262> comparison_value,
                     Field<detail::FieldType::T263> budget_value) noexcept
        : kind(kind_value),
          id(id_value),
          point(point_value),
          quantity(quantity_value),
          basis(basis_value),
          comparison(comparison_value),
          budget(budget_value) {}
    Field<detail::FieldType::T258> kind; ///< Kind; choose one of point.
    Field<detail::FieldType::T259> id; ///< Id in the Observation identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T142>
        point; ///< Point payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T260>
        quantity; ///< Quantity; choose one of position,displacement,velocity,acceleration,force.
    Field<detail::FieldType::T261>
        basis; ///< Basis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T262>
        comparison; ///< Comparison payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T263>
        budget; ///< Budget payload; record() selects its checked concrete structure.
};
/// Checked Observation stress payload; all fields retain one immutable case owner.
struct ObservationStress {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ObservationStress(Field<detail::FieldType::T264> kind_value,
                      Field<detail::FieldType::T259> id_value,
                      Field<detail::FieldType::T142> point_value,
                      Field<detail::FieldType::T265> layer_value,
                      Field<detail::FieldType::T266> interface_value,
                      Field<detail::FieldType::T267> measure_value,
                      Field<detail::FieldType::T268> contribution_value,
                      Field<detail::FieldType::T261> basis_value,
                      Field<detail::FieldType::T262> comparison_value,
                      Field<detail::FieldType::T263> budget_value) noexcept
        : kind(kind_value),
          id(id_value),
          point(point_value),
          layer(layer_value),
          interface(interface_value),
          measure(measure_value),
          contribution(contribution_value),
          basis(basis_value),
          comparison(comparison_value),
          budget(budget_value) {}
    Field<detail::FieldType::T264> kind; ///< Kind; choose one of stress.
    Field<detail::FieldType::T259> id; ///< Id in the Observation identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T142>
        point; ///< Point payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T265>
        layer; ///< Layer referring to a checked Layer definition in this case.
    Field<detail::FieldType::T266>
        interface; ///< Interface; choose one of interior,lower,upper,both.
    Field<detail::FieldType::T267>
        measure; ///< Measure; choose one of green,second_pk,first_pk,cauchy.
    Field<detail::FieldType::T268>
        contribution; ///< Contribution; choose one of elastic,viscous,total.
    Field<detail::FieldType::T261>
        basis; ///< Basis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T262>
        comparison; ///< Comparison payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T263>
        budget; ///< Budget payload; record() selects its checked concrete structure.
};
/// Checked Observation section payload; all fields retain one immutable case owner.
struct ObservationSection {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ObservationSection(Field<detail::FieldType::T269> kind_value,
                       Field<detail::FieldType::T259> id_value,
                       Field<detail::FieldType::T133> region_value,
                       Field<detail::FieldType::T26> origin_value,
                       Field<detail::FieldType::T261> basis_value,
                       Field<detail::FieldType::T270> quantity_value,
                       Field<detail::FieldType::T262> comparison_value,
                       Field<detail::FieldType::T263> budget_value) noexcept
        : kind(kind_value),
          id(id_value),
          region(region_value),
          origin(origin_value),
          basis(basis_value),
          quantity(quantity_value),
          comparison(comparison_value),
          budget(budget_value) {}
    Field<detail::FieldType::T269> kind; ///< Kind; choose one of section.
    Field<detail::FieldType::T259> id; ///< Id in the Observation identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T133>
        region; ///< Region payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T26>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T261>
        basis; ///< Basis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T270> quantity; ///< Quantity; choose one of force,moment.
    Field<detail::FieldType::T262>
        comparison; ///< Comparison payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T263>
        budget; ///< Budget payload; record() selects its checked concrete structure.
};
/// Checked Observation reaction payload; all fields retain one immutable case owner.
struct ObservationReaction {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ObservationReaction(Field<detail::FieldType::T271> kind_value,
                        Field<detail::FieldType::T259> id_value,
                        Field<detail::FieldType::T272> support_value,
                        Field<detail::FieldType::T26> origin_value,
                        Field<detail::FieldType::T261> basis_value,
                        Field<detail::FieldType::T262> comparison_value,
                        Field<detail::FieldType::T263> budget_value) noexcept
        : kind(kind_value),
          id(id_value),
          support(support_value),
          origin(origin_value),
          basis(basis_value),
          comparison(comparison_value),
          budget(budget_value) {}
    Field<detail::FieldType::T271> kind; ///< Kind; choose one of reaction.
    Field<detail::FieldType::T259> id; ///< Id in the Observation identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T272>
        support; ///< Support referring to a checked Support definition in this case.
    Field<detail::FieldType::T26>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T261>
        basis; ///< Basis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T262>
        comparison; ///< Comparison payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T263>
        budget; ///< Budget payload; record() selects its checked concrete structure.
};
/// Checked Observation aero payload; all fields retain one immutable case owner.
struct ObservationAero {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ObservationAero(Field<detail::FieldType::T273> kind_value,
                    Field<detail::FieldType::T259> id_value,
                    Field<detail::FieldType::T274> surfaces_value,
                    Field<detail::FieldType::T26> origin_value,
                    Field<detail::FieldType::T261> basis_value,
                    Field<detail::FieldType::T275> quantity_value,
                    Field<detail::FieldType::T262> comparison_value,
                    Field<detail::FieldType::T263> budget_value) noexcept
        : kind(kind_value),
          id(id_value),
          surfaces(surfaces_value),
          origin(origin_value),
          basis(basis_value),
          quantity(quantity_value),
          comparison(comparison_value),
          budget(budget_value) {}
    Field<detail::FieldType::T273> kind; ///< Kind; choose one of aero.
    Field<detail::FieldType::T259> id; ///< Id in the Observation identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T274> surfaces; ///< Ordered surfaces entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T26>
        origin; ///< Coordinate or time origin in the field's declared case units.
    Field<detail::FieldType::T261>
        basis; ///< Basis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T275> quantity; ///< Quantity; choose one of force,moment,coefficient.
    Field<detail::FieldType::T262>
        comparison; ///< Comparison payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T263>
        budget; ///< Budget payload; record() selects its checked concrete structure.
};
/// Checked Observation mode payload; all fields retain one immutable case owner.
struct ObservationMode {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ObservationMode(Field<detail::FieldType::T276> kind_value,
                    Field<detail::FieldType::T259> id_value,
                    Field<detail::FieldType::T277> targets_value,
                    Field<detail::FieldType::T278> points_value,
                    Field<detail::FieldType::T261> basis_value,
                    Field<detail::FieldType::T279> quantity_value,
                    Field<detail::FieldType::T262> comparison_value,
                    Field<detail::FieldType::T263> budget_value) noexcept
        : kind(kind_value),
          id(id_value),
          targets(targets_value),
          points(points_value),
          basis(basis_value),
          quantity(quantity_value),
          comparison(comparison_value),
          budget(budget_value) {}
    Field<detail::FieldType::T276> kind; ///< Kind; choose one of mode.
    Field<detail::FieldType::T259> id; ///< Id in the Observation identifier namespace; equal IDs in
                                       ///< other namespaces remain distinct.
    Field<detail::FieldType::T277> targets; ///< Ordered targets entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T278> points;  ///< Ordered points entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T261>
        basis; ///< Basis payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T279> quantity; ///< Quantity; choose one of frequency,growth,subspace.
    Field<detail::FieldType::T262>
        comparison; ///< Comparison payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T263>
        budget; ///< Budget payload; record() selects its checked concrete structure.
};
/// Checked Basis world/material payload; all fields retain one immutable case owner.
struct BasisWorldMaterial {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BasisWorldMaterial(Field<detail::FieldType::T280> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T280>
        kind; ///< Kind; choose one of world,reference_material,current_material.
};
/// Checked Basis explicit payload; all fields retain one immutable case owner.
struct BasisExplicit {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BasisExplicit(Field<detail::FieldType::T57> kind_value,
                  Field<detail::FieldType::T27> rotation_value) noexcept
        : kind(kind_value), rotation(rotation_value) {}
    Field<detail::FieldType::T57> kind;     ///< Kind; choose one of explicit.
    Field<detail::FieldType::T27> rotation; ///< Rotation matrix with ordered rows; frame regularity
                                            ///< is checked during preparation.
};
/// Checked Comparison payload; all fields retain one immutable case owner.
struct Comparison {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Comparison(Field<detail::FieldType::T46> window_value,
               Field<detail::FieldType::T281> statistic_value,
               Field<detail::FieldType::T282> side_value,
               Field<detail::FieldType::T283> reference_value) noexcept
        : window(window_value),
          statistic(statistic_value),
          side(side_value),
          reference(reference_value) {}
    Field<detail::FieldType::T46>
        window; ///< Ordered window components; each numeric element exposes its physical dimensions
                ///< and exact input value.
    Field<detail::FieldType::T281>
        statistic; ///< Statistic; choose one of instantaneous,history,peak,rms.
    Field<detail::FieldType::T282> side; ///< Side; choose one of left,right,both.
    Field<detail::FieldType::T283>
        reference; ///< Reference payload; record() selects its checked concrete structure.
};
/// Checked ComparisonReference none payload; all fields retain one immutable case owner.
struct ComparisonReferenceNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ComparisonReferenceNone(Field<detail::FieldType::T179> kind_value) noexcept
        : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked ComparisonReference file payload; all fields retain one immutable case owner.
struct ComparisonReferenceFile {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ComparisonReferenceFile(Field<detail::FieldType::T284> kind_value,
                            Field<detail::FieldType::T285> source_value,
                            Field<detail::FieldType::T25> columns_value,
                            Field<detail::FieldType::T171> uncertainty_value) noexcept
        : kind(kind_value),
          source(source_value),
          columns(columns_value),
          uncertainty(uncertainty_value) {}
    Field<detail::FieldType::T284> kind;   ///< Kind; choose one of file.
    Field<detail::FieldType::T285> source; ///< Source source reference, including owned source ID,
                                           ///< path, digest and expected role.
    Field<detail::FieldType::T25> columns; ///< Ordered columns entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T171> uncertainty; ///< Ordered uncertainty entries; size() and at()
                                                ///< expose the checked sequence without copying.
};
/// Checked ObservationBudget payload; all fields retain one immutable case owner.
struct ObservationBudget {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ObservationBudget(Field<detail::FieldType::T171> absolute_value,
                      Field<detail::FieldType::T171> relative_value,
                      Field<detail::FieldType::T171> scale_value,
                      Field<detail::FieldType::T286> groups_value) noexcept
        : absolute(absolute_value),
          relative(relative_value),
          scale(scale_value),
          groups(groups_value) {}
    Field<detail::FieldType::T171> absolute; ///< Ordered absolute entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T171> relative; ///< Ordered relative entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T171> scale;    ///< Ordered scale entries; size() and at() expose the
                                             ///< checked sequence without copying.
    Field<detail::FieldType::T286> groups;   ///< Ordered groups entries; size() and at() expose the
                                             ///< checked sequence without copying.
};
/// Checked ErrorGroup not_applicable payload; all fields retain one immutable case owner.
struct ErrorGroupNotApplicable {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ErrorGroupNotApplicable(Field<detail::FieldType::T287> status_value,
                            Field<detail::FieldType::T288> kind_value,
                            Field<detail::FieldType::T24> reason_value) noexcept
        : status(status_value), kind(kind_value), reason(reason_value) {}
    Field<detail::FieldType::T287> status; ///< Status; choose one of not_applicable.
    Field<detail::FieldType::T288>
        kind; ///< Kind; choose one of
              ///< geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping.
    Field<detail::FieldType::T24> reason; ///< Reason as decoded UTF-8 text; no normalization or
                                          ///< implicit replacement is performed.
};
/// Checked ErrorGroup covered_by payload; all fields retain one immutable case owner.
struct ErrorGroupCoveredBy {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ErrorGroupCoveredBy(Field<detail::FieldType::T289> status_value,
                        Field<detail::FieldType::T288> kind_value,
                        Field<detail::FieldType::T288> representative_value) noexcept
        : status(status_value), kind(kind_value), representative(representative_value) {}
    Field<detail::FieldType::T289> status; ///< Status; choose one of covered_by.
    Field<detail::FieldType::T288>
        kind; ///< Kind; choose one of
              ///< geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping.
    Field<detail::FieldType::T288>
        representative; ///< Representative; choose one of
                        ///< geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping.
};
/// Checked ErrorGroup required payload; all fields retain one immutable case owner.
struct ErrorGroupRequired {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ErrorGroupRequired(Field<detail::FieldType::T290> status_value,
                       Field<detail::FieldType::T288> kind_value,
                       Field<detail::FieldType::T291> covers_value,
                       Field<detail::FieldType::T171> absolute_value,
                       Field<detail::FieldType::T171> relative_value) noexcept
        : status(status_value),
          kind(kind_value),
          covers(covers_value),
          absolute(absolute_value),
          relative(relative_value) {}
    Field<detail::FieldType::T290> status; ///< Status; choose one of required.
    Field<detail::FieldType::T288>
        kind; ///< Kind; choose one of
              ///< geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping.
    Field<detail::FieldType::T291> covers;   ///< Ordered covers entries; size() and at() expose the
                                             ///< checked sequence without copying.
    Field<detail::FieldType::T171> absolute; ///< Ordered absolute entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T171> relative; ///< Ordered relative entries; size() and at() expose
                                             ///< the checked sequence without copying.
};
/// Checked Numerics payload; all fields retain one immutable case owner.
struct Numerics {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Numerics(
        Field<detail::FieldType::T292> patch_value, Field<detail::FieldType::T293> quadrature_value,
        Field<detail::FieldType::T294> linear_value, Field<detail::FieldType::T295> newton_value,
        Field<detail::FieldType::T296> precision_value, Field<detail::FieldType::T297> eigen_value,
        Field<detail::FieldType::T298> wake_rows_value, Field<detail::FieldType::T220> scales_value,
        Field<detail::FieldType::T286> error_groups_value) noexcept
        : patch(patch_value),
          quadrature(quadrature_value),
          linear(linear_value),
          newton(newton_value),
          precision(precision_value),
          eigen(eigen_value),
          wake_rows(wake_rows_value),
          scales(scales_value),
          error_groups(error_groups_value) {}
    Field<detail::FieldType::T292>
        patch; ///< Patch payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T293>
        quadrature; ///< Quadrature payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T294>
        linear; ///< Linear payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T295>
        newton; ///< Newton payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T296>
        precision; ///< Precision payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T297>
        eigen; ///< Eigen payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T298>
        wake_rows; ///< Wake rows payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T220> scales; ///< Ordered scales entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T286> error_groups; ///< Ordered error groups entries; size() and at()
                                                 ///< expose the checked sequence without copying.
};
/// Checked WakeCapacity none payload; all fields retain one immutable case owner.
struct WakeCapacityNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    WakeCapacityNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked WakeCapacity limited payload; all fields retain one immutable case owner.
struct WakeCapacityLimited {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    WakeCapacityLimited(Field<detail::FieldType::T299> kind_value,
                        Field<detail::FieldType::T219> maximum_rows_value) noexcept
        : kind(kind_value), maximum_rows(maximum_rows_value) {}
    Field<detail::FieldType::T299> kind; ///< Kind; choose one of limited.
    Field<detail::FieldType::T219>
        maximum_rows; ///< Maximum rows payload; record() selects its checked concrete structure.
};
/// Checked PatchPolicy none payload; all fields retain one immutable case owner.
struct PatchPolicyNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PatchPolicyNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked PatchPolicy active payload; all fields retain one immutable case owner.
struct PatchPolicyActive {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PatchPolicyActive(Field<detail::FieldType::T189> kind_value,
                      Field<detail::FieldType::T67> maximum_degree_value,
                      Field<detail::FieldType::T219> maximum_partition_value,
                      Field<detail::FieldType::T61> geometry_error_value,
                      Field<detail::FieldType::T61> trace_error_value,
                      Field<detail::FieldType::T300> rank_value) noexcept
        : kind(kind_value),
          maximum_degree(maximum_degree_value),
          maximum_partition(maximum_partition_value),
          geometry_error(geometry_error_value),
          trace_error(trace_error_value),
          rank(rank_value) {}
    Field<detail::FieldType::T189> kind;              ///< Kind; choose one of active.
    Field<detail::FieldType::T67> maximum_degree;     ///< Nonnegative count of maximum degree.
    Field<detail::FieldType::T219> maximum_partition; ///< Maximum partition payload; record()
                                                      ///< selects its checked concrete structure.
    Field<detail::FieldType::T61>
        geometry_error; ///< Exact length quantity for geometry error, in declared case units;
                        ///< number() retains its dimensions.
    Field<detail::FieldType::T61>
        trace_error; ///< Exact length quantity for trace error, in declared case units; number()
                     ///< retains its dimensions.
    Field<detail::FieldType::T300> rank; ///< Rank; choose one of exact_no_drop.
};
/// Checked Quadrature payload; all fields retain one immutable case owner.
struct Quadrature {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Quadrature(Field<detail::FieldType::T301> structure_value,
               Field<detail::FieldType::T302> aero_value) noexcept
        : structure(structure_value), aero(aero_value) {}
    Field<detail::FieldType::T301>
        structure; ///< Structure payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T302>
        aero; ///< Aero payload; record() selects its checked concrete structure.
};
/// Checked StructuralQuadrature none payload; all fields retain one immutable case owner.
struct StructuralQuadratureNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    StructuralQuadratureNone(Field<detail::FieldType::T179> kind_value) noexcept
        : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked StructuralQuadrature active payload; all fields retain one immutable case owner.
struct StructuralQuadratureActive {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    StructuralQuadratureActive(Field<detail::FieldType::T189> kind_value,
                               Field<detail::FieldType::T219> surface_degree_value,
                               Field<detail::FieldType::T219> thickness_points_value,
                               Field<detail::FieldType::T219> independent_surface_degree_value,
                               Field<detail::FieldType::T219> independent_thickness_points_value,
                               Field<detail::FieldType::T88> absolute_value,
                               Field<detail::FieldType::T88> relative_value,
                               Field<detail::FieldType::T67> maximum_subdivisions_value) noexcept
        : kind(kind_value),
          surface_degree(surface_degree_value),
          thickness_points(thickness_points_value),
          independent_surface_degree(independent_surface_degree_value),
          independent_thickness_points(independent_thickness_points_value),
          absolute(absolute_value),
          relative(relative_value),
          maximum_subdivisions(maximum_subdivisions_value) {}
    Field<detail::FieldType::T189> kind;           ///< Kind; choose one of active.
    Field<detail::FieldType::T219> surface_degree; ///< Surface degree payload; record() selects its
                                                   ///< checked concrete structure.
    Field<detail::FieldType::T219> thickness_points; ///< Thickness points payload; record() selects
                                                     ///< its checked concrete structure.
    Field<detail::FieldType::T219>
        independent_surface_degree; ///< Independent surface degree payload; record() selects its
                                    ///< checked concrete structure.
    Field<detail::FieldType::T219>
        independent_thickness_points; ///< Independent thickness points payload; record() selects
                                      ///< its checked concrete structure.
    Field<detail::FieldType::T88>
        absolute; ///< Exact dimensionless quantity for absolute, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T88>
        relative; ///< Exact dimensionless quantity for relative, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T67>
        maximum_subdivisions; ///< Nonnegative count of maximum subdivisions.
};
/// Checked AeroQuadrature none payload; all fields retain one immutable case owner.
struct AeroQuadratureNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AeroQuadratureNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked AeroQuadrature active payload; all fields retain one immutable case owner.
struct AeroQuadratureActive {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    AeroQuadratureActive(Field<detail::FieldType::T303> kind_value,
                         Field<detail::FieldType::T88> absolute_value,
                         Field<detail::FieldType::T88> relative_value) noexcept
        : kind(kind_value), absolute(absolute_value), relative(relative_value) {}
    Field<detail::FieldType::T303> kind; ///< Kind; choose one of wp1_induction_v1.
    Field<detail::FieldType::T88>
        absolute; ///< Exact dimensionless quantity for absolute, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T88>
        relative; ///< Exact dimensionless quantity for relative, in declared case units; number()
                  ///< retains its dimensions.
};
/// Checked Scale uniform_components payload; all fields retain one immutable case owner.
struct ScaleUniformComponents {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ScaleUniformComponents(Field<detail::FieldType::T304> kind_value,
                           Field<detail::FieldType::T305> field_value,
                           Field<detail::FieldType::T277> targets_value,
                           Field<detail::FieldType::T171> values_value) noexcept
        : kind(kind_value), field(field_value), targets(targets_value), values(values_value) {}
    Field<detail::FieldType::T304> kind; ///< Kind; choose one of uniform_components.
    Field<detail::FieldType::T305>
        field; ///< Field; choose one of
               ///< shell_position,beam_translation,beam_rotation,rigid_translation,rigid_rotation,structural_velocity,structural_angular_velocity,circulation,wake_position,wake_strength,constraint,force_balance,moment_balance,aero_residual,transport_residual,newton_correction,arclength_state.
    Field<detail::FieldType::T277> targets; ///< Ordered targets entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T171> values;  ///< Ordered values entries; size() and at() expose the
                                            ///< checked sequence without copying.
};
/// Checked Scale physical_layout_policy payload; all fields retain one immutable case owner.
struct ScalePhysicalLayoutPolicy {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ScalePhysicalLayoutPolicy(Field<detail::FieldType::T306> kind_value,
                              Field<detail::FieldType::T307> field_value,
                              Field<detail::FieldType::T308> constraint_value) noexcept
        : kind(kind_value), field(field_value), constraint(constraint_value) {}
    Field<detail::FieldType::T306> kind;  ///< Kind; choose one of physical_layout_policy.
    Field<detail::FieldType::T307> field; ///< Field; choose one of mixed_blocks.
    Field<detail::FieldType::T308>
        constraint; ///< Constraint payload; record() selects its checked concrete structure.
};
/// Checked ConstraintScalePolicy none payload; all fields retain one immutable case owner.
struct ConstraintScalePolicyNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ConstraintScalePolicyNone(Field<detail::FieldType::T179> kind_value) noexcept
        : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked ConstraintScalePolicy physical payload; all fields retain one immutable case owner.
struct ConstraintScalePolicyPhysical {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ConstraintScalePolicyPhysical(Field<detail::FieldType::T309> kind_value,
                                  Field<detail::FieldType::T310> position_value,
                                  Field<detail::FieldType::T310> orientation_value,
                                  Field<detail::FieldType::T311> multiplier_work_scale_value,
                                  Field<detail::FieldType::T312> reaction_force_value,
                                  Field<detail::FieldType::T313> reaction_moment_value) noexcept
        : kind(kind_value),
          position(position_value),
          orientation(orientation_value),
          multiplier_work_scale(multiplier_work_scale_value),
          reaction_force(reaction_force_value),
          reaction_moment(reaction_moment_value) {}
    Field<detail::FieldType::T309> kind;     ///< Kind; choose one of physical.
    Field<detail::FieldType::T310> position; ///< Position-constraint scales in the case length
                                             ///< unit, or none when no position constraints exist.
    Field<detail::FieldType::T310>
        orientation; ///< Orientation-constraint scales in the case angle unit, or none when no
                     ///< orientation constraints exist.
    Field<detail::FieldType::T311>
        multiplier_work_scale; ///< Exact energy quantity for multiplier work scale, in declared
                               ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T312>
        reaction_force; ///< Ordered reaction force components; each numeric element exposes its
                        ///< physical dimensions and exact input value.
    Field<detail::FieldType::T313>
        reaction_moment; ///< Ordered reaction moment components; each numeric element exposes its
                         ///< physical dimensions and exact input value.
};
/// Checked ComponentScale3 none payload; all fields retain one immutable case owner.
struct ComponentScale3None {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ComponentScale3None(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked ComponentScale3 values payload; all fields retain one immutable case owner.
struct ComponentScale3Values {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ComponentScale3Values(Field<detail::FieldType::T314> kind_value,
                          Field<detail::FieldType::T31> values_value) noexcept
        : kind(kind_value), values(values_value) {}
    Field<detail::FieldType::T314> kind; ///< Kind; choose one of values.
    Field<detail::FieldType::T31>
        values; ///< Ordered values components; each numeric element exposes its physical dimensions
                ///< and exact input value.
};
/// Checked Scale explicit_layout payload; all fields retain one immutable case owner.
struct ScaleExplicitLayout {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ScaleExplicitLayout(Field<detail::FieldType::T315> kind_value,
                        Field<detail::FieldType::T305> field_value,
                        Field<detail::FieldType::T277> targets_value,
                        Field<detail::FieldType::T240> layout_value,
                        Field<detail::FieldType::T316> dimensions_value,
                        Field<detail::FieldType::T171> values_value) noexcept
        : kind(kind_value),
          field(field_value),
          targets(targets_value),
          layout(layout_value),
          dimensions(dimensions_value),
          values(values_value) {}
    Field<detail::FieldType::T315> kind; ///< Kind; choose one of explicit_layout.
    Field<detail::FieldType::T305>
        field; ///< Field; choose one of
               ///< shell_position,beam_translation,beam_rotation,rigid_translation,rigid_rotation,structural_velocity,structural_angular_velocity,circulation,wake_position,wake_strength,constraint,force_balance,moment_balance,aero_residual,transport_residual,newton_correction,arclength_state.
    Field<detail::FieldType::T277> targets; ///< Ordered targets entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T240> layout;  ///< Layout as decoded UTF-8 text; no normalization or
                                            ///< implicit replacement is performed.
    Field<detail::FieldType::T316>
        dimensions; ///< Length, mass, time and angle exponents, in that order.
    Field<detail::FieldType::T171> values; ///< Ordered values entries; size() and at() expose the
                                           ///< checked sequence without copying.
};
/// Checked LinearPolicy direct payload; all fields retain one immutable case owner.
struct LinearPolicyDirect {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LinearPolicyDirect(Field<detail::FieldType::T317> backend_value,
                       Field<detail::FieldType::T318> accuracy_value,
                       Field<detail::FieldType::T67> maximum_refinements_value,
                       Field<detail::FieldType::T319> backend_memory_value) noexcept
        : backend(backend_value),
          accuracy(accuracy_value),
          maximum_refinements(maximum_refinements_value),
          backend_memory(backend_memory_value) {}
    Field<detail::FieldType::T317>
        backend; ///< Backend; choose one of dense_lu,sparse_lu,sparse_llt.
    Field<detail::FieldType::T318>
        accuracy; ///< Accuracy payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T67>
        maximum_refinements;                       ///< Nonnegative count of maximum refinements.
    Field<detail::FieldType::T319> backend_memory; ///< Backend memory payload; record() selects its
                                                   ///< checked concrete structure.
};
/// Checked LinearPolicy gmres payload; all fields retain one immutable case owner.
struct LinearPolicyGmres {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LinearPolicyGmres(Field<detail::FieldType::T320> backend_value,
                      Field<detail::FieldType::T318> accuracy_value,
                      Field<detail::FieldType::T219> maximum_iterations_value,
                      Field<detail::FieldType::T219> restart_value,
                      Field<detail::FieldType::T67> maximum_applications_value,
                      Field<detail::FieldType::T67> maximum_adjoint_value,
                      Field<detail::FieldType::T67> maximum_preconditioner_value,
                      Field<detail::FieldType::T67> maximum_true_residual_checks_value,
                      Field<detail::FieldType::T319> backend_memory_value) noexcept
        : backend(backend_value),
          accuracy(accuracy_value),
          maximum_iterations(maximum_iterations_value),
          restart(restart_value),
          maximum_applications(maximum_applications_value),
          maximum_adjoint(maximum_adjoint_value),
          maximum_preconditioner(maximum_preconditioner_value),
          maximum_true_residual_checks(maximum_true_residual_checks_value),
          backend_memory(backend_memory_value) {}
    Field<detail::FieldType::T320> backend; ///< Backend; choose one of gmres.
    Field<detail::FieldType::T318>
        accuracy; ///< Accuracy payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T219> maximum_iterations; ///< Maximum iterations payload; record()
                                                       ///< selects its checked concrete structure.
    Field<detail::FieldType::T219>
        restart; ///< Restart payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T67>
        maximum_applications;                      ///< Nonnegative count of maximum applications.
    Field<detail::FieldType::T67> maximum_adjoint; ///< Nonnegative count of maximum adjoint.
    Field<detail::FieldType::T67>
        maximum_preconditioner; ///< Nonnegative count of maximum preconditioner.
    Field<detail::FieldType::T67>
        maximum_true_residual_checks; ///< Nonnegative count of maximum true residual checks.
    Field<detail::FieldType::T319> backend_memory; ///< Backend memory payload; record() selects its
                                                   ///< checked concrete structure.
};
/// Checked LinearAccuracy payload; all fields retain one immutable case owner.
struct LinearAccuracy {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LinearAccuracy(Field<detail::FieldType::T88> absolute_value,
                   Field<detail::FieldType::T88> relative_value,
                   Field<detail::FieldType::T88> backward_value,
                   Field<detail::FieldType::T321> norm_value) noexcept
        : absolute(absolute_value),
          relative(relative_value),
          backward(backward_value),
          norm(norm_value) {}
    Field<detail::FieldType::T88>
        absolute; ///< Exact dimensionless quantity for absolute, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T88>
        relative; ///< Exact dimensionless quantity for relative, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T88>
        backward; ///< Exact dimensionless quantity for backward, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T321> norm; ///< Norm; choose one of infinity,l2.
};
/// Checked NewtonPolicy payload; all fields retain one immutable case owner.
struct NewtonPolicy {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    NewtonPolicy(Field<detail::FieldType::T219> maximum_iterations_value,
                 Field<detail::FieldType::T219> maximum_line_searches_value,
                 Field<detail::FieldType::T88> absolute_value,
                 Field<detail::FieldType::T88> relative_value,
                 Field<detail::FieldType::T88> correction_value,
                 Field<detail::FieldType::T88> minimum_step_value) noexcept
        : maximum_iterations(maximum_iterations_value),
          maximum_line_searches(maximum_line_searches_value),
          absolute(absolute_value),
          relative(relative_value),
          correction(correction_value),
          minimum_step(minimum_step_value) {}
    Field<detail::FieldType::T219> maximum_iterations; ///< Maximum iterations payload; record()
                                                       ///< selects its checked concrete structure.
    Field<detail::FieldType::T219>
        maximum_line_searches; ///< Maximum line searches payload; record() selects its checked
                               ///< concrete structure.
    Field<detail::FieldType::T88>
        absolute; ///< Exact dimensionless quantity for absolute, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T88>
        relative; ///< Exact dimensionless quantity for relative, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T88>
        correction; ///< Exact dimensionless quantity for correction, in declared case units;
                    ///< number() retains its dimensions.
    Field<detail::FieldType::T88>
        minimum_step; ///< Exact dimensionless quantity for minimum step, in declared case units;
                      ///< number() retains its dimensions.
};
/// Checked PrecisionPolicy payload; all fields retain one immutable case owner.
struct PrecisionPolicy {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PrecisionPolicy(Field<detail::FieldType::T322> input_value,
                    Field<detail::FieldType::T323> retained_value,
                    Field<detail::FieldType::T324> returned_value,
                    Field<detail::FieldType::T88> projection_absolute_value,
                    Field<detail::FieldType::T88> projection_relative_value) noexcept
        : input(input_value),
          retained(retained_value),
          returned(returned_value),
          projection_absolute(projection_absolute_value),
          projection_relative(projection_relative_value) {}
    Field<detail::FieldType::T322> input; ///< Input; choose one of preserve.
    Field<detail::FieldType::T323>
        retained; ///< Retained; choose one of adaptive_binary64_two_word.
    Field<detail::FieldType::T324>
        returned; ///< Returned; choose one of retained,binary64_projection.
    Field<detail::FieldType::T88>
        projection_absolute; ///< Exact dimensionless quantity for projection absolute, in declared
                             ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        projection_relative; ///< Exact dimensionless quantity for projection relative, in declared
                             ///< case units; number() retains its dimensions.
};
/// Checked BackendMemory payload; all fields retain one immutable case owner.
struct BackendMemory {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BackendMemory(Field<detail::FieldType::T325> admission_value,
                  Field<detail::FieldType::T67> maximum_bytes_value) noexcept
        : admission(admission_value), maximum_bytes(maximum_bytes_value) {}
    Field<detail::FieldType::T325>
        admission; ///< Admission; choose one of estimate,require_proven_bound.
    Field<detail::FieldType::T67> maximum_bytes; ///< Nonnegative count of maximum bytes.
};
/// Checked EigenPolicy none payload; all fields retain one immutable case owner.
struct EigenPolicyNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenPolicyNone(Field<detail::FieldType::T179> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T179> kind; ///< Kind; choose one of none.
};
/// Checked EigenPolicy full payload; all fields retain one immutable case owner.
struct EigenPolicyFull {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenPolicyFull(Field<detail::FieldType::T326> kind_value,
                    Field<detail::FieldType::T327> target_value,
                    Field<detail::FieldType::T328> eigenvalue_scale_value,
                    Field<detail::FieldType::T329> packet_value,
                    Field<detail::FieldType::T330> requirements_value,
                    Field<detail::FieldType::T331> limits_value) noexcept
        : kind(kind_value),
          target(target_value),
          eigenvalue_scale(eigenvalue_scale_value),
          packet(packet_value),
          requirements(requirements_value),
          limits(limits_value) {}
    Field<detail::FieldType::T326> kind;             ///< Kind; choose one of full.
    Field<detail::FieldType::T327> target;           ///< Target; choose one of step_map.
    Field<detail::FieldType::T328> eigenvalue_scale; ///< Exact numeric eigenvalue scale; original
                                                     ///< decimal digits and exponent are retained.
    Field<detail::FieldType::T329> packet; ///< Packet; choose one of spectrum,pairs,subspace.
    Field<detail::FieldType::T330>
        requirements; ///< Requirements payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T331>
        limits; ///< Limits payload; record() selects its checked concrete structure.
};
/// Checked EigenPolicy generalized payload; all fields retain one immutable case owner.
struct EigenPolicyGeneralized {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenPolicyGeneralized(Field<detail::FieldType::T332> kind_value,
                           Field<detail::FieldType::T333> target_value,
                           Field<detail::FieldType::T334> eigenvalue_scale_value,
                           Field<detail::FieldType::T330> requirements_value,
                           Field<detail::FieldType::T331> limits_value) noexcept
        : kind(kind_value),
          target(target_value),
          eigenvalue_scale(eigenvalue_scale_value),
          requirements(requirements_value),
          limits(limits_value) {}
    Field<detail::FieldType::T332> kind;   ///< Kind; choose one of generalized.
    Field<detail::FieldType::T333> target; ///< Target; choose one of structural_modes.
    Field<detail::FieldType::T334>
        eigenvalue_scale; ///< Exact 1/T² quantity for eigenvalue scale, in declared case units;
                          ///< number() retains its dimensions.
    Field<detail::FieldType::T330>
        requirements; ///< Requirements payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T331>
        limits; ///< Limits payload; record() selects its checked concrete structure.
};
/// Checked EigenPolicy arnoldi payload; all fields retain one immutable case owner.
struct EigenPolicyArnoldi {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenPolicyArnoldi(Field<detail::FieldType::T335> kind_value,
                       Field<detail::FieldType::T327> target_value,
                       Field<detail::FieldType::T328> eigenvalue_scale_value,
                       Field<detail::FieldType::T336> packet_value,
                       Field<detail::FieldType::T337> selector_value,
                       Field<detail::FieldType::T219> nev_value,
                       Field<detail::FieldType::T219> ncv_value,
                       Field<detail::FieldType::T338> seed_value,
                       Field<detail::FieldType::T88> backend_tolerance_value,
                       Field<detail::FieldType::T330> requirements_value,
                       Field<detail::FieldType::T331> limits_value) noexcept
        : kind(kind_value),
          target(target_value),
          eigenvalue_scale(eigenvalue_scale_value),
          packet(packet_value),
          selector(selector_value),
          nev(nev_value),
          ncv(ncv_value),
          seed(seed_value),
          backend_tolerance(backend_tolerance_value),
          requirements(requirements_value),
          limits(limits_value) {}
    Field<detail::FieldType::T335> kind;             ///< Kind; choose one of arnoldi.
    Field<detail::FieldType::T327> target;           ///< Target; choose one of step_map.
    Field<detail::FieldType::T328> eigenvalue_scale; ///< Exact numeric eigenvalue scale; original
                                                     ///< decimal digits and exponent are retained.
    Field<detail::FieldType::T336> packet;           ///< Packet; choose one of pairs,subspace.
    Field<detail::FieldType::T337>
        selector; ///< Selector payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T219>
        nev; ///< Nev payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T219>
        ncv; ///< Ncv payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T338>
        seed; ///< Seed payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T88>
        backend_tolerance; ///< Exact dimensionless quantity for backend tolerance, in declared case
                           ///< units; number() retains its dimensions.
    Field<detail::FieldType::T330>
        requirements; ///< Requirements payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T331>
        limits; ///< Limits payload; record() selects its checked concrete structure.
};
/// Checked Complex payload; all fields retain one immutable case owner.
struct Complex {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Complex(Field<detail::FieldType::T328> real_value,
            Field<detail::FieldType::T328> imaginary_value) noexcept
        : real(real_value), imaginary(imaginary_value) {}
    Field<detail::FieldType::T328>
        real; ///< Exact numeric real; original decimal digits and exponent are retained.
    Field<detail::FieldType::T328>
        imaginary; ///< Exact numeric imaginary; original decimal digits and exponent are retained.
};
/// Checked ArnoldiSeed uniform payload; all fields retain one immutable case owner.
struct ArnoldiSeedUniform {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ArnoldiSeedUniform(Field<detail::FieldType::T339> kind_value,
                       Field<detail::FieldType::T340> value_value) noexcept
        : kind(kind_value), value(value_value) {}
    Field<detail::FieldType::T339> kind; ///< Kind; choose one of uniform_nonzero.
    Field<detail::FieldType::T340>
        value; ///< Value payload; record() selects its checked concrete structure.
};
/// Checked ArnoldiSeed explicit payload; all fields retain one immutable case owner.
struct ArnoldiSeedExplicit {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ArnoldiSeedExplicit(Field<detail::FieldType::T315> kind_value,
                        Field<detail::FieldType::T240> layout_value,
                        Field<detail::FieldType::T341> values_value) noexcept
        : kind(kind_value), layout(layout_value), values(values_value) {}
    Field<detail::FieldType::T315> kind;   ///< Kind; choose one of explicit_layout.
    Field<detail::FieldType::T240> layout; ///< Layout as decoded UTF-8 text; no normalization or
                                           ///< implicit replacement is performed.
    Field<detail::FieldType::T341> values; ///< Ordered values entries; size() and at() expose the
                                           ///< checked sequence without copying.
};
/// Checked EigenSelector largest payload; all fields retain one immutable case owner.
struct EigenSelectorLargest {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenSelectorLargest(Field<detail::FieldType::T342> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T342> kind; ///< Kind; choose one of largest_magnitude.
};
/// Checked EigenSelector shift payload; all fields retain one immutable case owner.
struct EigenSelectorShift {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenSelectorShift(Field<detail::FieldType::T343> kind_value,
                       Field<detail::FieldType::T340> shift_value) noexcept
        : kind(kind_value), shift(shift_value) {}
    Field<detail::FieldType::T343> kind; ///< Kind; choose one of shift_invert.
    Field<detail::FieldType::T340>
        shift; ///< Shift payload; record() selects its checked concrete structure.
};
/// Checked EigenRequirements payload; all fields retain one immutable case owner.
struct EigenRequirements {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenRequirements(Field<detail::FieldType::T88> residual_value,
                      Field<detail::FieldType::T344> residual_observation_value,
                      Field<detail::FieldType::T88> normalization_value,
                      Field<detail::FieldType::T344> normalization_observation_value,
                      Field<detail::FieldType::T88> orthogonality_value,
                      Field<detail::FieldType::T344> orthogonality_observation_value,
                      Field<detail::FieldType::T345> backward_value,
                      Field<detail::FieldType::T88> association_absolute_value,
                      Field<detail::FieldType::T88> association_relative_value,
                      Field<detail::FieldType::T135> left_value,
                      Field<detail::FieldType::T346> clusters_value,
                      Field<detail::FieldType::T135> physical_coverage_value) noexcept
        : residual(residual_value),
          residual_observation(residual_observation_value),
          normalization(normalization_value),
          normalization_observation(normalization_observation_value),
          orthogonality(orthogonality_value),
          orthogonality_observation(orthogonality_observation_value),
          backward(backward_value),
          association_absolute(association_absolute_value),
          association_relative(association_relative_value),
          left(left_value),
          clusters(clusters_value),
          physical_coverage(physical_coverage_value) {}
    Field<detail::FieldType::T88>
        residual; ///< Exact dimensionless quantity for residual, in declared case units; number()
                  ///< retains its dimensions.
    Field<detail::FieldType::T344>
        residual_observation; ///< Residual observation referring to a checked Observation
                              ///< definition in this case.
    Field<detail::FieldType::T88>
        normalization; ///< Exact dimensionless quantity for normalization, in declared case units;
                       ///< number() retains its dimensions.
    Field<detail::FieldType::T344>
        normalization_observation; ///< Normalization observation referring to a checked Observation
                                   ///< definition in this case.
    Field<detail::FieldType::T88>
        orthogonality; ///< Exact dimensionless quantity for orthogonality, in declared case units;
                       ///< number() retains its dimensions.
    Field<detail::FieldType::T344>
        orthogonality_observation; ///< Orthogonality observation referring to a checked Observation
                                   ///< definition in this case.
    Field<detail::FieldType::T345>
        backward; ///< Backward payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T88>
        association_absolute; ///< Exact dimensionless quantity for association absolute, in
                              ///< declared case units; number() retains its dimensions.
    Field<detail::FieldType::T88>
        association_relative; ///< Exact dimensionless quantity for association relative, in
                              ///< declared case units; number() retains its dimensions.
    Field<detail::FieldType::T135> left;     ///< Whether left is enabled.
    Field<detail::FieldType::T346> clusters; ///< Ordered clusters entries; size() and at() expose
                                             ///< the checked sequence without copying.
    Field<detail::FieldType::T135> physical_coverage; ///< Whether physical coverage is enabled.
};
/// Checked Backward none payload; all fields retain one immutable case owner.
struct BackwardNone {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BackwardNone(Field<detail::FieldType::T347> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T347> kind; ///< Kind; choose one of not_requested.
};
/// Checked Backward required payload; all fields retain one immutable case owner.
struct BackwardRequired {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    BackwardRequired(Field<detail::FieldType::T290> kind_value,
                     Field<detail::FieldType::T88> budget_value,
                     Field<detail::FieldType::T344> observation_value) noexcept
        : kind(kind_value), budget(budget_value), observation(observation_value) {}
    Field<detail::FieldType::T290> kind;  ///< Kind; choose one of required.
    Field<detail::FieldType::T88> budget; ///< Exact dimensionless quantity for budget, in declared
                                          ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T344>
        observation; ///< Observation referring to a checked Observation definition in this case.
};
/// Checked EigenCluster payload; all fields retain one immutable case owner.
struct EigenCluster {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenCluster(Field<detail::FieldType::T143> members_value,
                 Field<detail::FieldType::T219> expected_dimension_value) noexcept
        : members(members_value), expected_dimension(expected_dimension_value) {}
    Field<detail::FieldType::T143> members; ///< Ordered members entries; size() and at() expose the
                                            ///< checked sequence without copying.
    Field<detail::FieldType::T219> expected_dimension; ///< Expected dimension payload; record()
                                                       ///< selects its checked concrete structure.
};
/// Checked EigenLimits payload; all fields retain one immutable case owner.
struct EigenLimits {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    EigenLimits(Field<detail::FieldType::T67> maximum_order_value,
                Field<detail::FieldType::T67> maximum_results_value,
                Field<detail::FieldType::T67> maximum_ncv_value,
                Field<detail::FieldType::T67> maximum_backend_iterations_value,
                Field<detail::FieldType::T67> maximum_attempts_value,
                Field<detail::FieldType::T67> maximum_forward_value,
                Field<detail::FieldType::T67> maximum_adjoint_value,
                Field<detail::FieldType::T67> maximum_inverse_value,
                Field<detail::FieldType::T319> backend_memory_value) noexcept
        : maximum_order(maximum_order_value),
          maximum_results(maximum_results_value),
          maximum_ncv(maximum_ncv_value),
          maximum_backend_iterations(maximum_backend_iterations_value),
          maximum_attempts(maximum_attempts_value),
          maximum_forward(maximum_forward_value),
          maximum_adjoint(maximum_adjoint_value),
          maximum_inverse(maximum_inverse_value),
          backend_memory(backend_memory_value) {}
    Field<detail::FieldType::T67> maximum_order;   ///< Nonnegative count of maximum order.
    Field<detail::FieldType::T67> maximum_results; ///< Nonnegative count of maximum results.
    Field<detail::FieldType::T67> maximum_ncv;     ///< Nonnegative count of maximum ncv.
    Field<detail::FieldType::T67>
        maximum_backend_iterations; ///< Nonnegative count of maximum backend iterations.
    Field<detail::FieldType::T67> maximum_attempts; ///< Nonnegative count of maximum attempts.
    Field<detail::FieldType::T67> maximum_forward;  ///< Nonnegative count of maximum forward.
    Field<detail::FieldType::T67> maximum_adjoint;  ///< Nonnegative count of maximum adjoint.
    Field<detail::FieldType::T67> maximum_inverse;  ///< Nonnegative count of maximum inverse.
    Field<detail::FieldType::T319> backend_memory; ///< Backend memory payload; record() selects its
                                                   ///< checked concrete structure.
};
/// Checked Resources payload; all fields retain one immutable case owner.
struct Resources {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Resources(Field<detail::FieldType::T348> profile_value,
              Field<detail::FieldType::T219> threads_value,
              Field<detail::FieldType::T67> maximum_memory_bytes_value,
              Field<detail::FieldType::T44> maximum_wall_seconds_value,
              Field<detail::FieldType::T67> maximum_steps_value,
              Field<detail::FieldType::T67> maximum_total_iterations_value,
              Field<detail::FieldType::T67> checkpoint_temporary_bytes_value,
              Field<detail::FieldType::T349> input_value) noexcept
        : profile(profile_value),
          threads(threads_value),
          maximum_memory_bytes(maximum_memory_bytes_value),
          maximum_wall_seconds(maximum_wall_seconds_value),
          maximum_steps(maximum_steps_value),
          maximum_total_iterations(maximum_total_iterations_value),
          checkpoint_temporary_bytes(checkpoint_temporary_bytes_value),
          input(input_value) {}
    Field<detail::FieldType::T348> profile; ///< Profile; choose one of 24gib,48gib,custom.
    Field<detail::FieldType::T219>
        threads; ///< Threads payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T67>
        maximum_memory_bytes; ///< Nonnegative count of maximum memory bytes.
    Field<detail::FieldType::T44>
        maximum_wall_seconds; ///< Exact time quantity for maximum wall seconds, in declared case
                              ///< units; number() retains its dimensions.
    Field<detail::FieldType::T67> maximum_steps; ///< Nonnegative count of maximum steps.
    Field<detail::FieldType::T67>
        maximum_total_iterations; ///< Nonnegative count of maximum total iterations.
    Field<detail::FieldType::T67>
        checkpoint_temporary_bytes; ///< Nonnegative count of checkpoint temporary bytes.
    Field<detail::FieldType::T349>
        input; ///< Input payload; record() selects its checked concrete structure.
};
/// Checked InputQuotas payload; all fields retain one immutable case owner.
struct InputQuotas {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InputQuotas(Field<detail::FieldType::T67> maximum_source_bytes_value,
                Field<detail::FieldType::T67> maximum_depth_value,
                Field<detail::FieldType::T67> maximum_nodes_value,
                Field<detail::FieldType::T67> maximum_members_value,
                Field<detail::FieldType::T67> maximum_decoded_string_bytes_value,
                Field<detail::FieldType::T67> maximum_decimal_token_bytes_value,
                Field<detail::FieldType::T67> maximum_total_text_bytes_value,
                Field<detail::FieldType::T67> maximum_typed_bytes_value,
                Field<detail::FieldType::T67> maximum_entities_value,
                Field<detail::FieldType::T67> maximum_references_value,
                Field<detail::FieldType::T67> maximum_scalars_value,
                Field<detail::FieldType::T67> maximum_polynomial_terms_value,
                Field<detail::FieldType::T67> maximum_table_samples_value,
                Field<detail::FieldType::T67> maximum_expression_nodes_value,
                Field<detail::FieldType::T67> maximum_expression_edges_value,
                Field<detail::FieldType::T67> maximum_expression_depth_value,
                Field<detail::FieldType::T67> maximum_expression_power_value,
                Field<detail::FieldType::T67> maximum_expression_points_value,
                Field<detail::FieldType::T67> maximum_diagnostic_path_bytes_value) noexcept
        : maximum_source_bytes(maximum_source_bytes_value),
          maximum_depth(maximum_depth_value),
          maximum_nodes(maximum_nodes_value),
          maximum_members(maximum_members_value),
          maximum_decoded_string_bytes(maximum_decoded_string_bytes_value),
          maximum_decimal_token_bytes(maximum_decimal_token_bytes_value),
          maximum_total_text_bytes(maximum_total_text_bytes_value),
          maximum_typed_bytes(maximum_typed_bytes_value),
          maximum_entities(maximum_entities_value),
          maximum_references(maximum_references_value),
          maximum_scalars(maximum_scalars_value),
          maximum_polynomial_terms(maximum_polynomial_terms_value),
          maximum_table_samples(maximum_table_samples_value),
          maximum_expression_nodes(maximum_expression_nodes_value),
          maximum_expression_edges(maximum_expression_edges_value),
          maximum_expression_depth(maximum_expression_depth_value),
          maximum_expression_power(maximum_expression_power_value),
          maximum_expression_points(maximum_expression_points_value),
          maximum_diagnostic_path_bytes(maximum_diagnostic_path_bytes_value) {}
    Field<detail::FieldType::T67>
        maximum_source_bytes;                      ///< Nonnegative count of maximum source bytes.
    Field<detail::FieldType::T67> maximum_depth;   ///< Nonnegative count of maximum depth.
    Field<detail::FieldType::T67> maximum_nodes;   ///< Nonnegative count of maximum nodes.
    Field<detail::FieldType::T67> maximum_members; ///< Nonnegative count of maximum members.
    Field<detail::FieldType::T67>
        maximum_decoded_string_bytes; ///< Nonnegative count of maximum decoded string bytes.
    Field<detail::FieldType::T67>
        maximum_decimal_token_bytes; ///< Nonnegative count of maximum decimal token bytes.
    Field<detail::FieldType::T67>
        maximum_total_text_bytes; ///< Nonnegative count of maximum total text bytes.
    Field<detail::FieldType::T67>
        maximum_typed_bytes;                          ///< Nonnegative count of maximum typed bytes.
    Field<detail::FieldType::T67> maximum_entities;   ///< Nonnegative count of maximum entities.
    Field<detail::FieldType::T67> maximum_references; ///< Nonnegative count of maximum references.
    Field<detail::FieldType::T67> maximum_scalars;    ///< Nonnegative count of maximum scalars.
    Field<detail::FieldType::T67>
        maximum_polynomial_terms; ///< Nonnegative count of maximum polynomial terms.
    Field<detail::FieldType::T67>
        maximum_table_samples; ///< Nonnegative count of maximum table samples.
    Field<detail::FieldType::T67>
        maximum_expression_nodes; ///< Nonnegative count of maximum expression nodes.
    Field<detail::FieldType::T67>
        maximum_expression_edges; ///< Nonnegative count of maximum expression edges.
    Field<detail::FieldType::T67>
        maximum_expression_depth; ///< Nonnegative count of maximum expression depth.
    Field<detail::FieldType::T67>
        maximum_expression_power; ///< Nonnegative count of maximum expression power.
    Field<detail::FieldType::T67>
        maximum_expression_points; ///< Nonnegative count of maximum expression points.
    Field<detail::FieldType::T67>
        maximum_diagnostic_path_bytes; ///< Nonnegative count of maximum diagnostic path bytes.
};
/// Checked Output payload; all fields retain one immutable case owner.
struct Output {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Output(Field<detail::FieldType::T24> directory_value,
           Field<detail::FieldType::T350> observations_value,
           Field<detail::FieldType::T351> visualization_value,
           Field<detail::FieldType::T352> checkpoint_value,
           Field<detail::FieldType::T353> overwrite_value,
           Field<detail::FieldType::T354> resume_value) noexcept
        : directory(directory_value),
          observations(observations_value),
          visualization(visualization_value),
          checkpoint(checkpoint_value),
          overwrite(overwrite_value),
          resume(resume_value) {}
    Field<detail::FieldType::T24> directory; ///< Directory as decoded UTF-8 text; no normalization
                                             ///< or implicit replacement is performed.
    Field<detail::FieldType::T350>
        observations; ///< Observations payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T351>
        visualization; ///< Visualization payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T352>
        checkpoint; ///< Checkpoint payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T353> overwrite; ///< Overwrite; choose one of reject_existing.
    Field<detail::FieldType::T354> resume;    ///< Resume; choose one of explicit_only.
};
/// Checked OutputCadence payload; all fields retain one immutable case owner.
struct OutputCadence {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    OutputCadence(Field<detail::FieldType::T219> every_steps_value,
                  Field<detail::FieldType::T89> observations_value) noexcept
        : every_steps(every_steps_value), observations(observations_value) {}
    Field<detail::FieldType::T219>
        every_steps; ///< Every steps payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T89> observations; ///< Ordered observations entries; size() and at()
                                                ///< expose the checked sequence without copying.
};
/// Checked Visualization payload; all fields retain one immutable case owner.
struct Visualization {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Visualization(Field<detail::FieldType::T135> enabled_value,
                  Field<detail::FieldType::T219> every_steps_value,
                  Field<detail::FieldType::T355> layers_value,
                  Field<detail::FieldType::T356> interfaces_value,
                  Field<detail::FieldType::T135> reference_geometry_value,
                  Field<detail::FieldType::T135> current_geometry_value) noexcept
        : enabled(enabled_value),
          every_steps(every_steps_value),
          layers(layers_value),
          interfaces(interfaces_value),
          reference_geometry(reference_geometry_value),
          current_geometry(current_geometry_value) {}
    Field<detail::FieldType::T135> enabled; ///< Whether enabled is enabled.
    Field<detail::FieldType::T219>
        every_steps; ///< Every steps payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T355> layers; ///< Ordered layers entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T356> interfaces;         ///< Interfaces; choose one of both.
    Field<detail::FieldType::T135> reference_geometry; ///< Whether reference geometry is enabled.
    Field<detail::FieldType::T135> current_geometry;   ///< Whether current geometry is enabled.
};
/// Checked CheckpointOutput off payload; all fields retain one immutable case owner.
struct CheckpointOutputOff {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    CheckpointOutputOff(Field<detail::FieldType::T357> kind_value) noexcept : kind(kind_value) {}
    Field<detail::FieldType::T357> kind; ///< Kind; choose one of off.
};
/// Checked CheckpointOutput on payload; all fields retain one immutable case owner.
struct CheckpointOutputOn {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    CheckpointOutputOn(Field<detail::FieldType::T358> kind_value,
                       Field<detail::FieldType::T219> every_steps_value,
                       Field<detail::FieldType::T24> directory_value,
                       Field<detail::FieldType::T219> retain_value) noexcept
        : kind(kind_value),
          every_steps(every_steps_value),
          directory(directory_value),
          retain(retain_value) {}
    Field<detail::FieldType::T358> kind; ///< Kind; choose one of on.
    Field<detail::FieldType::T219>
        every_steps; ///< Every steps payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T24> directory; ///< Directory as decoded UTF-8 text; no normalization
                                             ///< or implicit replacement is performed.
    Field<detail::FieldType::T219>
        retain; ///< Retain payload; record() selects its checked concrete structure.
};
/// Checked ReferenceField<d,n> expression payload; all fields retain one immutable case owner.
struct ReferenceFieldExpression {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ReferenceFieldExpression(Field<detail::FieldType::T359> kind_value,
                             Field<detail::FieldType::T360> expression_value,
                             Field<detail::FieldType::T39> error_value) noexcept
        : kind(kind_value), expression(expression_value), error(error_value) {}
    Field<detail::FieldType::T359> kind; ///< Kind; choose one of expression.
    Field<detail::FieldType::T360>
        expression; ///< Expression payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T39> error; ///< Ordered error components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
};
/// Checked TimeSeries<n> expression payload; all fields retain one immutable case owner.
struct TimeSeriesExpression {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    TimeSeriesExpression(Field<detail::FieldType::T359> kind_value,
                         Field<detail::FieldType::T360> expression_value) noexcept
        : kind(kind_value), expression(expression_value) {}
    Field<detail::FieldType::T359> kind; ///< Kind; choose one of expression.
    Field<detail::FieldType::T360>
        expression; ///< Expression payload; record() selects its checked concrete structure.
};
/// Checked Load current_expression payload; all fields retain one immutable case owner.
struct LoadCurrentExpression {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    LoadCurrentExpression(Field<detail::FieldType::T361> kind_value,
                          Field<detail::FieldType::T157> id_value,
                          Field<detail::FieldType::T133> region_value,
                          Field<detail::FieldType::T362> expression_value) noexcept
        : kind(kind_value), id(id_value), region(region_value), expression(expression_value) {}
    Field<detail::FieldType::T361> kind; ///< Kind; choose one of current_pressure.
    Field<detail::FieldType::T157>
        id; ///< Id in the Load identifier namespace; equal IDs in other namespaces remain distinct.
    Field<detail::FieldType::T133>
        region; ///< Region payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T362>
        expression; ///< Expression payload; record() selects its checked concrete structure.
};
/// Checked Expression<n> payload; all fields retain one immutable case owner.
struct Expression {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Expression(Field<detail::FieldType::T363> domain_value,
               Field<detail::FieldType::T364> nodes_value,
               Field<detail::FieldType::T365> outputs_value) noexcept
        : domain(domain_value), nodes(nodes_value), outputs(outputs_value) {}
    Field<detail::FieldType::T363>
        domain; ///< Domain payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T364> nodes; ///< Ordered nodes entries; size() and at() expose the
                                          ///< checked sequence without copying.
    Field<detail::FieldType::T365>
        outputs; ///< Ordered result node indices; their count equals the consumer component count.
};
/// Checked ExpressionDomain reference payload; all fields retain one immutable case owner.
struct ExpressionDomainReference {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExpressionDomainReference(Field<detail::FieldType::T213> kind_value,
                              Field<detail::FieldType::T41> coordinates_value,
                              Field<detail::FieldType::T31> lower_value,
                              Field<detail::FieldType::T31> upper_value) noexcept
        : kind(kind_value),
          coordinates(coordinates_value),
          lower(lower_value),
          upper(upper_value) {}
    Field<detail::FieldType::T213> kind; ///< Kind; choose one of reference.
    Field<detail::FieldType::T41>
        coordinates; ///< Coordinate system used by this field; its choices are closed.
    Field<detail::FieldType::T31> lower; ///< Ordered lower components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
    Field<detail::FieldType::T31> upper; ///< Ordered upper components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
};
/// Checked ExpressionDomain time payload; all fields retain one immutable case owner.
struct ExpressionDomainTime {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExpressionDomainTime(Field<detail::FieldType::T366> kind_value,
                         Field<detail::FieldType::T44> lower_value,
                         Field<detail::FieldType::T44> upper_value,
                         Field<detail::FieldType::T48> events_value) noexcept
        : kind(kind_value), lower(lower_value), upper(upper_value), events(events_value) {}
    Field<detail::FieldType::T366> kind; ///< Kind; choose one of time.
    Field<detail::FieldType::T44> lower; ///< Exact time quantity for lower, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T44> upper; ///< Exact time quantity for upper, in declared case units;
                                         ///< number() retains its dimensions.
    Field<detail::FieldType::T48> events; ///< Ordered events entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked ExpressionDomain current_position_time payload; all fields retain one immutable case
/// owner.
struct ExpressionDomainCurrentPositionTime {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExpressionDomainCurrentPositionTime(Field<detail::FieldType::T367> kind_value,
                                        Field<detail::FieldType::T26> position_lower_value,
                                        Field<detail::FieldType::T26> position_upper_value,
                                        Field<detail::FieldType::T44> time_lower_value,
                                        Field<detail::FieldType::T44> time_upper_value,
                                        Field<detail::FieldType::T48> events_value) noexcept
        : kind(kind_value),
          position_lower(position_lower_value),
          position_upper(position_upper_value),
          time_lower(time_lower_value),
          time_upper(time_upper_value),
          events(events_value) {}
    Field<detail::FieldType::T367> kind; ///< Kind; choose one of current_position_time.
    Field<detail::FieldType::T26>
        position_lower; ///< Ordered position lower components; each numeric element exposes its
                        ///< physical dimensions and exact input value.
    Field<detail::FieldType::T26>
        position_upper; ///< Ordered position upper components; each numeric element exposes its
                        ///< physical dimensions and exact input value.
    Field<detail::FieldType::T44> time_lower; ///< Exact time quantity for time lower, in declared
                                              ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T44> time_upper; ///< Exact time quantity for time upper, in declared
                                              ///< case units; number() retains its dimensions.
    Field<detail::FieldType::T48> events; ///< Ordered events entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked ExprNode constant payload; all fields retain one immutable case owner.
struct ExprNodeConstant {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExprNodeConstant(Field<detail::FieldType::T38> op_value,
                     Field<detail::FieldType::T328> value_value,
                     Field<detail::FieldType::T368> dimensions_value) noexcept
        : op(op_value), value(value_value), dimensions(dimensions_value) {}
    Field<detail::FieldType::T38> op; ///< Op; choose one of constant.
    Field<detail::FieldType::T328>
        value; ///< Exact numeric value; original decimal digits and exponent are retained.
    Field<detail::FieldType::T368>
        dimensions; ///< Length, mass, time and angle exponents, in that order.
};
/// Checked ExprNode variable payload; all fields retain one immutable case owner.
struct ExprNodeVariable {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExprNodeVariable(Field<detail::FieldType::T369> op_value,
                     Field<detail::FieldType::T35> axis_value) noexcept
        : op(op_value), axis(axis_value) {}
    Field<detail::FieldType::T369> op;  ///< Op; choose one of variable.
    Field<detail::FieldType::T35> axis; ///< Zero-based index of axis.
};
/// Checked ExprNode binary payload; all fields retain one immutable case owner.
struct ExprNodeBinary {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExprNodeBinary(Field<detail::FieldType::T370> op_value,
                   Field<detail::FieldType::T35> first_value,
                   Field<detail::FieldType::T35> second_value) noexcept
        : op(op_value), first(first_value), second(second_value) {}
    Field<detail::FieldType::T370> op;    ///< Op; choose one of add,subtract,multiply,divide.
    Field<detail::FieldType::T35> first;  ///< First earlier operand node index.
    Field<detail::FieldType::T35> second; ///< Second earlier operand node index.
};
/// Checked ExprNode negate payload; all fields retain one immutable case owner.
struct ExprNodeNegate {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExprNodeNegate(Field<detail::FieldType::T371> op_value,
                   Field<detail::FieldType::T35> operand_value) noexcept
        : op(op_value), operand(operand_value) {}
    Field<detail::FieldType::T371> op; ///< Op; choose one of negate.
    Field<detail::FieldType::T35>
        operand; ///< Index of an earlier expression node; forward references are rejected.
};
/// Checked ExprNode power payload; all fields retain one immutable case owner.
struct ExprNodePower {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExprNodePower(Field<detail::FieldType::T372> op_value,
                  Field<detail::FieldType::T35> operand_value,
                  Field<detail::FieldType::T373> exponent_value) noexcept
        : op(op_value), operand(operand_value), exponent(exponent_value) {}
    Field<detail::FieldType::T372> op; ///< Op; choose one of integer_power.
    Field<detail::FieldType::T35>
        operand; ///< Index of an earlier expression node; forward references are rejected.
    Field<detail::FieldType::T373>
        exponent; ///< Exponent payload; record() selects its checked concrete structure.
};
/// Checked ExprNode analytic payload; all fields retain one immutable case owner.
struct ExprNodeAnalytic {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ExprNodeAnalytic(Field<detail::FieldType::T374> op_value,
                     Field<detail::FieldType::T35> operand_value) noexcept
        : op(op_value), operand(operand_value) {}
    Field<detail::FieldType::T374> op; ///< Op; choose one of sin,cos,exp,log,sqrt.
    Field<detail::FieldType::T35>
        operand; ///< Index of an earlier expression node; forward references are rejected.
};
/// Checked Perturbation physical payload; all fields retain one immutable case owner.
struct PerturbationPhysical {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PerturbationPhysical(Field<detail::FieldType::T309> kind_value,
                         Field<detail::FieldType::T375> shells_value,
                         Field<detail::FieldType::T238> beams_value,
                         Field<detail::FieldType::T239> rigids_value,
                         Field<detail::FieldType::T376> displacements_value,
                         Field<detail::FieldType::T376> velocities_value) noexcept
        : kind(kind_value),
          shells(shells_value),
          beams(beams_value),
          rigids(rigids_value),
          displacements(displacements_value),
          velocities(velocities_value) {}
    Field<detail::FieldType::T309> kind;   ///< Kind; choose one of physical.
    Field<detail::FieldType::T375> shells; ///< Ordered shells entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T238> beams;  ///< Ordered beams entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T239> rigids; ///< Ordered rigids entries; size() and at() expose the
                                           ///< checked sequence without copying.
    Field<detail::FieldType::T376>
        displacements; ///< Ordered displacements entries; size() and at() expose the checked
                       ///< sequence without copying.
    Field<detail::FieldType::T376> velocities; ///< Ordered velocities entries; size() and at()
                                               ///< expose the checked sequence without copying.
};
/// Checked ShellPhysicalState payload; all fields retain one immutable case owner.
struct ShellPhysicalState {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    ShellPhysicalState(Field<detail::FieldType::T177> shell_value,
                       Field<detail::FieldType::T377> displacement_value,
                       Field<detail::FieldType::T377> velocity_value) noexcept
        : shell(shell_value), displacement(displacement_value), velocity(velocity_value) {}
    Field<detail::FieldType::T177>
        shell; ///< Shell referring to a checked Shell definition in this case.
    Field<detail::FieldType::T377>
        displacement; ///< Displacement payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T377>
        velocity; ///< Velocity payload; record() selects its checked concrete structure.
};
/// Checked FullShellField global_polynomial payload; all fields retain one immutable case owner.
struct FullShellFieldGlobalPolynomial {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    FullShellFieldGlobalPolynomial(Field<detail::FieldType::T378> kind_value,
                                   Field<detail::FieldType::T379> r_value,
                                   Field<detail::FieldType::T379> p_value,
                                   Field<detail::FieldType::T379> c_value) noexcept
        : kind(kind_value), r(r_value), p(p_value), c(c_value) {}
    Field<detail::FieldType::T378> kind; ///< Kind; choose one of global_polynomial.
    Field<detail::FieldType::T379>
        r; ///< R payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T379>
        p; ///< P payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T379>
        c; ///< C payload; record() selects its checked concrete structure.
};
/// Checked FullShellField half_cells payload; all fields retain one immutable case owner.
struct FullShellFieldHalfCells {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    FullShellFieldHalfCells(Field<detail::FieldType::T380> kind_value,
                            Field<detail::FieldType::T381> cells_value) noexcept
        : kind(kind_value), cells(cells_value) {}
    Field<detail::FieldType::T380> kind;  ///< Kind; choose one of half_cells.
    Field<detail::FieldType::T381> cells; ///< Ordered cells entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked CellField payload; all fields retain one immutable case owner.
struct CellField {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    CellField(Field<detail::FieldType::T136> cell_value, Field<detail::FieldType::T379> r_value,
              Field<detail::FieldType::T379> p_value,
              Field<detail::FieldType::T379> c_value) noexcept
        : cell(cell_value), r(r_value), p(p_value), c(c_value) {}
    Field<detail::FieldType::T136>
        cell; ///< Two half-cell indices along the shell reference coordinates.
    Field<detail::FieldType::T379>
        r; ///< R payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T379>
        p; ///< P payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T379>
        c; ///< C payload; record() selects its checked concrete structure.
};
/// Checked VectorPolynomial payload; all fields retain one immutable case owner.
struct VectorPolynomial {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    VectorPolynomial(Field<detail::FieldType::T382> terms_value) noexcept : terms(terms_value) {}
    Field<detail::FieldType::T382> terms; ///< Ordered terms entries; size() and at() expose the
                                          ///< checked sequence without copying.
};
/// Checked VectorTerm payload; all fields retain one immutable case owner.
struct VectorTerm {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    VectorTerm(Field<detail::FieldType::T63> powers_value,
               Field<detail::FieldType::T383> value_value) noexcept
        : powers(powers_value), value(value_value) {}
    Field<detail::FieldType::T63> powers; ///< Nonnegative monomial exponents, in coordinate order.
    Field<detail::FieldType::T383>
        value; ///< Ordered value components; each numeric element exposes its physical dimensions
               ///< and exact input value.
};
/// Checked PointCheck payload; all fields retain one immutable case owner.
struct PointCheck {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    PointCheck(Field<detail::FieldType::T384> point_value,
               Field<detail::FieldType::T31> value_value,
               Field<detail::FieldType::T31> tolerance_value) noexcept
        : point(point_value), value(value_value), tolerance(tolerance_value) {}
    Field<detail::FieldType::T384>
        point; ///< Point payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T31> value; ///< Ordered value components; each numeric element exposes
                                         ///< its physical dimensions and exact input value.
    Field<detail::FieldType::T31>
        tolerance; ///< Ordered tolerance components; each numeric element exposes its physical
                   ///< dimensions and exact input value.
};
/// Checked InitialPoint shell payload; all fields retain one immutable case owner.
struct InitialPointShell {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialPointShell(Field<detail::FieldType::T385> target_value,
                      Field<detail::FieldType::T31> coordinates_value,
                      Field<detail::FieldType::T32> side_s_value,
                      Field<detail::FieldType::T32> side_t_value) noexcept
        : target(target_value),
          coordinates(coordinates_value),
          side_s(side_s_value),
          side_t(side_t_value) {}
    Field<detail::FieldType::T385>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T31> coordinates; ///< Shell reference coordinates s, t and thickness
                                               ///< coordinate zeta, all dimensionless.
    Field<detail::FieldType::T32>
        side_s; ///< One-sided trace selection along the first shell reference coordinate.
    Field<detail::FieldType::T32>
        side_t; ///< One-sided trace selection along the second shell reference coordinate.
};
/// Checked InitialPoint beam payload; all fields retain one immutable case owner.
struct InitialPointBeam {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialPointBeam(Field<detail::FieldType::T386> target_value,
                     Field<detail::FieldType::T31> coordinates_value) noexcept
        : target(target_value), coordinates(coordinates_value) {}
    Field<detail::FieldType::T386>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T31> coordinates; ///< Beam material coordinate in the case length unit
                                               ///< followed by two exact zeros.
};
/// Checked InitialPoint rigid payload; all fields retain one immutable case owner.
struct InitialPointRigid {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    InitialPointRigid(Field<detail::FieldType::T387> target_value,
                      Field<detail::FieldType::T31> coordinates_value) noexcept
        : target(target_value), coordinates(coordinates_value) {}
    Field<detail::FieldType::T387>
        target; ///< Target payload; record() selects its checked concrete structure.
    Field<detail::FieldType::T31>
        coordinates; ///< Rigid-body reference coordinates, ordered x, y, z in the case length unit.
};
/// Checked Rational payload; all fields retain one immutable case owner.
struct Rational {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    Rational(Field<detail::FieldType::T24> numerator_value,
             Field<detail::FieldType::T24> denominator_value) noexcept
        : numerator(numerator_value), denominator(denominator_value) {}
    Field<detail::FieldType::T24>
        numerator; ///< Exact signed decimal integer numerator; no binary floating conversion.
    Field<detail::FieldType::T24>
        denominator; ///< Exact positive decimal integer denominator; fractions need not be reduced.
};
/// Checked FileReference payload; all fields retain one immutable case owner.
struct FileReference {
    /// Bind every declared field explicitly; no missing value is silently defaulted.
    FileReference(Field<detail::FieldType::T388> source_id_value,
                  Field<detail::FieldType::T24> path_value,
                  Field<detail::FieldType::T240> sha256_value,
                  Field<detail::FieldType::T446> role_value) noexcept
        : source_id(source_id_value), path(path_value), sha256(sha256_value), role(role_value) {}
    Field<detail::FieldType::T388> source_id; ///< Source identity in the source namespace; repeated
                                              ///< IDs require identical path, digest and role.
    Field<detail::FieldType::T24>
        path; ///< Decoded UTF-8 relative source path; file access belongs to preparation.
    Field<detail::FieldType::T240>
        sha256; ///< Lowercase hexadecimal SHA-256 of the referenced source bytes.
    Field<detail::FieldType::T446>
        role; ///< Closed external-file purpose; it must match the referring consumer.
};
/// Checked MaterialMap payload; all fields retain one immutable case owner.
struct MaterialMap {
    explicit MaterialMap(Field<detail::FieldType::T448> pieces_value) noexcept
        : pieces(pieces_value) {}
    Field<detail::FieldType::T448> pieces;
};
/// Checked MaterialMapPiece payload; all fields retain one immutable case owner.
struct MaterialMapPiece {
    MaterialMapPiece(Field<detail::FieldType::T136> chord_cells_value,
                     Field<detail::FieldType::T136> span_cells_value,
                     Field<detail::FieldType::T265> layer_value,
                     Field<detail::FieldType::T32> side_value,
                     Field<detail::FieldType::T450> terms_value) noexcept
        : chord_cells(chord_cells_value),
          span_cells(span_cells_value),
          layer(layer_value),
          side(side_value),
          terms(terms_value) {}
    Field<detail::FieldType::T136> chord_cells;
    Field<detail::FieldType::T136> span_cells;
    Field<detail::FieldType::T265> layer;
    Field<detail::FieldType::T32> side;
    Field<detail::FieldType::T450> terms;
};
/// Checked MaterialMapTerm payload; all fields retain one immutable case owner.
struct MaterialMapTerm {
    MaterialMapTerm(Field<detail::FieldType::T63> powers_value,
                    Field<detail::FieldType::T452> coefficient_value) noexcept
        : powers(powers_value), coefficient(coefficient_value) {}
    Field<detail::FieldType::T63> powers;
    Field<detail::FieldType::T452> coefficient;
};
} // namespace schema
namespace detail {
template <>
struct RecordAccess<FieldType::T1> {
    using Value = schema::Units;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T2> {
    using Value = schema::Model;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T3> {
    using Value = std::variant<schema::FlowNone, schema::FlowKnown>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T4> {
    using Value =
        std::variant<schema::AnalysisStatic, schema::AnalysisTransient, schema::AnalysisFlutter>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T6> {
    using Value = schema::Numerics;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T7> {
    using Value = schema::Resources;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T8> {
    using Value = schema::Output;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T9> {
    using Value =
        std::variant<schema::InitialEquilibrium, schema::InitialCold, schema::InitialPrepared>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T10> {
    using Value = schema::Restart;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T11> {
    using Value = schema::Provenance;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T30> {
    using Value = schema::Target;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T49> {
    using Value = std::variant<schema::OutsideRejectHold, schema::OutsidePeriodic>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T52> {
    using Value = schema::Frame;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T53> {
    using Value = std::variant<schema::TimeSeriesConstant, schema::TimeSeriesPolynomial,
                               schema::TimeSeriesPchip, schema::TimeSeriesExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T55> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T56> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T58> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T62> {
    using Value = schema::Rational;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T66> {
    using Value = schema::Thickness;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T77> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T78> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T82> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T83> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T84> {
    using Value = std::variant<schema::ReferenceGeometrySupplied, schema::ReferenceGeometryFit>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T86> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T94> {
    using Value = std::variant<schema::ApplicabilityUnevaluated, schema::ApplicabilityDeclared>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T107> {
    using Value = std::variant<schema::AxisProjected, schema::AxisExplicit>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T108> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T115> {
    using Value = std::variant<schema::BeamCoordinateArclength, schema::BeamCoordinateScaled>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T116> {
    using Value =
        std::variant<schema::SectionMatrix, schema::SectionCircle, schema::SectionRectangle>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T133> {
    using Value = schema::Region;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T142> {
    using Value = schema::MaterialPoint;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T144> {
    using Value = schema::Motion;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T150> {
    using Value = std::variant<schema::LedgerAdd, schema::LedgerReplace>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T163> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T164> {
    using Value = std::variant<schema::TimeSeriesConstant, schema::TimeSeriesPolynomial,
                               schema::TimeSeriesPchip, schema::TimeSeriesExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T166> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T168> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T170> {
    using Value = std::variant<schema::AeroAttachmentFixed, schema::AeroAttachmentPrescribed,
                               schema::AeroAttachmentStructural>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T181> {
    using Value = std::variant<schema::GustNone, schema::GustSine>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T182> {
    using Value = schema::FlowReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T190> {
    using Value = std::variant<schema::ConvectionFree, schema::ConvectionPrescribed>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T191> {
    using Value = std::variant<schema::RetentionAll, schema::RetentionWindow>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T199> {
    using Value = schema::Branch;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T200> {
    using Value = std::variant<schema::AnalysisSpeedNone, schema::AnalysisSpeedValue>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T201> {
    using Value = std::variant<schema::WakeNone, schema::WakeActive>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T208> {
    using Value = schema::Search;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T211> {
    using Value = std::variant<schema::BranchSeedReference, schema::BranchSeedPrepared>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T212> {
    using Value = std::variant<schema::ContinuationFixed, schema::ContinuationSweep>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T215> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T218> {
    using Value = schema::ArcLength;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T225> {
    using Value = schema::Association;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T226> {
    using Value = schema::Coverage;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T230> {
    using Value = std::variant<schema::PerturbationNone, schema::PerturbationPhysical>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T232> {
    using Value = schema::InitialState;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T234> {
    using Value = schema::InitialWake;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T236> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T242> {
    using Value = std::variant<schema::PlanesBinary64, schema::PlanesTwoWord>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T248> {
    using Value = schema::FrameStatePayload;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T261> {
    using Value = std::variant<schema::BasisWorldMaterial, schema::BasisExplicit>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T262> {
    using Value = schema::Comparison;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T263> {
    using Value = schema::ObservationBudget;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T283> {
    using Value = std::variant<schema::ComparisonReferenceNone, schema::ComparisonReferenceFile>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T285> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T292> {
    using Value = std::variant<schema::PatchPolicyNone, schema::PatchPolicyActive>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T293> {
    using Value = schema::Quadrature;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T294> {
    using Value = std::variant<schema::LinearPolicyDirect, schema::LinearPolicyGmres>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T295> {
    using Value = schema::NewtonPolicy;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T296> {
    using Value = schema::PrecisionPolicy;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T297> {
    using Value = std::variant<schema::EigenPolicyNone, schema::EigenPolicyFull,
                               schema::EigenPolicyGeneralized, schema::EigenPolicyArnoldi>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T298> {
    using Value = std::variant<schema::WakeCapacityNone, schema::WakeCapacityLimited>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T301> {
    using Value =
        std::variant<schema::StructuralQuadratureNone, schema::StructuralQuadratureActive>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T302> {
    using Value = std::variant<schema::AeroQuadratureNone, schema::AeroQuadratureActive>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T308> {
    using Value =
        std::variant<schema::ConstraintScalePolicyNone, schema::ConstraintScalePolicyPhysical>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T310> {
    using Value = std::variant<schema::ComponentScale3None, schema::ComponentScale3Values>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T318> {
    using Value = schema::LinearAccuracy;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T319> {
    using Value = schema::BackendMemory;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T330> {
    using Value = schema::EigenRequirements;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T331> {
    using Value = schema::EigenLimits;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T337> {
    using Value = std::variant<schema::EigenSelectorLargest, schema::EigenSelectorShift>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T338> {
    using Value = std::variant<schema::ArnoldiSeedUniform, schema::ArnoldiSeedExplicit>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T340> {
    using Value = schema::Complex;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T345> {
    using Value = std::variant<schema::BackwardNone, schema::BackwardRequired>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T349> {
    using Value = schema::InputQuotas;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T350> {
    using Value = schema::OutputCadence;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T351> {
    using Value = schema::Visualization;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T352> {
    using Value = std::variant<schema::CheckpointOutputOff, schema::CheckpointOutputOn>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T360> {
    using Value = schema::Expression;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T362> {
    using Value = schema::Expression;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T363> {
    using Value = std::variant<schema::ExpressionDomainReference, schema::ExpressionDomainTime,
                               schema::ExpressionDomainCurrentPositionTime>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T377> {
    using Value =
        std::variant<schema::FullShellFieldGlobalPolynomial, schema::FullShellFieldHalfCells>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T379> {
    using Value = schema::VectorPolynomial;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T384> {
    using Value = std::variant<schema::InitialPointShell, schema::InitialPointBeam,
                               schema::InitialPointRigid>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T385> {
    using Value = schema::Target;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T386> {
    using Value = schema::Target;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T387> {
    using Value = schema::Target;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T389> {
    using Value = schema::Case;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T390> {
    using Value = schema::PolyTerm;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T391> {
    using Value = std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                               schema::ReferenceFieldExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T392> {
    using Value = std::variant<schema::TimeSeriesConstant, schema::TimeSeriesPolynomial,
                               schema::TimeSeriesPchip, schema::TimeSeriesExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T393> {
    using Value = std::variant<schema::GeometryPg, schema::GeometrySurfaceGraph,
                               schema::GeometryNastran, schema::GeometryGmsh>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T394> {
    using Value = schema::PgTerm;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T395> {
    using Value = schema::Boundary;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T396> {
    using Value = schema::Continuity;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T397> {
    using Value = std::variant<schema::MaterialIsotropic, schema::MaterialOrthotropic,
                               schema::MaterialEffective>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T398> {
    using Value = schema::Layer;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T399> {
    using Value = schema::Shell;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T400> {
    using Value = schema::Beam;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T401> {
    using Value =
        std::variant<schema::JointAdhesive, schema::JointRigidAttachment, schema::JointWeld>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T402> {
    using Value = schema::TraceMap;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T403> {
    using Value = schema::MapTerm;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T404> {
    using Value = std::variant<schema::SupportConstraint, schema::SupportTraceConstraint,
                               schema::SupportElastic>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T405> {
    using Value = std::variant<schema::MassParticle, schema::MassRigid>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T406> {
    using Value = std::variant<schema::LoadPointForce, schema::LoadCouple, schema::LoadGravity,
                               schema::LoadTraction, schema::LoadPressure, schema::LoadBody,
                               schema::LoadCurrentExpression>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T407> {
    using Value = schema::AeroSurface;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T408> {
    using Value = schema::AeroSharing;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T409> {
    using Value = schema::TimeInterval;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T410> {
    using Value = schema::SpeedInterval;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T411> {
    using Value = schema::CoefficientState;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T412> {
    using Value = schema::BeamState;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T413> {
    using Value = schema::RigidState;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T414> {
    using Value = schema::SurfaceWake;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T415> {
    using Value = std::variant<schema::ObservationPoint, schema::ObservationStress,
                               schema::ObservationSection, schema::ObservationReaction,
                               schema::ObservationAero, schema::ObservationMode>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T416> {
    using Value = std::variant<schema::ErrorGroupNotApplicable, schema::ErrorGroupCoveredBy,
                               schema::ErrorGroupRequired>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T417> {
    using Value = std::variant<schema::ScaleUniformComponents, schema::ScalePhysicalLayoutPolicy,
                               schema::ScaleExplicitLayout>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T418> {
    using Value = schema::EigenCluster;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T419> {
    using Value = schema::Expression;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T420> {
    using Value =
        std::variant<schema::ExprNodeConstant, schema::ExprNodeVariable, schema::ExprNodeBinary,
                     schema::ExprNodeNegate, schema::ExprNodePower, schema::ExprNodeAnalytic>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T421> {
    using Value = schema::ShellPhysicalState;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T422> {
    using Value = schema::CellField;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T423> {
    using Value = schema::VectorTerm;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T424> {
    using Value = schema::PointCheck;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T425> {
    using Value = schema::FileReference;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T444> {
    using Value = std::variant<SchemaNumberView, schema::Rational>;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T445> {
    using Value = schema::PolyTerm;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T447> {
    using Value = schema::MaterialMap;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T449> {
    using Value = schema::MaterialMapPiece;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
template <>
struct RecordAccess<FieldType::T451> {
    using Value = schema::MaterialMapTerm;
    static Value read(const SchemaData*, std::uint64_t) noexcept;
};
} // namespace detail

} // namespace uvlas::io
