#pragma once
#include "case_schema_detail.hpp"

#include <exception>
#include <variant>
namespace uvlas::io::detail {
struct SchemaData {
    core::PlannedBuffer<InputDocumentData> snapshot;
    core::PlannedBuffer<SchemaNode> nodes;
    core::PlannedBuffer<std::uint64_t> array_elements;
    std::array<std::uint64_t, 196> record_counts{};
    core::PlannedBuffer<std::optional<schema::Case>> records_Case;
    core::PlannedBuffer<std::optional<schema::Units>> records_Units;
    core::PlannedBuffer<std::optional<schema::Model>> records_Model;
    core::PlannedBuffer<std::optional<schema::Provenance>> records_Provenance;
    core::PlannedBuffer<std::optional<schema::Frame>> records_Frame;
    core::PlannedBuffer<std::optional<schema::Target>> records_Target;
    core::PlannedBuffer<std::optional<schema::MaterialPoint>> records_MaterialPoint;
    core::PlannedBuffer<std::optional<schema::Region>> records_Region;
    core::PlannedBuffer<std::optional<schema::PolyTerm>> records_PolyTerm;
    core::PlannedBuffer<std::optional<schema::ReferenceFieldConstant>>
        records_ReferenceFieldConstant;
    core::PlannedBuffer<std::optional<schema::ReferenceFieldPolynomial>>
        records_ReferenceFieldPolynomial;
    core::PlannedBuffer<std::optional<schema::TimeSeriesConstant>> records_TimeSeriesConstant;
    core::PlannedBuffer<std::optional<schema::TimeSeriesPolynomial>> records_TimeSeriesPolynomial;
    core::PlannedBuffer<std::optional<schema::TimeSeriesPchip>> records_TimeSeriesPchip;
    core::PlannedBuffer<std::optional<schema::OutsideRejectHold>> records_OutsideRejectHold;
    core::PlannedBuffer<std::optional<schema::OutsidePeriodic>> records_OutsidePeriodic;
    core::PlannedBuffer<std::optional<schema::Motion>> records_Motion;
    core::PlannedBuffer<std::optional<schema::AxisProjected>> records_AxisProjected;
    core::PlannedBuffer<std::optional<schema::AxisExplicit>> records_AxisExplicit;
    core::PlannedBuffer<std::optional<schema::GeometryPg>> records_GeometryPg;
    core::PlannedBuffer<std::optional<schema::PgTerm>> records_PgTerm;
    core::PlannedBuffer<std::optional<schema::Thickness>> records_Thickness;
    core::PlannedBuffer<std::optional<schema::Boundary>> records_Boundary;
    core::PlannedBuffer<std::optional<schema::Continuity>> records_Continuity;
    core::PlannedBuffer<std::optional<schema::GeometrySurfaceGraph>> records_GeometrySurfaceGraph;
    core::PlannedBuffer<std::optional<schema::GeometryNastran>> records_GeometryNastran;
    core::PlannedBuffer<std::optional<schema::GeometryGmsh>> records_GeometryGmsh;
    core::PlannedBuffer<std::optional<schema::ReferenceGeometrySupplied>>
        records_ReferenceGeometrySupplied;
    core::PlannedBuffer<std::optional<schema::ReferenceGeometryFit>> records_ReferenceGeometryFit;
    core::PlannedBuffer<std::optional<schema::MaterialIsotropic>> records_MaterialIsotropic;
    core::PlannedBuffer<std::optional<schema::MaterialOrthotropic>> records_MaterialOrthotropic;
    core::PlannedBuffer<std::optional<schema::MaterialEffective>> records_MaterialEffective;
    core::PlannedBuffer<std::optional<schema::ApplicabilityUnevaluated>>
        records_ApplicabilityUnevaluated;
    core::PlannedBuffer<std::optional<schema::ApplicabilityDeclared>> records_ApplicabilityDeclared;
    core::PlannedBuffer<std::optional<schema::Layer>> records_Layer;
    core::PlannedBuffer<std::optional<schema::Shell>> records_Shell;
    core::PlannedBuffer<std::optional<schema::Beam>> records_Beam;
    core::PlannedBuffer<std::optional<schema::BeamCoordinateArclength>>
        records_BeamCoordinateArclength;
    core::PlannedBuffer<std::optional<schema::BeamCoordinateScaled>> records_BeamCoordinateScaled;
    core::PlannedBuffer<std::optional<schema::SectionMatrix>> records_SectionMatrix;
    core::PlannedBuffer<std::optional<schema::SectionCircle>> records_SectionCircle;
    core::PlannedBuffer<std::optional<schema::SectionRectangle>> records_SectionRectangle;
    core::PlannedBuffer<std::optional<schema::JointAdhesive>> records_JointAdhesive;
    core::PlannedBuffer<std::optional<schema::TraceMap>> records_TraceMap;
    core::PlannedBuffer<std::optional<schema::MapTerm>> records_MapTerm;
    core::PlannedBuffer<std::optional<schema::JointRigidAttachment>> records_JointRigidAttachment;
    core::PlannedBuffer<std::optional<schema::JointWeld>> records_JointWeld;
    core::PlannedBuffer<std::optional<schema::SupportConstraint>> records_SupportConstraint;
    core::PlannedBuffer<std::optional<schema::SupportTraceConstraint>>
        records_SupportTraceConstraint;
    core::PlannedBuffer<std::optional<schema::SupportElastic>> records_SupportElastic;
    core::PlannedBuffer<std::optional<schema::MassParticle>> records_MassParticle;
    core::PlannedBuffer<std::optional<schema::MassRigid>> records_MassRigid;
    core::PlannedBuffer<std::optional<schema::LedgerAdd>> records_LedgerAdd;
    core::PlannedBuffer<std::optional<schema::LedgerReplace>> records_LedgerReplace;
    core::PlannedBuffer<std::optional<schema::LoadPointForce>> records_LoadPointForce;
    core::PlannedBuffer<std::optional<schema::LoadCouple>> records_LoadCouple;
    core::PlannedBuffer<std::optional<schema::LoadGravity>> records_LoadGravity;
    core::PlannedBuffer<std::optional<schema::LoadTraction>> records_LoadTraction;
    core::PlannedBuffer<std::optional<schema::LoadPressure>> records_LoadPressure;
    core::PlannedBuffer<std::optional<schema::LoadBody>> records_LoadBody;
    core::PlannedBuffer<std::optional<schema::AeroSurface>> records_AeroSurface;
    core::PlannedBuffer<std::optional<schema::AeroAttachmentFixed>> records_AeroAttachmentFixed;
    core::PlannedBuffer<std::optional<schema::AeroAttachmentPrescribed>>
        records_AeroAttachmentPrescribed;
    core::PlannedBuffer<std::optional<schema::AeroAttachmentStructural>>
        records_AeroAttachmentStructural;
    core::PlannedBuffer<std::optional<schema::AeroSharing>> records_AeroSharing;
    core::PlannedBuffer<std::optional<schema::FlowNone>> records_FlowNone;
    core::PlannedBuffer<std::optional<schema::FlowKnown>> records_FlowKnown;
    core::PlannedBuffer<std::optional<schema::GustNone>> records_GustNone;
    core::PlannedBuffer<std::optional<schema::GustSine>> records_GustSine;
    core::PlannedBuffer<std::optional<schema::FlowReference>> records_FlowReference;
    core::PlannedBuffer<std::optional<schema::WakeNone>> records_WakeNone;
    core::PlannedBuffer<std::optional<schema::WakeActive>> records_WakeActive;
    core::PlannedBuffer<std::optional<schema::ConvectionFree>> records_ConvectionFree;
    core::PlannedBuffer<std::optional<schema::ConvectionPrescribed>> records_ConvectionPrescribed;
    core::PlannedBuffer<std::optional<schema::RetentionAll>> records_RetentionAll;
    core::PlannedBuffer<std::optional<schema::RetentionWindow>> records_RetentionWindow;
    core::PlannedBuffer<std::optional<schema::AnalysisStatic>> records_AnalysisStatic;
    core::PlannedBuffer<std::optional<schema::AnalysisTransient>> records_AnalysisTransient;
    core::PlannedBuffer<std::optional<schema::AnalysisFlutter>> records_AnalysisFlutter;
    core::PlannedBuffer<std::optional<schema::AnalysisSpeedNone>> records_AnalysisSpeedNone;
    core::PlannedBuffer<std::optional<schema::AnalysisSpeedValue>> records_AnalysisSpeedValue;
    core::PlannedBuffer<std::optional<schema::Branch>> records_Branch;
    core::PlannedBuffer<std::optional<schema::BranchSeedReference>> records_BranchSeedReference;
    core::PlannedBuffer<std::optional<schema::BranchSeedPrepared>> records_BranchSeedPrepared;
    core::PlannedBuffer<std::optional<schema::ContinuationFixed>> records_ContinuationFixed;
    core::PlannedBuffer<std::optional<schema::ContinuationSweep>> records_ContinuationSweep;
    core::PlannedBuffer<std::optional<schema::ArcLength>> records_ArcLength;
    core::PlannedBuffer<std::optional<schema::TimeInterval>> records_TimeInterval;
    core::PlannedBuffer<std::optional<schema::SpeedInterval>> records_SpeedInterval;
    core::PlannedBuffer<std::optional<schema::Search>> records_Search;
    core::PlannedBuffer<std::optional<schema::Association>> records_Association;
    core::PlannedBuffer<std::optional<schema::Coverage>> records_Coverage;
    core::PlannedBuffer<std::optional<schema::InitialEquilibrium>> records_InitialEquilibrium;
    core::PlannedBuffer<std::optional<schema::InitialCold>> records_InitialCold;
    core::PlannedBuffer<std::optional<schema::InitialPrepared>> records_InitialPrepared;
    core::PlannedBuffer<std::optional<schema::PerturbationNone>> records_PerturbationNone;
    core::PlannedBuffer<std::optional<schema::InitialState>> records_InitialState;
    core::PlannedBuffer<std::optional<schema::CoefficientState>> records_CoefficientState;
    core::PlannedBuffer<std::optional<schema::PlanesBinary64>> records_PlanesBinary64;
    core::PlannedBuffer<std::optional<schema::PlanesTwoWord>> records_PlanesTwoWord;
    core::PlannedBuffer<std::optional<schema::BeamState>> records_BeamState;
    core::PlannedBuffer<std::optional<schema::RigidState>> records_RigidState;
    core::PlannedBuffer<std::optional<schema::FrameStatePayload>> records_FrameStatePayload;
    core::PlannedBuffer<std::optional<schema::InitialWake>> records_InitialWake;
    core::PlannedBuffer<std::optional<schema::SurfaceWake>> records_SurfaceWake;
    core::PlannedBuffer<std::optional<schema::Restart>> records_Restart;
    core::PlannedBuffer<std::optional<schema::ObservationPoint>> records_ObservationPoint;
    core::PlannedBuffer<std::optional<schema::ObservationStress>> records_ObservationStress;
    core::PlannedBuffer<std::optional<schema::ObservationSection>> records_ObservationSection;
    core::PlannedBuffer<std::optional<schema::ObservationReaction>> records_ObservationReaction;
    core::PlannedBuffer<std::optional<schema::ObservationAero>> records_ObservationAero;
    core::PlannedBuffer<std::optional<schema::ObservationMode>> records_ObservationMode;
    core::PlannedBuffer<std::optional<schema::BasisWorldMaterial>> records_BasisWorldMaterial;
    core::PlannedBuffer<std::optional<schema::BasisExplicit>> records_BasisExplicit;
    core::PlannedBuffer<std::optional<schema::Comparison>> records_Comparison;
    core::PlannedBuffer<std::optional<schema::ComparisonReferenceNone>>
        records_ComparisonReferenceNone;
    core::PlannedBuffer<std::optional<schema::ComparisonReferenceFile>>
        records_ComparisonReferenceFile;
    core::PlannedBuffer<std::optional<schema::ObservationBudget>> records_ObservationBudget;
    core::PlannedBuffer<std::optional<schema::ErrorGroupNotApplicable>>
        records_ErrorGroupNotApplicable;
    core::PlannedBuffer<std::optional<schema::ErrorGroupCoveredBy>> records_ErrorGroupCoveredBy;
    core::PlannedBuffer<std::optional<schema::ErrorGroupRequired>> records_ErrorGroupRequired;
    core::PlannedBuffer<std::optional<schema::Numerics>> records_Numerics;
    core::PlannedBuffer<std::optional<schema::WakeCapacityNone>> records_WakeCapacityNone;
    core::PlannedBuffer<std::optional<schema::WakeCapacityLimited>> records_WakeCapacityLimited;
    core::PlannedBuffer<std::optional<schema::PatchPolicyNone>> records_PatchPolicyNone;
    core::PlannedBuffer<std::optional<schema::PatchPolicyActive>> records_PatchPolicyActive;
    core::PlannedBuffer<std::optional<schema::Quadrature>> records_Quadrature;
    core::PlannedBuffer<std::optional<schema::StructuralQuadratureNone>>
        records_StructuralQuadratureNone;
    core::PlannedBuffer<std::optional<schema::StructuralQuadratureActive>>
        records_StructuralQuadratureActive;
    core::PlannedBuffer<std::optional<schema::AeroQuadratureNone>> records_AeroQuadratureNone;
    core::PlannedBuffer<std::optional<schema::AeroQuadratureActive>> records_AeroQuadratureActive;
    core::PlannedBuffer<std::optional<schema::ScaleUniformComponents>>
        records_ScaleUniformComponents;
    core::PlannedBuffer<std::optional<schema::ScalePhysicalLayoutPolicy>>
        records_ScalePhysicalLayoutPolicy;
    core::PlannedBuffer<std::optional<schema::ConstraintScalePolicyNone>>
        records_ConstraintScalePolicyNone;
    core::PlannedBuffer<std::optional<schema::ConstraintScalePolicyPhysical>>
        records_ConstraintScalePolicyPhysical;
    core::PlannedBuffer<std::optional<schema::ComponentScale3None>> records_ComponentScale3None;
    core::PlannedBuffer<std::optional<schema::ComponentScale3Values>> records_ComponentScale3Values;
    core::PlannedBuffer<std::optional<schema::ScaleExplicitLayout>> records_ScaleExplicitLayout;
    core::PlannedBuffer<std::optional<schema::LinearPolicyDirect>> records_LinearPolicyDirect;
    core::PlannedBuffer<std::optional<schema::LinearPolicyGmres>> records_LinearPolicyGmres;
    core::PlannedBuffer<std::optional<schema::LinearAccuracy>> records_LinearAccuracy;
    core::PlannedBuffer<std::optional<schema::NewtonPolicy>> records_NewtonPolicy;
    core::PlannedBuffer<std::optional<schema::PrecisionPolicy>> records_PrecisionPolicy;
    core::PlannedBuffer<std::optional<schema::BackendMemory>> records_BackendMemory;
    core::PlannedBuffer<std::optional<schema::EigenPolicyNone>> records_EigenPolicyNone;
    core::PlannedBuffer<std::optional<schema::EigenPolicyFull>> records_EigenPolicyFull;
    core::PlannedBuffer<std::optional<schema::EigenPolicyGeneralized>>
        records_EigenPolicyGeneralized;
    core::PlannedBuffer<std::optional<schema::EigenPolicyArnoldi>> records_EigenPolicyArnoldi;
    core::PlannedBuffer<std::optional<schema::Complex>> records_Complex;
    core::PlannedBuffer<std::optional<schema::ArnoldiSeedUniform>> records_ArnoldiSeedUniform;
    core::PlannedBuffer<std::optional<schema::ArnoldiSeedExplicit>> records_ArnoldiSeedExplicit;
    core::PlannedBuffer<std::optional<schema::EigenSelectorLargest>> records_EigenSelectorLargest;
    core::PlannedBuffer<std::optional<schema::EigenSelectorShift>> records_EigenSelectorShift;
    core::PlannedBuffer<std::optional<schema::EigenRequirements>> records_EigenRequirements;
    core::PlannedBuffer<std::optional<schema::BackwardNone>> records_BackwardNone;
    core::PlannedBuffer<std::optional<schema::BackwardRequired>> records_BackwardRequired;
    core::PlannedBuffer<std::optional<schema::EigenCluster>> records_EigenCluster;
    core::PlannedBuffer<std::optional<schema::EigenLimits>> records_EigenLimits;
    core::PlannedBuffer<std::optional<schema::Resources>> records_Resources;
    core::PlannedBuffer<std::optional<schema::InputQuotas>> records_InputQuotas;
    core::PlannedBuffer<std::optional<schema::Output>> records_Output;
    core::PlannedBuffer<std::optional<schema::OutputCadence>> records_OutputCadence;
    core::PlannedBuffer<std::optional<schema::Visualization>> records_Visualization;
    core::PlannedBuffer<std::optional<schema::CheckpointOutputOff>> records_CheckpointOutputOff;
    core::PlannedBuffer<std::optional<schema::CheckpointOutputOn>> records_CheckpointOutputOn;
    core::PlannedBuffer<std::optional<schema::ReferenceFieldExpression>>
        records_ReferenceFieldExpression;
    core::PlannedBuffer<std::optional<schema::TimeSeriesExpression>> records_TimeSeriesExpression;
    core::PlannedBuffer<std::optional<schema::LoadCurrentExpression>> records_LoadCurrentExpression;
    core::PlannedBuffer<std::optional<schema::Expression>> records_Expression;
    core::PlannedBuffer<std::optional<schema::ExpressionDomainReference>>
        records_ExpressionDomainReference;
    core::PlannedBuffer<std::optional<schema::ExpressionDomainTime>> records_ExpressionDomainTime;
    core::PlannedBuffer<std::optional<schema::ExpressionDomainCurrentPositionTime>>
        records_ExpressionDomainCurrentPositionTime;
    core::PlannedBuffer<std::optional<schema::ExprNodeConstant>> records_ExprNodeConstant;
    core::PlannedBuffer<std::optional<schema::ExprNodeVariable>> records_ExprNodeVariable;
    core::PlannedBuffer<std::optional<schema::ExprNodeBinary>> records_ExprNodeBinary;
    core::PlannedBuffer<std::optional<schema::ExprNodeNegate>> records_ExprNodeNegate;
    core::PlannedBuffer<std::optional<schema::ExprNodePower>> records_ExprNodePower;
    core::PlannedBuffer<std::optional<schema::ExprNodeAnalytic>> records_ExprNodeAnalytic;
    core::PlannedBuffer<std::optional<schema::PerturbationPhysical>> records_PerturbationPhysical;
    core::PlannedBuffer<std::optional<schema::ShellPhysicalState>> records_ShellPhysicalState;
    core::PlannedBuffer<std::optional<schema::FullShellFieldGlobalPolynomial>>
        records_FullShellFieldGlobalPolynomial;
    core::PlannedBuffer<std::optional<schema::FullShellFieldHalfCells>>
        records_FullShellFieldHalfCells;
    core::PlannedBuffer<std::optional<schema::CellField>> records_CellField;
    core::PlannedBuffer<std::optional<schema::VectorPolynomial>> records_VectorPolynomial;
    core::PlannedBuffer<std::optional<schema::VectorTerm>> records_VectorTerm;
    core::PlannedBuffer<std::optional<schema::PointCheck>> records_PointCheck;
    core::PlannedBuffer<std::optional<schema::InitialPointShell>> records_InitialPointShell;
    core::PlannedBuffer<std::optional<schema::InitialPointBeam>> records_InitialPointBeam;
    core::PlannedBuffer<std::optional<schema::InitialPointRigid>> records_InitialPointRigid;
    core::PlannedBuffer<std::optional<schema::Rational>> records_Rational;
    core::PlannedBuffer<std::optional<schema::FileReference>> records_FileReference;
    core::PlannedBuffer<std::optional<schema::MaterialMap>> records_MaterialMap;
    core::PlannedBuffer<std::optional<schema::MaterialMapPiece>> records_MaterialMapPiece;
    core::PlannedBuffer<std::optional<schema::MaterialMapTerm>> records_MaterialMapTerm;
    core::PlannedBuffer<std::optional<schema::EdgePairing>> records_EdgePairing;
    core::PlannedBuffer<std::optional<schema::EdgePairingKnot>> records_EdgePairingKnot;

    core::PlannedBuffer<char> digits;
    core::PlannedBuffer<SchemaIdentifier> identifiers;
    core::PlannedBuffer<SchemaReference> references;
    std::uint64_t digit_count = 0;
    std::uint64_t identifier_count = 0;
    std::uint64_t reference_count = 0;
    SchemaMemoryReport report{0, 0, 0, 0, 0, 0, 0};
};
template <class T>
const T& ready_record(const std::optional<T>& record) noexcept {
    // Only a fully populated factory-owned arena can enter a public SchemaCase.
    if (!record)
        std::terminate();
    return *record;
}
struct SchemaFactory {
    template <FieldType Type>
    static schema::Field<Type> field(const SchemaData* owner, std::uint64_t node) noexcept {
        return schema::Field<Type>{owner, node};
    }
    template <FieldType Type>
    static schema::Field<Type> required_field(SchemaContext& context, std::uint64_t node,
                                              std::string_view key) {
        const auto* owner = context.output;
        const auto child = schema_find(context, node, key);
        if (!child)
            throw core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel};
        return field<Type>(owner, *child);
    }
    template <FieldType Type>
    static std::optional<schema::Field<Type>>
    optional_field(SchemaContext& context, std::uint64_t node, std::string_view key) {
        const auto* owner = context.output;
        const auto child = schema_find(context, node, key);
        if (!child)
            return std::nullopt;
        return field<Type>(owner, *child);
    }
    static SchemaCase case_result(core::PlannedBuffer<SchemaData> data) noexcept {
        return SchemaCase{std::move(data)};
    }
    static void diagnostic(SchemaDiagnostic& slot, SchemaDiagnosticDetail fixed,
                           core::PlannedBuffer<char> path = {}) noexcept {
        slot.detail_ = fixed;
        slot.path_ = std::move(path);
    }
};
} // namespace uvlas::io::detail
