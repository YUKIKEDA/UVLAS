// Derived test-only public traversal from case_schema_snapshot.hpp SHA-256
// 0737379415d03c4c890dbb3f8253a6c4b866a86c17f125936e5c0648645f07f1. Expected numbers are never
// generated; this snapshots every exposed SI value for preservation. Generated public traversal; no
// product value is used as an oracle. Interface SHA-256:
// cbd38a49b742828fdfdafcfd1c50f172cdaa9144ad4ca9c30f041c729e8171a6 Generator SHA-256:
// 64780ef13564cd37ed97305f664f1e00ae92775bbedd33ee9170fa65ef0be527
#pragma once
#include "case_schema_snapshot.hpp"

#include <bit>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <uvlas/io/case_schema.hpp>
#include <uvlas/io/si_case.hpp>
#include <variant>
#include <vector>

namespace uvlas::test::si_ownership_snapshot {
struct Bytes : std::vector<std::uint8_t> {
    explicit Bytes(const io::SiCase& value) : owner(value) {}
    const io::SiCase& owner;
    std::vector<io::SiQuantityKey> keys;
};
struct Snapshot {
    std::vector<std::uint8_t> bytes;
    std::vector<io::SiQuantityKey> keys;
    bool operator==(const Snapshot&) const = default;
};
inline void integer(Bytes& out, std::uint64_t value) {
    for (unsigned shift = 0; shift < 64; shift += 8)
        out.push_back(static_cast<std::uint8_t>(value >> shift));
}
inline void text(Bytes& out, std::string_view value) {
    integer(out, value.size());
    for (const char c : value)
        out.push_back(static_cast<std::uint8_t>(c));
}
inline void append(Bytes& out, const io::SchemaNumberView& value) {
    integer(out, static_cast<std::uint64_t>(value.sign));
    text(out, value.digits);
    integer(out, static_cast<std::uint64_t>(value.exponent10));
    integer(out, static_cast<std::uint64_t>(value.original_kind));
    for (const auto exponent : value.dimensions.exponents)
        integer(out, static_cast<std::uint64_t>(exponent));
}
inline void append(Bytes&, const io::schema::Case&);
inline void append(Bytes&, const io::schema::Units&);
inline void append(Bytes&, const io::schema::Model&);
inline void append(Bytes&, const io::schema::Provenance&);
inline void append(Bytes&, const io::schema::Frame&);
inline void append(Bytes&, const io::schema::Target&);
inline void append(Bytes&, const io::schema::MaterialPoint&);
inline void append(Bytes&, const io::schema::Region&);
inline void append(Bytes&, const io::schema::PolyTerm&);
inline void append(Bytes&, const io::schema::ReferenceFieldConstant&);
inline void append(Bytes&, const io::schema::ReferenceFieldPolynomial&);
inline void append(Bytes&, const io::schema::TimeSeriesConstant&);
inline void append(Bytes&, const io::schema::TimeSeriesPolynomial&);
inline void append(Bytes&, const io::schema::TimeSeriesPchip&);
inline void append(Bytes&, const io::schema::OutsideRejectHold&);
inline void append(Bytes&, const io::schema::OutsidePeriodic&);
inline void append(Bytes&, const io::schema::Motion&);
inline void append(Bytes&, const io::schema::AxisProjected&);
inline void append(Bytes&, const io::schema::AxisExplicit&);
inline void append(Bytes&, const io::schema::GeometryPg&);
inline void append(Bytes&, const io::schema::PgTerm&);
inline void append(Bytes&, const io::schema::Thickness&);
inline void append(Bytes&, const io::schema::Boundary&);
inline void append(Bytes&, const io::schema::Continuity&);
inline void append(Bytes&, const io::schema::GeometrySurfaceGraph&);
inline void append(Bytes&, const io::schema::GeometryNastran&);
inline void append(Bytes&, const io::schema::GeometryGmsh&);
inline void append(Bytes&, const io::schema::ReferenceGeometrySupplied&);
inline void append(Bytes&, const io::schema::ReferenceGeometryFit&);
inline void append(Bytes&, const io::schema::MaterialIsotropic&);
inline void append(Bytes&, const io::schema::MaterialOrthotropic&);
inline void append(Bytes&, const io::schema::MaterialEffective&);
inline void append(Bytes&, const io::schema::ApplicabilityUnevaluated&);
inline void append(Bytes&, const io::schema::ApplicabilityDeclared&);
inline void append(Bytes&, const io::schema::Layer&);
inline void append(Bytes&, const io::schema::Shell&);
inline void append(Bytes&, const io::schema::Beam&);
inline void append(Bytes&, const io::schema::BeamCoordinateArclength&);
inline void append(Bytes&, const io::schema::BeamCoordinateScaled&);
inline void append(Bytes&, const io::schema::SectionMatrix&);
inline void append(Bytes&, const io::schema::SectionCircle&);
inline void append(Bytes&, const io::schema::SectionRectangle&);
inline void append(Bytes&, const io::schema::JointAdhesive&);
inline void append(Bytes&, const io::schema::TraceMap&);
inline void append(Bytes&, const io::schema::MapTerm&);
inline void append(Bytes&, const io::schema::JointRigidAttachment&);
inline void append(Bytes&, const io::schema::JointWeld&);
inline void append(Bytes&, const io::schema::SupportConstraint&);
inline void append(Bytes&, const io::schema::SupportTraceConstraint&);
inline void append(Bytes&, const io::schema::SupportElastic&);
inline void append(Bytes&, const io::schema::MassParticle&);
inline void append(Bytes&, const io::schema::MassRigid&);
inline void append(Bytes&, const io::schema::LedgerAdd&);
inline void append(Bytes&, const io::schema::LedgerReplace&);
inline void append(Bytes&, const io::schema::LoadPointForce&);
inline void append(Bytes&, const io::schema::LoadCouple&);
inline void append(Bytes&, const io::schema::LoadGravity&);
inline void append(Bytes&, const io::schema::LoadTraction&);
inline void append(Bytes&, const io::schema::LoadPressure&);
inline void append(Bytes&, const io::schema::LoadBody&);
inline void append(Bytes&, const io::schema::AeroSurface&);
inline void append(Bytes&, const io::schema::AeroAttachmentFixed&);
inline void append(Bytes&, const io::schema::AeroAttachmentPrescribed&);
inline void append(Bytes&, const io::schema::AeroAttachmentStructural&);
inline void append(Bytes&, const io::schema::AeroSharing&);
inline void append(Bytes&, const io::schema::FlowNone&);
inline void append(Bytes&, const io::schema::FlowKnown&);
inline void append(Bytes&, const io::schema::GustNone&);
inline void append(Bytes&, const io::schema::GustSine&);
inline void append(Bytes&, const io::schema::FlowReference&);
inline void append(Bytes&, const io::schema::WakeNone&);
inline void append(Bytes&, const io::schema::WakeActive&);
inline void append(Bytes&, const io::schema::ConvectionFree&);
inline void append(Bytes&, const io::schema::ConvectionPrescribed&);
inline void append(Bytes&, const io::schema::RetentionAll&);
inline void append(Bytes&, const io::schema::RetentionWindow&);
inline void append(Bytes&, const io::schema::AnalysisStatic&);
inline void append(Bytes&, const io::schema::AnalysisTransient&);
inline void append(Bytes&, const io::schema::AnalysisFlutter&);
inline void append(Bytes&, const io::schema::AnalysisSpeedNone&);
inline void append(Bytes&, const io::schema::AnalysisSpeedValue&);
inline void append(Bytes&, const io::schema::Branch&);
inline void append(Bytes&, const io::schema::BranchSeedReference&);
inline void append(Bytes&, const io::schema::BranchSeedPrepared&);
inline void append(Bytes&, const io::schema::ContinuationFixed&);
inline void append(Bytes&, const io::schema::ContinuationSweep&);
inline void append(Bytes&, const io::schema::ArcLength&);
inline void append(Bytes&, const io::schema::TimeInterval&);
inline void append(Bytes&, const io::schema::SpeedInterval&);
inline void append(Bytes&, const io::schema::Search&);
inline void append(Bytes&, const io::schema::Association&);
inline void append(Bytes&, const io::schema::Coverage&);
inline void append(Bytes&, const io::schema::InitialEquilibrium&);
inline void append(Bytes&, const io::schema::InitialCold&);
inline void append(Bytes&, const io::schema::InitialPrepared&);
inline void append(Bytes&, const io::schema::PerturbationNone&);
inline void append(Bytes&, const io::schema::InitialState&);
inline void append(Bytes&, const io::schema::CoefficientState&);
inline void append(Bytes&, const io::schema::PlanesBinary64&);
inline void append(Bytes&, const io::schema::PlanesTwoWord&);
inline void append(Bytes&, const io::schema::BeamState&);
inline void append(Bytes&, const io::schema::RigidState&);
inline void append(Bytes&, const io::schema::FrameStatePayload&);
inline void append(Bytes&, const io::schema::InitialWake&);
inline void append(Bytes&, const io::schema::SurfaceWake&);
inline void append(Bytes&, const io::schema::Restart&);
inline void append(Bytes&, const io::schema::ObservationPoint&);
inline void append(Bytes&, const io::schema::ObservationStress&);
inline void append(Bytes&, const io::schema::ObservationSection&);
inline void append(Bytes&, const io::schema::ObservationReaction&);
inline void append(Bytes&, const io::schema::ObservationAero&);
inline void append(Bytes&, const io::schema::ObservationMode&);
inline void append(Bytes&, const io::schema::BasisWorldMaterial&);
inline void append(Bytes&, const io::schema::BasisExplicit&);
inline void append(Bytes&, const io::schema::Comparison&);
inline void append(Bytes&, const io::schema::ComparisonReferenceNone&);
inline void append(Bytes&, const io::schema::ComparisonReferenceFile&);
inline void append(Bytes&, const io::schema::ObservationBudget&);
inline void append(Bytes&, const io::schema::ErrorGroupNotApplicable&);
inline void append(Bytes&, const io::schema::ErrorGroupCoveredBy&);
inline void append(Bytes&, const io::schema::ErrorGroupRequired&);
inline void append(Bytes&, const io::schema::Numerics&);
inline void append(Bytes&, const io::schema::WakeCapacityNone&);
inline void append(Bytes&, const io::schema::WakeCapacityLimited&);
inline void append(Bytes&, const io::schema::PatchPolicyNone&);
inline void append(Bytes&, const io::schema::PatchPolicyActive&);
inline void append(Bytes&, const io::schema::Quadrature&);
inline void append(Bytes&, const io::schema::StructuralQuadratureNone&);
inline void append(Bytes&, const io::schema::StructuralQuadratureActive&);
inline void append(Bytes&, const io::schema::AeroQuadratureNone&);
inline void append(Bytes&, const io::schema::AeroQuadratureActive&);
inline void append(Bytes&, const io::schema::ScaleUniformComponents&);
inline void append(Bytes&, const io::schema::ScalePhysicalLayoutPolicy&);
inline void append(Bytes&, const io::schema::ConstraintScalePolicyNone&);
inline void append(Bytes&, const io::schema::ConstraintScalePolicyPhysical&);
inline void append(Bytes&, const io::schema::ComponentScale3None&);
inline void append(Bytes&, const io::schema::ComponentScale3Values&);
inline void append(Bytes&, const io::schema::ScaleExplicitLayout&);
inline void append(Bytes&, const io::schema::LinearPolicyDirect&);
inline void append(Bytes&, const io::schema::LinearPolicyGmres&);
inline void append(Bytes&, const io::schema::LinearAccuracy&);
inline void append(Bytes&, const io::schema::NewtonPolicy&);
inline void append(Bytes&, const io::schema::PrecisionPolicy&);
inline void append(Bytes&, const io::schema::BackendMemory&);
inline void append(Bytes&, const io::schema::EigenPolicyNone&);
inline void append(Bytes&, const io::schema::EigenPolicyFull&);
inline void append(Bytes&, const io::schema::EigenPolicyGeneralized&);
inline void append(Bytes&, const io::schema::EigenPolicyArnoldi&);
inline void append(Bytes&, const io::schema::Complex&);
inline void append(Bytes&, const io::schema::ArnoldiSeedUniform&);
inline void append(Bytes&, const io::schema::ArnoldiSeedExplicit&);
inline void append(Bytes&, const io::schema::EigenSelectorLargest&);
inline void append(Bytes&, const io::schema::EigenSelectorShift&);
inline void append(Bytes&, const io::schema::EigenRequirements&);
inline void append(Bytes&, const io::schema::BackwardNone&);
inline void append(Bytes&, const io::schema::BackwardRequired&);
inline void append(Bytes&, const io::schema::EigenCluster&);
inline void append(Bytes&, const io::schema::EigenLimits&);
inline void append(Bytes&, const io::schema::Resources&);
inline void append(Bytes&, const io::schema::InputQuotas&);
inline void append(Bytes&, const io::schema::Output&);
inline void append(Bytes&, const io::schema::OutputCadence&);
inline void append(Bytes&, const io::schema::Visualization&);
inline void append(Bytes&, const io::schema::CheckpointOutputOff&);
inline void append(Bytes&, const io::schema::CheckpointOutputOn&);
inline void append(Bytes&, const io::schema::ReferenceFieldExpression&);
inline void append(Bytes&, const io::schema::TimeSeriesExpression&);
inline void append(Bytes&, const io::schema::LoadCurrentExpression&);
inline void append(Bytes&, const io::schema::Expression&);
inline void append(Bytes&, const io::schema::ExpressionDomainReference&);
inline void append(Bytes&, const io::schema::ExpressionDomainTime&);
inline void append(Bytes&, const io::schema::ExpressionDomainCurrentPositionTime&);
inline void append(Bytes&, const io::schema::ExprNodeConstant&);
inline void append(Bytes&, const io::schema::ExprNodeVariable&);
inline void append(Bytes&, const io::schema::ExprNodeBinary&);
inline void append(Bytes&, const io::schema::ExprNodeNegate&);
inline void append(Bytes&, const io::schema::ExprNodePower&);
inline void append(Bytes&, const io::schema::ExprNodeAnalytic&);
inline void append(Bytes&, const io::schema::PerturbationPhysical&);
inline void append(Bytes&, const io::schema::ShellPhysicalState&);
inline void append(Bytes&, const io::schema::FullShellFieldGlobalPolynomial&);
inline void append(Bytes&, const io::schema::FullShellFieldHalfCells&);
inline void append(Bytes&, const io::schema::CellField&);
inline void append(Bytes&, const io::schema::VectorPolynomial&);
inline void append(Bytes&, const io::schema::VectorTerm&);
inline void append(Bytes&, const io::schema::PointCheck&);
inline void append(Bytes&, const io::schema::InitialPointShell&);
inline void append(Bytes&, const io::schema::InitialPointBeam&);
inline void append(Bytes&, const io::schema::InitialPointRigid&);
inline void append(Bytes&, const io::schema::Rational&);
inline void append(Bytes&, const io::schema::FileReference&);

template <io::detail::FieldType Type>
void append(Bytes&, const io::schema::Field<Type>&);
template <class... T>
void append(Bytes& out, const std::variant<T...>& value) {
    integer(out, value.index());
    std::visit([&](const auto& selected) { append(out, selected); }, value);
}
template <class T>
void append(Bytes& out, const std::optional<T>& value) {
    integer(out, value.has_value() ? 1U : 0U);
    if (value)
        append(out, *value);
}
template <io::detail::FieldType Type>
void append(Bytes& out, const io::schema::Field<Type>& value) {
    // The static field type distinguishes units, ID namespaces, and selected schemas.
    integer(out, static_cast<std::uint64_t>(Type));
    constexpr auto kind = io::detail::field_kinds[static_cast<std::size_t>(Type)];
    if constexpr (kind == io::detail::FieldKind::Number ||
                  kind == io::detail::FieldKind::ExactQuantity ||
                  Type == io::detail::FieldType::T62) {
        const auto quantity = out.owner.quantity(value);
        if (!quantity)
            throw std::runtime_error("Owned SI quantity lookup failed during snapshot");
        out.keys.push_back(quantity->key);
        integer(out, quantity->ledger_origin);
        integer(out, static_cast<std::uint64_t>(quantity->recipe.pi_power));
        for (const auto e : quantity->recipe.dimensions.exponents)
            integer(out, static_cast<std::uint64_t>(e));
        const auto units = quantity->recipe.units;
        integer(out, static_cast<std::uint64_t>(units.length));
        integer(out, static_cast<std::uint64_t>(units.mass));
        integer(out, static_cast<std::uint64_t>(units.time));
        integer(out, static_cast<std::uint64_t>(units.angle));
        text(out, quantity->recipe.original);
        integer(out, quantity->recipe.value.index());
        const auto decimal = [&](const io::SiDecimalView& x) {
            integer(out, static_cast<std::uint64_t>(x.sign));
            text(out, x.digits);
            integer(out, static_cast<std::uint64_t>(x.exponent10));
        };
        std::visit(
            [&](const auto& original) {
                using Original = std::decay_t<decltype(original)>;
                if constexpr (std::is_same_v<Original, io::SiDecimalView>)
                    decimal(original);
                else {
                    decimal(original.numerator);
                    decimal(original.denominator);
                }
            },
            quantity->recipe.value);
        const auto request_exact = [&](const io::SiExactInput& input) {
            integer(out, input.index());
            std::visit(
                [&](const auto& x) {
                    using X = std::decay_t<decltype(x)>;
                    if constexpr (std::is_same_v<X, io::SiDecimalInput>)
                        text(out, x.token);
                    else {
                        text(out, x.numerator);
                        text(out, x.denominator);
                    }
                },
                input);
        };
        request_exact(quantity->requirement.scale.value);
        for (const auto e : quantity->requirement.scale.dimensions.exponents)
            integer(out, static_cast<std::uint64_t>(e));
        request_exact(quantity->requirement.absolute);
        request_exact(quantity->requirement.relative);
        integer(out, quantity->enclosure.index());
        std::visit(
            [&](const auto& x) {
                using X = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<X, io::ScaledBinary64>)
                    integer(out, std::bit_cast<std::uint64_t>(x.center));
                else {
                    integer(out, std::bit_cast<std::uint64_t>(x.high));
                    integer(out, std::bit_cast<std::uint64_t>(x.low));
                }
                integer(out, std::bit_cast<std::uint64_t>(x.radius));
                integer(out, static_cast<std::uint64_t>(x.exponent2));
            },
            quantity->enclosure);
    } else if constexpr (kind == io::detail::FieldKind::Integer) {
        if constexpr (io::detail::field_is_signed[static_cast<std::size_t>(Type)])
            integer(out, static_cast<std::uint64_t>(value.signed_value()));
        else
            integer(out, value.unsigned_value());
    } else if constexpr (kind == io::detail::FieldKind::Identifier ||
                         kind == io::detail::FieldKind::Reference)
        integer(out, value.unsigned_value());
    else if constexpr (kind == io::detail::FieldKind::Text)
        text(out, value.text());
    else if constexpr (kind == io::detail::FieldKind::Enumeration) {
        text(out, value.text());
        integer(out, static_cast<std::uint64_t>(value.selected()));
    } else if constexpr (kind == io::detail::FieldKind::Boolean)
        integer(out, value.boolean() ? 1U : 0U);
    else if constexpr (kind == io::detail::FieldKind::Array) {
        integer(out, value.size());
        for (std::uint64_t i = 0; i < value.size(); ++i) {
            const auto element = value.at(static_cast<core::Index>(i));
            if (!element)
                throw std::runtime_error("Checked schema array became invalid");
            append(out, *element);
        }
    } else
        append(out, value.record());
}

inline void append(Bytes& out, const io::schema::Case& value) {
    text(out, "Case");
    text(out, "schema");
    append(out, value.schema);
    text(out, "units");
    append(out, value.units);
    text(out, "model");
    append(out, value.model);
    text(out, "flow");
    append(out, value.flow);
    text(out, "analysis");
    append(out, value.analysis);
    text(out, "observations");
    append(out, value.observations);
    text(out, "numerics");
    append(out, value.numerics);
    text(out, "resources");
    append(out, value.resources);
    text(out, "output");
    append(out, value.output);
    text(out, "initial");
    append(out, value.initial);
    text(out, "restart");
    append(out, value.restart);
    text(out, "provenance");
    append(out, value.provenance);
}
inline void append(Bytes& out, const io::schema::Units& value) {
    text(out, "Units");
    text(out, "length");
    append(out, value.length);
    text(out, "mass");
    append(out, value.mass);
    text(out, "time");
    append(out, value.time);
    text(out, "angle");
    append(out, value.angle);
}
inline void append(Bytes& out, const io::schema::Model& value) {
    text(out, "Model");
    text(out, "geometry");
    append(out, value.geometry);
    text(out, "materials");
    append(out, value.materials);
    text(out, "shells");
    append(out, value.shells);
    text(out, "beams");
    append(out, value.beams);
    text(out, "joints");
    append(out, value.joints);
    text(out, "supports");
    append(out, value.supports);
    text(out, "masses");
    append(out, value.masses);
    text(out, "aero_surfaces");
    append(out, value.aero_surfaces);
}
inline void append(Bytes& out, const io::schema::Provenance& value) {
    text(out, "Provenance");
    text(out, "description");
    append(out, value.description);
    text(out, "references");
    append(out, value.references);
    text(out, "assumptions");
    append(out, value.assumptions);
}
inline void append(Bytes& out, const io::schema::Frame& value) {
    text(out, "Frame");
    text(out, "origin");
    append(out, value.origin);
    text(out, "rotation");
    append(out, value.rotation);
}
inline void append(Bytes& out, const io::schema::Target& value) {
    text(out, "Target");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
}
inline void append(Bytes& out, const io::schema::MaterialPoint& value) {
    text(out, "MaterialPoint");
    text(out, "target");
    append(out, value.target);
    text(out, "coordinates");
    append(out, value.coordinates);
    text(out, "side");
    append(out, value.side);
}
inline void append(Bytes& out, const io::schema::Region& value) {
    text(out, "Region");
    text(out, "target");
    append(out, value.target);
    text(out, "kind");
    append(out, value.kind);
    text(out, "geometry");
    append(out, value.geometry);
    text(out, "boundary");
    append(out, value.boundary);
    text(out, "side");
    append(out, value.side);
}
inline void append(Bytes& out, const io::schema::PolyTerm& value) {
    text(out, "PolyTerm");
    text(out, "powers");
    append(out, value.powers);
    text(out, "coefficient");
    append(out, value.coefficient);
}
inline void append(Bytes& out, const io::schema::ReferenceFieldConstant& value) {
    text(out, "ReferenceFieldConstant");
    text(out, "kind");
    append(out, value.kind);
    text(out, "value");
    append(out, value.value);
}
inline void append(Bytes& out, const io::schema::ReferenceFieldPolynomial& value) {
    text(out, "ReferenceFieldPolynomial");
    text(out, "kind");
    append(out, value.kind);
    text(out, "coordinates");
    append(out, value.coordinates);
    text(out, "origin");
    append(out, value.origin);
    text(out, "lengths");
    append(out, value.lengths);
    text(out, "components");
    append(out, value.components);
    text(out, "error");
    append(out, value.error);
}
inline void append(Bytes& out, const io::schema::TimeSeriesConstant& value) {
    text(out, "TimeSeriesConstant");
    text(out, "kind");
    append(out, value.kind);
    text(out, "value");
    append(out, value.value);
}
inline void append(Bytes& out, const io::schema::TimeSeriesPolynomial& value) {
    text(out, "TimeSeriesPolynomial");
    text(out, "kind");
    append(out, value.kind);
    text(out, "origin");
    append(out, value.origin);
    text(out, "scale");
    append(out, value.scale);
    text(out, "coefficients");
    append(out, value.coefficients);
    text(out, "domain");
    append(out, value.domain);
}
inline void append(Bytes& out, const io::schema::TimeSeriesPchip& value) {
    text(out, "TimeSeriesPchip");
    text(out, "kind");
    append(out, value.kind);
    text(out, "times");
    append(out, value.times);
    text(out, "values");
    append(out, value.values);
    text(out, "outside");
    append(out, value.outside);
}
inline void append(Bytes& out, const io::schema::OutsideRejectHold& value) {
    text(out, "OutsideRejectHold");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::OutsidePeriodic& value) {
    text(out, "OutsidePeriodic");
    text(out, "kind");
    append(out, value.kind);
    text(out, "period");
    append(out, value.period);
}
inline void append(Bytes& out, const io::schema::Motion& value) {
    text(out, "Motion");
    text(out, "origin");
    append(out, value.origin);
    text(out, "translation");
    append(out, value.translation);
    text(out, "rotation_vector");
    append(out, value.rotation_vector);
    text(out, "events");
    append(out, value.events);
}
inline void append(Bytes& out, const io::schema::AxisProjected& value) {
    text(out, "AxisProjected");
    text(out, "kind");
    append(out, value.kind);
    text(out, "direction");
    append(out, value.direction);
    text(out, "angle");
    append(out, value.angle);
}
inline void append(Bytes& out, const io::schema::AxisExplicit& value) {
    text(out, "AxisExplicit");
    text(out, "kind");
    append(out, value.kind);
    text(out, "columns");
    append(out, value.columns);
}
inline void append(Bytes& out, const io::schema::GeometryPg& value) {
    text(out, "GeometryPg");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "frame");
    append(out, value.frame);
    text(out, "reference_length");
    append(out, value.reference_length);
    text(out, "x0");
    append(out, value.x0);
    text(out, "width");
    append(out, value.width);
    text(out, "skew");
    append(out, value.skew);
    text(out, "degree");
    append(out, value.degree);
    text(out, "coefficients");
    append(out, value.coefficients);
    text(out, "domain");
    append(out, value.domain);
    text(out, "thickness");
    append(out, value.thickness);
    text(out, "n");
    append(out, value.n);
    text(out, "convention");
    append(out, value.convention);
    text(out, "boundaries");
    append(out, value.boundaries);
    text(out, "continuity");
    append(out, value.continuity);
}
inline void append(Bytes& out, const io::schema::PgTerm& value) {
    text(out, "PgTerm");
    text(out, "powers");
    append(out, value.powers);
    text(out, "coefficient");
    append(out, value.coefficient);
}
inline void append(Bytes& out, const io::schema::Thickness& value) {
    text(out, "Thickness");
    text(out, "vertices");
    append(out, value.vertices);
    text(out, "triangles");
    append(out, value.triangles);
    text(out, "values");
    append(out, value.values);
}
inline void append(Bytes& out, const io::schema::Boundary& value) {
    text(out, "Boundary");
    text(out, "name");
    append(out, value.name);
    text(out, "vertices");
    append(out, value.vertices);
}
inline void append(Bytes& out, const io::schema::Continuity& value) {
    text(out, "Continuity");
    text(out, "other");
    append(out, value.other);
    text(out, "kind");
    append(out, value.kind);
    text(out, "own_boundary");
    append(out, value.own_boundary);
    text(out, "other_boundary");
    append(out, value.other_boundary);
    text(out, "orientation");
    append(out, value.orientation);
}
inline void append(Bytes& out, const io::schema::GeometrySurfaceGraph& value) {
    text(out, "GeometrySurfaceGraph");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "frame");
    append(out, value.frame);
    text(out, "reference_length");
    append(out, value.reference_length);
    text(out, "x0");
    append(out, value.x0);
    text(out, "width");
    append(out, value.width);
    text(out, "skew");
    append(out, value.skew);
    text(out, "degree");
    append(out, value.degree);
    text(out, "coefficients");
    append(out, value.coefficients);
    text(out, "domain");
    append(out, value.domain);
    text(out, "boundaries");
    append(out, value.boundaries);
    text(out, "continuity");
    append(out, value.continuity);
}
inline void append(Bytes& out, const io::schema::GeometryNastran& value) {
    text(out, "GeometryNastran");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "deck");
    append(out, value.deck);
    text(out, "sidecar");
    append(out, value.sidecar);
    text(out, "dialect");
    append(out, value.dialect);
    text(out, "mode");
    append(out, value.mode);
    text(out, "excluded_cards");
    append(out, value.excluded_cards);
}
inline void append(Bytes& out, const io::schema::GeometryGmsh& value) {
    text(out, "GeometryGmsh");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "mesh");
    append(out, value.mesh);
    text(out, "sidecar");
    append(out, value.sidecar);
    text(out, "reference");
    append(out, value.reference);
}
inline void append(Bytes& out, const io::schema::ReferenceGeometrySupplied& value) {
    text(out, "ReferenceGeometrySupplied");
    text(out, "kind");
    append(out, value.kind);
    text(out, "file");
    append(out, value.file);
}
inline void append(Bytes& out, const io::schema::ReferenceGeometryFit& value) {
    text(out, "ReferenceGeometryFit");
    text(out, "kind");
    append(out, value.kind);
    text(out, "degree");
    append(out, value.degree);
    text(out, "distance_budget");
    append(out, value.distance_budget);
    text(out, "normal_budget");
    append(out, value.normal_budget);
    text(out, "conversion_observations");
    append(out, value.conversion_observations);
}
inline void append(Bytes& out, const io::schema::MaterialIsotropic& value) {
    text(out, "MaterialIsotropic");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "young");
    append(out, value.young);
    text(out, "poisson");
    append(out, value.poisson);
    text(out, "density");
    append(out, value.density);
    text(out, "provenance");
    append(out, value.provenance);
    text(out, "applicability");
    append(out, value.applicability);
}
inline void append(Bytes& out, const io::schema::MaterialOrthotropic& value) {
    text(out, "MaterialOrthotropic");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "young");
    append(out, value.young);
    text(out, "shear");
    append(out, value.shear);
    text(out, "poisson");
    append(out, value.poisson);
    text(out, "density");
    append(out, value.density);
    text(out, "provenance");
    append(out, value.provenance);
    text(out, "applicability");
    append(out, value.applicability);
}
inline void append(Bytes& out, const io::schema::MaterialEffective& value) {
    text(out, "MaterialEffective");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "stiffness");
    append(out, value.stiffness);
    text(out, "components");
    append(out, value.components);
    text(out, "density");
    append(out, value.density);
    text(out, "provenance");
    append(out, value.provenance);
    text(out, "applicability");
    append(out, value.applicability);
}
inline void append(Bytes& out, const io::schema::ApplicabilityUnevaluated& value) {
    text(out, "ApplicabilityUnevaluated");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::ApplicabilityDeclared& value) {
    text(out, "ApplicabilityDeclared");
    text(out, "kind");
    append(out, value.kind);
    text(out, "strain_lower");
    append(out, value.strain_lower);
    text(out, "strain_upper");
    append(out, value.strain_upper);
    text(out, "reference");
    append(out, value.reference);
}
inline void append(Bytes& out, const io::schema::Layer& value) {
    text(out, "Layer");
    text(out, "id");
    append(out, value.id);
    text(out, "interval");
    append(out, value.interval);
    text(out, "material");
    append(out, value.material);
    text(out, "axis");
    append(out, value.axis);
    text(out, "initial_green");
    append(out, value.initial_green);
    text(out, "viscosity");
    append(out, value.viscosity);
}
inline void append(Bytes& out, const io::schema::Shell& value) {
    text(out, "Shell");
    text(out, "id");
    append(out, value.id);
    text(out, "geometry");
    append(out, value.geometry);
    text(out, "layers");
    append(out, value.layers);
    text(out, "damping_components");
    append(out, value.damping_components);
    text(out, "loads");
    append(out, value.loads);
}
inline void append(Bytes& out, const io::schema::Beam& value) {
    text(out, "Beam");
    text(out, "id");
    append(out, value.id);
    text(out, "nodes");
    append(out, value.nodes);
    text(out, "coordinate");
    append(out, value.coordinate);
    text(out, "section");
    append(out, value.section);
    text(out, "initial_strain");
    append(out, value.initial_strain);
    text(out, "residual_energy");
    append(out, value.residual_energy);
    text(out, "viscosity");
    append(out, value.viscosity);
    text(out, "loads");
    append(out, value.loads);
}
inline void append(Bytes& out, const io::schema::BeamCoordinateArclength& value) {
    text(out, "BeamCoordinateArclength");
    text(out, "kind");
    append(out, value.kind);
    text(out, "length");
    append(out, value.length);
}
inline void append(Bytes& out, const io::schema::BeamCoordinateScaled& value) {
    text(out, "BeamCoordinateScaled");
    text(out, "kind");
    append(out, value.kind);
    text(out, "length");
    append(out, value.length);
    text(out, "arc_scale");
    append(out, value.arc_scale);
}
inline void append(Bytes& out, const io::schema::SectionMatrix& value) {
    text(out, "SectionMatrix");
    text(out, "kind");
    append(out, value.kind);
    text(out, "stiffness");
    append(out, value.stiffness);
    text(out, "line_density");
    append(out, value.line_density);
    text(out, "center");
    append(out, value.center);
    text(out, "inertia");
    append(out, value.inertia);
    text(out, "envelope");
    append(out, value.envelope);
    text(out, "warping");
    append(out, value.warping);
    text(out, "recovery");
    append(out, value.recovery);
    text(out, "source");
    append(out, value.source);
}
inline void append(Bytes& out, const io::schema::SectionCircle& value) {
    text(out, "SectionCircle");
    text(out, "kind");
    append(out, value.kind);
    text(out, "radius");
    append(out, value.radius);
    text(out, "young");
    append(out, value.young);
    text(out, "poisson");
    append(out, value.poisson);
    text(out, "density");
    append(out, value.density);
    text(out, "effective_shear");
    append(out, value.effective_shear);
    text(out, "warping");
    append(out, value.warping);
}
inline void append(Bytes& out, const io::schema::SectionRectangle& value) {
    text(out, "SectionRectangle");
    text(out, "kind");
    append(out, value.kind);
    text(out, "sides");
    append(out, value.sides);
    text(out, "young");
    append(out, value.young);
    text(out, "poisson");
    append(out, value.poisson);
    text(out, "density");
    append(out, value.density);
    text(out, "effective_shear");
    append(out, value.effective_shear);
    text(out, "torsion");
    append(out, value.torsion);
    text(out, "warping");
    append(out, value.warping);
}
inline void append(Bytes& out, const io::schema::JointAdhesive& value) {
    text(out, "JointAdhesive");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "first");
    append(out, value.first);
    text(out, "second");
    append(out, value.second);
    text(out, "mapping");
    append(out, value.mapping);
    text(out, "split_at_attachment");
    append(out, value.split_at_attachment);
}
inline void append(Bytes& out, const io::schema::TraceMap& value) {
    text(out, "TraceMap");
    text(out, "first_cell");
    append(out, value.first_cell);
    text(out, "second_cell");
    append(out, value.second_cell);
    text(out, "domain");
    append(out, value.domain);
    text(out, "first");
    append(out, value.first);
    text(out, "second");
    append(out, value.second);
    text(out, "orientation");
    append(out, value.orientation);
}
inline void append(Bytes& out, const io::schema::MapTerm& value) {
    text(out, "MapTerm");
    text(out, "powers");
    append(out, value.powers);
    text(out, "coefficient");
    append(out, value.coefficient);
}
inline void append(Bytes& out, const io::schema::JointRigidAttachment& value) {
    text(out, "JointRigidAttachment");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "attached");
    append(out, value.attached);
    text(out, "frame");
    append(out, value.frame);
    text(out, "offset");
    append(out, value.offset);
}
inline void append(Bytes& out, const io::schema::JointWeld& value) {
    text(out, "JointWeld");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "first");
    append(out, value.first);
    text(out, "second");
    append(out, value.second);
    text(out, "relative");
    append(out, value.relative);
}
inline void append(Bytes& out, const io::schema::SupportConstraint& value) {
    text(out, "SupportConstraint");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "target");
    append(out, value.target);
    text(out, "components");
    append(out, value.components);
    text(out, "motion");
    append(out, value.motion);
}
inline void append(Bytes& out, const io::schema::SupportTraceConstraint& value) {
    text(out, "SupportTraceConstraint");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "region");
    append(out, value.region);
    text(out, "components");
    append(out, value.components);
    text(out, "motion");
    append(out, value.motion);
}
inline void append(Bytes& out, const io::schema::SupportElastic& value) {
    text(out, "SupportElastic");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "target");
    append(out, value.target);
    text(out, "reference");
    append(out, value.reference);
    text(out, "stiffness");
    append(out, value.stiffness);
    text(out, "damping");
    append(out, value.damping);
    text(out, "frame");
    append(out, value.frame);
}
inline void append(Bytes& out, const io::schema::MassParticle& value) {
    text(out, "MassParticle");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "target");
    append(out, value.target);
    text(out, "mass");
    append(out, value.mass);
    text(out, "ledger");
    append(out, value.ledger);
}
inline void append(Bytes& out, const io::schema::MassRigid& value) {
    text(out, "MassRigid");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "target");
    append(out, value.target);
    text(out, "mass");
    append(out, value.mass);
    text(out, "center");
    append(out, value.center);
    text(out, "inertia");
    append(out, value.inertia);
    text(out, "frame");
    append(out, value.frame);
    text(out, "ledger");
    append(out, value.ledger);
}
inline void append(Bytes& out, const io::schema::LedgerAdd& value) {
    text(out, "LedgerAdd");
    text(out, "kind");
    append(out, value.kind);
    text(out, "source");
    append(out, value.source);
    text(out, "owner");
    append(out, value.owner);
}
inline void append(Bytes& out, const io::schema::LedgerReplace& value) {
    text(out, "LedgerReplace");
    text(out, "kind");
    append(out, value.kind);
    text(out, "source");
    append(out, value.source);
    text(out, "owner");
    append(out, value.owner);
    text(out, "replaced");
    append(out, value.replaced);
}
inline void append(Bytes& out, const io::schema::LoadPointForce& value) {
    text(out, "LoadPointForce");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "target");
    append(out, value.target);
    text(out, "frame");
    append(out, value.frame);
    text(out, "series");
    append(out, value.series);
}
inline void append(Bytes& out, const io::schema::LoadCouple& value) {
    text(out, "LoadCouple");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "target");
    append(out, value.target);
    text(out, "frame");
    append(out, value.frame);
    text(out, "series");
    append(out, value.series);
}
inline void append(Bytes& out, const io::schema::LoadGravity& value) {
    text(out, "LoadGravity");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "acceleration");
    append(out, value.acceleration);
}
inline void append(Bytes& out, const io::schema::LoadTraction& value) {
    text(out, "LoadTraction");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "region");
    append(out, value.region);
    text(out, "area");
    append(out, value.area);
    text(out, "frame");
    append(out, value.frame);
    text(out, "spatial");
    append(out, value.spatial);
    text(out, "multiplier");
    append(out, value.multiplier);
}
inline void append(Bytes& out, const io::schema::LoadPressure& value) {
    text(out, "LoadPressure");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "region");
    append(out, value.region);
    text(out, "spatial");
    append(out, value.spatial);
    text(out, "multiplier");
    append(out, value.multiplier);
}
inline void append(Bytes& out, const io::schema::LoadBody& value) {
    text(out, "LoadBody");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "region");
    append(out, value.region);
    text(out, "spatial");
    append(out, value.spatial);
    text(out, "multiplier");
    append(out, value.multiplier);
}
inline void append(Bytes& out, const io::schema::AeroSurface& value) {
    text(out, "AeroSurface");
    text(out, "id");
    append(out, value.id);
    text(out, "geometry");
    append(out, value.geometry);
    text(out, "attachment");
    append(out, value.attachment);
    text(out, "chord");
    append(out, value.chord);
    text(out, "span");
    append(out, value.span);
    text(out, "leading");
    append(out, value.leading);
    text(out, "trailing");
    append(out, value.trailing);
    text(out, "orientation");
    append(out, value.orientation);
    text(out, "sharing");
    append(out, value.sharing);
}
inline void append(Bytes& out, const io::schema::AeroAttachmentFixed& value) {
    text(out, "AeroAttachmentFixed");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::AeroAttachmentPrescribed& value) {
    text(out, "AeroAttachmentPrescribed");
    text(out, "kind");
    append(out, value.kind);
    text(out, "motion");
    append(out, value.motion);
}
inline void append(Bytes& out, const io::schema::AeroAttachmentStructural& value) {
    text(out, "AeroAttachmentStructural");
    text(out, "kind");
    append(out, value.kind);
    text(out, "host");
    append(out, value.host);
}
inline void append(Bytes& out, const io::schema::AeroSharing& value) {
    text(out, "AeroSharing");
    text(out, "other");
    append(out, value.other);
    text(out, "own_boundary");
    append(out, value.own_boundary);
    text(out, "other_boundary");
    append(out, value.other_boundary);
    text(out, "orientation");
    append(out, value.orientation);
}
inline void append(Bytes& out, const io::schema::FlowNone& value) {
    text(out, "FlowNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::FlowKnown& value) {
    text(out, "FlowKnown");
    text(out, "kind");
    append(out, value.kind);
    text(out, "density");
    append(out, value.density);
    text(out, "uniform");
    append(out, value.uniform);
    text(out, "gust");
    append(out, value.gust);
    text(out, "reference");
    append(out, value.reference);
    text(out, "events");
    append(out, value.events);
}
inline void append(Bytes& out, const io::schema::GustNone& value) {
    text(out, "GustNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::GustSine& value) {
    text(out, "GustSine");
    text(out, "kind");
    append(out, value.kind);
    text(out, "amplitude");
    append(out, value.amplitude);
    text(out, "direction");
    append(out, value.direction);
    text(out, "wavevector");
    append(out, value.wavevector);
    text(out, "angular_frequency");
    append(out, value.angular_frequency);
    text(out, "phase");
    append(out, value.phase);
    text(out, "origin");
    append(out, value.origin);
    text(out, "time_origin");
    append(out, value.time_origin);
}
inline void append(Bytes& out, const io::schema::FlowReference& value) {
    text(out, "FlowReference");
    text(out, "density");
    append(out, value.density);
    text(out, "speed");
    append(out, value.speed);
    text(out, "area");
    append(out, value.area);
    text(out, "chord");
    append(out, value.chord);
    text(out, "span");
    append(out, value.span);
    text(out, "moment_origin");
    append(out, value.moment_origin);
}
inline void append(Bytes& out, const io::schema::WakeNone& value) {
    text(out, "WakeNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::WakeActive& value) {
    text(out, "WakeActive");
    text(out, "kind");
    append(out, value.kind);
    text(out, "convection");
    append(out, value.convection);
    text(out, "retention");
    append(out, value.retention);
    text(out, "core_radius");
    append(out, value.core_radius);
    text(out, "model");
    append(out, value.model);
    text(out, "mapping");
    append(out, value.mapping);
}
inline void append(Bytes& out, const io::schema::ConvectionFree& value) {
    text(out, "ConvectionFree");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::ConvectionPrescribed& value) {
    text(out, "ConvectionPrescribed");
    text(out, "kind");
    append(out, value.kind);
    text(out, "velocity");
    append(out, value.velocity);
}
inline void append(Bytes& out, const io::schema::RetentionAll& value) {
    text(out, "RetentionAll");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::RetentionWindow& value) {
    text(out, "RetentionWindow");
    text(out, "kind");
    append(out, value.kind);
    text(out, "duration");
    append(out, value.duration);
    text(out, "distance");
    append(out, value.distance);
    text(out, "error_budget");
    append(out, value.error_budget);
}
inline void append(Bytes& out, const io::schema::AnalysisStatic& value) {
    text(out, "AnalysisStatic");
    text(out, "kind");
    append(out, value.kind);
    text(out, "physics");
    append(out, value.physics);
    text(out, "branch");
    append(out, value.branch);
    text(out, "speed");
    append(out, value.speed);
    text(out, "wake");
    append(out, value.wake);
    text(out, "loads");
    append(out, value.loads);
}
inline void append(Bytes& out, const io::schema::AnalysisTransient& value) {
    text(out, "AnalysisTransient");
    text(out, "kind");
    append(out, value.kind);
    text(out, "physics");
    append(out, value.physics);
    text(out, "branch");
    append(out, value.branch);
    text(out, "intervals");
    append(out, value.intervals);
    text(out, "wake");
    append(out, value.wake);
    text(out, "loads");
    append(out, value.loads);
}
inline void append(Bytes& out, const io::schema::AnalysisFlutter& value) {
    text(out, "AnalysisFlutter");
    text(out, "kind");
    append(out, value.kind);
    text(out, "physics");
    append(out, value.physics);
    text(out, "branches");
    append(out, value.branches);
    text(out, "intervals");
    append(out, value.intervals);
    text(out, "step");
    append(out, value.step);
    text(out, "wake");
    append(out, value.wake);
    text(out, "search");
    append(out, value.search);
    text(out, "loads");
    append(out, value.loads);
}
inline void append(Bytes& out, const io::schema::AnalysisSpeedNone& value) {
    text(out, "AnalysisSpeedNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::AnalysisSpeedValue& value) {
    text(out, "AnalysisSpeedValue");
    text(out, "kind");
    append(out, value.kind);
    text(out, "value");
    append(out, value.value);
}
inline void append(Bytes& out, const io::schema::Branch& value) {
    text(out, "Branch");
    text(out, "id");
    append(out, value.id);
    text(out, "seed");
    append(out, value.seed);
    text(out, "continuation");
    append(out, value.continuation);
}
inline void append(Bytes& out, const io::schema::BranchSeedReference& value) {
    text(out, "BranchSeedReference");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::BranchSeedPrepared& value) {
    text(out, "BranchSeedPrepared");
    text(out, "kind");
    append(out, value.kind);
    text(out, "history");
    append(out, value.history);
}
inline void append(Bytes& out, const io::schema::ContinuationFixed& value) {
    text(out, "ContinuationFixed");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::ContinuationSweep& value) {
    text(out, "ContinuationSweep");
    text(out, "kind");
    append(out, value.kind);
    text(out, "initial_increment");
    append(out, value.initial_increment);
    text(out, "minimum_increment");
    append(out, value.minimum_increment);
    text(out, "maximum_increment");
    append(out, value.maximum_increment);
    text(out, "arclength");
    append(out, value.arclength);
    text(out, "maximum_points");
    append(out, value.maximum_points);
    text(out, "maximum_retries");
    append(out, value.maximum_retries);
}
inline void append(Bytes& out, const io::schema::ArcLength& value) {
    text(out, "ArcLength");
    text(out, "state_scales");
    append(out, value.state_scales);
    text(out, "speed_scale");
    append(out, value.speed_scale);
    text(out, "orientation");
    append(out, value.orientation);
    text(out, "enable_at_limit");
    append(out, value.enable_at_limit);
}
inline void append(Bytes& out, const io::schema::TimeInterval& value) {
    text(out, "TimeInterval");
    text(out, "start");
    append(out, value.start);
    text(out, "end");
    append(out, value.end);
    text(out, "step");
    append(out, value.step);
    text(out, "restart");
    append(out, value.restart);
}
inline void append(Bytes& out, const io::schema::SpeedInterval& value) {
    text(out, "SpeedInterval");
    text(out, "lower");
    append(out, value.lower);
    text(out, "upper");
    append(out, value.upper);
}
inline void append(Bytes& out, const io::schema::Search& value) {
    text(out, "Search");
    text(out, "queue");
    append(out, value.queue);
    text(out, "speed_absolute");
    append(out, value.speed_absolute);
    text(out, "speed_relative");
    append(out, value.speed_relative);
    text(out, "growth_budget");
    append(out, value.growth_budget);
    text(out, "mode_association");
    append(out, value.mode_association);
    text(out, "maximum_points");
    append(out, value.maximum_points);
    text(out, "maximum_retries");
    append(out, value.maximum_retries);
    text(out, "coverage");
    append(out, value.coverage);
}
inline void append(Bytes& out, const io::schema::Association& value) {
    text(out, "Association");
    text(out, "physical_observations");
    append(out, value.physical_observations);
    text(out, "mac_minimum");
    append(out, value.mac_minimum);
    text(out, "subspace_angle_maximum");
    append(out, value.subspace_angle_maximum);
    text(out, "eigenvalue_absolute");
    append(out, value.eigenvalue_absolute);
    text(out, "eigenvalue_relative");
    append(out, value.eigenvalue_relative);
}
inline void append(Bytes& out, const io::schema::Coverage& value) {
    text(out, "Coverage");
    text(out, "branches");
    append(out, value.branches);
    text(out, "mode_observations");
    append(out, value.mode_observations);
    text(out, "require_physical");
    append(out, value.require_physical);
}
inline void append(Bytes& out, const io::schema::InitialEquilibrium& value) {
    text(out, "InitialEquilibrium");
    text(out, "kind");
    append(out, value.kind);
    text(out, "branch");
    append(out, value.branch);
    text(out, "perturbation");
    append(out, value.perturbation);
}
inline void append(Bytes& out, const io::schema::InitialCold& value) {
    text(out, "InitialCold");
    text(out, "kind");
    append(out, value.kind);
    text(out, "state");
    append(out, value.state);
    text(out, "circulation");
    append(out, value.circulation);
    text(out, "wake");
    append(out, value.wake);
}
inline void append(Bytes& out, const io::schema::InitialPrepared& value) {
    text(out, "InitialPrepared");
    text(out, "kind");
    append(out, value.kind);
    text(out, "file");
    append(out, value.file);
}
inline void append(Bytes& out, const io::schema::PerturbationNone& value) {
    text(out, "PerturbationNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::InitialState& value) {
    text(out, "InitialState");
    text(out, "shells");
    append(out, value.shells);
    text(out, "beams");
    append(out, value.beams);
    text(out, "rigids");
    append(out, value.rigids);
}
inline void append(Bytes& out, const io::schema::CoefficientState& value) {
    text(out, "CoefficientState");
    text(out, "shell");
    append(out, value.shell);
    text(out, "reference");
    append(out, value.reference);
    text(out, "representation");
    append(out, value.representation);
    text(out, "position");
    append(out, value.position);
    text(out, "velocity");
    append(out, value.velocity);
}
inline void append(Bytes& out, const io::schema::PlanesBinary64& value) {
    text(out, "PlanesBinary64");
    text(out, "kind");
    append(out, value.kind);
    text(out, "high");
    append(out, value.high);
}
inline void append(Bytes& out, const io::schema::PlanesTwoWord& value) {
    text(out, "PlanesTwoWord");
    text(out, "kind");
    append(out, value.kind);
    text(out, "high");
    append(out, value.high);
    text(out, "low");
    append(out, value.low);
}
inline void append(Bytes& out, const io::schema::BeamState& value) {
    text(out, "BeamState");
    text(out, "beam");
    append(out, value.beam);
    text(out, "nodes");
    append(out, value.nodes);
}
inline void append(Bytes& out, const io::schema::RigidState& value) {
    text(out, "RigidState");
    text(out, "rigid");
    append(out, value.rigid);
    text(out, "state");
    append(out, value.state);
}
inline void append(Bytes& out, const io::schema::FrameStatePayload& value) {
    text(out, "FrameStatePayload");
    text(out, "position");
    append(out, value.position);
    text(out, "rotation");
    append(out, value.rotation);
    text(out, "linear_velocity");
    append(out, value.linear_velocity);
    text(out, "angular_velocity");
    append(out, value.angular_velocity);
}
inline void append(Bytes& out, const io::schema::InitialWake& value) {
    text(out, "InitialWake");
    text(out, "surfaces");
    append(out, value.surfaces);
}
inline void append(Bytes& out, const io::schema::SurfaceWake& value) {
    text(out, "SurfaceWake");
    text(out, "surface");
    append(out, value.surface);
    text(out, "rows");
    append(out, value.rows);
    text(out, "points");
    append(out, value.points);
    text(out, "coefficients");
    append(out, value.coefficients);
    text(out, "row_ids");
    append(out, value.row_ids);
    text(out, "point_ids");
    append(out, value.point_ids);
    text(out, "band_ids");
    append(out, value.band_ids);
    text(out, "birth_times");
    append(out, value.birth_times);
}
inline void append(Bytes& out, const io::schema::Restart& value) {
    text(out, "Restart");
    text(out, "directory");
    append(out, value.directory);
    text(out, "manifest_sha256");
    append(out, value.manifest_sha256);
    text(out, "compatibility");
    append(out, value.compatibility);
    text(out, "budget");
    append(out, value.budget);
}
inline void append(Bytes& out, const io::schema::ObservationPoint& value) {
    text(out, "ObservationPoint");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "point");
    append(out, value.point);
    text(out, "quantity");
    append(out, value.quantity);
    text(out, "basis");
    append(out, value.basis);
    text(out, "comparison");
    append(out, value.comparison);
    text(out, "budget");
    append(out, value.budget);
}
inline void append(Bytes& out, const io::schema::ObservationStress& value) {
    text(out, "ObservationStress");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "point");
    append(out, value.point);
    text(out, "layer");
    append(out, value.layer);
    text(out, "interface");
    append(out, value.interface);
    text(out, "measure");
    append(out, value.measure);
    text(out, "contribution");
    append(out, value.contribution);
    text(out, "basis");
    append(out, value.basis);
    text(out, "comparison");
    append(out, value.comparison);
    text(out, "budget");
    append(out, value.budget);
}
inline void append(Bytes& out, const io::schema::ObservationSection& value) {
    text(out, "ObservationSection");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "region");
    append(out, value.region);
    text(out, "origin");
    append(out, value.origin);
    text(out, "basis");
    append(out, value.basis);
    text(out, "quantity");
    append(out, value.quantity);
    text(out, "comparison");
    append(out, value.comparison);
    text(out, "budget");
    append(out, value.budget);
}
inline void append(Bytes& out, const io::schema::ObservationReaction& value) {
    text(out, "ObservationReaction");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "support");
    append(out, value.support);
    text(out, "origin");
    append(out, value.origin);
    text(out, "basis");
    append(out, value.basis);
    text(out, "comparison");
    append(out, value.comparison);
    text(out, "budget");
    append(out, value.budget);
}
inline void append(Bytes& out, const io::schema::ObservationAero& value) {
    text(out, "ObservationAero");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "surfaces");
    append(out, value.surfaces);
    text(out, "origin");
    append(out, value.origin);
    text(out, "basis");
    append(out, value.basis);
    text(out, "quantity");
    append(out, value.quantity);
    text(out, "comparison");
    append(out, value.comparison);
    text(out, "budget");
    append(out, value.budget);
}
inline void append(Bytes& out, const io::schema::ObservationMode& value) {
    text(out, "ObservationMode");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "targets");
    append(out, value.targets);
    text(out, "points");
    append(out, value.points);
    text(out, "basis");
    append(out, value.basis);
    text(out, "quantity");
    append(out, value.quantity);
    text(out, "comparison");
    append(out, value.comparison);
    text(out, "budget");
    append(out, value.budget);
}
inline void append(Bytes& out, const io::schema::BasisWorldMaterial& value) {
    text(out, "BasisWorldMaterial");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::BasisExplicit& value) {
    text(out, "BasisExplicit");
    text(out, "kind");
    append(out, value.kind);
    text(out, "rotation");
    append(out, value.rotation);
}
inline void append(Bytes& out, const io::schema::Comparison& value) {
    text(out, "Comparison");
    text(out, "window");
    append(out, value.window);
    text(out, "statistic");
    append(out, value.statistic);
    text(out, "side");
    append(out, value.side);
    text(out, "reference");
    append(out, value.reference);
}
inline void append(Bytes& out, const io::schema::ComparisonReferenceNone& value) {
    text(out, "ComparisonReferenceNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::ComparisonReferenceFile& value) {
    text(out, "ComparisonReferenceFile");
    text(out, "kind");
    append(out, value.kind);
    text(out, "source");
    append(out, value.source);
    text(out, "columns");
    append(out, value.columns);
    text(out, "uncertainty");
    append(out, value.uncertainty);
}
inline void append(Bytes& out, const io::schema::ObservationBudget& value) {
    text(out, "ObservationBudget");
    text(out, "absolute");
    append(out, value.absolute);
    text(out, "relative");
    append(out, value.relative);
    text(out, "scale");
    append(out, value.scale);
    text(out, "groups");
    append(out, value.groups);
}
inline void append(Bytes& out, const io::schema::ErrorGroupNotApplicable& value) {
    text(out, "ErrorGroupNotApplicable");
    text(out, "status");
    append(out, value.status);
    text(out, "kind");
    append(out, value.kind);
    text(out, "reason");
    append(out, value.reason);
}
inline void append(Bytes& out, const io::schema::ErrorGroupCoveredBy& value) {
    text(out, "ErrorGroupCoveredBy");
    text(out, "status");
    append(out, value.status);
    text(out, "kind");
    append(out, value.kind);
    text(out, "representative");
    append(out, value.representative);
}
inline void append(Bytes& out, const io::schema::ErrorGroupRequired& value) {
    text(out, "ErrorGroupRequired");
    text(out, "status");
    append(out, value.status);
    text(out, "kind");
    append(out, value.kind);
    text(out, "covers");
    append(out, value.covers);
    text(out, "absolute");
    append(out, value.absolute);
    text(out, "relative");
    append(out, value.relative);
}
inline void append(Bytes& out, const io::schema::Numerics& value) {
    text(out, "Numerics");
    text(out, "patch");
    append(out, value.patch);
    text(out, "quadrature");
    append(out, value.quadrature);
    text(out, "linear");
    append(out, value.linear);
    text(out, "newton");
    append(out, value.newton);
    text(out, "precision");
    append(out, value.precision);
    text(out, "eigen");
    append(out, value.eigen);
    text(out, "wake_rows");
    append(out, value.wake_rows);
    text(out, "scales");
    append(out, value.scales);
    text(out, "error_groups");
    append(out, value.error_groups);
}
inline void append(Bytes& out, const io::schema::WakeCapacityNone& value) {
    text(out, "WakeCapacityNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::WakeCapacityLimited& value) {
    text(out, "WakeCapacityLimited");
    text(out, "kind");
    append(out, value.kind);
    text(out, "maximum_rows");
    append(out, value.maximum_rows);
}
inline void append(Bytes& out, const io::schema::PatchPolicyNone& value) {
    text(out, "PatchPolicyNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::PatchPolicyActive& value) {
    text(out, "PatchPolicyActive");
    text(out, "kind");
    append(out, value.kind);
    text(out, "maximum_degree");
    append(out, value.maximum_degree);
    text(out, "maximum_partition");
    append(out, value.maximum_partition);
    text(out, "geometry_error");
    append(out, value.geometry_error);
    text(out, "trace_error");
    append(out, value.trace_error);
    text(out, "rank");
    append(out, value.rank);
}
inline void append(Bytes& out, const io::schema::Quadrature& value) {
    text(out, "Quadrature");
    text(out, "structure");
    append(out, value.structure);
    text(out, "aero");
    append(out, value.aero);
}
inline void append(Bytes& out, const io::schema::StructuralQuadratureNone& value) {
    text(out, "StructuralQuadratureNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::StructuralQuadratureActive& value) {
    text(out, "StructuralQuadratureActive");
    text(out, "kind");
    append(out, value.kind);
    text(out, "surface_degree");
    append(out, value.surface_degree);
    text(out, "thickness_points");
    append(out, value.thickness_points);
    text(out, "independent_surface_degree");
    append(out, value.independent_surface_degree);
    text(out, "independent_thickness_points");
    append(out, value.independent_thickness_points);
    text(out, "absolute");
    append(out, value.absolute);
    text(out, "relative");
    append(out, value.relative);
    text(out, "maximum_subdivisions");
    append(out, value.maximum_subdivisions);
}
inline void append(Bytes& out, const io::schema::AeroQuadratureNone& value) {
    text(out, "AeroQuadratureNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::AeroQuadratureActive& value) {
    text(out, "AeroQuadratureActive");
    text(out, "kind");
    append(out, value.kind);
    text(out, "absolute");
    append(out, value.absolute);
    text(out, "relative");
    append(out, value.relative);
}
inline void append(Bytes& out, const io::schema::ScaleUniformComponents& value) {
    text(out, "ScaleUniformComponents");
    text(out, "kind");
    append(out, value.kind);
    text(out, "field");
    append(out, value.field);
    text(out, "targets");
    append(out, value.targets);
    text(out, "values");
    append(out, value.values);
}
inline void append(Bytes& out, const io::schema::ScalePhysicalLayoutPolicy& value) {
    text(out, "ScalePhysicalLayoutPolicy");
    text(out, "kind");
    append(out, value.kind);
    text(out, "field");
    append(out, value.field);
    text(out, "constraint");
    append(out, value.constraint);
}
inline void append(Bytes& out, const io::schema::ConstraintScalePolicyNone& value) {
    text(out, "ConstraintScalePolicyNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::ConstraintScalePolicyPhysical& value) {
    text(out, "ConstraintScalePolicyPhysical");
    text(out, "kind");
    append(out, value.kind);
    text(out, "position");
    append(out, value.position);
    text(out, "orientation");
    append(out, value.orientation);
    text(out, "multiplier_work_scale");
    append(out, value.multiplier_work_scale);
    text(out, "reaction_force");
    append(out, value.reaction_force);
    text(out, "reaction_moment");
    append(out, value.reaction_moment);
}
inline void append(Bytes& out, const io::schema::ComponentScale3None& value) {
    text(out, "ComponentScale3None");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::ComponentScale3Values& value) {
    text(out, "ComponentScale3Values");
    text(out, "kind");
    append(out, value.kind);
    text(out, "values");
    append(out, value.values);
}
inline void append(Bytes& out, const io::schema::ScaleExplicitLayout& value) {
    text(out, "ScaleExplicitLayout");
    text(out, "kind");
    append(out, value.kind);
    text(out, "field");
    append(out, value.field);
    text(out, "targets");
    append(out, value.targets);
    text(out, "layout");
    append(out, value.layout);
    text(out, "dimensions");
    append(out, value.dimensions);
    text(out, "values");
    append(out, value.values);
}
inline void append(Bytes& out, const io::schema::LinearPolicyDirect& value) {
    text(out, "LinearPolicyDirect");
    text(out, "backend");
    append(out, value.backend);
    text(out, "accuracy");
    append(out, value.accuracy);
    text(out, "maximum_refinements");
    append(out, value.maximum_refinements);
    text(out, "backend_memory");
    append(out, value.backend_memory);
}
inline void append(Bytes& out, const io::schema::LinearPolicyGmres& value) {
    text(out, "LinearPolicyGmres");
    text(out, "backend");
    append(out, value.backend);
    text(out, "accuracy");
    append(out, value.accuracy);
    text(out, "maximum_iterations");
    append(out, value.maximum_iterations);
    text(out, "restart");
    append(out, value.restart);
    text(out, "maximum_applications");
    append(out, value.maximum_applications);
    text(out, "maximum_adjoint");
    append(out, value.maximum_adjoint);
    text(out, "maximum_preconditioner");
    append(out, value.maximum_preconditioner);
    text(out, "maximum_true_residual_checks");
    append(out, value.maximum_true_residual_checks);
    text(out, "backend_memory");
    append(out, value.backend_memory);
}
inline void append(Bytes& out, const io::schema::LinearAccuracy& value) {
    text(out, "LinearAccuracy");
    text(out, "absolute");
    append(out, value.absolute);
    text(out, "relative");
    append(out, value.relative);
    text(out, "backward");
    append(out, value.backward);
    text(out, "norm");
    append(out, value.norm);
}
inline void append(Bytes& out, const io::schema::NewtonPolicy& value) {
    text(out, "NewtonPolicy");
    text(out, "maximum_iterations");
    append(out, value.maximum_iterations);
    text(out, "maximum_line_searches");
    append(out, value.maximum_line_searches);
    text(out, "absolute");
    append(out, value.absolute);
    text(out, "relative");
    append(out, value.relative);
    text(out, "correction");
    append(out, value.correction);
    text(out, "minimum_step");
    append(out, value.minimum_step);
}
inline void append(Bytes& out, const io::schema::PrecisionPolicy& value) {
    text(out, "PrecisionPolicy");
    text(out, "input");
    append(out, value.input);
    text(out, "retained");
    append(out, value.retained);
    text(out, "returned");
    append(out, value.returned);
    text(out, "projection_absolute");
    append(out, value.projection_absolute);
    text(out, "projection_relative");
    append(out, value.projection_relative);
}
inline void append(Bytes& out, const io::schema::BackendMemory& value) {
    text(out, "BackendMemory");
    text(out, "admission");
    append(out, value.admission);
    text(out, "maximum_bytes");
    append(out, value.maximum_bytes);
}
inline void append(Bytes& out, const io::schema::EigenPolicyNone& value) {
    text(out, "EigenPolicyNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::EigenPolicyFull& value) {
    text(out, "EigenPolicyFull");
    text(out, "kind");
    append(out, value.kind);
    text(out, "target");
    append(out, value.target);
    text(out, "eigenvalue_scale");
    append(out, value.eigenvalue_scale);
    text(out, "packet");
    append(out, value.packet);
    text(out, "requirements");
    append(out, value.requirements);
    text(out, "limits");
    append(out, value.limits);
}
inline void append(Bytes& out, const io::schema::EigenPolicyGeneralized& value) {
    text(out, "EigenPolicyGeneralized");
    text(out, "kind");
    append(out, value.kind);
    text(out, "target");
    append(out, value.target);
    text(out, "eigenvalue_scale");
    append(out, value.eigenvalue_scale);
    text(out, "requirements");
    append(out, value.requirements);
    text(out, "limits");
    append(out, value.limits);
}
inline void append(Bytes& out, const io::schema::EigenPolicyArnoldi& value) {
    text(out, "EigenPolicyArnoldi");
    text(out, "kind");
    append(out, value.kind);
    text(out, "target");
    append(out, value.target);
    text(out, "eigenvalue_scale");
    append(out, value.eigenvalue_scale);
    text(out, "packet");
    append(out, value.packet);
    text(out, "selector");
    append(out, value.selector);
    text(out, "nev");
    append(out, value.nev);
    text(out, "ncv");
    append(out, value.ncv);
    text(out, "seed");
    append(out, value.seed);
    text(out, "backend_tolerance");
    append(out, value.backend_tolerance);
    text(out, "requirements");
    append(out, value.requirements);
    text(out, "limits");
    append(out, value.limits);
}
inline void append(Bytes& out, const io::schema::Complex& value) {
    text(out, "Complex");
    text(out, "real");
    append(out, value.real);
    text(out, "imaginary");
    append(out, value.imaginary);
}
inline void append(Bytes& out, const io::schema::ArnoldiSeedUniform& value) {
    text(out, "ArnoldiSeedUniform");
    text(out, "kind");
    append(out, value.kind);
    text(out, "value");
    append(out, value.value);
}
inline void append(Bytes& out, const io::schema::ArnoldiSeedExplicit& value) {
    text(out, "ArnoldiSeedExplicit");
    text(out, "kind");
    append(out, value.kind);
    text(out, "layout");
    append(out, value.layout);
    text(out, "values");
    append(out, value.values);
}
inline void append(Bytes& out, const io::schema::EigenSelectorLargest& value) {
    text(out, "EigenSelectorLargest");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::EigenSelectorShift& value) {
    text(out, "EigenSelectorShift");
    text(out, "kind");
    append(out, value.kind);
    text(out, "shift");
    append(out, value.shift);
}
inline void append(Bytes& out, const io::schema::EigenRequirements& value) {
    text(out, "EigenRequirements");
    text(out, "residual");
    append(out, value.residual);
    text(out, "residual_observation");
    append(out, value.residual_observation);
    text(out, "normalization");
    append(out, value.normalization);
    text(out, "normalization_observation");
    append(out, value.normalization_observation);
    text(out, "orthogonality");
    append(out, value.orthogonality);
    text(out, "orthogonality_observation");
    append(out, value.orthogonality_observation);
    text(out, "backward");
    append(out, value.backward);
    text(out, "association_absolute");
    append(out, value.association_absolute);
    text(out, "association_relative");
    append(out, value.association_relative);
    text(out, "left");
    append(out, value.left);
    text(out, "clusters");
    append(out, value.clusters);
    text(out, "physical_coverage");
    append(out, value.physical_coverage);
}
inline void append(Bytes& out, const io::schema::BackwardNone& value) {
    text(out, "BackwardNone");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::BackwardRequired& value) {
    text(out, "BackwardRequired");
    text(out, "kind");
    append(out, value.kind);
    text(out, "budget");
    append(out, value.budget);
    text(out, "observation");
    append(out, value.observation);
}
inline void append(Bytes& out, const io::schema::EigenCluster& value) {
    text(out, "EigenCluster");
    text(out, "members");
    append(out, value.members);
    text(out, "expected_dimension");
    append(out, value.expected_dimension);
}
inline void append(Bytes& out, const io::schema::EigenLimits& value) {
    text(out, "EigenLimits");
    text(out, "maximum_order");
    append(out, value.maximum_order);
    text(out, "maximum_results");
    append(out, value.maximum_results);
    text(out, "maximum_ncv");
    append(out, value.maximum_ncv);
    text(out, "maximum_backend_iterations");
    append(out, value.maximum_backend_iterations);
    text(out, "maximum_attempts");
    append(out, value.maximum_attempts);
    text(out, "maximum_forward");
    append(out, value.maximum_forward);
    text(out, "maximum_adjoint");
    append(out, value.maximum_adjoint);
    text(out, "maximum_inverse");
    append(out, value.maximum_inverse);
    text(out, "backend_memory");
    append(out, value.backend_memory);
}
inline void append(Bytes& out, const io::schema::Resources& value) {
    text(out, "Resources");
    text(out, "profile");
    append(out, value.profile);
    text(out, "threads");
    append(out, value.threads);
    text(out, "maximum_memory_bytes");
    append(out, value.maximum_memory_bytes);
    text(out, "maximum_wall_seconds");
    append(out, value.maximum_wall_seconds);
    text(out, "maximum_steps");
    append(out, value.maximum_steps);
    text(out, "maximum_total_iterations");
    append(out, value.maximum_total_iterations);
    text(out, "checkpoint_temporary_bytes");
    append(out, value.checkpoint_temporary_bytes);
    text(out, "input");
    append(out, value.input);
}
inline void append(Bytes& out, const io::schema::InputQuotas& value) {
    text(out, "InputQuotas");
    text(out, "maximum_source_bytes");
    append(out, value.maximum_source_bytes);
    text(out, "maximum_depth");
    append(out, value.maximum_depth);
    text(out, "maximum_nodes");
    append(out, value.maximum_nodes);
    text(out, "maximum_members");
    append(out, value.maximum_members);
    text(out, "maximum_decoded_string_bytes");
    append(out, value.maximum_decoded_string_bytes);
    text(out, "maximum_decimal_token_bytes");
    append(out, value.maximum_decimal_token_bytes);
    text(out, "maximum_total_text_bytes");
    append(out, value.maximum_total_text_bytes);
    text(out, "maximum_typed_bytes");
    append(out, value.maximum_typed_bytes);
    text(out, "maximum_entities");
    append(out, value.maximum_entities);
    text(out, "maximum_references");
    append(out, value.maximum_references);
    text(out, "maximum_scalars");
    append(out, value.maximum_scalars);
    text(out, "maximum_polynomial_terms");
    append(out, value.maximum_polynomial_terms);
    text(out, "maximum_table_samples");
    append(out, value.maximum_table_samples);
    text(out, "maximum_expression_nodes");
    append(out, value.maximum_expression_nodes);
    text(out, "maximum_expression_edges");
    append(out, value.maximum_expression_edges);
    text(out, "maximum_expression_depth");
    append(out, value.maximum_expression_depth);
    text(out, "maximum_expression_power");
    append(out, value.maximum_expression_power);
    text(out, "maximum_expression_points");
    append(out, value.maximum_expression_points);
    text(out, "maximum_diagnostic_path_bytes");
    append(out, value.maximum_diagnostic_path_bytes);
}
inline void append(Bytes& out, const io::schema::Output& value) {
    text(out, "Output");
    text(out, "directory");
    append(out, value.directory);
    text(out, "observations");
    append(out, value.observations);
    text(out, "visualization");
    append(out, value.visualization);
    text(out, "checkpoint");
    append(out, value.checkpoint);
    text(out, "overwrite");
    append(out, value.overwrite);
    text(out, "resume");
    append(out, value.resume);
}
inline void append(Bytes& out, const io::schema::OutputCadence& value) {
    text(out, "OutputCadence");
    text(out, "every_steps");
    append(out, value.every_steps);
    text(out, "observations");
    append(out, value.observations);
}
inline void append(Bytes& out, const io::schema::Visualization& value) {
    text(out, "Visualization");
    text(out, "enabled");
    append(out, value.enabled);
    text(out, "every_steps");
    append(out, value.every_steps);
    text(out, "layers");
    append(out, value.layers);
    text(out, "interfaces");
    append(out, value.interfaces);
    text(out, "reference_geometry");
    append(out, value.reference_geometry);
    text(out, "current_geometry");
    append(out, value.current_geometry);
}
inline void append(Bytes& out, const io::schema::CheckpointOutputOff& value) {
    text(out, "CheckpointOutputOff");
    text(out, "kind");
    append(out, value.kind);
}
inline void append(Bytes& out, const io::schema::CheckpointOutputOn& value) {
    text(out, "CheckpointOutputOn");
    text(out, "kind");
    append(out, value.kind);
    text(out, "every_steps");
    append(out, value.every_steps);
    text(out, "directory");
    append(out, value.directory);
    text(out, "retain");
    append(out, value.retain);
}
inline void append(Bytes& out, const io::schema::ReferenceFieldExpression& value) {
    text(out, "ReferenceFieldExpression");
    text(out, "kind");
    append(out, value.kind);
    text(out, "expression");
    append(out, value.expression);
    text(out, "error");
    append(out, value.error);
}
inline void append(Bytes& out, const io::schema::TimeSeriesExpression& value) {
    text(out, "TimeSeriesExpression");
    text(out, "kind");
    append(out, value.kind);
    text(out, "expression");
    append(out, value.expression);
}
inline void append(Bytes& out, const io::schema::LoadCurrentExpression& value) {
    text(out, "LoadCurrentExpression");
    text(out, "kind");
    append(out, value.kind);
    text(out, "id");
    append(out, value.id);
    text(out, "region");
    append(out, value.region);
    text(out, "expression");
    append(out, value.expression);
}
inline void append(Bytes& out, const io::schema::Expression& value) {
    text(out, "Expression");
    text(out, "domain");
    append(out, value.domain);
    text(out, "nodes");
    append(out, value.nodes);
    text(out, "outputs");
    append(out, value.outputs);
}
inline void append(Bytes& out, const io::schema::ExpressionDomainReference& value) {
    text(out, "ExpressionDomainReference");
    text(out, "kind");
    append(out, value.kind);
    text(out, "coordinates");
    append(out, value.coordinates);
    text(out, "lower");
    append(out, value.lower);
    text(out, "upper");
    append(out, value.upper);
}
inline void append(Bytes& out, const io::schema::ExpressionDomainTime& value) {
    text(out, "ExpressionDomainTime");
    text(out, "kind");
    append(out, value.kind);
    text(out, "lower");
    append(out, value.lower);
    text(out, "upper");
    append(out, value.upper);
    text(out, "events");
    append(out, value.events);
}
inline void append(Bytes& out, const io::schema::ExpressionDomainCurrentPositionTime& value) {
    text(out, "ExpressionDomainCurrentPositionTime");
    text(out, "kind");
    append(out, value.kind);
    text(out, "position_lower");
    append(out, value.position_lower);
    text(out, "position_upper");
    append(out, value.position_upper);
    text(out, "time_lower");
    append(out, value.time_lower);
    text(out, "time_upper");
    append(out, value.time_upper);
    text(out, "events");
    append(out, value.events);
}
inline void append(Bytes& out, const io::schema::ExprNodeConstant& value) {
    text(out, "ExprNodeConstant");
    text(out, "op");
    append(out, value.op);
    text(out, "value");
    append(out, value.value);
    text(out, "dimensions");
    append(out, value.dimensions);
}
inline void append(Bytes& out, const io::schema::ExprNodeVariable& value) {
    text(out, "ExprNodeVariable");
    text(out, "op");
    append(out, value.op);
    text(out, "axis");
    append(out, value.axis);
}
inline void append(Bytes& out, const io::schema::ExprNodeBinary& value) {
    text(out, "ExprNodeBinary");
    text(out, "op");
    append(out, value.op);
    text(out, "first");
    append(out, value.first);
    text(out, "second");
    append(out, value.second);
}
inline void append(Bytes& out, const io::schema::ExprNodeNegate& value) {
    text(out, "ExprNodeNegate");
    text(out, "op");
    append(out, value.op);
    text(out, "operand");
    append(out, value.operand);
}
inline void append(Bytes& out, const io::schema::ExprNodePower& value) {
    text(out, "ExprNodePower");
    text(out, "op");
    append(out, value.op);
    text(out, "operand");
    append(out, value.operand);
    text(out, "exponent");
    append(out, value.exponent);
}
inline void append(Bytes& out, const io::schema::ExprNodeAnalytic& value) {
    text(out, "ExprNodeAnalytic");
    text(out, "op");
    append(out, value.op);
    text(out, "operand");
    append(out, value.operand);
}
inline void append(Bytes& out, const io::schema::PerturbationPhysical& value) {
    text(out, "PerturbationPhysical");
    text(out, "kind");
    append(out, value.kind);
    text(out, "shells");
    append(out, value.shells);
    text(out, "beams");
    append(out, value.beams);
    text(out, "rigids");
    append(out, value.rigids);
    text(out, "displacements");
    append(out, value.displacements);
    text(out, "velocities");
    append(out, value.velocities);
}
inline void append(Bytes& out, const io::schema::ShellPhysicalState& value) {
    text(out, "ShellPhysicalState");
    text(out, "shell");
    append(out, value.shell);
    text(out, "displacement");
    append(out, value.displacement);
    text(out, "velocity");
    append(out, value.velocity);
}
inline void append(Bytes& out, const io::schema::FullShellFieldGlobalPolynomial& value) {
    text(out, "FullShellFieldGlobalPolynomial");
    text(out, "kind");
    append(out, value.kind);
    text(out, "r");
    append(out, value.r);
    text(out, "p");
    append(out, value.p);
    text(out, "c");
    append(out, value.c);
}
inline void append(Bytes& out, const io::schema::FullShellFieldHalfCells& value) {
    text(out, "FullShellFieldHalfCells");
    text(out, "kind");
    append(out, value.kind);
    text(out, "cells");
    append(out, value.cells);
}
inline void append(Bytes& out, const io::schema::CellField& value) {
    text(out, "CellField");
    text(out, "cell");
    append(out, value.cell);
    text(out, "r");
    append(out, value.r);
    text(out, "p");
    append(out, value.p);
    text(out, "c");
    append(out, value.c);
}
inline void append(Bytes& out, const io::schema::VectorPolynomial& value) {
    text(out, "VectorPolynomial");
    text(out, "terms");
    append(out, value.terms);
}
inline void append(Bytes& out, const io::schema::VectorTerm& value) {
    text(out, "VectorTerm");
    text(out, "powers");
    append(out, value.powers);
    text(out, "value");
    append(out, value.value);
}
inline void append(Bytes& out, const io::schema::PointCheck& value) {
    text(out, "PointCheck");
    text(out, "point");
    append(out, value.point);
    text(out, "value");
    append(out, value.value);
    text(out, "tolerance");
    append(out, value.tolerance);
}
inline void append(Bytes& out, const io::schema::InitialPointShell& value) {
    text(out, "InitialPointShell");
    text(out, "target");
    append(out, value.target);
    text(out, "coordinates");
    append(out, value.coordinates);
    text(out, "side_s");
    append(out, value.side_s);
    text(out, "side_t");
    append(out, value.side_t);
}
inline void append(Bytes& out, const io::schema::InitialPointBeam& value) {
    text(out, "InitialPointBeam");
    text(out, "target");
    append(out, value.target);
    text(out, "coordinates");
    append(out, value.coordinates);
}
inline void append(Bytes& out, const io::schema::InitialPointRigid& value) {
    text(out, "InitialPointRigid");
    text(out, "target");
    append(out, value.target);
    text(out, "coordinates");
    append(out, value.coordinates);
}
inline void append(Bytes& out, const io::schema::Rational& value) {
    text(out, "Rational");
    text(out, "numerator");
    append(out, value.numerator);
    text(out, "denominator");
    append(out, value.denominator);
}
inline void append(Bytes& out, const io::schema::FileReference& value) {
    text(out, "FileReference");
    text(out, "source_id");
    append(out, value.source_id);
    text(out, "path");
    append(out, value.path);
    text(out, "sha256");
    append(out, value.sha256);
    text(out, "role");
    append(out, value.role);
}
inline Snapshot serialize(const io::SiCase& value) {
    const auto source = value.source_case();
    if (!source)
        throw std::runtime_error("SI owner lost its source case");
    const auto root = source->get().root();
    if (!root)
        throw std::runtime_error("SI source lost its root");
    Bytes out{value};
    const auto original = uvlas::test::schema_snapshot::serialize(*root);
    out.insert(out.end(), original.begin(), original.end());
    append(out, *root);
    const auto defaults = value.defaults();
    const auto memory = value.memory_report();
    if (!defaults || !memory)
        throw std::runtime_error("SI owner lost public metadata");
    integer(out, defaults->provenance_missing ? 1U : 0U);
    integer(out, defaults->initial_missing ? 1U : 0U);
    integer(out, defaults->restart_missing ? 1U : 0U);
    for (const auto count :
         {memory->owned_bytes, memory->scratch_bytes, memory->cumulative_bytes, memory->work,
          memory->planned_owned_bytes, memory->planned_scratch_bytes, memory->planned_work,
          memory->admission_bytes})
        integer(out, count);
    return {std::move(static_cast<std::vector<std::uint8_t>&>(out)), std::move(out.keys)};
}
} // namespace uvlas::test::si_ownership_snapshot
