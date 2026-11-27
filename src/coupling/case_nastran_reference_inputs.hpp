#pragma once

#include "case_external_pg.hpp"
#include "io/nastran_beam_si.hpp"
#include "io/nastran_lumped_mass_geometry.hpp"
#include "io/nastran_rigid_attachment.hpp"

#include <cstdint>
#include <optional>
#include <span>
#include <variant>

namespace uvlas::coupling::detail {

struct CaseNastranReferenceInputLimits {
    CaseNastranPgSetLimits shells;
    io::detail::NastranBeamLimits beam_cards{};
    io::detail::NastranPbeamLimits beam_properties{};
    io::detail::NastranLumpedMassLimits mass_cards{};
    io::detail::NastranRigidLinkLimits rigid_cards{};
    io::detail::NastranBeamSiLimits beam;
    io::detail::NastranMassGeometryLimits mass;
    io::detail::NastranRigidAttachmentLimits rigid;
    std::uint64_t maximum_bound_cards{};
};

struct CaseNastranCardSelection {
    io::SchemaNamespace target_space;
    std::uint64_t case_id;
    io::detail::NastranElementId element;
    core::Index declaration;
};

struct CaseNastranReferenceDiagnostic {
    CaseNastranPgSetDiagnostic shells{};
    std::optional<CaseNastranCardSelection> selection{};
    std::variant<std::monostate, io::detail::NastranBeamDiagnostic,
                 io::detail::NastranPbeamDiagnostic, io::detail::NastranLumpedMassDiagnostic,
                 io::detail::NastranRigidLinkDiagnostic, io::detail::NastranBeamSiDiagnostic,
                 io::detail::NastranMassGeometryDiagnostic, io::detail::NastranAttachmentDiagnostic>
        card{};
};

template <io::SchemaNamespace Space, class Physical>
struct CaseNastranReferenceCard {
    io::detail::NastranSidecarElementBinding<Space> binding;
    Physical physical;
};

/// One Case and authenticated Nastran sidecar, with fully prepared planar Shells and owned
/// native beam, mass and rigid attachment reference inputs for the corresponding Case IDs.
/// The latter remain source properties for later Case equivalence checks and assembly.
class CaseNastranReferenceInputs {
public:
    using Beam =
        CaseNastranReferenceCard<io::SchemaNamespace::Beam, io::detail::NastranSiBeamInput>;
    using Mass =
        CaseNastranReferenceCard<io::SchemaNamespace::Mass, io::detail::NastranBasicMassGeometry>;
    using Joint =
        CaseNastranReferenceCard<io::SchemaNamespace::Joint, io::detail::NastranRigidAttachment>;

    [[nodiscard]] static core::Result<CaseNastranReferenceInputs>
    prepare(const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&,
            const io::detail::NastranLaminateCards&, const CasePgReferenceLimits&,
            const io::detail::ExternalSidecarPgLimits&, const CaseNastranReferenceInputLimits&,
            core::Workspace&, core::ResourcePlan&,
            CaseNastranReferenceDiagnostic* = nullptr) noexcept;

    [[nodiscard]] const std::optional<CaseNastranPlanarPgSet>& shells() const& noexcept {
        return shells_;
    }
    [[nodiscard]] const std::optional<CaseNastranPlanarPgSet>& shells() const&& = delete;

    // Each published slot is engaged. Arrays preserve sidecar declaration order.
    [[nodiscard]] std::span<const std::optional<Beam>> beams() const& noexcept {
        return beams_.values();
    }
    [[nodiscard]] std::span<const std::optional<Beam>> beams() const&& = delete;
    [[nodiscard]] std::span<const std::optional<Mass>> masses() const& noexcept {
        return masses_.values();
    }
    [[nodiscard]] std::span<const std::optional<Mass>> masses() const&& = delete;
    [[nodiscard]] std::span<const std::optional<Joint>> joints() const& noexcept {
        return joints_.values();
    }
    [[nodiscard]] std::span<const std::optional<Joint>> joints() const&& = delete;
    [[nodiscard]] io::detail::RootSnapshotLocation sidecar_source() const noexcept {
        return sidecar_;
    }
    [[nodiscard]] io::SourceDigest sidecar_digest() const noexcept {
        return sidecar_digest_;
    }
    [[nodiscard]] io::detail::RootSnapshotLocation deck_source() const noexcept {
        return deck_;
    }
    [[nodiscard]] io::SourceDigest deck_digest() const noexcept {
        return deck_digest_;
    }

private:
    CaseNastranReferenceInputs(io::detail::RootSnapshotLocation sidecar,
                               io::SourceDigest sidecar_digest,
                               io::detail::RootSnapshotLocation deck,
                               io::SourceDigest deck_digest) noexcept
        : sidecar_(sidecar),
          sidecar_digest_(sidecar_digest),
          deck_(deck),
          deck_digest_(deck_digest) {}

    io::detail::RootSnapshotLocation sidecar_;
    io::SourceDigest sidecar_digest_;
    io::detail::RootSnapshotLocation deck_;
    io::SourceDigest deck_digest_;
    std::optional<CaseNastranPlanarPgSet> shells_;
    core::PlannedBuffer<std::optional<Beam>> beams_;
    core::PlannedBuffer<std::optional<Mass>> masses_;
    core::PlannedBuffer<std::optional<Joint>> joints_;
};

using CaseNastranMassId = io::detail::MeshCaseId<io::SchemaNamespace::Mass>;

enum class CaseNastranParticleFailure : std::uint8_t {
    Source,
    Selection,
    Kind,
    Mass,
    Position,
    Frame,
    Centre,
    Inertia,
};

struct CaseNastranParticleDiagnostic {
    CaseNastranParticleFailure reason = CaseNastranParticleFailure::Source;
    std::optional<io::detail::NastranElementId> element;
    std::optional<io::detail::NastranBulkPosition> field;
};

struct CaseNastranMassAttachment {
    io::detail::RootSnapshotLocation deck, sidecar;
    io::SourceDigest deck_digest, sidecar_digest;
    io::detail::NastranSidecarElementBinding<io::SchemaNamespace::Mass> binding;
    io::detail::MeshShellId shell;
    structure::PgReferenceIdentity patch;
    io::detail::NastranGridId grid;
    std::uint64_t face;
    structure::ExactRational mass_kg;
};

using CaseNastranParticleMass = CaseNastranMassAttachment;

struct CaseNastranRigidMass {
    CaseNastranMassAttachment attachment;
    std::array<structure::ExactRational, 3> centre_metres;
    // Centroidal tensor in basic axes, lower triangle in row order.
    std::array<structure::ExactRational, 6> centre_inertia_kg_m2;
};

// A particle requires one authenticated CONM2, zero lever and centre inertia, and
// the same original shell GRID at the Case material point on the midsurface.
[[nodiscard]] core::Result<CaseNastranParticleMass>
certify_case_nastran_particle_mass(const io::SiCase&, const CaseNastranReferenceInputs&,
                                   CaseNastranMassId, const structure::PgReferenceLimits&,
                                   core::ResourcePlan&,
                                   CaseNastranParticleDiagnostic* = nullptr) noexcept;

/// Compare the same shell GRID, the centre p + R c and the centroidal tensor
/// R I_c R^T in basic coordinates, using an authenticated rigid CONM2.
[[nodiscard]] core::Result<CaseNastranRigidMass>
certify_case_nastran_rigid_mass(const io::SiCase&, const CaseNastranReferenceInputs&,
                                CaseNastranMassId, const structure::PgReferenceLimits&,
                                core::ResourcePlan&,
                                CaseNastranParticleDiagnostic* = nullptr) noexcept;

} // namespace uvlas::coupling::detail
