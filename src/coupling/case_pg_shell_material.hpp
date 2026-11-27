#pragma once

#include <uvlas/io/si_case.hpp>
#include <uvlas/structure/pg_shell_material.hpp>

namespace uvlas::coupling::detail {

class CaseStaticPgPointPatch;

struct CasePgShellLayer {
    core::LayerId id{0};
    std::uint64_t material_id{};
    structure::ExactRational lower, upper;
    std::size_t law{};
};

/// Fixed coordinates of one reference point. Material and Cartesian
/// coordinates are SI lengths and use the same exact radical sqrt(g).
struct CasePgMaterialCoordinates {
    const std::array<structure::ExactRational, 3>& chart; ///< s, t, zeta.
    const structure::PgReferenceVector& material_metres;
    const structure::PgReferenceVector& world_reference_metres;
    const structure::ExactRational& g;
    bool has_material_metres = true;
};

struct CasePgInitialGreenValue {
    core::LayerId layer;
    structure::ExactRational g;
    std::array<structure::PgReferenceScalar, 6> mandel;
    std::array<structure::ExactRational, 6> declared_error;
    structure::PgReferenceCallInfo work;

    [[nodiscard]] structure::PgInitialGreenEnclosure enclosure() const& noexcept {
        return {g, mandel, declared_error};
    }
    [[nodiscard]] structure::PgInitialGreenEnclosure enclosure() const&& = delete;
};

struct CasePgProjectedAxisValue {
    structure::PgReferenceVector direction;
    std::array<structure::ExactRational, 3> direction_error;
    structure::PgMaterialEnclosure cosine, sine;
};

struct CasePgExplicitAxisValue {
    std::array<structure::PgReferenceScalar, 9> columns;
};

struct CasePgAxisValue {
    core::LayerId layer;
    std::variant<CasePgProjectedAxisValue, CasePgExplicitAxisValue> value;
    structure::PgReferenceCallInfo work;

    [[nodiscard]] structure::PgShellPointAxes axes() const& noexcept;
    [[nodiscard]] structure::PgShellPointAxes axes() const&& = delete;
};

/// Own the layer intervals and reusable normalized laws for a qualified Case
/// shell. Fields of the reference axis and initial strain still belong to the
/// typed original Case. Retrieve them only against the same source identity.
class CasePgShellMaterials {
public:
    CasePgShellMaterials(CasePgShellMaterials&&) noexcept = default;
    CasePgShellMaterials& operator=(CasePgShellMaterials&&) noexcept = default;
    CasePgShellMaterials(const CasePgShellMaterials&) = delete;
    CasePgShellMaterials& operator=(const CasePgShellMaterials&) = delete;

    [[nodiscard]] static core::Result<CasePgShellMaterials>
    prepare(const io::SiCase&, std::uint64_t shell_id, const structure::PgReferenceLimits&,
            core::ResourcePlan&) noexcept;

    [[nodiscard]] std::span<const CasePgShellLayer> layers() const& noexcept {
        return layers_.values();
    }
    [[nodiscard]] std::span<const CasePgShellLayer> layers() const&& = delete;
    [[nodiscard]] std::size_t law_count() const noexcept {
        return layers_.size() == 0 ? 0 : law_count_;
    }
    [[nodiscard]] structure::PgReferenceCallInfo call_info() const noexcept {
        return work_;
    }
    [[nodiscard]] core::Result<std::reference_wrapper<const structure::PgShellElasticLaw>>
    law_for(core::Index layer) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const structure::PgShellElasticLaw>>
        law_for(core::Index) const&& = delete;
    [[nodiscard]] core::Result<io::schema::Layer> source_layer(const io::SiCase&,
                                                               core::Index layer) const noexcept;
    [[nodiscard]] core::Result<io::schema::Layer> source_layer(const io::SiCase&&,
                                                               core::Index) const = delete;

    /// Evaluate the fixed initial field for this selected layer. At a shared
    /// interface, the caller selects the layer for the requested physical side.
    [[nodiscard]] core::Result<CasePgInitialGreenValue>
    sample_initial_green(const io::SiCase&, core::Index layer, const CasePgMaterialCoordinates&,
                         const structure::PgReferenceLimits&, core::ResourcePlan&) const noexcept;

    [[nodiscard]] core::Result<CasePgAxisValue> sample_axes(const io::SiCase&, core::Index layer,
                                                            const CasePgMaterialCoordinates&,
                                                            const structure::PgReferenceLimits&,
                                                            core::ResourcePlan&) const noexcept;

private:
    using Laws = core::PlannedBuffer<std::optional<structure::PgShellElasticLaw>>;
    CasePgShellMaterials(io::SchemaIdentity, core::Index shell,
                         core::PlannedBuffer<CasePgShellLayer>&&, Laws&&, std::size_t law_count,
                         structure::PgReferenceCallInfo) noexcept;

    io::SchemaIdentity source_;
    core::Index shell_;
    core::PlannedBuffer<CasePgShellLayer> layers_;
    Laws laws_;
    std::size_t law_count_;
    structure::PgReferenceCallInfo work_;
};

/// Own the rule and sampled fixed Case fields for one authenticated shell.
class CasePgShellVolume {
public:
    /// When used by a field, supply material coordinates in the rule's order.
    /// Chart and world coordinates come from the owning reference geometry.
    [[nodiscard]] static core::Result<CasePgShellVolume>
    prepare(const io::SiCase&, const CaseStaticPgPointPatch&,
            std::span<const structure::PgOwnedDiscreteVolumePoint>,
            std::span<const structure::PgReferenceVector> material_metres,
            const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<structure::PgShellElasticVolumeAction>
    elastic(const CaseStaticPgPointPatch&, const structure::PgShellPatchState&,
            const std::optional<structure::PgShellPatchDirection>&,
            const structure::PgShellVolumeResolution&, const structure::PgReferenceLimits&,
            core::ResourcePlan&) const noexcept;

private:
    struct Sample {
        structure::PgOwnedDiscreteVolumePoint quadrature;
        core::Index layer{};
        std::optional<CasePgInitialGreenValue> initial;
        std::optional<CasePgAxisValue> axes;
    };
    CasePgShellVolume(io::SchemaIdentity, structure::PgPreparedPatchIdentity, std::uint64_t,
                      CasePgShellMaterials&&, core::PlannedBuffer<Sample>&&) noexcept;

    io::SchemaIdentity source_;
    structure::PgPreparedPatchIdentity patch_;
    std::uint64_t shell_;
    CasePgShellMaterials materials_;
    core::PlannedBuffer<Sample> samples_;
};

} // namespace uvlas::coupling::detail
