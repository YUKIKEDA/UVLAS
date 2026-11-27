#pragma once

#include <uvlas/io/si_case.hpp>
#include <uvlas/structure/pg_shell_material.hpp>

namespace uvlas::coupling::detail {

struct CasePgShellLayer {
    core::LayerId id{0};
    std::uint64_t material_id{};
    structure::ExactRational lower, upper;
    std::size_t law{};
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
    [[nodiscard]] core::Result<std::reference_wrapper<const structure::PgShellElasticLaw>>
    law_for(core::Index layer) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const structure::PgShellElasticLaw>>
        law_for(core::Index) const&& = delete;
    [[nodiscard]] core::Result<io::schema::Layer> source_layer(const io::SiCase&,
                                                               core::Index layer) const noexcept;
    [[nodiscard]] core::Result<io::schema::Layer> source_layer(const io::SiCase&&,
                                                               core::Index) const = delete;

private:
    using Laws = core::PlannedBuffer<std::optional<structure::PgShellElasticLaw>>;
    CasePgShellMaterials(io::SchemaIdentity, core::Index shell,
                         core::PlannedBuffer<CasePgShellLayer>&&, Laws&&,
                         std::size_t law_count) noexcept;

    io::SchemaIdentity source_;
    core::Index shell_;
    core::PlannedBuffer<CasePgShellLayer> layers_;
    Laws laws_;
    std::size_t law_count_;
};

} // namespace uvlas::coupling::detail
