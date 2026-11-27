#pragma once

#include "case_material_rows.hpp"

#include <uvlas/structure/pg_rigid_support_rows.hpp>

namespace uvlas::coupling::detail {

/// Unrounded conditioned shape row of a Case material vertex, retaining its
/// actual patch association. Its coefficients share the owned PG radicand.
class CaseMaterialAliasRow {
public:
    CaseMaterialAliasRow(CaseMaterialAliasRow&&) noexcept = default;
    CaseMaterialAliasRow& operator=(CaseMaterialAliasRow&&) noexcept = default;
    CaseMaterialAliasRow(const CaseMaterialAliasRow&) = delete;
    CaseMaterialAliasRow& operator=(const CaseMaterialAliasRow&) = delete;

    [[nodiscard]] static core::Result<CaseMaterialAliasRow>
    prepare(const CaseMaterialChart&, const CaseReferenceSurface&,
            const structure::PgPreparedPatch&, CaseMaterialVertexSite,
            const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] CaseMaterialChartIdentity surface() const noexcept {
        return surface_;
    }
    [[nodiscard]] structure::PgPreparedPatchIdentity patch() const noexcept {
        return patch_;
    }
    [[nodiscard]] const structure::ExactRational& radicand() const& noexcept {
        return radicand_;
    }
    [[nodiscard]] const structure::ExactRational& radicand() const&& = delete;
    [[nodiscard]] std::span<const structure::PgConditionedReferenceColumn>
    coefficients() const& noexcept {
        return row_.values();
    }
    [[nodiscard]] structure::PgReferenceCallInfo call_info() const noexcept {
        return info_;
    }
    [[nodiscard]] std::span<const structure::PgConditionedReferenceColumn>
    coefficients() const&& = delete;

private:
    CaseMaterialAliasRow(CaseMaterialChartIdentity, structure::PgPreparedPatchIdentity,
                         structure::ExactRational&&,
                         core::PlannedBuffer<structure::PgConditionedReferenceColumn>&&,
                         structure::PgReferenceCallInfo) noexcept;

    structure::PgReferenceCallInfo info_;
    CaseMaterialChartIdentity surface_;
    structure::PgPreparedPatchIdentity patch_;
    structure::ExactRational radicand_;
    core::PlannedBuffer<structure::PgConditionedReferenceColumn> row_;
};

/// Compare actual rows in the same fixed affine chart, including each support
/// group's translation and rotation generators. This establishes equality of
/// position and every time derivative of the prescribed term.
[[nodiscard]] core::Result<void> qualify_case_material_constraint_alias(
    const CaseMaterialAliasRow&, const CaseMaterialAliasRow&, const structure::PgRigidSupportRows&,
    const structure::PgReferenceLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

} // namespace uvlas::coupling::detail
