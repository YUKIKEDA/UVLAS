#pragma once

#include "case_gmsh_support_boundary.hpp"
#include "case_nastran_static_point_support.hpp"
#include "case_pg_chart.hpp"

namespace uvlas::coupling::detail {

struct CaseGmshPgChartLimits {
    CasePgChartLimits chart;
    io::detail::ExternalSidecarPgLimits correspondence;
    CaseGmshPgSetLimits shells;
    CaseGmshSupportBoundaryLimits support;
    std::uint64_t maximum_support_mappings{};
    std::uint64_t maximum_joint_mappings{};
};

struct CaseGmshPgChartDiagnostic {
    CaseGmshPgSetDiagnostic shell;
    std::optional<CaseGmshSupportBoundaryDiagnostic> support;
    std::optional<CaseGmshJointBoundaryDiagnostic> joint;
};

/// One Case chart and its imported evidence. The chart owns the same
/// qualified PG references used to compare the original mesh and support lines.
struct CaseGmshPgChart {
    CasePgChart chart;
    core::PlannedBuffer<std::optional<io::detail::PreparedGmshSidecarPg>> shells;
    core::PlannedBuffer<std::optional<CaseGmshSupportBoundary>> supports;
    core::PlannedBuffer<std::optional<CaseGmshJointBoundary>> joints;
};

/// Shell evidence follows sidecar first-declaration order and support evidence
/// follows the support and joint declarations. Parametric Shells retain the Case path.
[[nodiscard]] core::Result<CaseGmshPgChart>
prepare_case_gmsh_pg_chart(const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&,
                           std::span<const CasePgRuleBinding>, const CaseGmshPgChartLimits&,
                           core::Workspace&, core::ResourcePlan&, const std::stop_token& = {},
                           CaseGmshPgChartDiagnostic* = nullptr) noexcept;

struct CaseNastranPgChartLimits {
    CasePgChartLimits chart;
    io::detail::ExternalSidecarPgLimits correspondence;
    CaseNastranPgSetLimits shells;
};

struct CaseNastranPgChart {
    CasePgChart chart;
    core::PlannedBuffer<std::optional<CaseNastranPgMatch>> shells;
};

/// Use the authenticated planar correspondence, published laminate and axes
/// and the existing Case support and adhesive contracts for structural rows.
[[nodiscard]] core::Result<CaseNastranPgChart>
prepare_case_nastran_pg_chart(const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&,
                              const io::detail::NastranLaminateCards&,
                              std::span<const CasePgRuleBinding>, const CaseNastranPgChartLimits&,
                              core::Workspace&, core::ResourcePlan&, const std::stop_token& = {},
                              CaseNastranPgSetDiagnostic* = nullptr) noexcept;

struct CaseNastranPointPgChartLimits {
    CaseNastranPgChartLimits base;
    CaseNastranPointSupportLimits support;
    std::uint64_t maximum_support_mappings{};
};

struct CaseNastranPointPgChartDiagnostic {
    CaseNastranPgSetDiagnostic shell;
    std::optional<CaseNastranPointSupportDiagnostic> support;
};

/// Own the qualified Shells, the static point-support evidence in sidecar order,
/// and the chart built from the same conditioned bases.
struct CaseNastranPointPgChart {
    CaseComponentPgChart structural;
    core::PlannedBuffer<std::optional<CaseNastranPgMatch>> shells;
    core::PlannedBuffer<std::optional<CaseNastranPointSupport>> supports;
};

[[nodiscard]] core::Result<CaseNastranPointPgChart> prepare_case_nastran_pg_chart(
    const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&,
    const io::detail::NastranLaminateCards&, std::span<const CasePgRuleBinding>,
    const CaseNastranPointPgChartLimits&, core::Workspace&, core::ResourcePlan&,
    const std::stop_token& = {}, CaseNastranPointPgChartDiagnostic* = nullptr) noexcept;

} // namespace uvlas::coupling::detail
