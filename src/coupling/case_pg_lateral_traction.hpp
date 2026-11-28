#pragma once

#include "case_pg_chart.hpp"

#include <uvlas/structure/pg_shell_material.hpp>

namespace uvlas::coupling::detail {

struct CasePgLateralTractionRule {
    std::span<const structure::PgLineQuadratureNode> edge, depth;
};

/// The two work limits in reference cover planning, points, root refinement and assembly.
struct CasePgLateralTractionLimits {
    structure::PgReferenceLimits reference;
    std::uint64_t maximum_loads;
    std::uint64_t maximum_atoms;
    std::uint64_t maximum_layers;
    std::uint64_t maximum_rule_nodes;
    std::uint64_t maximum_points;
};

struct CasePgLateralTractionAccuracy {
    std::uint32_t sum_fraction_bits;
    const structure::ExactRational& maximum_force_width_newtons;
};

struct CasePgLateralTractionDiagnostic {
    std::optional<std::uint64_t> load{};
    std::optional<std::uint64_t> boundary{};
    std::optional<core::Index> segment{}, atom{}, face{};
    std::optional<core::LayerId> layer{};
};

/// All indices address declared original or owned partition records.
struct CasePgLateralTractionPoint {
    std::uint64_t load{}, boundary{};
    core::Index segment{}, atom{}, face{};
    core::LayerId layer{0};
    structure::ExactRational source_along;  ///< Atom fraction in the original segment direction.
    structure::ExactRational face_along;    ///< Same position in the physical face direction.
    structure::ExactRational depth, weight; ///< Layer fraction and exact tensor weight.
};

struct CasePgLateralTractionForce {
    structure::PgPreparedPatchIdentity patch;
    numerics::PatchId shell;
    /// Positive applied force, original conditioned column by world XYZ. Subtract from residuals.
    /// This reference-world force has exactly zero derivative with respect to the original state.
    core::PlannedBuffer<structure::PgShellCartesianEnclosure> newtons;
    structure::PgReferenceCallInfo call_info;
};

/// FD-05 shell edge tractions in Case order with owned integration point provenance.
class CasePgLateralTraction {
public:
    CasePgLateralTraction(CasePgLateralTraction&& other) noexcept
        : source_(std::exchange(other.source_, std::nullopt)),
          shell_(other.shell_),
          patch_(other.patch_),
          loads_(std::move(other.loads_)),
          points_(std::move(other.points_)) {}

    [[nodiscard]] static core::Result<CasePgLateralTraction>
    prepare(const io::SiCase&, const CaseStaticPgPointPatch&, const CasePgLateralTractionRule&,
            const CasePgLateralTractionLimits&, core::ResourcePlan&,
            CasePgLateralTractionDiagnostic* = nullptr) noexcept;

    [[nodiscard]] core::Result<CasePgLateralTractionForce>
    evaluate(const CaseStaticPgPointPatch&, double time_seconds,
             const CasePgLateralTractionAccuracy&, const CasePgLateralTractionLimits&,
             core::ResourcePlan&, CasePgLateralTractionDiagnostic* = nullptr) const noexcept;

    [[nodiscard]] std::span<const CasePgLateralTractionPoint> points() const& noexcept {
        return points_.values();
    }
    [[nodiscard]] std::span<const CasePgLateralTractionPoint> points() const&& = delete;

private:
    struct Load {
        std::array<structure::ExactRational, 3> pascals;
        core::Index point_count{};
    };
    CasePgLateralTraction(io::SchemaIdentity source, numerics::PatchId shell,
                          structure::PgPreparedPatchIdentity patch,
                          core::PlannedBuffer<Load>&& loads,
                          core::PlannedBuffer<CasePgLateralTractionPoint>&& points) noexcept
        : source_(source),
          shell_(shell),
          patch_(patch),
          loads_(std::move(loads)),
          points_(std::move(points)) {}

    std::optional<io::SchemaIdentity> source_;
    numerics::PatchId shell_;
    structure::PgPreparedPatchIdentity patch_;
    core::PlannedBuffer<Load> loads_;
    core::PlannedBuffer<CasePgLateralTractionPoint> points_;
};

} // namespace uvlas::coupling::detail
