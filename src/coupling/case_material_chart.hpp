#pragma once

#include "case_pg_reference.hpp"

#include <array>

namespace uvlas::coupling::detail {

enum class CaseMaterialSide : std::uint8_t { Interior, Lower, Upper };

struct CaseMaterialChartLimits {
    CasePgReferenceLimits support;
    structure::PgReferenceLimits chart;
    std::uint64_t maximum_subdivisions; ///< Shared across the aerodynamic surface.
};

struct CaseMaterialChartIdentity {
    io::SchemaIdentity source;
    std::uint64_t surface_id;
    std::uint64_t shell_id;
    std::uint64_t geometry_id;
    bool reversed;
};

struct CaseMaterialChartPoint {
    std::array<structure::ExactRational, 2> aerodynamic;
    std::array<structure::ExactRational, 3> material, chord_tangent, span_tangent;
    core::LayerId layer{0};
    CaseMaterialSide side{};
};

/// Own the fixed Case material map and the physical reference selected from the same Case.
class CaseMaterialChart {
public:
    CaseMaterialChart(CaseMaterialChart&&) noexcept;
    CaseMaterialChart& operator=(CaseMaterialChart&&) noexcept;
    CaseMaterialChart(const CaseMaterialChart&) = delete;
    CaseMaterialChart& operator=(const CaseMaterialChart&) = delete;
    ~CaseMaterialChart();

    [[nodiscard]] static core::Result<CaseMaterialChart>
    prepare(const io::SiCase&, std::uint64_t surface_id, const CaseMaterialChartLimits&,
            core::Workspace&, core::ResourcePlan&) noexcept;
    [[nodiscard]] core::Result<CaseMaterialChartIdentity> identity() const noexcept;
    [[nodiscard]] core::Result<const structure::PgReferenceGeometry*> support() const& noexcept;
    [[nodiscard]] core::Result<const structure::PgReferenceGeometry*> support() const&& = delete;

    /// Coordinates at three-quarter chord and midspan of the selected physical panel.
    [[nodiscard]] core::Result<CaseMaterialChartPoint>
    collocation(core::Index chord_cell, core::Index span_cell, const structure::PgReferenceLimits&,
                core::ResourcePlan&) const noexcept;

    /// Quarter-chord vertex column or final trailing edge. An incident owning span panel is
    /// explicit, so each load contribution can retain its own original material side.
    [[nodiscard]] core::Result<CaseMaterialChartPoint>
    vertex(core::Index chord_column, core::Index span_knot, core::Index owning_span_cell,
           const structure::PgReferenceLimits&, core::ResourcePlan&) const noexcept;

private:
    struct Data;
    explicit CaseMaterialChart(core::PlannedBuffer<Data>) noexcept;
    core::PlannedBuffer<Data> data_;
};

} // namespace uvlas::coupling::detail
