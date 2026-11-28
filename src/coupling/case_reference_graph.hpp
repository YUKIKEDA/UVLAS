#pragma once

#include <array>
#include <uvlas/io/si_case.hpp>
#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::structure::detail {
struct ReferenceArithmetic;
}

namespace uvlas::coupling::detail {

struct CaseReferenceSurfaceSample {
    std::array<structure::ExactRational, 3> position;
    std::array<structure::ExactRational, 3> chord_tangent;
    std::array<structure::ExactRational, 3> span_tangent;
    structure::PgReferenceCallInfo call_info;
};

struct CaseReferenceSurfaceIdentity {
    io::SchemaIdentity source;
    std::uint64_t geometry_id;
};

/// Owned exact reference for the Case's independent, thin surface graph.
class CaseReferenceSurface {
public:
    CaseReferenceSurface(CaseReferenceSurface&&) noexcept;
    CaseReferenceSurface& operator=(CaseReferenceSurface&&) noexcept;
    CaseReferenceSurface(const CaseReferenceSurface&) = delete;
    CaseReferenceSurface& operator=(const CaseReferenceSurface&) = delete;
    ~CaseReferenceSurface();

    /// Select a Geometry ID from this SI case, then check the frame and chart domain.
    [[nodiscard]] static core::Result<CaseReferenceSurface>
    prepare(const io::SiCase&, std::uint64_t geometry_id, const structure::PgReferenceLimits&,
            core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<CaseReferenceSurfaceIdentity> identity() const noexcept;
    [[nodiscard]] structure::PgReferenceCallInfo preparation_work() const noexcept;

    /// Return Cartesian metres and analytic derivatives with respect to chord and span.
    /// Points on the chart polygon boundary use the same polynomial extension as its interior.
    [[nodiscard]] core::Result<CaseReferenceSurfaceSample>
    evaluate(const structure::ExactRational& chord, const structure::ExactRational& span,
             const structure::PgReferenceLimits&, core::ResourcePlan&) const noexcept;

    // For callers that already qualified a chart subdivision and own the arithmetic budget.
    [[nodiscard]] CaseReferenceSurfaceSample evaluate_chart(const structure::PgPartitionPoint&,
                                                            structure::detail::ReferenceArithmetic&,
                                                            bool tangents) const;
    [[nodiscard]] structure::PgPartitionPoint
    chart_projection(const std::array<structure::ExactRational, 3>&,
                     structure::detail::ReferenceArithmetic&) const;

private:
    struct Data;
    explicit CaseReferenceSurface(core::PlannedBuffer<Data>) noexcept;
    core::PlannedBuffer<Data> data_;
};

} // namespace uvlas::coupling::detail
