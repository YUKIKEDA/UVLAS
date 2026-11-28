#pragma once

#include "aero/geometry_snapshot.hpp"
#include "case_reference_graph.hpp"
#include "structure/pg_reference_arithmetic.hpp"

namespace uvlas::coupling::detail {
/// Original Case field and source indexes for a failed physical preparation.
enum class CaseWingPhase : std::uint8_t {
    Case,
    Clock,
    Reference,
    Boundaries,
    Chord,
    Span,
    Pairing,
    Coverage,
    Vortex,
    Collocation,
    PrescribedMotion,
    Sharing,
    CoreRadius,
    Geometry,
    TrailingEdge
};

struct CaseWingDiagnostic {
    CaseWingPhase phase = CaseWingPhase::Case;
    std::optional<std::uint64_t> surface, other_surface, boundary, other_boundary;
    std::optional<std::uint64_t> strip, other_strip, chord, span;
    std::optional<aero::detail::GeometryProblem> geometry_entity;
    std::optional<aero::detail::FacetId> other_facet;
};

namespace wing_chart {
using Rational = structure::ExactRational;
using Math = structure::detail::ReferenceArithmetic;
using Point = structure::PgPartitionPoint;

struct Knot {
    Rational eta;
    Point leading, trailing;
};

struct Site {
    Point chart;
    Point chord_tangent, span_tangent;
    std::size_t strip{};
};

struct BoundaryRange {
    std::size_t offset{}, count{};
};

struct Chart {
    CaseReferenceSurface graph;
    std::uint64_t surface_id{}, leading_index{}, trailing_index{};
    bool reversed{};
    core::PlannedBuffer<Rational> chord, span;
    core::PlannedBuffer<Knot> knots;
    std::size_t knot_count{};
    core::PlannedBuffer<BoundaryRange> boundaries;
    core::PlannedBuffer<Point> boundary_points;
    [[nodiscard]] std::span<const Point> boundary(std::size_t) const;
    [[nodiscard]] Site site(const Rational&, const Rational&, bool derivatives, Math&,
                            CaseWingDiagnostic&) const;
};

[[nodiscard]] Chart prepare_chart(const io::SiCase&, const io::schema::AeroSurface&,
                                  const io::schema::GeometrySurfaceGraph&, Math&,
                                  CaseWingDiagnostic&);
[[nodiscard]] Point blend(const Point&, const Point&, const Rational&, Math&);
[[nodiscard]] Rational blend(const Rational&, const Rational&, const Rational&, Math&);
[[nodiscard]] bool on_segment(const Point&, const Point&, const Point&, Math&);
[[nodiscard]] std::optional<Rational> boundary_parameter(std::span<const Point>, const Point&,
                                                         Math&);
[[nodiscard]] Rational integer(std::size_t, Math&);
} // namespace wing_chart
} // namespace uvlas::coupling::detail
