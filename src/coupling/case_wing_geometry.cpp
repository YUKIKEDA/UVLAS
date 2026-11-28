#include "case_wing_sharing_private.hpp"
#include "discrete_nonpenetration.hpp"
#include "io/si_rational_projection.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <memory_resource>
#include <ranges>
#include <type_traits>

namespace uvlas::coupling::detail {
namespace wing_geometry {
namespace a = aero::detail;
namespace p = a::pchip;
namespace pg = structure::detail;
namespace s = io::schema;
namespace si = io::detail;
using namespace wing_chart;
using enum core::FailureCode;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_mul;
using pg::partition_take;
std::size_t size(std::uint64_t count) {
    return static_cast<std::size_t>(partition_index(count));
}

constexpr auto operation = core::Operation::PrepareModel;
using XYZ = std::array<Rational, 3>;

struct StraightWake {
    XYZ displacement, velocity;
    io::SiEnclosure duration;
};

struct Owner {
    std::optional<io::SchemaIdentity> source;
    std::optional<a::GeometryContext> context;
    std::optional<CaseFlowClockQuery> address;
    std::optional<StraightWake> wake;
    core::PlannedBuffer<CaseWingWakeBand> wake_bands;
    core::PlannedBuffer<std::size_t> terminal_anchor;
    si::SiRationalWords instant;
    core::PlannedBuffer<Patch> patches;
    core::PlannedBuffer<CaseWingPoint> points;
    core::PlannedBuffer<CaseWingCell> cells;
    core::PlannedBuffer<CaseWingTrailingBand> bands;
    core::PlannedBuffer<CaseWingSharing> sharing;
    std::size_t band_count{};
    std::optional<a::GlobalContactState> global;
    std::optional<a::AliasKinematicsState> kinematics;
};

template <class T>
core::PlannedBuffer<T> allocate(Math& m, std::uint64_t count) {
    return partition_take(m.plan.allocate<T>(partition_index(count), operation));
}

template <class T, class E>
T take_aero(std::expected<T, E> result, CaseWingDiagnostic& diagnostic) {
    if (!result) {
        if constexpr (requires { result.error().cause; }) {
            diagnostic.geometry_entity = result.error().problem;
            throw pg::PartitionAbort{result.error().cause};
        } else {
            diagnostic.geometry_entity = result.error().origin.problem;
            diagnostic.other_facet = result.error().other_facet;
            throw pg::PartitionAbort{result.error().origin.cause};
        }
    }
    return std::move(*result);
}

template <class T>
class Inputs {
public:
    Inputs(Math& m, std::size_t count, const T& value)
        : payload_(allocate<std::byte>(
              m, partition_add(partition_mul(count, sizeof(T)), alignof(T) - 1))),
          arena_(payload_.values().data(), payload_.values().size_bytes(),
                 std::pmr::null_memory_resource()) {
        static_assert(std::is_trivially_destructible_v<T>);
        if (count == 0)
            return;
        // The planned arena covers one array. Debug vector iterators may need another allocation.
        auto* storage = std::pmr::polymorphic_allocator<T>{&arena_}.allocate(count);
        std::uninitialized_fill_n(storage, count, value);
        values_ = std::span<T>{storage, count};
    }

    [[nodiscard]] std::span<T> values() noexcept {
        return values_;
    }

private:
    core::PlannedBuffer<std::byte> payload_;
    std::pmr::monotonic_buffer_resource arena_;
    std::span<T> values_;
};

bool same_source(io::SchemaIdentity a, io::SchemaIdentity b) noexcept {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}

bool same_address(const std::optional<CaseFlowClockQuery>& a,
                  const std::optional<CaseFlowClockQuery>& b) noexcept {
    return a.has_value() == b.has_value() &&
           (!a || (a->stamp == b->stamp && a->interval == b->interval && a->step == b->step &&
                   a->side == b->side));
}

si::SiRationalWordsView si_view(p::RationalView words) {
    return {words.sign, words.numerator, words.denominator};
}

p::RationalView curve_view(const si::SiRationalWords& words) {
    const auto view = words.view();
    return {view.sign, view.numerator, view.denominator};
}

struct Projector {
    const CaseWingProjection& request;
    std::uint64_t remaining;

    io::SiEnclosure apply(const Rational& value, const CaseFlowComponentBudget& budget,
                          core::Precision representation, Math& m) {
        const auto owner =
            partition_take(si::copy_rational_words(value, m.plan, m.limits.cancellation));
        const auto projected = partition_take(si::project_rational_words(
            owner.view(),
            {si_view(budget.scale), si_view(budget.absolute), si_view(budget.relative)},
            representation, {request.exact.input_integer_bits, remaining}, m.limits.cancellation,
            {}));
        remaining -= projected.used_steps;
        return projected.value;
    }

    std::array<io::SiEnclosure, 3>
    vector(const XYZ& value, const std::array<CaseFlowComponentBudget, 3>& budgets, Math& m) {
        std::array<io::SiEnclosure, 3> result;
        for (std::size_t axis = 0; axis < result.size(); ++axis)
            result[axis] = apply(value[axis], budgets[axis], request.representation, m);
        return result;
    }
};

template <class Record, class Variant>
Record select(const Variant& value) {
    if (const auto* record = std::get_if<Record>(&value))
        return *record;
    partition_fail(Unsupported);
}

template <class Field>
Rational number(const io::SiCase& source, Field field, Math& m) {
    m.charge();
    const auto recipe = partition_take(source.quantity(field)).recipe;
    if (recipe.pi_power != 0)
        partition_fail(InvalidInput);
    return partition_take(io::materialize_rational(recipe, m.leaf_limits, m.plan));
}

XYZ constant_velocity(const io::SiCase& source, const s::TimeSeriesConstant& series, Math& m) {
    return {number(source, partition_take(series.value.at(0)), m),
            number(source, partition_take(series.value.at(1)), m),
            number(source, partition_take(series.value.at(2)), m)};
}

StraightWake prepare_straight_wake(const io::SiCase& source, const s::Case& root,
                                   const s::WakeActive& wake, const CaseWingStraightWake& request,
                                   Projector& project, Math& m) {
    const auto convection = select<s::ConvectionPrescribed>(wake.convection.record());
    const auto retention = select<s::RetentionWindow>(wake.retention.record());
    const auto flow = select<s::FlowKnown>(root.flow.record());
    (void)select<s::GustNone>(flow.gust.record());
    auto velocity =
        constant_velocity(source, select<s::TimeSeriesConstant>(convection.velocity.record()), m);
    if (velocity !=
        constant_velocity(source, select<s::TimeSeriesConstant>(flow.uniform.record()), m))
        partition_fail(Unsupported);
    const auto duration = number(source, retention.duration, m);
    const auto distance = number(source, retention.distance, m);
    if (duration.sign() <= 0 || distance.sign() <= 0)
        partition_fail(InvalidInput);
    XYZ advance;
    Rational squared;
    for (std::size_t axis = 0; axis < advance.size(); ++axis) {
        advance[axis] = m.mul(duration, velocity[axis]);
        if (request.expected_displacement_metres &&
            (*request.expected_displacement_metres)[axis] != advance[axis])
            partition_fail(StaleGeneration);
        squared = m.add(squared, m.square(advance[axis]));
    }
    if (squared != m.square(distance))
        partition_fail(InvalidInput);
    auto duration_projection =
        project.apply(duration, project.request.time, core::Precision::TwoWord, m);
    return {std::move(advance), std::move(velocity), std::move(duration_projection)};
}

double unscale(double word, std::int64_t exponent) {
    if (word == 0)
        return word;
    if (exponent < -2200 || exponent > 2200)
        partition_fail(PrecisionUnresolved);
    const auto shift = static_cast<int>(exponent);
    const auto value = std::scalbn(word, shift);
    if (!std::isfinite(value) || std::scalbn(value, -shift) != word)
        partition_fail(PrecisionUnresolved);
    return value;
}

a::ScalarWords retained(io::SiEnclosure value) {
    if (const auto* single = std::get_if<io::ScaledBinary64>(&value))
        return a::SingleWord{unscale(single->center, single->exponent2)};
    const auto pair = std::get<io::ScaledTwoWord>(value);
    return a::DoubleWord{unscale(pair.high, pair.exponent2), unscale(pair.low, pair.exponent2)};
}

XYZ tangent(const CaseReferenceSurfaceSample& graph, const Point& chart, Math& m) {
    XYZ result;
    for (std::size_t axis = 0; axis < result.size(); ++axis)
        result[axis] = m.add(m.mul(graph.chord_tangent[axis], chart.s),
                             m.mul(graph.span_tangent[axis], chart.t));
    return result;
}

void reverse_tangent(NonpenetrationVector& tangent) {
    for (auto& component : tangent) {
        if (auto* single = std::get_if<a::SingleWord>(&component)) {
            single->value = -single->value;
        } else {
            auto& pair = std::get<a::DoubleWord>(component);
            pair.high = -pair.high;
            pair.low = -pair.low;
        }
    }
}

s::GeometrySurfaceGraph find_graph(const s::Model& model, std::uint64_t id, Math& m) {
    for (std::size_t i = 0; i < model.geometry.size(); ++i) {
        m.visit(1);
        const auto geometry = partition_take(model.geometry.at(partition_index(i))).record();
        if (std::visit([](const auto& value) { return value.id.unsigned_value(); }, geometry) != id)
            continue;
        if (const auto* graph = std::get_if<s::GeometrySurfaceGraph>(&geometry))
            return *graph;
        partition_fail(Unsupported);
    }
    partition_fail(InvalidInput);
}

struct Stage {
    Rational instant;
    p::PchipInputStamp stamp;
    p::RequestedSide side;
    std::optional<case_clock_detail::Grid> grid;

    [[nodiscard]] CaseMotionTimeWindow window() const {
        if (!grid)
            return {instant, instant};
        const auto intervals = grid->intervals.values();
        return {intervals.front().begin, intervals.back().end};
    }
};

Stage stage(const io::SiCase& source, const a::GeometryContext& context,
            std::optional<CaseFlowClockQuery> address, bool transient, const CaseWingLimits& limits,
            Math& m) {
    if (transient != address.has_value())
        partition_fail(StaleGeneration);
    if (!address)
        return {{},
                {p::PchipIssuerId{context.model.value()},
                 p::PchipInputGeneration{context.current.value()}},
                p::RequestedSide::Automatic,
                {}};
    auto grid = partition_take(
        prepare_case_reference_clock(source, address->stamp, limits.clock, limits.clock_event_work,
                                     m.plan, m.limits.cancellation, context.diagnostics));
    auto sample = partition_take(evaluate_case_reference_clock(
        grid, grid, *address, limits.clock, m.plan, m.limits.cancellation, context.diagnostics));
    return {std::move(sample.seconds), address->stamp, address->side, std::move(grid)};
}

struct MotionPoints {
    using Words = std::array<si::SiRationalWords, 3>;
    std::array<CaseFlowComponentBudget, 6> positions, chord, span;
    core::PlannedBuffer<Words> words;
    core::PlannedBuffer<CaseMotionReferencePoint> requests;

    MotionPoints(const CaseWingProjection& projection, std::size_t count, Math& m)
        : words(allocate<Words>(m, count)), requests(allocate<CaseMotionReferencePoint>(m, count)) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            positions[axis] = projection.position[axis];
            chord[axis] = projection.chord_tangent[axis];
            span[axis] = projection.span_tangent[axis];
            positions[axis + 3] = chord[axis + 3] = span[axis + 3] = projection.velocity[axis];
        }
    }

    void set(std::size_t index, const XYZ& point, std::span<const CaseFlowComponentBudget> budgets,
             CaseMotionReferenceKind kind, Math& m) {
        auto& values = words.values()[index];
        auto& request = requests.values()[index];
        request.components = budgets;
        request.kind = kind;
        for (std::size_t axis = 0; axis < values.size(); ++axis) {
            values[axis] =
                partition_take(si::copy_rational_words(point[axis], m.plan, m.limits.cancellation));
            request.world_metres[axis] = curve_view(values[axis]);
        }
    }
};

std::array<io::SiEnclosure, 3> motion_value(const CaseMotionEvaluation& motion, std::size_t index,
                                            p::JetOrder order) {
    return {partition_take(motion.point(index, 0, order)),
            partition_take(motion.point(index, 1, order)),
            partition_take(motion.point(index, 2, order))};
}

void sample_patch(const io::SiCase& source, const s::AeroSurface& input, Patch& patch,
                  const Stage& time, const CaseWingLimits& limits, Projector& projector,
                  p::CurveWorkspace& curve, Owner& owner, Math& m, CaseWingDiagnostic& diagnostic) {
    const auto& chart = *patch.chart;
    diagnostic = {};
    diagnostic.phase = CaseWingPhase::PrescribedMotion;
    diagnostic.surface = chart.surface_id;
    const auto attachment = input.attachment.record();
    const auto* prescribed = std::get_if<s::AeroAttachmentPrescribed>(&attachment);
    if (!prescribed && !std::holds_alternative<s::AeroAttachmentFixed>(attachment))
        partition_fail(Unsupported);
    const auto point_count = (patch.m() + 1) * (patch.n() + 1);
    const auto cell_count = patch.m() * patch.n();
    const auto query_count = partition_add(point_count, partition_mul(3, cell_count));
    MotionPoints batch{projector.request, prescribed ? size(query_count) : 0U, m};
    const auto zero = Rational{};
    const auto one = m.constant("1");
    const auto four = m.constant("4");
    const auto quarter = m.div(one, four);
    const auto three_quarters = m.div(m.constant("3"), four);
    const auto half = m.div(one, m.constant("2"));
    const auto zero_velocity = projector.vector(XYZ{}, projector.request.velocity, m);
    diagnostic.phase = CaseWingPhase::Vortex;
    for (std::size_t i = 0; i <= patch.m(); ++i) {
        const auto xi = i == patch.m() ? m.copy(one)
                                       : blend(chart.chord.values()[i], chart.chord.values()[i + 1],
                                               quarter, m);
        for (std::size_t j = 0; j <= patch.n(); ++j) {
            m.visit(1);
            diagnostic.chord = i;
            diagnostic.span = j;
            auto site = chart.site(xi, chart.span.values()[j], false, m, diagnostic);
            const auto graph = chart.graph.evaluate_chart(site.chart, m, false);
            const auto index = patch.point(i, j);
            auto& point = owner.points.values()[index];
            point = {a::PointId{index + 1},
                     a::PointId{index + 1},
                     chart.surface_id,
                     i,
                     j,
                     std::min(i, patch.m() - 1),
                     j == 0 ? 0 : j - 1,
                     std::min(j, patch.n() - 1),
                     site.strip,
                     std::move(site.chart),
                     projector.vector(graph.position, projector.request.position, m),
                     {},
                     zero_velocity};
            point.position = point.reference;
            if (prescribed)
                batch.set(index - patch.points, graph.position, batch.positions,
                          CaseMotionReferenceKind::Position, m);
        }
    }
    diagnostic.phase = CaseWingPhase::Collocation;
    for (std::size_t i = 0; i < patch.m(); ++i) {
        const auto xi =
            blend(chart.chord.values()[i], chart.chord.values()[i + 1], three_quarters, m);
        for (std::size_t j = 0; j < patch.n(); ++j) {
            m.visit(1);
            const auto eta = blend(chart.span.values()[j], chart.span.values()[j + 1], half, m);
            diagnostic.chord = i;
            diagnostic.span = j;
            auto site = chart.site(xi, eta, true, m, diagnostic);
            const auto graph = chart.graph.evaluate_chart(site.chart, m, true);
            const auto chord = tangent(graph, site.chord_tangent, m);
            const auto span = tangent(graph, site.span_tangent, m);
            const auto index = patch.cells + i * patch.n() + j;
            auto& cell = owner.cells.values()[index];
            cell = {chart.surface_id,
                    i,
                    j,
                    site.strip,
                    chart.reversed,
                    a::PhysicalFaceId{index + 1},
                    {a::FacetId{2 * index + 1}, a::FacetId{2 * index + 2}},
                    std::move(site.chart),
                    {},
                    zero_velocity,
                    {},
                    {},
                    {},
                    {}};
            if (prescribed) {
                const auto start = point_count + 3 * (index - patch.cells);
                batch.set(start, graph.position, batch.positions, CaseMotionReferenceKind::Position,
                          m);
                batch.set(start + 1, chord, batch.chord, CaseMotionReferenceKind::Tangent, m);
                batch.set(start + 2, span, batch.span, CaseMotionReferenceKind::Tangent, m);
            } else {
                cell.position = projector.vector(graph.position, projector.request.position, m);
                cell.chord_tangent = projector.vector(chord, projector.request.chord_tangent, m);
                cell.span_tangent = projector.vector(span, projector.request.span_tangent, m);
            }
        }
    }
    if (!prescribed)
        return;
    diagnostic = {};
    diagnostic.phase = CaseWingPhase::PrescribedMotion;
    diagnostic.surface = chart.surface_id;
    const auto ready =
        partition_take(prepare_case_motion(source, prescribed->motion, time.window(), time.stamp,
                                           curve, limits.motion, m.plan, m.limits.cancellation));
    const auto motion = partition_take(evaluate_case_motion(
        ready, ready, {time.stamp, time.instant, p::JetOrder::First, time.side},
        {projector.request.representation, projector.request.motion_pose, batch.requests.values()},
        curve, limits.motion, m.plan, m.limits.cancellation));
    patch.motion = partition_take(
        snapshot_case_motion(ready, curve, limits.motion, m.plan, m.limits.cancellation));
    for (std::size_t i = 0; i < point_count; ++i) {
        m.visit(1);
        auto& point = owner.points.values()[patch.points + i];
        point.position = motion_value(motion, i, p::JetOrder::Value);
        point.velocity = motion_value(motion, i, p::JetOrder::First);
    }
    for (std::size_t i = 0; i < cell_count; ++i) {
        m.visit(1);
        auto& cell = owner.cells.values()[patch.cells + i];
        const auto index = point_count + 3 * i;
        cell.position = motion_value(motion, index, p::JetOrder::Value);
        cell.velocity = motion_value(motion, index, p::JetOrder::First);
        cell.chord_tangent = motion_value(motion, index + 1, p::JetOrder::Value);
        cell.span_tangent = motion_value(motion, index + 2, p::JetOrder::Value);
    }
}

a::PointId point_id(std::size_t index) {
    return a::PointId{index + 1};
}
a::EdgeId edge_id(std::size_t index) {
    return a::EdgeId{index + 1};
}

struct MeshCounts {
    std::uint64_t points{}, cells{}, edges{}, bands{};

    void append(Patch& patch) {
        patch.points = size(points);
        patch.cells = size(cells);
        patch.edges = size(edges);
        patch.bands = size(bands);
        const auto rows = partition_add(patch.m(), 1), cols = partition_add(patch.n(), 1);
        points = partition_add(points, partition_mul(rows, cols));
        const auto faces = partition_mul(patch.m(), patch.n());
        cells = partition_add(cells, faces);
        edges = partition_add(edges,
                              partition_add(faces, partition_add(partition_mul(rows, patch.n()),
                                                                 partition_mul(patch.m(), cols))));
        bands = partition_add(bands, patch.n());
        (void)partition_index(partition_mul(2, cells));
    }

    void admit(const a::GeometryCounts& bounds) const {
        if (points > bounds.points || cells > bounds.faces || cells > bounds.quads ||
            bands > bounds.bands || edges > bounds.edges)
            partition_fail(ResourceExceeded);
        (void)partition_index(partition_mul(3, points));
        (void)partition_index(edges);
    }
};

struct WakeCounts {
    std::uint64_t vertices{}, bands{};
};

WakeCounts count_wake(const Owner& owner, std::span<std::size_t> point_owners,
                      std::span<const bool> internal, std::span<std::size_t> terminal, Math& m) {
    WakeCounts count;
    for (const auto& patch : owner.patches.values())
        for (std::size_t j = 0; j < patch.n(); ++j) {
            m.visit(1);
            if (internal[patch.span_edge(patch.m(), j)])
                continue;
            ++count.bands;
            for (const auto end : {j, j + 1}) {
                auto& index = terminal[root(point_owners, patch.point(patch.m(), end), m)];
                if (index == 0)
                    index = size(++count.vertices);
            }
        }
    if (count.bands == 0)
        partition_fail(GeometryInvalid);
    return count;
}

class Mesh {
public:
    Mesh(const MeshCounts& count, core::Precision format, Math& m)
        : points_(m, size(count.points), a::PointId{0}),
          bindings_(m, size(count.points),
                    {a::ReferenceMapId{0}, a::MappingGeneration{0}, a::CanonicalPoint{}}),
          edges_(m, size(count.edges), {a::EdgeId{0}, {a::PointId{0}, a::PointId{0}}}),
          faces_(m, size(count.cells), {a::PhysicalFaceId{0}, a::BandId{0}}),
          bands_(m, size(count.bands), a::BandId{0}),
          quads_(m, size(count.cells),
                 {a::PhysicalFaceId{0},
                  {a::PointId{0}, a::PointId{0}, a::PointId{0}, a::PointId{0}},
                  {a::EdgeId{0}, a::EdgeId{0}, a::EdgeId{0}, a::EdgeId{0}},
                  a::EdgeId{0},
                  {a::FacetId{0}, a::FacetId{0}}}),
          roles_(m, size(count.cells), {a::PhysicalFaceId{0}, a::AeroFaceRole::Bound}),
          owners_(m, size(count.edges), {a::EdgeId{0}, a::CanonicalEdge{}}),
          trailing_(m, size(count.bands), {a::EdgeId{0}, {a::PointId{0}, a::PointId{0}}}),
          ref_high_(allocate<double>(m, partition_mul(count.points, 3))),
          ref_low_(allocate<double>(m, format == core::Precision::TwoWord ? ref_high_.size() : 0)),
          cur_high_(allocate<double>(m, ref_high_.size())),
          cur_low_(allocate<double>(m, ref_low_.size())) {}

    void place_points(Owner& owner, std::span<std::size_t> point_owners, Math& m) {
        const auto stamp = owner.context->mapping;
        auto points = owner.points.values();
        const auto place = [&](io::SiEnclosure value, std::size_t index, auto high, auto low) {
            const auto words = retained(value);
            if (const auto* two = std::get_if<a::DoubleWord>(&words)) {
                high[index] = two->high;
                low[index] = two->low;
            } else {
                high[index] = std::get<a::SingleWord>(words).value;
            }
        };
        for (std::size_t i = 0; i < points.size(); ++i) {
            m.visit(1);
            const auto id = point_id(i);
            const auto canonical = root(point_owners, i, m);
            points_.values()[i] = id;
            auto& point = points[i];
            const auto& chosen = points[canonical];
            point.canonical = chosen.original;
            bindings_.values()[i] = {a::ReferenceMapId{point.canonical.value()}, stamp,
                                     i == canonical
                                         ? a::PointOwnership{a::CanonicalPoint{}}
                                         : a::PointOwnership{a::AliasOf{point.canonical}}};
            // The exact reference curve and motion have already proved equality.
            // Reuse the one finite representation for both topological aliases.
            point.reference = chosen.reference;
            point.position = chosen.position;
            point.velocity = chosen.velocity;
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const auto coordinate = 3 * i + axis;
                place(point.reference[axis], coordinate, ref_high_.values(), ref_low_.values());
                place(point.position[axis], coordinate, cur_high_.values(), cur_low_.values());
            }
        }
    }

    void connect(const Patch& patch, Owner& owner, std::span<std::size_t> edge_owners,
                 std::span<const bool> internal, Math& m) {
        const auto point = [&](std::size_t i, std::size_t j) {
            return point_id(patch.point(i, j));
        };
        const auto edge = [&](std::size_t n, a::PointId first, a::PointId second) {
            m.visit(1);
            const auto id = edge_id(n);
            const auto primary = root(edge_owners, n, m);
            edges_.values()[n] = {id, {first, second}};
            owners_.values()[n] = {id, n == primary ? a::EdgeOwnership{a::CanonicalEdge{}}
                                                    : a::EdgeOwnership{edge_id(primary)}};
        };
        for (std::size_t i = 0; i < patch.m(); ++i)
            for (std::size_t j = 0; j <= patch.n(); ++j)
                edge(patch.chord_edge(i, j), point(i, j), point(i + 1, j));
        for (std::size_t i = 0; i <= patch.m(); ++i)
            for (std::size_t j = 0; j < patch.n(); ++j)
                edge(patch.span_edge(i, j), point(i, j), point(i, j + 1));
        for (std::size_t j = 0; j < patch.n(); ++j)
            bands_.values()[patch.bands + j] = a::BandId{patch.bands + j + 1};
        const bool reverse = patch.chart->reversed;
        for (std::size_t i = 0; i < patch.m(); ++i)
            for (std::size_t j = 0; j < patch.n(); ++j) {
                m.visit(1);
                const auto diagonal = patch.diagonal(i, j);
                std::array vertex{point(i, j), point(i, j + 1), point(i + 1, j + 1),
                                  point(i + 1, j)};
                std::array outline{
                    edge_id(patch.span_edge(i, j)), edge_id(patch.chord_edge(i, j + 1)),
                    edge_id(patch.span_edge(i + 1, j)), edge_id(patch.chord_edge(i, j))};
                edge(diagonal, vertex[0], vertex[2]);
                const auto f = patch.cells + i * patch.n() + j;
                const auto& cell = owner.cells.values()[f];
                faces_.values()[f] = {cell.face, bands_.values()[patch.bands + j]};
                roles_.values()[f] = {cell.face, a::AeroFaceRole::Bound};
                if (reverse) {
                    std::swap(vertex[1], vertex[3]);
                    std::ranges::reverse(outline);
                }
                quads_.values()[f] = {cell.face, vertex, outline, edge_id(diagonal), cell.facets};
            }
        for (std::size_t j = 0; j < patch.n(); ++j) {
            m.visit(1);
            const auto trailing = patch.span_edge(patch.m(), j);
            if (internal[trailing])
                continue;
            const auto& low = owner.points.values()[patch.point(patch.m(), j)];
            const auto& high = owner.points.values()[patch.point(patch.m(), j + 1)];
            const auto end = std::array{low.original, high.original};
            trailing_.values()[owner.band_count] = {edge_id(trailing),
                                                    reverse ? std::array{end[1], end[0]} : end};
            const auto& face = owner.cells.values()[patch.cells + (patch.m() - 1) * patch.n() + j];
            const auto boundary = patch.chart->boundary(size(patch.chart->trailing_index));
            auto begin = boundary_parameter(boundary, low.chart, m);
            auto end_param = boundary_parameter(boundary, high.chart, m);
            if (!begin || !end_param)
                partition_fail(GeometryInvalid);
            owner.bands.values()[owner.band_count++] = {patch.chart->surface_id,
                                                        patch.chart->trailing_index,
                                                        j,
                                                        {low.strip, high.strip},
                                                        {edge_id(trailing), face.face, 0, 0},
                                                        {std::move(*begin), std::move(*end_param)}};
        }
    }

    void extend(const MeshCounts& bound, std::span<const std::size_t> terminal, Owner& owner,
                Projector& projection, Math& m) {
        const auto& delta = owner.wake->displacement;
        const auto terminal_id = [&](std::size_t anchor) {
            return point_id(size(bound.points) + terminal[anchor] - 1);
        };
        const auto stream_edge = [&](std::size_t anchor) {
            return edge_id(size(bound.edges) + terminal[anchor] - 1);
        };
        const auto edge = [&](a::EdgeId id, a::PointId first, a::PointId last) {
            edges_.values()[size(id.value() - 1)] = {id, {first, last}};
            owners_.values()[size(id.value() - 1)] = {id, a::CanonicalEdge{}};
        };
        for (std::size_t i = 0; i < terminal.size(); ++i) {
            m.visit(1);
            if (terminal[i] == 0)
                continue;
            const auto id = terminal_id(i);
            const auto index = size(id.value() - 1);
            points_.values()[index] = id;
            bindings_.values()[index] = {a::ReferenceMapId{id.value()}, owner.context->mapping,
                                         a::CanonicalPoint{}};
            owner.terminal_anchor.values()[terminal[i] - 1] = i;
            const auto& source = owner.points.values()[i];
            const auto patch = std::ranges::find_if(owner.patches.values(), [&](const Patch& p) {
                m.visit(1);
                return p.chart->surface_id == source.surface;
            });
            if (patch == owner.patches.values().end())
                partition_fail(StaleGeneration);
            const auto front = patch->chart->graph.evaluate_chart(source.chart, m, false);
            for (std::size_t axis = 0; axis < delta.size(); ++axis) {
                // A transverse coordinate does not move, so its already qualified
                // front representation remains valid at the end of the wake.
                const auto projected = [&]() -> a::ScalarWords {
                    if (delta[axis].sign() == 0)
                        return retained(source.position[axis]);
                    const auto value = m.add(front.position[axis], delta[axis]);
                    return retained(projection.apply(value, projection.request.position[axis],
                                                     core::Precision::TwoWord, m));
                }();
                const auto pair = std::holds_alternative<a::DoubleWord>(projected)
                                      ? std::get<a::DoubleWord>(projected)
                                      : a::DoubleWord{std::get<a::SingleWord>(projected).value, 0};
                const auto slot = 3 * index + axis;
                ref_high_.values()[slot] = cur_high_.values()[slot] = pair.high;
                ref_low_.values()[slot] = cur_low_.values()[slot] = pair.low;
            }
            edge(stream_edge(i), point_id(i), id);
        }
        for (std::size_t i = 0; i < owner.band_count; ++i) {
            m.visit(1);
            const auto front = trailing_.values()[i];
            const auto& band = owner.bands.values()[i];
            const auto canonical = [&](a::PointId id) {
                return size(owner.points.values()[size(id.value() - 1)].canonical.value() - 1);
            };
            const auto p0 = canonical(front.points[0]), p1 = canonical(front.points[1]);
            const auto end0 = terminal_id(p0), end1 = terminal_id(p1);
            const auto tip_edge = edge_id(size(bound.edges) + owner.terminal_anchor.size() + 2 * i);
            const auto diagonal = edge_id(size(tip_edge.value()));
            edge(tip_edge, end0, end1);
            edge(diagonal, front.points[0], end1);
            const auto row = size(bound.cells) + i;
            const auto face = a::PhysicalFaceId{row + 1};
            const auto facet = a::FacetId{2 * row + 1};
            const auto bound_band =
                faces_.values()[size(band.transfer.bound_face.value() - 1)].band;
            faces_.values()[row] = {face, bound_band};
            roles_.values()[row] = {face, a::AeroFaceRole::Wake};
            quads_.values()[row] = {face,
                                    {front.points[0], front.points[1], end1, end0},
                                    {front.edge, stream_edge(p1), tip_edge, stream_edge(p0)},
                                    diagonal,
                                    {facet, a::FacetId{2 * row + 2}}};
            owner.wake_bands.values()[i] = {
                front.edge, band.transfer.bound_face, face, facet, front.points, {end0, end1}};
        }
    }

    void publish(Owner& owner, a::ScalarWords core_radius, const CaseWingLimits& limits,
                 Projector& projection, a::GeometryWorkspace& geometry_work,
                 a::GlobalWorkspace& global_work, Math& m, CaseWingDiagnostic& diagnostic,
                 const a::GlobalTopology* topology) {
        const auto position = [](std::span<const double> high,
                                 std::span<const double> low) -> a::PositionView {
            if (low.empty())
                return a::SinglePositions{high};
            return a::DoublePositions{high, low};
        };
        const a::GeometryInput input{points_.values(),
                                     bindings_.values(),
                                     position(ref_high_.values(), ref_low_.values()),
                                     position(cur_high_.values(), cur_low_.values()),
                                     edges_.values(),
                                     faces_.values(),
                                     bands_.values(),
                                     {},
                                     quads_.values(),
                                     core_radius,
                                     *owner.context,
                                     curve_view(owner.instant)};
        const auto stop = m.limits.cancellation;
        diagnostic = {};
        diagnostic.phase = CaseWingPhase::Geometry;
        const auto geometry = take_aero(
            a::prepare_owned_geometry(input, limits.geometry, geometry_work, m.plan, stop),
            diagnostic);
        std::optional<a::GlobalTopology> built;
        if (topology == nullptr) {
            built.emplace(
                take_aero(a::build_global_topology(geometry,
                                                   {roles_.values(), owners_.values(),
                                                    trailing_.values().first(owner.band_count)},
                                                   limits.global, global_work, m.plan, stop),
                          diagnostic));
            topology = &*built;
        }
        auto global = take_aero(
            a::check_global_geometry(*topology, geometry, limits.global, global_work, m.plan, stop),
            diagnostic);
        const a::SingleWord zero{0.0};
        Inputs<a::PointKinematicsInput> velocities{
            m, points_.values().size(), {a::PointId{0}, {}, {zero, zero, zero}}};
        for (std::size_t i = 0; i < velocities.values().size(); ++i) {
            m.visit(1);
            const auto anchor = i < owner.points.size()
                                    ? i
                                    : owner.terminal_anchor.values()[i - owner.points.size()];
            const auto& point = owner.points.values()[anchor];
            velocities.values()[i] = {point_id(i),
                                      {},
                                      {retained(point.velocity[0]), retained(point.velocity[1]),
                                       retained(point.velocity[2])}};
        }
        owner.kinematics = take_aero(
            a::check_alias_kinematics(
                geometry,
                {a::DofLayoutGeneration{owner.context->mapping.value()}, {}, velocities.values()},
                limits.global, global_work, m.plan, stop),
            diagnostic);
        // Edges and band drafts are emitted in increasing original edge order.
        diagnostic.phase = CaseWingPhase::TrailingEdge;
        const auto rows = topology->trailing_edges();
        if (rows.size() != owner.band_count)
            partition_fail(GeometryInvalid);
        for (std::size_t i = 0; i < rows.size(); ++i) {
            m.visit(1);
            auto& band = owner.bands.values()[i];
            diagnostic.surface = band.surface;
            diagnostic.boundary = band.source_boundary;
            diagnostic.span = band.span;
            if (rows[i].edge != band.transfer.edge ||
                rows[i].bound_face != band.transfer.bound_face)
                partition_fail(GeometryInvalid);
            band.transfer = rows[i];
        }
        const auto words = [](const auto& values) -> NonpenetrationVector {
            return {retained(values[0]), retained(values[1]), retained(values[2])};
        };
        diagnostic = {};
        diagnostic.phase = CaseWingPhase::Collocation;
        for (auto& cell : owner.cells.values()) {
            m.visit(1);
            diagnostic.surface = cell.surface;
            diagnostic.chord = cell.chord;
            diagnostic.span = cell.span;
            diagnostic.strip = cell.strip;
            const auto normal = partition_take(project_case_wing_normal(
                geometry, cell.facets, cell.chord_tangent, cell.span_tangent, cell.reversed,
                projection.request.representation, projection.request.normal,
                {projection.request.exact.input_integer_bits, projection.remaining}, stop));
            projection.remaining -= normal.used_steps;
            cell.normal = normal.unit;
            const auto retained_position = words(cell.position);
            const auto free_point = [&]() -> a::PointWords {
                if (projection.request.representation == core::Precision::Binary64)
                    return a::FreeSingle{{std::get<a::SingleWord>(retained_position[0]).value,
                                          std::get<a::SingleWord>(retained_position[1]).value,
                                          std::get<a::SingleWord>(retained_position[2]).value}};
                a::FreeDouble point{};
                for (std::size_t axis = 0; axis < point.high.size(); ++axis) {
                    const auto value = std::get<a::DoubleWord>(retained_position[axis]);
                    point.high[axis] = value.high;
                    point.low[axis] = value.low;
                }
                return point;
            };
            auto query =
                take_aero(a::resolve_query(
                              geometry, {a::QueryId{cell.face.value()}, a::FreePoint{free_point()}},
                              *owner.context, limits.geometry, geometry_work, m.plan, stop),
                          diagnostic);
            m.visit(query.used_steps());
            cell.collocation = std::move(query);
        }
        owner.global = std::move(global);
    }

private:
    Inputs<a::PointId> points_;
    Inputs<a::PointBinding> bindings_;
    Inputs<a::EdgeInput> edges_;
    Inputs<a::FaceInput> faces_;
    Inputs<a::BandId> bands_;
    Inputs<a::QuadInput> quads_;
    Inputs<a::FaceRoleInput> roles_;
    Inputs<a::EdgeBindingInput> owners_;
    Inputs<a::TrailingEdgeInput> trailing_;
    core::PlannedBuffer<double> ref_high_, ref_low_, cur_high_, cur_low_;
};

struct Build {
    static CaseWingGeometry run(const io::SiCase& source, const a::GeometryContext& context,
                                std::optional<CaseFlowClockQuery> address,
                                const CaseWingProjection& projection, const CaseWingLimits& limits,
                                a::GeometryWorkspace& geometry_work,
                                a::GlobalWorkspace& global_work, p::CurveWorkspace& curve,
                                core::ResourcePlan& plan, CaseWingDiagnostic& diagnostic,
                                std::optional<CaseWingStraightWake> straight, const Owner* anchor) {
        Math m{limits.chart, plan};
        m.cancel();
        const auto& schema = partition_take(source.source_case()).get();
        const auto identity = partition_take(schema.identity());
        std::optional<a::GlobalTopology> topology;
        if (anchor != nullptr) {
            diagnostic.phase = CaseWingPhase::Clock;
            if (!anchor->source || !same_source(*anchor->source, identity) || !anchor->global)
                partition_fail(StaleGeneration);
            if (anchor->wake) {
                diagnostic.phase = CaseWingPhase::TrailingEdge;
                partition_fail(Unsupported);
            }
            if (!anchor->address || !address || anchor->address->stamp != address->stamp ||
                (anchor->context->current == context.current &&
                 !same_address(anchor->address, address)))
                partition_fail(StaleGeneration);
            topology.emplace(take_aero(anchor->global->retain_reference_topology(), diagnostic));
        }
        const auto root = partition_take(schema.root());
        const auto analysis = root.analysis.record();
        const bool transient = std::holds_alternative<s::AnalysisTransient>(analysis);
        if (!transient && !std::holds_alternative<s::AnalysisStatic>(analysis))
            partition_fail(Unsupported);
        const auto wake = std::visit(
            [](const auto& item) -> s::WakeActive {
                if constexpr (requires { item.wake; }) {
                    const auto value = item.wake.record();
                    if (const auto* active = std::get_if<s::WakeActive>(&value))
                        return *active;
                }
                partition_fail(Unsupported);
            },
            analysis);
        auto carrier = allocate<Owner>(m, 1);
        auto& result = carrier.values().front();
        result.source = identity;
        result.context = context;
        result.address = address;
        diagnostic.phase = CaseWingPhase::Clock;
        const auto time = stage(source, context, address, transient, limits, m);
        result.instant =
            partition_take(si::copy_rational_words(time.instant, plan, m.limits.cancellation));
        Projector projector{projection, projection.exact.steps};
        if (straight) {
            if (transient)
                partition_fail(Unsupported);
            m.admit();
            result.wake.emplace(prepare_straight_wake(source, root, wake, *straight, projector, m));
        }
        if (std::get<a::SingleWord>(retained(projector.apply(time.instant, projection.time,
                                                             core::Precision::Binary64, m)))
                .value != context.time)
            partition_fail(StaleGeneration);

        const auto model = root.model.record();
        const auto surfaces = model.aero_surfaces;
        if (surfaces.size() == 0)
            partition_fail(InvalidInput);
        result.patches = allocate<Patch>(m, surfaces.size());
        auto patches = result.patches.values();
        MeshCounts counts;
        std::uint64_t sharings = 0;
        for (std::size_t i = 0; i < patches.size(); ++i) {
            m.visit(1);
            const auto surface = partition_take(surfaces.at(partition_index(i))).record();
            if (straight &&
                !std::holds_alternative<s::AeroAttachmentFixed>(surface.attachment.record()))
                partition_fail(Unsupported);
            auto& patch = patches[i];
            const auto id = surface.geometry.unsigned_value();
            diagnostic = {};
            diagnostic.phase = CaseWingPhase::Reference;
            diagnostic.surface = surface.id.unsigned_value();
            const auto graph = find_graph(model, id, m);
            patch.chart.emplace(prepare_chart(source, surface, graph, m, diagnostic));
            counts.append(patch);
            sharings = partition_add(sharings, surface.sharing.size());
        }
        diagnostic = {};
        diagnostic.phase = CaseWingPhase::Geometry;
        counts.admit(limits.geometry.counts);
        m.leaf_limits = limits.chart.algebra;
        m.admit(partition_add(
            pg::reference_live_rationals + (straight ? 6 : 0),
            partition_add(partition_mul(2, counts.points), partition_mul(2, counts.cells))));
        result.points = allocate<CaseWingPoint>(m, counts.points);
        result.cells = allocate<CaseWingCell>(m, counts.cells);
        result.bands = allocate<CaseWingTrailingBand>(m, counts.bands);
        result.sharing = allocate<CaseWingSharing>(m, sharings);
        for (std::size_t i = 0; i < patches.size(); ++i) {
            const auto surface = partition_take(surfaces.at(partition_index(i))).record();
            sample_patch(source, surface, patches[i], time, limits, projector, curve, result, m,
                         diagnostic);
        }

        auto point_owners = allocate<std::size_t>(m, counts.points);
        auto edge_owners = allocate<std::size_t>(m, counts.edges);
        auto internal = allocate<bool>(m, counts.edges);
        const auto initialize = [&](auto& buffer) {
            for (std::size_t i = 0; i < buffer.size(); ++i) {
                m.visit(1);
                buffer.values()[i] = i;
            }
        };
        initialize(point_owners);
        initialize(edge_owners);
        diagnostic = {};
        diagnostic.phase = CaseWingPhase::Sharing;
        std::size_t share_index = 0;
        for (std::size_t i = 0; i < patches.size(); ++i) {
            m.visit(1);
            const auto surface = partition_take(surfaces.at(partition_index(i))).record();
            for (std::size_t k = 0; k < surface.sharing.size(); ++k) {
                m.visit(1);
                const auto seam = partition_take(surface.sharing.at(partition_index(k))).record();
                const auto id = seam.other.unsigned_value();
                const auto found = std::ranges::find_if(patches, [&](const Patch& patch) {
                    m.visit(1);
                    return patch.chart->surface_id == id;
                });
                if (found == patches.end() || &*found == &patches[i])
                    partition_fail(InvalidInput);
                const auto own_edge = seam.own_boundary.unsigned_value();
                const auto other_edge = seam.other_boundary.unsigned_value();
                const auto reversed = seam.orientation.text() == "reversed";
                diagnostic.surface = patches[i].chart->surface_id;
                diagnostic.other_surface = id;
                diagnostic.boundary = own_edge;
                diagnostic.other_boundary = other_edge;
                result.sharing.values()[share_index++] = {patches[i].chart->surface_id, own_edge,
                                                          id, other_edge, reversed};
                share(patches[i], own_edge, *found, other_edge, reversed, result.points.values(),
                      point_owners.values(), edge_owners.values(), internal.values(), m,
                      diagnostic);
            }
        }

        diagnostic = {};
        diagnostic.phase = CaseWingPhase::CoreRadius;
        const auto value = partition_take(source.quantity(wake.core_radius));
        const auto radius =
            partition_take(io::materialize_rational(value.recipe, m.leaf_limits, plan));
        if (radius.sign() <= 0)
            partition_fail(InvalidInput);
        const auto projected =
            projector.apply(radius, projection.core_radius, projection.representation, m);
        diagnostic.phase = CaseWingPhase::Geometry;
        auto totals = counts;
        auto terminal = allocate<std::size_t>(m, straight ? counts.points : 0);
        if (straight) {
            const auto wake_counts =
                count_wake(result, point_owners.values(), internal.values(), terminal.values(), m);
            totals.points = partition_add(totals.points, wake_counts.vertices);
            totals.cells = partition_add(totals.cells, wake_counts.bands);
            totals.edges =
                partition_add(totals.edges, partition_add(wake_counts.vertices,
                                                          partition_mul(2, wake_counts.bands)));
            totals.admit(limits.geometry.counts);
            (void)partition_index(partition_mul(2, totals.cells));
            result.wake_bands = allocate<CaseWingWakeBand>(m, wake_counts.bands);
            result.terminal_anchor = allocate<std::size_t>(m, wake_counts.vertices);
        }
        Mesh mesh{totals, straight ? core::Precision::TwoWord : projection.representation, m};
        mesh.place_points(result, point_owners.values(), m);
        for (const auto& patch : patches)
            mesh.connect(patch, result, edge_owners.values(), internal.values(), m);
        if (straight)
            mesh.extend(counts, terminal.values(), result, projector, m);
        mesh.publish(result, retained(projected), limits, projector, geometry_work, global_work, m,
                     diagnostic, topology ? &*topology : nullptr);
        m.cancel();
        return CaseWingGeometry{std::move(carrier)};
    }
};
core::Result<CaseWingGeometry>
prepare(const io::SiCase& source, const aero::detail::GeometryContext& context,
        std::optional<CaseFlowClockQuery> address, const CaseWingProjection& projection,
        const CaseWingLimits& limits, aero::detail::GeometryWorkspace& geometry,
        aero::detail::GlobalWorkspace& global, aero::detail::pchip::CurveWorkspace& curve,
        core::ResourcePlan& plan, CaseWingDiagnostic* diagnostic,
        std::optional<CaseWingStraightWake> straight, const Owner* anchor) noexcept {
    CaseWingDiagnostic location;
    auto result = core::guard_entry(
        wing_geometry::operation, core::FailureCode::PrecisionUnresolved, context.diagnostics,
        [&]() -> core::Result<CaseWingGeometry> {
            try {
                return Build::run(source, context, address, projection, limits, geometry, global,
                                  curve, plan, location, straight, anchor);
            } catch (const structure::detail::PartitionAbort& rejected) {
                auto inner = rejected.failure;
                inner.operation = wing_geometry::operation;
                if (!inner.details.component && location.surface)
                    inner.details.component = core::ComponentId{*location.surface};
                return std::unexpected(inner);
            }
        });
    if (diagnostic)
        *diagnostic = result ? CaseWingDiagnostic{} : location;
    return result;
}
} // namespace wing_geometry

CaseWingGeometry::CaseWingGeometry(core::PlannedBuffer<wing_geometry::Owner> owner) noexcept
    : owner_(std::move(owner)) {}
CaseWingGeometry::CaseWingGeometry(CaseWingGeometry&&) noexcept = default;
CaseWingGeometry& CaseWingGeometry::operator=(CaseWingGeometry&&) noexcept = default;
CaseWingGeometry::~CaseWingGeometry() = default;

core::Result<CaseWingGeometry> CaseWingGeometry::prepare(
    const io::SiCase& source, const aero::detail::GeometryContext& context,
    std::optional<CaseFlowClockQuery> address, const CaseWingProjection& projection,
    const CaseWingLimits& limits, aero::detail::GeometryWorkspace& geometry,
    aero::detail::GlobalWorkspace& global, aero::detail::pchip::CurveWorkspace& curve,
    core::ResourcePlan& plan, CaseWingDiagnostic* diagnostic,
    std::optional<CaseWingStraightWake> straight) noexcept {
    return wing_geometry::prepare(source, context, address, projection, limits, geometry, global,
                                  curve, plan, diagnostic, straight, nullptr);
}

core::Result<CaseWingGeometry> CaseWingGeometry::evaluate_stage(
    const io::SiCase& source, const aero::detail::GeometryContext& context,
    CaseFlowClockQuery address, const CaseWingProjection& projection, const CaseWingLimits& limits,
    aero::detail::GeometryWorkspace& geometry, aero::detail::GlobalWorkspace& global,
    aero::detail::pchip::CurveWorkspace& curve, core::ResourcePlan& plan,
    CaseWingDiagnostic* diagnostic) const noexcept {
    if (owner_.size() != 1) {
        if (diagnostic)
            *diagnostic = {};
        const auto code = limits.chart.cancellation.stop_requested()
                              ? core::FailureCode::Cancelled
                              : core::FailureCode::StaleGeneration;
        return std::unexpected(core::Failure{code, wing_geometry::operation, context.diagnostics});
    }
    return wing_geometry::prepare(source, context, address, projection, limits, geometry, global,
                                  curve, plan, diagnostic, std::nullopt, &owner_.values().front());
}

core::Result<void>
CaseWingGeometry::validate_for(const io::SiCase& source,
                               const aero::detail::GeometryContext& context,
                               std::optional<CaseFlowClockQuery> address) const noexcept {
    const auto checked = source.source_case();
    const auto identity = checked ? checked->get().identity() : std::unexpected(checked.error());
    if (owner_.size() == 1 && identity) {
        const auto& value = owner_.values().front();
        if (wing_geometry::same_source(*value.source, *identity) && *value.context == context &&
            wing_geometry::same_address(value.address, address))
            return {};
    }
    return std::unexpected(core::Failure{core::FailureCode::StaleGeneration,
                                         wing_geometry::operation, context.diagnostics});
}

core::Result<io::SchemaIdentity> CaseWingGeometry::source() const noexcept {
    if (owner_.size() == 1)
        return *owner_.values().front().source;
    return std::unexpected(
        core::Failure{core::FailureCode::StaleGeneration, wing_geometry::operation});
}

std::optional<CaseFlowClockQuery> CaseWingGeometry::clock_address() const noexcept {
    return owner_.size() == 1 ? owner_.values().front().address : std::nullopt;
}

std::size_t CaseWingGeometry::surface_count() const noexcept {
    return owner_.size() == 1 ? owner_.values().front().patches.size() : 0;
}

core::Result<CaseWingSurfaceView> CaseWingGeometry::surface(std::size_t i) const& noexcept {
    if (i >= surface_count())
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, wing_geometry::operation});
    const auto& chart = *owner_.values().front().patches.values()[i].chart;
    const auto identity = chart.graph.identity();
    if (!identity)
        return std::unexpected(identity.error());
    const auto range = [&](std::uint64_t index) noexcept {
        const auto selected = chart.boundaries.values()[static_cast<std::size_t>(index)];
        return chart.boundary_points.values().subspan(selected.offset, selected.count);
    };
    return CaseWingSurfaceView{chart.surface_id,
                               identity->geometry_id,
                               chart.leading_index,
                               chart.trailing_index,
                               chart.reversed,
                               chart.chord.values(),
                               chart.span.values(),
                               chart.knots.values().first(chart.knot_count),
                               range(chart.leading_index),
                               range(chart.trailing_index)};
}

std::span<const CaseWingSharing> CaseWingGeometry::sharing() const& noexcept {
    return owner_.size() == 1 ? owner_.values().front().sharing.values()
                              : std::span<const CaseWingSharing>{};
}

std::span<const CaseWingPoint> CaseWingGeometry::points() const& noexcept {
    return owner_.size() == 1 ? owner_.values().front().points.values()
                              : std::span<const CaseWingPoint>{};
}

core::Result<DiscreteCollocation> CaseWingGeometry::collocation(std::size_t i) const& noexcept {
    return core::guard_entry(
        wing_geometry::operation, core::FailureCode::PrecisionUnresolved, {},
        [&]() -> core::Result<DiscreteCollocation> {
            try {
                const auto rows = cells();
                if (i >= rows.size())
                    wing_geometry::partition_fail(core::FailureCode::InvalidInput);
                const auto& row = rows[i];
                if (!row.collocation)
                    wing_geometry::partition_fail(core::FailureCode::StaleGeneration);
                const auto words = [](const auto& values) -> NonpenetrationVector {
                    using wing_geometry::retained;
                    return {retained(values[0]), retained(values[1]), retained(values[2])};
                };
                auto span = words(row.span_tangent);
                if (row.reversed)
                    wing_geometry::reverse_tangent(span);
                return DiscreteCollocation{row.facets[0],
                                           *row.collocation,
                                           {words(row.chord_tangent), span},
                                           words(row.velocity),
                                           aero::detail::SingleWord{0.0},
                                           aero::detail::SingleWord{0.0}};
            } catch (const structure::detail::PartitionAbort& error) {
                return std::unexpected(error.failure);
            }
        });
}

std::span<const CaseWingCell> CaseWingGeometry::cells() const& noexcept {
    return owner_.size() == 1 ? owner_.values().front().cells.values()
                              : std::span<const CaseWingCell>{};
}

std::span<const CaseWingTrailingBand> CaseWingGeometry::trailing_bands() const& noexcept {
    return owner_.size() == 1
               ? owner_.values().front().bands.values().first(owner_.values().front().band_count)
               : std::span<const CaseWingTrailingBand>{};
}

std::optional<CaseWingStraightWakeView> CaseWingGeometry::straight_wake() const& noexcept {
    if (owner_.size() != 1 || !owner_.values().front().wake)
        return std::nullopt;
    const auto& wake = *owner_.values().front().wake;
    return CaseWingStraightWakeView{wake.displacement, wake.velocity, wake.duration};
}

std::span<const CaseWingWakeBand> CaseWingGeometry::wake_bands() const& noexcept {
    return owner_.size() == 1 ? owner_.values().front().wake_bands.values()
                              : std::span<const CaseWingWakeBand>{};
}

core::Result<std::reference_wrapper<const aero::detail::GlobalContactState>>
CaseWingGeometry::global() const& noexcept {
    if (owner_.size() == 1)
        return std::cref(*owner_.values().front().global);
    return std::unexpected(
        core::Failure{core::FailureCode::StaleGeneration, wing_geometry::operation});
}

core::Result<std::reference_wrapper<const aero::detail::AliasKinematicsState>>
CaseWingGeometry::kinematics() const& noexcept {
    if (owner_.size() == 1)
        return std::cref(*owner_.values().front().kinematics);
    return std::unexpected(
        core::Failure{core::FailureCode::StaleGeneration, wing_geometry::operation});
}

} // namespace uvlas::coupling::detail
