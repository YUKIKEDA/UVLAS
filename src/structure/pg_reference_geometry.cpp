#include "pg_reference_domain.hpp"
#include "pg_reference_hash.hpp"

#include <array>
#include <span>
#include <utility>
#include <uvlas/structure/pg_reference_basis.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;

[[nodiscard]] std::array<R, 3> normal_numerator(ReferenceArithmetic& m, const PgGraphView& graph,
                                                const R& fs, const R& ft, bool derivative) {
    auto x = m.div(fs, graph.width);
    auto y = m.sub(m.mul(graph.skew, x), ft);
    return {m.neg(x), std::move(y), derivative ? R{} : m.constant("1")};
}

[[nodiscard]] R dot(ReferenceArithmetic& m, const std::array<R, 3>& a, const std::array<R, 3>& b) {
    R result;
    for (std::size_t i = 0; i < a.size(); ++i)
        result = m.add(result, m.mul(a[i], b[i]));
    return result;
}

struct ReferenceShape {
    PgReferenceScalar value;
    PgReferenceVector chart_derivatives;
};

template <std::size_t Orders>
[[nodiscard]] ReferenceShape reference_shape(ReferenceArithmetic& m,
                                             const std::array<std::array<R, Orders>, 3>& numerators,
                                             const PgReferenceValue& point,
                                             const std::array<R, 2>& g_gradient, const R& length) {
    static_assert(Orders >= 3);
    ReferenceShape result;
    const auto& [r, p, c] = numerators;
    const auto two = m.constant("2"), three = m.constant("3");
    const auto two_length = m.mul(two, length);
    const auto z_zeta = m.div(point.full_thickness_metres, two_length);
    const auto z = m.mul(z_zeta, point.zeta);
    const auto half_z_squared = m.div(m.square(z), two);
    const auto g_squared = m.square(point.g);
    const auto director = m.div(p[0], g_squared);
    result.value = {m.add(r[0], m.mul(half_z_squared, c[0])), m.mul(z, director)};
    for (std::size_t i = 0; i < g_gradient.size(); ++i) {
        const auto z_a = m.div(m.mul(point.thickness_gradient_metres[i], point.zeta), two_length);
        const auto numerator =
            m.sub(p[i + 1], m.div(m.mul(three, m.mul(p[0], g_gradient[i])), m.mul(two, point.g)));
        const auto d_a = m.div(numerator, g_squared);
        result.chart_derivatives[i] = {
            m.add(m.add(r[i + 1], m.mul(half_z_squared, c[i + 1])), m.mul(m.mul(z, z_a), c[0])),
            m.add(m.mul(z, d_a), m.mul(z_a, director)),
        };
    }
    result.chart_derivatives[2] = {m.mul(m.mul(z, z_zeta), c[0]), m.mul(z_zeta, director)};
    return result;
}

[[nodiscard]] PgReferenceScalar inverse_determinant(ReferenceArithmetic& m, const R& g,
                                                    const PgReferenceScalar& determinant) {
    const auto& [a, b] = determinant;
    const auto one = m.constant("1");
    if (b.sign() == 0)
        return {m.div(one, a), {}};
    if (a.sign() == 0)
        return {{}, m.div(one, m.mul(b, g))};
    const auto norm = m.sub(m.square(a), m.mul(g, m.square(b)));
    if (norm.sign() == 0) {
        // A nonzero determinant and zero conjugate imply b sqrt(g) = a.
        return {m.div(one, m.add(a, a)), {}};
    }
    return {m.div(a, norm), m.div(m.neg(b), norm)};
}

[[nodiscard]] PgReferenceVector physical_gradient(ReferenceArithmetic& m,
                                                  const PgReferenceValue& point,
                                                  const PgReferenceVector& chart_derivatives) {
    auto field = ReferenceRadicals{m, point.g};
    const auto& jacobian = point.jacobian_columns_metres;
    const auto inverse = inverse_determinant(m, point.g, point.determinant_metres_cubed);
    PgReferenceVector gradient;
    for (std::size_t d = 0; d < chart_derivatives.size(); ++d) {
        const auto dual = field.cross(jacobian[(d + 1) % 3], jacobian[(d + 2) % 3]);
        for (std::size_t j = 0; j < gradient.size(); ++j)
            gradient[j] = field.add(gradient[j], field.mul(chart_derivatives[d], dual[j]));
    }
    for (auto& value : gradient)
        value = field.mul(value, inverse);
    return gradient;
}

[[nodiscard]] PgReferenceColumnValue
column_at(ReferenceArithmetic& m, const PgCoefficientDictionary& dictionary,
          const PgReferenceValue& point, const std::array<R, 2>& g_gradient,
          const PgCoefficientColumn column, const R& axis_scale, const R& length) {
    const auto u = m.sub(m.mul(axis_scale, point.s), index_value(m, 2 * point.cell.first + 1));
    const auto v = m.sub(m.mul(axis_scale, point.t), index_value(m, 2 * point.cell.second + 1));
    const auto g_poly = partition_take(dictionary.g(point.cell));
    if (polynomial_derivative(m, *g_poly, {}, u, v) != point.g)
        partition_fail(core::FailureCode::StaleGeneration);

    constexpr std::array fields{PgCoefficientField::R, PgCoefficientField::P,
                                PgCoefficientField::C};
    const auto scale_squared = m.square(axis_scale);
    PgReferenceColumnValue result{column, {}, {}, {}};
    for (std::size_t field = 0; field < fields.size(); ++field) {
        const auto polynomial =
            partition_take(dictionary.coefficient(column, point.cell, fields[field]));
        for (std::size_t i = 0; i < result.numerators[field].size(); ++i) {
            const auto order = reference_orders[i];
            auto value = polynomial_derivative(m, *polynomial, order, u, v);
            if (i > 0)
                value = m.mul(value, i <= 2 ? axis_scale : scale_squared);
            result.numerators[field][i] = std::move(value);
        }
    }
    auto position = reference_shape(m, result.numerators, point, g_gradient, length);
    result.value = std::move(position.value);
    result.chart_derivatives = std::move(position.chart_derivatives);
    return result;
}

[[nodiscard]] PgReferenceValue
evaluate_reference(ReferenceArithmetic& m, const PgGraphView& graph, const PgReferenceInputs& input,
                   const PgReferenceFrame& frame, const PgReferencePoint& query,
                   std::optional<PgCoefficientColumn> column, const R* determinant_bound,
                   std::array<R, 2>* point_g_gradient = nullptr) {
    m.visit(1);
    const auto chart = std::visit(
        [&](const auto& point) {
            if constexpr (std::is_same_v<std::decay_t<decltype(point)>, PgReferenceLateralPoint>)
                return lateral_point(m, input.partition, point);
            else
                return triangle_point(m, input.partition, point);
        },
        query);
    const auto& piece = reference_at(partition_take(input.partition.pieces()), chart.piece);
    const auto plane =
        thickness_plane(m, partition_take(input.partition.thickness()), piece.thickness_triangle);

    PgReferenceValue out{};
    out.cell = piece.cell;
    out.piece = chart.piece;
    out.thickness_triangle = piece.thickness_triangle;
    out.partition_record = chart.record;
    out.layer = chart.layer.id;
    out.s = m.copy(chart.s);
    out.t = m.copy(chart.t);
    out.zeta = m.copy(chart.zeta);
    out.full_thickness_metres = plane.evaluate(m, chart.s, chart.t);
    if (out.full_thickness_metres.sign() <= 0)
        partition_fail(core::FailureCode::GeometryInvalid);
    for (std::size_t d = 0; d < out.thickness_gradient_metres.size(); ++d)
        out.thickness_gradient_metres[d] = m.copy(plane.gradient[d]);
    for (std::size_t i = 0; i < reference_orders.size(); ++i)
        out.height_chart_jet[i] =
            polynomial_derivative(m, graph.height, reference_orders[i], chart.s, chart.t);
    const auto& f = out.height_chart_jet;
    const auto n = normal_numerator(m, graph, f[1], f[2], false);
    const std::array n_derivatives{
        normal_numerator(m, graph, f[3], f[4], true),
        normal_numerator(m, graph, f[4], f[5], true),
    };
    out.g = dot(m, n, n);
    const auto two = m.constant("2");
    const auto two_g = m.mul(two, out.g);
    const auto rotated_n = physical_rotate(m, frame.rotation, n);
    std::array<R, 2> g_gradient;
    std::array<std::array<R, 3>, 2> normal_b_gradient;
    for (std::size_t d = 0; d < g_gradient.size(); ++d) {
        g_gradient[d] = m.mul(two, dot(m, n, n_derivatives[d]));
        const auto factor = m.div(g_gradient[d], two_g);
        const auto rotated_n_d = physical_rotate(m, frame.rotation, n_derivatives[d]);
        for (std::size_t i = 0; i < rotated_n.size(); ++i) {
            normal_b_gradient[d][i] =
                m.div(m.sub(rotated_n_d[i], m.mul(rotated_n[i], factor)), out.g);
        }
    }
    const auto position = physical_rotate(
        m, frame.rotation,
        {m.add(graph.x0, m.add(m.mul(graph.width, chart.s), m.mul(graph.skew, chart.t))),
         m.copy(chart.t), m.copy(f[0])});
    const std::array midsurface_tangents{
        physical_rotate(m, frame.rotation, {m.copy(graph.width), R{}, m.copy(f[1])}),
        physical_rotate(m, frame.rotation, {m.copy(graph.skew), m.constant("1"), m.copy(f[2])}),
    };
    const auto half_h = m.div(out.full_thickness_metres, two);
    const auto z = m.mul(half_h, chart.zeta);
    const auto half_zeta = m.div(chart.zeta, two);
    for (std::size_t i = 0; i < out.normal.size(); ++i) {
        const auto normal_b = m.div(rotated_n[i], out.g);
        out.normal[i] = {R{}, m.copy(normal_b)};
        out.position_metres[i] = {
            m.add(frame.origin_metres[i], m.mul(frame.length_metres, position[i])),
            m.mul(z, normal_b)};
        for (std::size_t d = 0; d < out.normal_chart_derivatives.size(); ++d) {
            out.normal_chart_derivatives[d][i] = {R{}, m.copy(normal_b_gradient[d][i])};
            out.jacobian_columns_metres[d][i] = {
                m.mul(frame.length_metres, midsurface_tangents[d][i]),
                m.add(m.mul(m.mul(half_zeta, plane.gradient[d]), normal_b),
                      m.mul(z, normal_b_gradient[d][i])),
            };
        }
        out.jacobian_columns_metres[2][i] = {R{}, m.mul(half_h, normal_b)};
    }

    auto field = ReferenceRadicals{m, out.g};
    const auto& j = out.jacobian_columns_metres;
    auto chart_area = field.cross(j[0], j[1]);
    out.determinant_metres_cubed = field.dot(chart_area, j[2]);
    if (field.sign(out.determinant_metres_cubed) <= 0)
        partition_fail(core::FailureCode::GeometryInvalid);
    if (determinant_bound != nullptr &&
        field.sign({m.sub(out.determinant_metres_cubed.rational, *determinant_bound),
                    m.copy(out.determinant_metres_cubed.radical)}) < 0)
        partition_fail(core::FailureCode::PrecisionUnresolved);

    if (std::holds_alternative<PgReferenceVolumePoint>(query)) {
        out.measure =
            PgReferenceVolumeMeasure{field.scale(out.determinant_metres_cubed, chart.map.scale)};
    } else {
        const bool lateral = std::holds_alternative<PgReferenceLateralPoint>(query);
        if (lateral) {
            PgReferenceVector tangent;
            for (std::size_t i = 0; i < tangent.size(); ++i)
                tangent[i] = field.add(field.scale(j[0][i], chart.map.tangent[0]),
                                       field.scale(j[1][i], chart.map.tangent[1]));
            chart_area = field.cross(tangent, j[2]);
        }
        const auto factor = chart.map.side < 0 ? m.neg(chart.map.scale) : m.copy(chart.map.scale);
        for (auto& x : chart_area)
            x = field.scale(x, factor);
        auto squared = field.dot(chart_area, chart_area);
        if (field.sign(squared) <= 0)
            partition_fail(core::FailureCode::GeometryInvalid);
        auto area = PgReferenceAreaMeasure{std::move(chart_area), std::move(squared)};
        if (lateral) {
            out.measure = PgReferenceLateralMeasure{std::move(area)};
        } else {
            const auto side =
                chart.map.side < 0 ? PgPartitionFaceSide::Lower : PgPartitionFaceSide::Upper;
            out.measure = PgReferenceHorizontalMeasure{side, std::move(area)};
        }
    }

    if (column)
        out.column = column_at(m, input.dictionary, out, g_gradient, *column,
                               index_value(m, 4 * graph.n), frame.length_metres);
    if (point_g_gradient != nullptr)
        *point_g_gradient = std::move(g_gradient);
    m.cancel();
    out.call_info = m.info();
    return out;
}

} // namespace

PgReferenceGeometry::PgReferenceGeometry(PgReferenceInputs&& source, std::array<R, 9> rotation,
                                         std::array<R, 3> origin, R length, R determinant_bound,
                                         PgReferenceIdentity identity,
                                         PgReferenceCallInfo info) noexcept
    : inputs_(std::move(source)),
      rotation_(std::move(rotation)),
      origin_(std::move(origin)),
      length_(std::move(length)),
      determinant_bound_(std::move(determinant_bound)),
      identity_(identity),
      info_(info) {}

core::Result<PgReferenceGeometry> PgReferenceGeometry::prepare(PgReferenceInputs&& inputs,
                                                               const PgReferenceFrame& frame,
                                                               const PgReferenceLimits& limits,
                                                               core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgReferenceGeometry> {
        auto math = ReferenceArithmetic{limits, plan};
        math.admit();
        const auto graph = validate_reference_inputs(math, inputs);
        validate_frame(math, frame);
        auto bound = normal_tube_bound(math, graph, inputs.partition, frame.length_metres);
        std::array<R, 9> rotation;
        std::array<R, 3> origin;
        for (std::size_t i = 0; i < rotation.size(); ++i)
            rotation[i] = math.copy(frame.rotation[i]);
        for (std::size_t i = 0; i < origin.size(); ++i)
            origin[i] = math.copy(frame.origin_metres[i]);
        auto length = math.copy(frame.length_metres);
        const auto identity = make_reference_identity(math, graph.identity,
                                                      partition_take(inputs.partition.identity()),
                                                      {rotation, origin, length});
        math.cancel();
        return PgReferenceGeometry{std::move(inputs), std::move(rotation), std::move(origin),
                                   std::move(length), std::move(bound),    identity,
                                   math.info()};
    });
}

core::Result<PgReferenceValue> PgReferenceGeometry::evaluate(
    const PgReferencePoint& query, std::optional<PgCoefficientColumn> column,
    const PgReferenceLimits& limits, core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgReferenceValue> {
        const auto graph = partition_take(inputs_.graph.view());
        auto math = ReferenceArithmetic{limits, plan};
        math.admit();
        auto output = evaluate_reference(math, graph, inputs_, {rotation_, origin_, length_}, query,
                                         column, &determinant_bound_);
        output.identity = identity_;
        return output;
    });
}

core::Result<void> PgReferenceGeometry::diagnose_point(const PgReferenceInputs& inputs,
                                                       const PgReferenceFrame& frame,
                                                       const PgReferencePoint& query,
                                                       const PgReferenceLimits& limits,
                                                       core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<void> {
        auto math = ReferenceArithmetic{limits, plan};
        math.admit();
        const auto graph = validate_reference_inputs(math, inputs);
        validate_frame(math, frame);
        static_cast<void>(
            evaluate_reference(math, graph, inputs, frame, query, std::nullopt, nullptr));
        return {};
    });
}

core::Result<PgReferenceIdentity> PgReferenceGeometry::identity() const noexcept {
    if (length_.sign() == 0)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, partition_operation});
    return identity_;
}

PgReferenceCallInfo PgReferenceGeometry::preparation() const noexcept {
    return info_;
}

core::Result<PgReferenceSource> PgReferenceGeometry::source() const& noexcept {
    if (length_.sign() == 0)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, partition_operation});
    return PgReferenceSource{inputs_.graph, inputs_.dictionary, inputs_.partition};
}

core::Result<const R*> PgReferenceGeometry::determinant_lower_bound() const& noexcept {
    if (length_.sign() == 0)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, partition_operation});
    return &determinant_bound_;
}

PgReferenceBasis::PgReferenceBasis(PgReferenceGeometry&& geometry,
                                   PgRankPreparation&& rank) noexcept
    : geometry_(std::move(geometry)), rank_(std::move(rank)) {}

core::Result<PgReferenceBasis> PgReferenceBasis::prepare(PgReferenceGeometry&& geometry,
                                                         PgRankPreparation&& rank,
                                                         const PgRankPreparationLimits& limits,
                                                         core::Workspace& workspace,
                                                         core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgReferenceBasis> {
        const auto source = partition_take(geometry.source());
        partition_take(
            rank.validate_for(source.dictionary, source.partition, limits, workspace, plan));
        return PgReferenceBasis{std::move(geometry), std::move(rank)};
    });
}

BasisReferencePoint PgReferenceBasis::locate(ReferenceArithmetic& m,
                                             const PgReferencePoint& query) const {
    const auto graph = partition_take(geometry_.inputs_.graph.view());
    const PgReferenceFrame frame{geometry_.rotation_, geometry_.origin_, geometry_.length_};
    std::array<R, 2> g_gradient;
    auto point = evaluate_reference(m, graph, geometry_.inputs_, frame, query, std::nullopt,
                                    &geometry_.determinant_bound_, &g_gradient);
    point.identity = geometry_.identity_;
    auto axis = index_value(m, 4 * graph.n);
    auto u = m.sub(m.mul(axis, point.s), index_value(m, 2 * point.cell.first + 1));
    auto v = m.sub(m.mul(axis, point.t), index_value(m, 2 * point.cell.second + 1));
    return {std::move(point), std::move(g_gradient), std::move(axis), std::move(u), std::move(v)};
}

std::array<std::array<ExactRational, 3>, 3>
PgReferenceBasis::numerator_jet(ReferenceArithmetic& m, PgRankBasisColumn column,
                                const BasisReferencePoint& point, bool derivatives,
                                PgRankMonomialPowers* powers) const {
    return rank_.basis_cell_jet(m, column, point.reference.cell, point.u, point.v, point.axis_scale,
                                derivatives, powers);
}

BasisReferenceSample
PgReferenceBasis::sample(ReferenceArithmetic& m, const BasisReferencePoint& point,
                         const std::array<std::array<ExactRational, 3>, 3>& numerators) const {
    auto shape =
        reference_shape(m, numerators, point.reference, point.g_gradient, geometry_.length_);
    auto gradient = physical_gradient(m, point.reference, shape.chart_derivatives);
    return {std::move(shape.value), std::move(shape.chart_derivatives), std::move(gradient)};
}

PgReferenceBasisValue
PgReferenceBasis::assemble(ReferenceArithmetic& m, BasisReferencePoint&& point,
                           PgRankBasisColumn column,
                           const std::array<std::array<ExactRational, 3>, 3>& numerators) const {
    auto shape = sample(m, point, numerators);
    m.cancel();
    point.reference.call_info = m.info();
    return {std::move(point.reference), column, std::move(shape.value),
            std::move(shape.chart_derivatives), std::move(shape.physical_gradient_per_metre)};
}

core::Result<PgReferenceBasisValue>
PgReferenceBasis::evaluate(const PgReferencePoint& query, PgRankBasisColumn column,
                           const PgReferenceLimits& limits,
                           core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgReferenceBasisValue> {
        auto m = ReferenceArithmetic{limits, plan};
        m.admit(basis_live_rationals);
        const auto rank = partition_take(rank_.preparation()).rank;
        if (column.value < 0 || static_cast<std::uint64_t>(column.value) >= rank)
            partition_fail(core::FailureCode::InvalidInput);
        auto point = locate(m, query);
        const auto numerators = numerator_jet(m, column, point);
        return assemble(m, std::move(point), column, numerators);
    });
}

core::Result<PgReferenceFrame> PgReferenceGeometry::frame() const& noexcept {
    if (const auto present = identity(); !present)
        return std::unexpected(present.error());
    return PgReferenceFrame{rotation_, origin_, length_};
}

core::Result<PgReferenceIdentity> PgReferenceBasis::reference_identity() const noexcept {
    return geometry_.identity();
}

core::Result<PgReferenceFrame> PgReferenceBasis::frame() const& noexcept {
    return geometry_.frame();
}

core::Result<std::uint64_t> PgReferenceBasis::size() const noexcept {
    const auto info = rank_.preparation();
    if (!info)
        return std::unexpected(info.error());
    return info->rank;
}

core::Result<const PgRankPreparation*> PgReferenceBasis::rank() const& noexcept {
    const auto info = rank_.preparation();
    if (!info)
        return std::unexpected(info.error());
    return &rank_;
}

const ExactRational& PgReferenceBasis::length_metres() const& noexcept {
    return geometry_.length_;
}

core::Result<PgReferenceSource> PgReferenceBasis::source() const& noexcept {
    return geometry_.source();
}

} // namespace uvlas::structure
