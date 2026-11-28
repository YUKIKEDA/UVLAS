#include "case_reference_graph.hpp"

#include "structure/pg_reference_polynomials.hpp"

#include <utility>
#include <variant>

namespace uvlas::coupling::detail {
namespace pg = structure::detail;
namespace {
using core::FailureCode;
using pg::partition_fail;
using pg::partition_take;
using structure::ExactRational;
using structure::PgPartitionPoint;

constexpr io::SchemaDimension dimensionless{{0, 0, 0, 0}};
constexpr io::SchemaDimension metres{{1, 0, 0, 0}};
constexpr std::uint64_t fixed_live_rationals = 64;

template <class Action>
auto run(core::Operation operation, Action action) noexcept -> std::invoke_result_t<Action> {
    return core::guard_entry(
        operation, FailureCode::PrecisionUnresolved, {}, [&]() -> std::invoke_result_t<Action> {
            try {
                return action();
            } catch (const pg::PartitionAbort& rejected) {
                const auto& inner = rejected.failure;
                return std::unexpected(core::Failure{inner.code, operation, inner.details});
            }
        });
}

template <class Field>
ExactRational exact(const io::SiCase& source, Field field, io::SchemaDimension expected,
                    pg::ReferenceArithmetic& arithmetic) {
    arithmetic.charge();
    const auto quantity = partition_take(source.quantity(field));
    if (quantity.recipe.dimensions.exponents != expected.exponents || quantity.recipe.pi_power != 0)
        partition_fail(FailureCode::InvalidInput);
    const auto input_digits = std::visit(
        [](const auto& value) -> std::uint64_t {
            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(value)>, io::SiDecimalView>)
                return value.digits.size();
            else
                return pg::partition_add(value.numerator.digits.size(),
                                         value.denominator.digits.size());
        },
        quantity.recipe.value);
    arithmetic.charge(input_digits);
    return partition_take(
        io::materialize_rational(quantity.recipe, arithmetic.leaf_limits, arithmetic.plan));
}

int orient(pg::ReferenceArithmetic& m, const PgPartitionPoint& a, const PgPartitionPoint& b,
           const PgPartitionPoint& p) {
    return m.sub(m.mul(m.sub(b.s, a.s), m.sub(p.t, a.t)), m.mul(m.sub(b.t, a.t), m.sub(p.s, a.s)))
        .sign();
}

bool on_axis(pg::ReferenceArithmetic& m, const ExactRational& a, const ExactRational& b,
             const ExactRational& x) {
    const auto side_a = m.compare(x, a), side_b = m.compare(x, b);
    return side_a == 0 || side_b == 0 || side_a != side_b;
}

bool on_segment(pg::ReferenceArithmetic& m, const PgPartitionPoint& a, const PgPartitionPoint& b,
                const PgPartitionPoint& p) {
    return on_axis(m, a.s, b.s, p.s) && on_axis(m, a.t, b.t, p.t);
}

bool intersects(pg::ReferenceArithmetic& m, const PgPartitionPoint& a, const PgPartitionPoint& b,
                const PgPartitionPoint& c, const PgPartitionPoint& d) {
    const auto abc = orient(m, a, b, c), abd = orient(m, a, b, d);
    const auto cda = orient(m, c, d, a), cdb = orient(m, c, d, b);
    return (abc * abd < 0 && cda * cdb < 0) || (abc == 0 && on_segment(m, a, b, c)) ||
           (abd == 0 && on_segment(m, a, b, d)) || (cda == 0 && on_segment(m, c, d, a)) ||
           (cdb == 0 && on_segment(m, c, d, b));
}

void validate_domain(pg::ReferenceArithmetic& m, std::span<const PgPartitionPoint> polygon) {
    if (polygon.size() < 3)
        partition_fail(FailureCode::GeometryInvalid);
    const auto one = m.constant("1");
    ExactRational twice_area;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        m.visit(1);
        const auto next = (i + 1) % polygon.size();
        const auto& a = polygon[i];
        const auto& b = polygon[next];
        if (a.s.sign() < 0 || a.t.sign() < 0 || m.compare(a.s, one) > 0 ||
            m.compare(a.t, one) > 0 || (a.s == b.s && a.t == b.t))
            partition_fail(FailureCode::GeometryInvalid);
        twice_area = m.add(twice_area, m.sub(m.mul(a.s, b.t), m.mul(a.t, b.s)));
        for (std::size_t j = i + 1; j < polygon.size(); ++j) {
            m.visit(1);
            const auto end = (j + 1) % polygon.size();
            const auto& c = polygon[j];
            const auto& d = polygon[end];
            const auto overlap = [&](const auto& first, const auto& joint, const auto& last) {
                return orient(m, first, joint, last) == 0 &&
                       (on_segment(m, first, joint, last) || on_segment(m, joint, last, first));
            };
            if (next == j ? overlap(a, b, d)
                          : (end == i ? overlap(c, d, b) : intersects(m, a, b, c, d)))
                partition_fail(FailureCode::GeometryInvalid);
        }
    }
    if (twice_area.sign() == 0)
        partition_fail(FailureCode::GeometryInvalid);
}

bool contains(pg::ReferenceArithmetic& m, std::span<const PgPartitionPoint> polygon,
              const PgPartitionPoint& p) {
    bool inside = false;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        m.visit(1);
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const auto side = orient(m, a, b, p);
        if (side == 0 && on_segment(m, a, b, p))
            return true;
        const bool below_a = m.compare(a.t, p.t) <= 0;
        const bool below_b = m.compare(b.t, p.t) <= 0;
        if (below_a != below_b && (side > 0) == below_a)
            inside = !inside;
    }
    return inside;
}

structure::ExactChartPolynomial height(const io::SiCase& source,
                                       const io::schema::GeometrySurfaceGraph& record,
                                       pg::ReferenceArithmetic& m) {
    constexpr auto maximum_terms = pg::reference_orders.size();
    const auto count = record.coefficients.size();
    if (count > maximum_terms)
        partition_fail(FailureCode::InvalidInput);
    std::array<structure::ExactPowers, maximum_terms> powers{};
    std::array<ExactRational, maximum_terms> coefficients;
    for (std::size_t i = 0; i < count; ++i) {
        const auto term = partition_take(record.coefficients.at(pg::partition_index(i))).record();
        powers[i] = {pg::partition_index(partition_take(term.powers.at(0)).unsigned_value()),
                     pg::partition_index(partition_take(term.powers.at(1)).unsigned_value())};
        coefficients[i] = exact(source, term.coefficient, dimensionless, m);
    }
    const auto terms = [&]<std::size_t... I>(std::index_sequence<I...>) {
        return std::array{structure::ExactInputTerm{powers[I], coefficients[I]}...};
    }(std::make_index_sequence<maximum_terms>{});
    m.charge(pg::partition_mul(count, count));
    return partition_take(structure::ExactChartPolynomial::make(
        std::span{terms}.first(static_cast<std::size_t>(count)), m.leaf_limits, m.plan));
}
} // namespace

struct CaseReferenceSurface::Data {
    CaseReferenceSurfaceIdentity identity{{{io::SourceId{0}, 0, {}}, 0, 0}, 0};
    ExactRational x0, width, skew, length;
    std::array<ExactRational, 3> origin;
    std::array<ExactRational, 9> rotation;
    structure::ExactChartPolynomial height;
    core::PlannedBuffer<PgPartitionPoint> domain;
    structure::PgReferenceCallInfo preparation;
};

CaseReferenceSurface::CaseReferenceSurface(core::PlannedBuffer<Data> data) noexcept
    : data_(std::move(data)) {}
CaseReferenceSurface::CaseReferenceSurface(CaseReferenceSurface&&) noexcept = default;
CaseReferenceSurface& CaseReferenceSurface::operator=(CaseReferenceSurface&&) noexcept = default;
CaseReferenceSurface::~CaseReferenceSurface() = default;

core::Result<CaseReferenceSurface>
CaseReferenceSurface::prepare(const io::SiCase& source, std::uint64_t geometry_id,
                              const structure::PgReferenceLimits& limits,
                              core::ResourcePlan& plan) noexcept {
    constexpr auto operation = core::Operation::PrepareModel;
    return run(operation, [&]() -> core::Result<CaseReferenceSurface> {
        pg::ReferenceArithmetic m{limits, plan};
        m.cancel();
        const auto& owner = partition_take(source.source_case()).get();
        const auto geometries = partition_take(owner.root()).model.record().geometry;
        for (std::uint64_t i = 0; i < geometries.size(); ++i) {
            m.visit(1);
            const auto definition = partition_take(geometries.at(pg::partition_index(i))).record();
            if (std::visit([](const auto& value) { return value.id.unsigned_value(); },
                           definition) != geometry_id)
                continue;
            const auto* graph = std::get_if<io::schema::GeometrySurfaceGraph>(&definition);
            if (!graph)
                partition_fail(FailureCode::Unsupported);
            m.admit(pg::partition_add(fixed_live_rationals,
                                      pg::partition_mul(2, graph->domain.size())));
            auto data = partition_take(plan.allocate<Data>(1, operation));
            auto& result = data.values().front();
            result.identity = {partition_take(owner.identity()), geometry_id};
            result.length = exact(source, graph->reference_length, metres, m);
            result.x0 = exact(source, graph->x0, dimensionless, m);
            result.width = exact(source, graph->width, dimensionless, m);
            result.skew = exact(source, graph->skew, dimensionless, m);
            if (result.width.sign() <= 0)
                partition_fail(FailureCode::GeometryInvalid);
            result.height = height(source, *graph, m);
            const auto frame = graph->frame.record();
            for (std::size_t j = 0; j < result.origin.size(); ++j) {
                const auto row_index = static_cast<core::Index>(j);
                result.origin[j] =
                    exact(source, partition_take(frame.origin.at(row_index)), metres, m);
                const auto row = partition_take(frame.rotation.at(row_index));
                for (std::size_t k = 0; k < result.origin.size(); ++k)
                    result.rotation[j * result.origin.size() + k] =
                        exact(source, partition_take(row.at(static_cast<core::Index>(k))),
                              dimensionless, m);
            }
            pg::validate_frame(m, {result.rotation, result.origin, result.length});
            result.domain = partition_take(plan.allocate<PgPartitionPoint>(
                pg::partition_index(graph->domain.size()), operation));
            for (std::size_t j = 0; j < result.domain.size(); ++j) {
                m.visit(1);
                const auto pair = partition_take(graph->domain.at(pg::partition_index(j)));
                result.domain.values()[j] = {
                    exact(source, partition_take(pair.at(0)), dimensionless, m),
                    exact(source, partition_take(pair.at(1)), dimensionless, m)};
            }
            validate_domain(m, result.domain.values());
            m.cancel();
            result.preparation = m.info();
            return CaseReferenceSurface{std::move(data)};
        }
        partition_fail(FailureCode::InvalidInput);
    });
}

core::Result<CaseReferenceSurfaceIdentity> CaseReferenceSurface::identity() const noexcept {
    if (data_.size() != 1)
        return std::unexpected(
            core::Failure{FailureCode::InvalidInput, core::Operation::PrepareModel});
    return data_.values().front().identity;
}

structure::PgReferenceCallInfo CaseReferenceSurface::preparation_work() const noexcept {
    return data_.size() == 1 ? data_.values().front().preparation
                             : structure::PgReferenceCallInfo{};
}

core::Result<CaseReferenceSurfaceSample>
CaseReferenceSurface::evaluate(const ExactRational& chord, const ExactRational& span,
                               const structure::PgReferenceLimits& limits,
                               core::ResourcePlan& plan) const noexcept {
    constexpr auto operation = core::Operation::EvaluateStep;
    return run(operation, [&]() -> core::Result<CaseReferenceSurfaceSample> {
        pg::ReferenceArithmetic m{limits, plan};
        m.cancel();
        if (data_.size() != 1)
            partition_fail(FailureCode::InvalidInput);
        const auto& data = data_.values().front();
        m.admit(pg::partition_add(fixed_live_rationals, pg::partition_mul(2, data.domain.size())));
        const PgPartitionPoint point{m.copy(chord), m.copy(span)};
        if (!contains(m, data.domain.values(), point))
            partition_fail(FailureCode::GeometryInvalid);
        auto result = evaluate_chart(point, m, true);
        m.cancel();
        result.call_info = m.info();
        return result;
    });
}

CaseReferenceSurfaceSample CaseReferenceSurface::evaluate_chart(const PgPartitionPoint& point,
                                                                pg::ReferenceArithmetic& m,
                                                                bool tangents) const {
    if (data_.size() != 1)
        partition_fail(FailureCode::InvalidInput);
    const auto& data = data_.values().front();
    const auto& [chord, span] = point;
    const std::array<ExactRational, 3> position{
        m.add(data.x0, m.add(m.mul(data.width, chord), m.mul(data.skew, span))), m.copy(span),
        pg::polynomial_derivative(m, data.height, {0, 0}, chord, span)};
    const auto world = [&](const auto& chart, bool translate) {
        auto result = pg::physical_rotate(m, data.rotation, chart);
        for (std::size_t axis = 0; axis < result.size(); ++axis) {
            result[axis] = m.mul(result[axis], data.length);
            if (translate)
                result[axis] = m.add(result[axis], data.origin[axis]);
        }
        return result;
    };
    CaseReferenceSurfaceSample result{world(position, true), {}, {}, {}};
    if (tangents) {
        result.chord_tangent =
            world(std::array<ExactRational, 3>{m.copy(data.width),
                                               {},
                                               pg::polynomial_derivative(m, data.height, {1, 0},
                                                                         chord, span)},
                  false);
        result.span_tangent =
            world(std::array<ExactRational, 3>{m.copy(data.skew), m.constant("1"),
                                               pg::polynomial_derivative(m, data.height, {0, 1},
                                                                         chord, span)},
                  false);
    }
    return result;
}

PgPartitionPoint CaseReferenceSurface::chart_projection(const std::array<ExactRational, 3>& world,
                                                        pg::ReferenceArithmetic& m) const {
    if (data_.size() != 1)
        partition_fail(FailureCode::InvalidInput);
    const auto& data = data_.values().front();
    std::array<ExactRational, 2> xy;
    for (std::size_t axis = 0; axis < world.size(); ++axis) {
        const auto translated = m.sub(world[axis], data.origin[axis]);
        for (std::size_t local = 0; local < xy.size(); ++local)
            xy[local] = m.add(xy[local], m.mul(data.rotation[axis * 3 + local], translated));
    }
    auto t = m.div(xy[1], data.length);
    auto s =
        m.div(m.sub(m.sub(m.div(xy[0], data.length), data.x0), m.mul(data.skew, t)), data.width);
    return {std::move(s), std::move(t)};
}

} // namespace uvlas::coupling::detail
