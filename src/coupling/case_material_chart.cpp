#include "case_material_chart.hpp"

#include "structure/pg_reference_polynomials.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <utility>

namespace uvlas::coupling::detail {
namespace pg = structure::detail;
namespace {
using core::FailureCode;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_mul;
using pg::partition_take;
using R = structure::ExactRational;
using Poly = structure::ExactChartPolynomial;
using Point = structure::PgPartitionPoint;

constexpr auto operation = core::Operation::PrepareModel;
constexpr std::size_t coordinate_degree = 4;
constexpr std::size_t determinant_degree = 6;

template <class F>
auto run(F action) noexcept -> std::invoke_result_t<F> {
    return core::guard_entry(operation, FailureCode::PrecisionUnresolved, {},
                             [&]() -> std::invoke_result_t<F> {
                                 try {
                                     return action();
                                 } catch (const pg::PartitionAbort& abort) {
                                     return std::unexpected(abort.failure);
                                 }
                             });
}

struct Range {
    core::Index first{}, past{};
    [[nodiscard]] bool contains(core::Index n) const noexcept {
        return first <= n && n < past;
    }
};
struct Piece {
    Range chord, span;
    std::array<Poly, 3> material;
    Poly determinant;
    std::size_t layer{};
    CaseMaterialSide side{};
};
struct Box {
    const R& xi0;
    const R& xi1;
    const R& eta0;
    const R& eta1;
};

constexpr std::uint64_t choose(std::size_t n, std::size_t k) noexcept {
    if (k > n)
        return 0;
    std::uint64_t result = 1;
    for (std::size_t i = 0; i < k; ++i)
        result = result * (n - i) / (i + 1);
    return result;
}

template <std::size_t Degree>
using Bernstein = std::array<std::array<R, Degree + 1>, Degree + 1>;

template <std::size_t Degree>
Bernstein<Degree> basis(pg::ReferenceArithmetic& m, const R& a, const R& b) {
    // Column i is the blossom evaluated at i copies of b and Degree-i copies of a.
    Bernstein<Degree> weight;
    std::array<R, Degree + 1> ap, bp;
    ap[0] = m.constant("1");
    bp[0] = m.constant("1");
    for (std::size_t p = 1; p <= Degree; ++p) {
        ap[p] = m.mul(ap[p - 1], a);
        bp[p] = m.mul(bp[p - 1], b);
    }
    for (std::size_t p = 0; p <= Degree; ++p) {
        const auto divisor = pg::index_value(m, partition_index(choose(Degree, p)));
        for (std::size_t i = 0; i <= Degree; ++i) {
            R sum;
            for (std::size_t k = 0; k <= p; ++k) {
                const auto count = choose(i, k) * choose(Degree - i, p - k);
                if (count != 0)
                    sum = m.add(sum, m.mul(pg::index_value(m, partition_index(count)),
                                           m.mul(ap[p - k], bp[k])));
            }
            weight[p][i] = m.div(sum, divisor);
        }
    }
    return weight;
}

template <std::size_t Degree>
Bernstein<Degree> coefficients(pg::ReferenceArithmetic& m, const Poly& poly,
                               const Bernstein<Degree>& x, const Bernstein<Degree>& y) {
    Bernstein<Degree> result;
    for (const auto& term : poly.terms()) {
        const auto a = static_cast<std::size_t>(term.powers.first);
        const auto b = static_cast<std::size_t>(term.powers.second);
        if (a > Degree || b > Degree)
            partition_fail(FailureCode::InvalidInput);
        for (std::size_t i = 0; i <= Degree; ++i)
            for (std::size_t j = 0; j <= Degree; ++j) {
                const auto& xi = x[a][i];
                const auto& eta = y[b][j];
                if (xi.sign() != 0 && eta.sign() != 0)
                    result[i][j] = m.add(result[i][j], m.mul(term.coefficient, m.mul(xi, eta)));
            }
    }
    return result;
}

R orientation(pg::ReferenceArithmetic& m, const Point& a, const Point& b, const R& s, const R& t) {
    return m.sub(m.mul(m.sub(b.s, a.s), m.sub(t, a.t)), m.mul(m.sub(b.t, a.t), m.sub(s, a.s)));
}

bool between(pg::ReferenceArithmetic& m, const R& x, const R& a, const R& b) {
    const auto first = m.compare(x, a), last = m.compare(x, b);
    return first == 0 || last == 0 || first != last;
}

struct Domain {
    std::span<const structure::PgPartitionSource> edges;
    std::span<const Point> points;

    bool contains(pg::ReferenceArithmetic& m, const R& s, const R& t, bool allow_boundary) const {
        bool inside = false;
        for (const auto& edge : edges) {
            m.visit(1);
            if (!std::holds_alternative<structure::PgPartitionDomainSource>(edge.descriptor))
                continue;
            const auto& a = points[static_cast<std::size_t>(edge.endpoints[0])];
            const auto& b = points[static_cast<std::size_t>(edge.endpoints[1])];
            const int orient = orientation(m, a, b, s, t).sign();
            if (orient == 0 && between(m, s, a.s, b.s) && between(m, t, a.t, b.t))
                return allow_boundary;
            const bool below_a = m.compare(a.t, t) <= 0;
            const bool below_b = m.compare(b.t, t) <= 0;
            if (below_a != below_b && (orient > 0) == below_a)
                inside = !inside;
        }
        return inside;
    }

    bool boundary_avoided(pg::ReferenceArithmetic& m, const Bernstein<coordinate_degree>& s,
                          const Bernstein<coordinate_degree>& t) const {
        // An open parameter box cannot meet an edge if its oriented edge function has
        // one sign there, or one coordinate stays strictly beyond both edge endpoints.
        for (const auto& edge : edges) {
            m.visit(1);
            if (!std::holds_alternative<structure::PgPartitionDomainSource>(edge.descriptor))
                continue;
            const auto& a = points[static_cast<std::size_t>(edge.endpoints[0])];
            const auto& b = points[static_cast<std::size_t>(edge.endpoints[1])];
            const auto dx = m.sub(b.s, a.s), dy = m.sub(b.t, a.t);
            const auto& x_min = dx.sign() >= 0 ? a.s : b.s;
            const auto& x_max = dx.sign() >= 0 ? b.s : a.s;
            const auto& y_min = dy.sign() >= 0 ? a.t : b.t;
            const auto& y_max = dy.sign() >= 0 ? b.t : a.t;
            bool positive = false, negative = false;
            bool left = true, right = true, down = true, up = true;
            for (std::size_t i = 0; i <= coordinate_degree; ++i)
                for (std::size_t j = 0; j <= coordinate_degree; ++j) {
                    const auto& x = s[i][j];
                    const auto& y = t[i][j];
                    const int side =
                        m.sub(m.mul(dx, m.sub(y, a.t)), m.mul(dy, m.sub(x, a.s))).sign();
                    positive = positive || side > 0;
                    negative = negative || side < 0;
                    left = left && m.compare(x, x_min) < 0;
                    right = right && m.compare(x, x_max) > 0;
                    down = down && m.compare(y, y_min) < 0;
                    up = up && m.compare(y, y_max) > 0;
                }
            if (positive == negative && !left && !right && !down && !up)
                return false;
        }
        return true;
    }
};

R middle(pg::ReferenceArithmetic& m, const R& a, const R& b) {
    return m.div(m.add(a, b), m.constant("2"));
}

bool layer_bound(pg::ReferenceArithmetic& m, const Bernstein<coordinate_degree>& zeta,
                 const R& bound, bool lower, const Box& box, const Box& panel) {
    // Nonnegative coefficients prove the closed-box bound. A positive coefficient
    // proves positivity on the open box, and similarly on every internal split edge.
    constexpr auto d = coordinate_degree;
    const std::array internal{box.xi0 != panel.xi0, box.xi1 != panel.xi1, box.eta0 != panel.eta0,
                              box.eta1 != panel.eta1};
    std::array<bool, 4> edge_positive{};
    bool any_positive = false, valid = true;
    for (std::size_t i = 0; i <= d; ++i)
        for (std::size_t j = 0; j <= d; ++j) {
            const auto cmp = lower ? m.compare(zeta[i][j], bound) : m.compare(bound, zeta[i][j]);
            const bool corner = (i == 0 || i == d) && (j == 0 || j == d);
            const bool interior_corner =
                corner && internal[i == 0 ? 0 : 1] && internal[j == 0 ? 2 : 3];
            if ((corner && cmp < 0) || (interior_corner && cmp == 0))
                partition_fail(FailureCode::GeometryInvalid);
            valid = valid && cmp >= 0;
            if (cmp <= 0)
                continue;
            any_positive = true;
            edge_positive[0] = edge_positive[0] || i == 0;
            edge_positive[1] = edge_positive[1] || i == d;
            edge_positive[2] = edge_positive[2] || j == 0;
            edge_positive[3] = edge_positive[3] || j == d;
        }
    for (std::size_t i = 0; i < internal.size(); ++i)
        if (internal[i] && !edge_positive[i])
            valid = false;
    return valid && any_positive;
}

R value(pg::ReferenceArithmetic& m, const Poly& polynomial, const R& xi, const R& eta) {
    return pg::polynomial_derivative(m, polynomial, {0, 0}, xi, eta);
}

bool certify_box(pg::ReferenceArithmetic& m, const Domain& domain, const Piece& piece,
                 const structure::PgPartitionOwnedLayer& layer, const Box& box, const Box& panel,
                 int& orientation_sign) {
    const auto xi = middle(m, box.xi0, box.xi1), eta = middle(m, box.eta0, box.eta1);
    const auto s = value(m, piece.material[0], xi, eta);
    const auto t = value(m, piece.material[1], xi, eta);
    const auto det = value(m, piece.determinant, xi, eta).sign();
    if (orientation_sign == 0)
        orientation_sign = det;
    if (det == 0 || orientation_sign != det || !domain.contains(m, s, t, false))
        partition_fail(FailureCode::GeometryInvalid);

    bool shape_known = true;
    {
        const auto bx = basis<determinant_degree>(m, box.xi0, box.xi1);
        const auto by = basis<determinant_degree>(m, box.eta0, box.eta1);
        const auto dets = coefficients<determinant_degree>(m, piece.determinant, bx, by);
        for (std::size_t i = 0; i <= determinant_degree; ++i)
            for (std::size_t j = 0; j <= determinant_degree; ++j)
                if (dets[i][j].sign() != orientation_sign) {
                    if ((i == 0 || i == determinant_degree) && (j == 0 || j == determinant_degree))
                        partition_fail(FailureCode::GeometryInvalid);
                    shape_known = false;
                }
    }
    const auto bx = basis<coordinate_degree>(m, box.xi0, box.xi1);
    const auto by = basis<coordinate_degree>(m, box.eta0, box.eta1);
    const auto zeta = coefficients<coordinate_degree>(m, piece.material[2], bx, by);
    if (piece.side == CaseMaterialSide::Interior) {
        const auto z = value(m, piece.material[2], xi, eta);
        if (m.compare(z, layer.lower) <= 0 || m.compare(z, layer.upper) >= 0)
            partition_fail(FailureCode::GeometryInvalid);
        const bool lower = layer_bound(m, zeta, layer.lower, true, box, panel);
        const bool upper = layer_bound(m, zeta, layer.upper, false, box, panel);
        shape_known = shape_known && lower && upper;
    } else {
        const auto& endpoint = piece.side == CaseMaterialSide::Lower ? layer.lower : layer.upper;
        for (const auto& row : zeta)
            for (const auto& z : row)
                if (m.compare(z, endpoint) != 0)
                    partition_fail(FailureCode::GeometryInvalid);
    }
    const auto bounded_s = coefficients<coordinate_degree>(m, piece.material[0], bx, by);
    const auto bounded_t = coefficients<coordinate_degree>(m, piece.material[1], bx, by);
    for (const auto i : {std::size_t{0}, coordinate_degree})
        for (const auto j : {std::size_t{0}, coordinate_degree})
            if (!domain.contains(m, bounded_s[i][j], bounded_t[i][j], true))
                partition_fail(FailureCode::GeometryInvalid);
    return domain.boundary_avoided(m, bounded_s, bounded_t) && shape_known;
}

struct PendingBox {
    R xi0, xi1, eta0, eta1;
    std::uint64_t depth{};
    [[nodiscard]] Box view() const noexcept {
        return {xi0, xi1, eta0, eta1};
    }
};

void certify_panel(pg::ReferenceArithmetic& m, const Domain& domain, const Piece& piece,
                   const structure::PgPartitionOwnedLayer& layer, const Box& panel,
                   std::uint64_t& subdivisions, int& orientation_sign) {
    if (certify_box(m, domain, piece, layer, panel, panel, orientation_sign))
        return;
    if (subdivisions == 0)
        partition_fail(FailureCode::PrecisionUnresolved);
    // At most one unvisited sibling per split stays live.
    auto pending = partition_take(
        m.plan.allocate<PendingBox>(partition_index(partition_add(subdivisions, 1)), operation));
    pending.values().front() = {m.copy(panel.xi0), m.copy(panel.xi1), m.copy(panel.eta0),
                                m.copy(panel.eta1), 0};
    std::size_t used = 1;
    while (used != 0) {
        auto current = std::move(pending.values()[--used]);
        if (current.depth != 0 &&
            certify_box(m, domain, piece, layer, current.view(), panel, orientation_sign))
            continue;
        if (subdivisions == 0)
            partition_fail(FailureCode::PrecisionUnresolved);
        --subdivisions;
        if (current.depth % 2 == 0) {
            auto mid = middle(m, current.xi0, current.xi1);
            pending.values()[used++] = {std::move(current.xi0), m.copy(mid), m.copy(current.eta0),
                                        m.copy(current.eta1), current.depth + 1};
            pending.values()[used++] = {std::move(mid), std::move(current.xi1),
                                        std::move(current.eta0), std::move(current.eta1),
                                        current.depth + 1};
        } else {
            auto mid = middle(m, current.eta0, current.eta1);
            pending.values()[used++] = {m.copy(current.xi0), m.copy(current.xi1),
                                        std::move(current.eta0), m.copy(mid), current.depth + 1};
            pending.values()[used++] = {std::move(current.xi0), std::move(current.xi1),
                                        std::move(mid), std::move(current.eta1), current.depth + 1};
        }
    }
}

template <class Field>
R exact(const io::SiCase& input, Field field, pg::ReferenceArithmetic& m,
        io::SchemaDimension expected = {{0, 0, 0, 0}}) {
    m.charge();
    const auto recipe = partition_take(input.quantity(field)).recipe;
    if (recipe.dimensions.exponents != expected.exponents || recipe.pi_power != 0)
        partition_fail(FailureCode::InvalidInput);
    const auto input_digits = std::visit(
        [](const auto& x) -> std::uint64_t {
            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(x)>, io::SiDecimalView>)
                return x.digits.size();
            else
                return partition_add(x.numerator.digits.size(), x.denominator.digits.size());
        },
        recipe.value);
    m.charge(input_digits);
    return partition_take(io::materialize_rational(recipe, m.leaf_limits, m.plan));
}

template <class Field>
Range cells(Field field) {
    return {partition_index(partition_take(field.at(0)).unsigned_value()),
            partition_index(partition_take(field.at(1)).unsigned_value())};
}

Piece read_piece(const io::SiCase& input, const io::schema::MaterialMapPiece& row,
                 std::span<const structure::PgPartitionOwnedLayer> layers,
                 pg::ReferenceArithmetic& m) {
    const auto layer = core::LayerId{row.layer.unsigned_value()};
    const auto found = std::ranges::find_if(layers, [&](const auto& x) {
        m.visit(1);
        return x.id == layer;
    });
    if (found == layers.end())
        partition_fail(FailureCode::GeometryInvalid);
    auto side = CaseMaterialSide::Interior;
    if (row.side.text() == "lower")
        side = CaseMaterialSide::Lower;
    else if (row.side.text() == "upper")
        side = CaseMaterialSide::Upper;
    Piece result{cells(row.chord_cells),
                 cells(row.span_cells),
                 {},
                 {},
                 static_cast<std::size_t>(found - layers.begin()),
                 side};

    constexpr auto max_terms = pg::reference_orders.size();
    const auto count = row.terms.size();
    if (count > max_terms)
        partition_fail(FailureCode::InvalidInput);
    std::array<structure::ExactPowers, max_terms> power{};
    std::array<std::array<R, max_terms>, 3> coefficient;
    for (std::size_t i = 0; i < count; ++i) {
        const auto term = partition_take(row.terms.at(partition_index(i))).record();
        const auto degree = cells(term.powers);
        power[i] = {degree.first, degree.past};
        for (std::size_t k = 0; k < coefficient.size(); ++k)
            coefficient[k][i] =
                exact(input, partition_take(term.coefficient.at(partition_index(k))), m);
    }
    for (std::size_t k = 0; k < result.material.size(); ++k) {
        const auto terms = [&]<std::size_t... I>(std::index_sequence<I...>) {
            return std::array{structure::ExactInputTerm{power[I], coefficient[k][I]}...};
        }(std::make_index_sequence<max_terms>{});
        m.charge(partition_mul(count, count));
        result.material[k] = partition_take(Poly::make(
            std::span{terms}.first(static_cast<std::size_t>(count)), m.leaf_limits, m.plan));
    }
    const auto derivative = [&](std::size_t axis, structure::ExactAxis direction) {
        m.charge(result.material[axis].terms().size());
        return partition_take(
            result.material[axis].differentiate(direction, 1, m.leaf_limits, m.plan));
    };
    const auto product = [&](const Poly& a, const Poly& b) {
        m.charge(partition_mul(a.terms().size(), b.terms().size()));
        return partition_take(a.multiply(b, m.leaf_limits, m.plan));
    };
    const auto ad = product(derivative(0, structure::ExactAxis::First),
                            derivative(1, structure::ExactAxis::Second));
    const auto bc = product(derivative(0, structure::ExactAxis::Second),
                            derivative(1, structure::ExactAxis::First));
    m.charge(partition_add(ad.terms().size(), bc.terms().size()));
    result.determinant = partition_take(ad.subtract(bc, m.leaf_limits, m.plan));
    return result;
}

const Piece& select(pg::ReferenceArithmetic& m, std::span<const Piece> pieces, core::Index chord,
                    core::Index span) {
    for (const auto& piece : pieces) {
        m.visit(1);
        if (piece.chord.contains(chord) && piece.span.contains(span))
            return piece;
    }
    partition_fail(FailureCode::InvalidInput);
}

R boundary_coefficient(pg::ReferenceArithmetic& m, const Poly& polynomial, bool fixed_chord,
                       const R& coordinate, core::Index varying_power) {
    R sum;
    for (const auto& term : polynomial.terms()) {
        m.charge();
        const auto varying = fixed_chord ? term.powers.second : term.powers.first;
        if (varying == varying_power) {
            const auto fixed = fixed_chord ? term.powers.first : term.powers.second;
            sum = m.add(sum, pg::times_monomial(m, m.copy(term.coefficient), fixed, 0, coordinate));
        }
    }
    return sum;
}

void interface(pg::ReferenceArithmetic& m, const Piece& a, const Piece& b, bool fixed_chord,
               const R& coordinate) {
    if (&a == &b)
        return;
    for (std::size_t k = 0; k < a.material.size(); ++k)
        for (core::Index d = 0; d <= static_cast<core::Index>(coordinate_degree); ++d)
            if (boundary_coefficient(m, a.material[k], fixed_chord, coordinate, d) !=
                boundary_coefficient(m, b.material[k], fixed_chord, coordinate, d))
                partition_fail(FailureCode::GeometryInvalid);
}

CaseMaterialChartPoint sample(pg::ReferenceArithmetic& m, const Piece& piece,
                              const structure::PgPartitionOwnedLayer& layer, R xi, R eta) {
    CaseMaterialChartPoint result;
    result.layer = layer.id;
    for (std::size_t k = 0; k < result.material.size(); ++k) {
        result.material[k] = value(m, piece.material[k], xi, eta);
        result.chord_tangent[k] = pg::polynomial_derivative(m, piece.material[k], {1, 0}, xi, eta);
        result.span_tangent[k] = pg::polynomial_derivative(m, piece.material[k], {0, 1}, xi, eta);
    }
    result.side = piece.side;
    if (result.side == CaseMaterialSide::Interior) {
        if (result.material[2] == layer.lower)
            result.side = CaseMaterialSide::Lower;
        else if (result.material[2] == layer.upper)
            result.side = CaseMaterialSide::Upper;
    }
    result.aerodynamic = {std::move(xi), std::move(eta)};
    return result;
}

R quarter(pg::ReferenceArithmetic& m, const R& a, const R& b, bool aft) {
    return m.add(a, m.div(m.mul(m.sub(b, a), m.constant(aft ? "3" : "1")), m.constant("4")));
}

io::schema::AeroSurface find_surface(const io::SchemaCase& source, std::uint64_t id,
                                     pg::ReferenceArithmetic& m) {
    const auto surfaces = partition_take(source.root()).model.record().aero_surfaces;
    for (std::uint64_t i = 0; i < surfaces.size(); ++i) {
        m.visit(1);
        const auto surface = partition_take(surfaces.at(partition_index(i))).record();
        if (surface.id.unsigned_value() == id)
            return surface;
    }
    partition_fail(FailureCode::InvalidInput);
}
} // namespace

struct CaseMaterialChart::Data {
    CaseMaterialChartIdentity identity{{{io::SourceId{0}, 0, {}}, 0, 0}, 0, 0, 0, false};
    std::optional<structure::PgReferenceGeometry> support;
    R geometry_error;
    core::PlannedBuffer<R> chord, span;
    core::PlannedBuffer<Piece> pieces;
};

CaseMaterialChart::CaseMaterialChart(core::PlannedBuffer<Data> data) noexcept
    : data_(std::move(data)) {}
CaseMaterialChart::CaseMaterialChart(CaseMaterialChart&&) noexcept = default;
CaseMaterialChart& CaseMaterialChart::operator=(CaseMaterialChart&&) noexcept = default;
CaseMaterialChart::~CaseMaterialChart() = default;

core::Result<CaseMaterialChart> CaseMaterialChart::prepare(const io::SiCase& source,
                                                           std::uint64_t id,
                                                           const CaseMaterialChartLimits& limits,
                                                           core::Workspace& workspace,
                                                           core::ResourcePlan& plan) noexcept {
    return run([&]() -> core::Result<CaseMaterialChart> {
        pg::ReferenceArithmetic m{limits.chart, plan};
        m.cancel();
        const auto& owner = partition_take(source.source_case()).get();
        const auto surface = find_surface(owner, id, m);
        const auto attachment = surface.attachment.record();
        const auto* structural = std::get_if<io::schema::AeroAttachmentStructural>(&attachment);
        if (!structural)
            partition_fail(FailureCode::InvalidInput);
        const auto map = structural->material_map;
        if (!map)
            partition_fail(FailureCode::InvalidInput);
        const auto pieces = map->record().pieces;
        const auto shell_id = structural->host.unsigned_value();
        auto support = partition_take(
            prepare_case_pg_reference(source, shell_id, limits.support, workspace, plan));
        const auto physical = partition_take(support.source());
        const auto layers = partition_take(physical.partition.layers());
        const Domain domain{partition_take(physical.partition.sources()),
                            partition_take(physical.partition.source_points())};

        constexpr std::uint64_t box_scratch = 512;
        constexpr std::uint64_t maximum_piece_terms = 3 * 15 + 28;
        const auto retained =
            partition_add(partition_add(surface.chord.size(), surface.span.size()),
                          partition_mul(pieces.size(), maximum_piece_terms));
        m.admit(partition_add(partition_add(box_scratch, retained),
                              partition_mul(4, limits.maximum_subdivisions)));
        auto data = partition_take(plan.allocate<Data>(1, operation));
        auto& result = data.values().front();
        result.identity = {partition_take(owner.identity()), id, shell_id,
                           surface.geometry.unsigned_value(),
                           surface.orientation.text() == "reversed"};
        const auto policy = partition_take(owner.root()).numerics.record().patch.record();
        const auto* active = std::get_if<io::schema::PatchPolicyActive>(&policy);
        if (!active)
            partition_fail(FailureCode::InvalidInput);
        result.geometry_error = exact(source, active->geometry_error, m, {{1, 0, 0, 0}});
        const auto knots = [&](auto field) {
            auto values =
                partition_take(plan.allocate<R>(partition_index(field.size()), operation));
            for (std::size_t i = 0; i < values.size(); ++i)
                values.values()[i] = exact(source, partition_take(field.at(partition_index(i))), m);
            return values;
        };
        result.chord = knots(surface.chord);
        result.span = knots(surface.span);
        result.pieces =
            partition_take(plan.allocate<Piece>(partition_index(pieces.size()), operation));
        for (std::size_t i = 0; i < result.pieces.size(); ++i) {
            m.visit(1);
            result.pieces.values()[i] = read_piece(
                source, partition_take(pieces.at(partition_index(i))).record(), layers, m);
        }
        const auto chord = result.chord.values(), span = result.span.values();
        const auto n_chord = chord.size() - 1, n_span = span.size() - 1;
        const auto piece_at = [&](std::size_t i, std::size_t j) -> const Piece& {
            return select(m, result.pieces.values(), partition_index(i), partition_index(j));
        };
        for (std::size_t i = 0; i < n_chord; ++i)
            for (std::size_t j = 0; j < n_span; ++j) {
                const auto& piece = piece_at(i, j);
                if (i != 0)
                    interface(m, piece_at(i - 1, j), piece, true, chord[i]);
                if (j != 0)
                    interface(m, piece_at(i, j - 1), piece, false, span[j]);
            }
        auto subdivisions = limits.maximum_subdivisions;
        int orientation_sign = 0;
        for (std::size_t i = 0; i < n_chord; ++i)
            for (std::size_t j = 0; j < n_span; ++j) {
                const auto& piece = piece_at(i, j);
                const Box panel{chord[i], chord[i + 1], span[j], span[j + 1]};
                certify_panel(m, domain, piece, layers[piece.layer], panel, subdivisions,
                              orientation_sign);
            }
        m.cancel();
        result.support.emplace(std::move(support));
        return CaseMaterialChart{std::move(data)};
    });
}

core::Result<CaseMaterialChartIdentity> CaseMaterialChart::identity() const noexcept {
    if (data_.size() != 1)
        return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
    return data_.values().front().identity;
}

core::Result<const R*> CaseMaterialChart::geometry_error_metres() const& noexcept {
    if (data_.size() != 1)
        return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
    return &data_.values().front().geometry_error;
}

core::Result<const structure::PgReferenceGeometry*> CaseMaterialChart::support() const& noexcept {
    const auto values = data_.values();
    if (values.size() == 1) {
        const auto& geometry = values.front().support;
        if (geometry)
            return &*geometry;
    }
    return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
}

core::Result<CaseMaterialChartPoint>
CaseMaterialChart::collocation(core::Index chord_cell, core::Index span_cell,
                               const structure::PgReferenceLimits& limits,
                               core::ResourcePlan& plan) const noexcept {
    return run([&]() -> core::Result<CaseMaterialChartPoint> {
        pg::ReferenceArithmetic m{limits, plan};
        const auto physical = partition_take(partition_take(support())->source());
        const auto layers = partition_take(physical.partition.layers());
        const auto& data = data_.values().front();
        const auto& piece = select(m, data.pieces.values(), chord_cell, span_cell);
        const auto i = static_cast<std::size_t>(chord_cell);
        const auto j = static_cast<std::size_t>(span_cell);
        const auto chord = data.chord.values(), span = data.span.values();
        m.admit();
        auto result =
            sample(m, piece, layers[piece.layer], quarter(m, chord[i], chord[i + 1], true),
                   middle(m, span[j], span[j + 1]));
        result.panel_extent = {m.sub(chord[i + 1], chord[i]), m.sub(span[j + 1], span[j])};
        result.call_info = m.info();
        return result;
    });
}

core::Result<CaseMaterialChartPoint>
CaseMaterialChart::vertex(core::Index chord_column, core::Index span_knot,
                          core::Index owning_span_cell, const structure::PgReferenceLimits& limits,
                          core::ResourcePlan& plan) const noexcept {
    return run([&]() -> core::Result<CaseMaterialChartPoint> {
        pg::ReferenceArithmetic m{limits, plan};
        m.cancel();
        const auto physical = partition_take(partition_take(support())->source());
        const auto layers = partition_take(physical.partition.layers());
        const auto& data = data_.values().front();
        const auto chord = data.chord.values(), span = data.span.values();
        const auto chord_last = partition_index(chord.size() - 1);
        const auto span_last = partition_index(span.size() - 1);
        if (chord_column < 0 || chord_column > chord_last || span_knot < 0 ||
            span_knot > span_last || owning_span_cell < 0 || owning_span_cell >= span_last ||
            (span_knot != owning_span_cell && span_knot != owning_span_cell + 1))
            partition_fail(FailureCode::InvalidInput);
        const auto& piece = select(m, data.pieces.values(), std::min(chord_column, chord_last - 1),
                                   owning_span_cell);
        const auto i = static_cast<std::size_t>(chord_column);
        m.admit();
        auto xi = chord_column == chord_last ? m.copy(chord[i])
                                             : quarter(m, chord[i], chord[i + 1], false);
        auto result = sample(m, piece, layers[piece.layer], std::move(xi),
                             m.copy(span[static_cast<std::size_t>(span_knot)]));
        result.call_info = m.info();
        return result;
    });
}

} // namespace uvlas::coupling::detail
