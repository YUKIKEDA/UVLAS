#include "../support/allocation_audit.hpp"
#include "pg_conditioning_test_support.hpp"
#include "pg_rank_storage.hpp"
#include "pg_reference_projection.hpp"
#include "pg_volume_position_arithmetic.hpp"

#include <bit>
#include <cmath>
#include <limits>
#include <uvlas/structure/pg_boundary_tangents.hpp>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <uvlas/structure/pg_volume_mass_reference.hpp>
#include <uvlas/structure/pg_volume_position.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;
PgReferenceBasis flat_basis(core::Workspace& ws, core::ResourcePlan& plan,
                            bool cyclic_frame = false) {
    Fixture f;
    for (auto& value : f.heights)
        value = f.r("2");
    const auto length = f.r("1");
    return make_basis(f, length, ws, plan, cyclic_frame);
}

Q moment(core::Index degree) {
    return degree % 2 == 0 ? Q{2} / (degree + 1) : Q{};
}

Q integrate_terms(const Term& a, const Term& b, const Q& h, const Q& measure) {
    // Independently integrate u, v and zeta; C contributes the factor 1/2.
    const auto i = static_cast<core::Index>(a.row.field);
    const auto j = static_cast<core::Index>(b.row.field);
    if (a.row.cell != b.row.cell || (i + j) % 2 != 0)
        return {};
    Q thickness_factor{1};
    for (core::Index k = 0; k < i + j; ++k)
        thickness_factor *= h;
    return a.value * b.value * moment(a.row.powers.first + b.row.powers.first) *
           moment(a.row.powers.second + b.row.powers.second) * moment(i + j) * measure *
           thickness_factor / ((i == 2 ? 2 : 1) * (j == 2 ? 2 : 1));
}

std::vector<Q> analytic_gram(const std::vector<std::vector<Term>>& basis, const Q& h = Q{1},
                             const Q& measure = Q("1/16")) {
    const auto count = basis.size();
    std::vector<Q> entries(slot(0, count));
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i)
            for (const auto& a : basis[i])
                for (const auto& b : basis[j])
                    entries[slot(i, j)] += integrate_terms(a, b, h, measure);
    return entries;
}

// The fixture cuts the chart by s+t <= 1+1/(4n). In local coordinates a
// boundary cell is a square minus one corner triangle or the opposite corner
// triangle alone. Integrate that triangle by successive exact integration:
// I(a,b) = (1/(a+1) - B(b+2,a+1))/(b+1).
Q fixture_corner_moment(core::Index a, core::Index b) {
    Q beta = Q{1} / (b + 2);
    for (core::Index k = 1; k <= a; ++k)
        beta *= Q{k} / (b + 2 + k);
    return (Q{1} / (a + 1) - beta) / (b + 1);
}

std::vector<Q> analytic_cut_gram(const std::vector<std::vector<Term>>& basis, core::Index n) {
    constexpr core::Index maximum_degree = 42;
    std::array<std::array<Q, maximum_degree + 1>, maximum_degree + 1> triangle{};
    for (core::Index a = 0; a <= maximum_degree; ++a)
        for (core::Index b = 0; b <= maximum_degree - a; ++b)
            triangle[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] =
                fixture_corner_moment(a, b);
    const auto count = basis.size();
    std::vector<Q> result(slot(0, count));
    const auto measure = Q{1} / (16 * n * n);
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i)
            for (const auto& a : basis[i])
                for (const auto& b : basis[j]) {
                    if (a.row.cell != b.row.cell)
                        continue;
                    const auto field_a = static_cast<core::Index>(a.row.field);
                    const auto field_b = static_cast<core::Index>(b.row.field);
                    const auto z = field_a + field_b;
                    if (z % 2 != 0)
                        continue;
                    const auto u = a.row.powers.first + b.row.powers.first;
                    const auto v = a.row.powers.second + b.row.powers.second;
                    const auto boundary = 4 * n - 1 - 2 * (a.row.cell.first + a.row.cell.second);
                    Q planar;
                    if (boundary >= 2) {
                        planar = moment(u) * moment(v);
                    } else {
                        if (boundary != -1 && boundary != 1)
                            throw std::runtime_error("Unexpected active cell in the cut oracle");
                        const auto& corner =
                            triangle[static_cast<std::size_t>(u)][static_cast<std::size_t>(v)];
                        if (boundary == 1)
                            planar = moment(u) * moment(v) - corner;
                        else
                            planar = (u + v) % 2 == 0 ? corner : -corner;
                    }
                    result[slot(i, j)] += a.value * b.value * planar * moment(z) * measure /
                                          ((field_a == 2 ? 2 : 1) * (field_b == 2 ? 2 : 1));
                }
    return result;
}

void make_sloping_cut(Fixture& f, const I& denominator) {
    f.domain.clear();
    f.domain.push_back(f.point("0", "0"));
    f.domain.push_back(f.point("1", "0"));
    f.domain.push_back({f.r("1"), rat(f, I{1}, denominator)});
    f.domain.push_back({rat(f, I{1}, denominator), f.r("1")});
    f.domain.push_back(f.point("0", "1"));
    for (auto& height : f.heights)
        height = f.r("2");
}

Q exact_congruence_row(const std::vector<Q>& gram, std::span<const ExactRational> transform,
                       std::size_t i, std::size_t count) {
    std::vector<Q> left(count);
    for (std::size_t k = 0; k <= i; ++k)
        for (std::size_t j = 0; j < count; ++j)
            left[j] += q(transform[slot(k, i)]) * gram[slot(std::min(k, j), std::max(k, j))];
    Q sum;
    for (std::size_t j = 0; j < count; ++j) {
        Q residual{i == j ? -1 : 0};
        for (std::size_t k = 0; k <= j; ++k)
            residual += left[k] * q(transform[slot(k, j)]);
        sum += residual < 0 ? -residual : residual;
    }
    return sum;
}

std::vector<PgGramEntry> enclose(Fixture& f, const std::vector<Q>& gram,
                                 std::uint32_t bits = dyadic_fractional_bits) {
    const I scale = I{1} << bits;
    std::vector<PgGramEntry> result;
    for (const auto& value : gram) {
        const Q scaled = value * scale;
        const I n = numerator(scaled), d = denominator(scaled);
        const I quotient = n / d;
        const bool fraction = n % d != 0;
        auto enclosure = PgGramEntry{rat(f, quotient - (fraction && n < 0 ? 1 : 0), scale),
                                     rat(f, quotient + (fraction && n > 0 ? 1 : 0), scale)};
        REQUIRE(q(enclosure.lower) <= value);
        REQUIRE(q(enclosure.upper) >= value);
        result.push_back(std::move(enclosure));
    }
    return result;
}

core::Result<PgConditioningCertificate> verify(PgReferenceBasis& basis,
                                               std::span<const PgGramEntry> gram,
                                               std::span<const ExactRational> candidate,
                                               const PgConditioningLimits& caps,
                                               core::Workspace& ws, core::ResourcePlan& plan) {
    return PgConditioningCertificate::certify(std::move(basis), {gram, dyadic_fractional_bits},
                                              {candidate, dyadic_fractional_bits}, caps, ws, plan);
}

struct AffineThicknessCase {
    core::Index n;
    Q s_min, s_max, t_min, t_max;
    bool certify;
};

ExactRational fixture_exact(Fixture& f, const Q& value) {
    return rat(f, numerator(value), denominator(value));
}

void set_affine_thickness_fixture(Fixture& f, const AffineThicknessCase& model) {
    f.domain.clear();
    for (const auto& s : {model.s_min, model.s_max}) {
        const auto forward = s == model.s_min;
        f.domain.push_back(
            {fixture_exact(f, s), fixture_exact(f, forward ? model.t_max : model.t_min)});
        f.domain.push_back(
            {fixture_exact(f, s), fixture_exact(f, forward ? model.t_min : model.t_max)});
    }
    f.thickness.clear();
    f.heights.clear();
    for (const Q& t : {Q{0}, Q{1}}) {
        for (const Q& s : {Q{0}, Q{"1/4"}, Q{1}}) {
            const Q h = s <= Q{"1/4"} ? Q{2 + 2 * s + t} : Q{1 + 6 * s + t};
            f.thickness.push_back({fixture_exact(f, s), fixture_exact(f, t)});
            f.heights.push_back(fixture_exact(f, h));
        }
    }
    f.triangles = {{0, 1, 4}, {0, 4, 3}, {1, 2, 5}, {1, 5, 4}};
    const auto width = f.r("4", "3"), skew = f.r("1", "9"), zero = f.r("0");
    f.graph = take(PgGraph::prepare({zero, width, skew, {}, model.n, PgGraphConvention::Version1},
                                    {exact_limits(), 64, {}}, f.input_plan));
}

constexpr std::array<std::size_t, 3> affine_oracle_degree{26, 42, 4};
using AffineOracleCell = std::array<std::vector<Q>, 3>;

std::array<Q, 48> rectangular_powers(const Q& low, const Q& high) {
    std::array<Q, 48> result{};
    Q a = low, b = high;
    for (std::size_t k = 0; k < result.size(); ++k) {
        result[k] = (b - a) / (k + 1);
        a *= low;
        b *= high;
    }
    return result;
}

Q choose(std::size_t n, std::size_t k) {
    Q result{1};
    for (std::size_t j = 1; j <= k; ++j)
        result *= Q{n + 1 - j} / j;
    return result;
}

void affine_rectangle_oracle(AffineOracleCell& result, const std::array<Q, 3>& h,
                             const std::array<Q, 48>& u, const std::array<Q, 48>& v) {
    for (std::size_t f = 0; f < result.size(); ++f) {
        const auto degree = affine_oracle_degree[f];
        auto& table = result[f];
        if (table.empty())
            table.resize((degree + 1) * (degree + 1));
        const auto power_of_h = 2 * f + 1;
        for (std::size_t i = 0; i <= power_of_h; ++i)
            for (std::size_t j = 0; j <= power_of_h - i; ++j) {
                const Q coefficient = choose(power_of_h, i) * choose(power_of_h - i, j) *
                                      power(h[0], static_cast<core::Index>(power_of_h - i - j)) *
                                      power(h[1], static_cast<core::Index>(i)) *
                                      power(h[2], static_cast<core::Index>(j));
                for (std::size_t a = 0; a <= degree; ++a)
                    for (std::size_t b = 0; b <= degree - a; ++b)
                        table[a * (degree + 1) + b] += coefficient * u[a + i] * v[b + j];
            }
    }
}

std::vector<AffineOracleCell> affine_oracle_tables(const AffineThicknessCase& model) {
    const auto axis = 2 * model.n, local = 4 * model.n;
    std::vector<AffineOracleCell> tables(static_cast<std::size_t>(axis * axis));
    const Q split{"1/4"};
    for (core::Index i = 0; i < axis; ++i)
        for (core::Index j = 0; j < axis; ++j) {
            auto& table = tables[static_cast<std::size_t>(i * axis + j)];
            const Q t_lo = std::max(model.t_min, Q{j} / axis);
            const Q t_hi = std::min(model.t_max, Q{j + 1} / axis);
            if (t_lo >= t_hi)
                continue;
            const auto v =
                rectangular_powers(local * t_lo - (2 * j + 1), local * t_hi - (2 * j + 1));
            for (const bool left : {true, false}) {
                const Q s_lo = std::max({model.s_min, Q{i} / axis, left ? Q{} : split});
                const Q s_hi = std::min({model.s_max, Q{i + 1} / axis, left ? split : Q{1}});
                if (s_lo >= s_hi)
                    continue;
                const auto u =
                    rectangular_powers(local * s_lo - (2 * i + 1), local * s_hi - (2 * i + 1));
                const auto slope = left ? 2 : 6, base = left ? 2 : 1;
                // Half thickness is H/(2L_ref) with L_ref = 3/2.
                const std::array<Q, 3> coefficients{
                    (Q{base} + Q{slope * (2 * i + 1) + 2 * j + 1} / local) / 3,
                    Q{slope} / (3 * local), Q{1} / (3 * local)};
                affine_rectangle_oracle(table, coefficients, u, v);
            }
        }
    return tables;
}

std::vector<Q> analytic_affine_gram(const std::vector<std::vector<Term>>& basis,
                                    const AffineThicknessCase& model) {
    const auto tables = affine_oracle_tables(model);
    const auto count = basis.size();
    std::vector<Q> result(slot(0, count));
    const auto axis = 2 * model.n;
    const Q measure = Q{"4/3"} / (16 * model.n * model.n);
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i)
            for (const auto& a : basis[i])
                for (const auto& b : basis[j]) {
                    if (a.row.cell != b.row.cell)
                        continue;
                    const auto field_a = static_cast<std::size_t>(a.row.field);
                    const auto field_b = static_cast<std::size_t>(b.row.field);
                    const auto z = field_a + field_b;
                    if (z % 2 != 0)
                        continue;
                    const auto index =
                        static_cast<std::size_t>(a.row.cell.first * axis + a.row.cell.second);
                    const auto& table = tables[index][z / 2];
                    const auto u =
                        static_cast<std::size_t>(a.row.powers.first + b.row.powers.first);
                    const auto v =
                        static_cast<std::size_t>(a.row.powers.second + b.row.powers.second);
                    const auto degree = affine_oracle_degree[z / 2];
                    const auto& planar = table.at(u * (degree + 1) + v);
                    result[slot(i, j)] += a.value * b.value * planar *
                                          moment(static_cast<core::Index>(z)) * measure /
                                          ((field_a == 2 ? 2 : 1) * (field_b == 2 ? 2 : 1));
                }
    return result;
}

std::array<Q, 4> conditioned_jet(const std::vector<std::vector<Term>>& source,
                                 std::span<const ExactRational> transform, std::size_t j,
                                 PgCellIndex cell, const Q& s, const Q& t, const Q& zeta) {
    const Q u = 4 * s - (2 * cell.first + 1);
    const Q v = 4 * t - (2 * cell.second + 1);
    std::array<Q, 4> out{};
    for (std::size_t i = 0; i <= j; ++i) {
        const auto& weight = q(transform[slot(i, j)]);
        if (weight == 0)
            continue;
        for (const auto& term : source[i]) {
            if (term.row.cell != cell)
                continue;
            const auto z = static_cast<core::Index>(term.row.field);
            const auto [a, b] = term.row.powers;
            const Q c = term.value * weight / (z == 2 ? 2 : 1);
            const auto h = power(zeta, z);
            out[0] += c * power(u, a) * power(v, b) * h;
            if (a > 0)
                out[1] += 4 * a * c * power(u, a - 1) * power(v, b) * h;
            if (b > 0)
                out[2] += 4 * b * c * power(u, a) * power(v, b - 1) * h;
            if (z > 0)
                out[3] += z * c * power(u, a) * power(v, b) * power(zeta, z - 1);
        }
    }
    return out;
}

void check_jet(const PgReferenceBasisValue& result, const std::array<Q, 4>& expected) {
    // This flat, unit-frame H=2 fixture has g=1 and physical J=I.
    REQUIRE(q(result.reference.g) == 1);
    const auto value = [](const PgReferenceScalar& x) -> Q { return q(x.rational) + q(x.radical); };
    CHECK(value(result.value) == expected[0]);
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(value(result.chart_derivatives[i]) == expected[i + 1]);
        CHECK(value(result.physical_gradient_per_metre[i]) == expected[i + 1]);
    }
}

Q flat_scalar(const PgReferenceScalar& x) {
    return q(x.rational) + q(x.radical);
}

std::array<Q, 3> flat_vector(const PgReferenceVector& vector) {
    return {flat_scalar(vector[0]), flat_scalar(vector[1]), flat_scalar(vector[2])};
}

void check_batch_limits(const PgConditioningCertificate& certificate, const PgReferencePoint& query,
                        const PgReferencePoint& invalid) {
    const auto one_point = std::span{&query, 1};
    const auto caps = query_limits();
    core::ResourcePlan control{profile};
    const auto reference = take(certificate.evaluate_reference(one_point, caps, control));
    const auto& used = reference.call_info;
    CHECK(used.allocated_owned_payload_bytes == control.allocated_bytes());
    auto exact = caps;
    exact.maximum_scalar_operations = used.scalar_operations;
    exact.maximum_partition_visits = used.partition_visits;
    core::ResourcePlan equal{used.allocated_owned_payload_bytes};
    const auto repeated = take(certificate.evaluate_reference(one_point, exact, equal));
    CHECK(repeated.call_info.allocated_owned_payload_bytes == equal.byte_limit());
    CHECK(repeated.points.values().front().identity == reference.points.values().front().identity);
    CHECK(flat_scalar(repeated.columns.values().back().value) ==
          flat_scalar(reference.columns.values().back().value));
    core::ResourcePlan short_bytes{used.allocated_owned_payload_bytes - 1};
    check_failure(certificate.evaluate_reference(one_point, exact, short_bytes),
                  core::FailureCode::ResourceExceeded);
    --exact.maximum_scalar_operations;
    check_failure(certificate.evaluate_reference(one_point, exact, control),
                  core::FailureCode::ResourceExceeded);
    ++exact.maximum_scalar_operations;
    --exact.maximum_partition_visits;
    check_failure(certificate.evaluate_reference(one_point, exact, control),
                  core::FailureCode::ResourceExceeded);
    const auto inline_bytes =
        sizeof(PgReferenceValue) + reference.column_count * sizeof(PgConditionedReferenceColumn);
    core::ResourcePlan short_inline{inline_bytes - 1};
    check_failure(certificate.evaluate_reference(one_point, caps, short_inline),
                  core::FailureCode::ResourceExceeded);
    CHECK(short_inline.allocated_bytes() == 0);
    core::ResourcePlan cancelled{0};
    exact.cancellation = requested_stop_token();
    check_failure(certificate.evaluate_reference(one_point, exact, cancelled),
                  core::FailureCode::Cancelled);
    check_failure(certificate.evaluate_reference({}, exact, cancelled),
                  core::FailureCode::Cancelled);
    CHECK(cancelled.allocated_bytes() == 0);
    auto short_external = caps;
    short_external.algebra.external_memory = ExactEstimateAdmission{0};
    core::ResourcePlan external{profile};
    check_failure(certificate.evaluate_reference(one_point, short_external, external),
                  core::FailureCode::ResourceExceeded);
    CHECK(external.allocated_bytes() == 0);
    const auto none = take(certificate.evaluate_reference({}, caps, cancelled));
    CHECK(none.column_count == reference.column_count);
    CHECK(none.points.size() == 0);
    CHECK(none.columns.size() == 0);
    CHECK(none.call_info.allocated_owned_payload_bytes == 0);
    CHECK(none.call_info.scalar_operations == 0);
    CHECK(none.call_info.partition_visits == 0);
    const std::array pair{query, invalid};
    check_failure(certificate.evaluate_reference(pair, caps, control),
                  core::FailureCode::InvalidInput);
    REQUIRE(certificate.evaluate_reference(one_point, caps, control));
}

void check_reference_batch(const PgConditioningCertificate& certificate,
                           const std::vector<std::vector<Term>>& terms,
                           std::span<const ExactRational> transform,
                           const PgIntegrationPartition& partition,
                           const std::array<PgReferencePoint, 5>& queries) {
    const auto caps = query_limits();
    const auto triangles = take(partition.triangles());
    const auto vertices = take(partition.vertices());
    const auto pieces = take(partition.pieces());
    const auto volumes = take(partition.volumes());
    const auto faces = take(partition.horizontal_faces());
    const Q s{"1/2"}, t{0}, zeta{"1/3"};
    const auto batch = [&] {
        core::ResourcePlan measured{profile};
        auto output = take(certificate.evaluate_reference(queries, caps, measured));
        CHECK(output.call_info.allocated_owned_payload_bytes == measured.allocated_bytes());
        return output;
    }();
    REQUIRE(batch.column_count == terms.size());
    REQUIRE(batch.points.size() == queries.size());
    REQUIRE(batch.columns.size() == queries.size() * terms.size());
    const auto points = batch.points.values();
    const auto lateral_face = take(partition.lateral_faces()).front();
    const auto edges = take(partition.edges());
    const auto& edge = edges[static_cast<std::size_t>(lateral_face.edge)];
    const auto first = lateral_face.canonical_direction ? 0U : 1U;
    const auto& edge_a = vertices[static_cast<std::size_t>(edge.vertices[first])];
    const auto& edge_b = vertices[static_cast<std::size_t>(edge.vertices[1 - first])];
    const Q edge_s = (q(edge_a.s) + q(edge_b.s)) / 2;
    const Q edge_t = (q(edge_a.t) + q(edge_b.t)) / 2;
    const auto inside_piece =
        (lateral_face.canonical_direction ? edge.left : edge.right).value_or(-1);
    REQUIRE(inside_piece >= 0);
    const auto edge_cell = pieces[static_cast<std::size_t>(inside_piece)].cell;
    const auto edge_layer = static_cast<std::size_t>(lateral_face.layer);
    const auto layers = take(partition.layers());
    const Q edge_zeta = (q(layers[edge_layer].lower) + q(layers[edge_layer].upper)) / 2;
    const Q interface_zeta{"-1/3"};
    const std::array<PgCellIndex, 5> cells{{{1, 0}, {0, 0}, edge_cell, {0, 0}, {0, 0}}};
    const std::array<std::array<Q, 3>, 5> positions{{{s, t, zeta},
                                                     {s, t, interface_zeta},
                                                     {edge_s, edge_t, edge_zeta},
                                                     {s, t, zeta},
                                                     {s, t, interface_zeta}}};
    const std::array<core::LayerId, 5> layer_ids{core::LayerId{29}, core::LayerId{11},
                                                 layers[edge_layer].id, core::LayerId{29},
                                                 core::LayerId{29}};
    const auto interface0 = std::get<PgReferenceHorizontalPoint>(queries[1]).face;
    const auto interface1 = std::get<PgReferenceHorizontalPoint>(queries[4]).face;
    const std::array<core::Index, 5> records{
        std::get<PgReferenceVolumePoint>(queries[0]).volume, interface0, 0,
        std::get<PgReferenceVolumePoint>(queries[3]).volume, interface1};
    const std::array<std::array<Q, 3>, 3> identity{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    std::uint64_t scalar_operations = 0;
    for (std::size_t p = 0; p < points.size(); ++p) {
        CAPTURE(p);
        const auto& point = points[p];
        CHECK(point.cell == cells[p]);
        CHECK(point.layer == layer_ids[p]);
        CHECK(point.partition_record == records[p]);
        CHECK(point.identity == points.front().identity);
        REQUIRE(q(point.g) == 1);
        CHECK(flat_vector(point.position_metres) == positions[p]);
        for (std::size_t axis = 0; axis < 3; ++axis)
            CHECK(flat_vector(point.jacobian_columns_metres[axis]) == identity[axis]);
        scalar_operations += point.call_info.scalar_operations;
        for (std::size_t i = 0; i < terms.size(); ++i) {
            CAPTURE(i);
            const auto& [x, y, z] = positions[p];
            const auto expected = conditioned_jet(terms, transform, i, cells[p], x, y, z);
            const auto& sample = batch.columns.values()[p * terms.size() + i];
            CHECK(flat_scalar(sample.value) == expected[0]);
            CHECK(flat_vector(sample.physical_gradient_per_metre) ==
                  std::array<Q, 3>{expected[1], expected[2], expected[3]});
        }
    }
    CHECK(batch.call_info.scalar_operations == scalar_operations);
    for (const auto p : {0U, 3U}) {
        REQUIRE(std::holds_alternative<PgReferenceVolumeMeasure>(points[p].measure));
        const auto& triangle = triangles[static_cast<std::size_t>(
            volumes[static_cast<std::size_t>(records[p])].triangle)];
        CHECK(flat_scalar(std::get<PgReferenceVolumeMeasure>(points[p].measure).jacobian) ==
              q(triangle.twice_area) * Q("4/3"));
    }
    const Q triangle_scale =
        q(triangles[static_cast<std::size_t>(faces[static_cast<std::size_t>(interface0)].triangle)]
              .twice_area);
    for (const auto p : {1U, 4U}) {
        REQUIRE(std::holds_alternative<PgReferenceHorizontalMeasure>(points[p].measure));
        const auto& face = std::get<PgReferenceHorizontalMeasure>(points[p].measure);
        const auto upper = p == 1;
        CHECK(face.side == (upper ? PgPartitionFaceSide::Upper : PgPartitionFaceSide::Lower));
        CHECK(flat_vector(face.area.outward) ==
              std::array<Q, 3>{Q{}, Q{}, upper ? triangle_scale : -triangle_scale});
        CHECK(flat_scalar(face.area.squared_norm) == triangle_scale * triangle_scale);
    }
    REQUIRE(std::holds_alternative<PgReferenceLateralMeasure>(points[2].measure));
    const auto& area = std::get<PgReferenceLateralMeasure>(points[2].measure).area;
    const Q depth = q(layers[edge_layer].upper) - q(layers[edge_layer].lower);
    const Q normal_x = depth * (q(edge_b.t) - q(edge_a.t));
    const Q normal_y = depth * (q(edge_a.s) - q(edge_b.s));
    CHECK(flat_vector(area.outward) == std::array<Q, 3>{normal_x, normal_y, Q{}});
    CHECK(flat_scalar(area.squared_norm) == normal_x * normal_x + normal_y * normal_y);
    const auto& lateral = std::get<PgReferenceLateralPoint>(queries[2]);
    check_batch_limits(certificate, queries.front(),
                       PgReferenceLateralPoint{-1, lateral.along, lateral.depth});
}

void check_reference_points(Fixture& f, const PgConditioningCertificate& certificate,
                            const std::vector<std::vector<Term>>& terms,
                            std::span<const ExactRational> transform) {
    const auto* basis = take(certificate.reference_basis());
    const auto& partition = take(basis->source()).partition;
    const auto triangles = take(partition.triangles());
    const auto vertices = take(partition.vertices());
    const auto pieces = take(partition.pieces());
    const auto volumes = take(partition.volumes());
    const auto faces = take(partition.horizontal_faces());
    const auto zero = f.r("0"), one = f.r("1"), half = f.r("1", "2");
    const Q s{"1/2"}, t{0}, zeta{"1/3"};
    struct VertexQuery {
        core::Index volume = -1, triangle = -1;
        std::size_t vertex = 0;
    };
    std::array<VertexQuery, 2> side{};
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        if (volumes[i].layer != 1)
            continue;
        const auto& triangle = triangles[static_cast<std::size_t>(volumes[i].triangle)];
        const auto cell = pieces[static_cast<std::size_t>(triangle.piece)].cell;
        if (cell.second != 0)
            continue;
        for (std::size_t k = 0; k < triangle.vertices.size(); ++k) {
            const auto& v = vertices[static_cast<std::size_t>(triangle.vertices[k])];
            if (q(v.s) == s && q(v.t) == t)
                side[static_cast<std::size_t>(cell.first)] = {static_cast<core::Index>(i),
                                                              volumes[i].triangle, k};
        }
    }
    REQUIRE(side[0].volume >= 0);
    REQUIRE(side[1].volume >= 0);
    std::size_t j = 0;
    for (; j < terms.size(); ++j) {
        const auto left = conditioned_jet(terms, transform, j, {0, 0}, s, t, zeta);
        const auto right = conditioned_jet(terms, transform, j, {1, 0}, s, t, zeta);
        const auto weights = transform.subspan(slot(0, j), j + 1);
        if (left[1] != right[1] &&
            std::ranges::count_if(weights, [](const auto& x) { return x.sign() != 0; }) > 1)
            break;
    }
    REQUIRE(j < terms.size());
    const PgRankBasisColumn column{static_cast<core::Index>(j)};
    auto caps = query_limits();
    std::array<Q, 2> values{};
    for (std::size_t i = 0; i < side.size(); ++i) {
        const auto& row = side[i];
        const PgReferenceVolumePoint query{row.volume, row.vertex == 1 ? one : zero,
                                           row.vertex == 2 ? one : zero, half};
        core::ResourcePlan evaluation{profile};
        const auto point = take(certificate.evaluate(query, column, caps, evaluation));
        const PgCellIndex cell{static_cast<core::Index>(i), 0};
        CHECK(point.reference.cell == cell);
        CHECK(point.reference.layer == core::LayerId{29});
        CHECK(q(point.reference.s) == s);
        CHECK(q(point.reference.t) == t);
        CHECK(q(point.reference.zeta) == zeta);
        CHECK(std::holds_alternative<PgReferenceVolumeMeasure>(point.reference.measure));
        const auto expected = conditioned_jet(terms, transform, j, cell, s, t, zeta);
        check_jet(point, expected);
        values[i] = expected[0];
        const auto& used = point.reference.call_info;
        CHECK(used.allocated_owned_payload_bytes == evaluation.allocated_bytes());
        auto exact = caps;
        exact.maximum_scalar_operations = used.scalar_operations;
        exact.maximum_partition_visits = used.partition_visits;
        core::ResourcePlan equal{used.allocated_owned_payload_bytes};
        REQUIRE(certificate.evaluate(query, column, exact, equal));
        core::ResourcePlan short_bytes{used.allocated_owned_payload_bytes - 1};
        check_failure(certificate.evaluate(query, column, exact, short_bytes),
                      core::FailureCode::ResourceExceeded);
        --exact.maximum_scalar_operations;
        check_failure(certificate.evaluate(query, column, exact, evaluation),
                      core::FailureCode::ResourceExceeded);
        ++exact.maximum_scalar_operations;
        --exact.maximum_partition_visits;
        check_failure(certificate.evaluate(query, column, exact, evaluation),
                      core::FailureCode::ResourceExceeded);
        core::ResourcePlan empty{0};
        check_failure(certificate.evaluate(query, {-1}, caps, empty),
                      core::FailureCode::InvalidInput);
        exact.cancellation = requested_stop_token();
        check_failure(certificate.evaluate(query, column, exact, empty),
                      core::FailureCode::Cancelled);
        CHECK(empty.allocated_bytes() == 0);
    }
    CHECK(values[0] == values[1]);

    const auto& left = side[0];
    std::array<core::Index, 2> interface{-1, -1};
    core::ResourcePlan faces_plan{profile};
    for (std::size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        if (face.triangle != left.triangle)
            continue;
        const auto direction =
            face.layer == 0 ? PgPartitionFaceSide::Upper : PgPartitionFaceSide::Lower;
        if (face.side != direction)
            continue;
        const PgReferenceHorizontalPoint query{static_cast<core::Index>(i),
                                               left.vertex == 1 ? one : zero,
                                               left.vertex == 2 ? one : zero};
        const auto point = take(certificate.evaluate(query, column, caps, faces_plan));
        CHECK(q(point.reference.zeta) == Q("-1/3"));
        CHECK(std::get<PgReferenceHorizontalMeasure>(point.reference.measure).side == direction);
        const auto layer = static_cast<std::size_t>(face.layer);
        CHECK(point.reference.layer == core::LayerId{layer == 0 ? 11U : 29U});
        interface[layer] = static_cast<core::Index>(i);
        check_jet(point, conditioned_jet(terms, transform, j, {0, 0}, s, t, Q("-1/3")));
    }
    REQUIRE(interface[0] >= 0);
    REQUIRE(interface[1] >= 0);

    const auto lateral = take(
        certificate.evaluate(PgReferenceLateralPoint{0, half, half}, column, caps, faces_plan));
    CHECK(std::holds_alternative<PgReferenceLateralMeasure>(lateral.reference.measure));
    const auto& r = lateral.reference;
    check_jet(lateral, conditioned_jet(terms, transform, j, r.cell, q(r.s), q(r.t), q(r.zeta)));
    check_failure(
        certificate.evaluate(PgReferenceLateralPoint{-1, half, half}, column, caps, faces_plan),
        core::FailureCode::InvalidInput);

    const auto volume_query = [&](const VertexQuery& row) {
        return PgReferenceVolumePoint{row.volume, row.vertex == 1 ? one : zero,
                                      row.vertex == 2 ? one : zero, half};
    };
    const auto interface_query = [&](core::Index face) {
        return PgReferenceHorizontalPoint{face, left.vertex == 1 ? one : zero,
                                          left.vertex == 2 ? one : zero};
    };
    const std::array<PgReferencePoint, 5> queries{
        volume_query(side[1]), interface_query(interface[0]),
        PgReferenceLateralPoint{0, half, half}, volume_query(side[0]),
        interface_query(interface[1])};
    check_reference_batch(certificate, terms, transform, partition, queries);
}

template <class Reference, class Apply>
void check_position_action(const Reference& reference, std::span<const Q> exact_shapes,
                           const Apply& apply, std::size_t rows_per_point = 1) {
    using namespace test_support;
    const auto points = rows_per_point * reference.points().size();
    const auto rank = reference.column_count();
    REQUIRE(exact_shapes.size() == points * rank);
    const auto n = 3 * rank;
    const auto m = 3 * points;
    const auto precision = reference.projection().representation;
    const bool paired = precision == core::Precision::TwoWord;
    const auto low_size = [paired](std::size_t size) { return paired ? size : 0; };
    const auto absolute_error = paired ? 0x1p-75 : 0x1p-32;
    const auto exact_error = independent_double(absolute_error);
    std::vector<double> q_hi(n), q_lo(low_size(n)), f_hi(m), f_lo(low_size(m));
    for (std::size_t i = 0; i < n; ++i)
        q_hi[i] = static_cast<double>(static_cast<int>(i % 7) - 3) / 4;
    for (std::size_t i = 0; i < m; ++i)
        f_hi[i] = static_cast<double>(i + 1) / 8;
    if (paired) {
        const auto active = std::ranges::find_if(exact_shapes, [](const Q& x) { return x != 0; });
        REQUIRE(active != exact_shapes.end());
        const auto column = static_cast<std::size_t>(active - exact_shapes.begin()) % rank;
        q_lo[3 * column] = 0x1p-50;
        f_lo[0] = -0x1p-51;
    }
    const PgCartesianWords q_input{precision, q_hi, q_lo}, f_input{precision, f_hi, f_lo};
    const auto forward = PgVolumePositionDirection::ToPoints;
    const auto transpose = PgVolumePositionDirection::TransposeToColumns;
    const auto forward_steps = q_hi.size() + q_lo.size() + m + m * rank;
    const auto transpose_steps = f_hi.size() + f_lo.size() + n + n * points;
    const PgVolumePositionLimits fwd_caps{forward_steps, absolute_error, {}};
    const PgVolumePositionLimits trans_caps{transpose_steps, absolute_error, {}};
    core::ResourcePlan plan{profile};
    auto output =
        take(PgVolumePositionResult::allocate(static_cast<core::Index>(n), precision, plan));
    const std::array<std::unique_ptr<core::Workspace>, 2> workspaces{
        take(core::Workspace::create(output.workspace_bytes(), plan)),
        take(core::Workspace::create(output.workspace_bytes() - 1, plan))};
    auto& workspace = *workspaces[0];
    CHECK(output.words().high.empty());

    const auto exact_word = [](PgCartesianWords words, std::size_t k) -> Q {
        auto value = independent_double(words.high[k]);
        if (!words.low.empty())
            value += independent_double(words.low[k]);
        return value;
    };
    const auto compare_run = [&](PgVolumePositionDirection direction, PgCartesianWords input,
                                 const PgVolumePositionLimits& caps) {
        AllocationAudit audit{cpp_backend, AuditThreadDomain::CallingThread};
        const auto applied = apply(reference, direction, input, caps, workspace, output);
        const auto allocations = audit.finish();
        REQUIRE(applied);
        CHECK(allocations.status == AllocationAuditStatus::Clean);
        CHECK(allocations.cpp_allocation_attempts == 0);
        return output.words();
    };
    const auto evaluate = [&](PgVolumePositionDirection direction, PgCartesianWords input,
                              const PgVolumePositionLimits& caps) {
        const auto words = compare_run(direction, input, caps);
        const bool to_points = direction == forward;
        const auto outer = to_points ? points : rank;
        const auto inner = to_points ? rank : points;
        REQUIRE(words.representation == precision);
        REQUIRE(words.high.size() == 3 * outer);
        REQUIRE(words.low.size() == low_size(3 * outer));
        std::vector<Q> retained(words.high.size());
        for (std::size_t i = 0; i < outer; ++i)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                CAPTURE(i, axis);
                Q expected;
                for (std::size_t j = 0; j < inner; ++j) {
                    const auto coefficient = to_points ? i * rank + j : j * rank + i;
                    expected += exact_shapes[coefficient] * exact_word(input, 3 * j + axis);
                }
                auto& value = retained[3 * i + axis];
                value = exact_word(words, 3 * i + axis);
                CHECK(abs(value - expected) <= exact_error);
            }
        return retained;
    };
    const auto positions = evaluate(forward, q_input, fwd_caps);
    const auto repeated = compare_run(forward, q_input, fwd_caps);
    for (std::size_t i = 0; i < positions.size(); ++i)
        CHECK(exact_word(repeated, i) == positions[i]);
    const auto forces = evaluate(transpose, f_input, trans_caps);

    Q forward_work{}, transpose_work{}, input_norm{};
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto scalar = exact_word(f_input, i);
        forward_work += positions[i] * scalar;
        input_norm += abs(scalar);
    }
    for (std::size_t i = 0; i < forces.size(); ++i) {
        const auto scalar = exact_word(q_input, i);
        transpose_work += forces[i] * scalar;
        input_norm += abs(scalar);
    }
    CHECK(abs(forward_work - transpose_work) <= exact_error * input_norm);

    auto caps = trans_caps;
    --caps.maximum_steps;
    const auto rejected = [&](PgCartesianWords input, const PgVolumePositionLimits& bounds,
                              core::Workspace& scratch, core::FailureCode failure) {
        check_failure(apply(reference, transpose, input, bounds, scratch, output), failure);
        REQUIRE(output.words().high.size() == forces.size());
        for (std::size_t i = 0; i < forces.size(); ++i)
            CHECK(exact_word(output.words(), i) == forces[i]);
    };
    rejected(f_input, caps, workspace, core::FailureCode::ResourceExceeded);
    rejected(f_input, trans_caps, *workspaces[1], core::FailureCode::ResourceExceeded);
    caps = trans_caps;
    caps.cancellation = requested_stop_token();
    rejected(f_input, caps, workspace, core::FailureCode::Cancelled);
    caps = trans_caps;
    caps.maximum_absolute_error = 0;
    rejected(f_input, caps, workspace, core::FailureCode::PrecisionUnresolved);
    auto bad_high = f_hi;
    bad_high.back() = std::numeric_limits<double>::infinity();
    rejected({precision, bad_high, f_lo}, trans_caps, workspace, core::FailureCode::Nonfinite);
    {
        const auto busy = take(workspace.acquire());
        static_cast<void>(busy);
        rejected(f_input, trans_caps, workspace, core::FailureCode::InvalidInput);
    }
    CHECK(evaluate(forward, q_input, fwd_caps) == positions);
    if (paired) {
        const PgCartesianWords mixed{core::Precision::Binary64, q_hi, {}};
        auto mixed_caps = fwd_caps;
        mixed_caps.maximum_steps -= q_lo.size();
        REQUIRE(evaluate(forward, mixed, mixed_caps) != positions);
    }
    // A single unit coefficient reproduces an already representable retained N.
    // The zero error limit still has to reject the nonzero source uncertainty.
    std::vector<double> unit(n);
    unit[1] = 1;
    const auto before = std::vector<double>(output.words().high.begin(), output.words().high.end());
    auto rigorous = fwd_caps;
    rigorous.maximum_absolute_error = 0;
    const PgCartesianWords unit_input{core::Precision::Binary64, unit, {}};
    AllocationAudit audit{cpp_backend, AuditThreadDomain::CallingThread};
    const auto refused = apply(reference, forward, unit_input, rigorous, workspace, output);
    const auto evidence = audit.finish();
    check_failure(refused, core::FailureCode::PrecisionUnresolved);
    CHECK(evidence.cpp_allocation_attempts == 0);
    CHECK(std::ranges::equal(output.words().high, before));
}

PgBoundaryTangentReference
check_tangent_preparation(PgBoundaryPositionReference&& position,
                          std::span<const PgBoundaryPositionPoint> queries,
                          std::span<const Q> derivatives) {
    // Keep the position and tangent results simultaneously usable. The position
    // and its certificate leave scope before the caller applies the tangent map.
    auto owner = std::move(position);
    const auto& certificate = owner.conditioning();
    const auto original = owner.projection();
    const PgBoundaryTangentProjection projection{original.representation,
                                                 original.maximum_shape_error};
    const auto caps = query_limits();
    const auto build = [&](std::span<const PgBoundaryPositionPoint> input,
                           const PgBoundaryTangentProjection& accuracy,
                           const PgReferenceLimits& limits, core::ResourcePlan& plan) {
        return PgBoundaryTangentReference::prepare(certificate, input, accuracy, limits, plan);
    };
    core::ResourcePlan none{0};
    auto stopped = caps;
    stopped.cancellation = requested_stop_token();
    check_failure(build(queries, projection, stopped, none), core::FailureCode::Cancelled);
    check_failure(build(queries, projection, caps, none), core::FailureCode::ResourceExceeded);
    const auto empty = take(build({}, projection, caps, none));
    REQUIRE(empty.points().empty());
    CHECK(empty.column_count() == owner.column_count());
    CHECK(none.allocated_bytes() == 0);
    {
        core::ResourcePlan plan{profile};
        auto strict = projection;
        strict.maximum_derivative_error = 0;
        check_failure(build(queries, strict, caps, plan), core::FailureCode::PrecisionUnresolved);
        auto limited = caps;
        limited.maximum_partition_visits = 1;
        check_failure(build(queries, projection, limited, plan),
                      core::FailureCode::ResourceExceeded);
        REQUIRE(certificate.info());
    }
    core::ResourcePlan plan{profile};
    auto result = take(build(queries, projection, caps, plan));
    CHECK(result.call_info().allocated_owned_payload_bytes == plan.allocated_bytes());
    REQUIRE(result.column_count() == owner.column_count());
    REQUIRE(result.points().size() == owner.points().size());
    for (std::size_t p = 0; p < result.points().size(); ++p) {
        CAPTURE(p);
        const auto& reference = owner.points()[p];
        const auto& derivative = result.points()[p];
        CHECK(derivative.identity == reference.identity);
        CHECK(derivative.cell == reference.cell);
        CHECK(derivative.layer == reference.layer);
        CHECK(derivative.partition_record == reference.partition_record);
        CHECK(derivative.side == reference.side);
    }

    const auto words = result.words();
    const bool paired = projection.representation == core::Precision::TwoWord;
    REQUIRE(words.high.size() == derivatives.size());
    REQUIRE(words.low.size() == (paired ? words.high.size() : 0));
    const auto tolerance = independent_double(projection.maximum_derivative_error);
    for (std::size_t i = 0; i < derivatives.size(); ++i) {
        const auto hi = independent_double(words.high[i]);
        const auto value = paired ? hi + independent_double(words.low[i]) : hi;
        CHECK(abs(value - derivatives[i]) <= tolerance);
    }

    const auto check_moved = [&] {
        auto output = take(PgVolumePositionResult::allocate(0, projection.representation, none));
        core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
        check_failure(apply_pg_boundary_tangents(result, PgVolumePositionDirection::ToPoints,
                                                 {projection.representation, {}, {}}, {0, 0, {}},
                                                 scratch, output),
                      core::FailureCode::InvalidInput);
    };
    auto moved = std::move(result);
    check_moved();
    return moved;
}

void check_boundary_position(Fixture& f, PgConditioningCertificate& certificate,
                             const std::vector<std::vector<Term>>& terms,
                             std::span<const ExactRational> transform,
                             PgBoundaryPositionProjection projection) {
    const auto* basis = take(certificate.reference_basis());
    const auto& partition = take(basis->source()).partition;
    const auto triangles = take(partition.triangles());
    const auto vertices = take(partition.vertices());
    const auto pieces = take(partition.pieces());
    const auto horizontal = take(partition.horizontal_faces());
    const auto layers = take(partition.layers());
    const auto upper = std::ranges::find_if(horizontal, [](const auto& face) {
        return face.layer == 0 && face.side == PgPartitionFaceSide::Upper;
    });
    REQUIRE(upper != horizontal.end());
    const auto lower = std::ranges::find_if(horizontal, [&](const auto& face) {
        return face.layer == 1 && face.triangle == upper->triangle &&
               face.side == PgPartitionFaceSide::Lower;
    });
    REQUIRE(lower != horizontal.end());
    const auto upper_id = static_cast<core::Index>(upper - horizontal.begin());
    const auto lower_id = static_cast<core::Index>(lower - horizontal.begin());
    const auto& triangle = triangles[static_cast<std::size_t>(upper->triangle)];
    const auto face_cell = pieces[static_cast<std::size_t>(triangle.piece)].cell;
    const auto face_coordinate = [&](auto member) -> Q {
        const auto a = q(vertices[static_cast<std::size_t>(triangle.vertices[0])].*member);
        const auto b = q(vertices[static_cast<std::size_t>(triangle.vertices[1])].*member);
        const auto c = q(vertices[static_cast<std::size_t>(triangle.vertices[2])].*member);
        return (2 * a + b + c) / 4;
    };
    const Q face_zeta = q(layers[0].upper);
    REQUIRE(face_zeta == q(layers[1].lower));

    // Choose a record whose CCW direction reverses its stored edge.
    const auto lateral = take(partition.lateral_faces());
    const auto reverse =
        std::ranges::find_if(lateral, [](const auto& face) { return !face.canonical_direction; });
    REQUIRE(reverse != lateral.end());
    const auto lateral_id = static_cast<core::Index>(reverse - lateral.begin());
    const auto edge = take(partition.edges())[static_cast<std::size_t>(reverse->edge)];
    const auto piece = edge.right.value_or(-1);
    REQUIRE(piece >= 0);
    const auto edge_cell = pieces[static_cast<std::size_t>(piece)].cell;
    const auto a = static_cast<std::size_t>(edge.vertices[0]);
    const auto b = static_cast<std::size_t>(edge.vertices[1]);
    const auto lateral_layer = static_cast<std::size_t>(reverse->layer);
    const auto& layer = layers[lateral_layer];
    const std::array<Q, 3> edge_position{(q(vertices[a].s) + q(vertices[b].s)) / 2,
                                         (q(vertices[a].t) + q(vertices[b].t)) / 2,
                                         (q(layer.lower) + q(layer.upper)) / 2};
    using Vertex = std::remove_cvref_t<decltype(vertices.front())>;
    const std::array<Q, 3> face_position{face_coordinate(&Vertex::s), face_coordinate(&Vertex::t),
                                         face_zeta};
    const auto& root = vertices[static_cast<std::size_t>(triangle.vertices[0])];
    const auto& e1 = vertices[static_cast<std::size_t>(triangle.vertices[1])];
    const auto& e2 = vertices[static_cast<std::size_t>(triangle.vertices[2])];
    const std::array<std::array<Q, 3>, 2> triangle_steps{
        {{q(e1.s) - q(root.s), q(e1.t) - q(root.t), Q{}},
         {q(e2.s) - q(root.s), q(e2.t) - q(root.t), Q{}}}};
    const std::array<std::array<Q, 3>, 2> lateral_steps{
        {{q(vertices[a].s) - q(vertices[b].s), q(vertices[a].t) - q(vertices[b].t), Q{}},
         {Q{}, Q{}, q(layer.upper) - q(layer.lower)}}};
    const auto quarter = f.r("1", "4"), half = f.r("1", "2");
    const PgReferenceHorizontalPoint lower_query{lower_id, quarter, quarter};
    const std::array<PgBoundaryPositionPoint, 4> queries{
        lower_query, PgReferenceLateralPoint{lateral_id, half, half},
        PgReferenceHorizontalPoint{upper_id, quarter, quarter}, lower_query};
    const std::array<PgCellIndex, 4> cells{face_cell, edge_cell, face_cell, face_cell};
    const std::array<core::Index, 4> records{lower_id, lateral_id, upper_id, lower_id};
    const std::array<core::LayerId, 4> expected_layers{layers[1].id, layer.id, layers[0].id,
                                                       layers[1].id};
    const std::array<PgBoundaryPositionSide, 4> sides{
        PgBoundaryPositionSide::Lower, PgBoundaryPositionSide::Lateral,
        PgBoundaryPositionSide::Upper, PgBoundaryPositionSide::Lower};
    const auto caps = query_limits();
    const auto prepare = [&](std::span<const PgBoundaryPositionPoint> points,
                             const PgBoundaryPositionProjection& requested,
                             const PgReferenceLimits& limits, core::ResourcePlan& plan) {
        return PgBoundaryPositionReference::prepare(std::move(certificate), points, requested,
                                                    limits, plan);
    };
    {
        core::ResourcePlan none{0};
        auto stopped = caps;
        stopped.cancellation = requested_stop_token();
        check_failure(prepare(queries, projection, stopped, none), core::FailureCode::Cancelled);
        CHECK(none.allocated_bytes() == 0);
        REQUIRE(certificate.info());
    }
    {
        core::ResourcePlan plan{profile};
        const std::array<PgBoundaryPositionPoint, 2> invalid{
            lower_query, PgReferenceLateralPoint{-1, half, half}};
        check_failure(prepare(invalid, projection, caps, plan), core::FailureCode::InvalidInput);
        REQUIRE(certificate.info());
        auto strict = projection;
        strict.maximum_shape_error = 0;
        check_failure(prepare(queries, strict, caps, plan), core::FailureCode::PrecisionUnresolved);
        REQUIRE(certificate.info());
    }
    auto result = [&] {
        core::ResourcePlan plan{profile};
        auto reference = take(prepare(queries, projection, caps, plan));
        CHECK(reference.call_info().allocated_owned_payload_bytes == plan.allocated_bytes());
        return reference;
    }();
    CHECK_FALSE(certificate.info());
    REQUIRE(result.conditioning().info());
    REQUIRE(result.column_count() == terms.size());
    REQUIRE(result.points().size() == queries.size());
    const auto words = result.words();
    REQUIRE(words.high.size() == terms.size() * queries.size());
    const bool paired = projection.representation == core::Precision::TwoWord;
    REQUIRE(words.low.size() == (paired ? words.high.size() : 0));
    const auto tolerance = independent_double(projection.maximum_shape_error);
    std::vector<Q> exact_shapes, exact_derivatives;
    for (std::size_t p = 0; p < queries.size(); ++p) {
        CAPTURE(p);
        const auto& row = result.points()[p];
        CHECK(row.identity == result.points().front().identity);
        CHECK(row.cell == cells[p]);
        CHECK(row.partition_record == records[p]);
        CHECK(row.layer == expected_layers[p]);
        CHECK(row.side == sides[p]);
        const auto& [s, t, z] = p == 1 ? edge_position : face_position;
        const auto& directions = p == 1 ? lateral_steps : triangle_steps;
        const auto first = 2 * p * terms.size();
        exact_derivatives.resize(first + 2 * terms.size());
        for (std::size_t j = 0; j < terms.size(); ++j) {
            const auto jet = conditioned_jet(terms, transform, j, cells[p], s, t, z);
            const auto& expected = jet[0];
            exact_shapes.push_back(expected);
            for (std::size_t d = 0; d < directions.size(); ++d) {
                const auto& step = directions[d];
                exact_derivatives[first + d * terms.size() + j] =
                    jet[1] * step[0] + jet[2] * step[1] + jet[3] * step[2];
            }
            const auto k = p * terms.size() + j;
            const auto hi = independent_double(words.high[k]);
            const auto value = paired ? hi + independent_double(words.low[k]) : hi;
            CHECK(abs(value - expected) <= tolerance);
        }
    }
    check_position_action(result, exact_shapes, apply_pg_boundary_position);
    const auto tangent = check_tangent_preparation(std::move(result), queries, exact_derivatives);
    check_position_action(tangent, exact_derivatives, apply_pg_boundary_tangents, 2);
}

void check_volume_mass_projection(Fixture& f, PgConditioningCertificate& certificate,
                                  const std::vector<std::vector<Term>>& terms,
                                  std::span<const ExactRational> transform,
                                  PgVolumeMassProjection projection, bool negative_controls) {
    const auto* basis = take(certificate.reference_basis());
    const auto& partition = take(basis->source()).partition;
    const auto volumes = take(partition.volumes());
    const auto vertices = take(partition.vertices());
    const auto pieces = take(partition.pieces());
    const auto triangles = take(partition.triangles());
    const auto layers = take(partition.layers());
    const auto first = volumes.front();
    const auto same_triangle = std::ranges::find_if(volumes, [&](const auto& x) {
        return x.triangle == first.triangle && x.layer != first.layer;
    });
    REQUIRE(same_triangle != volumes.end());
    const auto second_index = static_cast<core::Index>(same_triangle - volumes.begin());
    const auto& triangle = triangles[static_cast<std::size_t>(first.triangle)];
    const auto cell = pieces[static_cast<std::size_t>(triangle.piece)].cell;
    const auto third = f.r("1", "3");
    const std::array queries{PgReferenceVolumePoint{second_index, third, third, third},
                             PgReferenceVolumePoint{0, third, third, third}};
    const auto request_limits = query_limits();
    const auto prepare = [&](const PgVolumeMassProjection& policy, const PgReferenceLimits& caps,
                             core::ResourcePlan& plan,
                             std::span<const PgReferenceVolumePoint> points) {
        return PgVolumeMassReference::prepare(std::move(certificate), points, policy, caps, plan);
    };
    if (negative_controls) {
        core::ResourcePlan none{0};
        auto invalid = projection;
        invalid.representation = static_cast<core::Precision>(255);
        check_failure(prepare(invalid, request_limits, none, queries),
                      core::FailureCode::PrecisionUnsupported);
        invalid = projection;
        invalid.maximum_shape_error = -0x1p-1074;
        check_failure(prepare(invalid, request_limits, none, queries),
                      core::FailureCode::InvalidInput);
        invalid.maximum_shape_error = std::numeric_limits<double>::infinity();
        check_failure(prepare(invalid, request_limits, none, queries),
                      core::FailureCode::Nonfinite);
        auto stopped = request_limits;
        stopped.cancellation = requested_stop_token();
        check_failure(prepare(projection, stopped, none, queries), core::FailureCode::Cancelled);
        check_failure(prepare(projection, request_limits, none, queries),
                      core::FailureCode::ResourceExceeded);
        CHECK(none.allocated_bytes() == 0);
        auto zero_work = request_limits;
        zero_work.maximum_partition_visits = 0;
        core::ResourcePlan work_reject{profile};
        check_failure(prepare(projection, zero_work, work_reject, queries),
                      core::FailureCode::ResourceExceeded);
        CHECK(work_reject.allocated_bytes() == 0);
        auto strict = projection;
        strict.maximum_volume_jacobian_error_metres_cubed = 0;
        core::ResourcePlan reject{profile};
        check_failure(prepare(strict, request_limits, reject, queries),
                      core::FailureCode::PrecisionUnresolved);
        const std::array invalid_query{queries.front(),
                                       PgReferenceVolumePoint{-1, third, third, third}};
        check_failure(prepare(projection, request_limits, reject, invalid_query),
                      core::FailureCode::InvalidInput);
        REQUIRE(certificate.info());
    }
    // Calculate all oracle coordinates before the certificate owner is transferred.
    Q s{}, t{};
    for (const auto vertex : triangle.vertices) {
        const auto& value = vertices[static_cast<std::size_t>(vertex)];
        s += q(value.s) / 3;
        t += q(value.t) / 3;
    }
    const std::array selected_layers{static_cast<std::size_t>(same_triangle->layer),
                                     static_cast<std::size_t>(first.layer)};
    const Q chart_area = q(triangle.twice_area);
    std::array<Q, 2> expected_jacobian, expected_zeta;
    std::array<core::LayerId, 2> expected_layer_id{core::LayerId{0}, core::LayerId{0}};
    for (std::size_t row = 0; row < queries.size(); ++row) {
        const auto& layer = layers[selected_layers[row]];
        const Q thickness = q(layer.upper) - q(layer.lower);
        expected_jacobian[row] = chart_area * thickness;
        expected_zeta[row] = q(layer.lower) + thickness / 3;
        expected_layer_id[row] = layer.id;
    }
    const auto result = [&] {
        core::ResourcePlan evaluation{profile};
        auto retained = take(prepare(projection, request_limits, evaluation, queries));
        CHECK(retained.call_info().allocated_owned_payload_bytes == evaluation.allocated_bytes());
        return retained;
    }();
    CHECK_FALSE(certificate.info());
    CHECK(take(result.conditioning().info()).rank == terms.size());
    REQUIRE(result.column_count() == terms.size());
    REQUIRE(result.points().size() == queries.size());
    CHECK(result.projection().representation == projection.representation);
    const auto words = result.words();
    const bool two_word = projection.representation == core::Precision::TwoWord;
    const auto stride = terms.size() + 1;
    REQUIRE(words.high.size() == queries.size() * stride);
    REQUIRE(words.low.size() == (two_word ? words.high.size() : 0));
    const auto stored = [&](std::size_t i) -> Q {
        const auto high = independent_double(words.high[i]);
        return two_word ? high + independent_double(words.low[i]) : high;
    };
    const auto budget_n = independent_double(projection.maximum_shape_error);
    const auto budget_j = independent_double(projection.maximum_volume_jacobian_error_metres_cubed);
    std::vector<Q> exact_shapes;
    exact_shapes.reserve(terms.size() * queries.size());
    for (std::size_t row = 0; row < queries.size(); ++row) {
        const auto& point = result.points()[row];
        CHECK(point.cell == cell);
        CHECK(point.partition_record == queries[row].volume);
        CHECK(point.layer == expected_layer_id[row]);
        CHECK(point.identity == result.points().front().identity);
        CHECK(abs(stored(row * stride) - expected_jacobian[row]) <= budget_j);
        for (std::size_t column = 0; column < terms.size(); ++column) {
            CAPTURE(row, column);
            const auto expected =
                conditioned_jet(terms, transform, column, cell, s, t, expected_zeta[row])[0];
            CHECK(abs(stored(row * stride + column + 1) - expected) <= budget_n);
            exact_shapes.push_back(expected);
        }
    }
    if (two_word) {
        CHECK(std::ranges::any_of(words.low, [](double low) { return low != 0; }));
        CHECK(abs(independent_double(words.high.front()) - expected_jacobian.front()) > budget_j);
    }
    check_position_action(result, exact_shapes, apply_pg_volume_position);
}

} // namespace

TEST_CASE("PG full flat Gram certifies the transform and preserves all mandatory fields",
          "[SH-D02][SH-D03][SH-D05][N4-02][N4-03][N4-06][N4-08][N4-09][verification]") {
    Fixture f;
    core::Workspace ws{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan plan{profile};
    auto basis = flat_basis(ws, plan);
    const auto* rank = take(basis.rank());
    const auto count = static_cast<std::size_t>(take(basis.size()));
    const auto protected_count =
        static_cast<std::size_t>(take(rank->preparation()).independent_required);
    REQUIRE(count == 109);
    REQUIRE(protected_count > 0);
    const auto terms = basis_terms(*rank, count, plan);
    const auto gram = analytic_gram(terms);
    const auto enclosure = enclose(f, gram);
    const auto transform = proposal(f, gram, count);
    const auto limits = query_limits();
    std::vector<Q> original_required(33 * count);
    for (std::size_t t = 0; t < 33; ++t)
        for (std::size_t i = 0; i < count; ++i)
            original_required[t * count + i] =
                q(take(rank->required({static_cast<core::Index>(i)}, {static_cast<core::Index>(t)},
                                      limits.algebra, plan)));
    Q last{};
    const auto canonical_rows = take(rank->rows());
    core::Index checked_row = 0;
    for (std::size_t row = 0; row < canonical_rows.size(); ++row) {
        checked_row = static_cast<core::Index>(row);
        for (std::size_t i = 0; i < count; ++i)
            last += q(take(rank->basis({static_cast<core::Index>(i)}, checked_row, limits.algebra,
                                       plan))) *
                    q(transform[slot(i, count - 1)]);
        if (last != 0)
            break;
    }
    REQUIRE(last != 0);
    const PgFlatGramLimits flat_caps{query_limits(), steps};
    auto supplied =
        take(PgFlatGram::prepare(std::move(basis), dyadic_fractional_bits, flat_caps, ws, plan));
    const auto owned = take(supplied.enclosure());
    REQUIRE(owned.upper_triangle.size() == enclosure.size());
    CHECK(owned.fractional_bits == dyadic_fractional_bits);
    for (std::size_t i = 0; i < enclosure.size(); ++i) {
        CHECK(owned.upper_triangle[i].lower == enclosure[i].lower);
        CHECK(owned.upper_triangle[i].upper == enclosure[i].upper);
    }
    const auto certify_owned = [&](PgFlatGram& source, core::ResourcePlan& budget) {
        return std::move(source).certify({transform, dyadic_fractional_bits}, {profile, steps, {}},
                                         ws, budget);
    };
    core::ResourcePlan refused{0};
    check_failure(certify_owned(supplied, refused), core::FailureCode::ResourceExceeded);
    REQUIRE(supplied.enclosure());
    const auto prior = plan.allocated_bytes();
    auto certificate = take(certify_owned(supplied, plan));
    const auto info = take(certificate.info());
    CHECK(info.rank == count);
    CHECK(info.payload_bytes == plan.allocated_bytes() - prior);
    CHECK(take(take(certificate.reference_basis())->size()) == count);
    CHECK(take(certificate.transform_coefficient({1}, {0}, limits.algebra, plan)).sign() == 0);
    CHECK(q(take(certificate.transform_coefficient({0}, {0}, limits.algebra, plan))) ==
          q(transform[slot(0, 0)]));
    CHECK(info.protected_rank == protected_count);
    const Q bound = q(take(certificate.congruence_bound(limits.algebra, plan)));
    CHECK(bound > 0);
    CHECK(bound <= Q("1/100000000000000000000"));
    for (const auto row : std::array<std::size_t, 3>{0, count / 2, count - 1})
        CHECK(exact_congruence_row(gram, transform, row, count) <= bound);
    CHECK(info.steps > 0);
    check_reference_points(f, certificate, terms, transform);
    CHECK(q(take(certificate.coefficient({static_cast<core::Index>(count - 1)}, checked_row, limits,
                                         plan))) == last);
    for (std::size_t t = 0; t < 33; ++t) {
        const auto coords =
            take(certificate.required_coordinates({static_cast<core::Index>(t)}, limits, ws, plan));
        REQUIRE(coords.size() == count);
        for (std::size_t i = 0; i < count; ++i) {
            Q reconstructed{};
            for (std::size_t j = i; j < protected_count; ++j)
                reconstructed += q(transform[slot(i, j)]) * q(coords.values()[j]);
            if (i >= protected_count)
                CHECK(coords.values()[i].sign() == 0);
            CHECK(reconstructed == original_required[t * count + i]);
        }
    }
    check_volume_mass_projection(f, certificate, terms, transform,
                                 {core::Precision::TwoWord, 0x1p-90, 0x1p-100}, true);
    auto second = flat_basis(ws, plan);
    auto binary = take(verify(second, enclosure, transform, {profile, steps, {}}, ws, plan));
    check_volume_mass_projection(f, binary, terms, transform,
                                 {core::Precision::Binary64, 0x1p-43, 0x1p-53}, false);
    for (const auto boundary_projection : std::array<PgBoundaryPositionProjection, 2>{
             {{core::Precision::Binary64, 0x1p-43}, {core::Precision::TwoWord, 0x1p-90}}}) {
        const auto cyclic_frame = boundary_projection.representation == core::Precision::TwoWord;
        auto owned_basis = flat_basis(ws, plan, cyclic_frame);
        auto boundary =
            take(verify(owned_basis, enclosure, transform, {profile, steps, {}}, ws, plan));
        check_boundary_position(f, boundary, terms, transform, boundary_projection);
    }
}

TEST_CASE("PG conditioning rejects row-sum excess and honors exact work and byte limits",
          "[SH-D03][N4-02][N4-08][negative][verification]") {
    Fixture f;
    core::Workspace ws{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan home{profile}, reject{profile};
    auto basis = flat_basis(ws, home);
    const auto count = static_cast<std::size_t>(take(basis.size()));
    std::vector<PgGramEntry> gram(slot(0, count));
    std::vector<ExactRational> transform(gram.size());
    for (std::size_t i = 0; i < count; ++i) {
        gram[slot(i, i)] = {f.r("1"), f.r("1")};
        transform[slot(i, i)] = f.r("1");
    }
    auto caps = PgConditioningLimits{profile, steps, {}};
    const I dyadic = I{1} << 67;
    gram[slot(0, 1)].upper = rat(f, I{1}, dyadic);
    gram[slot(0, 2)].upper = rat(f, I{1}, dyadic);
    check_failure(verify(basis, gram, transform, caps, ws, reject),
                  core::FailureCode::PrecisionUnresolved);
    REQUIRE(basis.size());

    gram[slot(0, 2)].upper = {};
    core::ResourcePlan empty{0};
    check_failure(verify(basis, std::span{gram}.first(1), transform, caps, ws, empty),
                  core::FailureCode::InvalidInput);
    caps.cancellation = requested_stop_token();
    check_failure(verify(basis, gram, transform, caps, ws, empty), core::FailureCode::Cancelled);
    CHECK(empty.allocated_bytes() == 0);
    caps.cancellation = {};
    gram[slot(0, 1)].lower = f.r("1");
    check_failure(verify(basis, gram, transform, caps, ws, reject),
                  core::FailureCode::InvalidInput);
    gram[slot(0, 1)].lower = {};
    transform[slot(0, 0)] = f.r("1", "3");
    check_failure(verify(basis, gram, transform, caps, ws, reject),
                  core::FailureCode::InvalidInput);
    transform[slot(0, 0)] = {};
    check_failure(verify(basis, gram, transform, caps, ws, reject),
                  core::FailureCode::RankUnresolved);
    transform[slot(0, 0)] = f.r("1");

    auto other = flat_basis(ws, home);
    core::ResourcePlan control{profile};
    const auto accepted = take(verify(other, gram, transform, caps, ws, control));
    const auto info = take(accepted.info());
    CHECK(q(take(accepted.congruence_bound(query_limits().algebra, home))) == Q(I{1}) / dyadic);
    caps.maximum_payload_bytes = info.payload_bytes;
    caps.maximum_steps = info.steps - 1;
    core::ResourcePlan short_work{info.payload_bytes};
    check_failure(verify(basis, gram, transform, caps, ws, short_work),
                  core::FailureCode::ResourceExceeded);
    caps.maximum_steps = info.steps;
    core::ResourcePlan short_bytes{info.payload_bytes - 1};
    check_failure(verify(basis, gram, transform, caps, ws, short_bytes),
                  core::FailureCode::ResourceExceeded);
    CHECK(short_bytes.allocated_bytes() == 0);
    core::ResourcePlan equal{info.payload_bytes};
    CHECK(verify(basis, gram, transform, caps, ws, equal));
}

TEST_CASE("PG analytic Gram rejects unqualified shapes before consuming the reference",
          "[SH-D03][SH-D07][N4-08][negative][verification]") {
    Fixture f;
    const auto length = f.r("1");
    auto expected = core::FailureCode::Unsupported;
    bool cancelled = false;
    SECTION("varying height") {
        const auto zero = f.r("0"), one = f.r("1"), height = f.r("1", "32");
        const std::array<ExactInputTerm, 1> terms{{{{2, 0}, height}}};
        const auto graph = take(ExactChartPolynomial::make(terms, exact_limits(), f.input_plan));
        f.graph = take(PgGraph::prepare({zero, one, zero, graph, 1, PgGraphConvention::Version1},
                                        {exact_limits(), 64, {}}, f.input_plan));
    }
    bool no_steps = false;
    SECTION("zero integer work budget") {
        no_steps = true;
        expected = core::FailureCode::ResourceExceeded;
    }
    SECTION("cancellation before input visits") {
        cancelled = true;
        expected = core::FailureCode::Cancelled;
    }
    core::ResourcePlan setup{profile};
    core::Workspace ws{core::PlannedBuffer<std::byte>{}};
    auto basis = make_basis(f, length, ws, setup);
    auto caps = PgFlatGramLimits{query_limits(), no_steps ? 0 : steps};
    if (cancelled)
        caps.arithmetic.cancellation = requested_stop_token();
    core::ResourcePlan refused{profile};
    const auto prepare = [](PgReferenceBasis& source, const PgFlatGramLimits& limits,
                            core::Workspace& workspace, core::ResourcePlan& plan) {
        return PgFlatGram::prepare(std::move(source), dyadic_fractional_bits, limits, workspace,
                                   plan);
    };
    check_failure(prepare(basis, caps, ws, refused), expected);
    REQUIRE(basis.source());
    if (cancelled)
        CHECK(refused.allocated_bytes() == 0);
}

TEST_CASE("PG analytic Gram retains physical chart scales on adopted larger grids",
          "[SH-D03][N4-02][verification]") {
    for (const core::Index n : {2, 4}) {
        CAPTURE(n);
        Fixture f{n};
        const auto length = f.r("3", "2"), width = f.r("4", "3"), skew = f.r("1", "9");
        for (auto& value : f.heights)
            value = f.r("5", "4");
        const auto zero = f.r("0"), constant = f.r("1", "7");
        const std::array<ExactInputTerm, 1> terms{{{{0, 0}, constant}}};
        const auto height = take(ExactChartPolynomial::make(terms, exact_limits(), f.input_plan));
        f.graph = take(PgGraph::prepare({zero, width, skew, height, n, PgGraphConvention::Version1},
                                        {exact_limits(), 64, {}}, f.input_plan));
        core::Workspace ws{core::PlannedBuffer<std::byte>{}};
        core::ResourcePlan plan{profile};
        auto basis = make_basis(f, length, ws, plan);
        const auto count = static_cast<std::size_t>(take(basis.size()));
        const auto* rank = take(basis.rank());
        // Include early R/P combinations and a last C column at both grids.
        const auto selected = std::array<std::size_t, 4>{0, 1, count / 2, count - 1};
        const auto rows = take(rank->rows());
        const auto algebra = query_limits().algebra;
        std::vector<std::vector<Term>> independent(selected.size());
        for (std::size_t k = 0; k < selected.size(); ++k)
            for (std::size_t i = 0; i < rows.size(); ++i) {
                const auto value = take(rank->basis({static_cast<core::Index>(selected[k])},
                                                    static_cast<core::Index>(i), algebra, plan));
                if (value.sign() != 0)
                    independent[k].push_back({rows[i], q(value)});
            }
        const Q h("5/12");
        const auto expected =
            enclose(f, analytic_gram(independent, h, h * Q("4/3") / (16 * n * n)));
        const auto flat = take(PgFlatGram::prepare(std::move(basis), dyadic_fractional_bits,
                                                   {query_limits(), steps}, ws, plan));
        const auto supplied = take(flat.enclosure());
        REQUIRE(supplied.upper_triangle.size() == slot(0, count));
        for (std::size_t j = 0; j < selected.size(); ++j)
            for (std::size_t i = 0; i <= j; ++i) {
                const auto& actual = supplied.upper_triangle[slot(selected[i], selected[j])];
                CHECK(actual.lower == expected[slot(i, j)].lower);
                CHECK(actual.upper == expected[slot(i, j)].upper);
            }
    }
}

TEST_CASE("PG LCM meters both limbs of its division before the final product",
          "[SH-D03][N4-08][negative][verification]") {
    using namespace uvlas::structure::detail;
    const RankWide base = RankWide{1} << std::numeric_limits<RankWord>::digits;
    const RankWide left = 6 * base, right = 9 * base;

    // Three Euclidean remainders (6P, 3P, 0) each read two limbs from
    // each operand: 12. The following 6P / 3P costs four, and 2 * 9P costs two.
    constexpr std::uint64_t gcd_steps = 12, division_steps = 4, product_steps = 2;
    constexpr auto expected = gcd_steps + division_steps + product_steps;
    const auto attempt = [&](std::uint64_t maximum_steps) {
        const PgRankPreparationLimits caps{0, maximum_steps, 0, {}};
        auto ledger = RankBudget{caps};
        auto result = partition_entry(
            [&]() -> core::Result<RankWide> { return rank_lcm(ledger, left, right); });
        return std::pair{std::move(result), ledger.steps};
    };

    const auto [equal, equal_used] = attempt(expected);
    REQUIRE(equal);
    CHECK(*equal == 18 * base);
    CHECK(equal_used == expected);
    const auto [short_product, product_used] = attempt(expected - 1);
    check_failure(short_product, core::FailureCode::ResourceExceeded);
    CHECK(product_used == gcd_steps + division_steps);
    const auto [short_division, division_used] = attempt(gcd_steps + division_steps - 1);
    check_failure(short_division, core::FailureCode::ResourceExceeded);
    CHECK(division_used == gcd_steps);
}

TEST_CASE("PG projection rounds the actual exact radical and enforces the requested error",
          "[SH-D03][N4-02][N4-06][N4-08][negative][verification]") {
    Fixture f;
    core::ResourcePlan plan{profile};
    const auto limits = query_limits();
    auto math = detail::ReferenceArithmetic{limits, plan};
    const auto round = [&](const PgReferenceScalar& value, const ExactRational& g, double budget,
                           core::Precision precision) {
        return detail::partition_entry([&]() -> core::Result<detail::PgProjectedWord> {
            const auto tolerance = detail::reference_binary64(math, budget);
            return detail::pg_project_scalar(math, value, g, tolerance, precision);
        });
    };
    const auto one = f.r("1"), two = f.r("2"), four = f.r("4");
    const auto radical = PgReferenceScalar{{}, f.r("1")};
    const auto binary = core::Precision::Binary64;
    const auto two_word = core::Precision::TwoWord;
    const auto root = take(round(radical, two, 0x1p-104, two_word));
    // Independently rounded from a 140 decimal digit evaluation of sqrt(2).
    CHECK(root.high == 0x1.6a09e667f3bcdp+0);
    CHECK(root.low == -0x1.bdd3413b26456p-54);
    check_failure(round(radical, two, 0x1p-104, binary), core::FailureCode::PrecisionUnresolved);
    const auto zero = take(round({f.r("-2"), f.r("1")}, four, 0, two_word));
    CHECK(zero.high == 0);
    CHECK(zero.low == 0);

    const auto tie = PgReferenceScalar{f.r("9007199254740993", "9007199254740992"), {}};
    const auto at_bound = take(round(tie, one, 0x1p-53, binary));
    CHECK(at_bound.high == 1);
    check_failure(round(tie, one, std::nextafter(0x1p-53, 0.0), binary),
                  core::FailureCode::PrecisionUnresolved);
    const auto kept = take(round(tie, one, 0, two_word));
    CHECK(kept.high == 1);
    CHECK(kept.low == 0x1p-53);
    const auto negative_tie = PgReferenceScalar{f.r("-9007199254740995", "9007199254740992"), {}};
    CHECK(take(round(negative_tie, one, 0x1p-53, binary)).high == -0x1.0000000000002p+0);

    const I half_denominator = I{1} << 1075U;
    const auto neg_half = PgReferenceScalar{rat(f, I{-1}, half_denominator), {}};
    const auto odd_half = PgReferenceScalar{rat(f, I{3}, half_denominator), {}};
    CHECK(std::bit_cast<std::uint64_t>(take(round(neg_half, one, 0x1p-1074, binary)).high) ==
          0x8000000000000000ULL);
    CHECK(std::bit_cast<std::uint64_t>(take(round(odd_half, one, 0x1p-1074, binary)).high) == 2);
    check_failure(round(neg_half, one, 0, two_word), core::FailureCode::PrecisionUnresolved);

    const I overflow = (I{1} << 1024U) - (I{1} << 970U);
    const auto overflow_tie = PgReferenceScalar{rat(f, overflow, I{1}), {}};
    const auto below_tie = PgReferenceScalar{rat(f, overflow - 1, I{1}), {}};
    check_failure(round(overflow_tie, one, 0x1p970, binary),
                  core::FailureCode::PrecisionUnresolved);
    CHECK(take(round(below_tie, one, 0x1p970, binary)).high == std::numeric_limits<double>::max());
}

TEST_CASE("PG dyadic dot preserves cancellation and subnormal rounding",
          "[SH-D02][SH-D05][N4-06][verification]") {
    using namespace detail;
    const auto max = volume_decode(std::numeric_limits<double>::max());
    const auto min = volume_decode(std::numeric_limits<double>::denorm_min());
    const auto one = volume_decode(1);
    const auto large_product = volume_product(max, max);
    const auto subnormal = volume_product(min, one);
    const auto cancelled = large_product + subnormal - large_product;
    CHECK(take(volume_round(cancelled)) == std::numeric_limits<double>::denorm_min());
    check_failure(volume_round(large_product), core::FailureCode::PrecisionUnresolved);
    const auto half_subnormal = volume_product(min, volume_decode(0.5));
    CHECK(take(volume_round(half_subnormal)) == 0);
    CHECK(take(volume_round(subnormal + half_subnormal)) ==
          2 * std::numeric_limits<double>::denorm_min());
    const auto tail = volume_lift(volume_decode(0x1p-53));
    const auto tie = volume_lift(one) + tail;
    const auto paired = take(volume_finalize(tie, {}, {}, true));
    CHECK(paired.high == 1);
    CHECK(paired.low == 0x1p-53);
    check_failure(volume_finalize(tie, {}, {}, false), core::FailureCode::PrecisionUnresolved);
}

TEST_CASE("PG analytic Gram integrates a clipped constant-thickness plane",
          "[SH-D03][SH-D07][N4-02][N4-08][verification]") {
    for (const core::Index n : {1, 2, 4}) {
        CAPTURE(n);
        Fixture f{n};
        make_sloping_cut(f, I{4 * n});
        if (n == 1) {
            f.cut_points.push_back({f.r("0"), f.r("1", "3")});
            f.cut_points.push_back({f.r("1"), f.r("3", "7")});
            f.refinements.push_back({23, 4, f.cut_points[0], f.cut_points[1]});
        }
        core::ResourcePlan plan{profile};
        core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
        auto basis = make_basis(f, f.r("1"), workspace, plan);
        const auto count = static_cast<std::size_t>(take(basis.size()));
        const auto* rank = take(basis.rank());
        const auto rows = take(rank->rows());
        std::vector<std::size_t> selected;
        if (n == 1) {
            selected.resize(count);
            for (std::size_t j = 0; j < count; ++j)
                selected[j] = j;
        } else {
            selected = {0, 1, count / 2, count - 1};
        }
        std::vector<std::vector<Term>> independent(selected.size());
        const auto arithmetic = query_limits();
        for (std::size_t j = 0; j < selected.size(); ++j) {
            const PgRankBasisColumn col{static_cast<core::Index>(selected[j])};
            for (std::size_t row = 0; row < rows.size(); ++row) {
                const auto x =
                    take(rank->basis(col, static_cast<core::Index>(row), arithmetic.algebra, plan));
                if (x.sign() != 0)
                    independent[j].push_back({rows[row], q(x)});
            }
        }
        const auto exact = analytic_cut_gram(independent, n);
        const auto expected = enclose(f, exact);
        const auto caps = PgFlatGramLimits{arithmetic, steps};
        const auto prepare = [&](core::ResourcePlan& resources) {
            return PgFlatGram::prepare(std::move(basis), dyadic_fractional_bits, caps, workspace,
                                       resources);
        };
        if (n == 1) {
            core::ResourcePlan short_payload{1024};
            check_failure(prepare(short_payload), core::FailureCode::ResourceExceeded);
        }
        auto flat = take(prepare(plan));
        const auto supplied = take(flat.enclosure());
        REQUIRE(supplied.upper_triangle.size() == slot(0, count));
        for (std::size_t j = 0; j < selected.size(); ++j)
            for (std::size_t i = 0; i <= j; ++i) {
                const auto& actual = supplied.upper_triangle[slot(selected[i], selected[j])];
                CHECK(actual.lower == expected[slot(i, j)].lower);
                CHECK(actual.upper == expected[slot(i, j)].upper);
            }
        if (n == 1) {
            const auto transform = proposal(f, exact, count);
            const auto accepted = take(std::move(flat).certify(
                {transform, dyadic_fractional_bits}, {profile, steps, {}}, workspace, plan));
            CHECK(take(accepted.info()).rank == count);
            CHECK(q(take(accepted.congruence_bound(arithmetic.algebra, plan))) <=
                  Q("1/100000000000000000000"));
            REQUIRE(accepted.reference_basis());
        }
    }
}

TEST_CASE("PG analytic Gram rejects a cut outside the fixed-width moment scale",
          "[SH-D03][SH-D07][N4-02][N4-08][negative][verification]") {
    Fixture f;
    make_sloping_cut(f, I{1} << 190U);
    core::ResourcePlan plan{profile};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto basis = make_basis(f, f.r("1"), workspace, plan);
    const auto project = [&](std::uint64_t limit) {
        return PgFlatGram::prepare(std::move(basis), dyadic_fractional_bits,
                                   {query_limits(), limit}, workspace, plan);
    };
    check_failure(project(steps), core::FailureCode::ResourceExceeded);
    REQUIRE(basis.source());
}

TEST_CASE("PG analytic Gram uses rectangular moments on a cell-aligned subdomain",
          "[SH-D03][SH-D07][N4-02][verification]") {
    Fixture f;
    f.domain.clear();
    f.domain.push_back(f.point("0", "0"));
    f.domain.push_back({f.r("1", "2"), f.r("0")});
    f.domain.push_back({f.r("1", "2"), f.r("1", "2")});
    f.domain.push_back({f.r("0"), f.r("1", "2")});
    for (auto& height : f.heights)
        height = f.r("2");
    core::ResourcePlan plan{profile};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto basis = make_basis(f, f.r("1"), workspace, plan);
    const auto count = static_cast<std::size_t>(take(basis.size()));
    const auto independent = basis_terms(*take(basis.rank()), count, plan);
    for (const auto& column : independent)
        for (const auto& term : column)
            REQUIRE(term.row.cell == PgCellIndex{0, 0});
    const auto expected = enclose(f, analytic_gram(independent));
    const auto flat = take(PgFlatGram::prepare(std::move(basis), dyadic_fractional_bits,
                                               {query_limits(), steps}, workspace, plan));
    const auto actual = take(flat.enclosure()).upper_triangle;
    REQUIRE(actual.size() == expected.size());
    for (std::size_t k = 0; k < actual.size(); ++k) {
        CHECK(actual[k].lower == expected[k].lower);
        CHECK(actual[k].upper == expected[k].upper);
    }
}

TEST_CASE("PG analytic Gram integrates certified piecewise-affine thickness",
          "[SH-D03][SH-D07][N4-02][N4-08][verification]") {
    constexpr std::uint32_t enclosure_bits = 256;
    const std::array models{AffineThicknessCase{1, Q{}, Q{1}, Q{}, Q{1}, false},
                            AffineThicknessCase{1, Q{"1/8"}, Q{"7/8"}, Q{"1/16"}, Q{"15/16"}, true},
                            AffineThicknessCase{2, Q{"1/8"}, Q{"3/8"}, Q{"1/8"}, Q{"3/8"}, false},
                            AffineThicknessCase{4, Q{"1/8"}, Q{"3/8"}, Q{"1/8"}, Q{"3/8"}, false}};
    for (const auto& model : models) {
        CAPTURE(model.n, model.certify);
        Fixture f{model.n};
        set_affine_thickness_fixture(f, model);
        core::ResourcePlan plan{profile};
        core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
        auto basis = make_basis(f, f.r("3", "2"), workspace, plan);
        const auto count = static_cast<std::size_t>(take(basis.size()));
        const auto* rank = take(basis.rank());
        const auto rows = take(rank->rows());
        std::vector<std::size_t> selected;
        if (model.certify) {
            selected.resize(count);
            for (std::size_t j = 0; j < count; ++j)
                selected[j] = j;
        } else {
            selected = {0, 1, count / 2, count - 1};
        }
        std::vector<std::vector<Term>> independent(selected.size());
        const auto caps = PgFlatGramLimits{query_limits(), steps};
        for (std::size_t j = 0; j < selected.size(); ++j) {
            const PgRankBasisColumn col{static_cast<core::Index>(selected[j])};
            for (std::size_t row = 0; row < rows.size(); ++row) {
                const auto x = take(
                    rank->basis(col, static_cast<core::Index>(row), caps.arithmetic.algebra, plan));
                if (x.sign() != 0)
                    independent[j].push_back({rows[row], q(x)});
            }
        }
        const auto exact = analytic_affine_gram(independent, model);
        const auto expected = enclose(f, exact, enclosure_bits);
        const auto prepare = [&](core::ResourcePlan& resources) {
            return PgFlatGram::prepare(std::move(basis), enclosure_bits, caps, workspace,
                                       resources);
        };
        if (model.certify) {
            core::ResourcePlan short_payload{1024};
            check_failure(prepare(short_payload), core::FailureCode::ResourceExceeded);
        }
        auto flat = take(prepare(plan));
        const auto supplied = take(flat.enclosure());
        REQUIRE(supplied.upper_triangle.size() == slot(0, count));
        for (std::size_t j = 0; j < selected.size(); ++j)
            for (std::size_t i = 0; i <= j; ++i) {
                const auto& actual = supplied.upper_triangle[slot(selected[i], selected[j])];
                CHECK(actual.lower == expected[slot(i, j)].lower);
                CHECK(actual.upper == expected[slot(i, j)].upper);
            }
        if (model.certify) {
            std::vector<Q> midpoint;
            midpoint.reserve(expected.size());
            for (const auto& interval : expected)
                midpoint.emplace_back((q(interval.lower) + q(interval.upper)) / 2);
            const auto independent_transform = proposal(f, midpoint, count, enclosure_bits);
            const PgConditioningGeneration grids{512, enclosure_bits};
            const auto certify = [&](const PgConditioningLimits& limits) {
                return std::move(flat).generate_and_certify(grids, limits, workspace, plan);
            };
            check_failure(certify({0, steps, {}}), core::FailureCode::ResourceExceeded);
            REQUIRE(flat.enclosure());
            const auto accepted = take(certify({profile, steps, {}}));
            CHECK(take(accepted.info()).rank == count);
            CHECK(q(take(accepted.congruence_bound(caps.arithmetic.algebra, plan))) <=
                  Q("1/100000000000000000000"));
            check_generated_required(accepted, enclosure_bits, independent_transform, workspace,
                                     plan);
        }
    }
}

TEST_CASE("PG analytic Gram rejects affine thickness outside the fixed moment width",
          "[SH-D03][SH-D07][N4-02][N4-08][negative][verification]") {
    Fixture f;
    const I d = I{1} << 1700U;
    f.heights.front() = rat(f, d + 1, d);
    core::ResourcePlan plan{profile};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto basis = make_basis(f, f.r("1"), workspace, plan);
    const auto prepare = [&] {
        return PgFlatGram::prepare(std::move(basis), dyadic_fractional_bits,
                                   {query_limits(), steps}, workspace, plan);
    };
    check_failure(prepare(), core::FailureCode::ResourceExceeded);
    REQUIRE(basis.source());
}

TEST_CASE("PG generated factor certifies the original interval and uses one budget",
          "[SH-D03][N4-02][N4-08][negative][verification]") {
    // Synthetic dyadic unit control for R = [[2,3/4],[0,1/2]], padded
    // with identity. This exercises the algebra, not geometric enclosure provenance.
    Fixture f;
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan setup{profile};
    auto basis = flat_basis(workspace, setup);
    const auto rank = static_cast<std::size_t>(take(basis.size()));
    std::vector<PgGramEntry> gram(slot(0, rank));
    const auto exact_interval = [&](std::size_t k, std::string_view n, std::string_view d = "1") {
        gram[k] = {f.r(n, d), f.r(n, d)};
    };
    for (std::size_t j = 0; j < rank; ++j)
        exact_interval(slot(j, j), "1");
    exact_interval(0, "4");
    exact_interval(1, "3", "2");
    exact_interval(2, "13", "16");
    const PgConditioningGeneration grids{64, 16};
    const PgConditioningLimits full{profile, steps, {}};
    const auto certify = [&](PgConditioningGeneration generation,
                             const PgConditioningLimits& limits, core::ResourcePlan& plan,
                             std::span<const PgGramEntry> source) {
        return PgConditioningCertificate::generate_and_certify(std::move(basis), {source, 16},
                                                               generation, limits, workspace, plan);
    };
    core::ResourcePlan none{0}, refused{profile};
    check_failure(certify(grids, full, none, std::span{gram}.first(1)),
                  core::FailureCode::InvalidInput);
    check_failure(certify({1536, 16}, full, none, gram), core::FailureCode::ResourceExceeded);
    check_failure(certify(grids, {profile, steps, requested_stop_token()}, none, gram),
                  core::FailureCode::Cancelled);
    CHECK(none.allocated_bytes() == 0);

    gram.front() = {f.r("3"), f.r("5")};
    check_failure(certify(grids, full, refused, gram), core::FailureCode::PrecisionUnresolved);
    gram.front() = {f.r("5"), f.r("3")};
    check_failure(certify(grids, full, refused, gram), core::FailureCode::InvalidInput);
    exact_interval(0, "1", "3");
    check_failure(certify(grids, full, refused, gram), core::FailureCode::InvalidInput);
    exact_interval(0, "1", "16");
    check_failure(certify({0, 16}, full, refused, gram), core::FailureCode::PrecisionUnresolved);
    exact_interval(0, "4");
    REQUIRE(basis.size());

    core::ResourcePlan observed{profile};
    const auto accepted = take(certify(grids, full, observed, gram));
    const auto used = take(accepted.info());
    CHECK(used.payload_bytes == observed.allocated_bytes());
    const auto exact = query_limits().algebra;
    CHECK(take(accepted.congruence_bound(exact, setup)).sign() == 0);
    CHECK(q(take(accepted.transform_coefficient({0}, {0}, exact, setup))) == Q("1/2"));
    CHECK(q(take(accepted.transform_coefficient({0}, {1}, exact, setup))) == Q("-3/4"));
    CHECK(q(take(accepted.transform_coefficient({1}, {1}, exact, setup))) == 2);
    CHECK(take(accepted.transform_coefficient({1}, {0}, exact, setup)).sign() == 0);

    basis = flat_basis(workspace, setup);
    core::ResourcePlan one_byte_short{used.payload_bytes - 1};
    check_failure(certify(grids, {used.payload_bytes, used.steps, {}}, one_byte_short, gram),
                  core::FailureCode::ResourceExceeded);
    CHECK(one_byte_short.allocated_bytes() == 0);
    check_failure(certify(grids, {used.payload_bytes, used.steps - 1, {}}, refused, gram),
                  core::FailureCode::ResourceExceeded);
    core::ResourcePlan equality{used.payload_bytes};
    const auto replay = take(certify(grids, {used.payload_bytes, used.steps, {}}, equality, gram));
    CHECK(take(replay.info()).steps == used.steps);
    CHECK(take(replay.info()).payload_bytes == used.payload_bytes);
}
