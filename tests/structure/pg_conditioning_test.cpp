#include "exact_polynomial_detail.hpp"
#include "pg_partition_test_support.hpp"
#include "pg_rank_storage.hpp"

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <uvlas/structure/pg_flat_gram.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using Q = detail::Rational;
using I = detail::Integer;
using Dec = boost::multiprecision::cpp_dec_float_100;
constexpr auto profile = std::uint64_t{48} << 30U;
constexpr auto steps = 10'000'000'000'000ULL;
constexpr std::uint32_t dyadic_fractional_bits = 128;

PgReferenceLimits query_limits() {
    auto a = exact_limits();
    a.maximum_value_bits = 8192;
    a.maximum_work_bits = 32768;
    return {a, 20'000'000, 10'000'000, {}};
}

PgReferenceBasis make_basis(Fixture& f, const ExactRational& length, core::Workspace& ws,
                            core::ResourcePlan& plan) {
    auto dictionary = take(PgCoefficientDictionary::prepare(
        f.chart_graph(), {exact_limits(), 2000, 20000000, 100000, {}}, ws, plan));
    auto partition = take(PgIntegrationPartition::prepare(f.request(), limits(), ws, plan));
    PgReferenceInputs input{std::move(f.chart_graph()), std::move(dictionary),
                            std::move(partition)};
    std::array<ExactRational, 9> rotation;
    rotation[0] = f.r("1");
    rotation[4] = f.r("1");
    rotation[8] = f.r("1");
    const std::array<ExactRational, 3> origin{};
    auto geometry = take(PgReferenceGeometry::prepare(std::move(input), {rotation, origin, length},
                                                      query_limits(), plan));
    const PgRankPreparationLimits rank_limits{profile, steps, 128, {}};
    const auto source = take(geometry.source());
    auto rank = take(
        PgRankPreparation::prepare(source.dictionary, source.partition, rank_limits, ws, plan));
    return take(
        PgReferenceBasis::prepare(std::move(geometry), std::move(rank), rank_limits, ws, plan));
}

PgReferenceBasis flat_basis(core::Workspace& ws, core::ResourcePlan& plan) {
    Fixture f;
    for (auto& value : f.heights)
        value = f.r("2");
    const auto length = f.r("1");
    return make_basis(f, length, ws, plan);
}

const Q& q(const ExactRational& value) {
    return detail::ExactFactory::value(value);
}
ExactRational rat(Fixture& f, const I& n, const I& d) {
    return f.r(n.str(), d.str());
}
std::size_t slot(std::size_t i, std::size_t j) {
    return j * (j + 1) / 2 + i;
}

struct Term {
    PgRankRow row;
    Q value;
};
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

std::vector<std::vector<Term>> basis_terms(const PgRankPreparation& rank, std::size_t count,
                                           core::ResourcePlan& plan) {
    const auto rows = take(rank.rows());
    const auto algebra = query_limits().algebra;
    std::vector<std::vector<Term>> basis(count);
    for (std::size_t k = 0; k < count; ++k)
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto x = take(rank.basis({static_cast<core::Index>(k)},
                                           static_cast<core::Index>(i), algebra, plan));
            if (x.sign() != 0)
                basis[k].push_back({rows[i], q(x)});
        }
    return basis;
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

std::vector<PgGramEntry> enclose(Fixture& f, const std::vector<Q>& gram) {
    const I scale = I{1} << dyadic_fractional_bits;
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

std::vector<ExactRational> proposal(Fixture& f, const std::vector<Q>& gram, std::size_t count) {
    // Approximate test-side generator, independent of the checker.
    std::vector<Dec> factor(gram.size()), inverse(gram.size());
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i) {
            const auto& x = gram[slot(i, j)];
            Dec residual = numerator(x).convert_to<Dec>() / denominator(x).convert_to<Dec>();
            for (std::size_t k = 0; k < i; ++k)
                residual -= factor[slot(k, i)] * factor[slot(k, j)];
            if (i == j) {
                REQUIRE(residual > 0);
                factor[slot(i, j)] = sqrt(residual);
            } else {
                factor[slot(i, j)] = residual / factor[slot(i, i)];
            }
        }
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t k = j + 1; k > 0; --k) {
            const auto i = k - 1;
            Dec residual = i == j ? 1 : 0;
            for (std::size_t l = i + 1; l <= j; ++l)
                residual -= factor[slot(i, l)] * inverse[slot(l, j)];
            inverse[slot(i, j)] = residual / factor[slot(i, i)];
        }
    const I scale = I{1} << dyadic_fractional_bits;
    const auto decimal_scale = scale.convert_to<Dec>();
    std::vector<ExactRational> result;
    for (const auto& value : inverse) {
        const Dec integer = floor(value * decimal_scale + Dec{"0.5"});
        // Boost's decimal-to-integer conversion negates ilogb(0)'s sentinel.
        const auto units = integer == 0 ? I{} : integer.convert_to<I>();
        result.push_back(rat(f, units, scale));
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

template <class T>
void check_failure(const core::Result<T>& failure, core::FailureCode expected) {
    REQUIRE_FALSE(failure);
    CHECK(failure.error().code == expected);
}
Q power(const Q& x, core::Index exponent) {
    Q result{1};
    for (core::Index i = 0; i < exponent; ++i)
        result *= x;
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

} // namespace

TEST_CASE("PG full flat Gram certifies the transform and preserves all mandatory fields",
          "[SH-D03][N4-02][N4-08][verification]") {
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
    const auto certificate = take(certify_owned(supplied, plan));
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
    SECTION("nonconstant thickness") {
        f.heights.front() = f.r("2");
    }
    SECTION("partial cell coverage") {
        f.domain.clear();
        f.domain.push_back(f.point("0", "0"));
        f.domain.push_back({f.r("1", "2"), f.r("0")});
        f.domain.push_back({f.r("1", "2"), f.r("1", "2")});
        f.domain.push_back({f.r("0"), f.r("1", "2")});
    }
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
