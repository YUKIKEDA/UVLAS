#include "exact_polynomial_detail.hpp"
#include "pg_partition_test_support.hpp"

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <uvlas/structure/pg_conditioning.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using Q = detail::Rational;
using I = detail::Integer;
using Dec = boost::multiprecision::cpp_dec_float_100;
constexpr auto profile = std::uint64_t{48} << 30U;
constexpr auto steps = 10'000'000'000'000ULL;
constexpr std::uint32_t bits = 128;

PgReferenceLimits query_limits() {
    auto a = exact_limits();
    a.maximum_value_bits = 8192;
    a.maximum_work_bits = 32768;
    return {a, 20'000'000, 10'000'000, {}};
}

PgReferenceBasis flat_basis(core::Workspace& ws, core::ResourcePlan& plan) {
    Fixture f;
    for (auto& value : f.heights)
        value = f.r("2");
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
    const auto length = f.r("1");
    auto geometry = take(PgReferenceGeometry::prepare(std::move(input), {rotation, origin, length},
                                                      query_limits(), plan));
    const PgRankPreparationLimits rank_limits{profile, steps, 128, {}};
    const auto source = take(geometry.source());
    auto rank = take(
        PgRankPreparation::prepare(source.dictionary, source.partition, rank_limits, ws, plan));
    return take(
        PgReferenceBasis::prepare(std::move(geometry), std::move(rank), rank_limits, ws, plan));
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

Q integrate_terms(const Term& a, const Term& b) {
    // For H=2, L=1, a flat unit chart has zhat=zeta and det J=1.
    // The local half-cell coordinates introduce 1/16 and C carries 1/2.
    const auto i = static_cast<core::Index>(a.row.field);
    const auto j = static_cast<core::Index>(b.row.field);
    if (a.row.cell != b.row.cell || (i + j) % 2 != 0)
        return {};
    return a.value * b.value * moment(a.row.powers.first + b.row.powers.first) *
           moment(a.row.powers.second + b.row.powers.second) * moment(i + j) /
           (16 * (i == 2 ? 2 : 1) * (j == 2 ? 2 : 1));
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

std::vector<Q> analytic_gram(const std::vector<std::vector<Term>>& basis) {
    const auto count = basis.size();
    std::vector<Q> entries(slot(0, count));
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i)
            for (const auto& a : basis[i])
                for (const auto& b : basis[j])
                    entries[slot(i, j)] += integrate_terms(a, b);
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
    const I scale = I{1} << bits;
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
    return PgConditioningCertificate::certify(std::move(basis), {gram, bits}, {candidate, bits},
                                              caps, ws, plan);
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
    std::array<bool, 2> layer_seen{};
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
        layer_seen[layer] = true;
        check_jet(point, conditioned_jet(terms, transform, j, {0, 0}, s, t, Q("-1/3")));
    }
    REQUIRE(layer_seen[0]);
    REQUIRE(layer_seen[1]);

    const auto lateral = take(
        certificate.evaluate(PgReferenceLateralPoint{0, half, half}, column, caps, faces_plan));
    CHECK(std::holds_alternative<PgReferenceLateralMeasure>(lateral.reference.measure));
    const auto& r = lateral.reference;
    check_jet(lateral, conditioned_jet(terms, transform, j, r.cell, q(r.s), q(r.t), q(r.zeta)));
    check_failure(
        certificate.evaluate(PgReferenceLateralPoint{-1, half, half}, column, caps, faces_plan),
        core::FailureCode::InvalidInput);
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
    const auto prior = plan.allocated_bytes();
    const auto certificate =
        take(verify(basis, enclosure, transform, {profile, steps, {}}, ws, plan));
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
