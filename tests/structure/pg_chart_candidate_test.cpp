#include "../numerics/sparse_test_support.hpp"
#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <string>
#include <uvlas/structure/pg_chart_candidate.hpp>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <uvlas/structure/pg_full_trace.hpp>
#include <uvlas/structure/pg_global_constraints.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;

PgPreparedPatch patch(Fixture& f, const PgReferenceFrame& frame, core::Workspace& ws,
                      core::ResourcePlan& plan) {
    const auto limits = query_limits();
    auto basis = make_basis(f, frame, ws, plan);
    auto gram = take(PgFlatGram::prepare(std::move(basis), 256, {limits, steps}, ws, plan));
    auto accepted =
        take(std::move(gram).generate_and_certify({320, 246}, {profile, steps, {}}, ws, plan));
    return take(PgPreparedPatch::prepare(std::move(accepted), limits, ws, plan));
}

const PgIntegrationPartition& partition(const PgPreparedPatch& patch) {
    return take(take(patch.conditioning().reference_basis())->source()).partition;
}

struct EdgeRange {
    core::Index edge{};
    ExactRational first, last;
    bool increasing = true;
};

std::vector<EdgeRange> vertical_edges(Fixture& fixture, const PgPreparedPatch& patch, const Q& s,
                                      const Q& scale) {
    const auto& part = partition(patch);
    const auto vertices = take(part.vertices());
    const auto edges = take(part.edges());
    std::vector<EdgeRange> result;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& edge = edges[i];
        if (!edge.external)
            continue;
        const auto& a = vertices[static_cast<std::size_t>(edge.vertices[0])];
        const auto& b = vertices[static_cast<std::size_t>(edge.vertices[1])];
        if (q(a.s) != s || q(b.s) != s)
            continue;
        const Q x = q(a.t) * scale, y = q(b.t) * scale;
        const auto lo = std::min(x, y), hi = std::max(x, y);
        result.push_back({static_cast<core::Index>(i), rat(fixture, numerator(lo), denominator(lo)),
                          rat(fixture, numerator(hi), denominator(hi)), x < y});
    }
    std::ranges::sort(result, [](const auto& a, const auto& b) { return q(a.first) < q(b.first); });
    return result;
}

std::vector<PgTraceEdgeSpan> edge_spans(std::span<const EdgeRange> ranges,
                                        const ExactRational& zero, const ExactRational& one,
                                        std::optional<core::Index> layer = {}) {
    std::vector<PgTraceEdgeSpan> result;
    result.reserve(ranges.size());
    for (const auto& edge : ranges)
        result.push_back({edge.edge, edge.first, edge.last, edge.increasing ? zero : one,
                          edge.increasing ? one : zero, layer});
    return result;
}

} // namespace

namespace {
using namespace numerics;
namespace linear = numerics::linear;
using Buffer = PrecisionBuffer<double>;

Buffer words(std::size_t rows, std::uint64_t generation, int shift, bool nonzero, bool two_word,
             core::ResourcePlan& plan) {
    const auto size = static_cast<core::Index>(rows);
    const DenseShape shape{size, 3, size};
    std::vector<double> hi(rows * 3), lo(hi.size());
    if (nonzero)
        for (std::size_t i = 0; i < hi.size(); ++i) {
            hi[i] = static_cast<double>(static_cast<int>(i % 5U) - 2 + shift);
            lo[i] = (i % 2U == 0 ? 0x1p-54 : -0x1p-54);
        }
    auto high = take(DenseReal::copy(shape, hi, plan, core::Operation::EvaluateStructure));
    if (!two_word)
        return take(Buffer::binary64(std::move(high), core::NumericalGeneration{generation}));
    auto low = take(DenseReal::copy(shape, lo, plan, core::Operation::EvaluateStructure));
    return take(
        Buffer::two_word(std::move(high), std::move(low), core::NumericalGeneration{generation}));
}

Q qwords(Fixture& f, const Buffer& x, std::size_t i) {
    const auto as_q = [&](double v) {
        const auto exact = take(ExactRational::from_binary64(v, exact_limits(), f.input_plan));
        return q(exact);
    };
    if (const auto* binary = std::get_if<Binary64<double>>(&x.storage()))
        return as_q(binary->values.values()[i]);
    const auto& pair = std::get<TwoWord<double>>(x.storage());
    return as_q(pair.hi.values()[i]) + as_q(pair.lo.values()[i]);
}

std::vector<PgCartesianExactPosition> exact_words(Fixture& f, const Buffer& value) {
    const auto rows = static_cast<std::size_t>(value.shape().rows);
    std::vector<PgCartesianExactPosition> result(rows);
    for (std::size_t axis = 0; axis < 3; ++axis)
        for (std::size_t row = 0; row < rows; ++row) {
            const auto x = qwords(f, value, axis * rows + row);
            result[row][axis] = rat(f, numerator(x), denominator(x));
        }
    return result;
}
} // namespace

TEST_CASE("PG numerical chart keeps the exact sparse pattern and checks original work",
          "[SH-D05][SH-D06][N4-04][N4-05][N4-06][verification]") {
    Fixture a, b;
    const auto zero = a.r("0"), one = a.r("1");
    std::array<ExactRational, 9> rotation;
    for (const auto i : {0U, 4U, 8U})
        rotation[i] = a.r("1");
    const std::array<ExactRational, 3> oa{};
    const std::array<ExactRational, 3> ob{b.r("1"), {}, {}};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan plan{profile};
    const auto left = patch(a, {rotation, oa, one}, workspace, plan);
    const auto right = patch(b, {rotation, ob, one}, workspace, plan);
    const auto a0 = vertical_edges(a, left, Q{0}, Q{1});
    const auto a1 = vertical_edges(a, left, Q{1}, Q{1});
    const auto b0 = vertical_edges(b, right, Q{0}, Q{1});
    const auto ar = edge_spans(a1, zero, one), bl = edge_spans(b0, zero, one);
    const std::vector<PgCartesianExactPosition> fixed_coordinates(
        left.reference_coordinates().size());
    auto arithmetic = query_limits();
    arithmetic.maximum_scalar_operations = 300'000'000;
    arithmetic.maximum_partition_visits = 100'000'000;
    const auto connection =
        take(PgFullTrace::prepare(left, ar, right, bl, {arithmetic, 30000}, workspace, plan));
    const auto support = take(
        PgFullTrace::prepare_support(left, edge_spans(a0, zero, one), take(left.info()).identity,
                                     fixed_coordinates, {arithmetic, 30000}, workspace, plan));
    const std::array patches{PgConstraintPatchBinding{numerics::PatchId{17}, left},
                             PgConstraintPatchBinding{numerics::PatchId{72}, right}};
    const std::array traces{PgConstraintTraceBinding{PgConstraintId{11}, connection, 0, 1},
                            PgConstraintTraceBinding{PgConstraintId{40}, support, 0, {}}};
    auto original = take(PgGlobalConstraints::prepare(
        patches, traces, {arithmetic, 2, 2, 4096, 30000}, workspace, plan));
    auto exact_limits = arithmetic;
    exact_limits.algebra.maximum_value_bits = 65536;
    exact_limits.algebra.maximum_work_bits = 262144;
    auto reference = take(PgSparseChart::prepare(std::move(original),
                                                 {exact_limits, 320, 480, 30000}, workspace, plan));
    const auto id = take(reference.identity());
    const auto info = take(reference.original().info());
    const auto rank = static_cast<core::Index>(reference.pivots().size());
    REQUIRE(rank > 0);
    sparse_test::Fixture scales{rank}, transposed_scales{rank};
    auto transposed = sparse_test::context();
    transposed.unknown_layout = core::ScopeId{601};
    transposed.residual_layout = core::ScopeId{602};
    transposed_scales.unknown[0] = 2.0;
    const auto fwd = scales.scales(), rev = transposed_scales.scales(transposed);
    const auto solver_limits = sparse_test::limits(rank, 3);
    constexpr double coefficient_error = 0x1p-30;
    const auto prepare = [&](PgSparseChart& chart, double epsilon) {
        return PgChartCandidate::prepare(std::move(chart), fwd, rev,
                                         {arithmetic, epsilon, solver_limits}, plan);
    };
    check_failure(prepare(reference, 0.0), core::FailureCode::PrecisionUnresolved);
    CHECK(take(reference.identity()) == id);
    auto candidate = take(prepare(reference, coefficient_error));
    const auto independent = candidate.exact().free_columns().size();
    const auto r = words(independent, 71, 0, true, true, plan);
    const auto c = words(static_cast<std::size_t>(info.rows), 72, 0, false, false, plan);
    const auto f = words(static_cast<std::size_t>(info.columns), 73, 1, true, true, plan);
    const auto r_exact = exact_words(a, r);
    const auto f_exact = exact_words(a, f);
    const auto exact_lift =
        take(candidate.exact().lift({id, r_exact}, exact_limits, workspace, plan));
    const auto exact_dual =
        take(candidate.exact().transpose({id, f_exact}, exact_limits, workspace, plan));
    const auto request = [&](const linear::Scales& scale) {
        return PgChartCandidateSolve{scale,
                                     scales.precision(),
                                     sparse_test::accuracy(),
                                     solver_limits,
                                     core::RhsGeneration{88},
                                     core::NumericalGeneration{90},
                                     core::NumericalGeneration{110},
                                     1e-5,
                                     2'000'000,
                                     {}};
    };
    auto primal = request(fwd);
    auto dual_request = request(rev);
    const auto lift = take(candidate.lift({id, r}, {id, c}, primal, workspace, plan));
    const auto dual = take(candidate.transpose({id, f}, dual_request, workspace, plan));
    CHECK(lift.coordinates.precision() == core::Precision::TwoWord);
    CHECK(dual.independent_force.precision() == core::Precision::TwoWord);
    CHECK(lift.original_residual_upper <= primal.maximum_absolute_original_defect);
    CHECK(dual.original_balance_upper <= dual_request.maximum_absolute_original_defect);
    constexpr double work_budget = 1e-6;
    const auto work =
        take(candidate.work_defect({id, f}, {id, r}, {id, c}, lift, dual, work_budget, 10'000, {}));
    CHECK(work <= work_budget);

    const Q comparison_error(1, 1'000'000);
    const auto compare = [&](const Buffer& value,
                             std::span<const PgCartesianExactPosition> expected) {
        for (std::size_t axis = 0; axis < 3; ++axis)
            for (std::size_t i = 0; i < expected.size(); ++i)
                CHECK(abs(qwords(a, value, axis * expected.size() + i) - q(expected[i][axis])) <=
                      comparison_error);
    };
    compare(lift.coordinates, exact_lift.values());
    compare(dual.independent_force, exact_dual.independent_force.values());
    compare(dual.original_multipliers, exact_dual.original_multipliers.values());
    for (std::size_t i = 0; i < independent; ++i) {
        const auto j = static_cast<std::size_t>(candidate.exact().free_columns()[i]);
        CHECK(qwords(a, lift.coordinates, j) == qwords(a, r, i));
    }
    const auto exact_bound = [&](double upper) {
        const auto value =
            take(ExactRational::from_binary64(upper, partition_test::exact_limits(), a.input_plan));
        return q(value);
    };
    const auto q_candidate = exact_words(a, lift.coordinates);
    const auto lambda_candidate = exact_words(a, dual.original_multipliers);
    const auto h_candidate = exact_words(a, dual.independent_force);
    std::vector<std::array<Q, 3>> original_balance(f_exact.size());
    for (std::size_t col = 0; col < f_exact.size(); ++col)
        for (std::size_t axis = 0; axis < 3; ++axis)
            original_balance[col][axis] = q(f_exact[col][axis]);
    Q original_defect = 0;
    for (std::size_t row = 0; row < lambda_candidate.size(); ++row) {
        std::array<Q, 3> sum{};
        for (const auto& term :
             take(candidate.exact().original().row(static_cast<core::Index>(row)))) {
            const auto column = static_cast<std::size_t>(term.column);
            const auto coefficient = q(term.coefficient);
            for (std::size_t axis = 0; axis < 3; ++axis) {
                sum[axis] += coefficient * q(q_candidate[column][axis]);
                original_balance[column][axis] -= coefficient * q(lambda_candidate[row][axis]);
            }
        }
        for (const auto& value : sum) {
            const Q magnitude = abs(value);
            original_defect = std::max(original_defect, magnitude);
        }
    }
    CHECK(original_defect <= exact_bound(lift.original_residual_upper));
    for (std::size_t col = 0; col < independent; ++col)
        for (std::size_t axis = 0; axis < 3; ++axis)
            original_balance[static_cast<std::size_t>(candidate.exact().free_columns()[col])]
                            [axis] -= q(h_candidate[col][axis]);
    Q original_work = 0;
    for (std::size_t col = 0; col < f_exact.size(); ++col)
        for (std::size_t axis = 0; axis < 3; ++axis) {
            CHECK(abs(original_balance[col][axis]) <= exact_bound(dual.original_balance_upper));
            original_work += q(f_exact[col][axis]) * q(q_candidate[col][axis]);
        }
    for (std::size_t col = 0; col < independent; ++col)
        for (std::size_t axis = 0; axis < 3; ++axis)
            original_work -= q(r_exact[col][axis]) * q(h_candidate[col][axis]);
    CHECK(abs(original_work) <= exact_bound(work));
    primal.maximum_absolute_original_defect = 0;
    check_failure(candidate.lift({id, r}, {id, c}, primal, workspace, plan),
                  core::FailureCode::PrecisionUnresolved);
    auto stale = id;
    stale.bytes[0] ^= std::byte{1};
    check_failure(candidate.transpose({stale, f}, dual_request, workspace, plan),
                  core::FailureCode::StaleGeneration);
    dual_request.maximum_visits = 1;
    check_failure(candidate.transpose({id, f}, dual_request, workspace, plan),
                  core::FailureCode::ResourceExceeded);
    dual_request.cancellation = requested_stop_token();
    check_failure(candidate.transpose({id, f}, dual_request, workspace, plan),
                  core::FailureCode::Cancelled);
}
