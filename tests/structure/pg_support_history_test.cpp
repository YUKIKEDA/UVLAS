#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <string>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <uvlas/structure/pg_full_trace.hpp>
#include <uvlas/structure/pg_global_constraints.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>
#include <uvlas/structure/pg_support_history.hpp>

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

TEST_CASE("PG support history qualifies complete common intervals of dependent rows",
          "[SH-D05][SH-D06][N4-04][N4-08][verification]") {
    Fixture f;
    const auto zero = f.r("0"), one = f.r("1");
    std::array<ExactRational, 9> rotation;
    for (const auto i : {0U, 4U, 8U})
        rotation[i] = f.r("1");
    const std::array<ExactRational, 3> origin{};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan plan{profile};
    const auto model = patch(f, {rotation, origin, one}, workspace, plan);
    const auto range = vertical_edges(f, model, Q{0}, Q{1});
    const auto spans = edge_spans(range, zero, one);
    const auto columns = model.reference_coordinates().size();
    const std::vector<PgCartesianExactPosition> fixed_coordinates(columns);
    auto normal = query_limits();
    normal.maximum_scalar_operations = 300'000'000;
    normal.maximum_partition_visits = 100'000'000;
    const auto support =
        take(PgFullTrace::prepare_support(model, spans, take(model.info()).identity,
                                          fixed_coordinates, {normal, 30'000}, workspace, plan));
    const std::array patches{PgConstraintPatchBinding{numerics::PatchId{7}, model}};
    constexpr auto primary = PgConstraintId{58}, dependent = PgConstraintId{69};
    const std::array sources{PgConstraintTraceBinding{primary, support, 0, {}},
                             PgConstraintTraceBinding{dependent, support, 0, {}}};
    auto exact = normal;
    exact.algebra.maximum_value_bits = 65536;
    exact.algebra.maximum_work_bits = 262144;
    const auto prepare_chart = [&] {
        auto original = take(PgGlobalConstraints::prepare(
            patches, sources, {normal, 1, 2, 4096, 30'000}, workspace, plan));
        return take(PgSparseChart::prepare(std::move(original), {exact, 140, 140, 12'000},
                                           workspace, plan));
    };
    auto chart = prepare_chart();
    const auto rank = chart.pivots().size();
    REQUIRE(rank > 0);
    const auto& original = chart.original();
    const auto rows = static_cast<std::size_t>(original.sources()[0].rows);
    REQUIRE(original.original_rows().size() == 2 * rows);
    const auto selected = chart.pivots().front().dependent_column;
    std::vector<Q> spatial(rows);
    for (std::size_t row = 0; row < rows; ++row) {
        const auto terms = take(original.row(static_cast<core::Index>(row)));
        for (const auto& term : terms)
            if (term.column == selected)
                spatial[row] = q(term.coefficient);
    }
    REQUIRE(std::ranges::any_of(spatial, [](const auto& x) { return x != 0; }));
    const auto left = f.r("1", "4"), right = f.r("7", "4");
    const auto acts = [](const Q& t, bool after) -> std::array<Q, 3> {
        const Q d = t - 1;
        if (after)
            return {1 + 3 * d + 5 * d * d, 3 + 10 * d, Q{10}};
        return {t * t * t, 3 * t * t, 6 * t};
    };
    const auto make_track = [&](std::span<const Q> times, std::optional<std::size_t> mutate) {
        struct Sample {
            ExactRational t;
            std::vector<PgCartesianExactPosition> value, speed;
        };
        std::vector<Sample> data;
        data.reserve(times.size());
        for (std::size_t i = 0; i < times.size(); ++i) {
            const auto& t = times[i];
            const auto derivatives = acts(t, t > 1);
            Sample sample{rat(f, numerator(t), denominator(t)),
                          std::vector<PgCartesianExactPosition>(rows),
                          std::vector<PgCartesianExactPosition>(rows)};
            for (std::size_t row = 0; row < rows; ++row) {
                const Q value = spatial[row] * derivatives[0];
                Q speed = spatial[row] * derivatives[1];
                if (mutate == i)
                    speed += spatial[row];
                sample.value[row][0] = rat(f, numerator(value), denominator(value));
                sample.speed[row][0] = rat(f, numerator(speed), denominator(speed));
            }
            data.push_back(std::move(sample));
        }
        return data;
    };
    const auto views = [](const auto& samples) {
        std::vector<PgSupportTimeKnot> result;
        result.reserve(samples.size());
        for (const auto& s : samples)
            result.push_back({s.t, s.value, s.speed});
        return result;
    };
    const std::array a_times{Q{0}, Q{1}, Q{2}};
    const std::array b_times{Q{0}, Q{1} / 2, Q{1}, Q{3} / 2, Q{2}};
    const auto caps = PgSupportHistoryLimits{exact, exact, 8, 4 * (2 * rows) * 4};
    const auto prepare_history = [&](PgSparseChart& input,
                                     std::span<const PgSupportTimeTrack> tracks,
                                     const PgSupportHistoryLimits& budget) {
        return PgSupportHistory::prepare(std::move(input), {left, right}, tracks, budget, workspace,
                                         plan);
    };
    const auto prepared = [&] {
        const auto a = make_track(a_times, {}), b = make_track(b_times, {});
        const auto invalid = make_track(b_times, 3);
        const auto ak = views(a), bk = views(b), bad = views(invalid);
        const std::array tracks{PgSupportTimeTrack{dependent, bk}, PgSupportTimeTrack{primary, ak}};
        const std::array mismatch{PgSupportTimeTrack{dependent, bad},
                                  PgSupportTimeTrack{primary, ak}};
        const auto identity = take(chart.identity());
        check_failure(prepare_history(chart, mismatch, caps), core::FailureCode::InvalidInput);
        const std::array duplicate{tracks[1], tracks[1]};
        check_failure(prepare_history(chart, duplicate, caps), core::FailureCode::InvalidInput);
        const std::array uncovered{PgSupportTimeTrack{dependent, std::span{bk}.subspan(1)},
                                   tracks[1]};
        check_failure(prepare_history(chart, uncovered, caps), core::FailureCode::InvalidInput);
        const std::array repeated{ak[0], ak[0], ak[2]};
        const std::array decreasing{ak[0], ak[2], ak[1]};
        const std::array incomplete{
            ak[0], PgSupportTimeKnot{ak[1].seconds, ak[1].value.first(rows - 1), ak[1].rate},
            ak[2]};
        for (const auto candidate :
             {std::span{repeated}, std::span{decreasing}, std::span{incomplete}}) {
            const std::array invalid_track{PgSupportTimeTrack{primary, candidate}, tracks[0]};
            check_failure(prepare_history(chart, invalid_track, caps),
                          core::FailureCode::InvalidInput);
        }
        const std::array unknown{PgSupportTimeTrack{PgConstraintId{99}, ak}};
        check_failure(prepare_history(chart, unknown, caps), core::FailureCode::InvalidInput);
        auto short_coefficients = caps;
        --short_coefficients.maximum_coefficient_triples;
        check_failure(prepare_history(chart, tracks, short_coefficients),
                      core::FailureCode::ResourceExceeded);
        auto short_knots = caps;
        --short_knots.maximum_input_knots;
        check_failure(prepare_history(chart, tracks, short_knots),
                      core::FailureCode::ResourceExceeded);
        auto stopped = caps;
        stopped.arithmetic.cancellation = requested_stop_token();
        check_failure(prepare_history(chart, tracks, stopped), core::FailureCode::Cancelled);
        CHECK(take(chart.identity()) == identity);
        // Independently derive the two live sets and leaf allowance for this
        // full-width product trial. Working bytes at the two stages overlap.
        constexpr std::uint64_t word_bytes = sizeof(boost::multiprecision::limb_type);
        const auto limbs = caps.arithmetic.algebra.maximum_work_bits / (8 * word_bytes);
        constexpr std::uint64_t helper_rationals = 256;
        constexpr std::uint64_t event_slots = 10;
        const auto held = [&](std::uint64_t rationals) {
            return (2 * rationals + 1) * 4 * limbs * word_bytes;
        };
        const auto event_bytes = held(helper_rationals + event_slots);
        const auto combined_bytes =
            held(helper_rationals + event_slots + 3 * caps.maximum_coefficient_triples);
        const auto leaf_bytes = (51 * limbs + 6) * word_bytes;
        const auto minimum = combined_bytes + leaf_bytes;
        auto short_external = caps;
        short_external.arithmetic.algebra.external_memory = ExactEstimateAdmission{minimum - 1};
        check_failure(prepare_history(chart, tracks, short_external),
                      core::FailureCode::ResourceExceeded);
        CHECK(take(chart.identity()) == identity);

        auto combined_only = caps;
        combined_only.arithmetic.algebra.external_memory = ExactEstimateAdmission{minimum};
        auto double_subtraction_boundary = caps;
        double_subtraction_boundary.arithmetic.algebra.external_memory =
            ExactEstimateAdmission{combined_bytes + event_bytes};
        {
            auto repeat = prepare_chart();
            const auto alternative =
                take(prepare_history(repeat, tracks, double_subtraction_boundary));
            CHECK(alternative.events_seconds().size() == 5);
        }
        return take(prepare_history(chart, tracks, combined_only));
    }();
    const auto events = prepared.events_seconds();
    const std::array expected_events{Q{1} / 4, Q{1} / 2, Q{1}, Q{3} / 2, Q{7} / 4};
    REQUIRE(events.size() == expected_events.size());
    for (std::size_t i = 0; i < events.size(); ++i)
        CHECK(q(events[i]) == expected_events[i]);
    using Side = PgSupportTimeSide;
    using Order = PgSupportTimeOrder;
    const std::array orders{Order::Value, Order::Rate, Order::Acceleration};
    const auto compare = [&](const ExactRational& time, Side side) {
        const Q t = q(time);
        const auto reference = acts(t, t > 1 || (t == 1 && side == Side::Right));
        for (std::size_t derivative = 0; derivative < orders.size(); ++derivative) {
            const auto actual =
                take(prepared.evaluate(time, orders[derivative], side, exact, workspace, plan));
            REQUIRE(actual.size() == 2 * rows);
            for (std::size_t row = 0; row < rows; ++row) {
                CHECK(q(actual.values()[row][0]) == spatial[row] * reference[derivative]);
                CHECK(q(actual.values()[rows + row][0]) == spatial[row] * reference[derivative]);
                for (const auto offset : {std::size_t{0}, rows}) {
                    CHECK(actual.values()[offset + row][1].sign() == 0);
                    CHECK(actual.values()[offset + row][2].sign() == 0);
                }
            }
        }
    };
    compare(left, Side::Right);
    compare(f.r("3", "8"), Side::Interior);
    compare(f.r("1", "2"), Side::Left);
    compare(f.r("1", "2"), Side::Right);
    compare(one, Side::Left);
    compare(one, Side::Right);
    compare(right, Side::Left);
    const auto evaluate = [&](const ExactRational& t, Side side) {
        return prepared.evaluate(t, Order::Acceleration, side, exact, workspace, plan);
    };
    check_failure(evaluate(one, Side::Interior), core::FailureCode::InvalidInput);
    check_failure(evaluate(left, Side::Left), core::FailureCode::InvalidInput);
    check_failure(evaluate(right, Side::Right), core::FailureCode::InvalidInput);
    check_failure(evaluate(zero, Side::Interior), core::FailureCode::InvalidInput);
    const auto prescribed = take(evaluate(one, Side::Right));
    const auto id = take(prepared.chart().identity());
    const std::vector<PgCartesianExactPosition> independent(prepared.chart().free_columns().size());
    const auto lifted = take(prepared.chart().lift_rhs({id, independent}, {id, prescribed.values()},
                                                       exact, workspace, plan));
    const auto& full = prepared.chart().original();
    const auto recovered =
        take(full.apply({full.layout(), lifted.values()}, normal, workspace, plan));
    CHECK(std::ranges::equal(recovered.values(), prescribed.values()));
}
