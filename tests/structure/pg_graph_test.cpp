#include "pg_graph_detail.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <limits>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>
#include <uvlas/structure/pg_graph.hpp>
namespace {
using namespace uvlas;
using namespace uvlas::structure;
constexpr auto budget = std::uint64_t{64} * 1024 * 1024;
PgGraphLimits limits() {
    return {{10000, 1000, 10000, 100, 200, 4096, 1000000, ExactEstimateAdmission{budget},
             ExactTimingMode::Disabled},
            64,
            {}};
}
template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw std::runtime_error("Unexpected graph failure " +
                                 std::to_string(static_cast<unsigned>(result.error().code)));
    return std::move(*result);
}
ExactRational r(std::string_view p, std::string_view q, core::ResourcePlan& plan) {
    return take(ExactRational::from_integer_pair(p, q, limits().algebra, plan));
}
void fraction(const ExactRational& value, std::string_view p, std::string_view q) {
    const auto& backend = detail::ExactFactory::value(value).backend();
    CHECK(backend.num().str(0, std::ios_base::dec) == p);
    CHECK(backend.denom().str(0, std::ios_base::dec) == q);
}
PgGraph graph(core::Index n, core::ResourcePlan& plan) {
    const auto x = r("1", "7", plan), w = r("4", "5", plan), a = r("3", "20", plan);
    const auto one = r("1", "1", plan), mixed = r("2", "3", plan), fourth = r("-2", "5", plan),
               constant = r("1", "11", plan);
    const std::array terms{ExactInputTerm{{4, 0}, one}, ExactInputTerm{{2, 1}, mixed},
                           ExactInputTerm{{0, 4}, fourth}, ExactInputTerm{{0, 0}, constant}};
    const auto f = take(ExactChartPolynomial::make(terms, limits().algebra, plan));
    return take(PgGraph::prepare({x, w, a, f, n, PgGraphConvention::Version1}, limits(), plan));
}
template <class T>
concept RvalueView = requires(const T&& t) { std::move(t).view(); };
static_assert(!RvalueView<PgGraph> && !RvalueView<PgCellPoint> && !RvalueView<PgCellPolynomials> &&
              !RvalueView<PgGraphValue>);
static_assert(!std::is_default_constructible_v<QualifiedRectangularActivityMap>);
static_assert(!std::is_copy_constructible_v<PgGraph>);
} // namespace
TEST_CASE("PG graph control uses exact graph and half-cell anchor literals",
          "[SH-D03][SH-D07][N4-02][N4-08][verification]") {
    core::ResourcePlan plan{budget};
    const auto g = graph(1, plan);
    const ExactRational zero;
    const auto p =
        take(g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, limits(), plan));
    const auto point = take(p.view());
    fraction(point.s, "1", "4");
    fraction(point.t, "1", "4");
    const auto v = take(g.evaluate(p, {0, 0}, PgDerivativeFrame::ChartST, limits(), plan));
    const auto value = take(v.view());
    REQUIRE(value.components.size() == 3);
    fraction(value.components[0], "213", "560");
    fraction(value.components[1], "1", "4");
    fraction(value.components[2], "4379", "42240");
    const auto d = take(g.evaluate(p, {4, 0}, PgDerivativeFrame::ChartST, limits(), plan));
    fraction(take(d.view()).components[2], "24", "1");
    const auto u = take(g.evaluate(p, {4, 0}, PgDerivativeFrame::CellUV, limits(), plan));
    fraction(take(u.view()).components[2], "3", "32");
    const auto anchor = take(g.anchor_at(point.t, PgSide::Interior, 0, limits(), plan));
    const auto av = take(anchor.view());
    fraction(av.components[0], "101", "560");
    fraction(av.components[2], "629", "7040");
    const auto coefficients = take(g.cell_polynomials({0, 0}, limits(), plan));
    const auto polys = take(coefficients.view());
    for (std::size_t i = 0; i < 3; ++i) {
        const auto actual = take(polys.components[i].evaluate(zero, zero, limits().algebra, plan));
        CHECK(actual == value.components[i]);
    }
    const auto activity = take(QualifiedRectangularActivityMap::prepare(g, limits(), plan));
    REQUIRE(activity.validate_for(g));
    fraction(take(activity.cell_area(limits(), plan)), "1", "4");
    fraction(take(activity.total_area(limits(), plan)), "1", "1");
}

namespace {
using Ref = boost::multiprecision::cpp_rational;
Ref ref(std::int64_t numerator, std::int64_t denominator = 1) {
    return Ref{numerator} / Ref{denominator};
}
Ref power(const Ref& x, unsigned n) {
    Ref result = 1;
    for (unsigned i = 0; i < n; ++i)
        result *= x;
    return result;
}
void reference_is(const ExactRational& actual, const Ref& expected) {
    fraction(actual, boost::multiprecision::numerator(expected).str(0, std::ios_base::dec),
             boost::multiprecision::denominator(expected).str(0, std::ios_base::dec));
}
void position_is(const PgGraphValue& actual, const Ref& s, const Ref& t) {
    const auto values = take(actual.view());
    reference_is(values.components[0], ref(1, 7) + ref(4, 5) * s + ref(3, 20) * t);
    reference_is(values.components[1], t);
    reference_is(values.components[2],
                 power(s, 4) + ref(2, 3) * s * s * t - ref(2, 5) * power(t, 4) + ref(1, 11));
}
std::string digest(const PgGraphContentIdentity& id) {
    constexpr auto alphabet = "0123456789abcdef";
    std::string result;
    result.reserve(64);
    for (const auto byte : id.bytes) {
        const auto b = std::to_integer<unsigned char>(byte);
        result.push_back(alphabet[b >> 4U]);
        result.push_back(alphabet[b & 15U]);
    }
    return result;
}
PgGraph flat_graph(core::Index n, core::ResourcePlan& plan) {
    const auto x = r("1", "7", plan), w = r("4", "5", plan), a = r("3", "20", plan);
    const ExactChartPolynomial f;
    return take(PgGraph::prepare({x, w, a, f, n, PgGraphConvention::Version1}, limits(), plan));
}
} // namespace
TEST_CASE("Every adopted half-cell has independent rational geometry and derivative scales",
          "[SH-D03][SH-D07][N4-02][verification]") {
    for (const core::Index n : {1, 2, 4}) {
        core::ResourcePlan plan{budget};
        const auto g = graph(n, plan);
        const auto definition = take(g.view());
        REQUIRE(definition.cells.size() == static_cast<std::size_t>(4 * n * n));
        const ExactRational zero;
        for (core::Index a = 0; a < 2 * n; ++a)
            for (core::Index b = 0; b < 2 * n; ++b) {
                CAPTURE(n, a, b);
                const PgCellIndex cell{a, b};
                CHECK(definition.cells[static_cast<std::size_t>(a * (2 * n) + b)] == cell);
                const auto p = take(g.from_local(cell, zero, zero, PgSide::Interior,
                                                 PgSide::Interior, limits(), plan));
                const auto pv = take(p.view());
                const auto s = ref(2 * a + 1, 4 * n), t = ref(2 * b + 1, 4 * n);
                reference_is(pv.s, s);
                reference_is(pv.t, t);
                const auto inverse =
                    take(g.locate(pv.s, pv.t, PgSide::Interior, PgSide::Interior, limits(), plan));
                const auto iv = take(inverse.view());
                CHECK(iv.cell == cell);
                fraction(iv.u, "0", "1");
                fraction(iv.v, "0", "1");
                const auto value =
                    take(g.evaluate(p, {0, 0}, PgDerivativeFrame::ChartST, limits(), plan));
                position_is(value, s, t);
                const auto first =
                    take(g.evaluate(p, {1, 0}, PgDerivativeFrame::ChartST, limits(), plan));
                const auto fv = take(first.view());
                fraction(fv.components[0], "4", "5");
                fraction(fv.components[1], "0", "1");
                reference_is(fv.components[2], 4 * s * s * s + ref(4, 3) * s * t);
                const auto first_t =
                    take(g.evaluate(p, {0, 1}, PgDerivativeFrame::ChartST, limits(), plan));
                const auto ft = take(first_t.view());
                const auto height_t = ref(2, 3) * s * s - ref(8, 5) * t * t * t;
                reference_is(ft.components[0], ref(3, 20));
                reference_is(ft.components[1], ref(1));
                reference_is(ft.components[2], height_t);
                const auto local_t =
                    take(g.evaluate(p, {0, 1}, PgDerivativeFrame::CellUV, limits(), plan));
                const auto lt = take(local_t.view());
                reference_is(lt.components[0], ref(3, 20) / ref(4 * n));
                reference_is(lt.components[1], ref(1, 4 * n));
                reference_is(lt.components[2], height_t / ref(4 * n));
                const auto mixed =
                    take(g.evaluate(p, {2, 1}, PgDerivativeFrame::CellUV, limits(), plan));
                const auto mv = take(mixed.view());
                reference_is(mv.components[2], ref(4, 3) / power(ref(4 * n), 3));
                const auto fourth =
                    take(g.evaluate(p, {0, 4}, PgDerivativeFrame::CellUV, limits(), plan));
                reference_is(take(fourth.view()).components[2], ref(-48, 5) / power(ref(4 * n), 4));
                const auto coefficients = take(g.cell_polynomials(cell, limits(), plan));
                const auto cv = take(coefficients.view());
                const auto vv = take(value.view());
                for (std::size_t i = 0; i < 3; ++i) {
                    const auto scalar =
                        take(cv.components[i].evaluate(zero, zero, limits().algebra, plan));
                    CHECK(scalar == vv.components[i]);
                }
                const auto activity =
                    take(QualifiedRectangularActivityMap::prepare(g, limits(), plan));
                REQUIRE(activity.validate_for(g));
                reference_is(take(activity.cell_area(limits(), plan)), ref(1, 4 * n * n));
                fraction(take(activity.total_area(limits(), plan)), "1", "1");
            }
    }
}
TEST_CASE("Half-cell boundary choices are exact and never average neighboring coordinates",
          "[SH-D03][N4-02][negative][verification]") {
    for (const core::Index n : {1, 2, 4}) {
        core::ResourcePlan plan{budget};
        const auto g = graph(n, plan);
        const auto minus = r("-1", "1", plan), plus = r("1", "1", plan);
        for (core::Index a = 0; a < 2 * n; ++a)
            for (core::Index b = 0; b < 2 * n; ++b)
                for (const auto* u : {&minus, &plus})
                    for (const auto* v : {&minus, &plus}) {
                        const auto su = u->sign() < 0 ? PgSide::Upper : PgSide::Lower,
                                   st = v->sign() < 0 ? PgSide::Upper : PgSide::Lower;
                        const auto point =
                            take(g.from_local({a, b}, *u, *v, su, st, limits(), plan));
                        const auto pv = take(point.view());
                        const auto located = take(g.locate(pv.s, pv.t, su, st, limits(), plan));
                        const auto lv = take(located.view());
                        CHECK(lv.cell == PgCellIndex{a, b});
                        CHECK(lv.u == *u);
                        CHECK(lv.v == *v);
                    }
        const auto t = r("1", std::to_string(4 * n), plan);
        for (core::Index j = 1; j < 2 * n; ++j) {
            const auto s = r(std::to_string(j), std::to_string(2 * n), plan);
            const auto lower =
                take(g.locate(s, t, PgSide::Lower, PgSide::Interior, limits(), plan));
            const auto upper =
                take(g.locate(s, t, PgSide::Upper, PgSide::Interior, limits(), plan));
            const auto lo = take(lower.view()), hi = take(upper.view());
            CHECK(lo.cell.first == j - 1);
            CHECK(hi.cell.first == j);
            fraction(lo.u, "1", "1");
            fraction(hi.u, "-1", "1");
            const auto lvalue =
                take(g.evaluate(lower, {0, 0}, PgDerivativeFrame::ChartST, limits(), plan));
            position_is(lvalue, ref(j, 2 * n), ref(1, 4 * n));
            core::ResourcePlan rejected{budget};
            const auto unspecified =
                g.locate(s, t, PgSide::Interior, PgSide::Interior, limits(), rejected);
            REQUIRE_FALSE(unspecified);
            CHECK(unspecified.error().code == core::FailureCode::InvalidInput);
            CHECK(rejected.allocated_bytes() == 0);
        }
        const auto before =
            r("633825300114114700748351602687", "1267650600228229401496703205376", plan);
        const auto after =
            r("633825300114114700748351602689", "1267650600228229401496703205376", plan);
        const auto left =
            take(g.locate(before, t, PgSide::Interior, PgSide::Interior, limits(), plan));
        const auto right =
            take(g.locate(after, t, PgSide::Interior, PgSide::Interior, limits(), plan));
        CHECK(take(left.view()).cell.first == n - 1);
        CHECK(take(right.view()).cell.first == n);
        const auto tiny = take(ExactRational::from_binary64(
            std::numeric_limits<double>::denorm_min(), limits().algebra, plan));
        const auto tiny_point =
            take(g.locate(tiny, t, PgSide::Interior, PgSide::Interior, limits(), plan));
        const auto tp = take(tiny_point.view());
        CHECK(tp.cell.first == 0);
        CHECK(tp.s == tiny);
        const Ref minsub = Ref{1} / Ref{boost::multiprecision::cpp_int{1} << 1074};
        reference_is(tp.u, Ref{-1} + Ref{4 * n} * minsub);
    }
}
TEST_CASE("Graph canonical identities match independent versioned byte fixtures",
          "[SH-D03][SH-D07][N4-02][verification]") {
    std::ifstream stream(UVLAS_PG_GRAPH_FIXTURE);
    REQUIRE(stream.good());
    nlohmann::json fixture;
    stream >> fixture;
    REQUIRE(fixture.at("schema") == "uvlas.pg-graph-reference/1");
    REQUIRE(fixture.at("id") == "pg-graph-canonical-v1");
    for (const auto& row : fixture.at("identities")) {
        core::ResourcePlan plan{budget};
        const auto n = row.at("n").get<core::Index>();
        const auto kind = row.at("kind").get<std::string>();
        REQUIRE((kind == "flat" || kind == "quartic"));
        const auto g = kind == "flat" ? flat_graph(n, plan) : graph(n, plan);
        const auto v = take(g.view());
        CHECK(digest(v.identity) == row.at("sha256").get<std::string>());
        const auto copied = take(g.clone(limits(), plan));
        CHECK(take(copied.view()).identity == v.identity);
    }
    core::ResourcePlan plan{budget};
    const auto g = graph(1, plan);
    const auto base = take(g.view());
    const auto equivalent_x = r("2", "14", plan), equivalent_w = r("8", "10", plan),
               equivalent_a = r("6", "40", plan);
    const auto same = take(PgGraph::prepare(
        {equivalent_x, equivalent_w, equivalent_a, base.height, 1, PgGraphConvention::Version1},
        limits(), plan));
    CHECK(take(same.view()).identity == base.identity);
    const auto zero = r("0", "1", plan), one = r("1", "1", plan);
    const auto point =
        take(g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, limits(), plan));
    const auto activity = take(QualifiedRectangularActivityMap::prepare(g, limits(), plan));
    const std::array constant_term{ExactInputTerm{{0, 0}, one}};
    const auto constant = take(ExactChartPolynomial::make(constant_term, limits().algebra, plan));
    const auto changed_f = take(base.height.add(constant, limits().algebra, plan));
    for (unsigned field = 0; field < 5; ++field) {
        CAPTURE(field);
        const auto changed = take(
            PgGraph::prepare({field == 0 ? zero : base.x0, field == 1 ? one : base.width,
                              field == 2 ? zero : base.skew, field == 3 ? changed_f : base.height,
                              field == 4 ? 2 : 1, PgGraphConvention::Version1},
                             limits(), plan));
        CHECK_FALSE(take(changed.view()).identity == base.identity);
        const auto rejected =
            changed.evaluate(point, {0, 0}, PgDerivativeFrame::ChartST, limits(), plan);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == core::FailureCode::StaleGeneration);
        const auto wrong_map = activity.validate_for(changed);
        REQUIRE_FALSE(wrong_map);
        CHECK(wrong_map.error().code == core::FailureCode::StaleGeneration);
    }
}
TEST_CASE("Graph outputs own their source and high derivatives need no giant scale",
          "[SH-D03][N4-02][N4-08][property][verification]") {
    const auto detached = []() {
        core::ResourcePlan input{budget};
        const auto g = graph(4, input);
        const ExactRational zero;
        return take(
            g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, limits(), input));
    }();
    core::ResourcePlan plan{budget};
    const auto g = graph(4, plan);
    const auto value =
        take(g.evaluate(detached, {0, 0}, PgDerivativeFrame::ChartST, limits(), plan));
    position_is(value, ref(1, 16), ref(1, 16));
    auto restricted = limits();
    restricted.algebra.maximum_work_bits = 0;
    restricted.algebra.external_memory = ExactEstimateAdmission{0};
    const auto derivative = take(g.evaluate(
        detached,
        {std::numeric_limits<core::Index>::max(), std::numeric_limits<core::Index>::max()},
        PgDerivativeFrame::CellUV, restricted, plan));
    for (const auto& item : take(derivative.view()).components)
        fraction(item, "0", "1");
    CHECK(derivative.preparation().arithmetic.algebraic_work_bits == 0);
    CHECK(derivative.preparation().arithmetic.planned_work_bits == 0);
    const auto flat = flat_graph(4, plan);
    const auto t = r("1", "3", plan);
    const auto flat_point =
        take(flat.locate(t, t, PgSide::Interior, PgSide::Interior, limits(), plan));
    const auto zero_jet =
        take(flat.evaluate(flat_point, {2, 0}, PgDerivativeFrame::CellUV, restricted, plan));
    for (const auto& item : take(zero_jet.view()).components)
        fraction(item, "0", "1");
    const auto old = take(g.view()).identity;
    auto moved = take(g.clone(limits(), plan));
    auto owner = std::move(moved);
    // NOLINTNEXTLINE(bugprone-use-after-move): checks the required moved-from error.
    CHECK_FALSE(moved.view());
    CHECK(take(owner.view()).identity == old);
    std::array<bool, 2> answers{};
    std::array<std::jthread, 2> workers;
    for (std::size_t i = 0; i < workers.size(); ++i)
        workers[i] = std::jthread([&, i]() {
            try {
                core::ResourcePlan local{budget};
                const auto clone = take(g.clone(limits(), local));
                answers[i] = take(clone.view()).identity == old;
            } catch (...) {
                answers[i] = false;
            }
        });
    for (auto& worker : workers)
        worker.join();
    CHECK(answers[0]);
    CHECK(answers[1]);
}
TEST_CASE("Graph resource preflight and invalid requests preserve input and old results",
          "[SH-D03][N4-02][N4-08][negative][verification]") {
    core::ResourcePlan input{budget};
    const auto g = graph(1, input);
    const auto definition = take(g.view());
    const ExactRational zero;
    const auto one = r("1", "1", input);
    for (const core::Index n : {0, 3, -1, 8}) {
        core::ResourcePlan output{0};
        const auto bad = PgGraph::prepare({definition.x0, definition.width, definition.skew,
                                           definition.height, n, PgGraphConvention::Version1},
                                          limits(), output);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::InvalidInput);
        CHECK(output.allocated_bytes() == 0);
    }
    {
        core::ResourcePlan output{0};
        const auto bad = PgGraph::prepare({definition.x0, zero, definition.skew, definition.height,
                                           1, PgGraphConvention::Version1},
                                          limits(), output);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::InvalidInput);
        CHECK(output.allocated_bytes() == 0);
    }
    {
        core::ResourcePlan output{0};
        const auto bad =
            PgGraph::prepare({definition.x0, definition.width, definition.skew, definition.height,
                              1, static_cast<PgGraphConvention>(255)},
                             limits(), output);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::InvalidInput);
    }
    for (const auto side : {PgSide::Interior, PgSide::Lower, static_cast<PgSide>(255)}) {
        core::ResourcePlan output{0};
        const auto bad = g.locate(zero, zero, side, PgSide::Upper, limits(), output);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::InvalidInput);
        CHECK(output.allocated_bytes() == 0);
    }
    {
        core::ResourcePlan output{0};
        const auto bad = g.locate(one, one, PgSide::Upper, PgSide::Lower, limits(), output);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::InvalidInput);
    }
    const auto probe = take(g.clone(limits(), input));
    const auto info = probe.preparation();
    auto exact = limits();
    exact.algebra.external_memory =
        ExactEstimateAdmission{info.arithmetic.estimated_external_payload_bytes};
    exact.algebra.maximum_work_bits = info.arithmetic.planned_work_bits;
    {
        core::ResourcePlan output{info.planned_owned_payload_upper_bytes};
        const auto ok = g.clone(exact, output);
        REQUIRE(ok);
        CHECK(ok->preparation().allocated_owned_payload_bytes == output.allocated_bytes());
    }
    for (unsigned field = 0; field < 5; ++field) {
        core::ResourcePlan output{field == 0 ? info.planned_owned_payload_upper_bytes - 1 : budget};
        auto insufficient = exact;
        if (field == 1)
            insufficient.algebra.external_memory =
                ExactEstimateAdmission{info.arithmetic.estimated_external_payload_bytes - 1};
        if (field == 2)
            insufficient.algebra.maximum_work_bits = info.arithmetic.planned_work_bits - 1;
        if (field == 3)
            insufficient.maximum_cells = 3;
        if (field == 4)
            insufficient.algebra.maximum_raw_terms = 3;
        const auto bad = g.clone(insufficient, output);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::ResourceExceeded);
        CHECK(output.allocated_bytes() == 0);
    }
    auto strict = limits();
    strict.algebra.external_memory = ExactRequireProvenBound{budget};
    core::ResourcePlan output{budget};
    const auto unproved = g.clone(strict, output);
    REQUIRE_FALSE(unproved);
    CHECK(unproved.error().code == core::FailureCode::Unsupported);
    CHECK(output.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop is non-const; libstdc++
    // declares it const.
    std::stop_source cancellation;
    CHECK(cancellation.request_stop());
    auto stopped = limits();
    stopped.cancellation = cancellation.get_token();
    const auto cancelled = g.clone(stopped, output);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK(output.allocated_bytes() == 0);
    CHECK(take(g.view()).identity == definition.identity);
}
TEST_CASE("A full quartic exercises the 129-raw-term cumulative allocation envelope",
          "[SH-D03][N4-02][N4-08][verification]") {
    core::ResourcePlan input{budget};
    const auto one = r("1", "1", input), x = r("1", "7", input), w = r("4", "5", input),
               a = r("3", "20", input);
    std::array<ExactInputTerm, 15> terms{{{{0, 0}, one},
                                          {{0, 1}, one},
                                          {{0, 2}, one},
                                          {{0, 3}, one},
                                          {{0, 4}, one},
                                          {{1, 0}, one},
                                          {{1, 1}, one},
                                          {{1, 2}, one},
                                          {{1, 3}, one},
                                          {{2, 0}, one},
                                          {{2, 1}, one},
                                          {{2, 2}, one},
                                          {{3, 0}, one},
                                          {{3, 1}, one},
                                          {{4, 0}, one}}};
    const auto f = take(ExactChartPolynomial::make(terms, limits().algebra, input));
    const auto g =
        take(PgGraph::prepare({x, w, a, f, 4, PgGraphConvention::Version1}, limits(), input));
    core::ResourcePlan output{budget};
    const auto values = take(g.cell_polynomials({0, 0}, limits(), output));
    const auto info = values.preparation();
    CHECK(info.arithmetic.planned_raw_terms == 129);
    CHECK(info.allocated_owned_payload_bytes == info.planned_owned_payload_upper_bytes);
    CHECK(output.allocated_bytes() == info.allocated_owned_payload_bytes);
    auto cap = limits();
    cap.algebra.maximum_raw_terms = 128;
    core::ResourcePlan none{budget};
    const auto bad = g.cell_polynomials({0, 0}, cap, none);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().code == core::FailureCode::ResourceExceeded);
    CHECK(none.allocated_bytes() == 0);
}

TEST_CASE("Known graph-call resource and shape failures occur before payload allocation",
          "[SH-D03][N4-02][N4-08][negative][verification]") {
    core::ResourcePlan input{budget};
    const auto g = flat_graph(4, input);
    const ExactRational zero;
    const auto point =
        take(g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, limits(), input));
    CHECK(point.preparation().arithmetic.algebraic_work_bits == 9);
    for (unsigned field = 0; field < 4; ++field) {
        auto restricted = limits();
        if (field == 0)
            restricted.algebra.maximum_axis_degree = 0;
        if (field == 1)
            restricted.algebra.maximum_total_degree = 0;
        if (field == 2)
            restricted.algebra.maximum_terms = 1;
        if (field == 3)
            restricted.algebra.maximum_input_characters = 1;
        core::ResourcePlan output{budget};
        const auto bad = g.cell_polynomials({0, 0}, restricted, output);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::ResourceExceeded);
        CHECK(output.allocated_bytes() == 0);
    }
    auto short_native = limits();
    short_native.algebra.maximum_input_characters = 1;
    core::ResourcePlan output{budget};
    const auto invalid_native =
        g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, short_native, output);
    REQUIRE_FALSE(invalid_native);
    CHECK(invalid_native.error().code == core::FailureCode::ResourceExceeded);
    CHECK(output.allocated_bytes() == 0);
    for (const auto orders : {PgDerivativeOrders{-1, 0}, PgDerivativeOrders{0, -1}}) {
        const auto invalid =
            g.evaluate(point, orders, PgDerivativeFrame::ChartST, limits(), output);
        REQUIRE_FALSE(invalid);
        CHECK(invalid.error().code == core::FailureCode::InvalidInput);
        CHECK(output.allocated_bytes() == 0);
    }
    const auto bad_frame =
        g.evaluate(point, {0, 0}, static_cast<PgDerivativeFrame>(255), limits(), output);
    REQUIRE_FALSE(bad_frame);
    CHECK(bad_frame.error().code == core::FailureCode::InvalidInput);
    CHECK(output.allocated_bytes() == 0);
    const auto x = r("1" + std::string(200, '0'), "3", input), one = r("1", "1", input);
    const ExactChartPolynomial empty;
    const auto large = take(
        PgGraph::prepare({x, one, zero, empty, 1, PgGraphConvention::Version1}, limits(), input));
    const auto map = take(QualifiedRectangularActivityMap::prepare(large, limits(), input));
    auto no_external = limits();
    no_external.algebra.external_memory = ExactEstimateAdmission{0};
    const auto area = map.cell_area(no_external, output);
    REQUIRE_FALSE(area);
    CHECK(area.error().code == core::FailureCode::ResourceExceeded);
    CHECK(output.allocated_bytes() == 0);
    const auto total = map.total_area(no_external, output);
    REQUIRE_FALSE(total);
    CHECK(total.error().code == core::FailureCode::ResourceExceeded);
    CHECK(output.allocated_bytes() == 0);
    const auto before = input.allocated_bytes();
    const auto normal = take(large.clone(limits(), input));
    CHECK(normal.preparation().allocated_owned_payload_bytes == input.allocated_bytes() - before);
}
static_assert(
    !std::is_constructible_v<uvlas::structure::PgGraph,
                             uvlas::core::PlannedBuffer<uvlas::structure::detail::PgGraphData>,
                             uvlas::structure::PgGraphCallInfo>);
static_assert(!std::is_constructible_v<
              uvlas::structure::QualifiedRectangularActivityMap,
              uvlas::core::PlannedBuffer<uvlas::structure::detail::PgRectangleActivityData>,
              uvlas::structure::PgGraphCallInfo>);

TEST_CASE("Every owned graph operation preserves an explicit cancellation request",
          "[SH-D03][N4-08][negative]") {
    core::ResourcePlan input{budget};
    const auto g = graph(1, input);
    const auto definition = take(g.view());
    const auto t = r("1", "3", input);
    const ExactRational zero;
    const auto p =
        take(g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, limits(), input));
    const auto activity = take(QualifiedRectangularActivityMap::prepare(g, limits(), input));
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop is non-const.
    std::stop_source cancellation;
    CHECK(cancellation.request_stop());
    auto stopped = limits();
    stopped.cancellation = cancellation.get_token();
    core::ResourcePlan output{0};
    const auto cancelled = [&](const auto& result) {
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::Cancelled);
        CHECK(result.error().operation == core::Operation::PrepareModel);
        CHECK(output.allocated_bytes() == 0);
    };
    cancelled(PgGraph::prepare({definition.x0, definition.width, definition.skew, definition.height,
                                1, PgGraphConvention::Version1},
                               stopped, output));
    cancelled(g.clone(stopped, output));
    cancelled(g.locate(t, t, PgSide::Interior, PgSide::Interior, stopped, output));
    cancelled(
        g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, stopped, output));
    cancelled(g.cell_polynomials({0, 0}, stopped, output));
    cancelled(g.evaluate(p, {0, 0}, PgDerivativeFrame::ChartST, stopped, output));
    cancelled(g.anchor_at(t, PgSide::Interior, 0, stopped, output));
    cancelled(QualifiedRectangularActivityMap::prepare(g, stopped, output));
    cancelled(activity.cell_area(stopped, output));
    cancelled(activity.total_area(stopped, output));
    CHECK(take(g.view()).identity == definition.identity);
}
TEST_CASE("Optional graph preparation clocks preserve exact values and identities",
          "[SH-D03][N4-02][property]") {
    core::ResourcePlan plan{budget};
    const auto g = graph(1, plan);
    auto timed = limits();
    timed.algebra.timing = ExactTimingMode::Enabled;
    const auto on = take(g.clone(timed, plan));
    CHECK(take(on.view()).identity == take(g.view()).identity);
    const auto& graph_timing = on.preparation().arithmetic.timing;
    if (!graph_timing) {
        FAIL_CHECK("Missing graph timing");
        return;
    }
    CHECK(graph_timing->calls == 1);
    CHECK_FALSE(g.preparation().arithmetic.timing);
    const ExactRational zero;
    const auto p =
        take(g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, limits(), plan));
    const auto tp =
        take(g.from_local({0, 0}, zero, zero, PgSide::Interior, PgSide::Interior, timed, plan));
    const auto a = take(p.view()), b = take(tp.view());
    CHECK(a.s == b.s);
    CHECK(a.t == b.t);
    CHECK(a.u == b.u);
    CHECK(a.v == b.v);
    CHECK(a.identity == b.identity);
    const auto value = take(g.evaluate(p, {2, 1}, PgDerivativeFrame::CellUV, limits(), plan));
    const auto measured = take(g.evaluate(tp, {2, 1}, PgDerivativeFrame::CellUV, timed, plan));
    const auto original = take(value.view()), observed = take(measured.view());
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(original.components[i] == observed.components[i]);
    CHECK(original.identity == observed.identity);
    CHECK_FALSE(value.preparation().arithmetic.timing);
    const auto& timing = measured.preparation().arithmetic.timing;
    if (!timing) {
        FAIL_CHECK("Missing evaluation timing");
        return;
    }
    CHECK(timing->calls == 1);
}
