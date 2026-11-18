#include "pg_coefficients_plan.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdint>
#include <future>
#include <ios>
#include <stop_token>
#include <string_view>
#include <type_traits>
#include <utility>
#include <uvlas/structure/pg_coefficients.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace uvlas::structure::detail;
constexpr auto payload_budget = std::uint64_t{1024} * 1024 * 1024;
PgCoefficientLimits coefficient_limits() {
    return {{64, 1000, 100000, 32, 32, 65536, 10000000, ExactEstimateAdmission{payload_budget},
             ExactTimingMode::Disabled},
            2000,
            20000000,
            100000,
            {}};
}
void fraction(const ExactRational& value, std::string_view numerator,
              std::string_view denominator) {
    const auto& b = ExactFactory::value(value).backend();
    CHECK(b.num().str(0, std::ios_base::dec) == numerator);
    CHECK(b.denom().str(0, std::ios_base::dec) == denominator);
}
} // namespace
TEST_CASE("PG quadratic and linear pieces satisfy independent Bernstein literals",
          "[SH-D03][SH-D07][N4-02][verification]") {
    core::ResourcePlan plan{payload_budget};
    const auto caps = coefficient_limits();
    PgCoefficientAlgebra op{caps, plan};
    const auto zero = op.number(0);
    const auto quadratic = pg_basis(op, 1, 0, ExactAxis::First, PgUnivariateKind::QuadraticC1);
    REQUIRE(quadratic.count == 3);
    fraction(op.e(quadratic.values[0], zero, zero), "9", "16");
    fraction(op.e(quadratic.values[1], zero, zero), "3", "8");
    fraction(op.e(quadratic.values[2], zero, zero), "1", "16");
    const auto hat = pg_basis(op, 1, 0, ExactAxis::First, PgUnivariateKind::LinearC0);
    fraction(op.e(hat.values[0], zero, zero), "3", "4");
    fraction(op.e(hat.values[1], zero, zero), "1", "4");
    const auto unit = op.constant(1);
    for (const core::Index n : {1, 2, 4})
        for (core::Index cell = 0; cell < 2 * n; ++cell)
            for (const auto axis : {ExactAxis::First, ExactAxis::Second})
                for (const auto kind :
                     {PgUnivariateKind::QuadraticC1, PgUnivariateKind::QuadraticC0,
                      PgUnivariateKind::LinearC0}) {
                    CAPTURE(n, cell, axis, kind);
                    const auto basis = pg_basis(op, n, cell, axis, kind);
                    PgPoly sum;
                    for (core::Index j = 0; j < basis.count; ++j) {
                        const auto& poly = basis.values[static_cast<std::size_t>(j)];
                        CHECK(op.e(poly, zero, zero).sign() >= 0);
                        sum = op.a(sum, poly);
                    }
                    CHECK(pg_equal(sum, unit));
                }
}
TEST_CASE("PG flat lifted fields retain full moment carry and physical-plane signs",
          "[SH-D03][SH-D05][N4-02][verification]") {
    core::ResourcePlan plan{payload_budget};
    const auto caps = coefficient_limits();
    PgCoefficientAlgebra op{caps, plan};
    const auto one = op.number(1), quarter = op.number(1, 4), three_quarters = op.number(3, 4),
               width_inverse = op.number(5, 4), skew = op.number(3, 20), scale = op.number(4);
    const auto legendre = pg_legendre(op, ExactAxis::First);
    REQUIRE(legendre[7].terms().size() == 4);
    fraction(legendre[7].terms()[0].coefficient, "-35", "16");
    fraction(legendre[7].terms()[3].coefficient, "429", "16");
    const auto moments = pg_moments(op, 1, 0, legendre[0], ExactAxis::First);
    fraction(op.e(moments.active[0], op.zero, op.zero), "1", "4");
    fraction(op.e(moments.active[1], op.zero, op.zero), "1", "32");
    fraction(op.e(moments.completed[0], op.zero, op.zero), "1", "2");
    fraction(op.e(moments.completed[1], op.zero, op.zero), "1", "8");
    const ExactChartPolynomial flat;
    const PgPoly zero;
    const auto s = op.affine(quarter, quarter, ExactAxis::First);
    const auto t = op.affine(quarter, quarter, ExactAxis::Second);
    const auto transverse = pg_powers(op, t);
    const auto v = pg_lifting(op, flat, zero, s, transverse, moments.active, PgPairedAxis::S,
                              width_inverse, skew);
    const auto unit = op.constant(one);
    const auto p = pg_director(op, v, zero, zero, unit, unit, zero, scale, width_inverse, skew);
    fraction(op.e(v[2], op.zero, op.zero), "1", "32");
    fraction(op.e(p[0], op.zero, op.zero), "-5", "16");
    fraction(op.e(p[1], op.zero, op.zero), "3", "64");
    CHECK(p[2].terms().empty());
    const auto tail_s = op.affine(three_quarters, quarter, ExactAxis::First);
    const auto tail = pg_lifting(op, flat, zero, tail_s, transverse, moments.completed,
                                 PgPairedAxis::S, width_inverse, skew);
    fraction(op.e(tail[2], op.zero, op.zero), "1", "4");
    const auto t_legendre = pg_legendre(op, ExactAxis::Second);
    const auto t_moments = pg_moments(op, 1, 0, t_legendre[0], ExactAxis::Second);
    const auto s_powers = pg_powers(op, s);
    const auto vt = pg_lifting(op, flat, zero, t, s_powers, t_moments.active, PgPairedAxis::T,
                               width_inverse, skew);
    const auto pt = pg_director(op, vt, zero, zero, unit, unit, zero, scale, width_inverse, skew);
    CHECK(pt[0].terms().empty());
    fraction(op.e(pt[1], op.zero, op.zero), "-1", "4");
}

namespace {
constexpr auto whole_budget = std::uint64_t{64} * 1024 * 1024 * 1024;
PgCoefficientLimits whole_limits() {
    return {{64, 1000, 4096, 32, 32, 65536, 2000000,
             ExactEstimateAdmission{std::uint64_t{16} * 1024 * 1024 * 1024},
             ExactTimingMode::Disabled},
            1000,
            20000000,
            2000,
            {}};
}
PgGraph flat_graph(core::ResourcePlan& plan, core::Index n = 1) {
    const auto limits = whole_limits();
    PgCoefficientAlgebra op{limits, plan};
    const auto x0 = op.number(2, 7), w = op.number(4, 5), a = op.number(3, 20);
    const ExactChartPolynomial flat;
    return pg_take(PgGraph::prepare({x0, w, a, flat, n, PgGraphConvention::Version1},
                                    {op.limits, 64, {}}, plan));
}
} // namespace
TEST_CASE("PG full preparation admits before allocating and owns its graph and exact values",
          "[SH-D03][SH-D05][SH-D07][N4-08][negative][property]") {
    core::ResourcePlan input{whole_budget};
    auto source = flat_graph(input);
    const auto original = pg_take(source.view()).identity;
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto limits = whole_limits();
    core::ResourcePlan pure{whole_budget};
    const auto forecast = pg_coefficient_plan(pg_take(source.view()), limits, pure);
    CHECK(pure.allocated_bytes() == 0);
    // Equality is accepted; each independent decrement is rejected before any payload.
    limits.algebra.maximum_work_bits = forecast.report.arithmetic.planned_work_bits;
    limits.algebra.maximum_raw_terms = forecast.report.arithmetic.planned_raw_terms;
    limits.algebra.external_memory =
        ExactEstimateAdmission{forecast.report.arithmetic.estimated_external_payload_bytes};
    limits.maximum_columns = forecast.full_columns;
    limits.maximum_map_entries = forecast.map_capacity + forecast.c_map_capacity;
    for (unsigned fault = 0; fault < 9; ++fault) {
        auto changed = limits;
        core::ResourcePlan output{fault == 0 ? forecast.report.planned_owned_payload_upper_bytes - 1
                                             : whole_budget};
        switch (fault) {
        case 0:
            break;
        case 1:
            --changed.algebra.maximum_work_bits;
            break;
        case 2:
            --changed.algebra.maximum_raw_terms;
            break;
        case 3:
            --changed.maximum_columns;
            break;
        case 4:
            --changed.maximum_map_entries;
            break;
        case 5:
            changed.maximum_stored_terms = 0;
            break;
        case 6:
            changed.algebra.external_memory = ExactEstimateAdmission{
                forecast.report.arithmetic.estimated_external_payload_bytes - 1};
            break;
        case 7:
            changed.algebra.external_memory = ExactRequireProvenBound{whole_budget};
            break;
        case 8:
            changed.algebra.timing = static_cast<ExactTimingMode>(255);
            break;
        default:
            FAIL("Unexpected fault selector");
        }
        const auto result = PgCoefficientDictionary::prepare(source, changed, workspace, output);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == (fault == 7   ? core::FailureCode::Unsupported
                                      : fault == 8 ? core::FailureCode::InvalidInput
                                                   : core::FailureCode::ResourceExceeded));
        CHECK(output.allocated_bytes() == 0);
        CHECK(pg_take(source.view()).identity == original);
    }
    {
        const auto lease = pg_take(workspace.acquire());
        core::ResourcePlan output{whole_budget};
        const auto busy = PgCoefficientDictionary::prepare(source, limits, workspace, output);
        REQUIRE_FALSE(busy);
        CHECK(busy.error().code == core::FailureCode::InvalidInput);
        CHECK(busy.error().operation == core::Operation::PrepareModel);
        CHECK(output.allocated_bytes() == 0);
    }
    // MSVC 14.51 request_stop is non-const; libstdc++ declares it const.
    // NOLINTNEXTLINE(misc-const-correctness): preserve the portable non-const handle.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    auto cancelled_limits = limits;
    cancelled_limits.cancellation = stop.get_token();
    core::ResourcePlan cancelled_plan{whole_budget};
    const auto cancelled =
        PgCoefficientDictionary::prepare(source, cancelled_limits, workspace, cancelled_plan);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK(cancelled_plan.allocated_bytes() == 0);
    core::ResourcePlan output{forecast.report.planned_owned_payload_upper_bytes};
    auto value = pg_take(PgCoefficientDictionary::prepare(source, limits, workspace, output));
    CHECK(value.preparation().observed_columns == 161);
    CHECK(value.preparation().allocated_owned_payload_bytes <=
          forecast.report.planned_owned_payload_upper_bytes);
    CHECK_FALSE(value.preparation().timings);
    CHECK(value.validate_for(source));
    const auto& coefficient = *pg_take(value.coefficient({0}, {0, 0}, PgCoefficientField::R));
    const auto& before = coefficient.terms()[0].coefficient;
    // A different actual snapshot cannot borrow this dictionary's qualification.
    auto different = flat_graph(input, 2);
    CHECK_FALSE(value.validate_for(different));
    source = std::move(different);
    CHECK(pg_take(value.graph()).identity == original);
    CHECK(coefficient.terms()[0].coefficient == before);
    auto moved = std::move(value);
    // NOLINTBEGIN(bugprone-use-after-move): the public contract explicitly rejects every moved-from
    // numerical query.
    CHECK_FALSE(value.columns());
    CHECK_FALSE(value.records());
    CHECK_FALSE(value.required_map());
    CHECK_FALSE(value.anchor_c_map());
    CHECK_FALSE(value.original_columns());
    CHECK_FALSE(value.coefficient({0}, {0, 0}, PgCoefficientField::R));
    // NOLINTEND(bugprone-use-after-move)
    CHECK(pg_take(moved.graph()).identity == original);
    CHECK_FALSE(moved.coefficient({-1}, {0, 0}, PgCoefficientField::R));
    CHECK_FALSE(moved.coefficient({0}, {2, 0}, PgCoefficientField::R));
    CHECK_FALSE(moved.coefficient({0}, {0, 0}, static_cast<PgCoefficientField>(255)));
    static_assert(!std::is_copy_constructible_v<PgCoefficientDictionary>);
}

TEST_CASE(
    "PG independent concurrent workspaces preserve canonical coefficients with timing enabled",
    "[SH-D03][N4-08][property]") {
    core::ResourcePlan inputs{whole_budget};
    const auto graph = flat_graph(inputs);
    const auto run = [&](ExactTimingMode timing) {
        auto limits = whole_limits();
        limits.algebra.timing = timing;
        core::ResourcePlan output{whole_budget};
        core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
        return PgCoefficientDictionary::prepare(graph, limits, workspace, output);
    };
    auto pending = std::async(std::launch::async, run, ExactTimingMode::Enabled);
    const auto disabled = run(ExactTimingMode::Disabled);
    const auto enabled = pending.get();
    REQUIRE(disabled);
    REQUIRE(enabled);
    const auto& timings = enabled->preparation().timings;
    if (!timings.has_value()) {
        FAIL_CHECK("Enabled product timing must be present");
        return;
    }
    CHECK(timings->coefficients.calls == 1);
    CHECK(timings->anchor.calls == 1);
    CHECK(timings->verification.calls == 1);
    const auto first = pg_take(disabled->records()), second = pg_take(enabled->records());
    REQUIRE(first.size() == second.size());
    for (std::size_t j = 0; j < first.size(); ++j) {
        CHECK(first[j].column == second[j].column);
        CHECK(first[j].cell == second[j].cell);
        CHECK(first[j].field == second[j].field);
        CHECK(pg_equal(first[j].polynomial, second[j].polynomial));
    }
    const auto first_map = pg_take(disabled->required_map()),
               second_map = pg_take(enabled->required_map());
    REQUIRE(first_map.size() == second_map.size());
    for (std::size_t j = 0; j < first_map.size(); ++j) {
        CHECK(first_map[j].source == second_map[j].source);
        CHECK(first_map[j].target == second_map[j].target);
        CHECK(first_map[j].coefficient == second_map[j].coefficient);
    }
}

namespace {
using OracleRational = boost::multiprecision::cpp_rational;
OracleRational oracle_power(const OracleRational& x, unsigned p) {
    OracleRational r = 1;
    for (unsigned j = 0; j < p; ++j)
        r *= x;
    return r;
}
struct OracleJet {
    OracleRational value, ds, dt;
};
OracleJet operator-(const OracleJet& a, const OracleJet& b) {
    return {a.value - b.value, a.ds - b.ds, a.dt - b.dt};
}
OracleJet operator*(const OracleJet& a, const OracleJet& b) {
    return {a.value * b.value, a.ds * b.value + a.value * b.ds, a.dt * b.value + a.value * b.dt};
}
OracleJet oracle_scale(const OracleJet& a, const OracleRational& b) {
    return {a.value * b, a.ds * b, a.dt * b};
}
OracleJet oracle_power(const OracleJet& x, unsigned p) {
    OracleJet r{1, 0, 0};
    for (unsigned j = 0; j < p; ++j)
        r = r * x;
    return r;
}
OracleRational oracle_value(const ExactRational& value) {
    const auto& r = ExactFactory::value(value).backend();
    return OracleRational{r.num().str(0, std::ios_base::dec)} /
           OracleRational{r.denom().str(0, std::ios_base::dec)};
}
OracleRational oracle_at_midpoint(const PgPoly& value) {
    for (const auto& term : value.terms())
        if (term.powers.first == 0 && term.powers.second == 0)
            return oracle_value(term.coefficient);
    return 0;
}
OracleRational oracle_chart_derivative_at_midpoint(const PgPoly& value, ExactAxis axis,
                                                   core::Index n) {
    const ExactPowers selected = axis == ExactAxis::First ? ExactPowers{1, 0} : ExactPowers{0, 1};
    for (const auto& term : value.terms())
        if (term.powers.first == selected.first && term.powers.second == selected.second)
            return oracle_value(term.coefficient) * (4 * n);
    return 0;
}
// Explicit classical Legendre coefficients, independent of the product recurrence.
constexpr std::array<std::array<int, 8>, 8> oracle_legendre{{{{1, 0, 0, 0, 0, 0, 0, 0}},
                                                             {{0, 1, 0, 0, 0, 0, 0, 0}},
                                                             {{-1, 0, 3, 0, 0, 0, 0, 0}},
                                                             {{0, -3, 0, 5, 0, 0, 0, 0}},
                                                             {{3, 0, -30, 0, 35, 0, 0, 0}},
                                                             {{0, 15, 0, -70, 0, 63, 0, 0}},
                                                             {{-5, 0, 105, 0, -315, 0, 231, 0}},
                                                             {{0, -35, 0, 315, 0, -693, 0, 429}}}};
constexpr std::array<int, 8> oracle_legendre_den{1, 1, 2, 2, 8, 8, 16, 16};
unsigned oracle_choose(unsigned n, unsigned k) {
    unsigned r = 1;
    for (unsigned j = 1; j <= k; ++j)
        r = r * (n + 1 - j) / j;
    return r;
}
OracleRational oracle_moment(unsigned m, unsigned k, core::Index n, bool completed) {
    OracleRational sum = 0;
    const OracleRational upper = completed ? 1 : 0;
    // Origin is the first half-cell: sigma=(1+x)/(4n). Integrate in x
    // directly; no local antiderivative or product evaluation generates expected values.
    for (unsigned a = 0; a <= m; ++a)
        for (unsigned j = 0; j <= k; ++j) {
            const unsigned power = a + j + 1;
            const OracleRational coefficient = OracleRational{oracle_choose(m, a)} *
                                               oracle_legendre[k][j] / oracle_legendre_den[k];
            sum += coefficient *
                   (oracle_power(upper, power) - oracle_power(OracleRational{-1}, power)) / power;
        }
    return sum / oracle_power(OracleRational{4 * n}, m + 1);
}
} // namespace
TEST_CASE("PG all quartic monomials and both full P7 families satisfy independent exact jets",
          "[SH-D03][SH-D05][N4-02][verification]") {
    // Unit monomials isolate every adopted height coefficient. At each n both
    // active and completed intervals are checked; no random sample or fitted tolerance.
    for (const core::Index n : {1, 2, 4})
        for (const auto axis : {PgPairedAxis::S, PgPairedAxis::T})
            for (unsigned degree = 0; degree < 8; ++degree)
                for (core::Index along_cell = 0; along_cell < 2 * n; ++along_cell)
                    for (core::Index transverse_cell = 0; transverse_cell < 2 * n;
                         ++transverse_cell) {
                        const bool completed = along_cell != 0;
                        core::ResourcePlan plan{whole_budget};
                        const auto limits = whole_limits();
                        PgCoefficientAlgebra op{limits, plan};
                        const auto along_axis =
                            axis == PgPairedAxis::S ? ExactAxis::First : ExactAxis::Second;
                        const auto legs = pg_legendre(op, along_axis);
                        const auto row = pg_moments(op, n, 0, legs[degree], along_axis);
                        const auto& moments = completed ? row.completed : row.active;
                        const auto cs = axis == PgPairedAxis::S ? along_cell : transverse_cell;
                        const auto ct = axis == PgPairedAxis::T ? along_cell : transverse_cell;
                        const auto s = op.affine(op.number(2 * cs + 1, 4 * n), op.number(1, 4 * n),
                                                 ExactAxis::First);
                        const auto t = op.affine(op.number(2 * ct + 1, 4 * n), op.number(1, 4 * n),
                                                 ExactAxis::Second);
                        const auto sp = pg_powers(op, s), tp = pg_powers(op, t);
                        const auto inverse = op.number(5, 4), skew = op.number(3, 20),
                                   scale = op.number(4 * n);
                        const auto one = op.constant(1);
                        const OracleJet sj{OracleRational{2 * cs + 1} / (4 * n), 1, 0};
                        const OracleJet tj{OracleRational{2 * ct + 1} / (4 * n), 0, 1};
                        const auto& q = axis == PgPairedAxis::S ? sj : tj;
                        const auto& z = axis == PgPairedAxis::S ? tj : sj;
                        std::array<OracleJet, 5> integral;
                        for (unsigned m = 0; m < 5; ++m) {
                            const OracleRational derivative =
                                completed ? OracleRational{0}
                                          : oracle_power(q.value, m) * oracle_legendre[degree][0] /
                                                oracle_legendre_den[degree];
                            integral[m] = {oracle_moment(m, degree, n, completed),
                                           axis == PgPairedAxis::S ? derivative : OracleRational{0},
                                           axis == PgPairedAxis::T ? derivative
                                                                   : OracleRational{0}};
                        }
                        for (unsigned i = 0; i <= 4; ++i)
                            for (unsigned j = 0; j <= 4 - i; ++j) {
                                CAPTURE(n, axis, degree, cs, ct, i, j);
                                const auto coefficient = op.number(1);
                                const std::array input{ExactInputTerm{{i, j}, coefficient}};
                                const auto global =
                                    pg_take(ExactChartPolynomial::make(input, op.limits, plan));
                                const auto local = op.m(sp[i], tp[j]);
                                const auto fx =
                                    op.s(op.d(local, ExactAxis::First), op.number(4 * n * 5, 4));
                                const auto fy = op.sub(op.s(op.d(local, ExactAxis::Second), scale),
                                                       op.s(fx, skew));
                                const auto U = op.a(one, op.m(fx, fx)), V = op.a(one, op.m(fy, fy)),
                                           xy = op.m(fx, fy);
                                const auto& along = axis == PgPairedAxis::S ? s : t;
                                const auto& transverse = axis == PgPairedAxis::S ? tp : sp;
                                const auto v = pg_lifting(op, global, local, along, transverse,
                                                          moments, axis, inverse, skew);
                                const auto p =
                                    pg_director(op, v, fx, fy, U, V, xy, scale, inverse, skew);
                                const unsigned a = axis == PgPairedAxis::S ? i : j,
                                               b = axis == PgPairedAxis::S ? j : i;
                                const auto f = oracle_power(sj, i) * oracle_power(tj, j);
                                const auto J = oracle_power(z, b) * integral[a];
                                const auto Jz =
                                    b == 0 ? OracleJet{0, 0, 0}
                                           : oracle_scale(oracle_power(z, b - 1) * integral[a], b);
                                const auto Jqz =
                                    b == 0
                                        ? OracleJet{0, 0, 0}
                                        : oracle_scale(oracle_power(z, b - 1) * integral[a + 1], b);
                                auto vr = std::array{J - f * integral[0], Jqz - q * Jz,
                                                     q * integral[0] - integral[1]};
                                if (axis == PgPairedAxis::T)
                                    std::swap(vr[0], vr[1]);
                                vr[0] = oracle_scale(vr[0], OracleRational{5} / 4);
                                vr[1] = vr[1] - oracle_scale(vr[0], OracleRational{3} / 20);
                                const OracleRational fxr = f.ds * 5 / 4, fyr = f.dt - fxr * 3 / 20;
                                std::array<OracleRational, 3> dx, dy, actual_ds, actual_dt;
                                for (std::size_t c = 0; c < 3; ++c) {
                                    CHECK(oracle_at_midpoint(v[c]) == vr[c].value);
                                    actual_ds[c] = oracle_chart_derivative_at_midpoint(
                                        v[c], ExactAxis::First, n);
                                    actual_dt[c] = oracle_chart_derivative_at_midpoint(
                                        v[c], ExactAxis::Second, n);
                                    CHECK(actual_ds[c] == vr[c].ds);
                                    CHECK(actual_dt[c] == vr[c].dt);
                                    dx[c] = vr[c].ds * 5 / 4;
                                    dy[c] = vr[c].dt - dx[c] * 3 / 20;
                                }
                                const OracleRational tx = dx[2] - fxr * dx[0] - fyr * dx[1];
                                const OracleRational ty = dy[2] - fxr * dy[0] - fyr * dy[1];
                                const std::array<OracleRational, 3> pr{
                                    -(1 + fyr * fyr) * tx + fxr * fyr * ty,
                                    fxr * fyr * tx - (1 + fxr * fxr) * ty, -fxr * tx - fyr * ty};
                                std::array<OracleRational, 3> actual_p;
                                for (std::size_t c = 0; c < 3; ++c) {
                                    actual_p[c] = oracle_at_midpoint(p[c]);
                                    CHECK(actual_p[c] == pr[c]);
                                    CHECK(pg_equal(
                                        op.d(op.d(v[c], ExactAxis::First), ExactAxis::Second),
                                        op.d(op.d(v[c], ExactAxis::Second), ExactAxis::First)));
                                }
                                CHECK(oracle_at_midpoint(fx) == fxr);
                                CHECK(oracle_at_midpoint(fy) == fyr);
                                const OracleRational actual_tx =
                                    (actual_ds[2] - fxr * actual_ds[0] - fyr * actual_ds[1]) * 5 /
                                    4;
                                const OracleRational actual_ty = actual_dt[2] - fxr * actual_dt[0] -
                                                                 fyr * actual_dt[1] -
                                                                 actual_tx * 3 / 20;
                                const OracleRational g = 1 + fxr * fxr + fyr * fyr;
                                CHECK(actual_p[0] + fxr * actual_p[2] + g * actual_tx == 0);
                                CHECK(actual_p[1] + fyr * actual_p[2] + g * actual_ty == 0);
                                CHECK(-fxr * actual_p[0] - fyr * actual_p[1] + actual_p[2] == 0);
                                // Only the along-family membrane component and shear
                                // vanish. The remaining component is generally nonzero.
                                const OracleRational ess =
                                    OracleRational{4} * actual_ds[0] / 5 + f.ds * actual_ds[2];
                                const OracleRational ett = actual_dt[1] +
                                                           OracleRational{3} * actual_dt[0] / 20 +
                                                           f.dt * actual_dt[2];
                                const OracleRational shear =
                                    OracleRational{4} * actual_dt[0] / 5 + actual_ds[1] +
                                    OracleRational{3} * actual_ds[0] / 20 + f.ds * actual_dt[2] +
                                    f.dt * actual_ds[2];
                                CHECK(shear == 0);
                                CHECK((axis == PgPairedAxis::S ? ess : ett) == 0);
                                const OracleRational remaining =
                                    b < 2
                                        ? OracleRational{0}
                                        : OracleRational{b * (b - 1)} *
                                              oracle_power(z.value, b - 2) *
                                              (integral[a + 1].value - q.value * integral[a].value);
                                CHECK((axis == PgPairedAxis::S ? ett : ess) == remaining);
                                if (n == 1 && axis == PgPairedAxis::S && degree == 0 &&
                                    !completed && transverse_cell == 0 && i == 0 && j == 2)
                                    CHECK(remaining == OracleRational{-1} / 16);
                            }
                    }
}
