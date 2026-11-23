#include "pg_partition_test_support.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <future>
#include <stop_token>
#include <string>
#include <string_view>
#include <utility>
#include <uvlas/structure/pg_rank_preparation.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using partition_test::Fixture;
using partition_test::take;
constexpr std::uint64_t gib = std::uint64_t{1} << 30U;
constexpr auto rank_steps = std::uint64_t{10000000000000};

struct Producers {
    core::ResourcePlan plan{48 * gib};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    PgCoefficientDictionary dictionary;
    PgIntegrationPartition partition;

    explicit Producers(Fixture& f)
        : dictionary(take(PgCoefficientDictionary::prepare(
              *f.graph, {partition_test::exact_limits(), 2000, 20000000, 100000, {}}, workspace,
              plan))),
          partition(take(PgIntegrationPartition::prepare(f.request(), partition_test::limits(),
                                                         workspace, plan))) {}
};

[[nodiscard]] PgRankPreparationLimits limits() {
    return {48 * gib, rank_steps, 128, {}};
}

[[nodiscard]] PgRankPreparation make_owned_flat() {
    Fixture f;
    Producers inputs(f);
    return take(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, limits(),
                                           inputs.workspace, inputs.plan));
}

template <class T>
void check_failure(const core::Result<T>& value, core::FailureCode code) {
    if (value)
        FAIL("The rank operation accepted a negative control");
    else
        CHECK(value.error().code == code);
}

void make_graph(Fixture& f, const ExactRational& x0, const ExactChartPolynomial& height) {
    const auto one = f.r("1");
    const auto zero = f.r("0");
    const auto algebra = partition_test::exact_limits();
    f.graph = take(PgGraph::prepare({x0, one, zero, height, 1, PgGraphConvention::Version1},
                                    {algebra, 64, {}}, f.input_plan));
}
} // namespace

TEST_CASE("PG local rank owns the full unconditioned basis and source maps",
          "[SH-D03][SH-D05][SH-D07][N4-02][verification]") {
    const auto rank = make_owned_flat();
    const auto info = take(rank.preparation());
    const auto selected = take(rank.selected());
    const auto rows = take(rank.rows());
    REQUIRE(info.rank == 109);
    CHECK(info.rank_c == 4);
    CHECK(info.prime == 2147483647);
    CHECK(info.proof_rows == info.active_rows * info.candidates);
    CHECK(info.verified_anchor_rows > 0);
    CHECK(info.allocated_payload_bytes <= info.preflight_bytes);
    CHECK(rows.size() == info.active_rows);
    CHECK(selected.size() == info.rank);

    Fixture literal;
    const auto one = literal.r("1");
    const auto zero = literal.r("0");
    core::ResourcePlan reads{gib};
    const auto algebra = partition_test::exact_limits();
    for (std::size_t k = 0; k < selected.size(); ++k) {
        const auto basis = PgRankBasisColumn{static_cast<core::Index>(k)};
        if (k != 0)
            REQUIRE(selected[k - 1] < selected[k]);
        const auto j = static_cast<core::Index>(selected[k]);
        if (j < 33) {
            CHECK(take(rank.required(basis, {j}, algebra, reads)) == one);
        } else {
            const auto source = PgRankOriginalColumn{j - 33};
            CHECK(take(rank.original(basis, source, algebra, reads)) == one);
            CHECK(take(rank.source_map(basis, source, algebra, reads)) == one);
        }
        const auto neighbor =
            PgRankBasisColumn{static_cast<core::Index>((k + 1) % selected.size())};
        if (j < 33)
            CHECK(take(rank.required(neighbor, {j}, algebra, reads)) == zero);
        else
            CHECK(take(rank.original(neighbor, {j - 33}, algebra, reads)) == zero);
    }
    REQUIRE(selected[0] == 0);
    const auto quarter = literal.r("1", "4");
    std::uint32_t own_terms = 0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& key = rows[i];
        if (key.field != PgCoefficientField::R || key.cell != PgCellIndex{0, 0})
            continue;
        const bool nonzero =
            key.powers.second == 0 && (key.powers.first == 0 || key.powers.first == 1);
        const auto value = take(rank.basis({0}, static_cast<core::Index>(i), algebra, reads));
        CHECK(value == (nonzero ? quarter : zero));
        if (nonzero) {
            ++own_terms;
            core::ResourcePlan no_storage{0};
            check_failure(rank.basis({0}, static_cast<core::Index>(i), algebra, no_storage),
                          core::FailureCode::ResourceExceeded);
            CHECK(no_storage.allocated_bytes() == 0);
        }
    }
    CHECK(own_terms == 2);

    const auto witness = take(rank.witness());
    CHECK(witness.n == 1);
    CHECK(witness.height_terms == 0);
    CHECK(witness.domain_vertices == 4);
    CHECK(witness.activity_cells == 4);
    const auto chart = take(rank.source_chart(algebra, reads));
    CHECK(chart.x0 == zero);
    CHECK(chart.width == one);
    CHECK(chart.skew == zero);
    const auto first = take(rank.source_domain(0, algebra, reads));
    const auto third = take(rank.source_domain(2, algebra, reads));
    CHECK(first.s == zero);
    CHECK(first.t == zero);
    CHECK(third.s == one);
    CHECK(third.t == one);
    for (core::Index i = 0; i < 4; ++i) {
        const auto cell = take(rank.source_activity(i, algebra, reads));
        CHECK(cell.cell == PgCellIndex{i / 2, i % 2});
        CHECK(cell.area == quarter);
        CHECK(cell.active);
    }
    check_failure(rank.source_domain(4, algebra, reads), core::FailureCode::InvalidInput);
    check_failure(rank.basis({0}, -1, algebra, reads), core::FailureCode::InvalidInput);
}

TEST_CASE("PG rank preparation rejects short budgets without changing sources",
          "[SH-D03][N4-08][negative][verification]") {
    Fixture f;
    Producers inputs(f);
    const auto input_identity = take(inputs.partition.identity());
    const auto initial = inputs.plan.allocated_bytes();
    const auto rank = take(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, limits(),
                                                      inputs.workspace, inputs.plan));
    const auto info = take(rank.preparation());
    CHECK(inputs.plan.allocated_bytes() - initial == info.allocated_payload_bytes);

    auto cap = limits();
    cap.maximum_payload_bytes = info.preflight_bytes - 1;
    core::ResourcePlan byte_short{48 * gib};
    check_failure(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, cap,
                                             inputs.workspace, byte_short),
                  core::FailureCode::ResourceExceeded);
    CHECK(byte_short.allocated_bytes() == 0);
    cap.maximum_payload_bytes = info.preflight_bytes;
    core::ResourcePlan plan_short{info.preflight_bytes - 1};
    check_failure(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, cap,
                                             inputs.workspace, plan_short),
                  core::FailureCode::ResourceExceeded);
    CHECK(plan_short.allocated_bytes() == 0);

    cap.maximum_steps = info.steps;
    core::ResourcePlan equal{info.preflight_bytes};
    const auto retry = take(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, cap,
                                                       inputs.workspace, equal));
    CHECK(take(retry.preparation()).steps == info.steps);
    CHECK(take(retry.preparation()).rank == info.rank);
    --cap.maximum_steps;
    core::ResourcePlan work_short{info.preflight_bytes};
    check_failure(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, cap,
                                             inputs.workspace, work_short),
                  core::FailureCode::ResourceExceeded);
    CHECK(take(inputs.partition.identity()) == input_identity);

    std::stop_source source;
    REQUIRE(source.request_stop());
    cap = limits();
    cap.cancellation = source.get_token();
    core::ResourcePlan cancelled{48 * gib};
    check_failure(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, cap,
                                             inputs.workspace, cancelled),
                  core::FailureCode::Cancelled);
    CHECK(cancelled.allocated_bytes() == 0);
    const auto lease = take(inputs.workspace.acquire());
    check_failure(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, limits(),
                                             inputs.workspace, cancelled),
                  core::FailureCode::InvalidInput);
    CHECK(cancelled.allocated_bytes() == 0);
    static_cast<void>(lease);
}

TEST_CASE("PG bounded rank retries a prime dividing a real input denominator",
          "[SH-D03][N4-02][negative][verification]") {
    Fixture f;
    const auto x0 = f.r("1", "2147483647");
    const ExactChartPolynomial flat;
    make_graph(f, x0, flat);
    Producers inputs(f);
    const auto rank = take(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, limits(),
                                                      inputs.workspace, inputs.plan));
    const auto info = take(rank.preparation());
    CHECK(info.prime == 2147483629);
    CHECK(info.rank == 109);
    CHECK(info.proof_rows == info.active_rows * info.candidates);
}

TEST_CASE("PG bounded rank retains tiny nonzero curvature and rejects insufficient lifting",
          "[SH-D03][N4-02][N4-08][negative][verification]") {
    const auto rank_trial = [](std::string_view denominator, std::uint32_t digits) {
        Fixture f;
        const auto origin = f.r("0");
        const auto value = f.r("1", denominator);
        const std::array<ExactInputTerm, 1> terms{{{{2, 1}, value}}};
        const auto algebra = partition_test::exact_limits();
        const auto height = take(ExactChartPolynomial::make(terms, algebra, f.input_plan));
        make_graph(f, origin, height);
        Producers inputs(f);
        auto cap = limits();
        cap.maximum_digits = digits;
        return PgRankPreparation::prepare(inputs.dictionary, inputs.partition, cap,
                                          inputs.workspace, inputs.plan);
    };
    constexpr std::string_view nonzero = "1267650600228229401496703205376"; // 2^100
    const auto tiny = take(rank_trial(nonzero, 128));
    const auto base = take(rank_trial("1", 128));
    CHECK(take(tiny.preparation()).rank == take(base.preparation()).rank);
    CHECK(take(tiny.preparation()).rank > 109);
    const auto too_wide = rank_trial(std::string("1") + std::string(600, '0'), 128);
    check_failure(too_wide, core::FailureCode::ResourceExceeded);
    const auto short_run = rank_trial(nonzero, 1);
    check_failure(short_run, core::FailureCode::RankUnresolved);
}

TEST_CASE("Independent PG rank workspaces preserve the same canonical exact result",
          "[SH-D03][N4-08][verification]") {
    Fixture f;
    const Producers inputs(f);
    const auto run = [&] {
        core::ResourcePlan plan{48 * gib};
        core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
        return take(PgRankPreparation::prepare(inputs.dictionary, inputs.partition, limits(),
                                               workspace, plan));
    };
    auto first = std::async(std::launch::async, run);
    auto second = std::async(std::launch::async, run);
    const auto a = first.get();
    const auto b = second.get();
    const auto info = take(a.preparation());
    CHECK(take(b.preparation()).steps == info.steps);
    CHECK(take(b.preparation()).allocated_payload_bytes == info.allocated_payload_bytes);
    REQUIRE(std::ranges::equal(take(a.selected()), take(b.selected())));
    core::ResourcePlan query{gib};
    const auto algebra = partition_test::exact_limits();
    const auto originals = take(inputs.dictionary.original_columns());
    const auto columns = static_cast<core::Index>(take(inputs.dictionary.columns()).size());
    const auto compare = [&](auto read, auto first_column, core::Index count) {
        for (core::Index j = 0; j < count; ++j)
            for (core::Index k = 0; k < static_cast<core::Index>(info.rank); ++k) {
                const auto basis = PgRankBasisColumn{k};
                first_column.value = j;
                CHECK(take((a.*read)(basis, first_column, algebra, query)) ==
                      take((b.*read)(basis, first_column, algebra, query)));
            }
    };
    compare(&PgRankPreparation::source_map, PgRankOriginalColumn{}, columns);
    compare(&PgRankPreparation::original, PgRankOriginalColumn{}, columns);
    compare(&PgRankPreparation::required, PgRankRequiredColumn{}, 33);
    compare(&PgRankPreparation::anchor, PgRankAnchorColumn{}, originals);
    for (core::Index k = 0; k < static_cast<core::Index>(info.rank); ++k)
        for (core::Index i = 0; i < static_cast<core::Index>(info.active_rows); ++i)
            CHECK(take(a.basis({k}, i, algebra, query)) == take(b.basis({k}, i, algebra, query)));
}
