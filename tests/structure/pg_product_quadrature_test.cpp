#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <uvlas/structure/pg_discrete_gram.hpp>

using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;

TEST_CASE("PG product rules preserve the triangle Jacobian on cut curved volumes",
          "[SH-D03][SH-D07][N4-02][N4-08][verification][negative]") {
    Fixture fixture;
    fixture.domain.clear();
    fixture.thickness.clear();
    fixture.heights.clear();
    const auto add = [&](std::string_view s, std::string_view tn, std::string_view td,
                         std::string_view hn, std::string_view hd) {
        fixture.domain.push_back({fixture.r(s), fixture.r(tn, td)});
        fixture.thickness.push_back({fixture.r(s), fixture.r(tn, td)});
        fixture.heights.push_back(fixture.r(hn, hd));
    };
    add("0", "0", "1", "1", "1");
    add("1", "1", "4", "3", "2");
    add("1", "7", "8", "3", "2");
    add("0", "1", "1", "1", "1");

    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan setup{profile}, low_plan{profile}, high_plan{profile};
    const auto make_curve_basis = [&] {
        const auto zero = fixture.r("0"), one = fixture.r("1"), length = fixture.r("2");
        const auto quadratic = fixture.r("1", "16");
        const std::array<ExactInputTerm, 1> terms{{{{2, 0}, quadratic}}};
        const auto curve =
            take(ExactChartPolynomial::make(terms, exact_limits(), fixture.input_plan));
        fixture.graph =
            take(PgGraph::prepare({zero, one, zero, curve, 1, PgGraphConvention::Version1},
                                  {exact_limits(), 64, {}}, fixture.input_plan));
        return make_basis(fixture, length, scratch, setup);
    };
    auto low_basis = make_curve_basis();
    auto high_basis = make_curve_basis();
    const auto size = static_cast<std::size_t>(take(high_basis.size()));
    const auto& partition = take(high_basis.source()).partition;
    const auto volumes = take(partition.volumes()).size();
    REQUIRE(volumes > 4);
    const auto* rank = take(high_basis.rank());
    std::vector<Q> translation;
    translation.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        auto coefficient = take(
            rank->required({static_cast<core::Index>(i)}, {30}, query_limits().algebra, setup));
        translation.push_back(q(coefficient));
    }

    const auto zero = fixture.r("0"), half = fixture.r("1", "2"), one = fixture.r("1");
    const auto quarter = fixture.r("1", "4"), three_quarters = fixture.r("3", "4");
    const std::array depth{PgLineQuadratureNode{half, one}};
    const std::array low{PgLineQuadratureNode{quarter, half},
                         PgLineQuadratureNode{three_quarters, half}};
    const PgDiscreteGramGrid grid{260, 192};
    const PgDiscreteGramLimits limits{query_limits(), steps};
    const std::array wrong_moment{PgLineQuadratureNode{quarter, one}};
    const std::array duplicate{PgLineQuadratureNode{half, half}, PgLineQuadratureNode{half, half}};
    const std::array zero_weight{PgLineQuadratureNode{half, zero}};
    core::ResourcePlan rejected{profile};
    for (const auto& input : std::array{PgVolumeProductRule{wrong_moment, low, depth},
                                        PgVolumeProductRule{low, duplicate, depth},
                                        PgVolumeProductRule{low, low, zero_weight},
                                        PgVolumeProductRule{low, low, {}}}) {
        check_failure(PgDiscreteGram::prepare_product(std::move(low_basis), input, grid, limits,
                                                      scratch, rejected),
                      core::FailureCode::InvalidInput);
        CHECK(take(low_basis.size()) == size);
    }
    core::ResourcePlan no_storage{0};
    check_failure(PgDiscreteGram::prepare_product(std::move(low_basis), {low, low, depth}, grid,
                                                  limits, scratch, no_storage),
                  core::FailureCode::ResourceExceeded);
    CHECK(no_storage.allocated_bytes() == 0);
    auto low_gram = take(PgDiscreteGram::prepare_product(std::move(low_basis), {low, low, depth},
                                                         grid, limits, scratch, low_plan));
    REQUIRE(low_gram.rule().size() == 4 * volumes);
    CHECK(low_gram.info().points == low_gram.rule().size());
    const std::array<Q, 4> u{Q("1/4"), Q("1/4"), Q("3/4"), Q("3/4")};
    const std::array<Q, 4> v{Q("3/16"), Q("9/16"), Q("1/16"), Q("3/16")};
    for (std::size_t i = 0; i < volumes; ++i) {
        Q weight_sum;
        for (std::size_t j = 0; j < 4; ++j) {
            const auto& point = low_gram.rule()[4 * i + j];
            CHECK(point.volume == static_cast<core::Index>(i));
            CHECK(q(point.u) == u[j]);
            CHECK(q(point.v) == v[j]);
            CHECK(q(point.depth) == Q("1/2"));
            CHECK(q(point.weight) == (1 - u[j]) / 4);
            weight_sum += q(point.weight);
        }
        CHECK(weight_sum == Q("1/2"));
    }

    // A symmetric five-point Gauss approximation to 42 decimal places.
    // The symmetric coordinates and reconstructed center weight keep the
    // mass and first moment exact for the actual finite rational rule.
    const I scale = pow(I{10}, 42);
    const I first("46910077030668003601186560850303517437174");
    const I second("230765344947158454481842789649895597516357");
    const I w1("118463442528094543757132020359958681321630");
    const I w2("239314335249683234020645757417819096456148");
    const auto x1 = rat(fixture, first, scale), x2 = rat(fixture, second, scale);
    const auto x4 = rat(fixture, scale - second, scale), x5 = rat(fixture, scale - first, scale);
    const auto a = rat(fixture, w1, scale), b = rat(fixture, w2, scale);
    const auto c = rat(fixture, scale - 2 * (w1 + w2), scale);
    const std::array high{PgLineQuadratureNode{x1, a}, PgLineQuadratureNode{x2, b},
                          PgLineQuadratureNode{half, c}, PgLineQuadratureNode{x4, b},
                          PgLineQuadratureNode{x5, a}};
    auto high_gram = take(PgDiscreteGram::prepare_product(
        std::move(high_basis), {high, high, depth}, grid, limits, scratch, high_plan));
    REQUIRE(high_gram.rule().size() == 25 * volumes);
    CHECK(high_gram.info().points == high_gram.rule().size());

    const auto volume_bounds = [&](const PgDiscreteGram& gram) {
        const auto entries = take(gram.enclosure()).upper_triangle;
        std::array<Q, 2> bound{};
        for (std::size_t j = 0; j < size; ++j) {
            for (std::size_t i = 0; i <= j; ++i) {
                const Q weight = translation[i] * translation[j] * (i == j ? 1 : 2);
                if (weight == 0)
                    continue;
                const auto& entry = entries[slot(i, j)];
                const Q x = weight * q(entry.lower), y = weight * q(entry.upper);
                bound[0] += std::min(x, y);
                bound[1] += std::max(x, y);
            }
        }
        return bound;
    };
    // Analytic normalized volume. Integrate the affine depth term in the
    // offset Jacobian exactly through the full thickness, leaving
    // integral_0^1 (1/2 + s/4)(1 - 3s/8) sqrt(1 + s^2/64) ds.
    const auto polynomial_integral = [](int order) -> Q {
        return Q{1, 2 * (order + 1)} + Q{1, 16 * (order + 2)} - Q{3, 32 * (order + 3)};
    };
    Q reference, term{1};
    for (int k = 0; k < 32; ++k) {
        reference += term * polynomial_integral(2 * k);
        term *= (Q("1/2") - k) / (64 * (k + 1));
    }
    const Q tail = abs(term); // alternating, decreasing terms for s in [0,1].
    REQUIRE(tail < Q("1/1000000000000000000000000000000000000000000000000000000"));
    const Q tolerance("1/1000000000000");
    const auto fine = volume_bounds(high_gram);
    CHECK(abs(fine[0] - reference) + tail < tolerance);
    CHECK(abs(fine[1] - reference) + tail < tolerance);
    const auto coarse = volume_bounds(low_gram);
    CHECK((coarse[1] < reference - tail - tolerance || coarse[0] > reference + tail + tolerance));
}
