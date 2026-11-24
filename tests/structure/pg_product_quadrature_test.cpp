#include "pg_conditioning_test_support.hpp"
#include "pg_rank_storage.hpp"
#include "pg_reference_arithmetic.hpp"

#include <algorithm>
#include <uvlas/structure/pg_discrete_gram.hpp>

using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;

TEST_CASE("PG common denominators admit guarded work before evaluating any term",
          "[SH-D07][N4-08][verification][negative]") {
    using Reuse = detail::PgRankMonomialPowers::Reuse;
    constexpr auto limb_bits = std::numeric_limits<detail::RankWord>::digits;
    const I denominator = I{1} << (4 * limb_bits);
    Fixture fixture;
    const auto u = rat(fixture, I{1}, denominator);
    const auto v = fixture.r("0");
    const std::array rows{PgRankRow{PgCoefficientField::R, {0, 0}, {1, 0}},
                          PgRankRow{PgCoefficientField::R, {0, 0}, {2, 0}}};
    const std::array<std::uint32_t, 2> indices{0, 1};
    using Coefficient = detail::RankFraction<detail::RankCell>;
    const std::array coefficients{Coefficient{1, 1}, Coefficient{1, 1}};
    core::ResourcePlan plan{profile};
    auto limits = query_limits();
    limits.algebra.maximum_value_bits = 11 * limb_bits;

    // The bound is 8 limbs plus 5 bits. The next limb and both guard limbs
    // must fit, even though the final reduced fraction has a smaller width.
    limits.algebra.maximum_work_bits = 10 * limb_bits;
    detail::ReferenceArithmetic rejected{limits, plan};
    detail::PgRankMonomialPowers powers{u, v, Reuse::CommonDenominator};
    CHECK_FALSE(powers.common_sum(rejected, rows, indices, coefficients));
    CHECK(rejected.info().scalar_operations == 0);
    CHECK(plan.allocated_bytes() == 0);

    limits.algebra.maximum_work_bits = 11 * limb_bits;
    detail::ReferenceArithmetic admitted{limits, plan};
    const auto result = powers.common_sum(admitted, rows, indices, coefficients);
    REQUIRE(result);
    CHECK(q(*result) == Q{denominator + 1} / (denominator * denominator));
}

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
    const auto five_sixths = fixture.r("5", "6");
    const auto depth_left_weight = fixture.r("4", "7"), depth_right_weight = fixture.r("3", "7");
    const std::array depth_pair{PgLineQuadratureNode{quarter, depth_left_weight},
                                PgLineQuadratureNode{five_sixths, depth_right_weight}};
    const std::array low{PgLineQuadratureNode{quarter, half},
                         PgLineQuadratureNode{three_quarters, half}};
    const PgDiscreteGramGrid grid{260, 192};
    const PgDiscreteGramLimits limits{query_limits(), steps};
    const std::array wrong_moment{PgLineQuadratureNode{quarter, one}};
    const std::array duplicate{PgLineQuadratureNode{half, half}, PgLineQuadratureNode{half, half}};
    const std::array zero_weight{PgLineQuadratureNode{half, zero}};
    const auto prepare = [&](PgReferenceBasis& basis, PgVolumeProductRule input,
                             core::ResourcePlan& plan) {
        return PgDiscreteGram::prepare_product(std::move(basis), input, grid, limits, scratch,
                                               plan);
    };
    core::ResourcePlan rejected{profile};
    for (const auto& input : std::array{PgVolumeProductRule{wrong_moment, low, depth},
                                        PgVolumeProductRule{low, duplicate, depth},
                                        PgVolumeProductRule{low, low, zero_weight},
                                        PgVolumeProductRule{low, low, {}}}) {
        check_failure(prepare(low_basis, input, rejected), core::FailureCode::InvalidInput);
        CHECK(take(low_basis.size()) == size);
    }
    core::ResourcePlan no_storage{0};
    check_failure(prepare(low_basis, {low, low, depth_pair}, no_storage),
                  core::FailureCode::ResourceExceeded);
    // Payload admission precedes exact verification of the radial moment.
    check_failure(prepare(low_basis, {wrong_moment, low, depth}, no_storage),
                  core::FailureCode::ResourceExceeded);
    CHECK(no_storage.allocated_bytes() == 0);
    auto low_gram = take(prepare(low_basis, {low, low, depth_pair}, low_plan));
    REQUIRE(low_gram.rule().size() == 8 * volumes);
    CHECK(low_gram.info().points == low_gram.rule().size());
    const std::array<Q, 4> u{Q("1/4"), Q("1/4"), Q("3/4"), Q("3/4")};
    const std::array<Q, 4> v{Q("3/16"), Q("9/16"), Q("1/16"), Q("3/16")};
    const std::array<Q, 2> d{Q("1/4"), Q("5/6")};
    const std::array<Q, 2> depth_weights{Q("4/7"), Q("3/7")};
    for (std::size_t i = 0; i < volumes; ++i) {
        Q weight_sum, depth_first, depth_second;
        for (std::size_t j = 0; j < 4; ++j) {
            for (std::size_t k = 0; k < d.size(); ++k) {
                const auto& point = low_gram.rule()[8 * i + 2 * j + k];
                CHECK(point.volume == static_cast<core::Index>(i));
                CHECK(q(point.u) == u[j]);
                CHECK(q(point.v) == v[j]);
                CHECK(q(point.depth) == d[k]);
                CHECK(q(point.weight) == (1 - u[j]) * depth_weights[k] / 4);
                weight_sum += q(point.weight);
                depth_first += q(point.depth) * q(point.weight);
                depth_second += q(point.depth) * q(point.depth) * q(point.weight);
            }
        }
        CHECK(weight_sum == Q("1/2"));
        CHECK(depth_first == Q("1/4"));
        CHECK(depth_second == Q("1/6"));
    }

    // Recompute every value from the same nodes through the explicit point path.
    std::vector<PgDiscreteVolumePoint> explicit_points;
    explicit_points.reserve(low_gram.rule().size());
    for (const auto& point : low_gram.rule()) {
        explicit_points.push_back({{point.volume, point.u, point.v, point.depth}, point.weight});
    }
    auto explicit_basis = make_curve_basis();
    core::ResourcePlan explicit_plan{profile};
    const auto explicit_gram = take(PgDiscreteGram::prepare(
        std::move(explicit_basis), explicit_points, grid, limits, scratch, explicit_plan));
    const auto product_entries = take(low_gram.enclosure()).upper_triangle;
    const auto explicit_entries = take(explicit_gram.enclosure()).upper_triangle;
    REQUIRE(product_entries.size() == explicit_entries.size());
    CHECK(low_gram.info().arithmetic.scalar_operations <
          explicit_gram.info().arithmetic.scalar_operations);
    for (std::size_t i = 0; i < product_entries.size(); ++i) {
        CHECK(product_entries[i].lower == explicit_entries[i].lower);
        CHECK(product_entries[i].upper == explicit_entries[i].upper);
    }

    // Admit mixed powers while excluding the integer powers and scratch used
    // by the common-denominator sum.
    auto axis_limits = limits;
    const auto held = detail::basis_live_rationals + 16 + 4 * low_gram.rule().size() +
                      2 * product_entries.size() + 3 * size +
                      detail::PgRankMonomialPowers::retained_rationals;
    const auto bytes_per_integer = 4 * (axis_limits.arithmetic.algebra.maximum_work_bits / 8);
    auto mixed_limits = limits;
    const auto mixed_held = held + detail::PgRankMonomialPowers::retained_mixed;
    mixed_limits.arithmetic.algebra.external_memory =
        ExactEstimateAdmission{(2 * (mixed_held + 32) + 1) * bytes_per_integer};
    auto mixed_basis = make_curve_basis();
    core::ResourcePlan mixed_plan{profile};
    const auto mixed_gram = take(PgDiscreteGram::prepare_product(
        std::move(mixed_basis), {low, low, depth_pair}, grid, mixed_limits, scratch, mixed_plan));
    const auto mixed_entries = take(mixed_gram.enclosure()).upper_triangle;
    REQUIRE(mixed_entries.size() == product_entries.size());
    CHECK(low_gram.info().arithmetic.allocated_owned_payload_bytes <
          mixed_gram.info().arithmetic.allocated_owned_payload_bytes);
    for (std::size_t i = 0; i < product_entries.size(); ++i) {
        CHECK(product_entries[i].lower == mixed_entries[i].lower);
        CHECK(product_entries[i].upper == mixed_entries[i].upper);
    }

    // Admit axis powers while excluding the additional array of mixed powers.
    axis_limits.arithmetic.algebra.external_memory =
        ExactEstimateAdmission{(2 * (held + 64) + 1) * bytes_per_integer};
    auto axis_basis = make_curve_basis();
    core::ResourcePlan axis_plan{profile};
    const auto axis_gram = take(PgDiscreteGram::prepare_product(
        std::move(axis_basis), {low, low, depth_pair}, grid, axis_limits, scratch, axis_plan));
    const auto axis_entries = take(axis_gram.enclosure()).upper_triangle;
    REQUIRE(axis_entries.size() == product_entries.size());
    CHECK(low_gram.info().arithmetic.scalar_operations <
          axis_gram.info().arithmetic.scalar_operations);
    for (std::size_t i = 0; i < product_entries.size(); ++i) {
        CHECK(product_entries[i].lower == axis_entries[i].lower);
        CHECK(product_entries[i].upper == axis_entries[i].upper);
    }

    // A smaller exact value limit leaves room for this geometry and finite rule,
    // but not for every cached power product allowed by the rank dictionary.
    auto narrow_basis = make_curve_basis();
    auto narrow_limits = limits;
    narrow_limits.arithmetic.algebra.maximum_value_bits = 1536;
    core::ResourcePlan narrow_plan{profile};
    const auto narrow_gram =
        take(PgDiscreteGram::prepare_product(std::move(narrow_basis), {low, low, depth_pair}, grid,
                                             narrow_limits, scratch, narrow_plan));
    const auto narrow_entries = take(narrow_gram.enclosure()).upper_triangle;
    REQUIRE(narrow_entries.size() == product_entries.size());
    CHECK(low_gram.info().arithmetic.scalar_operations <
          narrow_gram.info().arithmetic.scalar_operations);
    CHECK(axis_gram.info().arithmetic.scalar_operations <
          narrow_gram.info().arithmetic.scalar_operations);
    for (std::size_t i = 0; i < product_entries.size(); ++i) {
        CHECK(product_entries[i].lower == narrow_entries[i].lower);
        CHECK(product_entries[i].upper == narrow_entries[i].upper);
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
    auto high_gram = take(prepare(high_basis, {high, high, depth}, high_plan));
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

TEST_CASE("PG common denominator Horner preserves sparse signed polynomial values",
          "[SH-D07][N4-08][verification][negative]") {
    using Reuse = detail::PgRankMonomialPowers::Reuse;
    using Coefficient = detail::RankFraction<detail::RankCell>;
    constexpr auto field = PgCoefficientField::R;
    const std::array rows{PgRankRow{field, {0, 0}, {1, 1}}, PgRankRow{field, {0, 0}, {1, 4}},
                          PgRankRow{field, {0, 0}, {3, 0}}, PgRankRow{field, {0, 0}, {6, 2}},
                          PgRankRow{field, {0, 0}, {6, 3}}};
    const std::array<std::uint32_t, rows.size()> indices{0, 1, 2, 3, 4};
    const std::array coefficients{Coefficient{-3, 5}, Coefficient{7, 6}, Coefficient{9, 10},
                                  Coefficient{-11, 12}, Coefficient{13, 30}};
    Fixture fixture;
    const auto u = fixture.r("-2", "3");
    const auto v = fixture.r("-3", "7");
    const Q x{-2, 3}, y{-3, 7};
    const Q x3 = x * x * x, y2 = y * y;
    // Expanded rational sum independent of the integer recurrence and its LCM.
    const Q expected = Q{-3, 5} * x * y + Q{7, 6} * x * y2 * y2 + Q{9, 10} * x3 +
                       Q{-11, 12} * x3 * x3 * y2 + Q{13, 30} * x3 * x3 * y2 * y;
    auto limits = query_limits();
    core::ResourcePlan plan{profile};
    detail::ReferenceArithmetic math{limits, plan};
    detail::PgRankMonomialPowers powers{u, v, Reuse::CommonDenominator};
    const auto value = powers.common_sum(math, rows, indices, coefficients);
    CHECK(value.transform([](const auto& exact) { return q(exact); }) == expected);

    const ExactRational zero;
    detail::PgRankMonomialPowers zero_coordinate{u, zero, Reuse::CommonDenominator};
    const auto only_zero_power = zero_coordinate.common_sum(math, rows, indices, coefficients);
    const Q expected_with_zero = Q{9, 10} * x3;
    CHECK(only_zero_power.transform([](const auto& exact) { return q(exact); }) ==
          expected_with_zero);

    const std::array cancel_rows{PgRankRow{field, {0, 0}, {0, 0}},
                                 PgRankRow{field, {0, 0}, {1, 0}}};
    const std::array<std::uint32_t, 2> cancel_indices{0, 1};
    const std::array cancel_coefficients{Coefficient{10, 3}, Coefficient{5, 1}};
    const auto cancelled =
        powers.common_sum(math, cancel_rows, cancel_indices, cancel_coefficients);
    CHECK(cancelled.transform([](const auto& exact) { return exact.sign(); }) == 0);

    const std::array<std::uint32_t, rows.size()> reversed{4, 3, 2, 1, 0};
    const auto before = math.info();
    CHECK_FALSE(powers.common_sum(math, rows, reversed, coefficients).has_value());
    CHECK(math.info().scalar_operations == before.scalar_operations);
    CHECK(math.info().allocated_owned_payload_bytes == before.allocated_owned_payload_bytes);
}
