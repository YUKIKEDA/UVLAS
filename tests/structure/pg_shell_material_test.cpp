#include <catch2/catch_test_macros.hpp>

#include <uvlas/structure/pg_shell_material.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;

struct Fixture {
    static constexpr std::uint64_t bytes = 64ULL << 20U;
    core::ResourcePlan memory{bytes};
    PgReferenceLimits limits{{256, 256, 256, 8, 8, 8192, 32768, ExactEstimateAdmission{bytes},
                              ExactTimingMode::Disabled},
                             100'000,
                             100'000,
                             {}};

    ExactRational q(std::string_view numerator, std::string_view denominator = "1") {
        auto value =
            ExactRational::from_integer_pair(numerator, denominator, limits.algebra, memory);
        REQUIRE(value);
        return std::move(*value);
    }

    std::array<ExactRational, 36> identity() {
        std::array<ExactRational, 36> result;
        for (std::size_t i = 0; i < 6; ++i)
            result[7 * i] = q("1");
        return result;
    }
};

void rejected(const core::Result<PgShellMaterial>& result, core::FailureCode code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
}

TEST_CASE("Isotropic and reciprocal orthotropic inputs yield the full 3D Mandel law",
          "[MAT-D04][MAT-D08][V-MAT][verification]") {
    Fixture f;
    const std::array<ExactRational, 6> zero;
    const PgInitialGreenStrain initial{zero, PgMaterialComponents::Mandel};
    const auto rho = f.q("4"), young = f.q("100"), poisson = f.q("1", "4");
    const auto iso = PgShellMaterial::prepare(PgIsotropicElasticity{young, poisson}, rho, initial,
                                              f.limits, f.memory);
    REQUIRE(iso);
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            const auto& value = iso->law().stiffness()[i * 6 + j];
            const auto expected = i == j ? (i < 3 ? "120" : "80") : (i < 3 && j < 3 ? "40" : "0");
            CHECK(value.rational == f.q(expected));
            CHECK(value.radical.sign() == 0);
        }
    }
    CHECK(iso->law().density_kg_m3() == rho);
    CHECK(iso->law().symmetry_correction_squared().sign() == 0);

    const std::array axes{f.q("120"), f.q("90"), f.q("60")};
    const std::array ratios{f.q("1", "4"), f.q("1", "5"), f.q("1", "3")};
    const std::array shear{f.q("7"), f.q("8"), f.q("9")};
    const auto ortho = PgShellMaterial::prepare(PgOrthotropicElasticity{axes, ratios, shear}, rho,
                                                initial, f.limits, f.memory);
    REQUIRE(ortho);
    // Exact row elimination of the reciprocal normal compliance:
    // [1/120 -1/480 -1/600; -1/480 1/90 -1/270; -1/600 -1/270 1/60].
    const std::array numerators{"4800000", "1144800", "734400", "1144800", "3810240",
                                "961200",  "734400",  "961200", "2470500"};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            const auto& value = ortho->law().stiffness()[i * 6 + j];
            CHECK(value.radical.sign() == 0);
            if (i < 3 && j < 3)
                CHECK(value.rational == f.q(numerators[3 * i + j], "36391"));
            else
                CHECK(value.rational == f.q(i == j ? std::array{"14", "16", "18"}[i - 3] : "0"));
        }
    }
}

TEST_CASE("Effective material preserves normal-shear coupling and initial strain convention",
          "[MAT-D04][MAT-D08][MAT-D09][V-MAT][verification]") {
    Fixture f;
    auto values = f.identity();
    values[3] = f.q("1", "8");
    values[18] = f.q("1", "8");
    values[22] = f.q("1", "16");
    values[27] = f.q("1", "16");
    std::array<ExactRational, 6> e0;
    e0[3] = f.q("1", "8");
    e0[4] = f.q("-1", "16");
    const auto rho = f.q("5");
    const auto effective = PgShellMaterial::prepare(
        PgEffectiveElasticity{values, PgMaterialComponents::Engineering}, rho,
        PgInitialGreenStrain{e0, PgMaterialComponents::Engineering}, f.limits, f.memory);
    REQUIRE(effective);
    CHECK(effective->law().stiffness()[3].radical == f.q("1", "8"));
    CHECK(effective->law().stiffness()[18].radical == f.q("1", "8"));
    CHECK(effective->law().stiffness()[3].rational.sign() == 0);
    CHECK(effective->law().stiffness()[21].rational == f.q("2"));
    CHECK(effective->law().stiffness()[22].rational == f.q("1", "8"));
    CHECK(effective->initial_strain()[3].radical == f.q("1", "16"));
    CHECK(effective->initial_strain()[4].radical == f.q("-1", "32"));
    CHECK(effective->law().symmetry_correction_squared().sign() == 0);

    values[1] = f.q("1", "10000000000000");
    const auto rounded = PgShellMaterial::prepare(
        PgEffectiveElasticity{values, PgMaterialComponents::Engineering}, rho,
        PgInitialGreenStrain{e0, PgMaterialComponents::Mandel}, f.limits, f.memory);
    REQUIRE(rounded);
    CHECK(rounded->law().stiffness()[1].rational == f.q("1", "20000000000000"));
    CHECK(rounded->law().stiffness()[6].rational == f.q("1", "20000000000000"));
    CHECK(rounded->initial_strain()[3].rational == f.q("1", "8"));
    CHECK(rounded->law().symmetry_correction_squared().sign() > 0);
    CHECK(effective->law().stiffness()[1].rational.sign() == 0);
}

TEST_CASE("Materials reject invalid metrics, unresolved stiffness and resource exhaustion",
          "[MAT-D08][MAT-D09][V-MAT][negative][verification]") {
    Fixture f;
    const auto rho = f.q("2");
    auto values = f.identity();
    std::array<ExactRational, 6> e0;
    const auto prepare = [&](const PgReferenceLimits& limits) {
        return PgShellMaterial::prepare(
            PgEffectiveElasticity{values, PgMaterialComponents::Mandel}, rho,
            PgInitialGreenStrain{e0, PgMaterialComponents::Engineering}, limits, f.memory);
    };
    values[1] = f.q("999999999999", "1000000000000");
    values[6] = f.q("999999999999", "1000000000000");
    // This 2-by-2 block has positive spectrum but lambda_min equals the strict cut.
    rejected(prepare(f.limits), core::FailureCode::PrecisionUnresolved);
    values[1] = f.q("999999999998", "1000000000000");
    values[6] = f.q("999999999998", "1000000000000");
    CHECK(prepare(f.limits));
    values[1] = f.q("1", "8");
    rejected(prepare(f.limits), core::FailureCode::InvalidInput);
    values = f.identity();
    values[14] = f.q("-1");
    rejected(prepare(f.limits), core::FailureCode::InvalidInput);
    values = f.identity();
    // Every leading 1- and 2-dimensional principal minor is positive.
    for (std::size_t i = 3; i < e0.size(); ++i)
        e0[i] = f.q("-3", "4");
    rejected(prepare(f.limits), core::FailureCode::InvalidInput);
    for (std::size_t i = 3; i < e0.size(); ++i)
        e0[i] = f.q("3", "4");
    CHECK(prepare(f.limits));

    auto limits = f.limits;
    limits.maximum_scalar_operations = 1;
    rejected(prepare(limits), core::FailureCode::ResourceExceeded);
    limits = f.limits;
    limits.algebra.external_memory = ExactEstimateAdmission{0};
    rejected(prepare(limits), core::FailureCode::ResourceExceeded);
    std::stop_source stopped;
    REQUIRE(stopped.request_stop());
    limits = f.limits;
    limits.cancellation = stopped.get_token();
    rejected(prepare(limits), core::FailureCode::Cancelled);
}

TEST_CASE("Material point energy and stress are conjugate in all Mandel components",
          "[MAT-D04][MAT-D08][MAT-D09][V-MAT][verification]") {
    Fixture f;
    auto values = f.identity();
    values[3] = f.q("1", "8");
    values[18] = f.q("1", "8");
    std::array<ExactRational, 6> initial;
    initial[3] = f.q("1", "8");
    const auto rho = f.q("5");
    const auto material = PgShellMaterial::prepare(
        PgEffectiveElasticity{values, PgMaterialComponents::Engineering}, rho,
        PgInitialGreenStrain{initial, PgMaterialComponents::Engineering}, f.limits, f.memory);
    REQUIRE(material);
    std::array<PgReferenceScalar, 6> green;
    green[0].rational = f.q("1", "8");
    green[3].radical = f.q("1", "8");
    const auto result = evaluate_pg_shell_material(*material, green, f.limits, f.memory);
    REQUIRE(result);
    // Independently evaluate engineering work on the two nonzero elastic
    // components (1/8,1/8). Both engineering stresses are 9/64 Pa.
    CHECK(result->energy_joules_per_cubic_metre.rational == f.q("9", "512"));
    CHECK(result->energy_joules_per_cubic_metre.radical.sign() == 0);
    for (std::size_t i = 0; i < green.size(); ++i) {
        const auto& stress = result->second_pk_mandel_pascal[i];
        CHECK(stress.rational == (i == 0 ? f.q("9", "64") : f.q("0")));
        CHECK(stress.radical == (i == 3 ? f.q("9", "64") : f.q("0")));
    }
    const auto unloaded =
        evaluate_pg_shell_material(*material, material->initial_strain(), f.limits, f.memory);
    REQUIRE(unloaded);
    CHECK(unloaded->energy_joules_per_cubic_metre.rational.sign() == 0);
    CHECK(unloaded->energy_joules_per_cubic_metre.radical.sign() == 0);
    for (const auto& stress : unloaded->second_pk_mandel_pascal) {
        CHECK(stress.rational.sign() == 0);
        CHECK(stress.radical.sign() == 0);
    }
}
TEST_CASE("Shell material bounds retain independent radicals and initial field errors",
          "[MAT-D04][MAT-D08][V-MAT][verification][negative]") {
    Fixture f;
    auto stiffness = f.identity();
    stiffness[3] = f.q("1", "8");
    stiffness[18] = f.q("1", "8");
    const auto law = PgShellElasticLaw::prepare(
        PgEffectiveElasticity{stiffness, PgMaterialComponents::Engineering}, f.q("1"), f.limits,
        f.memory);
    REQUIRE(law);
    const auto g = f.q("3");
    std::array<PgReferenceScalar, 6> initial;
    initial[0].radical = f.q("1", "8");
    std::array<ExactRational, 6> uncertainty;
    const std::array<PgMaterialEnclosure, 6> green;
    const PgInitialGreenEnclosure input{g, initial, uncertainty};
    const auto evaluated = enclose_pg_shell_material(*law, input, green, f.limits, f.memory);
    REQUIRE(evaluated);
    const auto compare = [&](const ExactRational& x, const ExactRational& y) {
        const auto relation = x.compare(y, f.limits.algebra, f.memory);
        REQUIRE(relation);
        return *relation;
    };
    // S_yz in Mandel form is -sqrt(6)/64. Compare through the exact square
    // so the test does not depend on a binary approximation of either radical.
    const auto& shear = evaluated->second_pk_mandel_pascal[3];
    REQUIRE(shear.upper.sign() < 0);
    const auto lower_square = shear.lower.multiply(shear.lower, f.limits.algebra, f.memory);
    const auto upper_square = shear.upper.multiply(shear.upper, f.limits.algebra, f.memory);
    REQUIRE(lower_square);
    REQUIRE(upper_square);
    const auto expected_square = f.q("3", "2048");
    CHECK(compare(*lower_square, expected_square) > 0);
    CHECK(compare(*upper_square, expected_square) < 0);
    const auto width = shear.upper.subtract(shear.lower, f.limits.algebra, f.memory);
    REQUIRE(width);
    CHECK(compare(*width, f.q("1", "100000000000000000000000000000")) < 0);
    const auto& energy = evaluated->energy_joules_per_cubic_metre;
    const auto expected_energy = f.q("3", "128");
    CHECK(compare(energy.lower, expected_energy) < 0);
    CHECK(compare(energy.upper, expected_energy) > 0);

    initial[0] = {f.q("-49", "100"), {}};
    uncertainty[0] = f.q("1", "50");
    const auto unresolved = enclose_pg_shell_material(*law, input, green, f.limits, f.memory);
    REQUIRE_FALSE(unresolved);
    CHECK(unresolved.error().code == core::FailureCode::PrecisionUnresolved);
    initial[0].rational = f.q("-3", "4");
    const auto invalid = enclose_pg_shell_material(*law, input, green, f.limits, f.memory);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
}
} // namespace
