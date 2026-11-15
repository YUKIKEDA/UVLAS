#include <boost/hash2/digest.hpp>
#include <boost/hash2/sha2.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <expected>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <type_traits>
#include <uvlas/numerics/dense_solve.hpp>
#include <vector>

namespace {
using namespace uvlas;
using namespace numerics;
using namespace numerics::linear;
using Big = boost::multiprecision::cpp_bin_float_100;
using Json = nlohmann::json;
constexpr core::Index order = 24;
constexpr auto operation = core::Operation::SolveLinear;
struct Fixture {
    std::vector<double> matrix, rhs;
    std::vector<Big> reference;
};
[[nodiscard]] std::string read(const std::filesystem::path& path) {
    std::ifstream stream{path, std::ios::binary};
    return {std::istreambuf_iterator<char>{stream}, std::istreambuf_iterator<char>{}};
}
[[nodiscard]] std::string sha256(const std::string& bytes) {
    boost::hash2::sha2_256 hash;
    hash.update(bytes.data(), bytes.size());
    return boost::hash2::to_string(hash.result());
}
[[nodiscard]] std::expected<Fixture, core::FailureCode> decode(const Json& data) {
    if (!data.is_object() || data.value("schema", "") != "uvlas.foundation-reference/1" ||
        data.value("id", "") != "j1-pg-free-history-24-v1") {
        return std::unexpected(core::FailureCode::InvalidInput);
    }
    const auto& layout = data.at("layout");
    if (layout != Json{{"rows", order},
                       {"columns", order},
                       {"order", "row-major"},
                       {"state_order", "logical-index-0-23"},
                       {"units", "abstract-dimensionless-replay-only"}} ||
        data.at("operator") != "rounded-identity-minus-stored-P" ||
        data.at("relative_solution_tolerance") != 1e-12 ||
        data.at("normalized_residual_tolerance") != 1e-12) {
        return std::unexpected(core::FailureCode::InvalidInput);
    }
    Fixture result;
    const auto load = [](const Json& values, std::size_t size,
                         std::vector<double>& target) -> std::expected<void, core::FailureCode> {
        if (!values.is_array() || values.size() != size) {
            return std::unexpected(core::FailureCode::InvalidInput);
        }
        for (const auto& value : values) {
            if (!value.is_number()) {
                return std::unexpected(core::FailureCode::InvalidInput);
            }
            const double scalar = value.get<double>();
            if (!std::isfinite(scalar)) {
                return std::unexpected(core::FailureCode::Nonfinite);
            }
            target.push_back(scalar);
        }
        return {};
    };
    auto loaded = load(data.at("matrix"), std::size_t{24} * 24, result.matrix);
    if (!loaded) {
        return std::unexpected(loaded.error());
    }
    loaded = load(data.at("rhs"), 24U, result.rhs);
    if (!loaded) {
        return std::unexpected(loaded.error());
    }
    const auto& reference = data.at("reference_solution");
    if (!reference.is_array() || reference.size() != 24) {
        return std::unexpected(core::FailureCode::InvalidInput);
    }
    for (const auto& value : reference) {
        if (!value.is_string()) {
            return std::unexpected(core::FailureCode::InvalidInput);
        }
        const Big scalar{value.get<std::string>()};
        if (!std::isfinite(scalar.convert_to<double>())) {
            return std::unexpected(core::FailureCode::Nonfinite);
        }
        result.reference.push_back(scalar);
    }
    return result;
}
[[nodiscard]] Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            0.25};
}
class NoPhysicalChecks final : public ConsumerEvaluator<double> {
public:
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>&, const Context&,
           const ConsumerCheck& check, core::CertificateContext expected) const override {
        FAIL_CHECK("The algebraic replay must not claim a physical observation");
        return core::ErrorCertificate::unevaluated(expected, check.budget);
    }
};
[[nodiscard]] Json document() {
    return Json::parse(
        read(std::filesystem::path{UVLAS_FOUNDATION_REFERENCE_ROOT} / "fixture.json"));
}
[[nodiscard]] std::vector<Big> actual(const PrecisionBuffer<double>& values) {
    return std::visit(
        [](const auto& storage) {
            std::vector<Big> result;
            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(storage)>,
                                         Binary64<double>>) {
                for (const double value : storage.values.values()) {
                    result.emplace_back(value);
                }
            } else {
                for (std::size_t i = 0; i < storage.hi.values().size(); ++i) {
                    result.emplace_back(Big{storage.hi.values()[i]} + Big{storage.lo.values()[i]});
                }
            }
            return result;
        },
        values.storage());
}
} // namespace

TEST_CASE(
    "Frozen PG matrix replay exercises the product solver against an independent Decimal oracle",
    "[N4-10][N4-05][N4-06]") {
    const auto data = document();
    auto fixture = decode(data);
    REQUIRE(fixture);
    core::ResourcePlan plan{std::uint64_t{16} * 1024 * 1024};
    std::vector<double> column_major(std::size_t{24} * 24), scales(24, 1);
    std::vector<UnitDimensions> units(24);
    // Only the numerical coefficients are replayed. This is not reconstruction of
    // physical PG/free-wake state units or of the differentiated step itself.
    for (std::size_t i = 0; i < 24; ++i) {
        for (std::size_t j = 0; j < 24; ++j) {
            column_major[j * 24 + i] = fixture->matrix[i * 24 + j];
        }
    }
    auto matrix = DenseReal::copy({order, order, order}, column_major, plan, operation);
    auto rhs_values = DenseReal::copy({order, 1, order}, fixture->rhs, plan, operation);
    REQUIRE(matrix);
    REQUIRE(rhs_values);
    auto buffer =
        PrecisionBuffer<double>::binary64(std::move(*rhs_values), core::NumericalGeneration{30});
    REQUIRE(buffer);
    const DenseOperator op{*matrix, context()};
    const Scales scale{scales, scales, units, units, context()};
    const Limits limits{order, 1, 0, 5, EstimateAdmission{std::uint64_t{16} * 1024 * 1024}};
    auto prepared = DensePrepared::prepare(op, scale, limits, plan);
    REQUIRE(prepared);
    const NoPhysicalChecks evaluator;
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    const Accuracy accuracy{Norm::Infinity,         1e-12, 1e-12, 1e-12, core::ObservationId{20},
                            core::ObservationId{21}};
    const Rhs<double> rhs{*buffer, context(), core::RhsGeneration{40}};
    const auto original_rhs = actual(*buffer);
    auto solved = solve_linear(op, rhs, scale, precision, accuracy, limits, *prepared, plan);
    REQUIRE(solved);
    const auto values = actual(solved->values());
    for (std::size_t i = 0; i < 24; ++i) {
        const Big magnitude = abs(fixture->reference[i]);
        const Big denominator = magnitude > 1 ? magnitude : Big{1};
        CHECK(abs(values[i] - fixture->reference[i]) / denominator <= Big{"1e-12"});
        Big residual{fixture->rhs[i]};
        for (std::size_t j = 0; j < 24; ++j) {
            residual -= Big{fixture->matrix[i * 24 + j]} * values[j];
        }
        CHECK(abs(residual) <= Big{"1e-12"} * (1 + abs(Big{fixture->rhs[i]})));
    }
    CHECK(solved->certificate().columns()[0].normalized_original_residual_upper <= 1);
    core::ResourcePlan no_output{0};
    auto failed = solve_linear(op, rhs, scale, precision, accuracy, limits, *prepared, no_output);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
    CHECK(actual(solved->values()) == values);
    CHECK(actual(*buffer) == original_rhs);
    CHECK(buffer->generation() == core::NumericalGeneration{30});
    CHECK(std::ranges::equal(matrix->values(), column_major));
}

TEST_CASE("Frozen fixture adapter rejects invalid layouts and nonfinite coefficients",
          "[N4-10][N4-08]") {
    const auto original = document();
    for (const int fault : {0, 1, 2, 3, 4, 5}) {
        auto bad = original;
        switch (fault) {
        case 0:
            bad["layout"]["rows"] = 23;
            break;
        case 1:
            bad["matrix"].erase(bad["matrix"].begin());
            break;
        case 2:
            bad["rhs"][0] = true;
            break;
        case 3:
            bad["matrix"][0] = std::numeric_limits<double>::infinity();
            break;
        case 4:
            bad["relative_solution_tolerance"] = 1e-6;
            break;
        default:
            bad["layout"]["order"] = "column-major";
            break;
        }
        auto rejected = decode(bad);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error() ==
              (fault == 3 ? core::FailureCode::Nonfinite : core::FailureCode::InvalidInput));
    }
}

TEST_CASE("Frozen PG source artifacts preserve the indexed original hashes", "[N4-10]") {
    const auto data = document();
    const std::filesystem::path root{UVLAS_FOUNDATION_REFERENCE_ROOT};
    for (const auto& artifact : data.at("source_artifacts")) {
        const auto name = std::filesystem::path{artifact.at("hint").get<std::string>()}.filename();
        CHECK(sha256(read(root / "source" / name)) == artifact.at("sha256").get<std::string>());
    }
    // This checks only the internal consistency of an imported full-spectrum
    // oracle. No production eigensolver or physical mode coverage is claimed.
    Big trace = 0, eigen_sum = 0, imaginary_sum = 0;
    for (std::size_t i = 0; i < 24; ++i) {
        trace += Big{data.at("step_values").at(i * 24 + i).get<double>()};
        eigen_sum += Big{data.at("eigenvalues").at(i).at(0).get<double>()};
        imaginary_sum += Big{data.at("eigenvalues").at(i).at(1).get<double>()};
    }
    CHECK(abs(trace - eigen_sum) <= Big{"1e-12"});
    CHECK(abs(imaginary_sum) <= Big{"1e-12"});
}
