// Checks that the headers reaching the compiler are the versions named in
// vcpkg.json and that the build keeps the floating-point
// model the adopted numerics rely on (N4-01, N4-06). This is a toolchain check,
// not a verification of any solver path.

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/Util/Version.h>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/version.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_version_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace {
// Override version schemes: https://learn.microsoft.com/en-us/vcpkg/reference/vcpkg-json#overrides
// This reads an exact pin, not a version range or a full manifest schema.
[[nodiscard]] std::optional<std::string> pinned_version(const nlohmann::json& item) {
    if (!item.is_object()) {
        return std::nullopt;
    }
    std::optional<std::string> result;
    for (const auto* field : {"version", "version-semver", "version-date", "version-string"}) {
        const auto value = item.find(field);
        if (value == item.end()) {
            continue;
        }
        if (result.has_value() || !value->is_string()) {
            return std::nullopt;
        }
        result = value->get<std::string>();
        if (result->empty()) {
            return std::nullopt;
        }
    }
    return result;
}
} // namespace

TEST_CASE("Override pin fields are unique nonempty strings", "[toolchain][negative][N4-01]") {
    constexpr std::array<std::pair<std::string_view, std::string_view>, 4> fields{
        {{"version", "1.2.3"},
         {"version-semver", "2.3.4"},
         {"version-date", "2026-11-15"},
         {"version-string", "release-tag"}}};
    // All sixteen subsets cover missing, each valid field, and every ambiguity.
    for (unsigned mask = 0; mask < 16; ++mask) {
        nlohmann::json item = {{"name", "fixture"}};
        std::string expected;
        for (std::size_t i = 0; i < fields.size(); ++i) {
            if ((mask & (1U << i)) != 0) {
                item[fields[i].first] = fields[i].second;
                expected = fields[i].second;
            }
        }
        const auto actual = pinned_version(item);
        REQUIRE(actual.has_value() == std::has_single_bit(mask));
        if (actual) {
            CHECK(*actual == expected);
        }
    }
    for (const auto& [field, value] : fields) {
        (void)value;
        for (const auto& invalid : nlohmann::json::array(
                 {nullptr, false, 1, 1.5, "", nlohmann::json::array(), nlohmann::json::object()})) {
            CHECK_FALSE(pinned_version(nlohmann::json{{field, invalid}}));
        }
    }
    CHECK_FALSE(pinned_version(nullptr));
    CHECK_FALSE(pinned_version(nlohmann::json::array()));
}

TEST_CASE("pinned dependency versions reach the compiler", "[toolchain][N4-01]") {
    std::ifstream input(UVLAS_MANIFEST_PATH);
    REQUIRE(input.is_open());
    const auto manifest = nlohmann::json::parse(input);
    std::map<std::string, std::string> versions;
    for (const auto& item : manifest.at("overrides")) {
        const auto pin = pinned_version(item);
        if (!pin.has_value()) {
            FAIL("Override must contain one nonempty version pin");
            return;
        }
        versions.emplace(item.at("name").get<std::string>(), *pin);
        if (item.at("name") == "nlohmann-json") {
            CHECK(item.at("port-version").get<int>() == 3);
        }
    }
    const auto version = [](int major, int minor, int patch) {
        return std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(patch);
    };
    CHECK(version(EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, EIGEN_PATCH_VERSION) ==
          versions.at("eigen3"));
    CHECK(version(SPECTRA_MAJOR_VERSION, SPECTRA_MINOR_VERSION, SPECTRA_PATCH_VERSION) ==
          versions.at("spectra"));
    const auto boost_version =
        version(BOOST_VERSION / 100000, BOOST_VERSION / 100 % 1000, BOOST_VERSION % 100);
    CHECK(boost_version == versions.at("boost-multiprecision"));
    CHECK(boost_version == versions.at("boost-hash2"));
    CHECK(version(NLOHMANN_JSON_VERSION_MAJOR, NLOHMANN_JSON_VERSION_MINOR,
                  NLOHMANN_JSON_VERSION_PATCH) == versions.at("nlohmann-json"));
    CHECK(version(CATCH_VERSION_MAJOR, CATCH_VERSION_MINOR, CATCH_VERSION_PATCH) ==
          versions.at("catch2"));
    // Eigen 5 keeps EIGEN_WORLD_VERSION at 3; the release is MAJOR.MINOR.PATCH.
    CHECK(EIGEN_MAJOR_VERSION == 5);
    CHECK(EIGEN_MINOR_VERSION == 0);
    CHECK(EIGEN_PATCH_VERSION == 1);
    CHECK(SPECTRA_MAJOR_VERSION == 1);
    CHECK(SPECTRA_MINOR_VERSION == 2);
    CHECK(SPECTRA_PATCH_VERSION == 0);
    CHECK(BOOST_VERSION == 109200);
    CHECK(NLOHMANN_JSON_VERSION_MAJOR == 3);
    CHECK(NLOHMANN_JSON_VERSION_MINOR == 12);
    CHECK(NLOHMANN_JSON_VERSION_PATCH == 0);
#ifdef UVLAS_NLOHMANN_JSON_ALLOCATION_REPAIR
    CHECK(UVLAS_NLOHMANN_JSON_ALLOCATION_REPAIR == 1);
#else
    FAIL("The JSON package must provide the UVLAS allocation-exception repair");
#endif
}

TEST_CASE("numerical build policy is in effect", "[toolchain][N4-01][N4-02]") {
#ifndef EIGEN_MPL2_ONLY
    FAIL("EIGEN_MPL2_ONLY is not defined");
#endif
#ifndef EIGEN_DONT_PARALLELIZE
    FAIL("EIGEN_DONT_PARALLELIZE is not defined");
#endif
    // N4-02: indices are signed 64-bit on every supported platform.
    STATIC_REQUIRE(sizeof(Eigen::Index) == 8);
    STATIC_REQUIRE(std::numeric_limits<Eigen::Index>::is_signed);
    STATIC_REQUIRE(std::numeric_limits<double>::is_iec559);
}

TEST_CASE("fma gives the exact product error", "[toolchain][N4-06]") {
    // Error-free two-product: p + e == a * b exactly. A build that contracts or
    // reorders floating-point expressions breaks the compensated residuals.
    const volatile double a = 1.0 + std::ldexp(1.0, -30);
    const volatile double b = 1.0 - std::ldexp(1.0, -31);
    const double p = a * b;
    const double e = std::fma(a, b, -p);

    using Big = boost::multiprecision::cpp_bin_float_100;
    const Big exact = Big(static_cast<double>(a)) * Big(static_cast<double>(b));
    CHECK(Big(p) + Big(e) == exact);
    CHECK(e != 0.0);
}
