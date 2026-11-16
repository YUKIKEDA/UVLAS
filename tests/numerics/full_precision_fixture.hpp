#pragma once

#include <boost/hash2/sha2.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace uvlas::test_support::full_precision {
using FixtureComplex = std::complex<double>;
struct Mode {
    std::string name;
    double h = 0, alpha = 0, xp = 0;
    std::array<double, 3> scales{};
};
struct Observation {
    std::string id, mode, kind;
    std::size_t column = 0;
    FixtureComplex value{};
    double scale = 0, lower = 0, upper = 0;
};
struct Fixture {
    std::int64_t m = 0;
    double k = 0;
    std::vector<double> xv, xc, fc, fu, mc0, mu0;
    std::vector<FixtureComplex> update, rhs, reference_hi, reference_lo;
    std::array<Mode, 3> modes;
    std::array<Observation, 39> observations;
    std::string base_sha256, effective_sha256, rhs_sha256;
};

inline double fixture_number(const nlohmann::json& value) {
    if (!value.is_number())
        throw std::runtime_error("Non-numeric fixture field");
    const double result = value.get<double>();
    if (!std::isfinite(result))
        throw std::runtime_error("Nonfinite fixture field");
    return result;
}
inline unsigned fixture_nibble(char value) {
    if (value >= '0' && value <= '9')
        return static_cast<unsigned>(value - '0');
    if (value >= 'a' && value <= 'f')
        return static_cast<unsigned>(value - 'a' + 10);
    throw std::runtime_error("Noncanonical fixture hex digit");
}
inline std::vector<unsigned char>
fixture_payload(const nlohmann::json& value, std::span<const std::size_t> shape, bool complex) {
    if (value.at("dtype") != (complex ? "le-c128" : "le-f64") || value.at("order") != "C")
        throw std::runtime_error("Unexpected fixture array representation");
    const auto& dimensions = value.at("shape");
    if (!dimensions.is_array() || dimensions.size() != shape.size())
        throw std::runtime_error("Unexpected fixture array rank");
    std::size_t count = 1;
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (!dimensions[i].is_number_integer() || dimensions[i].get<std::int64_t>() < 0 ||
            dimensions[i].get<std::uint64_t>() != shape[i])
            throw std::runtime_error("Unexpected fixture array shape");
        // Only fixed 4096/4097 and three-column shapes reach this test helper.
        count *= shape[i];
    }
    const std::size_t bytes = count * (complex ? 16U : 8U);
    const auto hex = value.at("payload_hex").get<std::string>();
    if (hex.size() != bytes * 2)
        throw std::runtime_error("Truncated fixture array payload");
    std::vector<unsigned char> result(bytes);
    for (std::size_t i = 0; i < bytes; ++i)
        result[i] = static_cast<unsigned char>((fixture_nibble(hex[2 * i]) << 4U) |
                                               fixture_nibble(hex[2 * i + 1]));
    boost::hash2::sha2_256 hash;
    hash.update(result.data(), result.size());
    if (boost::hash2::to_string(hash.result()) != value.at("sha256").get<std::string>())
        throw std::runtime_error("Fixture array hash mismatch");
    return result;
}
inline double fixture_word(std::span<const unsigned char> bytes, std::size_t offset) {
    std::uint64_t bits = 0;
    for (unsigned i = 0; i < 8; ++i)
        bits |= static_cast<std::uint64_t>(bytes[offset + i]) << (8U * i);
    const double result = std::bit_cast<double>(bits);
    if (!std::isfinite(result))
        throw std::runtime_error("Nonfinite fixture array word");
    return result;
}
inline std::vector<double> fixture_reals(const nlohmann::json& value, std::size_t count) {
    const auto bytes = fixture_payload(value, std::array{count}, false);
    std::vector<double> result(count);
    for (std::size_t i = 0; i < count; ++i)
        result[i] = fixture_word(bytes, i * 8);
    return result;
}
inline std::vector<FixtureComplex> fixture_complex(const nlohmann::json& value,
                                                   std::span<const std::size_t> shape) {
    const auto bytes = fixture_payload(value, shape, true);
    std::vector<FixtureComplex> result(bytes.size() / 16);
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = {fixture_word(bytes, i * 16), fixture_word(bytes, i * 16 + 8)};
    return result;
}
inline nlohmann::json checked_fixture_json(std::string_view bytes) {
    boost::hash2::sha2_256 hash;
    hash.update(bytes.data(), bytes.size());
    // Bind observation scales/values, coefficient recipes and both reference
    // words to reviewed bytes, rather than trusting hashes stored in those bytes.
    if (boost::hash2::to_string(hash.result()) !=
        "9ce2ec50d12a321f0c64dabb052a1a54358cb4c27422751c17b5a2c173f8e1ee")
        throw std::runtime_error("Frozen full precision fixture bytes changed");
    return nlohmann::json::parse(bytes);
}
inline Fixture load_fixture(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
        throw std::runtime_error("Cannot open full precision fixture");
    const std::string bytes{std::istreambuf_iterator<char>{stream}, {}};
    const auto data = checked_fixture_json(bytes);
    if (data.at("schema") != "uvlas.precision-system/1" ||
        data.at("id") != "precision-4096-system-v1" || !data.at("m").is_number_integer() ||
        data.at("m").get<std::int64_t>() != 4096 || fixture_number(data.at("k")) != .1)
        throw std::runtime_error("Unexpected full precision fixture identity");
    const auto& original = data.at("original");
    Fixture result;
    result.m = 4096;
    result.k = .1;
    result.base_sha256 = original.at("base_sha256").get<std::string>();
    result.effective_sha256 = original.at("effective_sha256").get<std::string>();
    result.rhs_sha256 = original.at("rhs_sha256").get<std::string>();
    if (original.at("id") != "uvlm-column-update-result" ||
        original.at("sha256") !=
            "ddd7649056e4253e5cafdf993322bbe1f8bf0c4ed82ec79c1ec7a9fe0ccc2fc9" ||
        result.base_sha256 != "c37f5c1a7d30c15777a84494917796d458b2f802a833de539be73a9f6ccfe5fb" ||
        result.effective_sha256 !=
            "8ab8b7f97bc415b3f69ae4c6a70a4bfcd29a10931d441605302c68bed4988cbd" ||
        result.rhs_sha256 != "de478c4580657fb266760e41169a13270d7f799114ed58f1e7ffba69fa0aad6d")
        throw std::runtime_error("Frozen input identity changed");
    const auto& arrays = data.at("arrays");
    result.xv = fixture_reals(arrays.at("xv"), 4097);
    result.xc = fixture_reals(arrays.at("xc"), 4096);
    result.fc = fixture_reals(arrays.at("fc"), 4096);
    result.fu = fixture_reals(arrays.at("fu"), 4096);
    result.mc0 = fixture_reals(arrays.at("mc0"), 4096);
    result.mu0 = fixture_reals(arrays.at("mu0"), 4096);
    result.update = fixture_complex(arrays.at("v"), std::array<std::size_t, 1>{4096});
    result.rhs = fixture_complex(arrays.at("b"), std::array<std::size_t, 2>{4096, 3});
    result.reference_hi = fixture_complex(arrays.at("hi"), std::array<std::size_t, 2>{4096, 3});
    result.reference_lo = fixture_complex(arrays.at("lo"), std::array<std::size_t, 2>{4096, 3});
    if (arrays.at("b").at("sha256") != result.rhs_sha256)
        throw std::runtime_error("Frozen RHS payload changed");
    const auto& modes = data.at("modes");
    if (!modes.is_array() || modes.size() != 3)
        throw std::runtime_error("Unexpected frozen mode count");
    constexpr std::array names{"heave-quarter", "pitch-quarter", "pitch-mid"};
    for (std::size_t i = 0; i < 3; ++i) {
        auto& mode = result.modes[i];
        const auto& entry = modes[i];
        mode.name = entry.at("name").get<std::string>();
        if (mode.name != names[i] || entry.at("scales").size() != 3)
            throw std::runtime_error("Unexpected frozen mode identity");
        mode.h = fixture_number(entry.at("H"));
        mode.alpha = fixture_number(entry.at("alpha"));
        mode.xp = fixture_number(entry.at("xp"));
        for (std::size_t j = 0; j < 3; ++j) {
            mode.scales[j] = fixture_number(entry.at("scales")[j]);
            if (mode.scales[j] <= 0)
                throw std::runtime_error("Nonpositive frozen mode scale");
        }
    }
    const auto& observations = data.at("observations");
    if (!observations.is_array() || observations.size() != 39)
        throw std::runtime_error("Unexpected frozen observation count");
    for (std::size_t i = 0; i < 39; ++i) {
        auto& target = result.observations[i];
        const auto& entry = observations[i];
        target.id = entry.at("id").get<std::string>();
        target.mode = entry.at("mode").get<std::string>();
        target.kind = entry.at("kind").get<std::string>();
        if (!entry.at("column").is_number_integer())
            throw std::runtime_error("Unexpected observation column type");
        target.column = entry.at("column").get<std::size_t>();
        if (target.column != i / 13 || target.mode != names[i / 13] ||
            target.id != std::string(names[i / 13]) + "/" + std::to_string(i % 13) ||
            target.kind != (i % 13 == 0   ? "force"
                            : i % 13 == 1 ? "moment"
                                          : "pressure"))
            throw std::runtime_error("Frozen observation order changed");
        if (!entry.at("value").is_array() || entry.at("value").size() != 2)
            throw std::runtime_error("Unexpected observation value shape");
        target.value = {fixture_number(entry.at("value")[0]), fixture_number(entry.at("value")[1])};
        target.scale = fixture_number(entry.at("scale"));
        if (target.scale <= 0)
            throw std::runtime_error("Nonpositive observation scale");
        if (target.kind == "pressure") {
            target.lower = fixture_number(entry.at("lower"));
            target.upper = fixture_number(entry.at("upper"));
            if (target.lower < 0 || target.upper > 1 || target.lower >= target.upper)
                throw std::runtime_error("Invalid frozen pressure window");
        }
    }
    return result;
}

inline void fixture_append_word(std::vector<unsigned char>& bytes, double value) {
    const auto bits = std::bit_cast<std::uint64_t>(value);
    for (unsigned i = 0; i < 8; ++i)
        bytes.push_back(static_cast<unsigned char>((bits >> (8U * i)) & 255U));
}
inline std::string base_hash(std::span<const double> base, std::size_t n) {
    boost::hash2::sha2_256 hash;
    std::vector<unsigned char> row;
    row.reserve(n * 8);
    for (std::size_t i = 0; i < n; ++i) {
        row.clear();
        for (std::size_t j = 0; j < n; ++j)
            fixture_append_word(row, base[j * n + i]);
        hash.update(row.data(), row.size());
    }
    return boost::hash2::to_string(hash.result());
}
inline std::string effective_hash(std::span<const double> base,
                                  std::span<const FixtureComplex> update, std::size_t n) {
    boost::hash2::sha2_256 hash;
    std::vector<unsigned char> row;
    row.reserve(n * 16);
    for (std::size_t i = 0; i < n; ++i) {
        row.clear();
        for (std::size_t j = 0; j < n; ++j) {
            const auto value = j + 1 == n ? FixtureComplex{base[j * n + i], 0} + update[i]
                                          : FixtureComplex{base[j * n + i], 0};
            fixture_append_word(row, value.real());
            fixture_append_word(row, value.imag());
        }
        hash.update(row.data(), row.size());
    }
    return boost::hash2::to_string(hash.result());
}
inline std::vector<double> reconstruct_base(const Fixture& fixture) {
    if (fixture.m != 4096 || fixture.xv.size() != 4097 || fixture.xc.size() != 4096 ||
        fixture.update.size() != 4096 || std::fegetround() != FE_TONEAREST)
        throw std::runtime_error("Invalid coefficient reconstruction precondition");
    constexpr std::size_t n = 4096;
    constexpr double pi = 0x1.921fb54442d18p+1;
    const auto signed_h = [](double s, double x) {
        const double r = s - x;
        if (r == 0)
            throw std::runtime_error("Fixture observation coincides with vortex");
        // Match the frozen numpy ufunc order, with contraction disabled by the
        // target policy. Stored coordinates avoid platform sin/cos differences.
        const double squared = r * r;
        const double root = std::sqrt(squared + 25000000.);
        const double quotient = root / r;
        const double numerator = quotient - 1.;
        const double denominator = (2. * pi) * 5000.;
        return numerator / denominator;
    };
    std::vector<double> result(n * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j + 1 < n; ++j)
            result[j * n + i] =
                signed_h(fixture.xv[j], fixture.xc[i]) - signed_h(fixture.xv[j + 1], fixture.xc[i]);
        result[(n - 1) * n + i] =
            signed_h(fixture.xv[n - 1], fixture.xc[i]) - signed_h(161., fixture.xc[i]);
    }
    if (base_hash(result, n) != fixture.base_sha256 ||
        effective_hash(result, fixture.update, n) != fixture.effective_sha256)
        throw std::runtime_error("Reconstructed coefficients differ from frozen input hashes");
    return result;
}
} // namespace uvlas::test_support::full_precision
