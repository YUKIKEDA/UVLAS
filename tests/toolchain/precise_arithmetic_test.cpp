#include <boost/hash2/sha2.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cstdint>
#include <string_view>

TEST_CASE("multiprecision retains a unit lost by binary64 at two to the power 200",
          "[toolchain][N4-01]") {
    using Big = boost::multiprecision::cpp_bin_float_100;
    const Big large = ldexp(Big(1), 200);
    const Big retained = (large + 1) - large;
    CHECK(retained == 1);
    const double rounded = (large + 1).convert_to<double>();
    // This is an exact rounding/representation check, not a tolerance test.
    CHECK(std::bit_cast<std::uint64_t>(rounded) ==
          std::bit_cast<std::uint64_t>(large.convert_to<double>()));
}

TEST_CASE("exact rational arithmetic preserves one third plus one sixth", "[toolchain][N4-01]") {
    using Rational = boost::multiprecision::cpp_rational;
    const Rational sum = Rational(1) / 3 + Rational(1) / 6;
    CHECK(sum == Rational(1) / 2);
    CHECK(sum * 6 == 3);
}

TEST_CASE("Hash2 SHA256 matches fixed standard vectors with incremental input",
          "[toolchain][N4-01]") {
    // SHA-256 known-answer vectors; abc and the 56-byte NIST example.
    constexpr std::array<std::string_view, 3> inputs{
        "", "abc", "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"};
    constexpr std::array<std::string_view, 3> expected{
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"};
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        boost::hash2::sha2_256 hash;
        for (const char character : inputs[i]) {
            hash.update(&character, 1);
        }
        CHECK(boost::hash2::to_string(hash.result()) == expected[i]);
    }
}
