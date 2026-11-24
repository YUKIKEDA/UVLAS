#include "nastran_number_fields.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <limits>
#include <string>
#include <utility>

namespace rd = uvlas::io::detail;
namespace core = uvlas::core;
using namespace std::string_view_literals;

TEST_CASE("MSC Bulk Data real fields preserve the exact decimal before conversion",
          "[IO4-03][unit][negative]") {
    constexpr std::array controls{
        std::pair{"7."sv, "7.0"sv},
        std::pair{"   +.7D+1  "sv, "0.7e+1"sv},
        std::pair{"-003.00-02"sv, "-3.00e-02"sv},
        std::pair{"0.125+d2"sv, ""sv},
        std::pair{".0002+13"sv, "0.0002e+13"sv},
        std::pair{"0002.d-3"sv, "2.0e-3"sv},
        std::pair{"-0.E0"sv, "-0.0e0"sv},
        std::pair{"1.1234567890123456789012345D1200"sv, "1.1234567890123456789012345e1200"sv},
    };
    for (const auto& [source, expected] : controls) {
        CAPTURE(source);
        const auto value = rd::read_nastran_real_field(source);
        if (expected.empty()) {
            CHECK_FALSE(value);
        } else if (value) {
            CHECK(value->decimal() == expected);
        } else {
            FAIL_CHECK("A real field was rejected");
        }
    }
    for (const auto source : {""sv, "   "sv, "7"sv, "-."sv, "1. E7"sv, "1.D"sv, "1.--2"sv, "nan"sv,
                              "2.0E+1junk"sv, "\t2.0"sv, "2.0\0"sv}) {
        CAPTURE(source);
        CHECK_FALSE(rd::read_nastran_real_field(source));
    }
    const auto widest = std::string(76, '8') + ".+10";
    const auto admitted = rd::read_nastran_real_field(widest);
    if (admitted) {
        CHECK(admitted->decimal() == std::string(76, '8') + ".0e+10");
    } else {
        FAIL_CHECK("The 80-byte field should fit");
    }
    const auto too_long = rd::read_nastran_real_field(widest + "0");
    if (too_long) {
        FAIL_CHECK("The field exceeded 80 bytes");
    } else {
        CHECK(too_long.error().code == core::FailureCode::ResourceExceeded);
    }
}

TEST_CASE("MSC Bulk Data integers require a complete signed integer field",
          "[IO4-03][unit][negative]") {
    CHECK(rd::read_nastran_integer_field("  +00027 ") == 27);
    CHECK(rd::read_nastran_integer_field("-9223372036854775808") ==
          std::numeric_limits<std::int64_t>::min());
    CHECK(rd::read_nastran_integer_field("9223372036854775807") ==
          std::numeric_limits<std::int64_t>::max());
    for (const auto source : {""sv, "+"sv, "+-1"sv, "++1"sv, "1.0"sv, "1-2"sv, "1e0"sv, "1 2"sv,
                              "9223372036854775808"sv, "-9223372036854775809"sv}) {
        CAPTURE(source);
        CHECK_FALSE(rd::read_nastran_integer_field(source));
    }
}
