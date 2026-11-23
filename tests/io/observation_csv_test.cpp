#include "observation_csv.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <charconv>
#include <limits>

namespace {

using uvlas::core::FailureCode;
using uvlas::core::LayerId;
using uvlas::core::ObservationId;
using uvlas::core::Operation;
using uvlas::core::ResourcePlan;
using uvlas::io::ScaledBinary64;
using uvlas::io::ScaledTwoWord;
using uvlas::io::SiEnclosure;
using uvlas::io::detail::encode_observation_csv;
using uvlas::io::detail::ObservationCsvColumn;
using uvlas::io::detail::ObservationCsvLimits;
using uvlas::io::detail::ObservationCsvSample;
using uvlas::io::detail::ObservationSide;

constexpr auto max_id = std::numeric_limits<std::uint64_t>::max();
constexpr ObservationCsvLimits normal_limits{8, 8, 128, 4096};
constexpr std::string_view expected_header =
    "time_hi_s,time_lo_s,checkpoint,segment,observation,component,quantity,unit,basis,"
    "reference_point,layer,side,representation,value_hi,value_lo,radius,exponent2\r\n";

constexpr ObservationCsvColumn plain_column{
    ObservationId{0}, 0, "displacement", "m", "global", "origin", {}, ObservationSide::None};

TEST_CASE("Observation CSV retains the accepted time and each physical component",
          "[unit][verification][IO4-05]") {
    // METADATA contains a multibyte character, a delimiter, a quote and a record newline.
    constexpr std::array columns{
        ObservationCsvColumn{ObservationId{max_id}, 0, "stress", "Pa", "local,\"A\"",
                             "P\n\xe5\xb1\xa4", LayerId{max_id}, ObservationSide::Lower},
        ObservationCsvColumn{ObservationId{max_id}, 0, "stress", "Pa", "local,\"A\"",
                             "P\n\xe5\xb1\xa4", LayerId{max_id}, ObservationSide::Upper},
        plain_column,
    };
    constexpr std::array<SiEnclosure, 3> values{ScaledTwoWord{0.5, -0.125, 0.0625, 2048},
                                                ScaledBinary64{-0.5, 0, -18},
                                                ScaledBinary64{-0.0, 0, 0}};
    const std::array samples{ObservationCsvSample{{2, 0.125}, max_id, max_id - 1, values}};
    ResourcePlan memory{8192};
    const auto result = encode_observation_csv(columns, samples, normal_limits, memory);
    REQUIRE(result);
    const std::string_view bytes{result->values()};
    REQUIRE(bytes.starts_with(expected_header));
    CHECK(bytes.substr(expected_header.size()) ==
          "2,0.125,18446744073709551615,18446744073709551614,18446744073709551615,0,"
          "stress,Pa,\"local,\"\"A\"\"\",\"P\n\xe5\xb1\xa4\",18446744073709551615,lower,two_word,"
          "0.5,-0.125,0.0625,2048\r\n"
          "2,0.125,18446744073709551615,18446744073709551614,18446744073709551615,0,"
          "stress,Pa,\"local,\"\"A\"\"\",\"P\n\xe5\xb1\xa4\",18446744073709551615,upper,binary64,"
          "-0.5,,0,-18\r\n"
          "2,0.125,18446744073709551615,18446744073709551614,0,0,displacement,m,global,origin,"
          ",none,binary64,-0,,0,0\r\n");
}

TEST_CASE("Observation CSV recovers binary64 endpoint words without narrowing",
          "[property][verification][IO4-05]") {
    const std::array columns{plain_column};
    constexpr auto min_word = std::numeric_limits<double>::denorm_min();
    constexpr auto max_word = std::numeric_limits<double>::max();
    const std::array<SiEnclosure, 1> values{ScaledTwoWord{max_word, -min_word, min_word, 0}};
    const std::array samples{ObservationCsvSample{{1, min_word}, 0, 0, values}};
    ResourcePlan memory{8192};
    const auto result = encode_observation_csv(columns, samples, normal_limits, memory);
    REQUIRE(result);
    auto field = std::string_view{result->values()}.substr(expected_header.size());
    const auto word_at = [&](std::size_t index) {
        auto tail = field;
        for (std::size_t i = 0; i < index; ++i) {
            const auto delimiter = tail.find(',');
            REQUIRE(delimiter != tail.npos);
            tail.remove_prefix(delimiter + 1);
        }
        const auto token = tail.substr(0, tail.find(','));
        double word{};
        const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), word);
        REQUIRE(error == std::errc{});
        REQUIRE(end == token.data() + token.size());
        return word;
    };
    CHECK(word_at(1) == min_word);
    CHECK(word_at(13) == max_word);
    CHECK(word_at(14) == -min_word);
    CHECK(word_at(15) == min_word);
}

TEST_CASE("Observation CSV rejects malformed metadata, samples and numerical words",
          "[negative][verification][IO4-05]") {
    std::array columns{plain_column};
    std::array<SiEnclosure, 1> values{ScaledTwoWord{0.5, 0, 0, 0}};
    std::array samples{ObservationCsvSample{{1, 0}, 0, 0, values}};
    std::span<const ObservationCsvColumn> column_view = columns;
    std::span<const ObservationCsvSample> sample_view = samples;
    const std::array duplicates{plain_column, plain_column};
    const std::array repeated_sample{samples[0], samples[0]};
    auto code = FailureCode::InvalidInput;
    SECTION("duplicate component identity") {
        column_view = duplicates;
    }
    SECTION("non-increasing accepted checkpoint") {
        sample_view = repeated_sample;
    }
    SECTION("incomplete sample") {
        samples[0].values = {};
    }
    SECTION("history without observable columns") {
        column_view = {};
        samples[0].values = {};
    }
    SECTION("invalid layer association") {
        columns[0].side = ObservationSide::Upper;
    }
    SECTION("invalid UTF-8") {
        columns[0].quantity = "\xed\xa0\x80";
    }
    SECTION("embedded NUL") {
        columns[0].unit = std::string_view{"m\0s", 3};
    }
    SECTION("nonfinite time word") {
        samples[0].time_seconds.lo = std::numeric_limits<double>::infinity();
        code = FailureCode::Nonfinite;
    }
    SECTION("nonfinite measurement word") {
        std::get<ScaledTwoWord>(values[0]).low = std::numeric_limits<double>::quiet_NaN();
        code = FailureCode::Nonfinite;
    }
    SECTION("negative outward radius") {
        values[0] = ScaledBinary64{0.5, -1, 0};
    }
    ResourcePlan memory{8192};
    const auto result = encode_observation_csv(column_view, sample_view, normal_limits, memory);
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    CHECK(result.error().operation == Operation::WriteResults);
}

TEST_CASE("Observation CSV admits exact budgets and returns no bytes on refusal",
          "[negative][verification][IO4-05][N4-02]") {
    const std::array columns{plain_column};
    const std::array<SiEnclosure, 1> values{ScaledBinary64{0.5, 0, 0}};
    const std::array samples{ObservationCsvSample{{1, 0}, 0, 0, values}};
    ResourcePlan baseline_memory{8192};
    const auto baseline = encode_observation_csv(columns, samples, normal_limits, baseline_memory);
    REQUIRE(baseline);
    auto budget = baseline_memory.allocated_bytes();
    ObservationCsvLimits limits{1, 1, plain_column.quantity.size(), baseline->size()};
    std::stop_source cancel;
    auto code = FailureCode::ResourceExceeded;
    bool succeeds = false;
    SECTION("exact admission") {
        succeeds = true;
    }
    SECTION("column count") {
        --limits.maximum_columns;
    }
    SECTION("sample count") {
        --limits.maximum_samples;
    }
    SECTION("label bytes") {
        --limits.maximum_label_bytes;
    }
    SECTION("output bytes") {
        --limits.maximum_output_bytes;
    }
    SECTION("shared allocation budget") {
        --budget;
    }
    SECTION("cancelled operation") {
        REQUIRE(cancel.request_stop());
        code = FailureCode::Cancelled;
    }
    ResourcePlan memory{budget};
    const auto result =
        encode_observation_csv(columns, samples, limits, memory, cancel.get_token());
    if (succeeds) {
        REQUIRE(result);
        CHECK(std::string_view{result->values()} == std::string_view{baseline->values()});
    } else {
        REQUIRE_FALSE(result);
        CHECK(result.error().code == code);
        CHECK(result.error().operation == Operation::WriteResults);
        if (code == FailureCode::Cancelled)
            CHECK(memory.allocated_bytes() == 0);
    }
}

TEST_CASE("Observation CSV has an empty history header", "[unit][verification][IO4-05]") {
    ResourcePlan memory{expected_header.size()};
    const auto result = encode_observation_csv({}, {}, {0, 0, 0, expected_header.size()}, memory);
    REQUIRE(result);
    CHECK(std::string_view{result->values()} == expected_header);
}

} // namespace
