#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <type_traits>
#include <uvlas/core/checked_dimensions.hpp>

namespace {
using namespace uvlas::core;
constexpr Operation operation = Operation::PrepareModel;

void check_failure(const auto& result, FailureCode code) {
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == code);
    CHECK(result.error().operation == operation);
}
} // namespace

TEST_CASE("Extent products preserve exact integers and empty arrays", "[unit][N4-02]") {
    STATIC_REQUIRE(std::is_same_v<Index, std::int64_t>);
    STATIC_REQUIRE(sizeof(Index) == 8);
    CHECK(checked_extent_product(std::span<const Index>{}, operation) == Index{1});
    CHECK(checked_extent_product(std::array<Index, 3>{2, 3, 7}, operation) == Index{42});
    CHECK(checked_extent_product(std::array<Index, 2>{0, 0}, operation) == Index{0});
    const Index maximum = std::numeric_limits<Index>::max();
    CHECK(checked_extent_product(std::array{maximum, Index{1}}, operation) == maximum);
    // The mathematical product is zero regardless of where the zero occurs.
    for (const auto& shape :
         {std::array{Index{0}, maximum, maximum}, std::array{maximum, Index{0}, maximum},
          std::array{maximum, maximum, Index{0}}}) {
        CHECK(checked_extent_product(shape, operation) == Index{0});
    }
}

TEST_CASE("Negative dimensions cannot hide behind zero and overflow is rejected",
          "[negative][N4-02]") {
    const Index maximum = std::numeric_limits<Index>::max();
    for (const auto& shape : {std::array{Index{-1}, Index{0}}, std::array{Index{0}, Index{-1}},
                              std::array{std::numeric_limits<Index>::min(), maximum}}) {
        check_failure(checked_extent_product(shape, operation), FailureCode::InvalidInput);
    }
    for (const auto& shape : {std::array{maximum, Index{2}}, std::array{Index{2}, maximum},
                              std::array{Index{3037000500}, Index{3037000500}}}) {
        check_failure(checked_extent_product(shape, operation), FailureCode::ResourceExceeded);
    }
    CHECK(checked_extent_product(std::array{Index{3037000499}, Index{3037000499}}, operation) ==
          Index{9223372030926249001});
}

TEST_CASE("Payload byte counts check the host allocation range", "[negative][N4-02]") {
    CHECK(checked_byte_count(0, sizeof(double), operation) == std::size_t{0});
    CHECK(checked_byte_count(3, sizeof(double), operation) == std::size_t{24});
    check_failure(checked_byte_count(-1, sizeof(double), operation), FailureCode::InvalidInput);
    check_failure(checked_byte_count(0, 0, operation), FailureCode::InvalidInput);
    const std::size_t maximum = std::numeric_limits<std::size_t>::max();
    CHECK(checked_byte_count(1, maximum, operation) == maximum);
    check_failure(checked_byte_count(2, maximum, operation), FailureCode::ResourceExceeded);
    // These calls only compute counts, never request large memory from the OS.
    const Index largest_double_count = static_cast<Index>(maximum / sizeof(double));
    CHECK(checked_byte_count(largest_double_count, sizeof(double), operation) ==
          maximum - (maximum % sizeof(double)));
    check_failure(checked_byte_count(largest_double_count + 1, sizeof(double), operation),
                  FailureCode::ResourceExceeded);
}

TEST_CASE("External signed 32-bit arguments are checked before narrowing", "[negative][N4-02]") {
    constexpr Index maximum = 2147483647;
    CHECK(checked_api32(0, operation) == std::int32_t{0});
    CHECK(checked_api32(maximum, operation) == std::int32_t{2147483647});
    check_failure(checked_api32(-1, operation), FailureCode::InvalidInput);
    check_failure(checked_api32(maximum + 1, operation), FailureCode::ResourceExceeded);
    CHECK(checked_api32_product(std::array{Index{7}, Index{11}}, operation) == std::int32_t{77});
    CHECK(checked_api32_product(std::array{Index{46340}, Index{46340}}, operation) ==
          std::int32_t{2147395600});
    check_failure(checked_api32_product(std::array{Index{46341}, Index{46341}}, operation),
                  FailureCode::ResourceExceeded);
    check_failure(
        checked_api32_product(std::array{std::numeric_limits<Index>::max(), Index{2}}, operation),
        FailureCode::ResourceExceeded);
    CHECK(checked_api32_product(std::array{maximum + 1, Index{0}}, operation) == std::int32_t{0});
    check_failure(checked_api32_product(std::array{Index{0}, Index{-1}}, operation),
                  FailureCode::InvalidInput);
}
