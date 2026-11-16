#include "full_precision_fixture.hpp"
#include "full_precision_oracle.hpp"

#include <catch2/catch_test_macros.hpp>

#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace reference = uvlas::test_support::full_precision;

TEST_CASE("Full precision oracle checks independent small exact systems",
          "[N4-05][N4-06][N4-10][verification]") {
    CHECK_NOTHROW(reference::self_test());
    std::ifstream input(UVLAS_FULL_PRECISION_FIXTURE, std::ios::binary);
    REQUIRE(input.good());
    const std::string bytes{std::istreambuf_iterator<char>{input},
                            std::istreambuf_iterator<char>{}};
    CHECK_NOTHROW(reference::checked_fixture_json(bytes));
    const auto original = reference::checked_fixture_json(bytes);
    auto altered = original;
    altered["observations"][0]["scale"] = 2.0;
    CHECK_THROWS_AS(reference::checked_fixture_json(altered.dump()), std::runtime_error);
    altered = original;
    altered["observations"][0]["value"][0] = 2.0;
    CHECK_THROWS_AS(reference::checked_fixture_json(altered.dump()), std::runtime_error);
    for (const auto* name : {"fc", "lo"}) {
        altered = original;
        auto payload = altered["arrays"][name]["payload_hex"].get<std::string>();
        payload.front() = payload.front() == '0' ? '1' : '0';
        altered["arrays"][name]["payload_hex"] = payload;
        // Even a replacement self-declared payload hash cannot rebind the
        // independently fixed digest of the complete fixture document.
        altered["arrays"][name]["sha256"] = std::string(64, '0');
        CHECK_THROWS_AS(reference::checked_fixture_json(altered.dump()), std::runtime_error);
    }
}
