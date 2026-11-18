#include "exact_polynomial_detail.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <uvlas/structure/pg_coefficients.hpp>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using Json = nlohmann::json;

// This is an explicit admission/ledger cap, not a reservation or an RSS bound.
// The fixture is small; the whole-call exact planner accounts cumulative scratch.
constexpr auto reference_budget = std::uint64_t{64} << 30U;
constexpr auto external_budget = std::uint64_t{16} << 30U;
constexpr std::size_t reference_case_count = 6;
constexpr std::size_t selections_per_case = 166;

[[nodiscard]] PgCoefficientLimits reference_limits() {
    return {{64, 1000, 4096, 32, 32, 65536, 2000000, ExactEstimateAdmission{external_budget},
             ExactTimingMode::Disabled},
            1000,
            20000000,
            2000,
            {}};
}

template <class T>
[[nodiscard]] T take_reference(core::Result<T> result) {
    if (!result) {
        throw std::runtime_error("Unexpected coefficient reference failure " +
                                 std::to_string(static_cast<unsigned>(result.error().code)));
    }
    return std::move(*result);
}

[[nodiscard]] core::Index reference_index(const Json& value) {
    REQUIRE(value.is_number_integer());
    if (value.is_number_unsigned()) {
        const auto integer = value.get<std::uint64_t>();
        REQUIRE(integer <= static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()));
        return static_cast<core::Index>(integer);
    }
    return value.get<core::Index>();
}

[[nodiscard]] ExactRational input_fraction(const Json& pair, const ExactLimits& limits,
                                           core::ResourcePlan& plan) {
    REQUIRE(pair.is_array());
    REQUIRE(pair.size() == 2);
    REQUIRE(pair[0].is_string());
    REQUIRE(pair[1].is_string());
    return take_reference(ExactRational::from_integer_pair(pair[0].get_ref<const std::string&>(),
                                                           pair[1].get_ref<const std::string&>(),
                                                           limits, plan));
}

[[nodiscard]] PgGraph reference_graph(const Json& entry, const ExactLimits& limits,
                                      core::ResourcePlan& plan) {
    const auto x0 = input_fraction(entry.at("x0"), limits, plan);
    const auto width = input_fraction(entry.at("width"), limits, plan);
    const auto skew = input_fraction(entry.at("skew"), limits, plan);
    const auto& height = entry.at("height");
    REQUIRE(height.is_array());
    std::vector<ExactRational> coefficients;
    coefficients.reserve(height.size());
    for (const auto& term : height) {
        REQUIRE(term.is_array());
        REQUIRE(term.size() == 4);
        REQUIRE(term[2].is_string());
        REQUIRE(term[3].is_string());
        coefficients.push_back(take_reference(
            ExactRational::from_integer_pair(term[2].get_ref<const std::string&>(),
                                             term[3].get_ref<const std::string&>(), limits, plan)));
    }
    // Build references only after the owning vector has reached its final size.
    std::vector<ExactInputTerm> terms;
    terms.reserve(height.size());
    for (std::size_t i = 0; i < height.size(); ++i) {
        terms.push_back(
            {{reference_index(height[i][0]), reference_index(height[i][1])}, coefficients[i]});
    }
    const auto polynomial = take_reference(ExactChartPolynomial::make(terms, limits, plan));
    const PgGraphLimits graph_limits{limits, 64, {}};
    return take_reference(PgGraph::prepare(
        {x0, width, skew, polynomial, reference_index(entry.at("n")), PgGraphConvention::Version1},
        graph_limits, plan));
}

void exact_reference_terms(const ExactCellPolynomial& actual, const Json& expected) {
    REQUIRE(expected.is_array());
    const auto terms = actual.terms();
    REQUIRE(terms.size() == expected.size());
    for (std::size_t i = 0; i < terms.size(); ++i) {
        CAPTURE(i);
        const auto& literal = expected[i];
        REQUIRE(literal.is_array());
        REQUIRE(literal.size() == 4);
        CHECK(terms[i].powers.first == reference_index(literal[0]));
        CHECK(terms[i].powers.second == reference_index(literal[1]));
        REQUIRE(literal[2].is_string());
        REQUIRE(literal[3].is_string());
        // Read canonical p/q directly. Do not regenerate expected polynomials or
        // normalize the oracle through the product's own rational arithmetic.
        const auto& rational = detail::ExactFactory::value(terms[i].coefficient).backend();
        CHECK(rational.num().str(0, std::ios_base::dec) ==
              literal[2].get_ref<const std::string&>());
        CHECK(rational.denom().str(0, std::ios_base::dec) ==
              literal[3].get_ref<const std::string&>());
    }
}

[[nodiscard]] PgCoefficientField reference_field(const Json& value) {
    REQUIRE(value.is_string());
    const auto& field = value.get_ref<const std::string&>();
    if (field == "R") {
        return PgCoefficientField::R;
    }
    if (field == "P") {
        return PgCoefficientField::P;
    }
    REQUIRE(field == "C");
    return PgCoefficientField::C;
}

void anchor_reference(const std::span<const PgCoefficientSource> columns, std::size_t& next,
                      const Json& expected, PgCoefficientKind kind) {
    REQUIRE(expected.is_array());
    for (const auto& original : expected) {
        REQUIRE(next < columns.size());
        CAPTURE(next);
        CHECK(columns[next].kind == kind);
        CHECK(columns[next].original.value == reference_index(original));
        ++next;
    }
}
} // namespace

TEST_CASE("PG complete dictionaries match the independent canonical coefficient fixture",
          "[SH-D03][SH-D05][SH-D07][N4-02][N4-08][verification]") {
    std::ifstream stream{UVLAS_PG_COEFFICIENTS_FIXTURE};
    REQUIRE(stream.good());
    const auto fixture = Json::parse(stream);
    REQUIRE(fixture.at("schema") == "uvlas.pg-coefficient-reference/1");
    const auto& cases = fixture.at("cases");
    REQUIRE(cases.is_array());
    REQUIRE(cases.size() == reference_case_count);
    constexpr std::array expected_ids{"flat-n1",          "mixed-quartic-n1", "flat-n2",
                                      "mixed-quartic-n2", "flat-n4",          "mixed-quartic-n4"};
    std::size_t checked_selections = 0;
    for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
        const auto& entry = cases[case_index];
        REQUIRE(entry.at("id").is_string());
        REQUIRE(entry.at("id").get_ref<const std::string&>() == expected_ids[case_index]);
        CAPTURE(case_index, expected_ids[case_index]);
        core::ResourcePlan plan{reference_budget};
        const auto limits = reference_limits();
        const auto graph = reference_graph(entry, limits.algebra, plan);
        auto workspace = take_reference(core::Workspace::create(0, plan));
        // Exactly one complete preparation per case; all selections share it.
        const auto dictionary =
            take_reference(PgCoefficientDictionary::prepare(graph, limits, *workspace, plan));
        const auto original = take_reference(dictionary.original_columns());
        REQUIRE(original == reference_index(entry.at("original_columns")));
        REQUIRE(original >= 0);
        const auto columns = take_reference(dictionary.columns());
        const auto raw_count = static_cast<std::size_t>(original);
        const auto& anchor_r = entry.at("anchor_r_sources");
        const auto& anchor_p = entry.at("anchor_p_sources");
        REQUIRE(anchor_r.is_array());
        REQUIRE(anchor_p.is_array());
        REQUIRE(columns.size() == raw_count + anchor_r.size() + anchor_p.size());
        for (std::size_t i = 0; i < raw_count; ++i) {
            CHECK(columns[i].kind != PgCoefficientKind::AnchorR);
            CHECK(columns[i].kind != PgCoefficientKind::AnchorP);
        }
        std::size_t next = raw_count;
        anchor_reference(columns, next, anchor_r, PgCoefficientKind::AnchorR);
        anchor_reference(columns, next, anchor_p, PgCoefficientKind::AnchorP);
        REQUIRE(next == columns.size());

        const auto& selected = entry.at("selected");
        REQUIRE(selected.is_array());
        REQUIRE(selected.size() == selections_per_case);
        for (const auto& expected : selected) {
            const auto column = reference_index(expected.at("column"));
            const auto& cell_value = expected.at("cell");
            REQUIRE(cell_value.is_array());
            REQUIRE(cell_value.size() == 2);
            const PgCellIndex cell{reference_index(cell_value[0]), reference_index(cell_value[1])};
            CAPTURE(column, cell.first, cell.second, expected.at("field"));
            const auto& field = expected.at("field");
            const auto actual =
                field == "g" ? dictionary.g(cell)
                             : dictionary.coefficient({column}, cell, reference_field(field));
            REQUIRE(actual.has_value());
            REQUIRE(*actual != nullptr);
            exact_reference_terms(**actual, expected.at("terms"));
            ++checked_selections;
        }
    }
    CHECK(checked_selections == 996);
}
