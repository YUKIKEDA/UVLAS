#pragma once
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <uvlas/structure/pg_integration_partition.hpp>
#include <vector>

namespace partition_test {
using namespace uvlas;
using namespace uvlas::structure;
inline constexpr auto cap = std::uint64_t{1} << 60U;
inline ExactLimits exact_limits() {
    return {10000,
            1000,
            10000,
            32,
            32,
            1000000000,
            2000000000,
            ExactEstimateAdmission{cap},
            ExactTimingMode::Disabled};
}
inline PgPartitionLimits limits() {
    return {exact_limits(), 1000000, 1000000, 1000000, 10000000, 10000000,
            10000000,       1000000, cap,     cap,     1000000,  {}};
}
template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::runtime_error("Unexpected partition failure " +
                                 std::to_string(static_cast<unsigned>(value.error().code)));
    return std::move(*value);
}
struct Fixture {
    core::ResourcePlan input_plan{cap};
    std::optional<PgGraph> graph;
    std::vector<PgPartitionPoint> domain, thickness, named;
    std::vector<ExactRational> heights;
    std::vector<std::array<core::Index, 3>> triangles;
    std::array<ExactRational, 3> zeta;
    std::vector<PgPartitionLayer> layers;
    std::vector<PgPartitionBoundary> names;
    std::vector<PgPartitionRefinement> refinements;
    explicit Fixture(core::Index n = 1) {
        const auto zero = r("0"), one = r("1");
        const ExactChartPolynomial f;
        graph = take(PgGraph::prepare({zero, one, zero, f, n, PgGraphConvention::Version1},
                                      {exact_limits(), 64, {}}, input_plan));
        for (const auto& pair : std::array<std::array<std::string_view, 2>, 4>{
                 {{"0", "0"}, {"1", "0"}, {"1", "1"}, {"0", "1"}}}) {
            domain.push_back(point(pair[0], pair[1]));
            thickness.push_back(point(pair[0], pair[1]));
            heights.push_back(r("1"));
        }
        triangles = {{0, 1, 2}, {0, 2, 3}};
        zeta[0] = r("-1");
        zeta[1] = r("-1", "3");
        zeta[2] = r("1");
        layers.push_back({core::LayerId{11}, zeta[0], zeta[1]});
        layers.push_back({core::LayerId{29}, zeta[1], zeta[2]});
    }
    ExactRational r(std::string_view numerator, std::string_view denominator = "1") {
        return take(
            ExactRational::from_integer_pair(numerator, denominator, exact_limits(), input_plan));
    }
    PgPartitionPoint point(std::string_view s, std::string_view t) {
        return {r(s), r(t)};
    }
    [[nodiscard]] PgPartitionRequest request() const {
        return {*graph, domain, {thickness, triangles, heights}, layers, names, refinements};
    }
};
} // namespace partition_test
