#pragma once
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <stop_token>
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
/// Pass the source by value so a token keeps the requested state on either STL.
[[nodiscard]] inline std::stop_token requested_stop_token(std::stop_source source = {}) {
    REQUIRE(source.request_stop());
    return source.get_token();
}

struct Fixture {
    core::ResourcePlan input_plan{cap};
    std::optional<PgGraph> graph;
    std::vector<PgPartitionPoint> domain, thickness, named, cut_points;
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
    [[nodiscard]] PgGraph& chart_graph() & {
        if (!graph)
            throw std::runtime_error("Partition fixture graph is absent");
        return *graph;
    }
    [[nodiscard]] PgPartitionRequest request() const {
        if (!graph)
            throw std::runtime_error("Partition fixture graph is absent");
        return {*graph, domain, {thickness, triangles, heights}, layers, names, refinements};
    }
};
/// Published incidence input; f owns every borrowed endpoint.
inline void incidence_fixture(Fixture& f) {
    f.domain.clear();
    for (const auto p : std::array<std::array<int, 2>, 8>{
             {{0, 0}, {8, 0}, {8, 1}, {1, 1}, {1, 3}, {8, 3}, {8, 4}, {0, 4}}})
        f.domain.push_back({f.r(std::to_string(p[0]), "8"), f.r(std::to_string(p[1]), "8")});
    for (std::size_t i = 1; i < f.heights.size(); ++i)
        f.heights[i] = f.r("2");
    f.named.push_back({f.r("1", "4"), f.r("0")});
    f.named.push_back({f.r("3", "4"), f.r("0")});
    f.names.push_back({"partial", f.named});
    f.cut_points.push_back({f.r("1", "10"), f.r("4", "5")});
    f.cut_points.push_back({f.r("9", "10"), f.r("1", "5")});
    f.refinements.push_back({23, 4, f.cut_points[0], f.cut_points[1]});
}

} // namespace partition_test
