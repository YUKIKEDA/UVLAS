#pragma once

#include "case_wing_geometry.hpp"

namespace uvlas::coupling::detail::wing_geometry {

struct Patch {
    std::optional<wing_chart::Chart> chart;
    std::optional<CaseMotionSnapshot> motion;
    std::size_t points{}, cells{}, edges{}, bands{};
    [[nodiscard]] std::size_t m() const {
        return chart->chord.size() - 1;
    }
    [[nodiscard]] std::size_t n() const {
        return chart->span.size() - 1;
    }
    [[nodiscard]] std::size_t point(std::size_t i, std::size_t j) const {
        return points + i * (n() + 1) + j;
    }
    [[nodiscard]] std::size_t chord_edge(std::size_t i, std::size_t j) const {
        return edges + i * (n() + 1) + j;
    }
    [[nodiscard]] std::size_t span_edge(std::size_t i, std::size_t j) const {
        return edges + m() * (n() + 1) + i * n() + j;
    }
    [[nodiscard]] std::size_t diagonal(std::size_t i, std::size_t j) const {
        return edges + m() * (n() + 1) + (m() + 1) * n() + i * n() + j;
    }
};

[[nodiscard]] std::size_t root(std::span<std::size_t>, std::size_t, wing_chart::Math&);
void share(const Patch&, std::uint64_t own_boundary, const Patch&, std::uint64_t other_boundary,
           bool reverse, std::span<const CaseWingPoint>, std::span<std::size_t> point_owners,
           std::span<std::size_t> edge_owners, std::span<bool> internal, wing_chart::Math&,
           CaseWingDiagnostic&);

} // namespace uvlas::coupling::detail::wing_geometry
