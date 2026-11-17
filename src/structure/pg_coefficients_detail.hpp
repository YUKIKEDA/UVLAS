#pragma once

#include "exact_polynomial_detail.hpp"
#include "pg_graph_detail.hpp"

#include <array>
#include <optional>
#include <uvlas/structure/pg_coefficients.hpp>

namespace uvlas::structure::detail {
struct PgCoefficientData {
    std::optional<PgGraph> source;
    core::Index n = 0;
    core::Index original = 0;
    std::size_t column_count = 0, record_count = 0, required_count = 0, c_map_count = 0;
    core::PlannedBuffer<PgCoefficientSource> columns;
    core::PlannedBuffer<PgCoefficientRecord> records;
    core::PlannedBuffer<ExactCellPolynomial> denominator;
    core::PlannedBuffer<PgRequiredMapEntry> required;
    core::PlannedBuffer<PgAnchorCMapEntry> c_map;
    ExactCellPolynomial zero;
};
struct PgCoefficientFactory {
    static const PgCoefficientData* data(const PgCoefficientDictionary&) noexcept;
    static PgCoefficientDictionary make(core::PlannedBuffer<PgCoefficientData>,
                                        PgCoefficientPreparation) noexcept;
};
/// One complete factory admission, distinct from observed successful payload.
struct PgCoefficientPlan {
    PgCoefficientPreparation report;
    std::uint64_t cells = 0;
    std::uint64_t raw_columns = 0;
    std::uint64_t full_columns = 0;
    std::uint64_t record_capacity = 0;
    std::uint64_t map_capacity = 0;
    std::uint64_t c_map_capacity = 0;
};
// Source/planner audit units (not independent caps or dummy reservations):
// packet: one Data, full-column provenance, sparse record-capacity metadata,
//         one denominator polynomial per cell, required/C-map entry capacities.
// raw: graph clone; cell geometry; univariate bases; moment active/endpoint rows;
//      basic tensor products and every (j,axis,origin,k) lifting scalar column.
// anchor: exact first-s-cell restriction, distinct nonzero R/P extensions;
//         original zero columns remain in provenance; C adds map entries only.
// maps: 30 mandatory + 3 translations, exact reconstruction over all raw columns.
// Each exact leaf keeps its original ResourcePlan charge and admission. The
// caller aggregates their cumulative payload and maximum overlapping capacity.
} // namespace uvlas::structure::detail
