#include "pg_reference_domain.hpp"
#include "pg_trace_polynomial.hpp"

#include <bit>
#include <uvlas/structure/pg_full_trace.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;

struct TracePatch {
    const PgPreparedPatch& patch;
    const PgRankPreparation& rank;
    const PgIntegrationPartition& partition;
    PgReferenceFrame frame;
    PgGraphView graph;
    std::uint64_t columns = 0;
};

[[nodiscard]] TracePatch trace_patch(const PgPreparedPatch& patch) {
    const auto* basis = partition_take(patch.conditioning().reference_basis());
    const auto source = partition_take(basis->source());
    return {patch,
            *partition_take(basis->rank()),
            source.partition,
            partition_take(basis->frame()),
            partition_take(source.graph.view()),
            partition_take(basis->size())};
}

void validate_spans(ReferenceArithmetic& math, std::span<const PgTraceEdgeSpan> spans) {
    if (spans.empty())
        partition_fail(core::FailureCode::InvalidInput);
    const auto one = math.constant("1");
    const ExactRational* previous = nullptr;
    for (const auto& span : spans) {
        math.visit(1);
        if ((previous != nullptr && math.compare(*previous, span.first) != 0) ||
            math.compare(span.first, span.last) >= 0 ||
            math.compare(span.along_first, span.along_last) == 0)
            partition_fail(core::FailureCode::InvalidInput);
        check_fraction(math, span.along_first, one);
        check_fraction(math, span.along_last, one);
        previous = &span.last;
    }
}

struct RestrictedTrace {
    PgCellIndex cell;
    PgTracePolynomial s, t;
    std::array<PgTracePolynomial, 3> midsurface, normal;
    PgTraceDepth depth;
    bool with_outward_boundary_order = false;
};

[[nodiscard]] PgTraceDepth trace_depth(ReferenceArithmetic& math, const TracePatch& patch,
                                       const PgTraceEdgeSpan& span, core::Index thickness,
                                       const PgPartitionPoint& first,
                                       const PgPartitionPoint& last) {
    const auto plane =
        thickness_plane(math, partition_take(patch.partition.thickness()), thickness);
    const auto two = math.constant("2");
    const auto h1 = math.div(plane.evaluate(math, first.s, first.t), two);
    const auto h2 = math.div(plane.evaluate(math, last.s, last.t), two);
    if (h1.sign() <= 0 || h2.sign() <= 0)
        partition_fail(core::FailureCode::GeometryInvalid);
    if (!span.layer)
        return {math.neg(h1), math.copy(h1), math.neg(h2), math.copy(h2)};
    const auto& layer = reference_at(partition_take(patch.partition.layers()), *span.layer);
    return {math.mul(layer.lower, h1), math.mul(layer.upper, h1), math.mul(layer.lower, h2),
            math.mul(layer.upper, h2)};
}

[[nodiscard]] RestrictedTrace restrict_trace(ReferenceArithmetic& math, const TracePatch& patch,
                                             const PgTraceEdgeSpan& span,
                                             const ExactRational& first,
                                             const ExactRational& last) {
    math.visit(6); // Edge, piece, two vertices, thickness triangle and optional layer.
    const auto& edge = reference_at(partition_take(patch.partition.edges()), span.edge);
    if (!edge.external || edge.left.has_value() == edge.right.has_value())
        partition_fail(core::FailureCode::InvalidInput);
    const auto& piece = reference_at(partition_take(patch.partition.pieces()),
                                     edge.left.value_or(edge.right.value_or(-1)));
    const auto vertices = partition_take(patch.partition.vertices());
    const auto& a = reference_at(vertices, edge.vertices[0]);
    const auto& b = reference_at(vertices, edge.vertices[1]);
    const auto interval = math.sub(span.last, span.first);
    const auto point = [&](const ExactRational& eta) -> PgPartitionPoint {
        const auto fraction = math.div(math.sub(eta, span.first), interval);
        const auto along = interpolate(math, span.along_first, span.along_last, fraction);
        return {interpolate(math, a.s, b.s, along), interpolate(math, a.t, b.t, along)};
    };
    const auto x = point(first), y = point(last);
    auto s = trace_affine(math, x.s, y.s);
    auto t = trace_affine(math, x.t, y.t);
    const auto sp = trace_powers(math, s, 4);
    const auto tp = trace_powers(math, t, 4);
    const auto& graph = patch.graph;
    const auto one = math.constant("1"), inv_w = math.div(one, graph.width);
    const auto shifted_height = trace_restrict(math, graph.height, sp, tp);
    PgTracePolynomial x_chart;
    x_chart.coefficients[0] = math.copy(graph.x0);
    trace_add_scaled(math, x_chart, s, graph.width);
    trace_add_scaled(math, x_chart, t, graph.skew);
    const std::array<const PgTracePolynomial*, 3> chart{&x_chart, &t, &shifted_height};
    const auto fs = trace_restrict(math, graph.height, sp, tp, {1, 0});
    const auto ft = trace_restrict(math, graph.height, sp, tp, {0, 1});
    std::array<PgTracePolynomial, 3> k;
    trace_add_scaled(math, k[0], fs, math.neg(inv_w));
    trace_add_scaled(math, k[1], fs, math.mul(graph.skew, inv_w));
    trace_add_scaled(math, k[1], ft, math.constant("-1"));
    k[2].coefficients[0] = math.copy(one);
    std::array<PgTracePolynomial, 3> midsurface, normal;
    for (std::size_t i = 0; i < midsurface.size(); ++i) {
        midsurface[i].coefficients[0] = math.copy(patch.frame.origin_metres[i]);
        for (std::size_t j = 0; j < chart.size(); ++j) {
            const auto& rotation = patch.frame.rotation[i * chart.size() + j];
            trace_add_scaled(math, midsurface[i], *chart[j],
                             math.mul(patch.frame.length_metres, rotation));
            trace_add_scaled(math, normal[i], k[j], rotation);
        }
    }
    const auto oriented = math.compare(span.along_last, span.along_first) > 0;
    return {piece.cell,
            std::move(s),
            std::move(t),
            std::move(midsurface),
            std::move(normal),
            trace_depth(math, patch, span, piece.thickness_triangle, x, y),
            oriented == edge.left.has_value()};
}

[[nodiscard]] std::size_t compare_geometry(ReferenceArithmetic& math, const RestrictedTrace& a,
                                           const RestrictedTrace& b) {
    if (a.with_outward_boundary_order == b.with_outward_boundary_order)
        partition_fail(core::FailureCode::GeometryInvalid);
    std::size_t component = a.normal.size();
    ExactRational dot;
    for (std::size_t i = 0; i < a.normal.size(); ++i) {
        if (!trace_same(math, a.midsurface[i], b.midsurface[i]))
            partition_fail(core::FailureCode::GeometryInvalid);
        const auto j = (i + 1) % a.normal.size(), k = (i + 2) % a.normal.size();
        if (!trace_same(math, trace_product(math, a.normal[j], b.normal[k]),
                        trace_product(math, a.normal[k], b.normal[j])))
            partition_fail(core::FailureCode::GeometryInvalid);
        if (component == a.normal.size() && !a.normal[i].zero())
            component = i;
        dot = math.add(dot, math.mul(a.normal[i].coefficients[0], b.normal[i].coefficients[0]));
    }
    if (dot.sign() <= 0 || component == a.normal.size())
        partition_fail(core::FailureCode::GeometryInvalid);
    return component;
}

[[nodiscard]] bool clip_depth(ReferenceArithmetic& math, const PgTraceDepth& a,
                              const PgTraceDepth& b, ExactRational& first, ExactRational& last) {
    const auto clip = [&](const ExactRational& at_start, const ExactRational& at_end) {
        if (at_start.sign() <= 0 && at_end.sign() <= 0)
            return false;
        if (at_start.sign() <= 0 || at_end.sign() <= 0) {
            auto crossing = math.div(at_start, math.sub(at_start, at_end));
            if (at_start.sign() <= 0) {
                if (math.compare(first, crossing) < 0)
                    first = std::move(crossing);
            } else if (math.compare(last, crossing) > 0) {
                last = std::move(crossing);
            }
        }
        return true;
    };
    return clip(math.sub(a.upper_first, b.lower_first), math.sub(a.upper_last, b.lower_last)) &&
           clip(math.sub(b.upper_first, a.lower_first), math.sub(b.upper_last, a.lower_last)) &&
           math.compare(first, last) < 0;
}

void fill_trace_columns(ReferenceArithmetic& math, const TracePatch& side,
                        const RestrictedTrace& trace, const PgTracePolynomial& opposite_normal,
                        std::span<ExactRational> entries) {
    const auto one = math.constant("1");
    const auto size = index_value(math, 4 * side.graph.n);
    const auto local = [&](const PgTracePolynomial& polynomial, core::Index cell) {
        PgTracePolynomial result;
        result.coefficients[0] = math.neg(index_value(math, 2 * cell + 1));
        trace_add_scaled(math, result, polynomial, size);
        return result;
    };
    const auto u = trace_powers(math, local(trace.s, trace.cell.first), pg_trace_basis_max_degree);
    const auto v = trace_powers(math, local(trace.t, trace.cell.second), pg_trace_basis_max_degree);
    const auto director =
        trace_product(math, trace_product(math, opposite_normal, opposite_normal), opposite_normal);
    const auto inv_length = math.div(one, side.frame.length_metres);
    const auto bending_scale = math.div(math.square(inv_length), math.constant("2"));
    const std::array<const ExactRational*, 3> scales{&one, &inv_length, &bending_scale};
    const auto rows = partition_take(side.rank.rows());
    const auto rank = static_cast<std::size_t>(side.columns);
    for (std::size_t column = 0; column < rank; ++column) {
        math.visit(1);
        const PgRankBasisColumn id{partition_index(column)};
        const auto stored = partition_take(side.rank.basis_indices(id));
        std::array<PgTracePolynomial, 3> fields;
        for (const auto ordinal : stored) {
            math.visit(1);
            const auto& row = rows[ordinal];
            if (row.cell != trace.cell)
                continue;
            math.visit(static_cast<std::uint64_t>(std::bit_width(stored.size())));
            math.charge();
            const auto coefficient =
                partition_take(side.rank.basis(id, ordinal, math.leaf_limits, math.plan));
            trace_add_scaled(math, fields[static_cast<std::size_t>(row.field)],
                             trace_product(math, u[static_cast<std::size_t>(row.powers.first)],
                                           v[static_cast<std::size_t>(row.powers.second)]),
                             coefficient);
        }
        fields[1] = trace_product(math, fields[1], director);
        std::size_t first = 0;
        for (std::size_t field = 0; field < fields.size(); ++field) {
            const auto count = pg_trace_field_slots[field];
            if (fields[field].degree >= count)
                partition_fail(core::FailureCode::ResourceExceeded);
            for (std::size_t j = 0; j < count; ++j)
                entries[column * pg_trace_row_count + first + j] =
                    math.mul(fields[field].coefficients[j], *scales[field]);
            first += count;
        }
    }
    // Apply B*T backwards to retain all earlier untransformed columns.
    for (std::size_t end = rank; end > 0; --end) {
        const auto j = end - 1;
        std::array<ExactRational, pg_trace_row_count> output{};
        for (std::size_t i = 0; i <= j; ++i) {
            math.visit(1);
            math.charge();
            const auto factor = partition_take(side.patch.conditioning().transform_coefficient(
                {partition_index(i)}, {partition_index(j)}, math.leaf_limits, math.plan));
            if (factor.sign() == 0)
                continue;
            for (std::size_t k = 0; k < output.size(); ++k) {
                math.visit(1);
                const auto& source = entries[i * pg_trace_row_count + k];
                if (source.sign() != 0)
                    output[k] = math.add(output[k], math.mul(source, factor));
            }
        }
        for (std::size_t k = 0; k < output.size(); ++k)
            entries[j * pg_trace_row_count + k] = std::move(output[k]);
    }
}

struct TraceStorage {
    std::span<PgTraceOriginalRow> rows;
    std::span<std::size_t> offsets;
    std::span<PgTraceSparseEntry> terms;
    std::span<PgTraceSparseEntry> scratch;
    std::size_t row_count = 0, unique_count = 0, entry_count = 0;
};

[[nodiscard]] PgCartesianExactPosition
trace_row_value(ReferenceArithmetic& math, std::span<const PgTraceSparseEntry> terms,
                std::span<const PgCartesianExactPosition> left,
                std::span<const PgCartesianExactPosition> right) {
    PgCartesianExactPosition value;
    for (const auto& term : terms) {
        const auto column = static_cast<std::size_t>(term.column);
        const auto& xyz = column < left.size() ? left[column] : right[column - left.size()];
        math.visit(1);
        for (std::size_t k = 0; k < value.size(); ++k)
            value[k] = math.add(value[k], math.mul(term.coefficient, xyz[k]));
    }
    return value;
}

void retain_trace_row(ReferenceArithmetic& math, TraceStorage& storage, PgTraceOriginalRow row,
                      std::span<ExactRational> table, std::size_t k,
                      std::span<const PgCartesianExactPosition> left,
                      std::span<const PgCartesianExactPosition> right) {
    std::size_t count = 0;
    for (std::size_t col = 0; col < left.size() + right.size(); ++col) {
        math.visit(1);
        auto& value = table[col * pg_trace_row_count + k];
        if (value.sign() == 0)
            continue;
        storage.scratch[count++] = {partition_index(col),
                                    col < left.size() ? std::move(value) : math.neg(value)};
    }
    auto terms = storage.scratch.first(count);
    if (!terms.empty()) {
        row.prescribed = trace_row_value(math, terms, left, right);
        if (!right.empty())
            for (const auto& value : row.prescribed)
                if (value.sign() != 0)
                    partition_fail(core::FailureCode::GeometryInvalid);
        row.scale = math.copy(terms.front().coefficient);
        for (auto& term : terms) {
            math.visit(1);
            term.coefficient = math.div(term.coefficient, row.scale);
        }
        std::size_t unique = 0;
        for (; unique < storage.unique_count; ++unique) {
            math.visit(1);
            const auto first = storage.offsets[unique];
            const auto prior = storage.terms.subspan(first, storage.offsets[unique + 1] - first);
            const auto same = [&](const auto& x, const auto& y) {
                math.visit(1);
                return x.column == y.column && math.compare(x.coefficient, y.coefficient) == 0;
            };
            if (std::ranges::equal(prior, terms, same))
                break;
        }
        if (unique == storage.unique_count) {
            if (terms.size() > storage.terms.size() - storage.entry_count)
                partition_fail(core::FailureCode::ResourceExceeded);
            for (auto& term : terms) {
                math.visit(1);
                storage.terms[storage.entry_count++] = std::move(term);
            }
            storage.offsets[++storage.unique_count] = storage.entry_count;
        }
        row.unique_row = partition_index(unique);
    }
    storage.rows[storage.row_count++] = std::move(row);
}

struct TraceCapacity {
    std::uint64_t overlaps = 0, columns = 0, rows = 0, terms = 0, table = 0;
};

[[nodiscard]] TraceCapacity admit_trace(ReferenceArithmetic& math, std::uint64_t columns,
                                        std::uint64_t max_overlaps, std::uint64_t maximum_entries) {
    const auto rows_count = partition_mul(max_overlaps, pg_trace_row_count);
    const auto term_cap = std::min(maximum_entries, partition_mul(rows_count, columns));
    const auto table_count = partition_mul(columns, pg_trace_row_count);
    // Bound new exact values in the dense patch row-block, powers, original
    // affine right-hand sides and unique sparse rows before arithmetic starts.
    auto live = partition_add(reference_live_rationals + 2048, table_count);
    live = partition_add(live, partition_mul(max_overlaps, 12));
    live = partition_add(
        live, partition_add(partition_mul(4, rows_count), partition_add(term_cap, columns)));
    math.admit(live);
    return {max_overlaps, columns, rows_count, term_cap, table_count};
}

struct TraceBuffers {
    core::PlannedBuffer<PgTraceOverlap> overlaps;
    core::PlannedBuffer<PgTraceOriginalRow> rows;
    core::PlannedBuffer<std::size_t> offsets;
    core::PlannedBuffer<PgTraceSparseEntry> terms, scratch;
    core::PlannedBuffer<ExactRational> table;
};

[[nodiscard]] TraceBuffers allocate_trace(ReferenceArithmetic& math, const TraceCapacity& shape) {
    auto& plan = math.plan;
    const auto [max_overlaps, columns, rows_count, term_cap, table_count] = shape;
    const auto offset_count = partition_add(rows_count, 1);
    auto bytes = partition_add(partition_mul(max_overlaps, sizeof(PgTraceOverlap)),
                               partition_mul(rows_count, sizeof(PgTraceOriginalRow)));
    bytes = partition_add(bytes, partition_mul(offset_count, sizeof(std::size_t)));
    bytes = partition_add(
        bytes, partition_mul(partition_add(columns, term_cap), sizeof(PgTraceSparseEntry)));
    bytes = partition_add(bytes, partition_mul(table_count, sizeof(ExactRational)));
    if (bytes > plan.remaining_bytes())
        partition_fail(core::FailureCode::ResourceExceeded);
    math.visit(partition_add(table_count, partition_add(term_cap, offset_count)));
    auto overlaps = partition_take(
        plan.allocate<PgTraceOverlap>(partition_index(max_overlaps), partition_operation));
    auto rows = partition_take(
        plan.allocate<PgTraceOriginalRow>(partition_index(rows_count), partition_operation));
    auto offsets = partition_take(
        plan.allocate<std::size_t>(partition_index(offset_count), partition_operation));
    auto terms = partition_take(
        plan.allocate<PgTraceSparseEntry>(partition_index(term_cap), partition_operation));
    auto scratch = partition_take(
        plan.allocate<PgTraceSparseEntry>(partition_index(columns), partition_operation));
    auto table = partition_take(
        plan.allocate<ExactRational>(partition_index(table_count), partition_operation));
    return {std::move(overlaps), std::move(rows),    std::move(offsets),
            std::move(terms),    std::move(scratch), std::move(table)};
}

void retain_overlap_rows(ReferenceArithmetic& math, TraceStorage& store,
                         std::span<ExactRational> table, std::size_t overlap,
                         std::span<const PgCartesianExactPosition> left,
                         std::span<const PgCartesianExactPosition> right = {}) {
    constexpr std::array fields{PgCoefficientField::R, PgCoefficientField::P,
                                PgCoefficientField::C};
    std::size_t k = 0;
    for (std::size_t field = 0; field < fields.size(); ++field)
        for (std::size_t power = 0; power < pg_trace_field_slots[field]; ++power)
            retain_trace_row(math, store,
                             {partition_index(overlap),
                              fields[field],
                              static_cast<std::uint32_t>(power),
                              -1,
                              {},
                              {}},
                             table, k++, left, right);
}

} // namespace

PgFullTrace::PgFullTrace(core::PlannedBuffer<PgTraceOverlap> overlaps,
                         core::PlannedBuffer<PgTraceOriginalRow> rows,
                         core::PlannedBuffer<std::size_t> offsets,
                         core::PlannedBuffer<PgTraceSparseEntry> terms, PgTraceInfo info) noexcept
    : overlaps_(std::move(overlaps)),
      rows_(std::move(rows)),
      offsets_(std::move(offsets)),
      coefficients_(std::move(terms)),
      info_(info) {}

core::Result<PgFullTrace>
PgFullTrace::prepare(const PgPreparedPatch& left, std::span<const PgTraceEdgeSpan> left_spans,
                     const PgPreparedPatch& right, std::span<const PgTraceEdgeSpan> right_spans,
                     const PgTraceLimits& limits, core::Workspace& workspace,
                     core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgFullTrace> {
        auto math = ReferenceArithmetic{limits.arithmetic, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto left_info = partition_take(left.info()),
                   right_info = partition_take(right.info());
        const auto a = trace_patch(left), b = trace_patch(right);
        const auto columns = partition_add(a.columns, b.columns);
        const auto max_overlaps = partition_add(left_spans.size(), right_spans.size());
        const auto shape = admit_trace(math, columns, max_overlaps, limits.maximum_nonzeros);
        validate_spans(math, left_spans);
        validate_spans(math, right_spans);
        if (math.compare(left_spans.front().first, right_spans.front().first) != 0 ||
            math.compare(left_spans.back().last, right_spans.back().last) != 0)
            partition_fail(core::FailureCode::InvalidInput);

        auto buffers = allocate_trace(math, shape);
        TraceStorage store{buffers.rows.values(), buffers.offsets.values(), buffers.terms.values(),
                           buffers.scratch.values()};
        std::size_t l = 0, r = 0, n = 0;
        while (l < left_spans.size() && r < right_spans.size()) {
            math.visit(1);
            const auto& x = left_spans[l];
            const auto& y = right_spans[r];
            const auto& first = math.compare(x.first, y.first) < 0 ? y.first : x.first;
            const auto ending = math.compare(x.last, y.last);
            const auto& last = ending < 0 ? x.last : y.last;
            if (math.compare(first, last) < 0) {
                auto xa = restrict_trace(math, a, x, first, last);
                auto ya = restrict_trace(math, b, y, first, last);
                const auto component = compare_geometry(math, xa, ya);
                ExactRational open_first;
                auto open_last = math.constant("1");
                if (clip_depth(math, xa.depth, ya.depth, open_first, open_last)) {
                    buffers.overlaps.values()[n] = {x.edge,
                                                    y.edge,
                                                    x.layer,
                                                    y.layer,
                                                    math.copy(first),
                                                    math.copy(last),
                                                    std::move(open_first),
                                                    std::move(open_last),
                                                    std::move(xa.depth),
                                                    std::move(ya.depth)};
                    const auto split = static_cast<std::size_t>(a.columns) * pg_trace_row_count;
                    fill_trace_columns(math, a, xa, ya.normal[component],
                                       buffers.table.values().first(split));
                    fill_trace_columns(math, b, ya, xa.normal[component],
                                       buffers.table.values().subspan(split));
                    retain_overlap_rows(math, store, buffers.table.values(), n,
                                        a.patch.reference_coordinates(),
                                        b.patch.reference_coordinates());
                    ++n;
                }
            }
            if (ending <= 0)
                ++l;
            if (ending >= 0)
                ++r;
        }
        if (n == 0)
            partition_fail(core::FailureCode::GeometryInvalid);
        math.cancel();
        return PgFullTrace{std::move(buffers.overlaps),
                           std::move(buffers.rows),
                           std::move(buffers.offsets),
                           std::move(buffers.terms),
                           {PgTraceKind::Connection, left_info.identity, right_info.identity,
                            partition_index(a.columns), partition_index(b.columns), n,
                            store.row_count, store.unique_count, store.entry_count, math.info()}};
    });
}

core::Result<PgFullTrace> PgFullTrace::prepare_support(
    const PgPreparedPatch& patch, std::span<const PgTraceEdgeSpan> spans,
    PgPreparedPatchIdentity identity, std::span<const PgCartesianExactPosition> prescribed,
    const PgTraceLimits& limits, core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgFullTrace> {
        auto math = ReferenceArithmetic{limits.arithmetic, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto reference = partition_take(patch.info());
        if (reference.identity != identity)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto a = trace_patch(patch);
        if (prescribed.size() != a.columns)
            partition_fail(core::FailureCode::InvalidInput);
        const auto shape = admit_trace(math, a.columns, spans.size(), limits.maximum_nonzeros);
        validate_spans(math, spans);
        auto buffers = allocate_trace(math, shape);
        TraceStorage store{buffers.rows.values(), buffers.offsets.values(), buffers.terms.values(),
                           buffers.scratch.values()};
        PgTracePolynomial cleared;
        cleared.coefficients[0] = math.constant("1");
        for (std::size_t i = 0; i < spans.size(); ++i) {
            math.visit(1);
            const auto& span = spans[i];
            auto trace = restrict_trace(math, a, span, span.first, span.last);
            buffers.overlaps.values()[i] = {span.edge,
                                            -1,
                                            span.layer,
                                            {},
                                            math.copy(span.first),
                                            math.copy(span.last),
                                            {},
                                            math.constant("1"),
                                            std::move(trace.depth),
                                            {}};
            fill_trace_columns(math, a, trace, cleared, buffers.table.values());
            retain_overlap_rows(math, store, buffers.table.values(), i, prescribed);
        }
        math.cancel();
        return PgFullTrace{std::move(buffers.overlaps),
                           std::move(buffers.rows),
                           std::move(buffers.offsets),
                           std::move(buffers.terms),
                           {PgTraceKind::Support,
                            identity,
                            {},
                            partition_index(a.columns),
                            0,
                            spans.size(),
                            store.row_count,
                            store.unique_count,
                            store.entry_count,
                            math.info()}};
    });
}

core::Result<PgTraceInfo> PgFullTrace::info() const noexcept {
    if (offsets_.size() == 0)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, partition_operation});
    return info_;
}

std::span<const PgTraceOverlap> PgFullTrace::overlaps() const& noexcept {
    return offsets_.size() == 0 ? std::span<const PgTraceOverlap>{}
                                : overlaps_.values().first(info_.overlaps);
}

std::span<const PgTraceOriginalRow> PgFullTrace::original_rows() const& noexcept {
    return offsets_.size() == 0 ? std::span<const PgTraceOriginalRow>{}
                                : rows_.values().first(info_.original_rows);
}

core::Result<std::span<const PgTraceSparseEntry>>
PgFullTrace::unique_row(core::Index row) const& noexcept {
    return partition_entry([&]() -> core::Result<std::span<const PgTraceSparseEntry>> {
        const auto current = partition_take(info());
        if (row < 0 || static_cast<std::uint64_t>(row) >= current.unique_rows)
            partition_fail(core::FailureCode::InvalidInput);
        const auto i = static_cast<std::size_t>(row);
        const auto offsets = offsets_.values();
        return coefficients_.values().subspan(offsets[i], offsets[i + 1] - offsets[i]);
    });
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgFullTrace::apply_impl(PgTraceCoordinates input, bool subtract_prescribed,
                        const PgReferenceLimits& limits, core::Workspace& workspace,
                        core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<PgCartesianExactPosition>> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto current = partition_take(info());
        if (input.left != current.left || input.right != current.right)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto left = static_cast<std::size_t>(current.left_columns);
        if (input.left_values.size() != left ||
            input.right_values.size() != static_cast<std::size_t>(current.right_columns))
            partition_fail(core::FailureCode::InvalidInput);
        math.admit(
            partition_add(reference_live_rationals, partition_mul(3, current.original_rows)));
        math.visit(current.original_rows);
        auto result = partition_take(plan.allocate<PgCartesianExactPosition>(
            partition_index(current.original_rows), partition_operation));
        for (std::size_t i = 0; i < current.original_rows; ++i) {
            math.visit(1);
            const auto& row = original_rows()[i];
            if (row.unique_row < 0)
                continue;
            auto value = trace_row_value(math, partition_take(unique_row(row.unique_row)),
                                         input.left_values, input.right_values);
            for (std::size_t k = 0; k < value.size(); ++k) {
                value[k] = math.mul(value[k], row.scale);
                if (subtract_prescribed)
                    value[k] = math.sub(value[k], row.prescribed[k]);
            }
            result.values()[i] = std::move(value);
        }
        math.cancel();
        return result;
    });
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgFullTrace::apply(PgTraceCoordinates input, const PgReferenceLimits& limits,
                   core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return apply_impl(input, false, limits, workspace, plan);
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgFullTrace::residual(PgTraceCoordinates input, const PgReferenceLimits& limits,
                      core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return apply_impl(input, true, limits, workspace, plan);
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgFullTrace::transpose(std::span<const PgCartesianExactPosition> lambda,
                       const PgReferenceLimits& limits, core::Workspace& workspace,
                       core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<PgCartesianExactPosition>> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto current = partition_take(info());
        if (lambda.size() != current.original_rows)
            partition_fail(core::FailureCode::InvalidInput);
        const auto size = partition_add(static_cast<std::uint64_t>(current.left_columns),
                                        static_cast<std::uint64_t>(current.right_columns));
        math.admit(partition_add(reference_live_rationals, partition_mul(3, size)));
        math.visit(size);
        auto result = partition_take(
            plan.allocate<PgCartesianExactPosition>(partition_index(size), partition_operation));
        for (std::size_t i = 0; i < lambda.size(); ++i) {
            math.visit(1);
            const auto& row = original_rows()[i];
            if (row.unique_row < 0)
                continue;
            PgCartesianExactPosition scaled;
            for (std::size_t k = 0; k < scaled.size(); ++k)
                scaled[k] = math.mul(row.scale, lambda[i][k]);
            for (const auto& term : partition_take(unique_row(row.unique_row))) {
                math.visit(1);
                auto& out = result.values()[static_cast<std::size_t>(term.column)];
                for (std::size_t k = 0; k < out.size(); ++k)
                    out[k] = math.add(out[k], math.mul(term.coefficient, scaled[k]));
            }
        }
        math.cancel();
        return result;
    });
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgFullTrace::restore_dual(std::span<const PgCartesianExactPosition> lambda,
                          const PgReferenceLimits& limits, core::Workspace& workspace,
                          core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<PgCartesianExactPosition>> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto current = partition_take(info());
        if (lambda.size() != current.unique_rows)
            partition_fail(core::FailureCode::InvalidInput);
        math.admit(
            partition_add(reference_live_rationals, partition_mul(3, current.original_rows)));
        math.visit(current.original_rows);
        auto result = partition_take(plan.allocate<PgCartesianExactPosition>(
            partition_index(current.original_rows), partition_operation));
        std::size_t next = 0;
        for (std::size_t i = 0; i < current.original_rows; ++i) {
            math.visit(1);
            const auto& row = original_rows()[i];
            if (row.unique_row < 0 || static_cast<std::size_t>(row.unique_row) != next)
                continue;
            auto& output = result.values()[i];
            for (std::size_t k = 0; k < output.size(); ++k)
                output[k] = math.div(lambda[next][k], row.scale);
            ++next;
        }
        if (next != lambda.size())
            partition_fail(core::FailureCode::RankUnresolved);
        math.cancel();
        return result;
    });
}

} // namespace uvlas::structure
