#include "pg_reference_projection.hpp"
#include "pg_volume_position_arithmetic.hpp"

#include <algorithm>
#include <uvlas/structure/pg_chart_candidate.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;
using namespace numerics;
namespace linear = numerics::linear;
using Words = numerics::TwoWordValue;
using Buffer = numerics::PrecisionBuffer<double>;
using Solved = std::variant<Buffer, linear::Solution<double>>;

struct Entry {
    std::size_t column;
    double high, error;
};

struct VisitBudget {
    std::uint64_t remaining;
    std::stop_token cancellation;

    void visit() {
        if (cancellation.stop_requested())
            partition_fail(core::FailureCode::Cancelled);
        if (remaining == 0)
            partition_fail(core::FailureCode::ResourceExceeded);
        --remaining;
    }
};

template <class T>
auto allocate(core::ResourcePlan& plan, std::uint64_t n) {
    return partition_take(plan.allocate<T>(partition_index(n), volume_operation));
}

struct Planes {
    std::span<const double> high, low;
    std::size_t stride;

    [[nodiscard]] Words at(std::size_t i, std::size_t c) const noexcept {
        const auto k = c * stride + i;
        return {high[k], low.empty() ? 0.0 : low[k]};
    }
    [[nodiscard]] bool two_word() const noexcept {
        return !low.empty();
    }
};

Planes planes(PgSparseChartIdentity id, PgChartCandidateInput input, std::size_t rows) {
    if (input.identity != id)
        partition_fail(core::FailureCode::StaleGeneration);
    const auto shape = input.value.shape();
    if (shape.rows != partition_index(rows) || shape.columns != 3)
        partition_fail(core::FailureCode::InvalidInput);
    const auto stride = static_cast<std::size_t>(shape.leading_dimension);
    if (const auto* binary = std::get_if<Binary64<double>>(&input.value.storage()))
        return {binary->values.values(), {}, stride};
    const auto& pair = std::get<TwoWord<double>>(input.value.storage());
    return {pair.hi.values(), pair.lo.values(), stride};
}

const Buffer& solved_values(const Solved& result) noexcept {
    if (const auto* value = std::get_if<Buffer>(&result))
        return *value;
    return std::get<linear::Solution<double>>(result).values();
}
VolumeInteger lifted(Words a) {
    return volume_lift(volume_decode(a.hi)) + volume_lift(volume_decode(a.lo));
}
VolumeInteger product(double a, Words b) {
    const auto x = volume_decode(a);
    return volume_product(x, volume_decode(b.hi)) + volume_product(x, volume_decode(b.lo));
}
VolumeInteger dot(Words a, Words b) {
    return product(a.hi, b) + product(a.lo, b);
}

VolumeInteger allowance(double limit) {
    if (!volume_finite(limit))
        partition_fail(core::FailureCode::Nonfinite);
    auto bound = volume_lift(volume_decode(limit));
    if (bound < 0)
        partition_fail(core::FailureCode::InvalidInput);
    return bound;
}

double upper(const VolumeInteger& exact) {
    const auto word = partition_take(volume_round(exact));
    if (volume_lift(volume_decode(word)) >= exact)
        return word;
    const auto bumped = std::bit_cast<double>(std::bit_cast<std::uint64_t>(word) + 1);
    if (!volume_finite(bumped))
        partition_fail(core::FailureCode::PrecisionUnresolved);
    return bumped;
}

Words round_pair(const VolumeInteger& value) {
    const auto hi = partition_take(volume_round(value));
    const auto lo = partition_take(volume_round(value - volume_lift(volume_decode(hi))));
    return {hi, lo};
}

Buffer buffer(std::size_t rows, std::span<const Words> values, bool two_word,
              core::NumericalGeneration generation, core::ResourcePlan& plan) {
    const DenseShape shape{partition_index(rows), 3, partition_index(rows)};
    auto high = allocate<double>(plan, values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
        high.values()[i] = values[i].hi;
    auto hi = partition_take(DenseReal::copy(shape, high.values(), plan, volume_operation));
    if (!two_word)
        return partition_take(Buffer::binary64(std::move(hi), generation));
    auto low = allocate<double>(plan, values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
        low.values()[i] = values[i].lo;
    auto lo = partition_take(DenseReal::copy(shape, low.values(), plan, volume_operation));
    return partition_take(Buffer::two_word(std::move(hi), std::move(lo), generation));
}
} // namespace

struct PgChartCandidate::Data {
    PgSparseChartIdentity identity;
    std::optional<linear::Context> forward, backward;
    std::optional<linear::SparsePrepared> factor;
    SparseReal matrix;
    core::PlannedBuffer<std::size_t> offsets;
    core::PlannedBuffer<Entry> entries;
    core::PlannedBuffer<core::Index> dependent, free;
    std::size_t rows{}, columns{}, rank{}, independent{};

    [[nodiscard]] std::span<const Entry> row(std::size_t r) const noexcept {
        return entries.values().subspan(offsets.values()[r],
                                        offsets.values()[r + 1] - offsets.values()[r]);
    }

    [[nodiscard]] Solved solve(std::span<const Words> rhs, bool transpose,
                               const PgChartCandidateSolve& request, core::ResourcePlan& plan) {
        const auto expected = transpose ? backward : forward;
        if (!expected || request.scales.context != *expected)
            partition_fail(core::FailureCode::StaleGeneration);
        auto owned = buffer(rank, rhs, true, request.rhs_buffer_generation, plan);
        if (rank == 0)
            return Solved{std::move(owned)};
        const linear::SparseOperator a{matrix, *forward};
        const linear::Rhs<double> b{owned, *expected, request.rhs_generation};
        auto result = partition_take(
            transpose
                ? linear::solve_transposed_linear(a, b, request.scales, request.precision,
                                                  request.accuracy, request.solver, *factor, plan)
                : linear::solve_linear(a, b, request.scales, request.precision, request.accuracy,
                                       request.solver, *factor, plan));
        return Solved{std::move(result)};
    }
};

PgChartCandidate::PgChartCandidate(PgSparseChart&& exact, core::PlannedBuffer<Data> data) noexcept
    : exact_(std::move(exact)), data_(std::move(data)) {}
PgChartCandidate::PgChartCandidate(PgChartCandidate&&) noexcept = default;
PgChartCandidate& PgChartCandidate::operator=(PgChartCandidate&&) noexcept = default;
PgChartCandidate::~PgChartCandidate() = default;

core::Result<PgChartCandidate> PgChartCandidate::prepare(PgSparseChart&& chart,
                                                         const linear::Scales& forward,
                                                         const linear::Scales& transpose,
                                                         const PgChartCandidatePreparation& request,
                                                         core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgChartCandidate> {
        ReferenceArithmetic math{request.projection, plan};
        math.cancel();
        if (forward.context.cache != transpose.context.cache ||
            forward.context.trial != transpose.context.trial ||
            forward.context.time != transpose.context.time || !volume_finite(forward.context.time))
            partition_fail(core::FailureCode::StaleGeneration);
        static_cast<void>(allowance(request.maximum_absolute_coefficient_error));
        const auto identity = partition_take(chart.identity());
        const auto& original = chart.original();
        const auto info = partition_take(original.info());
        auto memory = allocate<Data>(plan, 1);
        auto& data = memory.values().front();
        data.identity = identity;
        data.rows = static_cast<std::size_t>(info.rows);
        data.columns = static_cast<std::size_t>(info.columns);
        data.rank = chart.pivots().size();
        data.independent = chart.free_columns().size();
        data.forward = forward.context;
        data.backward = transpose.context;
        if (request.solver.maximum_rhs_columns < 3 ||
            request.solver.maximum_order < partition_index(data.rank))
            partition_fail(core::FailureCode::ResourceExceeded);
        data.entries = allocate<Entry>(plan, static_cast<std::uint64_t>(info.nonzeros));
        data.offsets = allocate<std::size_t>(plan, partition_add(data.rows, 1));
        data.dependent = allocate<core::Index>(plan, data.columns);
        data.free = allocate<core::Index>(plan, data.columns);
        for (std::size_t i = 0; i < data.columns; ++i) {
            math.visit(1);
            data.dependent.values()[i] = -1;
            data.free.values()[i] = -1;
        }
        for (std::size_t i = 0; i < data.rank; ++i) {
            math.visit(1);
            const auto column = chart.pivots()[i].dependent_column;
            data.dependent.values()[static_cast<std::size_t>(column)] = partition_index(i);
        }
        for (std::size_t i = 0; i < data.independent; ++i) {
            math.visit(1);
            data.free.values()[static_cast<std::size_t>(chart.free_columns()[i])] =
                partition_index(i);
        }
        math.admit();
        const auto tolerance = reference_binary64(math, request.maximum_absolute_coefficient_error);
        const ExactRational zero;
        const PgProjectionComparator comparator{math, zero, 0};
        std::size_t used = 0;
        for (std::size_t i = 0; i < data.rows; ++i) {
            math.visit(1);
            for (const auto& term : partition_take(original.row(partition_index(i)))) {
                math.visit(1);
                const auto high = pg_round_word(comparator, term.coefficient);
                if (high == 0.0)
                    partition_fail(core::FailureCode::PrecisionUnresolved);
                const auto difference =
                    math.abs(math.sub(term.coefficient, reference_binary64(math, high)));
                if (math.compare(difference, tolerance) > 0)
                    partition_fail(core::FailureCode::PrecisionUnresolved);
                auto error = pg_round_word(comparator, difference);
                if (math.compare(reference_binary64(math, error), difference) < 0)
                    error = std::bit_cast<double>(std::bit_cast<std::uint64_t>(error) + 1);
                data.entries.values()[used++] = {static_cast<std::size_t>(term.column), high,
                                                 error};
            }
            data.offsets.values()[i + 1] = used;
        }

        std::uint64_t selected = 0;
        for (const auto& pivot : chart.pivots())
            for (const auto& term : data.row(static_cast<std::size_t>(pivot.original_row))) {
                math.visit(1);
                selected += data.dependent.values()[term.column] >= 0 ? 1U : 0U;
            }
        auto entries = allocate<Triplet>(plan, selected);
        used = 0;
        for (std::size_t i = 0; i < data.rank; ++i)
            for (const auto& term :
                 data.row(static_cast<std::size_t>(chart.pivots()[i].original_row))) {
                math.visit(1);
                const auto column = data.dependent.values()[term.column];
                if (column >= 0)
                    entries.values()[used++] = {partition_index(i), column, term.high};
            }
        if (data.rank != 0) {
            const auto n = partition_index(data.rank);
            data.matrix = partition_take(
                SparseReal::from_triplets(n, n, entries.values(), plan, volume_operation));
            const linear::SparseOperator op{data.matrix, forward.context};
            math.cancel();
            data.factor = partition_take(linear::SparsePrepared::prepare(
                op, forward, transpose, linear::SparseMethod::GeneralLu, request.solver, plan));
        } else if (!forward.unknown.empty() || !forward.residual.empty() ||
                   !forward.unknown_units.empty() || !forward.residual_units.empty() ||
                   !transpose.unknown.empty() || !transpose.residual.empty() ||
                   !transpose.unknown_units.empty() || !transpose.residual_units.empty()) {
            partition_fail(core::FailureCode::InvalidInput);
        }
        math.cancel();
        return PgChartCandidate{std::move(chart), std::move(memory)};
    });
}

const PgSparseChart& PgChartCandidate::exact() const& noexcept {
    return exact_;
}

core::Result<PgChartCandidateLift> PgChartCandidate::lift(PgChartCandidateInput independent,
                                                          PgChartCandidateInput original_rhs,
                                                          const PgChartCandidateSolve& request,
                                                          core::Workspace& workspace,
                                                          core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgChartCandidateLift> {
        VisitBudget budget{request.maximum_visits, request.cancellation};
        budget.visit();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        if (data_.size() == 0)
            partition_fail(core::FailureCode::StaleGeneration);
        auto& data = data_.values().front();
        const auto r = planes(data.identity, independent, data.independent);
        const auto c = planes(data.identity, original_rhs, data.rows);
        const auto allowed = allowance(request.maximum_absolute_original_defect);
        auto rhs = allocate<Words>(plan, partition_mul(3, data.rank));
        for (std::size_t axis = 0; axis < 3; ++axis)
            for (std::size_t i = 0; i < data.rank; ++i) {
                budget.visit();
                const auto row = static_cast<std::size_t>(exact_.pivots()[i].original_row);
                auto value = lifted(c.at(row, axis));
                for (const auto& term : data.row(row)) {
                    budget.visit();
                    const auto free = data.free.values()[term.column];
                    if (free >= 0)
                        value -= product(term.high, r.at(static_cast<std::size_t>(free), axis));
                }
                rhs.values()[axis * data.rank + i] = round_pair(value);
            }
        budget.visit();
        const auto solved = data.solve(rhs.values(), false, request, plan);
        const auto dependent =
            planes(data.identity, {data.identity, solved_values(solved)}, data.rank);
        auto output = allocate<Words>(plan, partition_mul(3, data.columns));
        VolumeInteger worst = 0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            auto lane = output.values().subspan(axis * data.columns, data.columns);
            for (std::size_t col = 0; col < lane.size(); ++col) {
                budget.visit();
                const auto p = data.dependent.values()[col];
                const auto f = data.free.values()[col];
                lane[col] = p >= 0 ? dependent.at(static_cast<std::size_t>(p), axis)
                                   : r.at(static_cast<std::size_t>(f), axis);
            }
            for (std::size_t row = 0; row < data.rows; ++row) {
                budget.visit();
                auto residual = lifted(c.at(row, axis));
                VolumeInteger error = 0;
                for (const auto& term : data.row(row)) {
                    budget.visit();
                    const auto x = lane[term.column];
                    residual -= product(term.high, x);
                    error += rank_abs(product(term.error, x));
                }
                worst = std::max(worst, rank_abs(residual) + error);
            }
        }
        if (worst > allowed)
            partition_fail(core::FailureCode::PrecisionUnresolved);
        auto result = buffer(data.columns, output.values(), r.two_word() || dependent.two_word(),
                             request.output_generation, plan);
        budget.visit();
        return PgChartCandidateLift{data.identity, std::move(result), upper(worst)};
    });
}

core::Result<PgChartCandidateDual> PgChartCandidate::transpose(PgChartCandidateInput force,
                                                               const PgChartCandidateSolve& request,
                                                               core::Workspace& workspace,
                                                               core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgChartCandidateDual> {
        VisitBudget budget{request.maximum_visits, request.cancellation};
        budget.visit();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        if (data_.size() == 0)
            partition_fail(core::FailureCode::StaleGeneration);
        auto& data = data_.values().front();
        const auto f = planes(data.identity, force, data.columns);
        const auto allowed = allowance(request.maximum_absolute_original_defect);
        auto rhs = allocate<Words>(plan, partition_mul(3, data.rank));
        for (std::size_t axis = 0; axis < 3; ++axis)
            for (std::size_t i = 0; i < data.rank; ++i) {
                budget.visit();
                rhs.values()[axis * data.rank + i] =
                    f.at(static_cast<std::size_t>(exact_.pivots()[i].dependent_column), axis);
            }
        budget.visit();
        const auto solved = data.solve(rhs.values(), true, request, plan);
        const auto multipliers =
            planes(data.identity, {data.identity, solved_values(solved)}, data.rank);
        auto balance = allocate<VolumeInteger>(plan, data.columns);
        auto errors = allocate<VolumeInteger>(plan, data.columns);
        auto lambda = allocate<Words>(plan, partition_mul(3, data.rows));
        auto h = allocate<Words>(plan, partition_mul(3, data.independent));
        VolumeInteger worst = 0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            auto reaction = lambda.values().subspan(axis * data.rows, data.rows);
            auto reduced = h.values().subspan(axis * data.independent, data.independent);
            for (std::size_t col = 0; col < data.columns; ++col) {
                budget.visit();
                balance.values()[col] = lifted(f.at(col, axis));
                errors.values()[col] = 0;
            }
            for (std::size_t i = 0; i < data.rank; ++i) {
                budget.visit();
                const auto row = static_cast<std::size_t>(exact_.pivots()[i].original_row);
                const auto x = multipliers.at(i, axis);
                reaction[row] = x;
                for (const auto& term : data.row(row)) {
                    budget.visit();
                    balance.values()[term.column] -= product(term.high, x);
                    errors.values()[term.column] += rank_abs(product(term.error, x));
                }
            }
            for (std::size_t col = 0; col < data.columns; ++col) {
                budget.visit();
                auto& residual = balance.values()[col];
                const auto free = data.free.values()[col];
                if (free >= 0) {
                    auto& x = reduced[static_cast<std::size_t>(free)];
                    x = round_pair(residual);
                    residual -= lifted(x);
                }
                worst = std::max(worst, rank_abs(residual) + errors.values()[col]);
            }
        }
        if (worst > allowed)
            partition_fail(core::FailureCode::PrecisionUnresolved);
        auto force_result =
            buffer(data.independent, h.values(), true, request.output_generation, plan);
        auto reaction_result = buffer(data.rows, lambda.values(), multipliers.two_word(),
                                      request.output_generation, plan);
        budget.visit();
        return PgChartCandidateDual{data.identity, std::move(force_result),
                                    std::move(reaction_result), upper(worst)};
    });
}

core::Result<double>
PgChartCandidate::work_defect(PgChartCandidateInput force, PgChartCandidateInput independent,
                              PgChartCandidateInput rhs, const PgChartCandidateLift& lift,
                              const PgChartCandidateDual& dual, double maximum_absolute_error,
                              std::uint64_t maximum_visits,
                              std::stop_token cancellation) const noexcept {
    return partition_entry([&]() -> core::Result<double> {
        VisitBudget budget{maximum_visits, cancellation};
        budget.visit();
        if (data_.size() == 0)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto& data = data_.values().front();
        const auto allowed = allowance(maximum_absolute_error);
        const auto f = planes(data.identity, force, data.columns);
        const auto q = planes(data.identity, {lift.identity, lift.coordinates}, data.columns);
        const auto r = planes(data.identity, independent, data.independent);
        const auto h =
            planes(data.identity, {dual.identity, dual.independent_force}, data.independent);
        const auto c = planes(data.identity, rhs, data.rows);
        const auto lambda =
            planes(data.identity, {dual.identity, dual.original_multipliers}, data.rows);
        VolumeInteger difference = 0;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            for (std::size_t i = 0; i < data.columns; ++i) {
                budget.visit();
                difference += dot(f.at(i, axis), q.at(i, axis));
            }
            for (std::size_t i = 0; i < data.independent; ++i) {
                budget.visit();
                difference -= dot(h.at(i, axis), r.at(i, axis));
            }
            for (std::size_t i = 0; i < data.rows; ++i) {
                budget.visit();
                difference -= dot(lambda.at(i, axis), c.at(i, axis));
            }
        }
        const auto magnitude = rank_abs(difference);
        if (magnitude > allowed)
            partition_fail(core::FailureCode::PrecisionUnresolved);
        return upper(magnitude);
    });
}

} // namespace uvlas::structure
