#include "pg_coefficients_plan.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <uvlas/core/instrumentation.hpp>

namespace uvlas::structure {
namespace {
constexpr auto operation = core::Operation::PrepareModel;
core::Failure invalid() noexcept {
    return {core::FailureCode::InvalidInput, operation};
}
const detail::PgCoefficientData* owned(const PgCoefficientDictionary& value) noexcept {
    return detail::PgCoefficientFactory::data(value);
}
bool same_graph(const PgGraphView& a, const PgGraphView& b) noexcept {
    if (a.identity != b.identity || a.n != b.n || a.convention != b.convention || !(a.x0 == b.x0) ||
        !(a.width == b.width) || !(a.skew == b.skew) || a.cells.size() != b.cells.size())
        return false;
    for (std::size_t i = 0; i < a.cells.size(); ++i)
        if (a.cells[i] != b.cells[i])
            return false;
    const auto f = a.height.terms(), g = b.height.terms();
    if (f.size() != g.size())
        return false;
    for (std::size_t i = 0; i < f.size(); ++i)
        if (f[i].powers.first != g[i].powers.first || f[i].powers.second != g[i].powers.second ||
            !(f[i].coefficient == g[i].coefficient))
            return false;
    return true;
}
core::Result<std::size_t> cell_ordinal(const detail::PgCoefficientData& data,
                                       PgCellIndex cell) noexcept {
    const auto axis = 2 * data.n;
    if (cell.first < 0 || cell.second < 0 || cell.first >= axis || cell.second >= axis)
        return std::unexpected(invalid());
    return static_cast<std::size_t>(cell.first * axis + cell.second);
}
} // namespace
namespace detail {
const PgCoefficientData* PgCoefficientFactory::data(const PgCoefficientDictionary& value) noexcept {
    return value.data_.size() == 1 ? value.data_.values().data() : nullptr;
}
PgCoefficientDictionary PgCoefficientFactory::make(core::PlannedBuffer<PgCoefficientData> data,
                                                   PgCoefficientPreparation info) noexcept {
    return PgCoefficientDictionary{std::move(data), info};
}
} // namespace detail
PgCoefficientDictionary::PgCoefficientDictionary(
    core::PlannedBuffer<detail::PgCoefficientData> data, PgCoefficientPreparation info) noexcept
    : data_(std::move(data)), info_(info) {}
PgCoefficientDictionary::PgCoefficientDictionary(PgCoefficientDictionary&&) noexcept = default;
PgCoefficientDictionary&
PgCoefficientDictionary::operator=(PgCoefficientDictionary&&) noexcept = default;
PgCoefficientDictionary::~PgCoefficientDictionary() = default;
core::Result<PgGraphView> PgCoefficientDictionary::graph() const& noexcept {
    const auto* data = owned(*this);
    if (data == nullptr || !data->source)
        return std::unexpected(invalid());
    return data->source->view();
}
core::Result<void>
PgCoefficientDictionary::validate_for(const PgGraph& graph_value) const noexcept {
    const auto actual = graph();
    if (!actual)
        return std::unexpected(actual.error());
    const auto requested = graph_value.view();
    if (!requested)
        return std::unexpected(requested.error());
    if (!same_graph(*actual, *requested))
        return std::unexpected(core::Failure{core::FailureCode::StaleGeneration, operation});
    return {};
}
core::Result<std::span<const PgCoefficientSource>>
PgCoefficientDictionary::columns() const& noexcept {
    const auto* data = owned(*this);
    if (data == nullptr)
        return std::unexpected(invalid());
    return data->columns.values().first(data->column_count);
}
core::Result<core::Index> PgCoefficientDictionary::original_columns() const noexcept {
    const auto* data = owned(*this);
    if (data == nullptr)
        return std::unexpected(invalid());
    return data->original;
}
core::Result<std::span<const PgCoefficientRecord>>
PgCoefficientDictionary::records() const& noexcept {
    const auto* data = owned(*this);
    if (data == nullptr)
        return std::unexpected(invalid());
    return data->records.values().first(data->record_count);
}
core::Result<const ExactCellPolynomial*>
PgCoefficientDictionary::coefficient(PgCoefficientColumn column, PgCellIndex cell,
                                     PgCoefficientField field) const& noexcept {
    const auto* data = owned(*this);
    if (data == nullptr || column.value < 0 ||
        static_cast<std::size_t>(column.value) >= data->column_count ||
        (field != PgCoefficientField::R && field != PgCoefficientField::P &&
         field != PgCoefficientField::C))
        return std::unexpected(invalid());
    const auto index = cell_ordinal(*data, cell);
    if (!index)
        return std::unexpected(index.error());
    const auto key = static_cast<std::uint64_t>(column.value) *
                         static_cast<std::uint64_t>(4 * data->n * data->n) * 3 +
                     static_cast<std::uint64_t>(*index) * 3 + static_cast<std::uint64_t>(field);
    const auto all = data->records.values().first(data->record_count);
    const auto found = std::lower_bound(
        all.begin(), all.end(), key, [&](const PgCoefficientRecord& entry, std::uint64_t target) {
            const auto cell_key = entry.cell.first * (2 * data->n) + entry.cell.second;
            const auto entry_key = static_cast<std::uint64_t>(entry.column.value) *
                                       static_cast<std::uint64_t>(4 * data->n * data->n) * 3 +
                                   static_cast<std::uint64_t>(cell_key) * 3 +
                                   static_cast<std::uint64_t>(entry.field);
            return entry_key < target;
        });
    if (found != all.end() && found->column == column && found->cell == cell &&
        found->field == field)
        return &found->polynomial;
    return &data->zero;
}
core::Result<const ExactCellPolynomial*>
PgCoefficientDictionary::g(PgCellIndex cell) const& noexcept {
    const auto* data = owned(*this);
    if (data == nullptr)
        return std::unexpected(invalid());
    const auto index = cell_ordinal(*data, cell);
    if (!index)
        return std::unexpected(index.error());
    return &data->denominator.values()[*index];
}
core::Result<std::span<const PgRequiredMapEntry>>
PgCoefficientDictionary::required_map() const& noexcept {
    const auto* data = owned(*this);
    if (data == nullptr)
        return std::unexpected(invalid());
    return data->required.values().first(data->required_count);
}
core::Result<std::span<const PgAnchorCMapEntry>>
PgCoefficientDictionary::anchor_c_map() const& noexcept {
    const auto* data = owned(*this);
    if (data == nullptr)
        return std::unexpected(invalid());
    return data->c_map.values().first(data->c_map_count);
}
const PgCoefficientPreparation& PgCoefficientDictionary::preparation() const& noexcept {
    return info_;
}
core::Result<PgCoefficientDictionary>
PgCoefficientDictionary::prepare(const PgGraph& graph_value, const PgCoefficientLimits& limits,
                                 core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::InvalidInput, {},
        [&]() -> core::Result<PgCoefficientDictionary> {
            try {
                using namespace detail;
                const auto view = pg_take(graph_value.view());
                auto acquired = workspace.acquire();
                if (!acquired) {
                    auto failure = acquired.error();
                    failure.operation = operation;
                    return std::unexpected(failure);
                }
                const auto lease = std::move(*acquired);
                const auto before = plan.allocated_bytes();
                auto forecast = pg_coefficient_plan(view, limits, plan);
                PgCoefficientAlgebra op{limits, plan};
                const bool measured = limits.algebra.timing == ExactTimingMode::Enabled;
                core::StageRecorder coefficients{measured, operation}, anchors{measured, operation},
                    verification{measured, operation};
                auto storage = pg_take(plan.allocate<PgCoefficientData>(1, operation));
                auto& data = storage.values()[0];
                data.n = view.n;
                data.source = pg_take(
                    graph_value.clone({op.limits, forecast.cells, limits.cancellation}, plan));
                data.columns = pg_take(plan.allocate<PgCoefficientSource>(
                    static_cast<core::Index>(forecast.full_columns), operation));
                data.records = pg_take(plan.allocate<PgCoefficientRecord>(
                    static_cast<core::Index>(forecast.record_capacity), operation));
                data.denominator = pg_take(
                    plan.allocate<PgPoly>(static_cast<core::Index>(forecast.cells), operation));
                data.required = pg_take(plan.allocate<PgRequiredMapEntry>(
                    static_cast<core::Index>(forecast.map_capacity), operation));
                data.c_map = pg_take(plan.allocate<PgAnchorCMapEntry>(
                    static_cast<core::Index>(forecast.c_map_capacity), operation));
                PgCoefficientGeneration generated{data, op, limits};
                const auto inverse = op.rd(op.number(1), view.width);
                const auto scale = op.number(4 * view.n);
                {
                    const core::StageScope phase{coefficients, core::MeasurementStage::Preparation};
                    generated.original_sources();
                    generated.geometry(view, inverse, scale);
                    generated.basic();
                    generated.paired(view, inverse, scale);
                }
                PgCoefficientCompletion completion{generated};
                {
                    const core::StageScope phase{anchors, core::MeasurementStage::Preparation};
                    completion.anchor();
                    completion.c_map();
                    completion.verify_anchor();
                }
                {
                    const core::StageScope phase{verification, core::MeasurementStage::Preparation};
                    PgCoefficientMapBuilder maps{generated};
                    maps.required(view);
                    maps.verify_basic(view);
                    maps.verify_curvatures(view, inverse, scale);
                }
                op.check();
                auto& info = forecast.report;
                const auto inspect_rational = [&](const ExactRational& value) {
                    info.observed_maximum_value_bits =
                        std::max(info.observed_maximum_value_bits,
                                 std::max(value.numerator_bits(), value.denominator_bits()));
                };
                const auto inspect_polynomial = [&](const PgPoly& value) {
                    info.observed_nonzero_terms =
                        pg_add(info.observed_nonzero_terms, value.terms().size());
                    for (const auto& term : value.terms())
                        inspect_rational(term.coefficient);
                };
                for (const auto& record : data.records.values().first(data.record_count))
                    inspect_polynomial(record.polynomial);
                for (const auto& denominator : data.denominator.values())
                    inspect_polynomial(denominator);
                for (const auto& entry : data.required.values().first(data.required_count))
                    inspect_rational(entry.coefficient);
                for (const auto& entry : data.c_map.values().first(data.c_map_count))
                    inspect_rational(entry.coefficient);
                info.observed_columns = data.column_count;
                info.allocated_owned_payload_bytes = plan.allocated_bytes() - before;
                // These checks detect disagreement between the pure recipe ledger
                // and real leaves. They never replace pre-allocation admission.
                if (info.allocated_owned_payload_bytes > info.planned_owned_payload_upper_bytes ||
                    info.observed_nonzero_terms > limits.maximum_stored_terms ||
                    op.leaf_maximum.algebraic_work_bits > info.arithmetic.algebraic_work_bits ||
                    op.leaf_maximum.planned_work_bits > info.arithmetic.planned_work_bits ||
                    op.leaf_maximum.planned_raw_terms > info.arithmetic.planned_raw_terms ||
                    op.leaf_maximum.estimated_external_payload_bytes >
                        info.arithmetic.estimated_external_payload_bytes)
                    PgCoefficientGeneration::exceeded();
                const auto coefficient_times = pg_take(coefficients.snapshot());
                const auto anchor_times = pg_take(anchors.snapshot());
                const auto verification_times = pg_take(verification.snapshot());
                if (measured) {
                    constexpr auto stage =
                        static_cast<std::size_t>(core::MeasurementStage::Preparation);
                    info.timings = PgCoefficientTimings{
                        coefficient_times[stage], anchor_times[stage], verification_times[stage]};
                }
                return PgCoefficientFactory::make(std::move(storage), info);
            } catch (const detail::PgCoefficientAbort& failure) {
                return std::unexpected(failure.failure);
            }
        });
}
} // namespace uvlas::structure
