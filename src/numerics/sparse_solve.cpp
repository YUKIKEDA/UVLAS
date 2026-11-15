#include "backend_memory_detail.hpp"
#include "compensated.hpp"
#include "linear_solve_detail.hpp"
#include "sparse_spd.hpp"
#include "sparse_triangular.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <limits>
#include <utility>
#include <uvlas/numerics/sparse_solve.hpp>

namespace uvlas::numerics::linear {
namespace {
constexpr auto operation = core::Operation::SolveLinear;
using Accumulator = uvlas::numerics::detail::CompensatedAccumulator;

auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}
core::Result<void> check_operator(const SparseOperator& input, const Limits& limits) noexcept {
    const auto valid = detail::check_context_inputs(input.context());
    if (!valid) {
        return valid;
    }
    if (input.context().cache.representation != core::Precision::Binary64) {
        return fail(core::FailureCode::PrecisionUnsupported);
    }
    if (input.values().rows() != input.values().columns()) {
        return fail(core::FailureCode::InvalidInput);
    }
    return detail::check_capacity_inputs(input.values().rows(), 0, 0, limits);
}
core::Result<BackendMemoryPlan> backend_plan(const SparseReal& matrix) noexcept {
    const auto n = static_cast<std::uint64_t>(matrix.rows());
    const auto nnz = static_cast<std::uint64_t>(matrix.values().size());
    if (n == 0) {
        return detail::BackendPlanFactory::make(0, 0);
    }
    // Deliberately a forecast: stock Eigen's data-dependent expansion is not
    // allocator-capped. Include ordering/copies/expansion scratch, never label
    // this multiplier a proved fill or whole-process upper bound.
    constexpr std::uint64_t per_entry = 1024, per_row = 512;
    if (nnz > std::numeric_limits<std::uint64_t>::max() / per_entry ||
        n > std::numeric_limits<std::uint64_t>::max() / per_row) {
        return fail(core::FailureCode::ResourceExceeded);
    }
    const auto entries = nnz * per_entry;
    const auto rows = n * per_row;
    if (entries > std::numeric_limits<std::uint64_t>::max() - rows) {
        return fail(core::FailureCode::ResourceExceeded);
    }
    return detail::BackendPlanFactory::make(entries + rows);
}
bool same_matrix(const SparseReal& a, const SparseReal& b) noexcept {
    return a.rows() == b.rows() && a.columns() == b.columns() &&
           std::ranges::equal(a.column_offsets(), b.column_offsets()) &&
           std::ranges::equal(a.row_indices(), b.row_indices()) &&
           std::ranges::equal(a.values(), b.values(), [](double x, double y) {
               return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y);
           });
}
struct WitnessAssessment {
    SparseSpdStatus status = SparseSpdStatus::Uncertified;
    std::optional<double> quadratic;
};
core::Result<WitnessAssessment> check_witness(const SparseReal& matrix,
                                              std::span<const double> witness,
                                              std::span<Accumulator> rows) noexcept {
    if (!std::ranges::any_of(witness, [](double value) { return value != 0; })) {
        return WitnessAssessment{};
    }
    std::fill(rows.begin(), rows.end(), Accumulator{});
    for (core::Index column = 0; column < matrix.columns(); ++column) {
        for (auto entry = matrix.column_offsets()[static_cast<std::size_t>(column)];
             entry < matrix.column_offsets()[static_cast<std::size_t>(column + 1)]; ++entry) {
            const auto index = static_cast<std::size_t>(entry);
            auto added = rows[static_cast<std::size_t>(matrix.row_indices()[index])].product(
                matrix.values()[index], witness[static_cast<std::size_t>(column)], 1);
            if (!added) {
                if (added.error().code == core::FailureCode::PrecisionUnsupported) {
                    return std::unexpected(added.error());
                }
                return WitnessAssessment{};
            }
        }
    }
    bool zero = true;
    Accumulator quadratic;
    double propagated = 0;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        const auto value = rows[row].result();
        if (!value) {
            return WitnessAssessment{};
        }
        zero &= value->value == 0 && value->error == 0;
        const auto added = quadratic.product(witness[row], value->value, 1);
        if (!added) {
            return WitnessAssessment{};
        }
        propagated = detail::sparse_up_add(
            propagated, detail::sparse_up_multiply(std::abs(witness[row]), value->error));
    }
    if (zero) {
        return WitnessAssessment{SparseSpdStatus::ZeroMode, std::nullopt};
    }
    const auto value = quadratic.result();
    if (!value) {
        return WitnessAssessment{};
    }
    const double error = detail::sparse_up_add(value->error, propagated);
    const double upper =
        error == 0 ? value->value
                   : std::nextafter(value->value + error, std::numeric_limits<double>::infinity());
    if (std::isfinite(upper) && upper < 0) {
        return WitnessAssessment{SparseSpdStatus::NegativeCurvature, upper};
    }
    return WitnessAssessment{};
}
} // namespace

namespace detail {
class SparsePreparedData {
public:
    SparsePreparedData(
        SparseReal source, SparseTriangularFactors factors, core::PlannedBuffer<double> work_in,
        core::PlannedBuffer<double> rhs_in, core::PlannedBuffer<double> rhs_low_in,
        core::PlannedBuffer<double> residual_in, core::PlannedBuffer<double> correction_in,
        core::PlannedBuffer<double> error_in, core::PlannedBuffer<double> denominator_in,
        core::PlannedBuffer<Accumulator> sums_in,
        core::PlannedBuffer<MagnitudeBounds> magnitudes_in, core::PlannedBuffer<double> last_in,
        BackendMemoryPlan plan_in, BackendAdmission admission_in, core::Index capacity_in,
        MetadataBuffers scales_in) noexcept
        : original(std::move(source)),
          triangular(std::move(factors)),
          work(std::move(work_in)),
          rhs(std::move(rhs_in)),
          rhs_low(std::move(rhs_low_in)),
          residual(std::move(residual_in)),
          correction(std::move(correction_in)),
          error(std::move(error_in)),
          denominator(std::move(denominator_in)),
          sums(std::move(sums_in)),
          magnitudes(std::move(magnitudes_in)),
          last_correction(std::move(last_in)),
          memory_plan(plan_in),
          admission(admission_in),
          capacity(capacity_in),
          scales(std::move(scales_in)) {}
    SparseReal original;
    SparseTriangularFactors triangular;
    core::PlannedBuffer<double> work, rhs, rhs_low, residual, correction, error, denominator;
    core::PlannedBuffer<Accumulator> sums;
    core::PlannedBuffer<MagnitudeBounds> magnitudes;
    core::PlannedBuffer<double> last_correction;
    BackendMemoryPlan memory_plan;
    BackendAdmission admission;
    core::Index capacity;
    MetadataBuffers scales;
    std::atomic_flag busy = ATOMIC_FLAG_INIT;
};
} // namespace detail

core::Result<SparseSpdDiagnostic> inspect_sparse_spd(const SparseOperator& original,
                                                     const Limits& limits,
                                                     core::ResourcePlan& plan) noexcept {
    auto outcome = core::guard_entry(
        operation, core::FailureCode::FactorizationFailed, {},
        [&]() -> core::Result<SparseSpdDiagnostic> {
            const auto checked = check_operator(original, limits);
            if (!checked) {
                return std::unexpected(checked.error());
            }
            const auto& matrix = original.values();
            const auto symmetric = detail::sparse_exact_symmetry(matrix);
            if (!symmetric) {
                return std::unexpected(symmetric.error());
            }
            if (!*symmetric) {
                return SparseSpdDiagnostic{
                    SparseSpdStatus::Asymmetric, original.context(), {}, {}, {}};
            }
            const auto n = matrix.rows();
            if (n == 0) {
                return SparseSpdDiagnostic{
                    SparseSpdStatus::PositiveDefinite, original.context(), 0.0, {}, {}};
            }
            auto direction = plan.allocate<double>(n, operation);
            if (!direction) {
                return std::unexpected(direction.error());
            }
            for (core::Index column = 0; column < n; ++column) {
                bool zero_column = true;
                double diagonal = 0;
                for (auto entry = matrix.column_offsets()[static_cast<std::size_t>(column)];
                     entry < matrix.column_offsets()[static_cast<std::size_t>(column + 1)];
                     ++entry) {
                    const auto index = static_cast<std::size_t>(entry);
                    zero_column &= matrix.values()[index] == 0;
                    if (matrix.row_indices()[index] == column) {
                        diagonal = matrix.values()[index];
                    }
                }
                if (zero_column || diagonal < 0) {
                    direction->values()[static_cast<std::size_t>(column)] = 1;
                    auto witness = DenseReal::copy({n, 1, n}, direction->values(), plan, operation);
                    if (!witness) {
                        return std::unexpected(witness.error());
                    }
                    return SparseSpdDiagnostic{zero_column ? SparseSpdStatus::ZeroMode
                                                           : SparseSpdStatus::NegativeCurvature,
                                               original.context(),
                                               {},
                                               zero_column ? std::optional<double>{}
                                                           : std::optional<double>{diagonal},
                                               std::move(*witness)};
                }
            }
            auto rows = plan.allocate<Accumulator>(n, operation);
            if (!rows) {
                return std::unexpected(rows.error());
            }
            for (unsigned candidate = 0; candidate < 2; ++candidate) {
                for (core::Index row = 0; row < n; ++row) {
                    direction->values()[static_cast<std::size_t>(row)] =
                        candidate == 0 || row % 2 == 0 ? 1.0 : -1.0;
                }
                const auto tested = check_witness(matrix, direction->values(), rows->values());
                if (!tested) {
                    return std::unexpected(tested.error());
                }
                if (tested->status != SparseSpdStatus::Uncertified) {
                    auto witness = DenseReal::copy({n, 1, n}, direction->values(), plan, operation);
                    if (!witness) {
                        return std::unexpected(witness.error());
                    }
                    return SparseSpdDiagnostic{tested->status,
                                               original.context(),
                                               {},
                                               tested->quadratic,
                                               std::move(*witness)};
                }
            }
            const auto forecast = backend_plan(matrix);
            if (!forecast) {
                return std::unexpected(forecast.error());
            }
            const auto admitted = admit_backend_memory(*forecast, limits.backend_memory);
            if (!admitted) {
                return std::unexpected(admitted.error());
            }
            const auto proof = detail::qualify_sparse_positive_definite(matrix, plan);
            if (!proof) {
                if (proof.error().code == core::FailureCode::RankUnresolved ||
                    proof.error().code == core::FailureCode::PrecisionUnresolved) {
                    return SparseSpdDiagnostic{
                        SparseSpdStatus::Uncertified, original.context(), {}, {}, {}};
                }
                return std::unexpected(proof.error());
            }
            return SparseSpdDiagnostic{proof->certified ? SparseSpdStatus::PositiveDefinite
                                                        : SparseSpdStatus::Uncertified,
                                       original.context(),
                                       proof->shift,
                                       {},
                                       {}};
        });
    if (!outcome) {
        // Preserve inner diagnostics and add the operator identities available here.
        auto& details = outcome.error().details;
        if (!details.model) {
            details.model = original.context().cache.model;
        }
        if (!details.trial) {
            details.trial = original.context().trial;
        }
        if (!details.numerical_generation) {
            details.numerical_generation = original.context().cache.matrix;
        }
    }
    return outcome;
}

SparsePrepared::SparsePrepared(std::unique_ptr<detail::SparsePreparedData> data, Context context,
                               SparseMethod method) noexcept
    : data_(std::move(data)), context_(context), method_(method) {}
SparsePrepared::SparsePrepared(SparsePrepared&&) noexcept = default;
SparsePrepared& SparsePrepared::operator=(SparsePrepared&&) noexcept = default;
SparsePrepared::~SparsePrepared() = default;
bool SparsePrepared::active() const noexcept {
    return data_ != nullptr;
}

core::Result<SparsePrepared> SparsePrepared::prepare(const SparseOperator& original,
                                                     const Scales& scales, SparseMethod method,
                                                     const Limits& limits,
                                                     core::ResourcePlan& plan) noexcept {
    auto outcome = core::guard_entry(
        operation, core::FailureCode::FactorizationFailed, {},
        [&]() -> core::Result<SparsePrepared> {
            const auto checked = check_operator(original, limits);
            if (!checked) {
                return std::unexpected(checked.error());
            }
            if (method != SparseMethod::GeneralLu && method != SparseMethod::PositiveDefiniteLlt) {
                return fail(core::FailureCode::InvalidInput);
            }
            const auto scale_check =
                detail::check_scale_inputs(original.values().rows(), scales, original.context());
            if (!scale_check) {
                return std::unexpected(scale_check.error());
            }
            const auto n = original.values().rows();
            const auto forecast = backend_plan(original.values());
            if (!forecast) {
                return std::unexpected(forecast.error());
            }
            const auto admitted = admit_backend_memory(*forecast, limits.backend_memory);
            if (!admitted) {
                return std::unexpected(admitted.error());
            }
            auto scale_snapshot = detail::SolutionFactory::prepare_metadata(
                original.values().rows(), 0, scales, {}, plan);
            if (!scale_snapshot) {
                return std::unexpected(scale_snapshot.error());
            }
            auto snapshot = SparseReal::copy_csc(n, n, original.values().column_offsets(),
                                                 original.values().row_indices(),
                                                 original.values().values(), plan, operation);
            if (!snapshot) {
                return std::unexpected(snapshot.error());
            }
            detail::SparseTriangularFactors factors;
            if (n != 0) {
                const auto mapped =
                    uvlas::numerics::detail::eigen_read_map(*snapshot, n, n, operation);
                if (!mapped) {
                    return std::unexpected(mapped.error());
                }
                if (method == SparseMethod::PositiveDefiniteLlt) {
                    const SparseOperator owned{*snapshot, original.context()};
                    const auto inspected = inspect_sparse_spd(owned, limits, plan);
                    if (!inspected) {
                        return std::unexpected(inspected.error());
                    }
                    if (inspected->status() != SparseSpdStatus::PositiveDefinite) {
                        return fail(core::FailureCode::RankUnresolved);
                    }
                    detail::SparseLltFactor factor;
                    factor.compute(*mapped);
                    auto extracted = detail::extract_sparse_llt(factor, plan);
                    if (!extracted) {
                        return std::unexpected(extracted.error());
                    }
                    factors = std::move(*extracted);
                } else {
                    detail::SparseLuFactor factor;
                    factor.compute(*mapped);
                    auto extracted = detail::extract_sparse_lu(factor, plan);
                    if (!extracted) {
                        return std::unexpected(extracted.error());
                    }
                    factors = std::move(*extracted);
                }
            }
            const auto count =
                core::checked_extent_product(std::array{n, limits.maximum_rhs_columns}, operation);
            if (!count) {
                return std::unexpected(count.error());
            }
            auto work = plan.allocate<double>(*count, operation);
            if (!work) {
                return std::unexpected(work.error());
            }
            auto rhs = plan.allocate<double>(*count, operation);
            if (!rhs) {
                return std::unexpected(rhs.error());
            }
            auto rhs_low = plan.allocate<double>(*count, operation);
            if (!rhs_low) {
                return std::unexpected(rhs_low.error());
            }
            auto residual = plan.allocate<double>(*count, operation);
            if (!residual) {
                return std::unexpected(residual.error());
            }
            auto correction = plan.allocate<double>(*count, operation);
            if (!correction) {
                return std::unexpected(correction.error());
            }
            auto error = plan.allocate<double>(n, operation);
            if (!error) {
                return std::unexpected(error.error());
            }
            auto denominator = plan.allocate<double>(n, operation);
            if (!denominator) {
                return std::unexpected(denominator.error());
            }
            auto sums = plan.allocate<Accumulator>(n, operation);
            if (!sums) {
                return std::unexpected(sums.error());
            }
            auto magnitudes = plan.allocate<detail::MagnitudeBounds>(n, operation);
            if (!magnitudes) {
                return std::unexpected(magnitudes.error());
            }
            auto last = plan.allocate<double>(n == 0 ? 0 : limits.maximum_rhs_columns, operation);
            if (!last) {
                return std::unexpected(last.error());
            }
            auto data = std::make_unique<detail::SparsePreparedData>(
                std::move(*snapshot), std::move(factors), std::move(*work), std::move(*rhs),
                std::move(*rhs_low), std::move(*residual), std::move(*correction),
                std::move(*error), std::move(*denominator), std::move(*sums),
                std::move(*magnitudes), std::move(*last), *forecast, *admitted,
                limits.maximum_rhs_columns, std::move(*scale_snapshot));
            return SparsePrepared{std::move(data), original.context(), method};
        });
    if (!outcome) {
        // Preserve inner diagnostics and add the operator identities available here.
        auto& details = outcome.error().details;
        if (!details.model) {
            details.model = original.context().cache.model;
        }
        if (!details.trial) {
            details.trial = original.context().trial;
        }
        if (!details.numerical_generation) {
            details.numerical_generation = original.context().cache.matrix;
        }
    }
    return outcome;
}

namespace {
class ScratchLease {
public:
    explicit ScratchLease(std::atomic_flag& flag) noexcept
        : flag_(flag), acquired_(!flag.test_and_set(std::memory_order_acquire)) {}
    ~ScratchLease() {
        if (acquired_) {
            flag_.clear(std::memory_order_release);
        }
    }
    ScratchLease(const ScratchLease&) = delete;
    ScratchLease& operator=(const ScratchLease&) = delete;
    [[nodiscard]] bool acquired() const noexcept {
        return acquired_;
    }

private:
    std::atomic_flag& flag_;
    bool acquired_;
};

core::Result<bool> evaluate_sparse_columns(detail::SparsePreparedData& data,
                                           std::span<const double> x, std::span<const double> x_low,
                                           core::Index columns, const Scales& scales,
                                           const Accuracy& accuracy,
                                           std::span<ColumnMetrics> metrics) noexcept {
    const auto n = data.original.rows();
    const auto& matrix = data.original;
    bool achieved = true;
    for (core::Index column = 0; column < columns; ++column) {
        const auto start = static_cast<std::size_t>(column * n);
        std::fill(data.sums.values().begin(), data.sums.values().end(), Accumulator{});
        for (core::Index row = 0; row < n; ++row) {
            const auto index = static_cast<std::size_t>(row);
            const double b = data.rhs.values()[start + index];
            auto added = data.sums.values()[index].add(b);
            if (!added) {
                return std::unexpected(added.error());
            }
            const double low = data.rhs_low.values()[start + index];
            added = data.sums.values()[index].add(low);
            if (!added) {
                return std::unexpected(added.error());
            }
            const auto magnitude = detail::actual_magnitude_bounds(b, low);
            if (!magnitude) {
                return std::unexpected(magnitude.error());
            }
            data.magnitudes.values()[index] = *magnitude;
        }
        for (core::Index coordinate = 0; coordinate < n; ++coordinate) {
            const double value = x[start + static_cast<std::size_t>(coordinate)];
            const double low =
                x_low.empty() ? 0 : x_low[start + static_cast<std::size_t>(coordinate)];
            const auto absolute = detail::actual_magnitude_bounds(value, low);
            if (!absolute) {
                return std::unexpected(absolute.error());
            }
            for (auto entry = matrix.column_offsets()[static_cast<std::size_t>(coordinate)];
                 entry < matrix.column_offsets()[static_cast<std::size_t>(coordinate + 1)];
                 ++entry) {
                const auto index = static_cast<std::size_t>(entry);
                const auto row = static_cast<std::size_t>(matrix.row_indices()[index]);
                const auto product =
                    data.sums.values()[row].product(matrix.values()[index], value, -1);
                if (!product) {
                    return std::unexpected(product.error());
                }
                const auto low_product =
                    data.sums.values()[row].product(matrix.values()[index], low, -1);
                if (!low_product) {
                    return std::unexpected(low_product.error());
                }
                const double a = std::abs(matrix.values()[index]);
                data.magnitudes.values()[row] = detail::add_magnitudes(
                    data.magnitudes.values()[row], detail::multiply_magnitudes({a, a}, *absolute));
            }
        }
        for (core::Index row = 0; row < n; ++row) {
            const auto index = static_cast<std::size_t>(row);
            const auto sum = data.sums.values()[index].result();
            if (!sum) {
                return std::unexpected(sum.error());
            }
            data.residual.values()[start + index] = sum->value;
            data.error.values()[index] = sum->error;
            data.denominator.values()[index] = data.magnitudes.values()[index].lower;
        }
        const auto metric = detail::evaluate_column_metrics<double>(
            n, data.residual.values().subspan(start, static_cast<std::size_t>(n)),
            data.error.values(), data.denominator.values(),
            data.rhs.values().subspan(start, static_cast<std::size_t>(n)),
            data.rhs_low.values().subspan(start, static_cast<std::size_t>(n)), scales.residual,
            accuracy);
        if (!metric) {
            return std::unexpected(metric.error());
        }
        metrics[static_cast<std::size_t>(column)] = *metric;
        metrics[static_cast<std::size_t>(column)].scaled_correction_upper =
            n == 0 ? 0 : data.last_correction.values()[static_cast<std::size_t>(column)];
        achieved &= detail::meets_accuracy(*metric, accuracy);
    }
    return achieved;
}
} // namespace

core::Result<Solution<double>>
solve_linear(const SparseOperator& original, const Rhs<double>& rhs, const Scales& scales,
             const Precision<double>& precision, const Accuracy& accuracy, const Limits& limits,
             SparsePrepared& prepared, core::ResourcePlan& result_plan) noexcept {
    auto outcome = core::guard_entry(
        operation, core::FailureCode::FactorizationFailed, {},
        [&]() -> core::Result<Solution<double>> {
            const auto checked = check_operator(original, limits);
            if (!checked) {
                return std::unexpected(checked.error());
            }
            if (!prepared.active()) {
                return fail(core::FailureCode::StaleGeneration);
            }
            auto& data = *prepared.data_;
            if (original.context() != prepared.context_ ||
                !same_matrix(original.values(), data.original)) {
                return fail(core::FailureCode::StaleGeneration);
            }
            const auto n = data.original.rows();
            const auto valid = detail::check_solve_inputs(n, 0, original.context(), rhs, scales,
                                                          precision, accuracy, limits);
            if (!valid) {
                return std::unexpected(valid.error());
            }
            // Scale labels alone cannot authorize replacement by different actual scales/units.
            for (std::size_t i = 0; i < scales.unknown.size(); ++i) {
                if (std::bit_cast<std::uint64_t>(scales.unknown[i]) !=
                        std::bit_cast<std::uint64_t>(data.scales.unknown_scales.values()[i]) ||
                    std::bit_cast<std::uint64_t>(scales.residual[i]) !=
                        std::bit_cast<std::uint64_t>(data.scales.residual_scales.values()[i]) ||
                    scales.unknown_units[i] != data.scales.unknown_units.values()[i] ||
                    scales.residual_units[i] != data.scales.residual_units.values()[i]) {
                    return fail(core::FailureCode::StaleGeneration);
                }
            }
            const auto* binary_source = std::get_if<Binary64<double>>(&rhs.values().storage());
            const auto* pair_source = std::get_if<TwoWord<double>>(&rhs.values().storage());
            const auto& source_hi = binary_source ? binary_source->values : pair_source->hi;
            const auto* source_lo = pair_source ? &pair_source->lo : nullptr;
            const auto shape = source_hi.shape();
            if (shape.columns > data.capacity) {
                return fail(core::FailureCode::ResourceExceeded);
            }
            const auto admitted = admit_backend_memory(data.memory_plan, limits.backend_memory);
            if (!admitted) {
                return std::unexpected(admitted.error());
            }
            const ScratchLease lease{data.busy};
            if (!lease.acquired()) {
                return fail(core::FailureCode::InvalidInput);
            }
            const auto extent =
                core::checked_extent_product(std::array{n, shape.columns}, operation);
            if (!extent) {
                return std::unexpected(extent.error());
            }
            const auto count = static_cast<std::size_t>(*extent);
            std::ranges::fill(data.work.values().first(count), 0);
            auto values = DenseReal::copy({n, shape.columns, n}, data.work.values().first(count),
                                          result_plan, operation);
            if (!values) {
                return std::unexpected(values.error());
            }
            std::optional<DenseReal> projected_values;
            if (std::holds_alternative<Binary64Projection>(precision.returned)) {
                auto output = DenseReal::copy(
                    {n, shape.columns, n}, data.work.values().first(count), result_plan, operation);
                if (!output) {
                    return std::unexpected(output.error());
                }
                projected_values.emplace(std::move(*output));
            }
            auto low_values = DenseReal::copy(
                {n, shape.columns, n}, data.work.values().first(count), result_plan, operation);
            if (!low_values) {
                return std::unexpected(low_values.error());
            }
            auto metadata = detail::SolutionFactory::prepare_metadata(
                n, shape.columns, scales, precision.required_checks, result_plan);
            if (!metadata) {
                return std::unexpected(metadata.error());
            }
            std::fill(data.last_correction.values().begin(), data.last_correction.values().end(),
                      0);
            for (core::Index column = 0; column < shape.columns; ++column) {
                for (core::Index row = 0; row < n; ++row) {
                    data.rhs.values()[static_cast<std::size_t>(column * n + row)] =
                        source_hi.values()[static_cast<std::size_t>(
                            column * shape.leading_dimension + row)];
                    data.rhs_low.values()[static_cast<std::size_t>(column * n + row)] =
                        source_lo ? source_lo->values()[static_cast<std::size_t>(
                                        column * shape.leading_dimension + row)]
                                  : 0;
                }
            }
            auto solved = detail::solve_sparse_triangular(
                data.triangular, data.rhs.values().first(count), shape.columns, data.work.values(),
                values->values());
            if (!solved) {
                return std::unexpected(solved.error());
            }
            bool two = pair_source != nullptr;
            RetentionPath retention =
                two ? RetentionPath::InitiallyTwoWord : RetentionPath::Binary64;
            if (two) {
                solved = detail::solve_sparse_triangular(
                    data.triangular, data.rhs_low.values().first(count), shape.columns,
                    data.work.values(), data.correction.values().first(count));
                if (!solved) {
                    return std::unexpected(solved.error());
                }
                for (std::size_t i = 0; i < count; ++i) {
                    const auto update = uvlas::numerics::detail::correct_two_word(
                        values->values()[i], low_values->values()[i], data.correction.values()[i]);
                    if (!update) {
                        return std::unexpected(update.error());
                    }
                    values->values()[i] = update->hi;
                    low_values->values()[i] = update->lo;
                }
            }
            unsigned corrections = 0;
            while (true) {
                const auto metrics = evaluate_sparse_columns(
                    data, values->values(), two ? low_values->values() : std::span<const double>{},
                    shape.columns, scales, accuracy, metadata->retained_columns.values());
                if (!metrics) {
                    return std::unexpected(metrics.error());
                }
                if (*metrics) {
                    break;
                }
                if (corrections == limits.maximum_refinements) {
                    return fail(core::FailureCode::LinearUnconverged);
                }
                if (!two && corrections > 0) {
                    for (std::size_t c = 0; c < metadata->retained_columns.size(); ++c) {
                        const auto now = metadata->retained_columns.values()[c],
                                   before = metadata->returned_columns.values()[c];
                        if (!detail::meets_accuracy(now, accuracy) &&
                            (now.normalized_original_residual_upper <= 1 ||
                             now.normalized_original_residual_upper >
                                 0.5 * before.normalized_original_residual_upper) &&
                            (now.componentwise_backward_error_upper <=
                                 accuracy.componentwise_backward_tolerance ||
                             now.componentwise_backward_error_upper >
                                 0.5 * before.componentwise_backward_error_upper)) {
                            two = true;
                            retention = RetentionPath::PromotedTwoWord;
                        }
                    }
                }
                std::ranges::copy(metadata->retained_columns.values(),
                                  metadata->returned_columns.values().begin());
                solved = detail::solve_sparse_triangular(
                    data.triangular, data.residual.values().first(count), shape.columns,
                    data.work.values(), data.correction.values().first(count));
                if (!solved) {
                    return std::unexpected(solved.error());
                }
                ++corrections;
                if (!two) {
                    bool stagnated = false;
                    for (core::Index column = 0; column < shape.columns; ++column) {
                        bool unchanged = true;
                        for (core::Index row = 0; row < n; ++row) {
                            const auto i = static_cast<std::size_t>(column * n + row);
                            unchanged &= values->values()[i] + data.correction.values()[i] ==
                                         values->values()[i];
                        }
                        stagnated |=
                            unchanged &&
                            !detail::meets_accuracy(metadata->retained_columns
                                                        .values()[static_cast<std::size_t>(column)],
                                                    accuracy);
                    }
                    if (stagnated) {
                        two = true;
                        retention = RetentionPath::PromotedTwoWord;
                    }
                }
                bool changed = false;
                for (core::Index column = 0; column < shape.columns; ++column) {
                    for (core::Index row = 0; row < n; ++row) {
                        const auto index = static_cast<std::size_t>(column * n + row);
                        const double correction = data.correction.values()[index];
                        const double updated = values->values()[index] + correction;
                        if (!std::isfinite(updated)) {
                            return fail(core::FailureCode::Nonfinite);
                        }
                        if (two) {
                            const auto update = uvlas::numerics::detail::correct_two_word(
                                values->values()[index], low_values->values()[index], correction);
                            if (!update) {
                                return std::unexpected(update.error());
                            }
                            changed |= update->hi != values->values()[index] ||
                                       update->lo != low_values->values()[index];
                            values->values()[index] = update->hi;
                            low_values->values()[index] = update->lo;
                        } else {
                            changed |= updated != values->values()[index];
                            values->values()[index] = updated;
                        }
                    }
                    std::ranges::fill(data.error.values(), 0);
                    std::ranges::fill(data.denominator.values(), 1);
                    auto correction_accuracy = accuracy;
                    correction_accuracy.absolute_tolerance = 1;
                    correction_accuracy.relative_tolerance = 0;
                    const auto offset = static_cast<std::size_t>(column * n),
                               size = static_cast<std::size_t>(n);
                    const auto correction_metric = detail::evaluate_column_metrics<double>(
                        n, data.correction.values().subspan(offset, size), data.error.values(),
                        data.denominator.values(), data.rhs.values().subspan(offset, size),
                        data.rhs_low.values().subspan(offset, size), scales.unknown,
                        correction_accuracy);
                    if (!correction_metric) {
                        return std::unexpected(correction_metric.error());
                    }
                    if (n != 0) {
                        data.last_correction.values()[static_cast<std::size_t>(column)] =
                            correction_metric->normalized_original_residual_upper;
                    }
                }
                if (!changed) {
                    return fail(core::FailureCode::PrecisionUnresolved);
                }
            }
            const auto generation = std::visit(
                [](const auto& request) {
                    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(request)>,
                                                 RetainedPrecision>) {
                        return request.generation;
                    } else {
                        return request.retained_generation;
                    }
                },
                precision.returned);
            auto retained = two ? PrecisionBuffer<double>::two_word(
                                      std::move(*values), std::move(*low_values), generation)
                                : PrecisionBuffer<double>::binary64(std::move(*values), generation);
            if (!retained) {
                return std::unexpected(retained.error());
            }
            std::optional<PrecisionBuffer<double>> projected;
            if (projected_values) {
                const auto* binary = std::get_if<Binary64<double>>(&retained->storage());
                const auto* pair = std::get_if<TwoWord<double>>(&retained->storage());
                for (std::size_t i = 0; i < count; ++i) {
                    const double projected_value =
                        binary ? binary->values.values()[i]
                               : pair->hi.values()[i] + pair->lo.values()[i];
                    if (!std::isfinite(projected_value)) {
                        return fail(core::FailureCode::Nonfinite);
                    }
                    projected_values->values()[i] = projected_value;
                }
                const auto checked_projection =
                    evaluate_sparse_columns(data, projected_values->values(), {}, shape.columns,
                                            scales, accuracy, metadata->returned_columns.values());
                if (!checked_projection) {
                    return std::unexpected(checked_projection.error());
                }
                if (!*checked_projection) {
                    return fail(core::FailureCode::PrecisionUnresolved);
                }
                auto output = PrecisionBuffer<double>::binary64(
                    std::move(*projected_values),
                    std::get<Binary64Projection>(precision.returned).target_generation);
                if (!output) {
                    return std::unexpected(output.error());
                }
                projected.emplace(std::move(*output));
            }
            const detail::Finalization context{original.context(),
                                               rhs.cache_request().rhs,
                                               rhs.values().generation(),
                                               accuracy,
                                               data.admission,
                                               *admitted,
                                               corrections,
                                               retention,
                                               NotEvaluated{}};
            return detail::SolutionFactory::finalize(std::move(*retained), std::move(projected),
                                                     precision, std::move(*metadata), context);
        });
    if (!outcome) {
        // Preserve inner diagnostics and add the operator identities available here.
        auto& details = outcome.error().details;
        if (!details.model) {
            details.model = original.context().cache.model;
        }
        if (!details.trial) {
            details.trial = original.context().trial;
        }
        if (!details.numerical_generation) {
            details.numerical_generation = original.context().cache.matrix;
        }
    }
    return outcome;
}

} // namespace uvlas::numerics::linear
