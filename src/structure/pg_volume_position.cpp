#include "pg_volume_position_arithmetic.hpp"

#include <array>
#include <cstring>
#include <optional>
#include <uvlas/structure/pg_boundary_tangents.hpp>
#include <uvlas/structure/pg_material_tangents.hpp>
#include <uvlas/structure/pg_prepared_volume.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;

[[nodiscard]] bool two_word(core::Precision precision) {
    if (precision == core::Precision::Binary64)
        return false;
    if (precision == core::Precision::TwoWord)
        return true;
    partition_fail(core::FailureCode::PrecisionUnsupported);
}

[[nodiscard]] bool valid_precision(core::Precision precision) noexcept {
    return precision == core::Precision::Binary64 || precision == core::Precision::TwoWord;
}

[[nodiscard]] std::unexpected<core::Failure> reject(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, volume_operation});
}

[[nodiscard]] double low_at(std::span<const double> plane, std::size_t i) noexcept {
    return plane.empty() ? 0.0 : plane[i];
}

void add_term(VolumeInteger& value, VolumeInteger& reference_error, double high, double low,
              PgCartesianWords input, std::size_t original, VolumeWord epsilon) {
    const auto n_hi = volume_decode(high);
    const auto n_lo = volume_decode(low);
    const auto q_hi = volume_decode(input.high[original]);
    const auto q_lo = volume_decode(low_at(input.low, original));
    value += volume_product(n_hi, q_hi);
    value += volume_product(n_lo, q_hi);
    value += volume_product(n_hi, q_lo);
    value += volume_product(n_lo, q_lo);
    reference_error += rank_abs(volume_product(epsilon, q_hi) + volume_product(epsilon, q_lo));
}

void stage_word(std::span<std::byte> bytes, std::size_t component, double value) noexcept {
    std::memcpy(bytes.data() + component * sizeof(double), &value, sizeof(value));
}

[[nodiscard]] core::Result<void> validate_plane(std::span<const double> plane,
                                                const PgVolumePositionLimits& limits) noexcept {
    for (const auto value : plane) {
        if (limits.cancellation.stop_requested())
            return reject(core::FailureCode::Cancelled);
        if (!volume_finite(value))
            return reject(core::FailureCode::Nonfinite);
    }
    return {};
}

[[nodiscard]] std::optional<std::uint64_t> checked_steps(std::uint64_t limit,
                                                         PgCartesianWords input, std::size_t inner,
                                                         std::size_t output) noexcept {
    auto remaining = limit;
    for (const auto count : {input.high.size(), input.low.size(), output}) {
        if (count > remaining)
            return {};
        remaining -= count;
    }
    if (inner != 0 && output > remaining / inner)
        return {};
    return (limit - remaining) + static_cast<std::uint64_t>(output) * inner;
}

} // namespace

PgVolumePositionResult::PgVolumePositionResult(core::Precision precision, core::Index bytes,
                                               core::PlannedBuffer<double> high,
                                               core::PlannedBuffer<double> low) noexcept
    : representation_(precision),
      workspace_bytes_(bytes),
      high_(std::move(high)),
      low_(std::move(low)) {}

core::Result<PgVolumePositionResult>
PgVolumePositionResult::allocate(core::Index capacity, core::Precision precision,
                                 core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgVolumePositionResult> {
        const auto paired = two_word(precision);
        const auto bytes = partition_take(core::checked_byte_count(
            capacity, sizeof(double) * (paired ? 2U : 1U), partition_operation));
        const auto scratch = partition_index(bytes);
        if (bytes > plan.remaining_bytes())
            partition_fail(core::FailureCode::ResourceExceeded);
        auto high = partition_take(plan.allocate<double>(capacity, partition_operation));
        auto low =
            partition_take(plan.allocate<double>(paired ? capacity : 0, partition_operation));
        return PgVolumePositionResult{precision, scratch, std::move(high), std::move(low)};
    });
}

core::Index PgVolumePositionResult::workspace_bytes() const noexcept {
    return high_.size() == 0 ? 0 : workspace_bytes_;
}

PgCartesianWords PgVolumePositionResult::words() const& noexcept {
    const auto size = std::min(size_, high_.size());
    return {representation_, high_.values().first(size),
            low_.values().empty() ? std::span<const double>{} : low_.values().first(size)};
}

core::Result<void> PgVolumePositionResult::apply(const MatrixView& reference,
                                                 PgVolumePositionDirection direction,
                                                 PgCartesianWords input,
                                                 const PgVolumePositionLimits& limits,
                                                 core::Workspace& workspace) noexcept {
    return core::guard_entry(
        volume_operation, core::FailureCode::InvalidInput, {}, [&]() -> core::Result<void> {
            if (limits.cancellation.stop_requested())
                return reject(core::FailureCode::Cancelled);
            if (reference.certificate != nullptr) {
                if (const auto accepted = reference.certificate->info(); !accepted)
                    return std::unexpected(accepted.error());
            } else if (reference.rank == 0) {
                return reject(core::FailureCode::InvalidInput);
            }
            if (!valid_precision(input.representation))
                return reject(core::FailureCode::PrecisionUnsupported);
            const auto input_paired = input.representation == core::Precision::TwoWord;
            const auto output_paired = representation_ == core::Precision::TwoWord;
            if (!volume_finite(limits.maximum_absolute_error))
                return reject(core::FailureCode::Nonfinite);
            const auto allowance = volume_lift(volume_decode(limits.maximum_absolute_error));
            if (allowance < 0)
                return reject(core::FailureCode::InvalidInput);
            const auto forward = direction == PgVolumePositionDirection::ToPoints;
            if (!forward && direction != PgVolumePositionDirection::TransposeToColumns)
                return reject(core::FailureCode::InvalidInput);
            const auto rank = reference.rank;
            const auto points = reference.rows;
            const auto inner = forward ? rank : points;
            const auto outer = forward ? points : rank;
            if (std::max(inner, outer) > std::numeric_limits<std::size_t>::max() / 3)
                return reject(core::FailureCode::ResourceExceeded);
            const auto input_size = 3 * inner;
            const auto output_size = 3 * outer;
            if (input.high.size() != input_size ||
                input.low.size() != (input_paired ? input_size : 0))
                return reject(core::FailureCode::InvalidInput);
            const auto steps = checked_steps(limits.maximum_steps, input, inner, output_size);
            if (output_size > high_.size() || !steps)
                return reject(core::FailureCode::ResourceExceeded);

            auto lease = workspace.acquire();
            if (!lease)
                return std::unexpected(lease.error());
            // Capacity was admitted in allocate(), including both word planes.
            const auto bytes = output_size * sizeof(double);
            auto scratch = lease->bytes();
            if (bytes * (output_paired ? 2U : 1U) > scratch.size())
                return reject(core::FailureCode::ResourceExceeded);
            const auto high = scratch.first(bytes);
            const auto low = output_paired ? scratch.subspan(bytes, bytes) : std::span<std::byte>{};
            if (const auto valid = validate_plane(input.high, limits); !valid)
                return valid;
            if (const auto valid = validate_plane(input.low, limits); !valid)
                return valid;

            const auto epsilon = volume_decode(reference.maximum_coefficient_error);
            const auto stride = rank * reference.column_stride + reference.prefix;
            for (std::size_t i = 0; i < outer; ++i)
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    VolumeInteger sum{}, reference_error{};
                    if (limits.cancellation.stop_requested())
                        return reject(core::FailureCode::Cancelled);
                    for (std::size_t k = 0; k < inner; ++k) {
                        if (limits.cancellation.stop_requested())
                            return reject(core::FailureCode::Cancelled);
                        const auto coefficient =
                            (forward ? i * stride + k * reference.column_stride
                                     : k * stride + i * reference.column_stride) +
                            reference.prefix;
                        add_term(sum, reference_error, reference.high[coefficient],
                                 low_at(reference.low, coefficient), input, 3 * k + axis, epsilon);
                    }
                    const auto value =
                        volume_finalize(std::move(sum), reference_error, allowance, output_paired);
                    if (!value)
                        return std::unexpected(value.error());
                    const auto component = 3 * i + axis;
                    stage_word(high, component, value->high);
                    if (output_paired)
                        stage_word(low, component, value->low);
                }
            if (limits.cancellation.stop_requested())
                return reject(core::FailureCode::Cancelled);
            if (bytes != 0) {
                std::memcpy(high_.values().data(), high.data(), bytes);
                if (output_paired)
                    std::memcpy(low_.values().data(), low.data(), bytes);
            }
            size_ = output_size;
            used_steps_ = *steps;
            return {};
        });
}

core::Result<void> apply_pg_volume_position(const PgVolumeMassReference& reference,
                                            PgVolumePositionDirection direction,
                                            PgCartesianWords input,
                                            const PgVolumePositionLimits& limits,
                                            core::Workspace& workspace,
                                            PgVolumePositionResult& output) noexcept {
    const auto words = reference.words();
    return output.apply({&reference.conditioning(), reference.column_count(),
                         reference.points().size(), words.high, words.low, 1,
                         reference.projection().maximum_shape_error},
                        direction, input, limits, workspace);
}

core::Result<void> apply_pg_boundary_position(const PgBoundaryPositionReference& reference,
                                              PgVolumePositionDirection direction,
                                              PgCartesianWords input,
                                              const PgVolumePositionLimits& limits,
                                              core::Workspace& workspace,
                                              PgVolumePositionResult& output) noexcept {
    const auto words = reference.words();
    return output.apply({&reference.conditioning(), reference.column_count(),
                         reference.points().size(), words.high, words.low, 0,
                         reference.projection().maximum_shape_error},
                        direction, input, limits, workspace);
}

core::Result<void> apply_pg_boundary_tangents(const PgBoundaryTangentReference& reference,
                                              PgVolumePositionDirection direction,
                                              PgCartesianWords input,
                                              const PgVolumePositionLimits& limits,
                                              core::Workspace& workspace,
                                              PgVolumePositionResult& output) noexcept {
    const auto words = reference.words();
    // Preparation admitted both rows per point before allocating their words.
    return output.apply({nullptr, reference.column_count(), 2 * reference.points().size(),
                         words.high, words.low, 0, reference.projection().maximum_derivative_error},
                        direction, input, limits, workspace);
}

core::Result<void> apply_pg_prepared_volume_position(const PgPreparedVolumeReference& reference,
                                                     PgPreparedPatchIdentity source,
                                                     PgVolumePositionDirection direction,
                                                     PgCartesianWords input,
                                                     const PgVolumePositionLimits& limits,
                                                     core::Workspace& workspace,
                                                     PgVolumePositionResult& output) noexcept {
    if (limits.cancellation.stop_requested())
        return reject(core::FailureCode::Cancelled);
    if (source != reference.source())
        return reject(core::FailureCode::StaleGeneration);
    const auto words = reference.words();
    return output.apply({nullptr, reference.column_count(), reference.points().size(), words.high,
                         words.low, 1, reference.projection().maximum_shape_error},
                        direction, input, limits, workspace);
}

core::Result<void> apply_pg_prepared_reference_position(const PgProjectedReferenceBatch& reference,
                                                        PgPreparedPatchIdentity source,
                                                        PgVolumePositionDirection direction,
                                                        PgCartesianWords input,
                                                        const PgVolumePositionLimits& limits,
                                                        core::Workspace& workspace,
                                                        PgVolumePositionResult& output) noexcept {
    if (limits.cancellation.stop_requested())
        return reject(core::FailureCode::Cancelled);
    if (source != reference.source())
        return reject(core::FailureCode::StaleGeneration);
    // X0[3] and measure precede N and the three gradient words per column.
    return output.apply({nullptr, reference.column_count(), reference.points().size(),
                         reference.high_words(), reference.low_words(), 4,
                         reference.projection().maximum_shape_error, 4},
                        direction, input, limits, workspace);
}

core::Result<void>
apply_pg_material_tangents(const PgMaterialTangentReference& reference,
                           PgPreparedPatchIdentity source, PgVolumePositionDirection direction,
                           PgCartesianWords input, const PgVolumePositionLimits& limits,
                           core::Workspace& workspace, PgVolumePositionResult& output) noexcept {
    if (limits.cancellation.stop_requested())
        return reject(core::FailureCode::Cancelled);
    if (source != reference.source_)
        return reject(core::FailureCode::StaleGeneration);
    return output.apply({nullptr, reference.column_count_, 2 * reference.points_.size(),
                         reference.high_.values(), reference.low_.values(), 0,
                         reference.projection_.maximum_derivative_error},
                        direction, input, limits, workspace);
}

} // namespace uvlas::structure
