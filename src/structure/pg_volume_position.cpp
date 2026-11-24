#include "pg_volume_position_arithmetic.hpp"

#include <array>
#include <cstring>
#include <uvlas/structure/pg_volume_position.hpp>

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

void add_term(VolumeInteger& value, VolumeInteger& reference_error, PgVolumeMassWords basis,
              std::size_t coefficient, PgCartesianWords input, std::size_t original,
              VolumeWord epsilon) {
    const auto n_hi = volume_decode(basis.high[coefficient]);
    const auto n_lo = volume_decode(low_at(basis.low, coefficient));
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

[[nodiscard]] bool fits_steps(std::uint64_t steps, PgCartesianWords input, std::size_t inner,
                              std::size_t output) noexcept {
    for (const auto count : {input.high.size(), input.low.size(), output}) {
        if (count > steps)
            return false;
        steps -= count;
    }
    return inner == 0 || output <= steps / inner;
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

core::Result<void> apply_pg_volume_position(const PgVolumeMassReference& reference,
                                            PgVolumePositionDirection direction,
                                            PgCartesianWords input,
                                            const PgVolumePositionLimits& limits,
                                            core::Workspace& workspace,
                                            PgVolumePositionResult& output) noexcept {
    return core::guard_entry(
        volume_operation, core::FailureCode::InvalidInput, {}, [&]() -> core::Result<void> {
            if (limits.cancellation.stop_requested())
                return reject(core::FailureCode::Cancelled);
            if (const auto accepted = reference.conditioning().info(); !accepted)
                return std::unexpected(accepted.error());
            if (!valid_precision(input.representation))
                return reject(core::FailureCode::PrecisionUnsupported);
            const auto input_paired = input.representation == core::Precision::TwoWord;
            const auto output_paired = output.representation_ == core::Precision::TwoWord;
            if (!volume_finite(limits.maximum_absolute_error))
                return reject(core::FailureCode::Nonfinite);
            const auto allowance = volume_lift(volume_decode(limits.maximum_absolute_error));
            if (allowance < 0)
                return reject(core::FailureCode::InvalidInput);
            const auto forward = direction == PgVolumePositionDirection::ToPoints;
            if (!forward && direction != PgVolumePositionDirection::TransposeToColumns)
                return reject(core::FailureCode::InvalidInput);
            const auto rank = reference.column_count();
            const auto points = reference.points().size();
            const auto inner = forward ? rank : points;
            const auto outer = forward ? points : rank;
            if (std::max(inner, outer) > std::numeric_limits<std::size_t>::max() / 3)
                return reject(core::FailureCode::ResourceExceeded);
            const auto input_size = 3 * inner;
            const auto output_size = 3 * outer;
            if (input.high.size() != input_size ||
                input.low.size() != (input_paired ? input_size : 0))
                return reject(core::FailureCode::InvalidInput);
            if (output_size > output.high_.size() ||
                !fits_steps(limits.maximum_steps, input, inner, output_size))
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

            const auto basis = reference.words();
            const auto epsilon = volume_decode(reference.projection().maximum_shape_error);
            const auto stride = rank + 1;
            for (std::size_t i = 0; i < outer; ++i)
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    VolumeInteger sum{}, reference_error{};
                    if (limits.cancellation.stop_requested())
                        return reject(core::FailureCode::Cancelled);
                    for (std::size_t k = 0; k < inner; ++k) {
                        if (limits.cancellation.stop_requested())
                            return reject(core::FailureCode::Cancelled);
                        const auto coefficient = forward ? i * stride + k + 1 : k * stride + i + 1;
                        add_term(sum, reference_error, basis, coefficient, input, 3 * k + axis,
                                 epsilon);
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
                std::memcpy(output.high_.values().data(), high.data(), bytes);
                if (output_paired)
                    std::memcpy(output.low_.values().data(), low.data(), bytes);
            }
            output.size_ = static_cast<std::size_t>(output_size);
            return {};
        });
}

} // namespace uvlas::structure
