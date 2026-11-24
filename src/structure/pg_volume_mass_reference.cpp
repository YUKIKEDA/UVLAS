#include "pg_reference_projection.hpp"

#include <uvlas/structure/pg_volume_mass_reference.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;

void validate_tolerance(double tolerance) {
    const auto bits = std::bit_cast<std::uint64_t>(tolerance);
    const auto magnitude = bits & ~pg_binary64_sign;
    if (magnitude >= pg_binary64_infinity)
        partition_fail(core::FailureCode::Nonfinite);
    if (magnitude != 0 && bits != magnitude)
        partition_fail(core::FailureCode::InvalidInput);
}

void project_mass_point(ReferenceArithmetic& math, const PgConditionedReferenceBatch& point,
                        const PgVolumeMassProjection& request, std::span<double> high,
                        std::span<double> low) {
    const auto& reference = point.points.values().front();
    const auto& jacobian = std::get<PgReferenceVolumeMeasure>(reference.measure).jacobian;
    const auto tolerance_measure =
        reference_binary64(math, request.maximum_volume_jacobian_error_metres_cubed);
    const auto measure =
        pg_project_scalar(math, jacobian, reference.g, tolerance_measure, request.representation);
    // A positive high word implies a positive retained value: the low word
    // rounds a residual of at most half an ulp of a finite rounded high word.
    const auto bits = std::bit_cast<std::uint64_t>(measure.high);
    if (bits == 0 || (bits & pg_binary64_sign) != 0)
        partition_fail(core::FailureCode::PrecisionUnresolved);
    high.front() = measure.high;
    if (!low.empty())
        low.front() = measure.low;
    const auto tolerance_shape = reference_binary64(math, request.maximum_shape_error);
    const auto columns = point.columns.values();
    for (std::size_t j = 0; j < columns.size(); ++j) {
        math.visit(1);
        const auto words = pg_project_scalar(math, columns[j].value, reference.g, tolerance_shape,
                                             request.representation);
        high[j + 1] = words.high;
        if (!low.empty())
            low[j + 1] = words.low;
    }
}
} // namespace

PgVolumeMassReference::PgVolumeMassReference(PgConditioningCertificate&& conditioning,
                                             std::size_t columns, PgVolumeMassProjection projection,
                                             PgReferenceCallInfo info,
                                             core::PlannedBuffer<PgVolumeMassPoint> points,
                                             core::PlannedBuffer<double> high,
                                             core::PlannedBuffer<double> low) noexcept
    : conditioning_(std::move(conditioning)),
      column_count_(columns),
      projection_(projection),
      call_info_(info),
      points_(std::move(points)),
      high_(std::move(high)),
      low_(std::move(low)) {}

core::Result<PgVolumeMassReference>
PgVolumeMassReference::prepare(PgConditioningCertificate&& conditioning,
                               std::span<const PgReferenceVolumePoint> queries,
                               const PgVolumeMassProjection& request,
                               const PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgVolumeMassReference> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto accepted = partition_take(conditioning.info());
        const auto two_word = request.representation == core::Precision::TwoWord;
        if (!two_word && request.representation != core::Precision::Binary64)
            partition_fail(core::FailureCode::PrecisionUnsupported);
        validate_tolerance(request.maximum_shape_error);
        validate_tolerance(request.maximum_volume_jacobian_error_metres_cubed);
        const auto rank = static_cast<std::size_t>(accepted.rank);
        if (queries.empty())
            return PgVolumeMassReference{
                std::move(conditioning), rank, request, math.info(), {}, {}, {}};
        const auto rows = partition_index(queries.size());
        const auto stride = partition_add(accepted.rank, 1);
        const auto word_count = partition_index(partition_mul(queries.size(), stride));
        const auto point_bytes = partition_take(
            core::checked_byte_count(rows, sizeof(PgVolumeMassPoint), partition_operation));
        const auto word_bytes = partition_take(
            core::checked_byte_count(word_count, sizeof(double), partition_operation));
        if (partition_add(point_bytes, partition_mul(two_word ? 2 : 1, word_bytes)) >
            plan.remaining_bytes())
            partition_fail(core::FailureCode::ResourceExceeded);
        // At most one exact reference row remains alive alongside the rounding
        // temporaries, using the same slot envelope as the exact point batch.
        math.admit(partition_add(reference_live_rationals + 74, partition_mul(8, accepted.rank)));
        const auto words = static_cast<std::uint64_t>(word_count);
        math.visit(partition_add(queries.size(), partition_mul(two_word ? 2 : 1, words)));
        auto points = partition_take(plan.allocate<PgVolumeMassPoint>(rows, partition_operation));
        auto high = partition_take(plan.allocate<double>(word_count, partition_operation));
        auto low =
            partition_take(plan.allocate<double>(two_word ? word_count : 0, partition_operation));
        auto high_row = high.values();
        auto low_row = low.values();
        const auto width = static_cast<std::size_t>(stride);
        for (std::size_t p = 0; p < queries.size(); ++p) {
            math.visit(1);
            auto remaining = limits;
            remaining.maximum_scalar_operations -= math.scalar_operations;
            remaining.maximum_partition_visits -= math.partition_visits;
            const PgReferencePoint query = queries[p];
            auto exact = partition_take(
                conditioning.evaluate_reference(std::span{&query, 1}, remaining, plan));
            math.charge(exact.call_info.scalar_operations);
            math.visit(exact.call_info.partition_visits);
            project_mass_point(math, exact, request, high_row.first(width),
                               two_word ? low_row.first(width) : std::span<double>{});
            const auto& geo = exact.points.values().front();
            points.values()[p] = {geo.identity, geo.cell, geo.partition_record, geo.layer};
            high_row = high_row.subspan(width);
            if (two_word)
                low_row = low_row.subspan(width);
        }
        math.cancel();
        return PgVolumeMassReference{
            std::move(conditioning), rank,          request, math.info(), std::move(points),
            std::move(high),         std::move(low)};
    });
}

} // namespace uvlas::structure
