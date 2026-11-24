#include "pg_reference_projection.hpp"

#include <uvlas/structure/pg_boundary_position.hpp>

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
}

void project_shapes(ReferenceArithmetic& math, const PgConditionedReferenceBatch& point,
                    double maximum_error, core::Precision precision, std::span<double> high,
                    std::span<double> low) {
    const auto tolerance = reference_binary64(math, maximum_error);
    const auto& reference = point.points.values().front();
    const auto columns = point.columns.values();
    for (std::size_t j = 0; j < columns.size(); ++j) {
        math.visit(1);
        const auto words =
            pg_project_scalar(math, columns[j].value, reference.g, tolerance, precision);
        high[j] = words.high;
        if (!low.empty())
            low[j] = words.low;
    }
}

[[nodiscard]] bool two_word_projection(core::Precision precision) {
    if (precision != core::Precision::TwoWord && precision != core::Precision::Binary64)
        partition_fail(core::FailureCode::PrecisionUnsupported);
    return precision == core::Precision::TwoWord;
}

[[nodiscard]] PgReferencePoint reference_query(const PgReferenceVolumePoint& query) {
    return query;
}

[[nodiscard]] PgReferencePoint reference_query(const PgBoundaryPositionPoint& query) {
    return std::visit([](const auto& point) -> PgReferencePoint { return point; }, query);
}

template <class Record>
struct ProjectedRows {
    core::PlannedBuffer<Record> points;
    core::PlannedBuffer<double> high, low;
};

template <class Record, class Query, class FillRow>
[[nodiscard]] ProjectedRows<Record>
project_rows(ReferenceArithmetic& math, const PgConditioningCertificate& conditioning,
             std::uint64_t rank, std::span<const Query> queries, bool two_word,
             std::uint64_t prefix, const FillRow& fill_row) {
    if (queries.empty())
        return {};
    auto& plan = math.plan;
    const auto rows = partition_index(queries.size());
    const auto stride = partition_add(rank, prefix);
    const auto word_count = partition_index(partition_mul(queries.size(), stride));
    const auto point_bytes =
        partition_take(core::checked_byte_count(rows, sizeof(Record), partition_operation));
    const auto word_bytes =
        partition_take(core::checked_byte_count(word_count, sizeof(double), partition_operation));
    if (partition_add(point_bytes, partition_mul(two_word ? 2 : 1, word_bytes)) >
        plan.remaining_bytes())
        partition_fail(core::FailureCode::ResourceExceeded);
    // Only one exact reference row remains alive while producing the rounded row.
    math.admit(partition_add(reference_live_rationals + 74, partition_mul(8, rank)));
    const auto words = static_cast<std::uint64_t>(word_count);
    math.visit(partition_add(queries.size(), partition_mul(two_word ? 2 : 1, words)));
    auto points = partition_take(plan.allocate<Record>(rows, partition_operation));
    auto high = partition_take(plan.allocate<double>(word_count, partition_operation));
    auto low =
        partition_take(plan.allocate<double>(two_word ? word_count : 0, partition_operation));
    auto high_row = high.values();
    auto low_row = low.values();
    const auto width = static_cast<std::size_t>(stride);
    for (std::size_t p = 0; p < queries.size(); ++p) {
        math.visit(1);
        auto remaining = math.limits;
        remaining.maximum_scalar_operations -= math.scalar_operations;
        remaining.maximum_partition_visits -= math.partition_visits;
        const auto query = reference_query(queries[p]);
        const auto exact =
            partition_take(conditioning.evaluate_reference(std::span{&query, 1}, remaining, plan));
        math.charge(exact.call_info.scalar_operations);
        math.visit(exact.call_info.partition_visits);
        points.values()[p] = fill_row(exact, high_row.first(width),
                                      two_word ? low_row.first(width) : std::span<double>{});
        high_row = high_row.subspan(width);
        if (two_word)
            low_row = low_row.subspan(width);
    }
    math.cancel();
    return {std::move(points), std::move(high), std::move(low)};
}

[[nodiscard]] PgBoundaryPositionSide boundary_side(const PgReferenceValue& point) {
    if (const auto* horizontal = std::get_if<PgReferenceHorizontalMeasure>(&point.measure))
        return horizontal->side == PgPartitionFaceSide::Upper ? PgBoundaryPositionSide::Upper
                                                              : PgBoundaryPositionSide::Lower;
    if (std::holds_alternative<PgReferenceLateralMeasure>(point.measure))
        return PgBoundaryPositionSide::Lateral;
    partition_fail(core::FailureCode::InvalidInput);
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
        const auto two_word = two_word_projection(request.representation);
        validate_tolerance(request.maximum_shape_error);
        validate_tolerance(request.maximum_volume_jacobian_error_metres_cubed);
        auto result = project_rows<PgVolumeMassPoint>(
            math, conditioning, accepted.rank, queries, two_word, 1,
            [&](const PgConditionedReferenceBatch& point, std::span<double> high,
                std::span<double> low) {
                project_mass_point(math, point, request, high, low);
                project_shapes(math, point, request.maximum_shape_error, request.representation,
                               high.subspan(1), two_word ? low.subspan(1) : low);
                const auto& geo = point.points.values().front();
                return PgVolumeMassPoint{geo.identity, geo.cell, geo.partition_record, geo.layer};
            });
        return PgVolumeMassReference{std::move(conditioning),
                                     static_cast<std::size_t>(accepted.rank),
                                     request,
                                     math.info(),
                                     std::move(result.points),
                                     std::move(result.high),
                                     std::move(result.low)};
    });
}

PgBoundaryPositionReference::PgBoundaryPositionReference(
    PgConditioningCertificate&& conditioning, std::size_t columns,
    PgBoundaryPositionProjection projection, PgReferenceCallInfo info,
    core::PlannedBuffer<PgBoundaryPositionRecord> points, core::PlannedBuffer<double> high,
    core::PlannedBuffer<double> low) noexcept
    : conditioning_(std::move(conditioning)),
      column_count_(columns),
      projection_(projection),
      call_info_(info),
      points_(std::move(points)),
      high_(std::move(high)),
      low_(std::move(low)) {}

core::Result<PgBoundaryPositionReference> PgBoundaryPositionReference::prepare(
    PgConditioningCertificate&& conditioning, std::span<const PgBoundaryPositionPoint> queries,
    const PgBoundaryPositionProjection& request, const PgReferenceLimits& limits,
    core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgBoundaryPositionReference> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto accepted = partition_take(conditioning.info());
        const auto two_word = two_word_projection(request.representation);
        validate_tolerance(request.maximum_shape_error);
        auto result = project_rows<PgBoundaryPositionRecord>(
            math, conditioning, accepted.rank, queries, two_word, 0,
            [&](const PgConditionedReferenceBatch& point, std::span<double> high,
                std::span<double> low) {
                project_shapes(math, point, request.maximum_shape_error, request.representation,
                               high, low);
                const auto& geo = point.points.values().front();
                return PgBoundaryPositionRecord{geo.identity, geo.cell, geo.partition_record,
                                                geo.layer, boundary_side(geo)};
            });
        return PgBoundaryPositionReference{std::move(conditioning),
                                           static_cast<std::size_t>(accepted.rank),
                                           request,
                                           math.info(),
                                           std::move(result.points),
                                           std::move(result.high),
                                           std::move(result.low)};
    });
}

} // namespace uvlas::structure
