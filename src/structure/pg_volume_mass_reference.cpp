#include "pg_reference_domain.hpp"
#include "pg_reference_projection.hpp"

#include <uvlas/structure/pg_boundary_tangents.hpp>
#include <uvlas/structure/pg_material_tangents.hpp>
#include <uvlas/structure/pg_prepared_volume.hpp>

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

[[nodiscard]] PgReferencePoint reference_query(const PgOwnedDiscreteVolumePoint& query) {
    return PgReferenceVolumePoint{query.volume, query.u, query.v, query.depth};
}

[[nodiscard]] PgReferencePoint reference_query(const PgReferencePoint& query) {
    return query;
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
             std::uint64_t prefix, std::uint64_t directions, const FillRow& fill_row,
             std::uint64_t retained_rationals_per_point = 0) {
    if (queries.empty())
        return {};
    auto& plan = math.plan;
    const auto rows = partition_index(queries.size());
    const auto stride = partition_add(partition_mul(rank, directions), prefix);
    const auto word_count = partition_index(partition_mul(queries.size(), stride));
    const auto point_bytes =
        partition_take(core::checked_byte_count(rows, sizeof(Record), partition_operation));
    const auto word_bytes =
        partition_take(core::checked_byte_count(word_count, sizeof(double), partition_operation));
    if (partition_add(point_bytes, partition_mul(two_word ? 2 : 1, word_bytes)) >
        plan.remaining_bytes())
        partition_fail(core::FailureCode::ResourceExceeded);
    // Only one exact reference row remains alive while producing the rounded row.
    // Multicomponent projections also need their Cartesian and area temporaries.
    const auto extra = directions > 1 ? 32U : 0U;
    const auto retained_rationals = partition_mul(queries.size(), retained_rationals_per_point);
    math.admit(
        partition_add(partition_add(reference_live_rationals + 74 + extra, retained_rationals),
                      partition_mul(8, rank)));
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
        math.reserve_nested_rationals(remaining.algebra,
                                      partition_mul(p, retained_rationals_per_point));
        const auto query = reference_query(queries[p]);
        const auto exact =
            partition_take(conditioning.evaluate_reference(std::span{&query, 1}, remaining, plan));
        math.charge(exact.call_info.scalar_operations);
        math.visit(exact.call_info.partition_visits);
        points.values()[p] = fill_row(p, exact, high_row.first(width),
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

[[nodiscard]] PgBoundaryPositionRecord boundary_record(const PgReferenceValue& point) {
    return {point.identity, point.cell, point.partition_record, point.layer, boundary_side(point)};
}

[[nodiscard]] PgReferenceVector planar_direction(ReferenceRadicals& field,
                                                 const PgReferenceValue& point,
                                                 const PgPartitionPoint& from,
                                                 const PgPartitionPoint& to) {
    const auto ds = field.math.sub(to.s, from.s), dt = field.math.sub(to.t, from.t);
    const auto& j = point.jacobian_columns_metres;
    PgReferenceVector direction;
    for (std::size_t axis = 0; axis < direction.size(); ++axis)
        direction[axis] = field.add(field.scale(j[0][axis], ds), field.scale(j[1][axis], dt));
    return direction;
}

[[nodiscard]] std::array<PgReferenceVector, 2>
face_directions(ReferenceRadicals& field, const PgReferenceValue& point,
                const PgIntegrationPartition& partition, PgBoundaryPositionSide side) {
    auto& math = field.math;
    math.visit(5); // Face plus triangle and three vertices, or edge, layer and two vertices.
    const auto vertices = partition_take(partition.vertices());
    if (side != PgBoundaryPositionSide::Lateral) {
        const auto& face =
            reference_at(partition_take(partition.horizontal_faces()), point.partition_record);
        const auto& triangle = reference_at(partition_take(partition.triangles()), face.triangle);
        const auto& a = reference_at(vertices, triangle.vertices[0]);
        const auto& b = reference_at(vertices, triangle.vertices[1]);
        const auto& c = reference_at(vertices, triangle.vertices[2]);
        return {planar_direction(field, point, a, b), planar_direction(field, point, a, c)};
    }
    const auto& face =
        reference_at(partition_take(partition.lateral_faces()), point.partition_record);
    const auto& edge = reference_at(partition_take(partition.edges()), face.edge);
    const auto& layer = reference_at(partition_take(partition.layers()), face.layer);
    const auto first = face.canonical_direction ? 0U : 1U;
    const auto& a = reference_at(vertices, edge.vertices[first]);
    const auto& b = reference_at(vertices, edge.vertices[1 - first]);
    const auto width = math.sub(layer.upper, layer.lower);
    PgReferenceVector depth;
    for (std::size_t axis = 0; axis < depth.size(); ++axis)
        depth[axis] = field.scale(point.jacobian_columns_metres[2][axis], width);
    return {planar_direction(field, point, a, b), std::move(depth)};
}

void project_tangents(ReferenceArithmetic& math, const PgConditionedReferenceBatch& batch,
                      const PgIntegrationPartition& partition, PgBoundaryPositionSide side,
                      const PgBoundaryTangentProjection& request, std::span<double> high,
                      std::span<double> low) {
    const auto& point = batch.points.values().front();
    ReferenceRadicals field{math, point.g};
    const auto directions = face_directions(field, point, partition, side);
    const auto tolerance = reference_binary64(math, request.maximum_derivative_error);
    const auto columns = batch.columns.values();
    std::size_t k = 0;
    for (const auto& direction : directions)
        for (const auto& column : columns) {
            math.visit(1);
            const auto derivative = field.dot(column.physical_gradient_per_metre, direction);
            const auto words =
                pg_project_scalar(math, derivative, point.g, tolerance, request.representation);
            high[k] = words.high;
            if (!low.empty())
                low[k] = words.low;
            ++k;
        }
}

[[nodiscard]] PgProjectedReferencePoint record(const PgReferenceValue& point) {
    auto domain = PgProjectedReferenceDomain::Volume;
    if (!std::holds_alternative<PgReferenceVolumeMeasure>(point.measure))
        switch (boundary_side(point)) {
        case PgBoundaryPositionSide::Lower:
            domain = PgProjectedReferenceDomain::Lower;
            break;
        case PgBoundaryPositionSide::Upper:
            domain = PgProjectedReferenceDomain::Upper;
            break;
        case PgBoundaryPositionSide::Lateral:
            domain = PgProjectedReferenceDomain::Lateral;
            break;
        }
    return {point.cell, point.partition_record, point.layer, domain};
}

void project_material_tangents(ReferenceArithmetic& math, const PgConditionedReferenceBatch& batch,
                               const PgMaterialDirections& motion,
                               const PgMaterialTangentProjection& request, std::span<double> high,
                               std::span<double> low) {
    const auto& point = batch.points.values().front();
    ReferenceRadicals field{math, point.g};
    const auto tolerance = reference_binary64(math, request.maximum_derivative_error);
    std::size_t slot = 0;
    for (const auto& change : motion) {
        PgReferenceVector world_direction;
        for (std::size_t axis = 0; axis < 3; ++axis)
            for (std::size_t chart = 0; chart < 3; ++chart)
                world_direction[axis] = field.add(
                    world_direction[axis],
                    field.scale(point.jacobian_columns_metres[chart][axis], change[chart]));
        for (const auto& column : batch.columns.values()) {
            math.visit(1);
            const auto d = field.dot(column.physical_gradient_per_metre, world_direction);
            const auto words =
                pg_project_scalar(math, d, point.g, tolerance, request.representation);
            high[slot] = words.high;
            if (!low.empty())
                low[slot] = words.low;
            ++slot;
        }
    }
}

[[nodiscard]] PgProjectedReferencePoint
project_reference_row(ReferenceArithmetic& math, const PgConditionedReferenceBatch& batch,
                      const PgPreparedReferenceProjection& request, std::span<double> high,
                      std::span<double> low) {
    const auto& point = batch.points.values().front();
    const auto precision = request.representation;
    const auto position_error = reference_binary64(math, request.maximum_position_error_metres);
    const auto shape_error = reference_binary64(math, request.maximum_shape_error);
    const auto gradient_error = reference_binary64(math, request.maximum_gradient_error_per_metre);
    std::size_t slot = 0;
    const auto write = [&](PgProjectedWord words) {
        high[slot] = words.high;
        if (!low.empty())
            low[slot] = words.low;
        ++slot;
    };
    const auto scalar = [&](const PgReferenceScalar& value, const ExactRational& error) {
        return pg_project_scalar(math, value, point.g, error, precision);
    };
    for (const auto& position : point.position_metres)
        write(scalar(position, position_error));

    PgProjectedReferenceDomain domain = PgProjectedReferenceDomain::Volume;
    PgProjectedWord measure;
    if (const auto* volume = std::get_if<PgReferenceVolumeMeasure>(&point.measure)) {
        const auto error = reference_binary64(math, request.maximum_volume_error_metres_cubed);
        measure = scalar(volume->jacobian, error);
    } else {
        const auto* side = std::get_if<PgReferenceHorizontalMeasure>(&point.measure);
        const auto& area =
            side ? side->area : std::get<PgReferenceLateralMeasure>(point.measure).area;
        domain =
            side ? (side->side == PgPartitionFaceSide::Lower ? PgProjectedReferenceDomain::Lower
                                                             : PgProjectedReferenceDomain::Upper)
                 : PgProjectedReferenceDomain::Lateral;
        const auto error = reference_binary64(math, request.maximum_area_error_metres_squared);
        measure = pg_project_positive_root(math, area.squared_norm, point.g, error, precision);
    }
    const auto measure_bits = std::bit_cast<std::uint64_t>(measure.high);
    if (measure_bits == 0 || (measure_bits & pg_binary64_sign) != 0)
        partition_fail(core::FailureCode::PrecisionUnresolved);
    write(measure);
    for (const auto& column : batch.columns.values()) {
        math.visit(1);
        write(scalar(column.value, shape_error));
        for (const auto& gradient : column.physical_gradient_per_metre)
            write(scalar(gradient, gradient_error));
    }
    return {point.cell, point.partition_record, point.layer, domain};
}

} // namespace

PgProjectedReferenceBatch::PgProjectedReferenceBatch(
    PgPreparedPatchIdentity source, std::size_t columns, PgPreparedReferenceProjection projection,
    PgReferenceCallInfo call_info, core::PlannedBuffer<PgProjectedReferencePoint> points,
    core::PlannedBuffer<double> high, core::PlannedBuffer<double> low,
    std::optional<std::size_t> first_rule_ordinal) noexcept
    : source_(source),
      columns_(columns),
      projection_(projection),
      call_info_(call_info),
      points_(std::move(points)),
      high_(std::move(high)),
      low_(std::move(low)),
      first_rule_ordinal_(first_rule_ordinal) {}

template <class Query>
core::Result<PgProjectedReferenceBatch>
PgPreparedPatch::project_points(std::span<const Query> queries,
                                const PgPreparedReferenceProjection& projection,
                                const PgReferenceLimits& limits, core::ResourcePlan& plan,
                                std::optional<std::size_t> first_rule_ordinal) const noexcept {
    return partition_entry([&]() -> core::Result<PgProjectedReferenceBatch> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto source = partition_take(info());
        const auto two_word = two_word_projection(projection.representation);
        for (const auto tolerance :
             {projection.maximum_position_error_metres, projection.maximum_shape_error,
              projection.maximum_gradient_error_per_metre,
              projection.maximum_volume_error_metres_cubed,
              projection.maximum_area_error_metres_squared})
            validate_tolerance(tolerance);
        auto rows = project_rows<PgProjectedReferencePoint>(
            math, conditioning(), source.basis.rank, queries, two_word, 4, 4,
            [&](std::size_t, const PgConditionedReferenceBatch& batch, std::span<double> high,
                std::span<double> low) {
                return project_reference_row(math, batch, projection, high, low);
            });
        return PgProjectedReferenceBatch{source.identity,
                                         static_cast<std::size_t>(source.basis.rank),
                                         projection,
                                         math.info(),
                                         std::move(rows.points),
                                         std::move(rows.high),
                                         std::move(rows.low),
                                         first_rule_ordinal};
    });
}

core::Result<PgProjectedReferenceBatch> PgPreparedPatch::project_reference(
    std::span<const PgReferencePoint> queries, const PgPreparedReferenceProjection& projection,
    const PgReferenceLimits& limits, core::ResourcePlan& plan) const noexcept {
    return project_points(queries, projection, limits, plan);
}

core::Result<PgProjectedReferenceBatch> PgPreparedPatch::project_fixed_reference(
    std::size_t first, std::size_t count, const PgPreparedReferenceProjection& projection,
    const PgReferenceLimits& limits, core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgProjectedReferenceBatch> {
        if (limits.cancellation.stop_requested())
            partition_fail(core::FailureCode::Cancelled);
        const auto rule = fixed_rule();
        if (!rule)
            partition_fail(core::FailureCode::Unsupported);
        if (count == 0 || first > rule->size() || count > rule->size() - first)
            partition_fail(core::FailureCode::InvalidInput);
        return project_points(rule->subspan(first, count), projection, limits, plan, first);
    });
}

core::Result<PgProjectedReferenceBatch>
PgPreparedPatch::project_boundary_reference(std::span<const PgBoundaryPositionPoint> queries,
                                            const PgPreparedReferenceProjection& projection,
                                            const PgReferenceLimits& limits,
                                            core::ResourcePlan& plan) const noexcept {
    return project_points(queries, projection, limits, plan);
}

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
            math, conditioning, accepted.rank, queries, two_word, 1, 1,
            [&](std::size_t, const PgConditionedReferenceBatch& point, std::span<double> high,
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
            math, conditioning, accepted.rank, queries, two_word, 0, 1,
            [&](std::size_t, const PgConditionedReferenceBatch& point, std::span<double> high,
                std::span<double> low) {
                project_shapes(math, point, request.maximum_shape_error, request.representation,
                               high, low);
                return boundary_record(point.points.values().front());
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

PgBoundaryTangentReference::PgBoundaryTangentReference(
    std::size_t columns, PgBoundaryTangentProjection projection, PgReferenceCallInfo info,
    core::PlannedBuffer<PgBoundaryPositionRecord> points, core::PlannedBuffer<double> high,
    core::PlannedBuffer<double> low) noexcept
    : column_count_(columns),
      projection_(projection),
      call_info_(info),
      points_(std::move(points)),
      high_(std::move(high)),
      low_(std::move(low)) {}

PgBoundaryTangentReference::PgBoundaryTangentReference(PgBoundaryTangentReference&& other) noexcept
    : PgBoundaryTangentReference(std::exchange(other.column_count_, 0), other.projection_,
                                 other.call_info_, std::move(other.points_), std::move(other.high_),
                                 std::move(other.low_)) {}

PgBoundaryTangentReference&
PgBoundaryTangentReference::operator=(PgBoundaryTangentReference&& other) noexcept {
    if (this != &other) {
        column_count_ = std::exchange(other.column_count_, 0);
        projection_ = other.projection_;
        call_info_ = other.call_info_;
        points_ = std::move(other.points_);
        high_ = std::move(other.high_);
        low_ = std::move(other.low_);
    }
    return *this;
}

core::Result<PgBoundaryTangentReference> PgBoundaryTangentReference::prepare(
    const PgConditioningCertificate& conditioning, std::span<const PgBoundaryPositionPoint> queries,
    const PgBoundaryTangentProjection& request, const PgReferenceLimits& limits,
    core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgBoundaryTangentReference> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto accepted = partition_take(conditioning.info());
        const auto two_word = two_word_projection(request.representation);
        validate_tolerance(request.maximum_derivative_error);
        const auto* basis = partition_take(conditioning.reference_basis());
        const auto& partition = partition_take(basis->source()).partition;
        auto result = project_rows<PgBoundaryPositionRecord>(
            math, conditioning, accepted.rank, queries, two_word, 0, 2,
            [&](std::size_t, const PgConditionedReferenceBatch& batch, std::span<double> high,
                std::span<double> low) {
                const auto record = boundary_record(batch.points.values().front());
                project_tangents(math, batch, partition, record.side, request, high, low);
                return record;
            });
        return PgBoundaryTangentReference{static_cast<std::size_t>(accepted.rank),
                                          request,
                                          math.info(),
                                          std::move(result.points),
                                          std::move(result.high),
                                          std::move(result.low)};
    });
}

PgMaterialTangentReference::PgMaterialTangentReference(
    PgPreparedPatchIdentity source, std::size_t columns, PgMaterialTangentProjection projection,
    PgReferenceCallInfo call_info, core::PlannedBuffer<PgProjectedReferencePoint> points,
    core::PlannedBuffer<double> high, core::PlannedBuffer<double> low) noexcept
    : source_(source),
      column_count_(columns),
      projection_(projection),
      call_info_(call_info),
      points_(std::move(points)),
      high_(std::move(high)),
      low_(std::move(low)) {}

PgMaterialTangentReference::PgMaterialTangentReference(PgMaterialTangentReference&& other) noexcept
    : PgMaterialTangentReference(other.source_, std::exchange(other.column_count_, 0),
                                 other.projection_, other.call_info_, std::move(other.points_),
                                 std::move(other.high_), std::move(other.low_)) {}

PgMaterialTangentReference&
PgMaterialTangentReference::operator=(PgMaterialTangentReference&& other) noexcept {
    if (this != &other) {
        source_ = other.source_;
        column_count_ = std::exchange(other.column_count_, 0);
        projection_ = other.projection_;
        call_info_ = other.call_info_;
        points_ = std::move(other.points_);
        high_ = std::move(other.high_);
        low_ = std::move(other.low_);
    }
    return *this;
}

core::Result<PgMaterialTangentReference> PgMaterialTangentReference::prepare(
    const PgPreparedPatch& patch, std::span<const PgReferencePoint> points,
    std::span<const PgMaterialDirections> directions, const PgMaterialTangentProjection& projection,
    const PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgMaterialTangentReference> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        if (points.size() != directions.size())
            partition_fail(core::FailureCode::InvalidInput);
        const auto source = partition_take(patch.info());
        const auto paired = two_word_projection(projection.representation);
        validate_tolerance(projection.maximum_derivative_error);
        auto rows = project_rows<PgProjectedReferencePoint>(
            math, patch.conditioning(), source.basis.rank, points, paired, 0, 2,
            [&](std::size_t i, const PgConditionedReferenceBatch& batch, std::span<double> high,
                std::span<double> low) {
                project_material_tangents(math, batch, directions[i], projection, high, low);
                return record(batch.points.values().front());
            });
        return PgMaterialTangentReference{source.identity,
                                          static_cast<std::size_t>(source.basis.rank),
                                          projection,
                                          math.info(),
                                          std::move(rows.points),
                                          std::move(rows.high),
                                          std::move(rows.low)};
    });
}

PgPreparedVolumeReference::PgPreparedVolumeReference(
    PgPreparedPatchIdentity source, std::size_t columns, PgVolumeMassProjection projection,
    PgReferenceCallInfo call_info, core::PlannedBuffer<PgPreparedVolumeSample> points,
    core::PlannedBuffer<double> high, core::PlannedBuffer<double> low) noexcept
    : source_(source),
      columns_(columns),
      projection_(projection),
      call_info_(call_info),
      points_(std::move(points)),
      high_(std::move(high)),
      low_(std::move(low)) {}

core::Result<PgPreparedVolumeReference>
PgPreparedVolumeReference::prepare(const PgPreparedPatch& patch, std::size_t first,
                                   std::size_t count, const PgVolumeMassProjection& projection,
                                   const PgReferenceLimits& limits,
                                   core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgPreparedVolumeReference> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto source = partition_take(patch.info());
        const auto rule = patch.fixed_rule();
        if (!rule)
            partition_fail(core::FailureCode::Unsupported);
        if (count == 0 || first > rule->size() || count > rule->size() - first)
            partition_fail(core::FailureCode::InvalidInput);
        const auto paired = two_word_projection(projection.representation);
        validate_tolerance(projection.maximum_shape_error);
        validate_tolerance(projection.maximum_volume_jacobian_error_metres_cubed);
        const auto queries = rule->subspan(first, count);
        auto ordinal = first;
        auto result = project_rows<PgPreparedVolumeSample>(
            math, patch.conditioning(), source.basis.rank, queries, paired, 1, 1,
            [&](std::size_t, const PgConditionedReferenceBatch& exact, std::span<double> high,
                std::span<double> low) {
                project_mass_point(math, exact, projection, high, low);
                project_shapes(math, exact, projection.maximum_shape_error,
                               projection.representation, high.subspan(1),
                               paired ? low.subspan(1) : low);
                const auto& reference = exact.points.values().front();
                const auto current = ordinal++;
                return PgPreparedVolumeSample{{reference.identity, reference.cell,
                                               reference.partition_record, reference.layer},
                                              current,
                                              math.copy((*rule)[current].weight)};
            },
            1); // The owned output retains one exact weight per completed point.
        return PgPreparedVolumeReference{source.identity,
                                         static_cast<std::size_t>(source.basis.rank),
                                         projection,
                                         math.info(),
                                         std::move(result.points),
                                         std::move(result.high),
                                         std::move(result.low)};
    });
}

} // namespace uvlas::structure
