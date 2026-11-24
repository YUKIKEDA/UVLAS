#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <uvlas/structure/pg_prepared_patch.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;

} // namespace

TEST_CASE("PG prepared patch reconstructs fixed physical position and translation",
          "[SH-D02][SH-D03][N4-03][N4-06][N4-08][verification]") {
    Fixture f;
    for (auto& h : f.heights)
        h = f.r("2");
    const auto length = f.r("2");
    std::array<ExactRational, 9> rotation;
    rotation[2] = f.r("1");
    rotation[3] = f.r("1");
    rotation[7] = f.r("1");
    const std::array origin{f.r("1"), f.r("-2"), f.r("3")};
    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan setup{profile};
    auto basis = make_basis(f, {rotation, origin, length}, scratch, setup);
    const auto initial_geometry = take(basis.reference_identity());
    const auto caps = query_limits();
    auto gram = take(PgFlatGram::prepare(std::move(basis), 256, {caps, steps}, scratch, setup));
    auto source = take(
        std::move(gram).generate_and_certify({320, 246}, {profile, steps, {}}, scratch, setup));
    const auto prepare = [&](const PgReferenceLimits& query, core::ResourcePlan& budget) {
        return PgPreparedPatch::prepare(std::move(source), query, scratch, budget);
    };
    core::ResourcePlan none{0}, measured{profile};
    auto cancelled = caps;
    cancelled.cancellation = requested_stop_token();
    check_failure(prepare(cancelled, none), core::FailureCode::Cancelled);
    check_failure(prepare(caps, none), core::FailureCode::ResourceExceeded);
    CHECK(none.allocated_bytes() == 0);
    auto short_visits = caps;
    short_visits.maximum_partition_visits = 1;
    check_failure(prepare(short_visits, measured), core::FailureCode::ResourceExceeded);
    REQUIRE(source.info());
    const auto consumed = [&] { CHECK_FALSE(source.info()); };
    auto patch = [&] {
        core::ResourcePlan local{profile};
        auto result = take(prepare(caps, local));
        CHECK(take(result.info()).preparation.allocated_owned_payload_bytes ==
              local.allocated_bytes());
        return result;
    }();
    consumed();
    const auto recorded = take(patch.info());
    CHECK(recorded.geometry == initial_geometry);
    CHECK_FALSE(patch.fixed_rule());
    CHECK(recorded.basis.rank == patch.reference_coordinates().size());
    CHECK(recorded.basis.rank == patch.translation().size());
    const auto* selected = take(patch.conditioning().reference_basis());
    const auto& partition = take(selected->source()).partition;
    const auto vertices = take(partition.vertices());
    const auto triangles = take(partition.triangles());
    const auto volumes = take(partition.volumes());
    const auto layers = take(partition.layers());
    const auto faces = take(partition.horizontal_faces());
    const auto u = f.r("1", "4"), v = f.r("1", "3"), depth = f.r("2", "5");
    const std::array<PgReferencePoint, 2> queries{PgReferenceVolumePoint{0, u, v, depth},
                                                  PgReferenceHorizontalPoint{0, u, v}};
    const auto batch = take(patch.evaluate_reference(queries, caps, setup));
    REQUIRE(batch.points.size() == queries.size());
    const auto cartesian = [&](auto member, const PgPartitionTriangle& tri) -> Q {
        const auto& a = vertices[static_cast<std::size_t>(tri.vertices[0])].*member;
        const auto& b = vertices[static_cast<std::size_t>(tri.vertices[1])].*member;
        const auto& c = vertices[static_cast<std::size_t>(tri.vertices[2])].*member;
        return (1 - q(u) - q(v)) * q(a) + q(u) * q(b) + q(v) * q(c);
    };
    const auto flat = [](const PgReferenceScalar& x) -> Q { return q(x.rational) + q(x.radical); };
    for (std::size_t p = 0; p < queries.size(); ++p) {
        CAPTURE(p);
        const auto triangle = p == 0 ? volumes.front().triangle : faces.front().triangle;
        const auto ordinal = p == 0 ? volumes.front().layer : faces.front().layer;
        const auto& tri = triangles[static_cast<std::size_t>(triangle)];
        const auto& layer = layers[static_cast<std::size_t>(ordinal)];
        const auto s = cartesian(&PgPartitionPoint::s, tri);
        const auto t = cartesian(&PgPartitionPoint::t, tri);
        const Q zeta = p == 0 ? q(layer.lower) + q(depth) * (q(layer.upper) - q(layer.lower))
                              : q(faces.front().side == PgPartitionFaceSide::Lower ? layer.lower
                                                                                   : layer.upper);
        const std::array<Q, 3> expected{1 + zeta, -2 + 2 * s, 3 + 2 * t};
        const auto& value = batch.points.values()[p];
        CHECK(value.identity == initial_geometry);
        REQUIRE(q(value.g) == 1);
        std::array<Q, 3> position{}, shift_gradient{};
        std::array<std::array<Q, 3>, 3> derivative{};
        Q shift{};
        const auto row = batch.columns.values().subspan(p * patch.translation().size(),
                                                        patch.translation().size());
        for (std::size_t i = 0; i < row.size(); ++i) {
            const Q n = flat(row[i].value);
            const auto& xyz = patch.reference_coordinates()[i];
            shift += n * q(patch.translation()[i]);
            for (std::size_t a = 0; a < expected.size(); ++a) {
                position[a] += n * q(xyz[a]);
                const Q g = flat(row[i].physical_gradient_per_metre[a]);
                shift_gradient[a] += g * q(patch.translation()[i]);
                for (std::size_t j = 0; j < expected.size(); ++j)
                    derivative[j][a] += g * q(xyz[j]);
            }
        }
        CHECK(position == expected);
        CHECK(shift == 1);
        for (std::size_t a = 0; a < expected.size(); ++a) {
            CHECK(flat(value.position_metres[a]) == expected[a]);
            CHECK(shift_gradient[a] == 0);
            for (std::size_t j = 0; j < expected.size(); ++j)
                CHECK(derivative[j][a] == (j == a ? 1 : 0));
        }
    }

    const auto interior = std::ranges::find_if(
        faces, [&](const auto& face) { return face.side != faces.front().side; });
    REQUIRE(interior != faces.end());
    const auto opposite_face = static_cast<core::Index>(interior - faces.begin());
    const std::array<PgReferencePoint, 4> mixed{queries[0], queries[1],
                                                PgReferenceHorizontalPoint{opposite_face, u, v},
                                                PgReferenceLateralPoint{0, u, depth}};
    const auto exact = take(patch.evaluate_reference(mixed, caps, setup));
    const auto first_side = faces.front().side == PgPartitionFaceSide::Lower
                                ? PgProjectedReferenceDomain::Lower
                                : PgProjectedReferenceDomain::Upper;
    const auto opposite_side = first_side == PgProjectedReferenceDomain::Lower
                                   ? PgProjectedReferenceDomain::Upper
                                   : PgProjectedReferenceDomain::Lower;
    const std::array domains{PgProjectedReferenceDomain::Volume, first_side, opposite_side,
                             PgProjectedReferenceDomain::Lateral};
    const auto stride = 4 * (static_cast<std::size_t>(exact.column_count) + 1);
    for (const auto precision : {core::Precision::Binary64, core::Precision::TwoWord}) {
        const auto paired = precision == core::Precision::TwoWord;
        const double requested = paired ? 0x1p-72 : 0x1p-32;
        const Q tolerance = independent_double(requested);
        const PgPreparedReferenceProjection projection{precision, requested, requested,
                                                       requested, requested, requested};
        const auto result = take(patch.project_reference(mixed, projection, caps, setup));
        CHECK(result.source() == recorded.identity);
        CHECK(result.column_count() == exact.column_count);
        REQUIRE(result.points().size() == mixed.size());
        REQUIRE(result.high_words().size() == stride * mixed.size());
        REQUIRE(result.low_words().size() == (paired ? stride * mixed.size() : 0));
        for (std::size_t p = 0; p < mixed.size(); ++p) {
            CAPTURE(p);
            const auto& record = result.points()[p];
            const auto& point = exact.points.values()[p];
            CHECK(record.domain == domains[p]);
            CHECK(record.layer == point.layer);
            CHECK(record.partition_record == point.partition_record);
            std::size_t slot = p * stride;
            const auto read = [&] {
                Q value = independent_double(result.high_words()[slot]);
                if (paired)
                    value += independent_double(result.low_words()[slot]);
                ++slot;
                return value;
            };
            const auto check_scalar = [&](const PgReferenceScalar& component) {
                CHECK(abs(read() - flat(component)) <= tolerance);
            };
            REQUIRE(q(point.g) == 1);
            for (const auto& position : point.position_metres)
                check_scalar(position);
            if (const auto* volume = std::get_if<PgReferenceVolumeMeasure>(&point.measure)) {
                check_scalar(volume->jacobian);
            } else {
                const auto* side = std::get_if<PgReferenceHorizontalMeasure>(&point.measure);
                const auto& area =
                    side ? side->area : std::get<PgReferenceLateralMeasure>(point.measure).area;
                const Q kept = read();
                CHECK(kept > tolerance);
                CHECK((kept - tolerance) * (kept - tolerance) <= flat(area.squared_norm));
                CHECK((kept + tolerance) * (kept + tolerance) >= flat(area.squared_norm));
            }
            const auto columns =
                exact.columns.values().subspan(p * result.column_count(), result.column_count());
            for (const auto& column : columns) {
                check_scalar(column.value);
                for (const auto& gradient : column.physical_gradient_per_metre)
                    check_scalar(gradient);
            }
            CHECK(slot == (p + 1) * stride);
        }
        check_failure(patch.project_reference(mixed, projection, caps, none),
                      core::FailureCode::ResourceExceeded);
        auto strict = projection;
        strict.maximum_position_error_metres = 0;
        check_failure(patch.project_reference(mixed, strict, caps, setup),
                      core::FailureCode::PrecisionUnresolved);
        CHECK(take(patch.info()).identity == recorded.identity);
    }
}
