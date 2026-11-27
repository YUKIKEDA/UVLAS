#include "local_geometry_packet.hpp"

#include "global_geometry.hpp"
#include "local_geometry_bounds_private.hpp"
#include "local_geometry_input_private.hpp"
#include "shared_p1_reuse_private.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <type_traits>

namespace uvlas::aero::detail {
static_assert(std::is_nothrow_copy_constructible_v<LocalAmbientDescriptor>);
static_assert(std::is_nothrow_copy_constructible_v<LocalQualification>);
static_assert(std::is_trivially_copyable_v<QueryRequest>);
namespace {
constexpr auto operation = core::Operation::EvaluateAero;
constexpr std::uint64_t enclosure_steps = 12 * (257 + 2);
OwnedGeometryFailure packet_failure(core::Failure failure, const GeometryContext& context,
                                    std::optional<GeometryProblem> problem = {}) noexcept {
    auto& d = failure.details;
    const auto& origin = context.diagnostics;
    if (!d.model)
        d.model = origin.model ? origin.model : std::optional{context.model};
    if (!d.parent)
        d.parent = origin.parent;
    if (!d.trial)
        d.trial = origin.trial;
    if (!d.numerical_generation)
        d.numerical_generation = origin.numerical_generation ? origin.numerical_generation
                                                             : std::optional{context.numerical};
    if (!d.last_accepted)
        d.last_accepted = origin.last_accepted;
    if (!d.resources)
        d.resources = origin.resources;
    if (!d.component)
        d.component = origin.component;
    if (!d.problem)
        d.problem = origin.problem;
    return {failure, problem, context};
}
OwnedGeometryFailure packet_failure(core::FailureCode code, const GeometryContext& context,
                                    std::optional<GeometryProblem> problem = {}) noexcept {
    return packet_failure(core::Failure{code, operation}, context, problem);
}
FacetScatter facet_scatter(const FacetStorage& facet) noexcept {
    if (const auto* quad = std::get_if<QuadOrigin>(&facet.origin)) {
        if (quad->half == 0)
            return QuadLRScatter{{1, 0, 0, 1, 1, 0}};
        return QuadLRScatter{{1, 0, 0, 1, 0, 1}};
    }
    return IndependentP1Scatter{{1, 0, 0, 0, 1, 0, 0, 0, 1}};
}
} // namespace

std::optional<core::FailureCode> LocalAccuracy::validation_error() const noexcept {
    if (!std::isfinite(tau))
        return core::FailureCode::Nonfinite;
    if (tau < 0 ||
        (representation != core::Precision::Binary64 && representation != core::Precision::TwoWord))
        return core::FailureCode::InvalidInput;
    return std::nullopt;
}

struct LocalPacketData {
    OwnedGeometry snapshot;
    LocalAmbientDescriptor ambient;
    QueryRequest expression;
    LocalAccuracy requested;
    GeometryLimits admitted;
    std::array<std::optional<ScaledEnclosure>, 12> values;
    std::optional<LocalQualification> qualification;
    std::optional<ClosestEntity> closest;
    FarRelation far;
    std::uint64_t steps;
    std::optional<LocalPacketTimings> timings;
    LocalPacketData(OwnedGeometry owner, LocalAmbientDescriptor layout, QueryRequest query,
                    LocalAccuracy accuracy, GeometryLimits limits) noexcept
        : snapshot(std::move(owner)),
          ambient(layout),
          expression(query),
          requested(accuracy),
          admitted(limits),
          far(FarRelation::Near),
          steps(0) {}
};
struct LocalPacketStorage {
    std::optional<LocalPacketData> payload;
};
class LocalPacketEngine {
    static bool midpoint(const QueryRequest& query, LocalRootScratch& scratch) {
        const auto* edge = std::get_if<SourceEdgePoint>(&query.definition);
        if (edge == nullptr)
            return false;
        if (const auto* single = std::get_if<SingleWord>(&edge->t))
            return single->value == .5;
        const auto* pair = std::get_if<DoubleWord>(&edge->t);
        if (pair == nullptr)
            return false;
        auto& sum = scratch.integers[0];
        auto& low = scratch.integers[1];
        auto& half = scratch.integers[2];
        local_decode(pair->high, sum);
        local_decode(pair->low, low);
        local_decode(.5, half);
        sum += low;
        return sum == half;
    }

    static AeroResult<LocalQualification> qualify(const LocalGeometryInput& input,
                                                  const NondegenerateTriangle& relation, bool bound,
                                                  std::optional<EdgeId> alias) noexcept {
        const auto& data = *input.owner;
        const auto& facet = data.facets.values()[static_cast<std::size_t>(input.facet)];
        const auto& support = input.query.support;
        constexpr std::array projection{ProjectionFeature::Edge01, ProjectionFeature::Edge12,
                                        ProjectionFeature::Edge20};
        const bool on =
            relation.plane == PlaneSide::On && relation.projection != ProjectionFeature::Outside;
        std::optional<std::variant<PointId, EdgeId, FacetId>> member;
        if (alias) {
            const auto* edge = std::get_if<EdgeSupportStorage>(&support);
            if (edge != nullptr && on) {
                for (std::size_t i = 0; i < facet.edges.size(); ++i) {
                    const auto id =
                        data.edges.values()[static_cast<std::size_t>(facet.edges[i])].id;
                    if (id == alias->value() && relation.projection == projection[i])
                        member = EdgeId{edge->id};
                }
            }
            if (!member)
                return std::unexpected(packet_failure(core::FailureCode::GeometryInvalid,
                                                      data.context, QueryId{input.query.id}));
        } else if (const auto* point = std::get_if<PointSupportStorage>(&support)) {
            for (const auto index : facet.points)
                if (data.points.values()[static_cast<std::size_t>(index)].id == point->id)
                    member = PointId{point->id};
        } else if (const auto* edge = std::get_if<EdgeSupportStorage>(&support)) {
            for (const auto index : facet.edges)
                if (data.edges.values()[static_cast<std::size_t>(index)].id == edge->id)
                    member = EdgeId{edge->id};
        } else if (const auto* face = std::get_if<FacetSupportStorage>(&support)) {
            if (facet.id == face->id)
                member = FacetId{face->id};
        }
        if (bound && on) {
            if (relation.projection == ProjectionFeature::Face)
                return BoundInterior{};
            const auto* origin = std::get_if<QuadOrigin>(&facet.origin);
            if (origin != nullptr) {
                const auto diagonal =
                    data.quads.values()[static_cast<std::size_t>(origin->source)].diagonal;
                for (std::size_t edge = 0; edge < projection.size(); ++edge)
                    if (relation.projection == projection[edge] &&
                        data.edges.values()[static_cast<std::size_t>(facet.edges[edge])].id ==
                            diagonal)
                        return BoundInterior{};
            }
        }
        if (member) {
            if (!on)
                return std::unexpected(packet_failure(core::FailureCode::GeometryInvalid,
                                                      data.context, QueryId{input.query.id}));
            return VerifiedIncidence{*member};
        }
        if (on)
            return std::unexpected(packet_failure(core::FailureCode::GeometryInvalid, data.context,
                                                  QueryId{input.query.id}));
        return Detached{};
    }
    static ClosestEntity nearest(const LocalGeometryInput& input,
                                 const NondegenerateTriangle& relation) noexcept {
        const auto& data = *input.owner;
        const auto& facet = data.facets.values()[static_cast<std::size_t>(input.facet)];
        std::optional<std::uint64_t> id;
        for (std::size_t i = 0; i < 3; ++i)
            if ((relation.minimizing_features & (1U << i)) != 0) {
                const auto next =
                    data.points.values()[static_cast<std::size_t>(facet.points[i])].id;
                if (!id || next < *id)
                    id = next;
            }
        if (id)
            return PointId{*id};
        for (std::size_t i = 0; i < 3; ++i)
            if ((relation.minimizing_features & (1U << (i + 3))) != 0) {
                const auto next = data.edges.values()[static_cast<std::size_t>(facet.edges[i])].id;
                if (!id || next < *id)
                    id = next;
            }
        if (id)
            return EdgeId{*id};
        return FacetId{facet.id};
    }

public:
    static AeroResult<LocalGeometryPacket>
    prepare(const OwnedGeometry& geometry, const ResolvedQuery& query, FacetId target,
            LocalAccuracy accuracy, const GeometryLimits& limits, LocalPacketWorkspace& workspace,
            core::ResourcePlan& plan, const std::stop_token& stop, bool bound = false,
            std::optional<EdgeId> alias = {}, std::uint64_t alias_steps = 0,
            SharedP1ReuseSlot<LocalGeometryPacket>* slot = nullptr) noexcept {
        const auto context = geometry.context();
        auto lease = workspace.acquire(context.diagnostics);
        if (!lease)
            return std::unexpected(packet_failure(lease.error(), context));
        try {
            if (stop.stop_requested())
                return std::unexpected(packet_failure(core::FailureCode::Cancelled, context));
            if (alias_steps > limits.max_steps || limits.max_steps - alias_steps < enclosure_steps)
                return std::unexpected(
                    packet_failure(core::FailureCode::ResourceExceeded, context));
            if (const auto error = accuracy.validation_error())
                return std::unexpected(packet_failure(*error, context));
            auto expression = query.definition();
            if (!expression)
                return std::unexpected(expression.error());
            auto& scratch = workspace.scratch_.values()[0];
            if (alias && !midpoint(*expression, scratch))
                return std::unexpected(
                    packet_failure(core::FailureCode::InvalidInput, context, expression->id));
            auto local_limits = limits;
            local_limits.max_steps -= alias_steps;
            LocalPacketTimings measured;
            auto input = [&] {
                const LocalClock clock{accuracy.collect_timings ? &measured.classification
                                                                : nullptr};
                return GeometryPacketAccess::load(geometry, query, target, local_limits,
                                                  *workspace.geometry_, *lease, plan, stop,
                                                  slot != nullptr ? slot->predicate_ : nullptr);
            }();
            if (!input) {
                auto failure = input.error();
                if (failure.cause.code == core::FailureCode::ResourceExceeded ||
                    failure.cause.code == core::FailureCode::Cancelled)
                    failure.problem.reset();
                return std::unexpected(failure);
            }
            const auto* outcome = input->outcome();
            if (!outcome)
                return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context));
            const auto* relation = std::get_if<NondegenerateTriangle>(outcome);
            if (!relation)
                return std::unexpected(
                    packet_failure(core::FailureCode::GeometryInvalid, context, target));
            auto qualification = [&] {
                const LocalClock clock{accuracy.collect_timings ? &measured.qualification
                                                                : nullptr};
                return qualify(*input, *relation, bound, alias);
            }();
            if (!qualification)
                return std::unexpected(qualification.error());
            auto& radius = scratch.integers[11];
            auto& temp = scratch.integers[12];
            const auto words = input->owner->radius;
            if (const auto* single = std::get_if<SingleWord>(&words))
                local_decode(single->value, radius);
            else if (const auto* pair = std::get_if<DoubleWord>(&words)) {
                local_decode(pair->high, radius);
                local_decode(pair->low, temp);
                radius += temp;
            } else
                return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context));
            if (std::holds_alternative<VerifiedIncidence>(*qualification) && radius == 0)
                return std::unexpected(
                    packet_failure(core::FailureCode::GeometryInvalid, context, target));
            if (stop.stop_requested())
                return std::unexpected(packet_failure(core::FailureCode::Cancelled, context));
            if (input->preparation_steps > local_limits.max_steps ||
                local_limits.max_steps - input->preparation_steps < enclosure_steps)
                return std::unexpected(
                    packet_failure(core::FailureCode::ResourceExceeded, context));
            const auto& facet =
                input->owner->facets.values()[static_cast<std::size_t>(input->facet)];
            core::PlannedBuffer<LocalPacketStorage> fresh;
            if (slot == nullptr || !slot->value_) {
                auto storage = plan.allocate<LocalPacketStorage>(1, operation);
                if (!storage)
                    return std::unexpected(packet_failure(storage.error(), context));
                if (slot != nullptr)
                    slot->value_.emplace(LocalGeometryPacket{std::move(*storage), context});
                else
                    fresh = std::move(*storage);
            }
            auto* storage = &fresh;
            if (slot != nullptr) {
                auto& current = slot->value_;
                if (!current)
                    return std::unexpected(
                        packet_failure(core::FailureCode::InvalidInput, context));
                storage = &current->storage_;
            }
            const auto& data = *input->owner;
            std::array<PointId, 3> vertices{PointId{0}, PointId{0}, PointId{0}};
            std::array<EdgeId, 3> edges{EdgeId{0}, EdgeId{0}, EdgeId{0}};
            for (std::size_t i = 0; i < 3; ++i) {
                vertices[i] =
                    PointId{data.points.values()[static_cast<std::size_t>(facet.points[i])].id};
                edges[i] = EdgeId{data.edges.values()[static_cast<std::size_t>(facet.edges[i])].id};
            }
            std::variant<std::array<PointId, 3>, std::array<PointId, 4>> physical{vertices};
            if (const auto* quad = std::get_if<QuadOrigin>(&facet.origin)) {
                const auto& original = data.quads.values()[static_cast<std::size_t>(quad->source)];
                physical = std::array{PointId{original.points[0]}, PointId{original.points[1]},
                                      PointId{original.points[2]}, PointId{original.points[3]}};
            }
            const auto& face = data.faces.values()[static_cast<std::size_t>(facet.face)];
            const LocalAmbientDescriptor layout{FacetId{facet.id},
                                                PhysicalFaceId{face.id},
                                                BandId{face.band},
                                                vertices,
                                                edges,
                                                facet.incidence,
                                                physical,
                                                facet_scatter(facet)};
            auto& payload = storage->values()[0].payload;
            std::optional<std::array<std::optional<ScaledEnclosure>, 7>> retained;
            if (payload && payload->steps != 0 && payload->snapshot.data_.get() == &data &&
                payload->ambient.facet == target &&
                payload->requested.representation == accuracy.representation &&
                payload->requested.tau == accuracy.tau) {
                // These slots depend only on the three immutable facet vertices.
                const auto constants = std::span{payload->values}.subspan<3, 7>();
                auto& cached = retained.emplace();
                std::ranges::copy(constants, cached.begin());
            }
            auto& packet = payload.emplace(GeometryPacketAccess::retain(*input), layout,
                                           *expression, accuracy, limits);
            packet.qualification = *qualification;
            packet.closest = nearest(*input, *relation);
            packet.far = relation->far;
            const auto& points = input->points.get();
            auto e = std::span{scratch.integers}.subspan<13, 3>();
            auto f = std::span{scratch.integers}.subspan<16, 3>();
            auto c = std::span{scratch.integers}.subspan<19, 3>();
            auto& j2 = scratch.integers[22];
            if (!retained) {
                for (std::size_t i = 0; i < 3; ++i) {
                    e[i] = points[1][i] - points[0][i];
                    f[i] = points[2][i] - points[0][i];
                }
                for (std::size_t i = 0; i < 3; ++i) {
                    const auto j = (i + 1) % 3, k = (i + 2) % 3;
                    c[i] = e[j] * f[k];
                    c[i] -= e[k] * f[j];
                }
                j2 = 0;
                for (const auto& component : c) {
                    temp = component * component;
                    j2 += temp;
                }
            }
            LocalWork work{alias_steps + input->preparation_steps, limits.max_steps, stop};
            {
                const LocalClock enclosure_clock{accuracy.collect_timings ? &measured.enclosure
                                                                          : nullptr};
                for (std::size_t i = 0; i < 12; ++i) {
                    if (retained && i >= 3 && i < 10) {
                        if (auto check = work.step(); !check)
                            return std::unexpected(packet_failure(check.error(), context));
                        packet.values[i] = (*retained)[i - 3];
                        continue;
                    }
                    core::Result<ScaledEnclosure> value = std::unexpected(
                        core::Failure{core::FailureCode::PrecisionUnresolved, operation});
                    if (i < 3)
                        value = LocalBounds::dyadic(points[3][i], -2148, accuracy.representation,
                                                    scratch);
                    else if (i < 6) {
                        radius = c[i - 3] * c[i - 3];
                        value = LocalBounds::root(
                            radius, j2, 0, c[i - 3] < 0, accuracy.representation, scratch, work,
                            accuracy.collect_timings ? &measured.root : nullptr);
                    } else if (i < 9)
                        value =
                            LocalBounds::dyadic(c[i - 6], -4297, accuracy.representation, scratch);
                    else if (i == 9) {
                        temp = 1;
                        value =
                            LocalBounds::root(relation->longest_edge_squared, temp, -2148, false,
                                              accuracy.representation, scratch, work,
                                              accuracy.collect_timings ? &measured.root : nullptr);
                    } else if (i == 10)
                        value =
                            LocalBounds::root(relation->distance_squared.numerator,
                                              relation->distance_squared.denominator, -2148, false,
                                              accuracy.representation, scratch, work,
                                              accuracy.collect_timings ? &measured.root : nullptr);
                    else {
                        if (const auto* single = std::get_if<SingleWord>(&words))
                            local_decode(single->value, radius);
                        else if (const auto* pair = std::get_if<DoubleWord>(&words)) {
                            local_decode(pair->high, radius);
                            local_decode(pair->low, temp);
                            radius += temp;
                        }
                        value =
                            LocalBounds::dyadic(radius, -1074, accuracy.representation, scratch);
                    }
                    if (!value)
                        return std::unexpected(packet_failure(value.error(), context, target));
                    if (auto check = work.step(); !check)
                        return std::unexpected(packet_failure(check.error(), context));
                    if (auto check = work.step(); !check)
                        return std::unexpected(packet_failure(check.error(), context));
                    const bool width_achieved = [&] {
                        const LocalClock clock{accuracy.collect_timings ? &measured.width
                                                                        : nullptr};
                        return LocalBounds::width(*value, accuracy.tau,
                                                  relation->longest_edge_squared, i >= 6 && i < 9,
                                                  i >= 3 && i < 6, scratch);
                    }();
                    if (!width_achieved)
                        return std::unexpected(packet_failure(
                            core::FailureCode::PrecisionUnresolved, context, target));
                    packet.values[i] = *value;
                }
            }
            if (accuracy.collect_timings)
                packet.timings = measured;
            packet.steps = work.used;
            return LocalGeometryPacket{std::move(*storage), context};
        } catch (const std::bad_alloc&) {
            return std::unexpected(packet_failure(core::FailureCode::ResourceExceeded, context));
        } catch (...) {
            return std::unexpected(packet_failure(core::FailureCode::PrecisionUnresolved, context));
        }
    }
};
LocalPacketWorkspace::LocalPacketWorkspace(Key, std::unique_ptr<GeometryWorkspace> geometry,
                                           core::PlannedBuffer<LocalRootScratch> scratch) noexcept
    : geometry_(std::move(geometry)), scratch_(std::move(scratch)) {}
LocalPacketWorkspace::~LocalPacketWorkspace() = default;
std::uint64_t local_packet_step_ceiling(std::uint64_t facets,
                                        const GeometryLimits& limits) noexcept {
    // Three source vertices plus a material query with up to three vertices.
    constexpr std::uint64_t decode_steps = 3 * (3 + 3);
    // Initial classification, projection, six distance candidates, optional
    // face, seven minimum comparisons and the far-boundary comparison.
    constexpr std::uint64_t predicate_steps = 2 + 6 + 1 + 7 + 1;
    // Check the alias call and its search, then examine all three facet edges.
    const auto alias_steps = 2 + static_cast<std::uint64_t>(std::bit_width(facets)) + 3;
    const auto total = enclosure_steps + decode_steps +
                       std::min(limits.predicate.max_steps, predicate_steps) + alias_steps;
    return std::min(total, limits.max_steps);
}
core::Result<std::unique_ptr<LocalPacketWorkspace>>
LocalPacketWorkspace::create(core::ResourcePlan& plan, core::FailureDetails details) noexcept {
    auto geometry = GeometryWorkspace::create(plan, details);
    if (!geometry)
        return std::unexpected(geometry.error());
    auto scratch = plan.allocate<LocalRootScratch>(1, operation);
    if (!scratch)
        return std::unexpected(scratch.error());
    try {
        return std::make_unique<LocalPacketWorkspace>(Key{}, std::move(*geometry),
                                                      std::move(*scratch));
    } catch (const std::bad_alloc&) {
        return std::unexpected(
            core::Failure{core::FailureCode::ResourceExceeded, operation, details});
    } catch (...) {
        return std::unexpected(
            core::Failure{core::FailureCode::PrecisionUnresolved, operation, details});
    }
}
core::Result<PredicateLease> LocalPacketWorkspace::acquire(core::FailureDetails details) noexcept {
    return geometry_->acquire(details);
}
std::uint64_t LocalPacketWorkspace::payload_bytes() noexcept {
    return GeometryWorkspace::payload_bytes() + sizeof(LocalRootScratch);
}
LocalGeometryPacket::LocalGeometryPacket(core::PlannedBuffer<LocalPacketStorage> storage,
                                         GeometryContext context) noexcept
    : storage_(std::move(storage)), context_(context) {}
LocalGeometryPacket::LocalGeometryPacket(LocalGeometryPacket&&) noexcept = default;
LocalGeometryPacket& LocalGeometryPacket::operator=(LocalGeometryPacket&&) noexcept = default;
LocalGeometryPacket::~LocalGeometryPacket() = default;
const LocalPacketData* LocalGeometryPacket::data() const noexcept {
    if (storage_.size() != 1)
        return nullptr;
    const auto& payload = storage_.values()[0].payload;
    if (!payload)
        return nullptr;
    return &*payload;
}
bool LocalGeometryPacket::active() const noexcept {
    return data() != nullptr;
}
AeroResult<std::reference_wrapper<const ScaledEnclosure>>
LocalGeometryPacket::value(LocalQuantity quantity) const& noexcept {
    const auto* owned = data();
    const auto index = static_cast<std::size_t>(quantity);
    if (owned == nullptr || index >= 12)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    const auto& value = owned->values[index];
    if (!value)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return std::cref(*value);
}
AeroResult<QueryRequest> LocalGeometryPacket::query() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->expression;
}
AeroResult<FacetScatter> LocalGeometryPacket::scatter() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->ambient.scatter;
}
AeroResult<FarRelation> LocalGeometryPacket::far() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->far;
}
AeroResult<LocalAccuracy> LocalGeometryPacket::accuracy() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->requested;
}
AeroResult<GeometryLimits> LocalGeometryPacket::admitted_limits() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->admitted;
}
AeroResult<std::uint64_t> LocalGeometryPacket::steps() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->steps;
}
AeroResult<LocalAmbientDescriptor> LocalGeometryPacket::ambient() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->ambient;
}
AeroResult<LocalQualification> LocalGeometryPacket::qualification() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    const auto& value = owned->qualification;
    if (!value)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return *value;
}
AeroResult<ClosestEntity> LocalGeometryPacket::closest() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    const auto& value = owned->closest;
    if (!value)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return *value;
}
AeroResult<std::array<std::byte, 32>> LocalGeometryPacket::snapshot_identity() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->snapshot.identity();
}
AeroResult<std::reference_wrapper<const OwnedGeometry>>
LocalGeometryPacket::geometry() const& noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return std::cref(owned->snapshot);
}
AeroResult<std::optional<LocalPacketTimings>> LocalGeometryPacket::timings() const noexcept {
    const auto* owned = data();
    if (owned == nullptr)
        return std::unexpected(packet_failure(core::FailureCode::InvalidInput, context_));
    return owned->timings;
}
std::uint64_t LocalGeometryPacket::payload_bytes() noexcept {
    return sizeof(LocalPacketStorage);
}
AeroResult<LocalGeometryPacket>
prepare_local_packet(const OwnedGeometry& geometry, const ResolvedQuery& query, FacetId target,
                     LocalAccuracy accuracy, const GeometryLimits& limits,
                     LocalPacketWorkspace& workspace, core::ResourcePlan& plan,
                     const std::stop_token& stop,
                     SharedP1ReuseSlot<LocalGeometryPacket>* slot) noexcept {
    return LocalPacketEngine::prepare(geometry, query, target, accuracy, limits, workspace, plan,
                                      stop, false, {}, 0, slot);
}
AeroResult<LocalGeometryPacket>
prepare_bound_local_packet(const OwnedGeometry& geometry, const ResolvedQuery& query,
                           FacetId target, LocalAccuracy accuracy, const GeometryLimits& limits,
                           LocalPacketWorkspace& workspace, core::ResourcePlan& plan,
                           const std::stop_token& stop,
                           SharedP1ReuseSlot<LocalGeometryPacket>* slot) noexcept {
    return LocalPacketEngine::prepare(geometry, query, target, accuracy, limits, workspace, plan,
                                      stop, true, {}, 0, slot);
}
AeroResult<LocalGeometryPacket>
prepare_aliased_midpoint_packet(const GlobalContactState& contact, const ResolvedQuery& query,
                                FacetId target, bool bound, LocalAccuracy accuracy,
                                const GeometryLimits& limits, LocalPacketWorkspace& workspace,
                                core::ResourcePlan& plan, const std::stop_token& stop,
                                SharedP1ReuseSlot<LocalGeometryPacket>* slot) noexcept {
    const auto checked = contact.geometry();
    if (!checked)
        return std::unexpected(checked.error().origin);
    const auto& geometry = checked->get();
    const auto context = geometry.context();
    if (stop.stop_requested())
        return std::unexpected(packet_failure(core::FailureCode::Cancelled, context));
    if (limits.max_steps == 0)
        return std::unexpected(packet_failure(core::FailureCode::ResourceExceeded, context));
    const auto expression = query.definition();
    if (!expression)
        return std::unexpected(expression.error());
    const auto* edge = std::get_if<SourceEdgePoint>(&expression->definition);
    if (edge == nullptr)
        return std::unexpected(
            packet_failure(core::FailureCode::InvalidInput, context, expression->id));
    const auto relationship =
        contact.facet_edge_alias(target, edge->edge, limits.max_steps - 1, stop);
    if (!relationship)
        return std::unexpected(relationship.error().origin);
    const auto [alias, used_steps] = *relationship;
    return LocalPacketEngine::prepare(geometry, query, target, accuracy, limits, workspace, plan,
                                      stop, bound, alias, 1 + used_steps, slot);
}
} // namespace uvlas::aero::detail
