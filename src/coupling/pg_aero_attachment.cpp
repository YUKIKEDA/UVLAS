#include "pg_aero_attachment.hpp"

#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_kernel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace s = structure;
namespace p = a::p1;

namespace {
using enum core::FailureCode;
constexpr auto operation = core::Operation::EvaluateAero;
constexpr std::size_t axes = 3;

[[noreturn]] void reject(core::FailureCode code) {
    throw p::ArithmeticFailure{code};
}
template <class Result>
auto take(Result result) {
    if (!result) {
        const auto& error = result.error();
        if constexpr (std::is_same_v<typename Result::error_type, core::Failure>)
            throw error;
        else if constexpr (std::is_same_v<typename Result::error_type, a::GlobalGeometryFailure>)
            throw error.origin.cause;
        else
            throw error.cause;
    }
    if constexpr (!std::is_void_v<typename Result::value_type>)
        return std::move(*result);
}

template <class F>
auto guarded(F action) noexcept {
    using Result = core::Result<decltype(action())>;
    return core::guard_entry(operation, InvalidInput, {}, [&]() -> Result {
        try {
            return action();
        } catch (const core::Failure& failure) {
            return std::unexpected(failure);
        } catch (const p::ArithmeticFailure& failure) {
            return std::unexpected(core::Failure{failure.code, operation});
        }
    });
}

std::uint64_t remaining(const p::Work& work) noexcept {
    return work.limit - work.used;
}
void account(p::Work& work, std::uint64_t steps) {
    if (steps > remaining(work))
        reject(ResourceExceeded);
    work.used += steps;
}
core::Index components(std::size_t count) {
    if (count > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) / axes)
        reject(ResourceExceeded);
    return static_cast<core::Index>(count * axes);
}
template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& plan, std::size_t count) {
    if (count > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        reject(ResourceExceeded);
    return take(plan.allocate<T>(static_cast<core::Index>(count), operation));
}
std::size_t vertex_row(std::span<const std::optional<a::PointId>> ids, a::PointId owner,
                       p::Work& work) {
    const auto it = std::ranges::lower_bound(ids, std::optional{owner}, [&](auto l, auto r) {
        work.step();
        return l < r;
    });
    return it == ids.end() || *it != owner ? ids.size()
                                           : static_cast<std::size_t>(it - ids.begin());
}

a::ScalarWords component(s::PgCartesianWords words, std::size_t i) {
    if (words.representation == core::Precision::Binary64)
        return a::SingleWord{words.high[i]};
    return a::DoubleWord{words.high[i], words.low[i]};
}
NonpenetrationVector vector(s::PgCartesianWords words, std::size_t row) {
    const auto i = row * axes;
    return {component(words, i), component(words, i + 1), component(words, i + 2)};
}
a::FreePoint point(s::PgCartesianWords words, std::size_t row) {
    const auto i = row * axes;
    const std::array<double, axes> high{words.high[i], words.high[i + 1], words.high[i + 2]};
    if (words.representation == core::Precision::Binary64)
        return a::FreePoint{a::FreeSingle{high}};
    return a::FreePoint{a::FreeDouble{high, {words.low[i], words.low[i + 1], words.low[i + 2]}}};
}

p::Interval exact_word(a::ScalarWords value, p::Work& work) {
    if (const auto* single = std::get_if<a::SingleWord>(&value))
        return a::decode_p1_word(single->value);
    const auto& pair = std::get<a::DoubleWord>(value);
    return p::add(a::decode_p1_word(pair.high), a::decode_p1_word(pair.low), work);
}
double unscale(double word, int exponent) {
    if (word == 0)
        return word;
    if (exponent < -2200 || exponent > 2200)
        reject(PrecisionUnresolved);
    const auto physical = std::scalbn(word, exponent);
    if (!std::isfinite(physical) || std::scalbn(physical, -exponent) != word)
        reject(PrecisionUnresolved);
    return physical;
}
NonpenetrationVector midpoint(s::PgCartesianWords words, const std::array<std::size_t, 2>& rows,
                              double error, double endpoint_error, p::Work& work) {
    a::LocalRootScratch scratch;
    const auto axis = [&](std::size_t k) -> a::ScalarWords {
        const auto first = exact_word(component(words, rows[0] * axes + k), work);
        const auto second = exact_word(component(words, rows[1] * axes + k), work);
        const auto average = p::scale(p::add(first, second, work), -1);
        const auto inherited = a::decode_p1_word(endpoint_error);
        const p::Interval qualified{p::subtract(average, inherited, work).lower,
                                    p::add(average, inherited, work).upper};
        const auto projected = a::enclose_p1_value(qualified, words.representation, scratch);
        const auto radius = p::scale(a::decode_p1_word(projected.radius()), projected.exponent());
        if (p::compare(radius.upper, a::decode_p1_word(error).lower) > 0)
            reject(PrecisionUnresolved);
        const auto center = projected.center();
        if (const auto* single = std::get_if<a::SingleWord>(&center))
            return a::SingleWord{unscale(single->value, projected.exponent())};
        const auto& pair = std::get<a::DoubleWord>(center);
        return a::DoubleWord{unscale(pair.high, projected.exponent()),
                             unscale(pair.low, projected.exponent())};
    };
    return {axis(0), axis(1), axis(2)};
}
} // namespace

struct PgAeroMap::Data {
    struct Row {
        a::FacetId facet{0};
        a::QueryId query{0};
        bool reverse_orientation = false;
    };
    s::PgProjectedReferenceBatch position;
    std::variant<s::PgBoundaryTangentReference, s::PgMaterialTangentReference> tangents;
    core::PlannedBuffer<std::optional<a::PointId>> vertices;
    core::PlannedBuffer<Row> rows;
};

PgAeroMap::PgAeroMap(std::shared_ptr<const Data> data) noexcept : data_(std::move(data)) {}

s::PgPreparedPatchIdentity PgAeroMap::source() const noexcept {
    return data_ == nullptr ? s::PgPreparedPatchIdentity{} : data_->position.source();
}

template <class Position, class Tangents>
core::Result<PgAeroMap> PgAeroMap::prepare_impl(
    std::span<const a::PointId> vertices, std::span<const PgAeroCollocation> rows,
    std::size_t material_count, const s::PgReferenceLimits& limits, core::ResourcePlan& plan,
    const Position& prepare_position, const Tangents& prepare_tangents) noexcept {
    return guarded([&] {
        if (limits.cancellation.stop_requested())
            reject(Cancelled);
        if (vertices.empty() || rows.empty() || vertices.size() > material_count ||
            rows.size() != material_count - vertices.size())
            reject(InvalidInput);
        if (material_count > limits.maximum_partition_visits)
            reject(ResourceExceeded);
        auto ids = allocate<std::optional<a::PointId>>(plan, vertices.size());
        auto rules = allocate<Data::Row>(plan, rows.size());
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (limits.cancellation.stop_requested())
                reject(Cancelled);
            if (i != 0 && vertices[i - 1] >= vertices[i])
                reject(InvalidInput);
            ids.values()[i] = vertices[i];
        }
        for (std::size_t i = 0; i < rows.size(); ++i) {
            if (limits.cancellation.stop_requested())
                reject(Cancelled);
            const auto& row = rows[i];
            rules.values()[i] = {row.facet, row.query, row.reverse_orientation};
        }
        auto budget = limits;
        budget.maximum_partition_visits -= material_count;
        auto position = take(prepare_position(budget));
        budget.maximum_partition_visits -= position.call_info().partition_visits;
        budget.maximum_scalar_operations -= position.call_info().scalar_operations;
        auto tangents = take(prepare_tangents(budget));
        return PgAeroMap{std::make_shared<const Data>(
            Data{std::move(position), std::move(tangents), std::move(ids), std::move(rules)})};
    });
}
core::Result<PgAeroMap> PgAeroMap::prepare(const s::PgPreparedPatch& patch,
                                           std::span<const a::PointId> vertices,
                                           std::span<const PgAeroCollocation> rows,
                                           std::span<const s::PgBoundaryPositionPoint> material,
                                           const s::PgPreparedReferenceProjection& projection,
                                           double tangent_error, const s::PgReferenceLimits& limits,
                                           core::ResourcePlan& plan) noexcept {
    return prepare_impl(
        vertices, rows, material.size(), limits, plan,
        [&](const s::PgReferenceLimits& budget) {
            return patch.project_boundary_reference(material, projection, budget, plan);
        },
        [&](const s::PgReferenceLimits& budget) {
            return s::PgBoundaryTangentReference::prepare(
                patch.conditioning(), material.subspan(vertices.size()),
                {projection.representation, tangent_error}, budget, plan);
        });
}

core::Result<PgAeroMap> PgAeroMap::prepare_material(
    const s::PgPreparedPatch& patch, std::span<const a::PointId> vertices,
    std::span<const PgAeroCollocation> rows, std::span<const s::PgReferencePoint> material,
    std::span<const s::PgMaterialDirections> directions,
    const s::PgPreparedReferenceProjection& projection, double tangent_error,
    const s::PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return prepare_impl(
        vertices, rows, material.size(), limits, plan,
        [&](const s::PgReferenceLimits& budget) {
            return patch.project_reference(material, projection, budget, plan);
        },
        [&](const s::PgReferenceLimits& budget) {
            return s::PgMaterialTangentReference::prepare(
                patch, material.subspan(vertices.size()), directions,
                {projection.representation, tangent_error}, budget, plan);
        });
}

struct PgAeroTrialRow {
    a::FacetId facet;
    a::ResolvedQuery query;
    ExplicitCollocationTangents tangents;
    NonpenetrationVector velocity;
};

struct PgAeroTrialData {
    std::shared_ptr<const PgAeroMap::Data> map;
    std::array<std::byte, 32> geometry{};
    a::P1ActionStamp stamp{};
    std::optional<a::GlobalContactState> contact;
    core::PlannedBuffer<std::optional<PgAeroTrialRow>> rows;
    core::PlannedBuffer<std::optional<EdgeMidpointVelocity>> midpoints;
};

PgAeroTrial::PgAeroTrial(core::PlannedBuffer<PgAeroTrialData>&& data) noexcept
    : data_(std::move(data)) {}
PgAeroTrial::PgAeroTrial(PgAeroTrial&&) noexcept = default;
PgAeroTrial& PgAeroTrial::operator=(PgAeroTrial&&) noexcept = default;
PgAeroTrial::~PgAeroTrial() = default;

core::Result<std::reference_wrapper<const a::GlobalContactState>>
PgAeroTrial::contact() const& noexcept {
    const auto trials = data_.values();
    if (!trials.empty())
        if (const auto& candidate = trials.front().contact)
            return std::cref(*candidate);
    return std::unexpected(core::Failure{InvalidInput, operation});
}
std::size_t PgAeroTrial::collocation_count() const noexcept {
    return data_.size() == 0 ? 0 : data_.values().front().rows.size();
}
core::Result<DiscreteCollocation> PgAeroTrial::collocation(std::size_t i) const& noexcept {
    if (i < collocation_count()) {
        if (const auto& row = data_.values().front().rows.values()[i]) {
            const a::ScalarWords zero = a::SingleWord{0};
            return DiscreteCollocation{row->facet,    row->query, row->tangents,
                                       row->velocity, zero,       zero};
        }
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}
std::size_t PgAeroTrial::midpoint_count() const noexcept {
    return data_.size() == 0 ? 0 : data_.values().front().midpoints.size();
}
core::Result<EdgeMidpointVelocity> PgAeroTrial::midpoint(std::size_t i) const noexcept {
    if (i < midpoint_count()) {
        if (const auto& value = data_.values().front().midpoints.values()[i])
            return *value;
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}

core::Result<core::PlannedBuffer<std::optional<DiscreteLoadVector>>>
PgAeroTrial::pullback_nodal_forces(const DiscreteCirculationLoads& forces,
                                   core::Precision representation, double maximum_error,
                                   std::uint64_t maximum_steps, core::ResourcePlan& memory,
                                   const std::stop_token& stop) const noexcept {
    return guarded([&] {
        p::Work work{0, maximum_steps, stop};
        work.step();
        if (data_.size() == 0 || !data_.values().front().map)
            reject(InvalidInput);
        if (representation != core::Precision::Binary64 &&
            representation != core::Precision::TwoWord)
            reject(PrecisionUnsupported);
        if (!std::isfinite(maximum_error))
            reject(Nonfinite);
        if (maximum_error < 0)
            reject(InvalidInput);
        const auto& trial = data_.values().front();
        const auto source = forces.source();
        const auto& geometry = take(take(contact()).get().geometry()).get();
        if (source.geometry != trial.geometry || source.context != geometry.context() ||
            source.stamp != trial.stamp)
            reject(StaleGeneration);
        const auto& map = *trial.map;
        const auto count = forces.point_count();
        if (count != map.vertices.size())
            reject(InvalidInput);

        const auto rank = map.position.column_count();
        const auto high = map.position.high_words();
        const auto low = map.position.low_words();
        const auto shape_error = a::decode_p1_word(map.position.projection().maximum_shape_error);
        using IntervalVector = std::array<p::Interval, axes>;
        auto loads = allocate<IntervalVector>(memory, count);
        for (std::size_t i = 0; i < count; ++i) {
            work.step();
            const auto point_load = take(forces.point(i));
            if (map.vertices.values()[i] != point_load.point)
                reject(InvalidInput);
            for (std::size_t k = 0; k < axes; ++k) {
                work.step();
                const auto& word = point_load.force[k];
                const auto center = exact_word(word.center(), work);
                const auto width = a::decode_p1_word(word.radius());
                loads.values()[i][k] = p::scale(
                    {p::subtract(center, width, work).lower, p::add(center, width, work).upper},
                    word.exponent());
            }
        }
        a::LocalRootScratch scratch;
        const auto allowance = a::decode_p1_word(maximum_error);
        auto output = allocate<std::optional<DiscreteLoadVector>>(memory, rank);
        for (std::size_t column = 0; column < rank; ++column) {
            IntervalVector total{};
            for (std::size_t i = 0; i < count; ++i) {
                work.step();
                // X0 and measure precede the column's shape and its three gradients.
                const auto index = 4 * (i * (rank + 1) + column + 1);
                const auto n = low.empty()
                                   ? a::decode_p1_word(high[index])
                                   : exact_word(a::DoubleWord{high[index], low[index]}, work);
                const p::Interval coefficient{p::subtract(n, shape_error, work).lower,
                                              p::add(n, shape_error, work).upper};
                for (std::size_t k = 0; k < axes; ++k)
                    total[k] = p::add(total[k],
                                      p::multiply(coefficient, loads.values()[i][k], work), work);
            }
            const auto project = [&](const p::Interval& value) {
                work.step();
                const auto word = a::enclose_p1_value(value, representation, scratch);
                const auto error = p::scale(a::decode_p1_word(word.radius()), word.exponent());
                if (p::compare(error.upper, allowance.lower) > 0)
                    reject(PrecisionUnresolved);
                return word;
            };
            output.values()[column].emplace(
                DiscreteLoadVector{project(total[0]), project(total[1]), project(total[2])});
        }
        if (stop.stop_requested())
            reject(Cancelled);
        return output;
    });
}

struct PgAeroBuilder {
    static PgAeroTrial evaluate(const PgAeroMap& map, const PgAeroState& state,
                                const PgAeroChart& chart, const PgAeroLimits& limits,
                                core::Workspace& scratch, a::GeometryWorkspace& coordinates,
                                a::GlobalWorkspace& global, core::ResourcePlan& memory,
                                const std::stop_token& stop) {
        p::Work work{0, limits.maximum_steps, stop};
        work.step();
        const auto& original = chart.geometry;
        if (map.data_ == nullptr)
            reject(InvalidInput);
        if (state.patch != map.source() || state.context != original.context)
            reject(StaleGeneration);
        if (map.data_->vertices.size() == 0 || map.data_->rows.size() == 0 ||
            chart.geometry.points.size() < map.data_->vertices.size())
            reject(InvalidInput);
        for (const auto error :
             {limits.maximum_midpoint_velocity_error, limits.maximum_reference_geometry_error}) {
            if (!std::isfinite(error))
                reject(Nonfinite);
            if (error < 0)
                reject(InvalidInput);
        }
        if (map.data_->rows.size() >
            static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) / 6)
            reject(ResourceExceeded);
        const auto format = map.data_->position.projection().representation;
        const auto size = components(map.data_->position.points().size());
        const auto tangent_size = components(map.data_->rows.size()) * 2;
        auto positions = take(s::PgVolumePositionResult::allocate(size, format, memory));
        auto velocities = take(s::PgVolumePositionResult::allocate(size, format, memory));
        auto tangents = take(s::PgVolumePositionResult::allocate(tangent_size, format, memory));
        const auto apply = [&](s::PgCartesianWords coefficients, s::PgVolumePositionResult& result,
                               s::PgVolumePositionLimits budget) {
            work.step();
            account(work, budget.maximum_steps);
            if (budget.cancellation.stop_requested())
                reject(Cancelled);
            budget.cancellation = stop;
            take(s::apply_pg_prepared_reference_position(map.data_->position, state.patch,
                                                         s::PgVolumePositionDirection::ToPoints,
                                                         coefficients, budget, scratch, result));
        };
        apply(state.position, positions, limits.position);
        apply(state.velocity, velocities, limits.velocity);
        account(work, limits.tangents.maximum_steps);
        if (limits.tangents.cancellation.stop_requested())
            reject(Cancelled);
        auto tangent_limits = limits.tangents;
        tangent_limits.cancellation = stop;
        std::visit(
            [&](const auto& map_tangents) {
                constexpr auto direction = s::PgVolumePositionDirection::ToPoints;
                using Owner = std::decay_t<decltype(map_tangents)>;
                if constexpr (std::is_same_v<Owner, s::PgBoundaryTangentReference>)
                    take(s::apply_pg_boundary_tangents(map_tangents, direction, state.position,
                                                       tangent_limits, scratch, tangents));
                else
                    take(s::apply_pg_material_tangents(map_tangents, state.patch, direction,
                                                       state.position, tangent_limits, scratch,
                                                       tangents));
            },
            map.data_->tangents);

        const auto& base = take(chart.topology.geometry()).get();
        const auto word_count = static_cast<std::size_t>(components(original.points.size()));
        const auto* single = std::get_if<a::SinglePositions>(&original.current);
        if ((single != nullptr) != (format == core::Precision::Binary64))
            reject(PrecisionUnsupported);
        const auto original_high =
            single != nullptr ? single->words : std::get<a::DoublePositions>(original.current).high;
        const auto original_low = single != nullptr
                                      ? std::span<const double>{}
                                      : std::get<a::DoublePositions>(original.current).low;
        if (original_high.size() != word_count ||
            original_low.size() != (single != nullptr ? 0 : word_count))
            reject(InvalidInput);
        auto high = allocate<double>(memory, word_count);
        auto low = allocate<double>(memory, original_low.size());
        const auto xs = positions.words();
        const auto vs = velocities.words();
        const auto ids = map.data_->vertices.values();
        auto used_vertices = allocate<std::uint8_t>(memory, ids.size());
        const auto reference_shape = map.data_->position.column_count();
        const auto projection = map.data_->position.projection();
        const s::PgCartesianWords reference_words{projection.representation,
                                                  map.data_->position.high_words(),
                                                  map.data_->position.low_words()};
        const auto radius = a::decode_p1_word(projection.maximum_position_error_metres);
        const auto tolerance = a::decode_p1_word(limits.maximum_reference_geometry_error);
        const auto threshold = p::multiply(tolerance, tolerance, work);
        for (std::size_t i = 0; i < ids.size(); ++i) {
            work.step();
            const auto& id = ids[i];
            if (!id || take(base.owner(*id)) != *id)
                reject(InvalidInput);
            const auto reference_xyz = take(base.point_words(*id, a::Configuration::Reference));
            p::Interval distance_sq{};
            for (std::size_t k = 0; k < axes; ++k) {
                const auto mapped =
                    exact_word(component(reference_words, 4 * (reference_shape + 1) * i + k), work);
                const auto difference =
                    p::subtract(mapped, exact_word(reference_xyz[k], work), work);
                const auto magnitude = p::absolute({p::subtract(difference, radius, work).lower,
                                                    p::add(difference, radius, work).upper});
                distance_sq = p::add(distance_sq, p::multiply(magnitude, magnitude, work), work);
            }
            if (p::compare(distance_sq.lower, threshold.upper) > 0)
                reject(GeometryInvalid);
            if (p::compare(distance_sq.upper, threshold.lower) > 0)
                reject(PrecisionUnresolved);
        }
        for (std::size_t i = 0; i < original.points.size(); ++i) {
            work.step();
            const auto owner = take(base.owner(original.points[i]));
            const auto row = vertex_row(ids, owner, work);
            const auto copy = [&](auto dst, auto raw, auto mapped) {
                if (dst.empty())
                    return;
                const auto words = row == ids.size() ? raw.subspan(i * axes, axes)
                                                     : mapped.subspan(row * axes, axes);
                std::ranges::copy(words, dst.subspan(i * axes, axes).begin());
            };
            copy(high.values(), original_high, xs.high);
            copy(low.values(), original_low, xs.low);
        }
        auto candidate = original;
        candidate.current = single != nullptr
                                ? a::PositionView{a::SinglePositions{high.values()}}
                                : a::PositionView{a::DoublePositions{high.values(), low.values()}};
        account(work, limits.geometry.max_steps);
        const auto geometry =
            take(a::prepare_owned_geometry(candidate, limits.geometry, coordinates, memory, stop));
        auto contact_limits = limits.contact;
        contact_limits.max_steps = std::min(contact_limits.max_steps, remaining(work));
        auto contact = take(a::check_global_geometry(chart.topology, geometry, contact_limits,
                                                     global, memory, stop));
        account(work, take(contact.work_plan()).executed_steps);
        auto result = allocate<PgAeroTrialData>(memory, 1);
        auto& out = result.values().front();
        out.rows = allocate<std::optional<PgAeroTrialRow>>(memory, map.data_->rows.size());
        out.midpoints =
            allocate<std::optional<EdgeMidpointVelocity>>(memory, chart.load_edges.size());
        const auto resolve = [&](a::QueryId id, a::QueryDefinition definition) {
            account(work, limits.geometry.max_steps);
            return take(a::resolve_query(geometry, {id, definition}, state.context, limits.geometry,
                                         coordinates, memory, stop));
        };
        for (std::size_t i = 0; i < map.data_->rows.size(); ++i) {
            work.step();
            const auto& row = map.data_->rows.values()[i];
            const auto panel = take(contact.uvlm_quad(row.facet, remaining(work), stop));
            account(work, panel.used_steps);
            if (panel.role != a::AeroFaceRole::Bound)
                reject(InvalidInput);
            for (const auto facet : panel.facets)
                for (const auto vertex : take(geometry.facet_points(facet))) {
                    const auto ordinal = vertex_row(ids, take(geometry.owner(vertex)), work);
                    if (ordinal == ids.size())
                        reject(InvalidInput);
                    used_vertices.values()[ordinal] = 1;
                }
            const auto position_row = ids.size() + i;
            auto query = resolve(row.query, point(xs, position_row));
            const auto derivatives = tangents.words();
            ExplicitCollocationTangents frame{vector(derivatives, 2 * i),
                                              vector(derivatives, 2 * i + 1)};
            if (row.reverse_orientation)
                std::swap(frame[0], frame[1]);
            out.rows.values()[i].emplace(
                PgAeroTrialRow{row.facet, std::move(query), frame, vector(vs, position_row)});
        }
        for (const auto used : used_vertices.values()) {
            work.step();
            if (used == 0)
                reject(InvalidInput);
        }
        for (std::size_t i = 0; i < chart.load_edges.size(); ++i) {
            work.step();
            const auto edge = chart.load_edges[i];
            if (i != 0 && chart.load_edges[i - 1] >= edge)
                reject(InvalidInput);
            const auto source = take(contact.uvlm_load_edge(edge, remaining(work), stop));
            account(work, source.used_steps);
            if (source.edge != edge || source.trailing)
                reject(InvalidInput);
            const std::array ordinal{
                vertex_row(ids, take(geometry.owner(source.endpoints[0])), work),
                vertex_row(ids, take(geometry.owner(source.endpoints[1])), work)};
            if (ordinal[0] == ids.size() || ordinal[1] == ids.size())
                reject(InvalidInput);
            out.midpoints.values()[i].emplace(
                edge, midpoint(vs, ordinal, limits.maximum_midpoint_velocity_error,
                               limits.velocity.maximum_absolute_error, work));
        }
        if (stop.stop_requested())
            reject(Cancelled);
        out.map = map.data_;
        out.geometry = take(geometry.identity());
        out.stamp = state.stamp;
        out.contact.emplace(std::move(contact));
        return PgAeroTrial{std::move(result)};
    }
};

core::Result<PgAeroTrial>
evaluate_pg_aero_trial(const PgAeroMap& map, const PgAeroState& state, const PgAeroChart& chart,
                       const PgAeroLimits& limits, core::Workspace& scratch,
                       a::GeometryWorkspace& coordinates, a::GlobalWorkspace& global,
                       core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return guarded([&] {
        return PgAeroBuilder::evaluate(map, state, chart, limits, scratch, coordinates, global,
                                       memory, stop);
    });
}

} // namespace uvlas::coupling::detail
