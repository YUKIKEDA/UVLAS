#include "wake_strength_transition.hpp"

#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_kernel.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <memory_resource>
#include <tuple>

namespace uvlas::coupling::detail {
namespace {
namespace a = aero::detail;
using Work = a::p1::Work;
using enum core::FailureCode;
constexpr auto operation = core::Operation::EvaluateAero;

[[noreturn]] void reject(core::FailureCode code) {
    throw core::Failure{code, operation};
}

template <class Result>
auto read(Result result) {
    if (!result) {
        if constexpr (std::is_same_v<typename Result::error_type, core::Failure>)
            throw result.error();
        else if constexpr (std::is_same_v<typename Result::error_type, a::GlobalGeometryFailure>)
            throw result.error().origin.cause;
        else
            throw result.error().cause;
    }
    return std::move(*result);
}

template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& memory, std::size_t size) {
    if (size > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        reject(ResourceExceeded);
    return read(memory.allocate<T>(static_cast<core::Index>(size), operation));
}

template <class Range, class Key, class Projection>
auto find(const Range& rows, Key key, Projection project, Work& work) {
    return std::ranges::lower_bound(
        rows, key,
        [&](const auto& lhs, const auto& rhs) {
            work.step();
            return lhs < rhs;
        },
        project);
}

std::array<double, 2> words(const a::ScalarWords& value, double sign = 1.) {
    if (const auto* one = std::get_if<a::SingleWord>(&value)) {
        if (!std::isfinite(one->value))
            reject(Nonfinite);
        return {sign * one->value, 0.};
    }
    const auto& two = std::get<a::DoubleWord>(value);
    if (!std::isfinite(two.high) || !std::isfinite(two.low))
        reject(Nonfinite);
    return {sign * two.high, sign * two.low};
}

struct WakeSpanBand {
    a::PointId from{0}, to{0};
    auto operator<=>(const WakeSpanBand&) const = default;
};

struct MaterialTraceKey {
    WakeSpanBand band;
    a::ReleaseEpoch epoch{0};
    auto operator<=>(const MaterialTraceKey&) const = default;
};

struct IndexedPoint {
    a::PointId owner{0};
    a::MaterialReleaseKey key{a::PointId{0}, a::ReleaseEpoch{0}};
};

struct Column {
    MaterialTraceKey key;
    a::P1StrengthId id;
    a::ScalarWords strength{a::SingleWord{0}};
    std::optional<MaterialTraceKey> leading_key;
    std::optional<a::P1StrengthId> leading_id;
    std::optional<WakeStrengthParent> parent;
};

struct AttachedColumn {
    a::EdgeId edge{0};
    WakeSpanBand band;
    std::array<double, 2> value{};
    WakeStrengthParent parent;
};

struct Section {
    core::PlannedBuffer<Column> frozen;
    core::PlannedBuffer<AttachedColumn> attached;
};

template <class Rows, class Projection>
void sort_unique(Rows rows, Projection project, Work& work) {
    std::ranges::sort(
        rows,
        [&](const auto& lhs, const auto& rhs) {
            work.step();
            return lhs < rhs;
        },
        project);
    for (std::size_t i = 1; i < rows.size(); ++i) {
        work.step();
        if (std::invoke(project, rows[i - 1]) == std::invoke(project, rows[i]))
            reject(InvalidInput);
    }
}

a::UvlmQuad panel(const a::GlobalContactState& contact, a::FacetId facet, Work& work) {
    const auto result = read(contact.uvlm_quad(facet, work.limit - work.used, work.stop));
    work.used += result.used_steps;
    return result;
}

Section describe(const WakeStrengthEndpoint& endpoint, core::ResourcePlan& memory, Work& work) {
    const auto& source = endpoint.source;
    const auto& geometry = read(source.contact.geometry()).get();
    if (geometry.context().diagnostics.parent != source.history.parent)
        reject(StaleGeneration);
    const auto counts = read(geometry.counts());
    if (counts.faces != counts.quads || counts.triangles != 0 ||
        source.collocations.size() > counts.faces ||
        source.wake.size() != counts.faces - source.collocations.size() ||
        source.wake.size() != source.history.columns.size())
        reject(InvalidInput);

    auto points = allocate<IndexedPoint>(memory, endpoint.material.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        work.step();
        const auto& row = endpoint.material[i];
        points.values()[i] = {read(geometry.owner(row.point)), row.release};
    }
    sort_unique(points.values(), &IndexedPoint::owner, work);
    const auto material = [&](a::PointId point) {
        const auto id = read(geometry.owner(point));
        const auto at = find(points.values(), id, &IndexedPoint::owner, work);
        if (at == points.values().end() || at->owner != id)
            reject(StaleGeneration);
        return at->key;
    };
    const auto key = [&](a::PointId left, a::PointId right) {
        const auto l = material(left);
        const auto r = material(right);
        if (l.epoch != r.epoch || l.point == r.point)
            reject(StaleGeneration);
        return MaterialTraceKey{{l.point, r.point}, l.epoch};
    };

    struct Bound {
        a::PhysicalFaceId face{0};
        a::ScalarWords gamma{a::SingleWord{0}};
    };
    auto bound = allocate<Bound>(memory, source.collocations.size());
    auto unique_faces = allocate<std::uint64_t>(memory, static_cast<std::size_t>(counts.faces));
    for (std::size_t i = 0; i < bound.size(); ++i) {
        const auto& row = source.collocations[i];
        const auto quad = panel(source.contact, row.facet, work);
        if (quad.role != a::AeroFaceRole::Bound)
            reject(InvalidInput);
        bound.values()[i] = {quad.face, row.circulation};
        unique_faces.values()[i] = quad.face.value();
    }
    sort_unique(bound.values(), &Bound::face, work);
    const auto history = source.history.columns;
    for (std::size_t i = 0; i < history.size(); ++i) {
        work.step();
        if (!history[i].id.wake_value() || (i != 0 && history[i - 1].id >= history[i].id))
            reject(InvalidInput);
    }

    auto traces = allocate<Column>(memory, source.wake.size());
    auto attached = allocate<AttachedColumn>(memory, source.contact.trailing_edge_count());
    auto identifiers = allocate<a::P1StrengthId>(memory, source.wake.size());
    std::size_t front = 0;
    for (std::size_t i = 0; i < source.wake.size(); ++i) {
        const auto& row = source.wake[i];
        const auto quad = panel(source.contact, row.facet, work);
        if (quad.role != a::AeroFaceRole::Wake)
            reject(InvalidInput);
        unique_faces.values()[bound.size() + i] = quad.face.value();
        const auto back = key(quad.corners[3], quad.corners[2]);
        const auto at = find(history, row.trailing, &FrozenWakeTrace::id, work);
        if (at == history.end() || at->id != row.trailing)
            reject(InvalidInput);
        auto& trace = traces.values()[i];
        (void)words(at->strength);
        trace = {back, row.trailing, at->strength, {}, {}, {}};
        identifiers.values()[i] = row.trailing;
        if (quad.attachment) {
            if (!std::holds_alternative<AttachedWakeFront>(row.leading) || front >= attached.size())
                reject(InvalidInput);
            const WakeSpanBand band{read(geometry.owner(quad.corners[0])),
                                    read(geometry.owner(quad.corners[1]))};
            if (band != back.band)
                reject(StaleGeneration);
            const auto& transfer = *quad.attachment;
            const auto gamma = find(bound.values(), transfer.bound_face, &Bound::face, work);
            if (gamma == bound.values().end() || gamma->face != transfer.bound_face)
                reject(InvalidInput);
            attached.values()[front++] = {
                transfer.edge,
                band,
                words(gamma->gamma, transfer.coefficient),
                {a::P1StrengthId::bound(transfer.bound_face), transfer.coefficient < 0}};
        } else {
            const auto* leading = std::get_if<a::P1StrengthId>(&row.leading);
            if (!leading)
                reject(InvalidInput);
            const auto upstream = key(quad.corners[0], quad.corners[1]);
            if (upstream.band != back.band || upstream.epoch <= back.epoch)
                reject(StaleGeneration);
            trace.leading_key = upstream;
            trace.leading_id = *leading;
        }
    }
    if (front != attached.size())
        reject(InvalidInput);
    sort_unique(unique_faces.values(), std::identity{}, work);
    sort_unique(identifiers.values(), std::identity{}, work);
    sort_unique(attached.values(), &AttachedColumn::edge, work);
    sort_unique(attached.values(), &AttachedColumn::band, work);
    sort_unique(traces.values(), &Column::key, work);

    // Every independent leading coefficient is the downstream slot of its actual
    // upstream quadrilateral, never a lookalike with the same value.
    const auto columns = traces.values();
    for (std::size_t i = 0; i < columns.size(); ++i) {
        work.step();
        const auto& row = columns[i];
        const bool has_upstream = i + 1 < columns.size() && columns[i + 1].key.band == row.key.band;
        if (has_upstream != row.leading_key.has_value() ||
            (has_upstream &&
             (columns[i + 1].key != *row.leading_key || columns[i + 1].id != row.leading_id)))
            reject(StaleGeneration);
    }
    return {std::move(traces), std::move(attached)};
}

} // namespace

core::Result<WakeStrengthTransition> prepare_wake_strength_transition(
    const io::SiCase& si, const WakeTransportClock& clock, const WakeStrengthEndpoint& parent,
    const WakeStrengthEndpoint& trial, a::ReleaseEpoch epoch, const WakeTransportLimits& limits,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    constexpr auto fallback = PrecisionUnresolved;
    return core::guard_entry(
        operation, fallback, {}, [&]() -> core::Result<WakeStrengthTransition> {
            try {
                auto material = read(prepare_wake_material_transition(
                    si, clock, {parent.source.contact, parent.material},
                    {trial.source.contact, trial.material}, epoch, limits, memory, stop));
                Work work{material.used_steps(), limits.maximum_steps, stop};
                const auto before = describe(parent, memory, work);
                auto after = describe(trial, memory, work);
                if (before.attached.size() != after.attached.size())
                    reject(StaleGeneration);
                for (std::size_t i = 0; i < before.attached.size(); ++i) {
                    work.step();
                    const auto& old = before.attached.values()[i];
                    const auto& fresh = after.attached.values()[i];
                    if (old.band != fresh.band || old.edge != fresh.edge)
                        reject(StaleGeneration);
                    const MaterialTraceKey newborn{old.band, epoch};
                    const auto row = find(after.frozen.values(), newborn, &Column::key, work);
                    if (row == after.frozen.values().end() || row->key != newborn ||
                        words(row->strength) != old.value)
                        reject(StaleGeneration);
                    row->parent = old.parent;
                }

                auto retired = allocate<std::optional<a::MaterialReleaseKey>>(
                    memory, material.retired_count());
                for (std::size_t i = 0; i < retired.size(); ++i) {
                    work.step();
                    retired.values()[i] = read(material.retired(i)).release;
                }
                const auto expired = [&](a::MaterialReleaseKey key) {
                    const auto at =
                        find(retired.values(), std::optional{key}, std::identity{}, work);
                    return at != retired.values().end() && *at == key;
                };
                const auto retired_column = [&](const Column& column) {
                    const auto& key = column.key;
                    return expired({key.band.from, key.epoch}) && expired({key.band.to, key.epoch});
                };
                std::size_t retained = 0;
                for (const auto& old : before.frozen.values()) {
                    work.step();
                    if (retired_column(old))
                        continue;
                    const auto row = find(after.frozen.values(), old.key, &Column::key, work);
                    if (row == after.frozen.values().end() || row->key != old.key ||
                        words(row->strength) != words(old.strength))
                        reject(StaleGeneration);
                    if (row->parent)
                        reject(StaleGeneration);
                    row->parent.emplace(old.id, false);
                    ++retained;
                }
                if (after.frozen.size() < retained ||
                    after.frozen.size() - retained != before.attached.size())
                    reject(StaleGeneration);

                const auto physical = [](const MaterialTraceKey& key) {
                    return WakeTraceMaterial{key.band.from, key.band.to, key.epoch};
                };
                auto retired_history =
                    allocate<RetiredWakeTrace>(memory, before.frozen.size() - retained);
                auto retired_at = retired_history.values().begin();
                for (const auto& old : before.frozen.values()) {
                    work.step();
                    if (retired_column(old))
                        *retired_at++ = {physical(old.key), {old.id, old.strength}};
                }
                auto columns = after.frozen.values();
                std::ranges::sort(columns, [&](const Column& lhs, const Column& rhs) {
                    work.step();
                    return lhs.id < rhs.id;
                });
                auto history = allocate<FrozenWakeTrace>(memory, columns.size());
                auto identities = allocate<WakeTraceMaterial>(memory, columns.size());
                auto parents = allocate<WakeStrengthParent>(memory, columns.size());
                for (std::size_t i = 0; i < columns.size(); ++i) {
                    work.step();
                    if (!columns[i].parent)
                        reject(StaleGeneration);
                    history.values()[i] = {columns[i].id, columns[i].strength};
                    identities.values()[i] = physical(columns[i].key);
                    parents.values()[i] = *columns[i].parent;
                }
                work.step();
                return WakeStrengthTransition{std::move(material), std::move(history),
                                              std::move(identities), std::move(retired_history),
                                              std::move(parents)};
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            } catch (const a::p1::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, operation});
            }
        });
}

core::Result<core::PlannedBuffer<FrozenWakeTrace>> apply_wake_strength_transfer(
    const WakeStrengthTransition& transition, std::span<const FrozenWakeTrace> parent_direction,
    std::uint64_t maximum_steps, core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return core::guard_entry(
        operation, PrecisionUnresolved, {},
        [&]() -> core::Result<core::PlannedBuffer<FrozenWakeTrace>> {
            try {
                Work work{0, maximum_steps, stop};
                work.step();
                for (std::size_t i = 0; i < parent_direction.size(); ++i) {
                    work.step();
                    if (i != 0 && parent_direction[i - 1].id >= parent_direction[i].id)
                        reject(InvalidInput);
                    (void)words(parent_direction[i].strength);
                }
                const auto columns = transition.history();
                const auto parents = transition.parents();
                auto output = allocate<FrozenWakeTrace>(memory, columns.size());
                for (std::size_t i = 0; i < columns.size(); ++i) {
                    work.step();
                    const auto& parent = parents[i];
                    const auto found =
                        find(parent_direction, parent.id, &FrozenWakeTrace::id, work);
                    a::ScalarWords variation{a::SingleWord{0.}};
                    if (found != parent_direction.end() && found->id == parent.id) {
                        const auto value = words(found->strength, parent.reversed ? -1. : 1.);
                        variation = std::holds_alternative<a::DoubleWord>(found->strength)
                                        ? a::ScalarWords{a::DoubleWord{value[0], value[1]}}
                                        : a::ScalarWords{a::SingleWord{value[0]}};
                    }
                    output.values()[i] = {columns[i].id, variation};
                }
                work.step();
                return output;
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            } catch (const a::p1::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, operation});
            }
        });
}

core::Result<FreeWakeStrengthCovector> WakeStrengthPullback::row(std::size_t i) const noexcept {
    const auto rows = rows_.values();
    if (i < rows.size()) {
        if (const auto& row = rows[i])
            return *row;
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}

core::Result<WakeStrengthPullback>
pullback_wake_strength_transfer(const WakeStrengthTransition& transition,
                                std::span<const FreeWakeStrengthCovector> seeds,
                                core::Precision format, std::uint64_t maximum_steps,
                                core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    namespace p = a::p1;
    return core::guard_entry(
        operation, PrecisionUnresolved, {}, [&]() -> core::Result<WakeStrengthPullback> {
            try {
                Work work{0, maximum_steps, stop};
                work.step();
                if (!a::p1_environment_supported() ||
                    (format != core::Precision::Binary64 && format != core::Precision::TwoWord))
                    reject(PrecisionUnsupported);
                const auto columns = transition.history();
                const auto parents = transition.parents();
                if (seeds.size() != columns.size())
                    reject(InvalidInput);
                struct Weight {
                    a::P1StrengthId parent;
                    std::size_t trial{};
                };
                auto ordered = allocate<Weight>(memory, columns.size());
                for (std::size_t i = 0; i < seeds.size(); ++i) {
                    work.step();
                    if (seeds[i].id != columns[i].id)
                        reject(StaleGeneration);
                    ordered.values()[i] = {parents[i].id, i};
                }
                const auto rows = ordered.values();
                std::ranges::sort(rows, [&](const Weight& lhs, const Weight& rhs) {
                    work.step();
                    return std::tie(lhs.parent, lhs.trial) < std::tie(rhs.parent, rhs.trial);
                });
                std::size_t count = 0;
                for (std::size_t i = 0; i < rows.size(); ++i) {
                    work.step();
                    if (i == 0 || rows[i - 1].parent != rows[i].parent)
                        ++count;
                }
                auto output = allocate<std::optional<FreeWakeStrengthCovector>>(memory, count);
                a::LocalRootScratch scratch;
                std::size_t next = 0;
                for (std::size_t i = 0; i < rows.size();) {
                    const auto id = rows[i].parent;
                    p::Interval sum{};
                    do {
                        work.step();
                        const auto j = rows[i++].trial;
                        const auto& source = seeds[j].value;
                        const auto center = words(source.center());
                        const auto value = p::add(a::decode_p1_word(center[0]),
                                                  a::decode_p1_word(center[1]), work);
                        const auto radius = a::decode_p1_word(source.radius());
                        const auto interval = p::scale({p::subtract(value, radius, work).lower,
                                                        p::add(value, radius, work).upper},
                                                       source.exponent());
                        sum = p::add(
                            sum, parents[j].reversed ? p::negate(interval, work) : interval, work);
                    } while (i < rows.size() && rows[i].parent == id);
                    output.values()[next++].emplace(id, a::enclose_p1_value(sum, format, scratch));
                }
                work.step();
                return WakeStrengthPullback{std::move(output)};
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            } catch (const p::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, operation});
            }
        });
}

core::Result<FreeWakeTrial> evaluate_free_wake_trial(
    const io::SiCase& si, const WakeTransportClock& clock, const FreeWakeTrialEndpoint& parent,
    const FreeWakeTrialEndpoint& trial, const FreeWakeTrialHistory& history, a::ReleaseEpoch epoch,
    std::span<const NonpenetrationVector> seeds, const FreeWakeTrialLimits& limits,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return core::guard_entry(
        operation, PrecisionUnresolved, {}, [&]() -> core::Result<FreeWakeTrial> {
            try {
                auto strength = read(prepare_wake_strength_transition(
                    si, clock, parent.wake, trial.wake, epoch, limits.material, memory, stop));
                const auto& links = strength.material();
                const auto count = links.size();
                Work dispatch{0, limits.transport.maximum_steps, stop};
                dispatch.step();
                if (count > limits.transport.maximum_particles)
                    reject(ResourceExceeded);
                if (!seeds.empty() && seeds.size() != count)
                    reject(InvalidInput);

                const bool restart = clock.destination.step == 1;
                std::size_t first_count = 0;
                for (std::size_t i = 0; i < count; ++i) {
                    dispatch.step();
                    if (restart || !read(links.row(i)).parent)
                        ++first_count;
                }
                const auto retained = count - first_count;
                const auto& previous_direction = history.retained_direction.previous;
                const auto& older_direction = history.retained_direction.older;
                const auto direction_count = trial.direction == nullptr ? 0 : retained;
                if (previous_direction.size() != direction_count ||
                    older_direction.size() != direction_count ||
                    (retained != 0 && history.older == nullptr) ||
                    (restart && history.older != nullptr) ||
                    (parent.direction != nullptr && trial.direction == nullptr) ||
                    (trial.direction != nullptr && !trial.direction->history.empty()))
                    reject(InvalidInput);

                const auto seed_bytes = seeds.empty() ? 0 : sizeof(NonpenetrationVector);
                const auto first_bytes = read(core::checked_byte_count(
                    static_cast<core::Index>(first_count),
                    sizeof(FirstFreeWakePoint) + 2 * sizeof(a::PointId) + seed_bytes, operation));
                const auto retained_bytes = read(core::checked_byte_count(
                    static_cast<core::Index>(retained),
                    3 * sizeof(WakeMaterialPoint) + sizeof(a::PointId) + seed_bytes, operation));
                if (retained_bytes > std::numeric_limits<std::size_t>::max() - first_bytes)
                    reject(ResourceExceeded);
                auto storage = allocate<std::byte>(memory, first_bytes + retained_bytes);
                std::pmr::monotonic_buffer_resource arena{storage.values().data(), storage.size(),
                                                          std::pmr::null_memory_resource()};
                const auto staging = [&]<class T>(std::in_place_type_t<T>, std::size_t size) {
                    static_assert(std::is_trivially_destructible_v<T>);
                    // Every array begins at the same alignment; the budget is exactly payload.
                    static_assert(alignof(T) == alignof(a::PointId));
                    return size == 0
                               ? std::span<T>{}
                               : std::span<T>{
                                     std::pmr::polymorphic_allocator<T>{&arena}.allocate(size),
                                     size};
                };
                const auto first_points =
                    staging(std::in_place_type<FirstFreeWakePoint>, first_count);
                const auto from = staging(std::in_place_type<a::PointId>, first_count);
                const auto at = staging(std::in_place_type<a::PointId>, first_count);
                const auto current = staging(std::in_place_type<WakeMaterialPoint>, retained);
                const auto previous = staging(std::in_place_type<WakeMaterialPoint>, retained);
                const auto older = staging(std::in_place_type<WakeMaterialPoint>, retained);
                const auto old_particles = staging(std::in_place_type<a::PointId>, retained);
                const auto first_seeds = staging(std::in_place_type<NonpenetrationVector>,
                                                 seeds.empty() ? 0 : first_count);
                const auto old_seeds =
                    staging(std::in_place_type<NonpenetrationVector>, seeds.empty() ? 0 : retained);
                const auto available_older = history.older == nullptr
                                                 ? std::span<const WakeMaterialPoint>{}
                                                 : history.older->points;
                for (std::size_t i = 1; i < available_older.size() && retained != 0; ++i) {
                    dispatch.step();
                    if (available_older[i - 1].release >= available_older[i].release)
                        reject(StaleGeneration);
                }
                std::size_t first_row = 0, retained_row = 0;
                for (std::size_t i = 0; i < count; ++i) {
                    dispatch.step();
                    const auto row = read(links.row(i));
                    if (restart || !row.parent) {
                        std::construct_at(first_points.data() + first_row, row);
                        std::construct_at(from.data() + first_row, row.parent
                                                                       ? row.parent->point
                                                                       : row.current.release.point);
                        std::construct_at(at.data() + first_row, row.current.point);
                        if (!seeds.empty())
                            std::construct_at(first_seeds.data() + first_row, seeds[i]);
                        ++first_row;
                        continue;
                    }
                    const auto key = row.current.release;
                    const auto found =
                        find(available_older, key, &WakeMaterialPoint::release, dispatch);
                    if (found == available_older.end() || found->release != key)
                        reject(StaleGeneration);
                    std::construct_at(current.data() + retained_row, row.current);
                    std::construct_at(previous.data() + retained_row, *row.parent);
                    std::construct_at(older.data() + retained_row, *found);
                    std::construct_at(old_particles.data() + retained_row, row.current.point);
                    if (!seeds.empty())
                        std::construct_at(old_seeds.data() + retained_row, seeds[i]);
                    ++retained_row;
                }
                core::PlannedBuffer<FrozenWakeTrace> trial_history_direction;
                std::optional<DiscreteFreeWakeDirection> current_direction;
                if (trial.direction != nullptr) {
                    core::PlannedBuffer<FrozenWakeTrace> parent_direction;
                    if (parent.direction != nullptr) {
                        const auto& before = parent.wake.source;
                        const auto old = parent.direction->history;
                        if (old.size() != before.history.columns.size())
                            reject(InvalidInput);
                        if (old.size() >
                            std::numeric_limits<std::size_t>::max() - before.collocations.size())
                            reject(ResourceExceeded);
                        parent_direction = allocate<FrozenWakeTrace>(
                            memory, old.size() + before.collocations.size());
                        auto values = parent_direction.values();
                        for (std::size_t i = 0; i < old.size(); ++i) {
                            dispatch.step();
                            if (old[i].id != before.history.columns[i].id)
                                reject(StaleGeneration);
                            values[i] = old[i];
                        }
                        for (std::size_t i = 0; i < before.collocations.size(); ++i) {
                            dispatch.step();
                            const auto& bound = before.collocations[i];
                            const auto face = panel(before.contact, bound.facet, dispatch).face;
                            values[old.size() + i] = {a::P1StrengthId::bound(face),
                                                      bound.circulation_direction};
                        }
                        sort_unique(values, &FrozenWakeTrace::id, dispatch);
                    }
                    trial_history_direction = read(
                        apply_wake_strength_transfer(strength, parent_direction.values(),
                                                     limits.material.maximum_steps, memory, stop));
                    current_direction = *trial.direction;
                    current_direction->history = trial_history_direction.values();
                }
                auto evaluation = limits.transport;
                evaluation.maximum_steps -= dispatch.used;
                const DiscreteNonpenetrationInput source{
                    trial.wake.source.contact,
                    trial.wake.source.collocations,
                    trial.wake.source.wake,
                    {trial.wake.source.history.parent, strength.history()},
                    trial.wake.source.trial_stamp,
                    trial.wake.source.flow_stamp,
                    trial.wake.source.flow_side};
                std::optional<WakeTransportResult> trapezoid, bdf2;
                if (first_count != 0) {
                    trapezoid.emplace(read(evaluate_first_free_wake_transport(
                        clock,
                        {source, trial.flow, at, current_direction ? &*current_direction : nullptr},
                        {parent.wake.source, parent.flow, from, parent.direction}, epoch,
                        first_points, first_seeds, evaluation, limits.velocity, shared, local,
                        kernel, memory, stop)));
                }
                if (retained != 0) {
                    const std::array<WakeMaterialSnapshot, 3> snapshots{
                        {{source.contact, current},
                         {parent.wake.source.contact, previous},
                         {history.older->contact, older, history.older->preceding_births}}};
                    bdf2.emplace(read(evaluate_fixed_free_wake_transport(
                        clock, snapshots, source, trial.flow, old_particles,
                        current_direction ? &*current_direction : nullptr,
                        history.retained_direction, old_seeds, evaluation, limits.velocity, shared,
                        local, kernel, memory, stop)));
                }
                dispatch.step();
                return FreeWakeTrial{std::move(strength), std::move(trapezoid), std::move(bdf2)};
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            } catch (const a::p1::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, operation});
            }
        });
}

} // namespace uvlas::coupling::detail
