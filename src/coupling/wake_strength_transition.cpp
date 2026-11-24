#include "wake_strength_transition.hpp"

#include "aero/p1_interval_private.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

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
    std::array<double, 2> value{};
    std::optional<MaterialTraceKey> leading_key;
    std::optional<a::P1StrengthId> leading_id;
};

struct AttachedColumn {
    a::EdgeId edge{0};
    WakeSpanBand band;
    std::array<double, 2> value{};
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
        trace = {back, row.trailing, words(at->strength), {}, {}};
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
            attached.values()[front++] = {transfer.edge, band,
                                          words(gamma->gamma, transfer.coefficient)};
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
                const auto after = describe(trial, memory, work);
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
                        row->value != old.value)
                        reject(StaleGeneration);
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
                std::size_t retained = 0;
                for (const auto& old : before.frozen.values()) {
                    work.step();
                    const auto& key = old.key;
                    if (expired({key.band.from, key.epoch}) && expired({key.band.to, key.epoch}))
                        continue;
                    const auto row = find(after.frozen.values(), key, &Column::key, work);
                    if (row == after.frozen.values().end() || row->key != key ||
                        row->value != old.value)
                        reject(StaleGeneration);
                    ++retained;
                }
                if (after.frozen.size() < retained ||
                    after.frozen.size() - retained != before.attached.size())
                    reject(StaleGeneration);

                const auto supplied = trial.source.history.columns;
                auto history = allocate<FrozenWakeTrace>(memory, supplied.size());
                for (std::size_t i = 0; i < supplied.size(); ++i) {
                    work.step();
                    history.values()[i] = supplied[i];
                }
                work.step();
                return WakeStrengthTransition{std::move(material), std::move(history)};
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            } catch (const a::p1::ArithmeticFailure& failure) {
                return std::unexpected(core::Failure{failure.code, operation});
            }
        });
}

} // namespace uvlas::coupling::detail
