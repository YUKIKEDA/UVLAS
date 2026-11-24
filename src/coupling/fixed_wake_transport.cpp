#include "fixed_wake_transport.hpp"

#include "aero/p1_interval_private.hpp"
#include "case_component_precision_private.hpp"
#include "io/si_case_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace m = a::p1;
namespace p = a::pchip;

namespace {
using enum core::FailureCode;
constexpr auto operation = core::Operation::EvaluateAero;
using Vector = std::array<m::Interval, 3>;

[[noreturn]] void reject(core::FailureCode reason) {
    throw m::ArithmeticFailure{reason};
}
template <class Result>
auto read(Result result) {
    if (!result) {
        const auto& error = result.error();
        if constexpr (std::is_same_v<typename Result::error_type, core::Failure>)
            throw error;
        else if constexpr (std::is_same_v<typename Result::error_type, a::GlobalGeometryFailure>)
            throw error.origin.cause;
        else
            throw error.cause;
    }
    return std::move(*result);
}
template <class F>
auto guarded(F action) noexcept {
    using Result = core::Result<decltype(action())>;
    return core::guard_entry(operation, PrecisionUnresolved, {}, [&]() -> Result {
        try {
            return action();
        } catch (const core::Failure& failure) {
            return std::unexpected(failure);
        } catch (const m::ArithmeticFailure& failure) {
            return std::unexpected(core::Failure{failure.code, operation});
        }
    });
}
template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& memory, std::size_t size) {
    if (size > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        reject(ResourceExceeded);
    return read(memory.allocate<T>(static_cast<core::Index>(size), operation));
}

m::Interval lift(double word) {
    if (!std::isfinite(word))
        reject(Nonfinite);
    return a::decode_p1_word(word);
}
m::Interval lift(const a::ScalarWords& words, m::Work& work) {
    work.step();
    if (const auto* single = std::get_if<a::SingleWord>(&words))
        return lift(single->value);
    const auto& paired = std::get<a::DoubleWord>(words);
    return m::add(lift(paired.high), lift(paired.low), work);
}
m::Interval lift(const a::ScaledEnclosure& words, m::Work& work) {
    const auto value = lift(words.center(), work);
    const auto width = lift(words.radius());
    return m::scale({m::subtract(value, width, work).lower, m::add(value, width, work).upper},
                    words.exponent());
}
m::Interval lift(const io::SiEnclosure& words, m::Work& work) {
    return std::visit(
        [&](const auto& value) {
            using Word = std::decay_t<decltype(value)>;
            const auto center = [&]() -> m::Interval {
                if constexpr (std::is_same_v<Word, io::ScaledBinary64>)
                    return lift(value.center);
                else
                    return m::add(lift(value.high), lift(value.low), work);
            }();
            const auto width = lift(value.radius);
            return m::scale(
                {m::subtract(center, width, work).lower, m::add(center, width, work).upper},
                m::exponent(value.exponent2));
        },
        words);
}
Vector lift(const NonpenetrationVector& v, m::Work& work) {
    return {lift(v[0], work), lift(v[1], work), lift(v[2], work)};
}
Vector lift(const DiscreteLoadVector& v, m::Work& work) {
    return {lift(v[0], work), lift(v[1], work), lift(v[2], work)};
}

p::RationalView view(const io::detail::SiRationalWords& words) {
    const auto parts = words.view();
    return {parts.sign, parts.numerator, parts.denominator};
}
bool same_time(p::RationalView x, p::RationalView y, m::Work& work) {
    work.step();
    const auto equal = [&](auto a, auto b) {
        work.step();
        return a == b;
    };
    return x.sign == y.sign && std::ranges::equal(x.numerator, y.numerator, equal) &&
           std::ranges::equal(x.denominator, y.denominator, equal);
}
bool same_model(const a::GeometryContext& x, const a::GeometryContext& y) {
    return x.fingerprint == y.fingerprint && x.model == y.model && x.reference == y.reference &&
           x.mapping == y.mapping && x.kernel == y.kernel && x.policy == y.policy &&
           x.radius == y.radius;
}

struct RetainedParticle {
    a::MaterialReleaseKey release{a::PointId{0}, a::ReleaseEpoch{0}};
    std::array<a::PointId, 3> point{a::PointId{0}, a::PointId{0}, a::PointId{0}};
    std::array<NonpenetrationVector, 3> positions{
        NonpenetrationVector{a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}},
        NonpenetrationVector{a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}},
        NonpenetrationVector{a::SingleWord{0}, a::SingleWord{0}, a::SingleWord{0}}};
};
struct TransportStencil {
    core::PlannedBuffer<RetainedParticle> rows;
    std::array<io::detail::SiRationalWords, 3> times;
    m::Interval twice_step{};
};

bool trailing_point(const a::PreparedTopologyOwner& topology, a::PointId point, m::Work& work) {
    // An upper bound for searching the sorted table indexed by a 64-bit ID.
    for (int i = 0; i <= std::numeric_limits<std::uint64_t>::digits; ++i)
        work.step();
    return read(topology.is_trailing_point(point));
}

a::PointId certify_particle(const WakeMaterialPoint& point, const a::OwnedGeometry& geometry,
                            const a::GlobalContactState& contact,
                            const a::PreparedTopologyOwner& topology, m::Work& work) {
    const auto canonical = read(geometry.owner(point.point));
    if (read(geometry.owner(point.release.point)) != point.release.point ||
        !trailing_point(topology, point.release.point, work) ||
        trailing_point(topology, canonical, work))
        reject(InvalidInput);
    if (!point.wake_facet)
        reject(InvalidInput);
    const auto wake = read(contact.uvlm_quad(*point.wake_facet, work.limit - work.used, work.stop));
    work.used += wake.used_steps;
    if (wake.role != a::AeroFaceRole::Wake)
        reject(InvalidInput);

    for (const auto facet : wake.facets)
        for (const auto id : read(geometry.facet_points(facet))) {
            work.step();
            if (read(geometry.owner(id)) == canonical)
                return canonical;
        }
    reject(InvalidInput);
}

TransportStencil prepare_stencil(const WakeTransportClock& clock,
                                 const std::array<WakeMaterialSnapshot, 3>& snapshots,
                                 const WakeTransportLimits& limits, core::ResourcePlan& memory,
                                 m::Work& work, bool allow_births = false) {
    work.step();
    const auto count = snapshots.front().points.size();
    if (count > limits.maximum_particles)
        reject(ResourceExceeded);
    if (count == 0 || clock.destination.step < 2)
        reject(InvalidInput);
    if (!a::p1_environment_supported() || (limits.representation != core::Precision::Binary64 &&
                                           limits.representation != core::Precision::TwoWord))
        reject(PrecisionUnsupported);
    auto retained = allocate<RetainedParticle>(memory, count);
    auto visited = allocate<std::uint64_t>(memory, count);
    std::array<io::detail::SiRationalWords, 3> times;
    const auto& base = read(snapshots.front().contact.geometry()).get();
    for (std::size_t k = 0; k < snapshots.size(); ++k) {
        work.step();
        const auto& snapshot = snapshots[k];
        if (snapshot.points.size() != count ||
            (snapshot.preceding_births && (!allow_births || k != 2)))
            reject(InvalidInput);
        auto address = clock.destination;
        if (k != 0) {
            address.step -= k;
            address.side =
                address.step == 0 ? p::RequestedSide::Right : p::RequestedSide::Automatic;
        }
        const auto instant = read(evaluate_case_reference_clock(clock.clock, clock.current, address,
                                                                limits.clock, memory, work.stop));
        times[k] = read(io::detail::copy_rational_words(instant.seconds, memory, work.stop));
        const auto& geometry = read(snapshot.contact.geometry()).get();
        const auto exact = read(geometry.exact_time());
        if (!exact || !same_model(base.context(), geometry.context()) ||
            !same_time(view(times[k]), *exact, work))
            reject(StaleGeneration);
        const auto topology = read(snapshot.contact.reference_topology());
        for (std::size_t i = 0; i < count; ++i) {
            work.step();
            const auto& point = snapshot.points[i];
            if ((i != 0 && snapshot.points[i - 1].release >= point.release) ||
                point.release != snapshots.front().points[i].release)
                reject(StaleGeneration);
            const auto canonical = [&] {
                auto birth = WakeBirthMatch::Missing;
                if (snapshot.preceding_births) {
                    const auto& detached = read(snapshots[1].contact.geometry()).get();
                    for (int bit = 0; bit <= std::numeric_limits<std::size_t>::digits; ++bit)
                        work.step();
                    birth = snapshot.preceding_births->match_birth(
                        point.release, clock, geometry, detached, retained.values()[i].point[1]);
                }
                if (birth == WakeBirthMatch::Stale ||
                    (birth == WakeBirthMatch::Verified && point.wake_facet))
                    reject(StaleGeneration);
                if (point.wake_facet)
                    return certify_particle(point, geometry, snapshot.contact, topology, work);
                if (birth != WakeBirthMatch::Verified || point.point != point.release.point ||
                    read(geometry.owner(point.point)) != point.point ||
                    !trailing_point(topology, point.point, work))
                    reject(InvalidInput);
                return point.point;
            }();
            const auto row = retained.values();
            row[i].release = point.release;
            row[i].point[k] = canonical;
            row[i].positions[k] = read(geometry.point_words(canonical, a::Configuration::Current));
            visited.values()[i] = canonical.value();
        }
        auto ordered = visited.values();
        std::ranges::sort(ordered, [&](auto x, auto y) {
            work.step();
            return x < y;
        });
        for (std::size_t i = 1; i < count; ++i) {
            work.step();
            if (ordered[i - 1] == ordered[i])
                reject(InvalidInput);
        }
    }
    const auto ordinal = static_cast<std::size_t>(clock.destination.interval);
    const auto& interval = clock.clock.intervals.values()[ordinal];
    const auto h = read(io::detail::copy_rational_words(interval.step, memory, work.stop));
    const auto step =
        m::scale(case_components::rational(view(h), limits.maximum_integer_bits, work), 1);
    if (m::compare(step.lower, {}) <= 0)
        reject(InvalidInput);
    return {std::move(retained), std::move(times), step};
}

Vector transport_rate(const Vector& x, const Vector& previous, const Vector& older,
                      const m::Interval& step, m::Work& work) {
    Vector result;
    for (std::size_t i = 0; i < result.size(); ++i) {
        const auto first = m::subtract(x[i], previous[i], work);
        const auto second = m::subtract(previous[i], older[i], work);
        result[i] = m::divide(m::subtract(m::multiply(m::integer(3), first, work), second, work),
                              step, work);
    }
    return result;
}
Vector difference(const Vector& lhs, const Vector& rhs, m::Work& work) {
    return {m::subtract(lhs[0], rhs[0], work), m::subtract(lhs[1], rhs[1], work),
            m::subtract(lhs[2], rhs[2], work)};
}
DiscreteLoadVector keep(const Vector& v, core::Precision format, a::LocalRootScratch& scratch,
                        m::Work& work) {
    const auto word = [&](std::size_t k) {
        work.step();
        return a::enclose_p1_value(v[k], format, scratch);
    };
    return {word(0), word(1), word(2)};
}
Vector prescribed_position(const Vector& previous, const Vector& older, const Vector& velocity,
                           const m::Interval& twice_step, m::Work& work) {
    Vector result;
    for (std::size_t k = 0; k < result.size(); ++k) {
        const auto move = m::add(m::subtract(previous[k], older[k], work),
                                 m::multiply(twice_step, velocity[k], work), work);
        result[k] = m::add(previous[k], m::divide(move, m::integer(3), work), work);
    }
    return result;
}
WakeTransportPositionDual stencil_dual(a::MaterialReleaseKey release,
                                       const NonpenetrationVector& seed,
                                       const m::Interval& twice_step, core::Precision format,
                                       a::LocalRootScratch& scratch, m::Work& work) {
    const auto input = lift(seed, work);
    const auto project = [&](int c) {
        Vector scaled;
        for (std::size_t k = 0; k < scaled.size(); ++k)
            scaled[k] = m::divide(m::multiply(m::integer(c), input[k], work), twice_step, work);
        return keep(scaled, format, scratch, work);
    };
    return {release, {project(3), project(-4), project(1)}};
}

Vector current_direction(const DiscreteFreeWakeDirection& direction, a::PointId owner,
                         m::Work& work) {
    const auto found = std::ranges::lower_bound(direction.points, owner, [&](auto lhs, auto rhs) {
        work.step();
        return lhs < rhs;
    });
    if (found == direction.points.end() || *found != owner)
        reject(InvalidInput);
    return lift(direction.displacements[static_cast<std::size_t>(found - direction.points.begin())],
                work);
}
} // namespace

struct WakeBirthHistory {
    io::SchemaIdentity source;
    p::PchipInputStamp stamp;
    a::GeometryContext origin_context, detached_context;
    std::array<std::byte, 32> origin_geometry, detached_geometry;
    core::PlannedBuffer<std::optional<a::PointId>> detached_points;
};

struct WakeTransportData {
    core::PlannedBuffer<std::optional<WakeTransportRow>> rows;
    core::PlannedBuffer<std::optional<WakeTransportPositionDual>> dual;
    std::optional<DiscreteFreeWakeVelocity> free;
    std::optional<DiscreteFreeWakeVelocity> parent_free;
    std::optional<WakeBirthHistory> birth_history;
    std::optional<a::ScaledEnclosure> prescribed_time;
};

WakeTransportResult::WakeTransportResult(core::PlannedBuffer<WakeTransportData>&& data) noexcept
    : data_(std::move(data)) {}
WakeTransportResult::WakeTransportResult(WakeTransportResult&&) noexcept = default;
WakeTransportResult& WakeTransportResult::operator=(WakeTransportResult&&) noexcept = default;
WakeTransportResult::~WakeTransportResult() = default;

std::size_t WakeTransportResult::size() const noexcept {
    return data_.size() == 0 ? 0 : data_.values().front().rows.size();
}
core::Result<WakeTransportRow> WakeTransportResult::row(std::size_t i) const noexcept {
    if (i < size()) {
        const auto rows = data_.values().front().rows.values();
        if (const auto& entry = rows[i])
            return *entry;
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}
WakeBirthMatch WakeTransportResult::match_birth(a::MaterialReleaseKey key,
                                                const WakeTransportClock& clock,
                                                const a::OwnedGeometry& origin,
                                                const a::OwnedGeometry& detached,
                                                a::PointId detached_point) const noexcept {
    if (data_.size() == 0)
        return WakeBirthMatch::Missing;
    const auto& data = data_.values().front();
    const auto& birth = data.birth_history;
    if (!birth)
        return WakeBirthMatch::Missing;

    const auto rows = data.rows.values();
    const auto at = std::ranges::lower_bound(
        rows, std::optional{key}, {}, [](const auto& row) -> std::optional<a::MaterialReleaseKey> {
            if (row)
                return row->release;
            return std::nullopt;
        });
    if (at == rows.end() || !*at || (*at)->release != key)
        return WakeBirthMatch::Missing;
    const auto& retained =
        birth->detached_points.values()[static_cast<std::size_t>(at - rows.begin())];
    if (!retained)
        return WakeBirthMatch::Missing;

    if (*retained != detached_point ||
        !io::detail::same_identity(birth->source, clock.clock.source) ||
        birth->stamp != clock.destination.stamp || birth->origin_context != origin.context() ||
        birth->detached_context != detached.context())
        return WakeBirthMatch::Stale;
    const auto origin_identity = origin.identity();
    const auto detached_identity = detached.identity();
    if (!origin_identity || !detached_identity || birth->origin_geometry != *origin_identity ||
        birth->detached_geometry != *detached_identity)
        return WakeBirthMatch::Stale;
    return WakeBirthMatch::Verified;
}

core::Result<WakeTransportPositionDual>
WakeTransportResult::position_dual(std::size_t i) const noexcept {
    if (data_.size() != 0) {
        const auto dual = data_.values().front().dual.values();
        if (i < dual.size()) {
            if (const auto& entry = dual[i])
                return *entry;
        }
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}
core::Result<std::reference_wrapper<const DiscreteFreeWakeVelocity>>
WakeTransportResult::free_velocity() const& noexcept {
    const auto owner = data_.values();
    if (!owner.empty()) {
        if (const auto& velocity = owner.front().free)
            return std::cref(*velocity);
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}

core::Result<std::reference_wrapper<const DiscreteFreeWakeVelocity>>
WakeTransportResult::parent_free_velocity() const& noexcept {
    const auto owner = data_.values();
    if (!owner.empty()) {
        if (const auto& velocity = owner.front().parent_free)
            return std::cref(*velocity);
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}

core::Result<a::ScaledEnclosure> WakeTransportResult::prescribed_time_dual() const noexcept {
    const auto owner = data_.values();
    if (!owner.empty()) {
        if (const auto& value = owner.front().prescribed_time)
            return *value;
    }
    return std::unexpected(core::Failure{InvalidInput, operation});
}

struct WakeTransportBuilder {
    static const a::OwnedGeometry& first_geometry(const WakeTransportClock& clock,
                                                  CaseFlowClockQuery address,
                                                  const FirstFreeWakeEndpoint& endpoint,
                                                  const WakeTransportLimits& limits,
                                                  core::ResourcePlan& memory, m::Work& work) {
        const auto instant = read(evaluate_case_reference_clock(clock.clock, clock.current, address,
                                                                limits.clock, memory, work.stop));
        const auto exact =
            read(io::detail::copy_rational_words(instant.seconds, memory, work.stop));
        const auto& geometry = read(endpoint.source.contact.geometry()).get();
        const auto geometry_time = read(geometry.exact_time());
        if (!geometry_time || !same_time(view(exact), *geometry_time, work) ||
            !same_time(view(exact), endpoint.flow.time.time, work) ||
            endpoint.flow.time.stamp != address.stamp || endpoint.flow.time.side != address.side ||
            !io::detail::same_identity(read(endpoint.flow.source.free_wake_case_identity()),
                                       clock.clock.source))
            reject(StaleGeneration);
        return geometry;
    }

    static Vector trapezoidal(const Vector& current, const Vector& parent,
                              const Vector& current_velocity, const Vector& parent_velocity,
                              const m::Interval& h, m::Work& work) {
        Vector result;
        for (std::size_t k = 0; k < result.size(); ++k) {
            const auto rate = m::divide(m::subtract(current[k], parent[k], work), h, work);
            const auto mean = m::scale(m::add(current_velocity[k], parent_velocity[k], work), -1);
            result[k] = m::subtract(rate, mean, work);
        }
        return result;
    }

    static WakeTransportResult
    first_free(const WakeTransportClock& clock, const FirstFreeWakeEndpoint& current,
               const FirstFreeWakeEndpoint& parent, std::optional<a::ReleaseEpoch> birth_epoch,
               std::span<const FirstFreeWakePoint> points,
               std::span<const NonpenetrationVector> seeds, const WakeTransportLimits& limits,
               const DiscreteFreeWakeLimits& velocity_limits, a::SharedP1Workspace& shared,
               a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
               const std::stop_token& stop) {
        m::Work work{0, limits.maximum_steps, stop};
        work.step();
        if (points.size() > limits.maximum_particles)
            reject(ResourceExceeded);
        if (points.empty() || current.particles.size() != points.size() ||
            parent.particles.size() != points.size() ||
            (!seeds.empty() && seeds.size() != points.size()) || clock.destination.step == 0)
            reject(InvalidInput);
        if (!a::p1_environment_supported() || (limits.representation != core::Precision::Binary64 &&
                                               limits.representation != core::Precision::TwoWord))
            reject(PrecisionUnsupported);
        auto parent_address = clock.destination;
        --parent_address.step;
        parent_address.side =
            parent_address.step == 0 ? p::RequestedSide::Right : p::RequestedSide::Automatic;
        const auto& x = first_geometry(clock, clock.destination, current, limits, memory, work);
        const auto& y = first_geometry(clock, parent_address, parent, limits, memory, work);
        if (!same_model(x.context(), y.context()))
            reject(StaleGeneration);

        const auto topology = read(current.source.contact.reference_topology());
        const auto parent_topology = read(parent.source.contact.reference_topology());
        const auto& interval =
            clock.clock.intervals.values()[static_cast<std::size_t>(clock.destination.interval)];
        const auto step = read(io::detail::copy_rational_words(interval.step, memory, stop));
        const auto h = case_components::rational(view(step), limits.maximum_integer_bits, work);
        if (m::compare(h.lower, {}) <= 0)
            reject(InvalidInput);

        for (std::size_t i = 0; i < points.size(); ++i) {
            work.step();
            const auto& point = points[i];
            const auto key = point.current.release;
            if (i != 0 && points[i - 1].current.release >= key)
                reject(InvalidInput);
            if (read(x.owner(current.particles[i])) !=
                certify_particle(point.current, x, current.source.contact, topology, work))
                reject(StaleGeneration);
            const auto parent_point = read(y.owner(parent.particles[i]));
            if (point.parent) {
                if (clock.destination.step != 1)
                    reject(InvalidInput);
                if (point.parent->release != key || (birth_epoch && key.epoch == *birth_epoch))
                    reject(StaleGeneration);
                if (parent_point != certify_particle(*point.parent, y, parent.source.contact,
                                                     parent_topology, work))
                    reject(StaleGeneration);
            } else {
                if (!birth_epoch || key.epoch != *birth_epoch || parent_point != key.point ||
                    read(y.owner(key.point)) != key.point ||
                    !trailing_point(parent_topology, key.point, work))
                    reject(StaleGeneration);
            }
        }

        const auto velocity = [&](const FirstFreeWakeEndpoint& endpoint) {
            return read(evaluate_discrete_free_wake_velocity(
                endpoint.source, endpoint.flow, endpoint.particles, endpoint.direction, seeds,
                velocity_limits, shared, local, kernel, memory, stop,
                FreeWakeSeedScale::NegativeHalf));
        };
        auto trial = velocity(current);
        auto accepted = velocity(parent);

        auto owner = allocate<WakeTransportData>(memory, 1);
        auto& out = owner.values().front();
        out.rows = allocate<std::optional<WakeTransportRow>>(memory, points.size());
        out.dual = allocate<std::optional<WakeTransportPositionDual>>(memory, seeds.size());
        if (birth_epoch)
            out.birth_history.emplace(clock.clock.source, clock.destination.stamp, y.context(),
                                      x.context(), read(y.identity()), read(x.identity()),
                                      allocate<std::optional<a::PointId>>(memory, points.size()));
        a::LocalRootScratch scratch;
        for (std::size_t i = 0; i < points.size(); ++i) {
            work.step();
            const auto now = read(trial.row(i));
            const auto last = read(accepted.row(i));
            const auto point = now.owner;
            const auto parent_point = last.owner;
            const auto residual = trapezoidal(
                lift(read(x.point_words(point, a::Configuration::Current)), work),
                lift(read(y.point_words(parent_point, a::Configuration::Current)), work),
                lift(now.velocity, work), lift(last.velocity, work), h, work);
            std::optional<DiscreteLoadVector> tangent;
            if (current.direction || parent.direction) {
                if ((current.direction != nullptr) != now.direction.has_value() ||
                    (parent.direction != nullptr) != last.direction.has_value())
                    reject(InvalidInput);
                const Vector dx = current.direction
                                      ? current_direction(*current.direction, point, work)
                                      : Vector{};
                const Vector dy = parent.direction
                                      ? current_direction(*parent.direction, parent_point, work)
                                      : Vector{};
                const Vector dv = now.direction ? lift(*now.direction, work) : Vector{};
                const Vector dp = last.direction ? lift(*last.direction, work) : Vector{};
                tangent.emplace(keep(trapezoidal(dx, dy, dv, dp, h, work), limits.representation,
                                     scratch, work));
            }
            const auto key = points[i].current.release;
            if (out.birth_history && !points[i].parent)
                out.birth_history->detached_points.values()[i] = point;
            out.rows.values()[i].emplace(WakeTransportRow{
                key, keep(residual, limits.representation, scratch, work), tangent, {}, {}});
            if (!seeds.empty()) {
                const auto seed = lift(seeds[i], work);
                Vector current_dual, parent_dual;
                for (std::size_t k = 0; k < seed.size(); ++k) {
                    current_dual[k] = m::divide(seed[k], h, work);
                    parent_dual[k] = m::subtract({}, current_dual[k], work);
                }
                out.dual.values()[i].emplace(WakeTransportPositionDual{
                    key,
                    {keep(current_dual, limits.representation, scratch, work),
                     keep(parent_dual, limits.representation, scratch, work),
                     keep({}, limits.representation, scratch, work)}});
            }
        }
        work.step();
        out.free.emplace(std::move(trial));
        out.parent_free.emplace(std::move(accepted));
        return WakeTransportResult{std::move(owner)};
    }

    static WakeTransportResult
    free(const WakeTransportClock& clock, const std::array<WakeMaterialSnapshot, 3>& snapshots,
         const DiscreteNonpenetrationInput& source, const DiscreteKnownFlowInput& flow,
         std::span<const a::PointId> particle, const DiscreteFreeWakeDirection* direction,
         const WakeTransportHistoryDirection& history, std::span<const NonpenetrationVector> seeds,
         const WakeTransportLimits& limits, const DiscreteFreeWakeLimits& velocity_limits,
         a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
         core::ResourcePlan& memory, const std::stop_token& stop) {
        m::Work work{0, limits.maximum_steps, stop};
        auto stencil = prepare_stencil(clock, snapshots, limits, memory, work, true);
        const auto size = stencil.rows.size();
        if ((direction == nullptr
                 ? !history.previous.empty() || !history.older.empty()
                 : history.previous.size() != size || history.older.size() != size) ||
            (!seeds.empty() && seeds.size() != size))
            reject(InvalidInput);
        const auto& x = read(snapshots.front().contact.geometry()).get();
        const auto& y = read(source.contact.geometry()).get();
        if (read(x.identity()) != read(y.identity()) || x.context() != y.context() ||
            !read(snapshots.front().contact.reference_topology())
                 .same_owner(read(source.contact.reference_topology())) ||
            !io::detail::same_identity(read(flow.source.free_wake_case_identity()),
                                       clock.clock.source) ||
            !same_time(view(stencil.times.front()), flow.time.time, work) ||
            flow.time.stamp != clock.destination.stamp || flow.time.side != clock.destination.side)
            reject(StaleGeneration);
        if (particle.size() != size)
            reject(InvalidInput);
        for (std::size_t i = 0; i < size; ++i) {
            work.step();
            if (read(y.owner(particle[i])) != stencil.rows.values()[i].point[0])
                reject(StaleGeneration);
        }
        auto samples = read(evaluate_discrete_free_wake_velocity(
            source, flow, particle, direction, seeds, velocity_limits, shared, local, kernel,
            memory, stop, FreeWakeSeedScale::Negative));
        auto owner = allocate<WakeTransportData>(memory, 1);
        const auto data = owner.values();
        data.front().rows = allocate<std::optional<WakeTransportRow>>(memory, size);
        data.front().dual =
            allocate<std::optional<WakeTransportPositionDual>>(memory, seeds.size());
        finish_free(data.front(), stencil, samples, direction, history, seeds,
                    limits.representation, work);
        data.front().free.emplace(std::move(samples));
        return WakeTransportResult{std::move(owner)};
    }

    static void finish_free(WakeTransportData& out, const TransportStencil& stencil,
                            const DiscreteFreeWakeVelocity& velocity,
                            const DiscreteFreeWakeDirection* direction,
                            const WakeTransportHistoryDirection& history,
                            std::span<const NonpenetrationVector> seeds, core::Precision format,
                            m::Work& work) {
        a::LocalRootScratch scratch;
        const auto rows = stencil.rows.values();
        for (std::size_t i = 0; i < rows.size(); ++i) {
            work.step();
            const auto& row = rows[i];
            const auto sample = read(velocity.row(i));
            if (sample.owner != row.point[0])
                reject(StaleGeneration);
            const auto lhs =
                transport_rate(lift(row.positions[0], work), lift(row.positions[1], work),
                               lift(row.positions[2], work), stencil.twice_step, work);
            const auto rhs = lift(sample.velocity, work);
            std::optional<DiscreteLoadVector> tangent;
            if (direction != nullptr) {
                if (!sample.direction)
                    reject(InvalidInput);
                const auto dx =
                    transport_rate(current_direction(*direction, row.point[0], work),
                                   lift(history.previous[i], work), lift(history.older[i], work),
                                   stencil.twice_step, work);
                tangent.emplace(keep(difference(dx, lift(*sample.direction, work), work), format,
                                     scratch, work));
            }
            out.rows.values()[i].emplace(
                WakeTransportRow{row.release,
                                 keep(difference(lhs, rhs, work), format, scratch, work),
                                 tangent,
                                 {},
                                 {}});
            if (!seeds.empty())
                out.dual.values()[i].emplace(
                    stencil_dual(row.release, seeds[i], stencil.twice_step, format, scratch, work));
        }
        work.step();
    }

    static WakeTransportResult prescribed(const WakeTransportClock& clock,
                                          const std::array<WakeMaterialSnapshot, 3>& snapshots,
                                          const PrescribedWakeClockSample& supplied,
                                          const PrescribedWakeTransportDirection* direction,
                                          std::span<const NonpenetrationVector> seeds,
                                          const WakeTransportLimits& limits,
                                          core::ResourcePlan& memory, const std::stop_token& stop) {
        m::Work work{0, limits.maximum_steps, stop};
        auto stencil = prepare_stencil(clock, snapshots, limits, memory, work);
        const auto size = stencil.rows.size();
        if ((direction != nullptr &&
             (direction->current.size() != size || direction->history.previous.size() != size ||
              direction->history.older.size() != size)) ||
            (!seeds.empty() && seeds.size() != size))
            reject(InvalidInput);
        const auto& field = std::visit(
            [&](const auto& owner) -> const ProjectedCasePrescribedWake& {
                const auto& sample = owner.get();
                const auto address = sample.address();
                if (!io::detail::same_identity(sample.source(), clock.clock.source) ||
                    address.stamp != clock.destination.stamp ||
                    address.interval != clock.destination.interval ||
                    address.step != clock.destination.step ||
                    address.side != clock.destination.side ||
                    !same_time(read(sample.time()), view(stencil.times.front()), work))
                    reject(StaleGeneration);
                return sample.prescribed();
            },
            supplied);
        auto owner = allocate<WakeTransportData>(memory, 1);
        const auto data = owner.values();
        data.front().rows = allocate<std::optional<WakeTransportRow>>(memory, size);
        data.front().dual =
            allocate<std::optional<WakeTransportPositionDual>>(memory, seeds.size());
        finish_prescribed(data.front(), stencil, field, direction, seeds, limits.representation,
                          work);
        return WakeTransportResult{std::move(owner)};
    }

    static void finish_prescribed(WakeTransportData& out, const TransportStencil& stencil,
                                  const ProjectedCasePrescribedWake& field,
                                  const PrescribedWakeTransportDirection* direction,
                                  std::span<const NonpenetrationVector> seeds,
                                  core::Precision format, m::Work& work) {
        Vector velocity, dt, dv;
        for (std::size_t k = 0; k < velocity.size(); ++k) {
            work.step();
            velocity[k] = lift(read(field.velocity(k)), work);
            if (direction != nullptr || !seeds.empty())
                dt[k] = lift(read(field.time_partial(k)), work);
            if (direction != nullptr)
                dv[k] = m::multiply(dt[k], lift(direction->flow_time, work), work);
        }
        a::LocalRootScratch scratch;
        m::Interval time_dual{};
        const auto rows = stencil.rows.values();
        for (std::size_t i = 0; i < rows.size(); ++i) {
            work.step();
            const auto& row = rows[i];
            const auto previous = lift(row.positions[1], work);
            const auto older = lift(row.positions[2], work);
            const auto lhs = transport_rate(lift(row.positions[0], work), previous, older,
                                            stencil.twice_step, work);
            const auto candidate =
                prescribed_position(previous, older, velocity, stencil.twice_step, work);
            std::optional<DiscreteLoadVector> tangent, candidate_direction;
            if (direction != nullptr) {
                const auto dp = lift(direction->history.previous[i], work);
                const auto old = lift(direction->history.older[i], work);
                const auto dx = transport_rate(lift(direction->current[i], work), dp, old,
                                               stencil.twice_step, work);
                tangent.emplace(keep(difference(dx, dv, work), format, scratch, work));
                candidate_direction.emplace(
                    keep(prescribed_position(dp, old, dv, stencil.twice_step, work), format,
                         scratch, work));
            }
            out.rows.values()[i].emplace(WakeTransportRow{
                row.release, keep(difference(lhs, velocity, work), format, scratch, work), tangent,
                keep(candidate, format, scratch, work), candidate_direction});
            if (!seeds.empty()) {
                out.dual.values()[i].emplace(
                    stencil_dual(row.release, seeds[i], stencil.twice_step, format, scratch, work));
                const auto weight = lift(seeds[i], work);
                for (std::size_t k = 0; k < weight.size(); ++k)
                    time_dual = m::subtract(time_dual, m::multiply(weight[k], dt[k], work), work);
            }
        }
        work.step();
        if (!seeds.empty())
            out.prescribed_time.emplace(a::enclose_p1_value(time_dual, format, scratch));
    }
};

core::Result<WakeTransportResult> evaluate_fixed_free_wake_transport(
    const WakeTransportClock& clock, const std::array<WakeMaterialSnapshot, 3>& snapshots,
    const DiscreteNonpenetrationInput& source, const DiscreteKnownFlowInput& flow,
    std::span<const a::PointId> particle, const DiscreteFreeWakeDirection* direction,
    const WakeTransportHistoryDirection& history, std::span<const NonpenetrationVector> seeds,
    const WakeTransportLimits& limits, const DiscreteFreeWakeLimits& velocity_limits,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return guarded([&] {
        return WakeTransportBuilder::free(clock, snapshots, source, flow, particle, direction,
                                          history, seeds, limits, velocity_limits, shared, local,
                                          kernel, memory, stop);
    });
}

core::Result<WakeTransportResult> evaluate_first_free_wake_transport(
    const WakeTransportClock& clock, const FirstFreeWakeEndpoint& current,
    const FirstFreeWakeEndpoint& parent, std::optional<a::ReleaseEpoch> birth_epoch,
    std::span<const FirstFreeWakePoint> points, std::span<const NonpenetrationVector> seeds,
    const WakeTransportLimits& limits, const DiscreteFreeWakeLimits& velocity_limits,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return guarded([&] {
        return WakeTransportBuilder::first_free(clock, current, parent, birth_epoch, points, seeds,
                                                limits, velocity_limits, shared, local, kernel,
                                                memory, stop);
    });
}

core::Result<WakeTransportResult> evaluate_fixed_prescribed_wake_transport(
    const WakeTransportClock& clock, const std::array<WakeMaterialSnapshot, 3>& snapshots,
    const PrescribedWakeClockSample& prescribed, const PrescribedWakeTransportDirection* direction,
    std::span<const NonpenetrationVector> seeds, const WakeTransportLimits& limits,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return guarded([&] {
        return WakeTransportBuilder::prescribed(clock, snapshots, prescribed, direction, seeds,
                                                limits, memory, stop);
    });
}

} // namespace uvlas::coupling::detail
