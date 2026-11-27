#include "closed_contour.hpp"

#include "certified_trigonometry_private.hpp"
#include "global_geometry_private.hpp"
#include "p1_quadrature_private.hpp"
#include "p1_shared_kernel.hpp"

#include <cmath>
#include <tuple>

namespace uvlas::aero::detail {
namespace {
using Interval = ContourInterval;
using Vector = ContourVector;
static_assert(std::numeric_limits<double>::is_iec559);
using global_detail::allocate;
using global_detail::bytes;
using global_detail::plus;
using global_detail::reject;
using global_detail::times;
constexpr auto invalid = core::FailureCode::InvalidInput;
constexpr auto unsupported = core::FailureCode::Unsupported;
constexpr auto unresolved = core::FailureCode::PrecisionUnresolved;
constexpr auto exhausted = core::FailureCode::ResourceExceeded;
constexpr Interval four_pi{0x1.921fb54442d18p+3, 0x1.921fb54442d19p+3};

template <class T, class F>
AeroResult<T> guard(GeometryContext context, F operation) noexcept {
    auto result = global_detail::guarded(context, [&]() -> GlobalResult<T> { return operation(); });
    if (!result)
        return std::unexpected(result.error().origin);
    return std::move(*result);
}
template <class T>
T take(GlobalResult<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}
template <class T>
T take(AeroResult<T> result) {
    if (!result)
        throw GlobalGeometryFailure{result.error(), {}};
    return std::move(*result);
}
struct Work {
    std::uint64_t limit, used{};
    std::stop_token stop;
    [[nodiscard]] std::uint64_t remaining() const noexcept {
        return limit - used;
    }
    void step(std::uint64_t amount = 1) {
        if (stop.stop_requested())
            reject(core::FailureCode::Cancelled);
        if (amount > remaining())
            reject(exhausted);
        used += amount;
    }
};

struct Lease {
    std::atomic_flag& flag;
    explicit Lease(std::atomic_flag& active) : flag(active) {
        if (flag.test_and_set(std::memory_order_acquire))
            reject(invalid);
    }
    ~Lease() {
        flag.clear(std::memory_order_release);
    }
};

bool valid(Interval x) {
    return std::isfinite(x.lower) && std::isfinite(x.upper) && x.lower <= x.upper;
}
bool zero(Interval x) {
    return x.lower == 0 && x.upper == 0;
}

template <bool Up>
double adjacent(double x) {
    constexpr auto sign = UINT64_C(1) << 63U;
    if (x == 0)
        return std::bit_cast<double>(Up ? UINT64_C(1) : sign | UINT64_C(1));
    const auto word = std::bit_cast<std::uint64_t>(x);
    return std::bit_cast<double>(Up == ((word & sign) == 0) ? word + 1 : word - 1);
}
Interval outward(double lower, double upper) {
    if (!valid({lower, upper}))
        reject(unresolved);
    const Interval result{adjacent<false>(lower), adjacent<true>(upper)};
    if (!valid(result))
        reject(unresolved);
    return result;
}
Interval operator-(Interval a) {
    return {-a.upper, -a.lower};
}
Interval operator+(Interval a, Interval b) {
    if (zero(a))
        return b;
    if (zero(b))
        return a;
    return outward(a.lower + b.lower, a.upper + b.upper);
}
Interval operator-(Interval a, Interval b) {
    if (a.lower == a.upper && b.lower == b.upper && a.lower == b.lower)
        return {};
    return a + -b;
}
template <class Apply>
Interval corners(Interval a, Interval b, Apply apply) {
    const auto [lower, upper] = std::minmax({apply(a.lower, b.lower), apply(a.lower, b.upper),
                                             apply(a.upper, b.lower), apply(a.upper, b.upper)});
    return outward(lower, upper);
}
Interval operator*(Interval a, Interval b) {
    if (zero(a) || zero(b))
        return {};
    return corners(a, b, [](double x, double y) { return x * y; });
}
Interval operator/(Interval a, Interval b) {
    if (b.lower <= 0 && b.upper >= 0)
        reject(unresolved);
    if (zero(a))
        return {};
    return corners(a, b, [](double x, double y) { return x / y; });
}
Interval square(Interval x) {
    if (zero(x))
        return {};
    const auto a = std::abs(x.lower), b = std::abs(x.upper);
    const auto maximum = std::max(a, b);
    const auto minimum = x.lower <= 0 && x.upper >= 0 ? 0. : std::min(a, b);
    auto result = outward(minimum * minimum, maximum * maximum);
    result.lower = std::max(0., result.lower);
    return result;
}
Interval root(Interval x) {
    if (x.upper < 0)
        reject(unresolved);
    auto result = outward(std::sqrt(std::max(x.lower, 0.)), std::sqrt(x.upper));
    result.lower = std::max(0., result.lower);
    return result;
}
Vector difference(const Vector& a, const Vector& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
Interval dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Interval norm2(const Vector& a) {
    return square(a[0]) + square(a[1]) + square(a[2]);
}
Vector cross(const Vector& a, const Vector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Interval from_words(ScalarWords words) {
    if (const auto* single = std::get_if<SingleWord>(&words))
        return {single->value, single->value};
    const auto& pair = std::get<DoubleWord>(words);
    return Interval{pair.high, pair.high} + Interval{pair.low, pair.low};
}

/// Binary carry reduction keeps roundoff accumulation logarithmic in the edge count.
class PairwiseSum {
public:
    void append(Interval value) {
        std::size_t level = 0;
        for (auto mask = count_; (mask & 1U) != 0; mask >>= 1U)
            value = levels_[level++] + value;
        levels_[level] = value;
        ++count_;
    }

    [[nodiscard]] Interval sum() const {
        Interval result{};
        for (std::size_t level = levels_.size(); level != 0; --level)
            if (((count_ >> (level - 1)) & 1U) != 0)
                result = result + levels_[level - 1];
        return result;
    }

private:
    // One append per edge in a buffer whose extent fits the signed core::Index.
    std::array<Interval, std::numeric_limits<std::size_t>::digits> levels_{};
    std::size_t count_{};
};

struct Term {
    std::size_t edge{}, column{};
    int sign{};
};
struct Edge {
    EdgeId id{0};
    std::array<std::size_t, 2> vertices{};
    Vector delta{};
    Interval length2{};
    std::size_t first{}, last{};
    bool used{};
};
} // namespace

struct ClosedContourEndpoint {
    Vector relative{};
    Interval distance{};
};
struct ClosedContourData {
    std::optional<ContourContactIndex> contacts;
    core::PlannedBuffer<Vector> vertices;
    core::PlannedBuffer<Edge> edges;
    core::PlannedBuffer<Term> terms;
    Interval core2{};
    std::size_t columns{};
    std::uint64_t steps{};
};

namespace {
template <class Scalar>
Scalar same_sign_integral(Scalar s, Scalar t, Scalar first_radius, Scalar second_radius) {
    return (s + t) / (first_radius * second_radius * (s * second_radius + t * first_radius));
}

Vector edge_value(const Edge& edge, std::span<const ClosedContourEndpoint> endpoints,
                  Interval core2) {
    const auto& a = endpoints[edge.vertices[0]];
    const auto& b = endpoints[edge.vertices[1]];
    const auto s = dot(edge.delta, a.relative);
    const auto z = s - edge.length2;
    const auto c = cross(edge.delta, a.relative);
    const bool same = (s.lower > 0 && z.lower > 0) || (s.upper < 0 && z.upper < 0);
    const auto integral =
        same ? same_sign_integral(s, z, a.distance, b.distance)
             : (s / a.distance - z / b.distance) / (norm2(c) + core2 * edge.length2);
    const auto factor = -integral / four_pi;
    return {c[0] * factor, c[1] * factor, c[2] * factor};
}

struct AgeDerivative {
    Interval value, derivative;
};
AgeDerivative operator+(AgeDerivative a, AgeDerivative b) {
    return {a.value + b.value, a.derivative + b.derivative};
}
AgeDerivative operator-(AgeDerivative a, AgeDerivative b) {
    return {a.value - b.value, a.derivative - b.derivative};
}
AgeDerivative operator*(AgeDerivative a, AgeDerivative b) {
    return {a.value * b.value, a.derivative * b.value + a.value * b.derivative};
}
AgeDerivative operator/(AgeDerivative a, AgeDerivative b) {
    const auto quotient = a.value / b.value;
    return {quotient, (a.derivative - quotient * b.derivative) / b.value};
}

Interval project_age_derivative(const Vector& delta, Interval length2, const Vector& velocity,
                                Interval core2, const Vector& first, const Vector& second,
                                const Vector& normal) {
    const auto first_radius = root(norm2(first) + core2);
    const auto second_radius = root(norm2(second) + core2);
    const auto cross_span = cross(delta, first);
    const auto cross_derivative = cross(delta, velocity);
    // This is the query derivative of the forward span edge. It is also the
    // translation derivative of the oppositely directed downstream edge.
    const AgeDerivative s{dot(delta, first), dot(delta, velocity)};
    const AgeDerivative t{s.value - length2, s.derivative};
    const AgeDerivative a{first_radius, dot(first, velocity) / first_radius};
    const AgeDerivative b{second_radius, dot(second, velocity) / second_radius};
    const bool same =
        (s.value.lower > 0 && t.value.lower > 0) || (s.value.upper < 0 && t.value.upper < 0);
    const auto integral =
        same ? same_sign_integral(s, t, a, b)
             : (s / a - t / b) / AgeDerivative{norm2(cross_span) + core2 * length2,
                                               Interval{2, 2} * dot(cross_span, cross_derivative)};
    const AgeDerivative projected{dot(cross_span, normal), dot(cross_derivative, normal)};
    const auto right =
        dot(cross(velocity, second), normal) / (square(second_radius) * second_radius);
    const auto left = dot(cross(velocity, first), normal) / (square(first_radius) * first_radius);
    return -((projected * integral).derivative + right - left) / four_pi;
}

template <class Visit>
std::uint64_t evaluate(const ClosedContourData& data, const ResolvedQuery& query,
                       std::optional<EdgeId> self, const GeometryLimits& query_limits,
                       GlobalLimits limits, GeometryWorkspace& geometry_work,
                       GlobalWorkspace& contact_work, std::span<ClosedContourEndpoint> endpoints,
                       const std::stop_token& stop, Visit visit) {
    Work work{limits.max_steps, 0, stop};
    work.step(0);
    if (!data.contacts || endpoints.size() != data.vertices.size())
        reject(invalid);
    if (!p1_environment_supported())
        reject(core::FailureCode::PrecisionUnsupported);
    work.step(take(data.contacts->qualify(query, self, query_limits, limits, geometry_work,
                                          contact_work, stop))
                  .steps);
    Vector location{};
    const auto expression = take(query.definition());
    if (self) {
        const auto edges = data.edges.values();
        const auto found = std::ranges::lower_bound(edges, *self, {}, &Edge::id);
        if (found == edges.end() || found->id != *self || !found->used)
            reject(unsupported, *self);
        const auto points = data.vertices.values();
        for (std::size_t axis = 0; axis < 3; ++axis)
            location[axis] = points[found->vertices[0]][axis] * Interval{.5, .5} +
                             points[found->vertices[1]][axis] * Interval{.5, .5};
    } else {
        const auto* point = std::get_if<FreePoint>(&expression.definition);
        if (point == nullptr)
            reject(unsupported, expression.id);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (const auto* single = std::get_if<FreeSingle>(&point->words))
                location[axis] = {single->words[axis], single->words[axis]};
            else if (const auto* pair = std::get_if<FreeDouble>(&point->words))
                location[axis] = from_words(DoubleWord{pair->high[axis], pair->low[axis]});
        }
    }
    const auto vertices = data.vertices.values();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        work.step();
        auto& endpoint = endpoints[i];
        endpoint.relative = difference(location, vertices[i]);
        endpoint.distance = root(norm2(endpoint.relative) + data.core2);
    }
    for (const auto& edge : data.edges.values()) {
        work.step();
        if (edge.first == edge.last || self == edge.id)
            continue;
        const auto velocity = edge_value(edge, endpoints, data.core2);
        const auto terms = data.terms.values().subspan(edge.first, edge.last - edge.first);
        work.step(terms.size());
        visit(velocity, terms);
    }
    work.step(0);
    return work.used;
}
} // namespace

AeroResult<StraightWakeAgeQuadrature>
prepare_straight_wake_age_quadrature(const GeometryContext& context,
                                     std::span<const double> boundaries, double frequency,
                                     std::uint64_t maximum_samples, std::uint64_t maximum_steps,
                                     core::ResourcePlan& plan, std::stop_token stop) noexcept {
    return guard<StraightWakeAgeQuadrature>(context, [&] {
        if (stop.stop_requested())
            reject(core::FailureCode::Cancelled);
        if (!p1_environment_supported())
            reject(core::FailureCode::PrecisionUnsupported);
        if (boundaries.size() < 2 || boundaries.front() != 0 || !std::isfinite(frequency))
            reject(invalid);
        constexpr auto order = p1::gauss_16.size();
        const auto count = times(boundaries.size() - 1, order);
        if (count > maximum_samples)
            reject(exhausted);
        p1::Work work{0, maximum_steps, stop};
        try {
            for (std::size_t i = 1; i < boundaries.size(); ++i) {
                work.step();
                if (!std::isfinite(boundaries[i]) || boundaries[i] <= boundaries[i - 1])
                    reject(invalid);
            }
            std::array<p1::Interval, order> nodes, weights;
            p1::build_rule(p1::gauss_16, std::span{nodes}, std::span{weights}, work);
            const auto half_frequency = p1::scale(decode_p1_word(frequency), -1);
            const auto pi = p1::pi_constant(work);
            const auto rounded = [](const p1::Bound& value, bool upper) {
                if (value.mantissa == 0)
                    return 0.;
                const auto top = static_cast<std::int64_t>(value.power) + p1::bits(value.mantissa);
                const auto quantum = p1::exponent(std::max<std::int64_t>(-1074, top - 53));
                const auto head = p1::at_power(value, quantum, upper);
                const auto result = std::scalbn(head.convert_to<double>(), quantum);
                if (!std::isfinite(result))
                    reject(unresolved);
                return result;
            };
            const auto enclosure = [&](const p1::Interval& value) {
                return Interval{rounded(value.lower, false), rounded(value.upper, true)};
            };
            core::PlannedBuffer<StraightWakeAgeSample> samples;
            allocate(samples, count, plan, context);
            std::size_t out = 0;
            auto left = decode_p1_word(boundaries.front());
            for (const auto hi : boundaries.subspan(1)) {
                const auto right = decode_p1_word(hi);
                const auto width = p1::subtract(right, left, work);
                for (std::size_t i = 0; i < order; ++i) {
                    const auto age = p1::add(left, p1::multiply(width, nodes[i], work), work);
                    const auto q = p1::multiply(width, weights[i], work);
                    const auto half_phase = p1::multiply(age, half_frequency, work);
                    const auto phase =
                        p1::sine_cosine(half_phase, p1::PhaseUnits::Radians, pi, work);
                    const auto weighted =
                        p1::negate(p1::scale(p1::multiply(q, phase.sine, work), 1), work);
                    samples.values()[out++] = {
                        enclosure(age), enclosure(p1::multiply(weighted, phase.sine, work)),
                        enclosure(p1::multiply(weighted, phase.cosine, work))};
                }
                left = right;
            }
            work.step();
            return StraightWakeAgeQuadrature{std::move(samples), work.used};
        } catch (const p1::ArithmeticFailure& failure) {
            reject(failure.code);
        }
    });
}

AeroResult<StraightWakeColumnUpdate> project_straight_wake_column_update(
    const GeometryContext& context, const ContourVector& first, const ContourVector& second,
    const ContourVector& velocity, ContourInterval core_radius, const ContourVector& point,
    const ContourVector& normal, std::span<const StraightWakeAgeSample> samples,
    std::uint64_t maximum_samples, std::stop_token stop) noexcept {
    return guard<StraightWakeColumnUpdate>(context, [&] {
        Work work{maximum_samples, 0, stop};
        work.step(0);
        if (!p1_environment_supported())
            reject(core::FailureCode::PrecisionUnsupported);
        if (!valid(core_radius) || core_radius.lower <= 0 || samples.empty())
            reject(invalid);
        for (const auto* vector : {&first, &second, &velocity, &point, &normal})
            if (!std::ranges::all_of(*vector, valid))
                reject(invalid);
        if (samples.size() > work.remaining())
            reject(exhausted);
        const auto delta = difference(second, first);
        const auto area2 = norm2(cross(delta, velocity));
        if (area2.upper == 0)
            reject(invalid);
        if (area2.lower <= 0)
            reject(unresolved);
        const auto length2 = norm2(delta);
        const auto core2 = square(core_radius);
        const auto first_offset = difference(point, first);
        const auto second_offset = difference(point, second);
        PairwiseSum real, imaginary;
        for (const auto& sample : samples) {
            work.step();
            if (!valid(sample.age_seconds) || sample.age_seconds.lower < 0 ||
                !valid(sample.real_weight_seconds) || !valid(sample.imaginary_weight_seconds))
                reject(invalid);
            if (zero(sample.real_weight_seconds) && zero(sample.imaginary_weight_seconds))
                continue;
            const Vector displacement{sample.age_seconds * velocity[0],
                                      sample.age_seconds * velocity[1],
                                      sample.age_seconds * velocity[2]};
            const auto derivative = project_age_derivative(
                delta, length2, velocity, core2, difference(first_offset, displacement),
                difference(second_offset, displacement), normal);
            real.append(derivative * sample.real_weight_seconds);
            imaginary.append(derivative * sample.imaginary_weight_seconds);
        }
        work.step(0);
        return StraightWakeColumnUpdate{real.sum(), imaginary.sum()};
    });
}

ClosedContourWorkspace::ClosedContourWorkspace(
    core::PlannedBuffer<ClosedContourEndpoint>&& endpoints) noexcept
    : endpoints_(std::move(endpoints)) {}
ClosedContourWorkspace::ClosedContourWorkspace(ClosedContourWorkspace&& other) noexcept
    : endpoints_(std::move(other.endpoints_)) {}
ClosedContourWorkspace& ClosedContourWorkspace::operator=(ClosedContourWorkspace&& other) noexcept {
    endpoints_ = std::move(other.endpoints_);
    return *this;
}
ClosedContourWorkspace::~ClosedContourWorkspace() = default;

ClosedContour::ClosedContour(core::PlannedBuffer<ClosedContourData>&& data,
                             GeometryContext context) noexcept
    : data_(std::move(data)), context_(context) {}
ClosedContour::ClosedContour(ClosedContour&&) noexcept = default;
ClosedContour& ClosedContour::operator=(ClosedContour&&) noexcept = default;
ClosedContour::~ClosedContour() = default;

AeroResult<ClosedContour> ClosedContour::prepare(const GlobalContactState& contact,
                                                 std::span<const P1FaceStrengths> maps,
                                                 std::span<const P1StrengthId> ids,
                                                 GlobalLimits limits, GlobalWorkspace& workspace,
                                                 core::ResourcePlan& plan,
                                                 std::stop_token stop) noexcept {
    const auto context = contact.context_;
    return guard<ClosedContour>(context, [&] {
        Work work{limits.max_steps, 0, stop};
        work.step(0);
        const auto& owned = take(contact.geometry()).get();
        const auto& geometry = *GlobalEngine::source(owned);
        const auto& topology = *contact.data_.values().front().topology;
        if (geometry.triangles.size() != 0 || maps.size() != geometry.faces.size() ||
            geometry.quads.size() != maps.size())
            reject(unsupported);
        for (std::size_t i = 1; i < ids.size(); ++i) {
            work.step();
            if (ids[i - 1] >= ids[i])
                reject(invalid);
        }
        for (std::size_t i = 1; i < maps.size(); ++i) {
            work.step();
            if (maps[i - 1].face >= maps[i].face)
                reject(invalid);
        }
        const auto term_count = times(geometry.facets.size(), 3);
        const auto payload =
            plus(bytes<ClosedContourData>(1),
                 plus(bytes<Vector>(geometry.points.size()),
                      plus(bytes<Edge>(geometry.edges.size()), bytes<Term>(term_count))));
        const auto index_payload =
            plus(bytes<GlobalContactData>(1), bytes<global_detail::Node>(topology.nodes.size()));
        global_detail::byte_admit(plus(payload, index_payload), plan, context);
        limits.max_steps = work.remaining();
        auto index = take(ContourContactIndex::prepare(contact, limits, workspace, plan, stop));
        work.step(index.preparation_steps());
        global_detail::byte_admit(payload, plan, context);
        core::PlannedBuffer<ClosedContourData> storage;
        allocate(storage, 1, plan, context);
        auto& data = storage.values().front();
        data.contacts.emplace(std::move(index));
        data.columns = ids.size();
        allocate(data.vertices, geometry.points.size(), plan, context);
        allocate(data.edges, geometry.edges.size(), plan, context);
        allocate(data.terms, term_count, plan, context);
        if (!p1_environment_supported())
            reject(core::FailureCode::PrecisionUnsupported);
        data.core2 = square(from_words(geometry.radius));
        const auto high = geometry.current_high.values();
        const auto low = geometry.current_low.values();
        auto points = data.vertices.values();
        for (std::size_t i = 0; i < points.size(); ++i) {
            work.step();
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const std::size_t coordinate = 3 * i + axis;
                points[i][axis] = low.empty()
                                      ? Interval{high[coordinate], high[coordinate]}
                                      : from_words(DoubleWord{high[coordinate], low[coordinate]});
            }
        }
        auto edges = data.edges.values();
        for (std::size_t i = 0; i < edges.size(); ++i) {
            work.step();
            const auto& source = geometry.edges.values()[i];
            edges[i].id = EdgeId{source.id};
            edges[i].vertices = {static_cast<std::size_t>(source.points[0]),
                                 static_cast<std::size_t>(source.points[1])};
        }
        auto terms = data.terms.values();
        std::size_t entered = 0;
        for (const auto& facet : geometry.facets.values()) {
            work.step();
            const auto index_of_face = static_cast<std::size_t>(facet.face);
            const auto face = PhysicalFaceId{geometry.faces.values()[index_of_face].id};
            const auto map = std::ranges::lower_bound(
                maps, face,
                [&](PhysicalFaceId a, PhysicalFaceId b) {
                    work.step();
                    return a < b;
                },
                &P1FaceStrengths::face);
            if (map == maps.end() || map->face != face)
                reject(invalid, face);
            const bool bound =
                topology.roles.values()[index_of_face] == GlobalTopologyData::bound_role;
            auto key = P1StrengthId::bound(face);
            if (bound) {
                if (!std::holds_alternative<P0BoundColumn>(map->columns))
                    reject(invalid, face);
            } else {
                const auto* constant = std::get_if<ConstantWakeColumn>(&map->columns);
                if (constant == nullptr)
                    reject(unsupported, face);
                key = constant->id;
            }
            const auto handle =
                std::ranges::lower_bound(ids, key, [&](P1StrengthId a, P1StrengthId b) {
                    work.step();
                    return a < b;
                });
            if (handle == ids.end() || *handle != key)
                reject(invalid, face);
            const auto column = static_cast<std::size_t>(handle - ids.begin());
            for (std::size_t j = 0; j < facet.edges.size(); ++j) {
                work.step();
                const auto raw = static_cast<std::size_t>(facet.edges[j]);
                const auto owner = static_cast<std::size_t>(topology.edges.values()[raw].owner);
                terms[entered++] = {owner, column, static_cast<int>(facet.incidence[j])};
            }
        }
        std::ranges::sort(terms, [&](const Term& a, const Term& b) {
            work.step();
            return std::tie(a.edge, a.column) < std::tie(b.edge, b.column);
        });
        std::size_t cursor = 0, written = 0;
        while (cursor < terms.size()) {
            const auto edge_index = terms[cursor].edge;
            auto& edge = edges[edge_index];
            edge.used = true;
            edge.first = written;
            do {
                const auto column = terms[cursor].column;
                std::int64_t sign = 0;
                do {
                    work.step();
                    sign += terms[cursor++].sign;
                } while (cursor < terms.size() && terms[cursor].edge == edge_index &&
                         terms[cursor].column == column);
                if (sign < -1 || sign > 1)
                    reject(unsupported, edge.id);
                if (sign != 0)
                    terms[written++] = {edge_index, column, static_cast<int>(sign)};
            } while (cursor < terms.size() && terms[cursor].edge == edge_index);
            edge.last = written;
            if (edge.first != edge.last) {
                edge.delta = difference(points[edge.vertices[1]], points[edge.vertices[0]]);
                edge.length2 = norm2(edge.delta);
            }
        }
        std::ranges::sort(edges, [&](const Edge& a, const Edge& b) {
            work.step();
            return a.id < b.id;
        });
        work.step(0);
        data.steps = work.used;
        return ClosedContour{std::move(storage), context};
    });
}

AeroResult<ClosedContourWorkspace>
ClosedContour::workspace(core::ResourcePlan& plan) const noexcept {
    if (data_.size() == 0)
        return std::unexpected(global_detail::failure(invalid, context_).origin);
    const auto& data = data_.values().front();
    return guard<ClosedContourWorkspace>(context_, [&] {
        core::PlannedBuffer<ClosedContourEndpoint> storage;
        allocate(storage, data.vertices.size(), plan, context_);
        return ClosedContourWorkspace{std::move(storage)};
    });
}

std::uint64_t ClosedContour::preparation_steps() const noexcept {
    return data_.size() == 0 ? 0 : data_.values().front().steps;
}

AeroResult<std::uint64_t>
ClosedContour::project(const ResolvedQuery& query, const ContourVector& normal,
                       std::span<ContourInterval> row, const GeometryLimits& query_limits,
                       GlobalLimits limits, GeometryWorkspace& geometry_work,
                       GlobalWorkspace& contact_work, ClosedContourWorkspace& work,
                       std::stop_token stop) const noexcept {
    return guard<std::uint64_t>(query.context(), [&] {
        const Lease lease{work.entered_};
        if (data_.size() == 0 || data_.values().front().columns != row.size() ||
            !std::ranges::all_of(normal, valid))
            reject(invalid);
        if (row.size() > limits.max_steps)
            reject(exhausted);
        if (stop.stop_requested())
            reject(core::FailureCode::Cancelled);
        std::ranges::fill(row, Interval{});
        limits.max_steps -= row.size();
        const auto used = evaluate(data_.values().front(), query, {}, query_limits, limits,
                                   geometry_work, contact_work, work.endpoints_.values(), stop,
                                   [&](const Vector& velocity, std::span<const Term> terms) {
                                       const auto value = dot(normal, velocity);
                                       for (const auto& term : terms)
                                           row[term.column] =
                                               row[term.column] + (term.sign == 1 ? value : -value);
                                   });
        return used + row.size();
    });
}

AeroResult<std::pair<ContourVector, std::uint64_t>>
ClosedContour::apply(const ResolvedQuery& query, std::optional<EdgeId> self,
                     std::span<const ContourInterval> strengths, const GeometryLimits& query_limits,
                     GlobalLimits limits, GeometryWorkspace& geometry_work,
                     GlobalWorkspace& contact_work, ClosedContourWorkspace& work,
                     std::stop_token stop) const noexcept {
    return guard<std::pair<ContourVector, std::uint64_t>>(query.context(), [&] {
        const Lease lease{work.entered_};
        if (data_.size() == 0 || data_.values().front().columns != strengths.size())
            reject(invalid);
        Work input_work{limits.max_steps, 0, stop};
        for (const auto value : strengths) {
            input_work.step();
            if (!valid(value))
                reject(invalid);
        }
        limits.max_steps = input_work.remaining();
        std::array<PairwiseSum, 3> sum{};
        const auto used = evaluate(data_.values().front(), query, self, query_limits, limits,
                                   geometry_work, contact_work, work.endpoints_.values(), stop,
                                   [&](const Vector& velocity, std::span<const Term> terms) {
                                       Interval weight{};
                                       for (const auto& term : terms) {
                                           const auto value = strengths[term.column];
                                           weight = weight + (term.sign == 1 ? value : -value);
                                       }
                                       for (std::size_t axis = 0; axis < sum.size(); ++axis)
                                           sum[axis].append(weight * velocity[axis]);
                                   });
        return std::pair{Vector{sum[0].sum(), sum[1].sum(), sum[2].sum()}, used + input_work.used};
    });
}

} // namespace uvlas::aero::detail
