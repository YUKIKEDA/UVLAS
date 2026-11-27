#include "p1_shared_action.hpp"

#include "geometry_snapshot_storage.hpp"
#include "global_geometry_storage.hpp"
#include "local_geometry_bounds_private.hpp"
#include "p1_interval_private.hpp"
#include "p1_shared_kernel.hpp"
#include "shared_p1_reuse_private.hpp"

#include <bit>
#include <tuple>

namespace uvlas::aero::detail {
namespace shared_p1 {
using p1::ArithmeticFailure;
using p1::Interval;
using p1::Work;

constexpr auto operation = core::Operation::EvaluateAero;
constexpr std::uint64_t allocator_estimate = 4 * sizeof(void*);
constexpr auto invalid = core::FailureCode::InvalidInput;
constexpr auto stale = core::FailureCode::StaleGeneration;
constexpr auto exceeded = core::FailureCode::ResourceExceeded;

[[nodiscard]] OwnedGeometryFailure failure(core::FailureCode code,
                                           GeometryContext context) noexcept {
    return {core::Failure{code, operation, context.diagnostics}, {}, context};
}
template <class T>
T needed(AeroResult<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}

[[nodiscard]] std::uint64_t sum(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        throw ArithmeticFailure{exceeded};
    return a + b;
}
[[nodiscard]] std::uint64_t product(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        throw ArithmeticFailure{exceeded};
    return a * b;
}
template <class T>
[[nodiscard]] std::uint64_t bytes(std::uint64_t count) {
    if (count > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        throw ArithmeticFailure{exceeded};
    return product(count, sizeof(T));
}
template <class T>
[[nodiscard]] core::PlannedBuffer<T> allocate(core::ResourcePlan& memory, std::uint64_t count,
                                              GeometryContext context,
                                              const std::stop_token& stop) {
    if (stop.stop_requested())
        throw ArithmeticFailure{core::FailureCode::Cancelled};
    (void)bytes<T>(count);
    auto created = memory.allocate<T>(static_cast<core::Index>(count), operation);
    if (!created)
        throw OwnedGeometryFailure{created.error(), {}, context};
    return std::move(*created);
}
template <class T, class Action>
[[nodiscard]] AeroResult<T> guarded(GeometryContext context, const Action& action) noexcept {
    try {
        return action();
    } catch (const OwnedGeometryFailure& reason) {
        return std::unexpected(reason);
    } catch (const ArithmeticFailure& reason) {
        return std::unexpected(failure(reason.code, context));
    } catch (const std::bad_alloc&) {
        return std::unexpected(failure(exceeded, context));
    } catch (...) {
        return std::unexpected(failure(core::FailureCode::PrecisionUnresolved, context));
    }
}

struct Lease {
    std::atomic_flag& flag;
    ~Lease() {
        flag.clear(std::memory_order_release);
    }
};

struct Complex {
    Interval real, imaginary;
};
struct MaterialPoint {
    std::array<std::size_t, 3> point{};
    std::array<Interval, 3> weight{};
    std::size_t count = 0;
};
struct EdgeRange {
    std::size_t first = 0, last = 0;
    bool positive = false;
};
struct Facet {
    std::array<std::size_t, 3> point{}, strength{};
    std::array<bool, 3> boundary{true, true, true};
    std::optional<P1Basis> basis;
    std::array<EdgeRange, 3> bound_edges{};
    bool far = false, bound = false;
};
struct EdgeOccurrence {
    std::size_t owner = 0;
    std::array<std::size_t, 2> trace{};
    std::size_t facet = 0, edge = 0;
    bool positive = false;
    bool retained_bound = false;
};
struct PointGradient {
    std::array<Interval, 3> value{};
};
struct CoordinateResult {
    std::array<std::optional<ScaledEnclosure>, 3> value;
};

[[nodiscard]] Interval words(ScalarWords original, Work& work) {
    work.step();
    if (const auto* single = std::get_if<SingleWord>(&original))
        return decode_p1_word(single->value);
    const auto* pair = std::get_if<DoubleWord>(&original);
    if (pair == nullptr)
        throw ArithmeticFailure{invalid};
    return p1::add(decode_p1_word(pair->high), decode_p1_word(pair->low), work);
}
[[nodiscard]] Complex words(const P1ComplexWords& original, Work& work) {
    return {words(original.real, work), words(original.imaginary, work)};
}
[[nodiscard]] Interval enclosure(const ScaledEnclosure& original, Work& work) {
    const auto center = words(original.center(), work);
    const auto radius = decode_p1_word(original.radius());
    const auto lower = p1::subtract(center, radius, work);
    const auto upper = p1::add(center, radius, work);
    return p1::scale(Interval{lower.lower, upper.upper}, original.exponent());
}
[[nodiscard]] Complex add(const Complex& a, const Complex& b, Work& work) {
    return {p1::add(a.real, b.real, work), p1::add(a.imaginary, b.imaginary, work)};
}
[[nodiscard]] Complex multiply(const Complex& a, const Interval& b, Work& work) {
    return {p1::multiply(a.real, b, work), p1::multiply(a.imaginary, b, work)};
}
[[nodiscard]] Interval dot(const Complex& left, const Complex& right, Work& work) {
    return p1::add(p1::multiply(left.real, right.real, work),
                   p1::multiply(left.imaginary, right.imaginary, work), work);
}
[[nodiscard]] P1ComplexEnclosure finish(const Complex& value, core::Precision precision,
                                        LocalRootScratch& scratch) {
    return {enclose_p1_value(value.real, precision, scratch),
            enclose_p1_value(value.imaginary, precision, scratch)};
}
[[nodiscard]] std::array<P1ComplexEnclosure, 3>
finish(const std::array<Complex, 3>& value, core::Precision precision, LocalRootScratch& scratch) {
    return {finish(value[0], precision, scratch), finish(value[1], precision, scratch),
            finish(value[2], precision, scratch)};
}

[[nodiscard]] bool same_word(double a, double b) noexcept {
    return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}
[[nodiscard]] bool same_words(const ScalarWords& a, const ScalarWords& b) noexcept {
    if (const auto* left = std::get_if<SingleWord>(&a)) {
        const auto* right = std::get_if<SingleWord>(&b);
        return right != nullptr && same_word(left->value, right->value);
    }
    const auto* left = std::get_if<DoubleWord>(&a);
    const auto* right = std::get_if<DoubleWord>(&b);
    return left != nullptr && right != nullptr && same_word(left->high, right->high) &&
           same_word(left->low, right->low);
}
[[nodiscard]] bool same_words(const PointWords& a, const PointWords& b) noexcept {
    if (const auto* left = std::get_if<FreeSingle>(&a)) {
        const auto* right = std::get_if<FreeSingle>(&b);
        return right != nullptr && std::ranges::equal(left->words, right->words, same_word);
    }
    const auto* left = std::get_if<FreeDouble>(&a);
    const auto* right = std::get_if<FreeDouble>(&b);
    return left != nullptr && right != nullptr &&
           std::ranges::equal(left->high, right->high, same_word) &&
           std::ranges::equal(left->low, right->low, same_word);
}
[[nodiscard]] bool same_query(const QueryRequest& a, const QueryRequest& b) noexcept {
    if (a.id != b.id)
        return false;
    if (const auto* left = std::get_if<FreePoint>(&a.definition)) {
        const auto* right = std::get_if<FreePoint>(&b.definition);
        return right != nullptr && same_words(left->words, right->words);
    }
    if (const auto* left = std::get_if<SourceVertex>(&a.definition)) {
        const auto* right = std::get_if<SourceVertex>(&b.definition);
        return right != nullptr && left->point == right->point;
    }
    if (const auto* left = std::get_if<SourceEdgePoint>(&a.definition)) {
        const auto* right = std::get_if<SourceEdgePoint>(&b.definition);
        return right != nullptr && left->edge == right->edge && same_words(left->t, right->t);
    }
    const auto* left = std::get_if<SourceFacePoint>(&a.definition);
    const auto* right = std::get_if<SourceFacePoint>(&b.definition);
    return left != nullptr && right != nullptr && left->facet == right->facet &&
           same_words(left->u, right->u) && same_words(left->v, right->v);
}

template <class T, class Id, class Key>
[[nodiscard]] std::size_t locate(std::span<const T> sorted, Id id, Key key, Work& work) {
    const auto position =
        std::lower_bound(sorted.begin(), sorted.end(), id, [&](const T& item, Id value) {
            work.step();
            return std::invoke(key, item) < value;
        });
    if (position == sorted.end() || std::invoke(key, *position) != id)
        throw ArithmeticFailure{invalid};
    return static_cast<std::size_t>(position - sorted.begin());
}

[[nodiscard]] std::size_t ordinal(std::span<const IdLookup> sorted, std::uint64_t id, Work& work) {
    return static_cast<std::size_t>(sorted[locate(sorted, id, &IdLookup::id, work)].ordinal);
}

[[nodiscard]] MaterialPoint material(const GeometryData& geometry, const QueryRequest& query,
                                     std::span<const std::size_t> columns, Work& work) {
    MaterialPoint result;
    if (const auto* vertex = std::get_if<SourceVertex>(&query.definition)) {
        const auto index = ordinal(geometry.point_lookup.values(), vertex->point.value(), work);
        result.point[0] = columns[static_cast<std::size_t>(geometry.points.values()[index].owner)];
        result.weight[0] = p1::integer(1);
        result.count = 1;
    }
    if (const auto* edge = std::get_if<SourceEdgePoint>(&query.definition)) {
        const auto index = ordinal(geometry.edge_lookup.values(), edge->edge.value(), work);
        const auto& endpoints = geometry.edges.values()[index].points;
        result.point[0] = columns[static_cast<std::size_t>(endpoints[0])];
        result.point[1] = columns[static_cast<std::size_t>(endpoints[1])];
        result.weight[1] = words(edge->t, work);
        result.weight[0] = p1::subtract(p1::integer(1), result.weight[1], work);
        result.count = 2;
    }
    if (const auto* facet = std::get_if<SourceFacePoint>(&query.definition)) {
        const auto index = ordinal(geometry.facet_lookup.values(), facet->facet.value(), work);
        const auto& points = geometry.facets.values()[index].points;
        for (std::size_t i = 0; i < points.size(); ++i)
            result.point[i] = columns[static_cast<std::size_t>(points[i])];
        result.weight[1] = words(facet->u, work);
        result.weight[2] = words(facet->v, work);
        result.weight[0] = p1::subtract(p1::subtract(p1::integer(1), result.weight[1], work),
                                        result.weight[2], work);
        result.count = 3;
    }
    return result;
}

[[nodiscard]] auto edge_key(const EdgeOccurrence& edge) noexcept {
    return std::tuple{edge.owner, edge.trace, edge.positive, edge.facet, edge.edge};
}
[[nodiscard]] std::span<EdgeOccurrence> sort_edges(std::span<EdgeOccurrence> source,
                                                   std::span<EdgeOccurrence> scratch, Work& work) {
    const auto size = source.size();
    for (std::size_t width = 1; width < size; width += std::min(width, size - width)) {
        for (std::size_t first = 0; first < size;) {
            const auto middle = first + std::min(width, size - first);
            const auto last = middle + std::min(width, size - middle);
            auto left = first;
            auto right = middle;
            for (std::size_t out = first; out < last; ++out) {
                work.step();
                if (right == last ||
                    (left < middle && edge_key(source[left]) < edge_key(source[right])))
                    scratch[out] = source[left++];
                else
                    scratch[out] = source[right++];
            }
            first = last;
        }
        std::swap(source, scratch);
    }
    return source;
}

[[nodiscard]] std::uint64_t cancel_pairs(std::span<const EdgeOccurrence> edges,
                                         std::span<Facet> facets, Work& work) {
    std::uint64_t cancellations = 0;
    for (std::size_t first = 0; first < edges.size();) {
        std::size_t last = first + 1;
        while (last < edges.size() && edges[last].owner == edges[first].owner &&
               edges[last].trace == edges[first].trace) {
            work.step();
            ++last;
        }
        auto positive = first;
        while (positive < last && !edges[positive].positive) {
            work.step();
            ++positive;
        }
        auto negative = first;
        const auto split = positive;
        while (negative < split && positive < last) {
            work.step();
            const auto& a = edges[negative];
            const auto& b = edges[positive];
            // Surface quadrature carries its whole facet. Mask only a pair of actual
            // analytic boundary terms from the same oriented trace.
            if (facets[a.facet].far) {
                ++negative;
            } else if (facets[b.facet].far) {
                ++positive;
            } else {
                facets[a.facet].boundary[a.edge] = false;
                facets[b.facet].boundary[b.edge] = false;
                ++negative;
                ++positive;
                ++cancellations;
            }
        }
        first = last;
    }
    return cancellations;
}
/// Pick one representative analytic edge per canonical bound-edge owner.
void group_bound_edges(std::span<EdgeOccurrence> edges, std::span<Facet> facets, Work& work) {
    for (std::size_t first = 0; first < edges.size();) {
        auto last = first;
        auto representative = edges.size();
        while (last < edges.size() && edges[last].owner == edges[first].owner) {
            work.step();
            auto& edge = edges[last];
            auto& facet = facets[edge.facet];
            edge.retained_bound = facet.bound && facet.boundary[edge.edge];
            if (edge.retained_bound) {
                if (representative == edges.size())
                    representative = last;
                else
                    facet.boundary[edge.edge] = false;
            }
            ++last;
        }
        if (representative != edges.size()) {
            const auto& edge = edges[representative];
            facets[edge.facet].bound_edges[edge.edge] = {first, last, edge.positive};
        }
        first = last;
    }
}
} // namespace shared_p1

struct BoundWakeFacetStorageData {
    core::PlannedBuffer<shared_p1::Facet> idle;
    core::PlannedBuffer<std::optional<LocalGeometryPacket>> packets;
    core::PlannedBuffer<TriangleOutcome> predicate;
    core::PlannedBuffer<std::uint64_t> points;
    core::PlannedBuffer<P1StrengthId> strengths;
    core::PlannedBuffer<std::size_t> point_columns;
    core::PlannedBuffer<shared_p1::EdgeOccurrence> edges;
    core::PlannedBuffer<shared_p1::EdgeOccurrence> edge_scratch;
    std::atomic_bool checked_out = false;
};

struct SharedP1Data {
    std::shared_ptr<const GlobalTopologyData> topology;
    std::optional<OwnedGeometry> geometry;
    QueryRequest query{QueryId{0}, FreePoint{FreeSingle{}}};
    P1ActionStamp stamp;
    shared_p1::MaterialPoint material;
    core::PlannedBuffer<std::uint64_t> points;
    core::PlannedBuffer<P1StrengthId> strengths;
    core::PlannedBuffer<shared_p1::Facet> facets;
    core::PlannedBuffer<shared_p1::EdgeOccurrence> occurrences;
    P1AssemblyReport report;
    std::size_t canonical_points = 0;
    core::Precision precision = core::Precision::Binary64;
    bool values_only{};
};

namespace shared_p1 {
[[nodiscard]] Complex oriented(const Complex& value, bool forward, Work& work) {
    if (forward)
        return value;
    return {p1::negate(value.real, work), p1::negate(value.imaginary, work)};
}

template <class Operation>
void visit_coefficients(const SharedP1Data& action, const Facet& facet, std::size_t column,
                        Work& work, Operation&& visit) {
    if (!facet.bound) {
        std::invoke(visit, facet.strength[column], true);
        return;
    }
    const auto& range = facet.bound_edges[column];
    for (std::size_t i = range.first; i < range.last; ++i) {
        work.step();
        const auto& edge = action.occurrences.values()[i];
        if (edge.retained_bound)
            std::invoke(visit, action.facets.values()[edge.facet].strength[edge.edge],
                        edge.positive == range.positive);
    }
}

[[nodiscard]] Complex coefficient(const SharedP1Data& action, const Facet& facet,
                                  std::size_t column, std::span<const P1ComplexWords> strengths,
                                  Work& work) {
    if (!facet.bound)
        return words(strengths[facet.strength[column]], work);
    Complex result{};
    visit_coefficients(action, facet, column, work, [&](std::size_t id, bool forward) {
        result = add(result, oriented(words(strengths[id], work), forward, work), work);
    });
    return result;
}
} // namespace shared_p1

struct SharedP1AdjointData {
    std::shared_ptr<const SharedP1Data> parent;
    core::PlannedBuffer<shared_p1::CoordinateResult> points;
    core::PlannedBuffer<std::optional<P1ComplexEnclosure>> strengths;
    shared_p1::CoordinateResult query;
    std::optional<ScaledEnclosure> radius;
    std::uint64_t steps = 0;
};

namespace {
using namespace shared_p1;

struct ActionHolder {
    explicit ActionHolder(core::PlannedBuffer<SharedP1Data> result) noexcept
        : payload(std::move(result)) {}
    core::PlannedBuffer<SharedP1Data> payload;
};

struct ReusingActionHolder : ActionHolder {
    using ActionHolder::ActionHolder;
    ~ReusingActionHolder() {
        if (reuse) {
            auto& action = payload.values()[0];
            reuse->idle = std::move(action.facets);
            reuse->packets = std::move(packets);
            if (action.points.size() != 0)
                reuse->points = std::move(action.points);
            if (action.strengths.size() != 0)
                reuse->strengths = std::move(action.strengths);
            if (action.occurrences.size() != 0) {
                auto& sorted_storage = sorted_into_scratch ? reuse->edge_scratch : reuse->edges;
                sorted_storage = std::move(action.occurrences);
            }
            reuse->checked_out.store(false, std::memory_order_release);
        }
    }
    std::shared_ptr<BoundWakeFacetStorageData> reuse;
    core::PlannedBuffer<std::optional<LocalGeometryPacket>> packets;
    bool sorted_into_scratch = false;
};

[[nodiscard]] std::uint64_t checked_limit(std::uint64_t needed, std::uint64_t available) {
    if (needed > available)
        throw ArithmeticFailure{exceeded};
    return needed;
}
void charge(Work& work, std::uint64_t steps) {
    if (work.stop.stop_requested())
        throw ArithmeticFailure{core::FailureCode::Cancelled};
    work.used += checked_limit(steps, work.limit - work.used);
}

[[nodiscard]] std::uint64_t logarithmic_visits(std::uint64_t count) noexcept {
    return static_cast<std::uint64_t>(std::bit_width(count));
}

[[nodiscard]] std::uint64_t apply_work(std::uint64_t facets, std::uint64_t points,
                                       std::uint64_t strengths, bool vjp) {
    // Per facet: 3 columns x 3 components x 13 geometry variables, with at most
    // six visits each to import and multiply the derivative enclosure, and up to
    // two visits per direction word. The remaining margin covers coefficient
    // terms, complex products/sums and output finishing.
    return sum(sum(product(facets, 1200), vjp ? 32 : 64),
               sum(product(points, 8), product(strengths, 5)));
}

void operator_environment() {
    if (!p1_environment_supported())
        throw ArithmeticFailure{core::FailureCode::PrecisionUnsupported};
}

} // namespace

class SharedP1Engine {
    template <class T>
    static bool same_owner(const std::shared_ptr<const T>& left,
                           const std::shared_ptr<const T>& right) noexcept {
        return left && right && !left.owner_before(right) && !right.owner_before(left);
    }

    static const GlobalContactData& source_data(const P1ActionSource& source) {
        if (source.contact.data_.size() != 1 || !source.query.active())
            throw ArithmeticFailure{invalid};
        const auto& verified = source.contact.data_.values()[0];
        if (!verified.geometry || !verified.topology || !source.query.geometry_)
            throw ArithmeticFailure{invalid};
        if (!same_owner(verified.geometry->data_, source.query.geometry_) ||
            verified.geometry->context() != source.query.context())
            throw ArithmeticFailure{stale};
        return verified;
    }

    static const SharedP1Data& matches(const SharedP1Action& action, const P1ActionSource& source) {
        if (!action.data_ || !action.data_->geometry)
            throw ArithmeticFailure{invalid};
        const auto& expected = source_data(source);
        if (!expected.geometry ||
            !same_owner(expected.geometry->data_, action.data_->geometry->data_) ||
            !same_owner(expected.topology, action.data_->topology) ||
            source.stamp != action.data_->stamp ||
            !same_query(needed(source.query.definition()), action.data_->query))
            throw ArithmeticFailure{stale};
        return *action.data_;
    }

    static void verify_shapes(const SharedP1Data& action, std::span<const P1ComplexWords> strengths,
                              const std::optional<P1ActionDirection>& direction) {
        if (strengths.size() != action.strengths.size())
            throw ArithmeticFailure{invalid};
        if (direction && (direction->points.size() != action.canonical_points ||
                          direction->strengths.size() != action.strengths.size() ||
                          direction->free_query.has_value() == (action.material.count != 0)))
            throw ArithmeticFailure{invalid};
    }

    static P1AssemblyPlan count(const P1ActionSource& source,
                                std::span<const P1FaceStrengths> faces,
                                std::span<const P1StrengthId> strengths,
                                const P1AssemblyLimits& limits, Work& scan, bool mixed) {
        const auto& verified = source_data(source);
        if (!verified.geometry || !verified.geometry->data_)
            throw ArithmeticFailure{invalid};
        const auto& geometry = *verified.geometry->data_;
        const auto& actual = geometry.counts;
        const auto& allowed = limits.geometry.counts;
        if (actual.points > allowed.points || actual.edges > allowed.edges ||
            actual.faces > allowed.faces || actual.bands > allowed.bands ||
            actual.triangles > allowed.triangles || actual.quads > allowed.quads)
            throw ArithmeticFailure{exceeded};
        if (faces.size() != geometry.faces.size())
            throw ArithmeticFailure{invalid};
        scan.step();
        const auto finite_outer = [&]() {
            for (std::size_t i = 0; i < strengths.size(); ++i) {
                scan.step();
                if (i != 0 && strengths[i - 1] >= strengths[i])
                    throw ArithmeticFailure{invalid};
                if (const auto face = strengths[i].bound_face()) {
                    const auto row = ordinal(geometry.face_lookup.values(), face->value(), scan);
                    if (!mixed ||
                        verified.topology->roles.values()[row] != GlobalTopologyData::bound_role)
                        throw ArithmeticFailure{invalid};
                }
            }
            for (std::size_t i = 0; i < faces.size(); ++i) {
                scan.step();
                if (i != 0 && faces[i - 1].face.value() >= faces[i].face.value())
                    throw ArithmeticFailure{invalid};
                const auto face = faces[i].face.value();
                const auto index = ordinal(geometry.face_lookup.values(), face, scan);
                const bool bound = std::holds_alternative<P0BoundColumn>(faces[i].columns);
                if (mixed ? bound != (verified.topology->roles.values()[index] ==
                                      GlobalTopologyData::bound_role)
                          : bound)
                    throw ArithmeticFailure{invalid};
                std::visit(
                    [&](const auto& keys) {
                        if constexpr (std::is_same_v<std::decay_t<decltype(keys)>, P0BoundColumn>)
                            (void)locate(strengths, P1StrengthId::bound(PhysicalFaceId{face}),
                                         std::identity{}, scan);
                        else {
                            for (const auto& key : keys)
                                (void)locate(strengths, key, std::identity{}, scan);
                        }
                    },
                    faces[i].columns);
            }
        };
        finite_outer();
        std::size_t bound_facets = 0;
        for (const auto& record : geometry.facets.values()) {
            scan.step();
            const auto face = geometry.faces.values()[static_cast<std::size_t>(record.face)].id;
            const auto row =
                locate(faces, face, [](const P1FaceStrengths& f) { return f.face.value(); }, scan);
            const auto& columns = faces[row].columns;
            const bool bound = std::holds_alternative<P0BoundColumn>(columns);
            if (std::holds_alternative<QuadOrigin>(record.origin) !=
                (std::holds_alternative<std::array<P1StrengthId, 2>>(columns) || bound))
                throw ArithmeticFailure{invalid};
            bound_facets += static_cast<std::size_t>(bound);
        }
        const auto n = geometry.facets.size();
        const auto e = product(n, 3);
        const auto p = geometry.points.size();
        const auto s = strengths.size();
        P1AssemblyPlan result;
        result.payload_bytes =
            sum(sum(bytes<SharedP1Data>(1), sum(bytes<std::uint64_t>(p), bytes<P1StrengthId>(s))),
                sum(sum(bytes<Facet>(n), bytes<std::optional<LocalGeometryPacket>>(n)),
                    sum(sum(bytes<std::size_t>(p), product(bytes<EdgeOccurrence>(e), 2)),
                        product(n, sum(ExactTriangleResult::payload_bytes(),
                                       sum(LocalGeometryPacket::payload_bytes(),
                                           P1Basis::payload_bytes()))))));
        result.estimated_metadata_bytes =
            sum(predicate_memory_plan().peak_payload_estimate_bytes(),
                sum(sizeof(ActionHolder), product(sum(9, product(3, n)), allocator_estimate)));
        const auto mode = limits.kernel.external_memory.selection();
        const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&mode);
        if (estimate == nullptr)
            throw ArithmeticFailure{core::FailureCode::Unsupported};
        (void)checked_limit(result.estimated_metadata_bytes, estimate->bytes);

        // Checked merge-sort performs one visit per occurrence in each merge pass.
        const auto sort = product(e, logarithmic_visits(e));
        const auto metadata = sum(sum(sum(product(e, mixed ? 9 : 8), sum(product(p, 3), s)), sort),
                                  sum(product(n, sum(30, sum(logarithmic_visits(faces.size()),
                                                             product(3, logarithmic_visits(s))))),
                                      64));
        result.preflight_steps = scan.used;
        // A caller may allow much more work than any local expression executes.
        // Reserve the ceiling the local evaluator actually enforces.
        const auto bound_steps = bound_facets == 0 ? 0 : P1Workspace::maximum_bound_steps();
        const auto wake_steps = bound_facets == n ? 0 : P1Workspace::maximum_steps();
        (void)checked_limit(std::max(bound_steps, wake_steps), limits.kernel.max_steps);
        const auto kernel_steps =
            sum(product(bound_facets, bound_steps), product(n - bound_facets, wake_steps));
        const auto packet_steps = local_packet_step_ceiling(n, limits.geometry);
        result.reserved_steps =
            sum(result.preflight_steps, sum(metadata, sum(product(n, packet_steps), kernel_steps)));
        (void)checked_limit(result.reserved_steps, limits.steps);
        return result;
    }

    static std::size_t populate_points(SharedP1Data& data, const GeometryData& geometry,
                                       std::span<std::size_t> point_columns, Work& work) {
        std::size_t count_points = 0;
        for (const auto& index : geometry.point_lookup.values()) {
            work.step();
            const auto raw = static_cast<std::size_t>(index.ordinal);
            if (raw == static_cast<std::size_t>(geometry.points.values()[raw].owner)) {
                data.points.values()[count_points] = index.id;
                point_columns[raw] = count_points++;
            }
        }
        return count_points;
    }

    static void map_facets(SharedP1Data& data, const GeometryData& geometry,
                           std::span<const P1FaceStrengths> faces,
                           std::span<const P1StrengthId> strengths,
                           std::span<const std::size_t> point_columns,
                           std::span<EdgeOccurrence> occurrences, Work& work) {
        if (!data.geometry || !data.topology)
            throw ArithmeticFailure{invalid};
        for (std::size_t i = 0; i < data.facets.size(); ++i) {
            work.step();
            auto& target = data.facets.values()[i];
            const auto& facet = geometry.facets.values()[i];
            const auto face = geometry.faces.values()[static_cast<std::size_t>(facet.face)].id;
            const auto& map =
                faces[locate(
                          faces, face, [](const P1FaceStrengths& f) { return f.face.value(); },
                          work)]
                    .columns;
            const auto scatter = needed(data.geometry->scatter(FacetId{facet.id}));
            const auto append = [&](const auto& keys, const auto& coefficients) {
                for (std::size_t local = 0; local < 3; ++local) {
                    for (std::size_t col = 0; col < keys.size(); ++col) {
                        work.step();
                        if (coefficients[local * keys.size() + col] == 1) {
                            target.strength[local] =
                                locate(strengths, keys[col], std::identity{}, work);
                            break;
                        }
                    }
                }
            };
            target.bound = std::holds_alternative<P0BoundColumn>(map);
            if (target.bound) {
                target.strength.fill(locate(strengths, P1StrengthId::bound(PhysicalFaceId{face}),
                                            std::identity{}, work));
            } else if (const auto* quad = std::get_if<QuadLRScatter>(&scatter)) {
                const auto* keys = std::get_if<std::array<P1StrengthId, 2>>(&map);
                if (keys == nullptr)
                    throw ArithmeticFailure{invalid};
                append(*keys, quad->coefficients);
            } else {
                const auto* tri = std::get_if<IndependentP1Scatter>(&scatter);
                const auto* keys = std::get_if<std::array<P1StrengthId, 3>>(&map);
                if (tri == nullptr || keys == nullptr)
                    throw ArithmeticFailure{invalid};
                append(*keys, tri->coefficients);
            }
            for (std::size_t edge = 0; edge < 3; ++edge) {
                work.step();
                target.point[edge] = point_columns[static_cast<std::size_t>(facet.points[edge])];
                const auto& binding =
                    data.topology->edges.values()[static_cast<std::size_t>(facet.edges[edge])];
                const auto next = (edge + 1) % 3;
                // Incidence from the geometry factory is already in canonical owner order.
                const bool positive = facet.incidence[edge] == 1;
                occurrences[3 * i + edge] = {
                    static_cast<std::size_t>(binding.owner),
                    positive ? std::array{target.strength[edge], target.strength[next]}
                             : std::array{target.strength[next], target.strength[edge]},
                    i, edge, positive};
            }
        }
    }

public:
    static AeroResult<P1AssemblyPlan>
    forecast(const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
             std::span<const P1StrengthId> strengths, const P1AssemblyLimits& limits,
             const std::stop_token& stop, bool mixed = false) noexcept {
        return guarded<P1AssemblyPlan>(source.query.context(), [&] {
            Work work{0, limits.steps, stop};
            return count(source, faces, strengths, limits, work, mixed);
        });
    }

    static std::uint64_t midpoint_steps(const GeometryData& geometry,
                                        const GlobalContactData& verified) {
        return sum(sum(logarithmic_visits(geometry.edge_lookup.size()),
                       logarithmic_visits(verified.topology->trailing.size())),
                   sum(1, product(geometry.facets.size(), 6)));
    }

    static AeroResult<P1SerialAssemblyPlan>
    forecast_serial(const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
                    std::span<const P1StrengthId> strengths, const P1AssemblyLimits& limits,
                    const std::stop_token& stop, bool midpoint = false) noexcept {
        return guarded<P1SerialAssemblyPlan>(source.query.context(), [&] {
            Work work{0, limits.steps, stop};
            auto ordinary = count(source, faces, strengths, limits, work, true);
            const auto& verified = source_data(source);
            const auto& geometry = *verified.geometry->data_;
            if (midpoint) {
                ordinary.reserved_steps =
                    sum(ordinary.reserved_steps, midpoint_steps(geometry, verified));
                (void)checked_limit(ordinary.reserved_steps, limits.steps);
            }
            const auto facets = geometry.facets.size();
            if (facets == 0)
                return P1SerialAssemblyPlan{ordinary, ordinary.payload_bytes,
                                            ordinary.payload_bytes, ordinary.reserved_steps, true};
            const auto first_payload = sum(
                ordinary.payload_bytes - product(facets - 1, ExactTriangleResult::payload_bytes()),
                bytes<BoundWakeFacetStorageData>(1));
            const auto repeated_payload = bytes<SharedP1Data>(1);
            const auto repeated_steps =
                sum(ordinary.reserved_steps, sum(product(2, facets), geometry.points.size()));
            const auto metadata =
                sum(ordinary.estimated_metadata_bytes,
                    sum(sizeof(ReusingActionHolder) - sizeof(ActionHolder),
                        sum(sizeof(core::PlannedBuffer<BoundWakeFacetStorageData>),
                            allocator_estimate)));
            const auto mode = limits.kernel.external_memory.selection();
            const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&mode);
            return P1SerialAssemblyPlan{ordinary, first_payload, repeated_payload, repeated_steps,
                                        estimate != nullptr && metadata <= estimate->bytes &&
                                            repeated_steps <= limits.steps};
        });
    }

    static AeroResult<SharedP1Action>
    prepare(const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
            std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy,
            P1Accuracy accuracy, const P1AssemblyLimits& limits, SharedP1Workspace& workspace,
            LocalPacketWorkspace& local, P1Workspace& kernel, core::ResourcePlan& memory,
            const std::stop_token& stop, bool mixed = false, std::optional<EdgeId> self = {},
            BoundWakeFacetStorage* reuse = nullptr, bool values_only = false) noexcept {
        const auto context = source.query.context();
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(invalid, context));
        const Lease lease{workspace.entered_};
        return guarded<SharedP1Action>(context, [&] {
            operator_environment();
            if (const auto error = local_accuracy.validation_error())
                throw ArithmeticFailure{*error};
            if (const auto error = accuracy.validation_error())
                throw ArithmeticFailure{*error};
            Work work{0, limits.steps, stop};
            auto plan = count(source, faces, strengths, limits, work, mixed);
            const auto& verified = source_data(source);
            if (!verified.geometry || !verified.geometry->data_)
                throw ArithmeticFailure{invalid};
            const auto& geometry = *verified.geometry->data_;
            const auto facet_count = geometry.facets.size();
            const auto point_slots = geometry.points.size();
            const auto edges = product(facet_count, 3);
            const bool use_storage = reuse != nullptr && facet_count != 0;
            if (reuse != nullptr && reuse->data_ &&
                reuse->data_->checked_out.load(std::memory_order_acquire))
                throw ArithmeticFailure{invalid};
            const bool has_idle =
                reuse != nullptr && reuse->data_ && reuse->data_->idle.size() != 0;
            const bool has_packets = has_idle && reuse->data_->packets.size() != 0;
            if (has_idle) {
                if (reuse->data_->idle.size() != facet_count ||
                    (has_packets && reuse->data_->packets.size() != facet_count))
                    throw ArithmeticFailure{invalid};
                plan.reserved_steps = sum(plan.reserved_steps, product(2, facet_count));
                (void)checked_limit(plan.reserved_steps, limits.steps);
                plan.payload_bytes -= bytes<Facet>(facet_count);
                if (has_packets)
                    plan.payload_bytes -= bytes<std::optional<LocalGeometryPacket>>(facet_count);
                for (std::size_t i = 0; i < facet_count; ++i) {
                    work.step();
                    if (has_packets && reuse->data_->packets.values()[i])
                        plan.payload_bytes -= LocalGeometryPacket::payload_bytes();
                    if (reuse->data_->idle.values()[i].basis)
                        plan.payload_bytes -= P1Basis::payload_bytes();
                }
            }
            const bool has_predicate =
                use_storage && reuse->data_ && reuse->data_->predicate.size() != 0;
            const bool has_point_columns =
                use_storage && reuse->data_ && reuse->data_->point_columns.size() == point_slots;
            const bool has_points =
                use_storage && reuse->data_ && reuse->data_->points.size() == point_slots;
            const bool has_strengths =
                use_storage && reuse->data_ && reuse->data_->strengths.size() == strengths.size();
            if (use_storage)
                plan.payload_bytes -= product(facet_count - (has_predicate ? 0U : 1U),
                                              ExactTriangleResult::payload_bytes());
            if (use_storage && reuse->data_) {
                const auto& stored = *reuse->data_;
                if (has_point_columns) {
                    plan.payload_bytes -= bytes<std::size_t>(point_slots);
                    plan.reserved_steps = sum(plan.reserved_steps, point_slots);
                    (void)checked_limit(plan.reserved_steps, limits.steps);
                }
                if (stored.edges.size() == edges)
                    plan.payload_bytes -= bytes<EdgeOccurrence>(edges);
                if (stored.edge_scratch.size() == edges)
                    plan.payload_bytes -= bytes<EdgeOccurrence>(edges);
                if (has_points)
                    plan.payload_bytes -= bytes<std::uint64_t>(point_slots);
                if (has_strengths)
                    plan.payload_bytes -= bytes<P1StrengthId>(strengths.size());
            }
            if (use_storage && !reuse->data_) {
                plan.payload_bytes = sum(plan.payload_bytes, bytes<BoundWakeFacetStorageData>(1));
                plan.estimated_metadata_bytes =
                    sum(plan.estimated_metadata_bytes,
                        sum(sizeof(core::PlannedBuffer<BoundWakeFacetStorageData>),
                            allocator_estimate));
            }
            if (use_storage) {
                plan.estimated_metadata_bytes =
                    sum(plan.estimated_metadata_bytes,
                        sizeof(ReusingActionHolder) - sizeof(ActionHolder));
                const auto mode = limits.kernel.external_memory.selection();
                const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&mode);
                if (estimate == nullptr)
                    throw ArithmeticFailure{core::FailureCode::Unsupported};
                (void)checked_limit(plan.estimated_metadata_bytes, estimate->bytes);
            }
            (void)checked_limit(plan.payload_bytes, memory.remaining_bytes());
            std::optional<std::size_t> excluded;
            if (self) {
                plan.reserved_steps = sum(plan.reserved_steps, midpoint_steps(geometry, verified));
                (void)checked_limit(plan.reserved_steps, limits.steps);
                work.step();
                const auto definition = needed(source.query.definition());
                const auto* edge = std::get_if<SourceEdgePoint>(&definition.definition);
                const auto* fraction =
                    edge != nullptr ? std::get_if<SingleWord>(&edge->t) : nullptr;
                if (edge == nullptr || fraction == nullptr || edge->edge != *self ||
                    fraction->value != .5)
                    throw ArithmeticFailure{invalid};
                const auto lookup = geometry.edge_lookup.values();
                const auto index = locate(lookup, self->value(), &IdLookup::id, work);
                const auto ordinal = static_cast<std::size_t>(lookup[index].ordinal);
                if (verified.topology->edges.values()[ordinal].owner != ordinal)
                    throw ArithmeticFailure{invalid};
                const auto trailing = verified.topology->trailing.values();
                const auto te = std::lower_bound(trailing.begin(), trailing.end(), *self,
                                                 [&](const TrailingEdgeRow& row, EdgeId id) {
                                                     work.step();
                                                     return row.edge < id;
                                                 });
                if (te != trailing.end() && te->edge == *self)
                    throw ArithmeticFailure{invalid};
                excluded = ordinal;
            }
            if (use_storage && !reuse->data_) {
                auto storage = allocate<BoundWakeFacetStorageData>(memory, 1, context, stop);
                auto pool = std::make_shared<core::PlannedBuffer<BoundWakeFacetStorageData>>(
                    std::move(storage));
                auto* state = pool->values().data();
                reuse->data_ = {std::move(pool), state};
            }
            if (use_storage && !has_predicate)
                reuse->data_->predicate = allocate<TriangleOutcome>(memory, 1, context, stop);
            auto* predicate = use_storage ? reuse->data_->predicate.values().data() : nullptr;
            auto storage = allocate<SharedP1Data>(memory, 1, context, stop);
            std::shared_ptr<ActionHolder> holder;
            ReusingActionHolder* reusable_holder = nullptr;
            if (use_storage) {
                auto reusable = std::make_shared<ReusingActionHolder>(std::move(storage));
                reusable_holder = reusable.get();
                holder = std::move(reusable);
            } else {
                holder = std::make_shared<ActionHolder>(std::move(storage));
            }
            auto& data = holder->payload.values()[0];
            data.facets = has_idle ? std::move(reuse->data_->idle)
                                   : allocate<Facet>(memory, facet_count, context, stop);
            if (reusable_holder) {
                reusable_holder->reuse = reuse->data_;
                reuse->data_->checked_out.store(true, std::memory_order_relaxed);
            }
            data.points = has_points ? std::move(reuse->data_->points)
                                     : allocate<std::uint64_t>(memory, point_slots, context, stop);
            data.strengths = has_strengths
                                 ? std::move(reuse->data_->strengths)
                                 : allocate<P1StrengthId>(memory, strengths.size(), context, stop);
            core::PlannedBuffer<std::optional<LocalGeometryPacket>> ordinary_packets;
            auto& packet_storage =
                reusable_holder != nullptr ? reusable_holder->packets : ordinary_packets;
            packet_storage = has_packets ? std::move(reuse->data_->packets)
                                         : allocate<std::optional<LocalGeometryPacket>>(
                                               memory, facet_count, context, stop);
            auto packets = packet_storage.values();
            core::PlannedBuffer<std::size_t> local_point_columns;
            core::PlannedBuffer<EdgeOccurrence> local_edges, local_edge_scratch;
            auto& point_columns = use_storage ? reuse->data_->point_columns : local_point_columns;
            auto& occurrences = use_storage ? reuse->data_->edges : local_edges;
            auto& sorting = use_storage ? reuse->data_->edge_scratch : local_edge_scratch;
            const auto fit = [&]<class T>(core::PlannedBuffer<T>& buffer, std::uint64_t count) {
                if (buffer.size() != count)
                    buffer = allocate<T>(memory, count, context, stop);
            };
            fit(point_columns, point_slots);
            fit(occurrences, edges);
            fit(sorting, edges);
            if (has_point_columns) {
                charge(work, point_slots);
                std::ranges::fill(point_columns.values(), std::size_t{0});
            }
            operator_environment();
            if (has_idle) {
                for (auto& facet : data.facets.values()) {
                    work.step();
                    auto basis = std::move(facet.basis);
                    facet = Facet{};
                    facet.basis = std::move(basis);
                }
            }
            data.query = needed(source.query.definition());
            data.stamp = source.stamp;
            data.topology = verified.topology;
            data.geometry.emplace(OwnedGeometry{verified.geometry->data_, context});
            data.precision = accuracy.representation;
            data.values_only = values_only;
            data.report.plan = plan;
            data.report.facets = facet_count;
            data.canonical_points = populate_points(data, geometry, point_columns.values(), work);
            data.material = material(geometry, data.query, point_columns.values(), work);
            for (std::size_t s = 0; s < strengths.size(); ++s) {
                work.step();
                data.strengths.values()[s] = strengths[s];
            }
            map_facets(data, geometry, faces, strengths, point_columns.values(),
                       occurrences.values(), work);
            auto packet_limits = limits.geometry;
            packet_limits.max_steps = local_packet_step_ceiling(facet_count, packet_limits);
            for (std::size_t i = 0; i < facet_count; ++i) {
                auto& facet = data.facets.values()[i];
                const auto prepare_packet =
                    facet.bound ? prepare_bound_local_packet : prepare_local_packet;
                const auto& original = geometry.facets.values()[i];
                bool has_alias = false;
                if (excluded) {
                    for (const auto edge : original.edges) {
                        work.step();
                        const auto raw = static_cast<std::size_t>(edge);
                        if (raw != *excluded &&
                            data.topology->edges.values()[raw].owner == *excluded)
                            has_alias = true;
                    }
                }
                const auto id = FacetId{original.id};
                SharedP1ReuseSlot<LocalGeometryPacket> slot{packets[i], predicate};
                auto* storage_slot = use_storage ? &slot : nullptr;
                auto packet = needed(
                    has_alias ? prepare_aliased_midpoint_packet(
                                    source.contact, source.query, id, facet.bound, local_accuracy,
                                    packet_limits, local, memory, stop, storage_slot)
                              : prepare_packet(*verified.geometry, source.query, id, local_accuracy,
                                               packet_limits, local, memory, stop, storage_slot));
                const auto steps = needed(packet.steps());
                const bool far = !facet.bound && needed(packet.far()) == FarRelation::Far;
                packets[i].emplace(std::move(packet));
                charge(work, steps);
                facet.far = far;
                if (far)
                    ++data.report.far_facets;
            }
            const auto ordered = sort_edges(occurrences.values(), sorting.values(), work);
            data.report.cancelled_analytic_pairs =
                cancel_pairs(ordered, data.facets.values(), work);
            if (excluded) {
                for (const auto& edge : ordered) {
                    work.step();
                    auto& facet = data.facets.values()[edge.facet];
                    if (edge.owner == *excluded && facet.bound) {
                        facet.boundary[edge.edge] = false;
                        ++data.report.excluded_bound_occurrences;
                    }
                }
                if (data.report.excluded_bound_occurrences == 0)
                    throw ArithmeticFailure{invalid};
            }
            if (mixed) {
                group_bound_edges(ordered, data.facets.values(), work);
                const bool sorted_into_scratch = ordered.data() != occurrences.values().data();
                data.occurrences =
                    sorted_into_scratch ? std::move(sorting) : std::move(occurrences);
                if (reusable_holder)
                    reusable_holder->sorted_into_scratch = sorted_into_scratch;
            }
            for (std::size_t i = 0; i < facet_count; ++i) {
                auto& entry = data.facets.values()[i];
                const auto& packet = packets[i];
                if (!packet)
                    throw ArithmeticFailure{invalid};
                SharedP1ReuseSlot<P1Basis> slot{entry.basis};
                auto* storage_slot = use_storage ? &slot : nullptr;
                auto basis = needed([&] {
                    if (values_only)
                        return evaluate_partition_values(*packet, entry.boundary, entry.bound,
                                                         accuracy, limits.kernel, kernel, memory,
                                                         stop, storage_slot);
                    if (entry.bound)
                        return evaluate_bound_partition(*packet, entry.boundary, accuracy,
                                                        limits.kernel, kernel, memory, stop,
                                                        storage_slot);
                    return evaluate_p1_partition(*packet, entry.boundary, accuracy, limits.kernel,
                                                 kernel, memory, stop, storage_slot);
                }());
                const auto steps = needed(basis.steps());
                entry.basis.emplace(std::move(basis));
                charge(work, steps);
            }
            work.step();
            data.report.used_steps = work.used;
            return SharedP1Action{std::shared_ptr<const SharedP1Data>{std::move(holder), &data},
                                  context};
        });
    }

    static AeroResult<SharedP1Action> prepare_reusing_facets(
        const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
        std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy, P1Accuracy accuracy,
        const P1AssemblyLimits& limits, SharedP1Workspace& workspace, LocalPacketWorkspace& local,
        P1Workspace& kernel, BoundWakeFacetStorage& reuse, const std::stop_token& stop,
        std::optional<EdgeId> self = {}, bool values_only = false) noexcept {
        if (reuse.preparing_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(invalid, source.query.context()));
        const Lease lease{reuse.preparing_};
        return prepare(source, faces, strengths, local_accuracy, accuracy, limits, workspace, local,
                       kernel, reuse.memory_, stop, true, self, &reuse, values_only);
    }

    static AeroResult<P1ActionValue> apply(const SharedP1Action& action,
                                           const P1ActionSource& source,
                                           std::span<const P1ComplexWords> strengths,
                                           const std::optional<P1ActionDirection>& direction,
                                           std::uint64_t steps, SharedP1Workspace& workspace,
                                           const std::stop_token& stop) noexcept {
        const auto context = source.query.context();
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(invalid, context));
        const Lease lease{workspace.entered_};
        return guarded<P1ActionValue>(context, [&] {
            operator_environment();
            const auto& data = matches(action, source);
            if (data.values_only && direction)
                throw ArithmeticFailure{invalid};
            verify_shapes(data, strengths, direction);
            (void)checked_limit(
                apply_work(data.facets.size(), data.canonical_points, data.strengths.size(), false),
                steps);
            Work work{0, steps, stop};
            validate_words(strengths, direction, work);
            std::array<Complex, 3> primal{}, jvp{};
            const auto delta_query = query_direction(data.material, direction, work);
            for (const auto& facet : data.facets.values()) {
                if (!facet.basis)
                    throw ArithmeticFailure{invalid};
                for (std::size_t k = 0; k < 3; ++k) {
                    if (facet.bound && !facet.boundary[k])
                        continue;
                    const auto mu = coefficient(data, facet, k, strengths, work);
                    const auto delta_mu =
                        direction ? coefficient(data, facet, k, direction->strengths, work)
                                  : Complex{};
                    for (std::size_t a = 0; a < 3; ++a) {
                        const auto value = enclosure(needed(facet.basis->basis(k, a)).get(), work);
                        primal[a] = add(primal[a], multiply(mu, value, work), work);
                        if (direction) {
                            const auto delta =
                                basis_direction(facet, k, a, delta_query, *direction, work);
                            const auto term = add(multiply(delta_mu, value, work),
                                                  multiply(mu, delta, work), work);
                            jvp[a] = add(jvp[a], term, work);
                        }
                    }
                }
            }
            LocalRootScratch rounding;
            charge(work, 6);
            const auto velocity = finish(primal, data.precision, rounding);
            std::optional<std::array<P1ComplexEnclosure, 3>> derivative;
            if (direction) {
                charge(work, 6);
                derivative.emplace(finish(jvp, data.precision, rounding));
            }
            work.step();
            return P1ActionValue{action.data_, velocity, derivative, work.used};
        });
    }

    static AeroResult<std::uint64_t>
    strength_transpose(const SharedP1Action& action, const P1ActionSource& source,
                       std::span<Interval> output, std::uint64_t steps,
                       SharedP1Workspace& workspace, const std::stop_token& stop) noexcept {
        const auto context = source.query.context();
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(invalid, context));
        const Lease lease{workspace.entered_};
        return guarded<std::uint64_t>(context, [&] {
            operator_environment();
            const auto& data = matches(action, source);
            constexpr std::size_t axes = 3;
            if (output.size() != product(data.strengths.size(), axes))
                throw ArithmeticFailure{invalid};
            Work work{0, steps, stop};
            for (auto& value : output) {
                work.step();
                value = {};
            }
            for (const auto& facet : data.facets.values()) {
                if (!facet.basis)
                    throw ArithmeticFailure{invalid};
                for (std::size_t k = 0; k < axes; ++k) {
                    if (facet.bound && !facet.boundary[k])
                        continue;
                    for (std::size_t axis = 0; axis < axes; ++axis) {
                        const auto value =
                            enclosure(needed(facet.basis->basis(k, axis)).get(), work);
                        visit_coefficients(data, facet, k, work, [&](std::size_t id, bool forward) {
                            auto& entry = output[axes * id + axis];
                            entry = forward ? p1::add(entry, value, work)
                                            : p1::subtract(entry, value, work);
                        });
                    }
                }
            }
            LocalRootScratch rounding;
            for (auto& value : output) {
                work.step();
                value = enclosure(enclose_p1_value(value, data.precision, rounding), work);
            }
            work.step();
            return work.used;
        });
    }

    static AeroResult<P1ActionAdjoint>
    adjoint(const SharedP1Action& action, const P1ActionSource& source,
            std::span<const P1ComplexWords> strengths, const std::array<P1ComplexWords, 3>& dual,
            std::uint64_t steps, SharedP1Workspace& workspace, core::ResourcePlan& memory,
            const std::stop_token& stop) noexcept {
        const auto context = source.query.context();
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(invalid, context));
        const Lease lease{workspace.entered_};
        return guarded<P1ActionAdjoint>(context, [&] {
            operator_environment();
            const auto& data = matches(action, source);
            if (data.values_only)
                throw ArithmeticFailure{invalid};
            verify_shapes(data, strengths, std::nullopt);
            const auto p = data.canonical_points;
            const auto s = data.strengths.size();
            (void)checked_limit(apply_work(data.facets.size(), p, s, true), steps);
            const auto payload =
                sum(bytes<SharedP1AdjointData>(1),
                    sum(sum(bytes<PointGradient>(p), bytes<CoordinateResult>(p)),
                        sum(bytes<Complex>(s), bytes<std::optional<P1ComplexEnclosure>>(s))));
            (void)checked_limit(payload, memory.remaining_bytes());
            Work work{0, steps, stop};
            validate_words(strengths, std::nullopt, work);
            const std::array<Complex, 3> lambda{words(dual[0], work), words(dual[1], work),
                                                words(dual[2], work)};
            auto owned = allocate<SharedP1AdjointData>(memory, 1, context, stop);
            auto point_gradient = allocate<PointGradient>(memory, p, context, stop);
            auto strength_gradient = allocate<Complex>(memory, s, context, stop);
            auto& result = owned.values()[0];
            result.parent = action.data_;
            result.points = allocate<CoordinateResult>(memory, p, context, stop);
            result.strengths =
                allocate<std::optional<P1ComplexEnclosure>>(memory, s, context, stop);
            operator_environment();
            std::array<Interval, 3> query{};
            Interval radius{};
            for (const auto& facet : data.facets.values()) {
                if (!facet.basis)
                    throw ArithmeticFailure{invalid};
                for (std::size_t k = 0; k < 3; ++k) {
                    if (facet.bound && !facet.boundary[k])
                        continue;
                    const auto mu = coefficient(data, facet, k, strengths, work);
                    for (std::size_t a = 0; a < 3; ++a) {
                        const auto value = enclosure(needed(facet.basis->basis(k, a)).get(), work);
                        const auto gradient = multiply(lambda[a], value, work);
                        visit_coefficients(data, facet, k, work, [&](std::size_t id, bool forward) {
                            auto& target = strength_gradient.values()[id];
                            target = add(target, oriented(gradient, forward, work), work);
                        });
                        const auto factor = dot(lambda[a], mu, work);
                        const auto scatter = [&](Interval& target, std::size_t variable) {
                            const auto derivative = enclosure(
                                needed(facet.basis->derivative(k, a, variable)).get(), work);
                            target = p1::add(target, p1::multiply(factor, derivative, work), work);
                        };
                        for (std::size_t axis = 0; axis < 3; ++axis) {
                            scatter(query[axis], axis);
                            for (std::size_t point = 0; point < 3; ++point)
                                scatter(point_gradient.values()[facet.point[point]].value[axis],
                                        3 + 3 * point + axis);
                        }
                        scatter(radius, 12);
                    }
                }
            }
            LocalRootScratch rounding;
            if (data.material.count == 0) {
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    work.step();
                    result.query.value[axis] =
                        enclose_p1_value(query[axis], data.precision, rounding);
                }
            } else {
                for (std::size_t i = 0; i < data.material.count; ++i) {
                    auto& target = point_gradient.values()[data.material.point[i]];
                    for (std::size_t axis = 0; axis < 3; ++axis)
                        target.value[axis] =
                            p1::add(target.value[axis],
                                    p1::multiply(data.material.weight[i], query[axis], work), work);
                }
            }
            for (std::size_t row = 0; row < p; ++row)
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    work.step();
                    result.points.values()[row].value[axis] = enclose_p1_value(
                        point_gradient.values()[row].value[axis], data.precision, rounding);
                }
            for (std::size_t i = 0; i < s; ++i) {
                work.step();
                result.strengths.values()[i] =
                    finish(strength_gradient.values()[i], data.precision, rounding);
            }
            result.radius = enclose_p1_value(radius, data.precision, rounding);
            work.step();
            result.steps = work.used;
            return P1ActionAdjoint{std::move(owned), context};
        });
    }

private:
    static void validate_words(std::span<const P1ComplexWords> strengths,
                               const std::optional<P1ActionDirection>& direction, Work& work) {
        const auto check_word = [&](const ScalarWords& w) {
            work.step();
            const auto finite = [](const auto& selected) {
                using T = std::decay_t<decltype(selected)>;
                if constexpr (std::is_same_v<T, SingleWord>)
                    return std::isfinite(selected.value);
                else
                    return std::isfinite(selected.high) && std::isfinite(selected.low);
            };
            if (!std::visit(finite, w))
                throw ArithmeticFailure{core::FailureCode::Nonfinite};
        };
        for (const auto& strength : strengths) {
            check_word(strength.real);
            check_word(strength.imaginary);
        }
        if (!direction)
            return;
        for (const auto& point : direction->points)
            for (const auto& axis : point)
                check_word(axis);
        for (const auto& strength : direction->strengths) {
            check_word(strength.real);
            check_word(strength.imaginary);
        }
        if (direction->free_query)
            for (const auto& axis : *direction->free_query)
                check_word(axis);
        check_word(direction->radius);
    }

    static std::array<Interval, 3>
    query_direction(const MaterialPoint& query, const std::optional<P1ActionDirection>& direction,
                    Work& work) {
        std::array<Interval, 3> result{};
        if (!direction)
            return result;
        if (direction->free_query) {
            for (std::size_t a = 0; a < 3; ++a)
                result[a] = words((*direction->free_query)[a], work);
            return result;
        }
        for (std::size_t i = 0; i < query.count; ++i)
            for (std::size_t a = 0; a < 3; ++a) {
                const auto value = words(direction->points[query.point[i]][a], work);
                result[a] = p1::add(result[a], p1::multiply(value, query.weight[i], work), work);
            }
        return result;
    }

    static Interval basis_direction(const Facet& facet, std::size_t basis, std::size_t axis,
                                    const std::array<Interval, 3>& query,
                                    const P1ActionDirection& direction, Work& work) {
        if (!facet.basis)
            throw ArithmeticFailure{invalid};
        const auto& owner = *facet.basis;
        Interval result{};
        const auto append = [&](std::size_t variable, const Interval& delta) {
            const auto jacobian =
                enclosure(needed(owner.derivative(basis, axis, variable)).get(), work);
            result = p1::add(result, p1::multiply(jacobian, delta, work), work);
        };
        for (std::size_t a = 0; a < 3; ++a) {
            append(a, query[a]);
            for (std::size_t point = 0; point < 3; ++point)
                append(3 + 3 * point + a, words(direction.points[facet.point[point]][a], work));
        }
        append(12, words(direction.radius, work));
        return result;
    }
};

bool SharedP1Action::active() const noexcept {
    return data_ != nullptr;
}
std::size_t SharedP1Action::point_count() const noexcept {
    return data_ ? data_->canonical_points : 0;
}
std::size_t SharedP1Action::strength_count() const noexcept {
    return data_ ? data_->strengths.size() : 0;
}
AeroResult<PointId> SharedP1Action::point(std::size_t index) const noexcept {
    if (!data_ || index >= data_->canonical_points)
        return std::unexpected(failure(invalid, context_));
    return PointId{data_->points.values()[index]};
}
AeroResult<P1StrengthId> SharedP1Action::strength(std::size_t index) const noexcept {
    if (!data_ || index >= data_->strengths.size())
        return std::unexpected(failure(invalid, context_));
    return data_->strengths.values()[index];
}
AeroResult<P1AssemblyReport> SharedP1Action::report() const noexcept {
    if (!data_)
        return std::unexpected(failure(invalid, context_));
    return data_->report;
}
AeroResult<P1EvaluationDiagnostic>
SharedP1Action::facet_diagnostic(std::size_t index) const noexcept {
    if (!data_ || index >= data_->facets.size())
        return std::unexpected(failure(invalid, context_));
    const auto& facet = data_->facets.values()[index];
    if (!facet.basis)
        return std::unexpected(failure(invalid, context_));
    return facet.basis->diagnostic();
}

P1ActionValue::P1ActionValue(std::shared_ptr<const SharedP1Data> owner,
                             const std::array<P1ComplexEnclosure, 3>& velocity,
                             const std::optional<std::array<P1ComplexEnclosure, 3>>& jvp,
                             std::uint64_t steps) noexcept
    : parent_(std::move(owner)), velocity_(velocity), jvp_(jvp), steps_(steps) {}
const std::array<P1ComplexEnclosure, 3>& P1ActionValue::velocity() const& noexcept {
    return velocity_;
}
const std::optional<std::array<P1ComplexEnclosure, 3>>& P1ActionValue::jvp() const& noexcept {
    return jvp_;
}
std::uint64_t P1ActionValue::used_steps() const noexcept {
    return steps_;
}

P1ActionAdjoint::P1ActionAdjoint(core::PlannedBuffer<SharedP1AdjointData> storage,
                                 GeometryContext context) noexcept
    : data_(std::move(storage)), context_(context) {}
P1ActionAdjoint::P1ActionAdjoint(P1ActionAdjoint&&) noexcept = default;
P1ActionAdjoint& P1ActionAdjoint::operator=(P1ActionAdjoint&&) noexcept = default;
P1ActionAdjoint::~P1ActionAdjoint() = default;
AeroResult<ScaledEnclosure> P1ActionAdjoint::point(std::size_t i, std::size_t axis) const noexcept {
    if (data_.size() != 1 || axis >= 3 || i >= data_.values()[0].points.size())
        return std::unexpected(failure(invalid, context_));
    const auto& value = data_.values()[0].points.values()[i].value[axis];
    if (!value)
        return std::unexpected(failure(invalid, context_));
    return *value;
}
AeroResult<P1ComplexEnclosure> P1ActionAdjoint::strength(std::size_t i) const noexcept {
    if (data_.size() != 1 || i >= data_.values()[0].strengths.size())
        return std::unexpected(failure(invalid, context_));
    const auto& value = data_.values()[0].strengths.values()[i];
    if (!value)
        return std::unexpected(failure(invalid, context_));
    return *value;
}
AeroResult<std::optional<ScaledEnclosure>>
P1ActionAdjoint::free_query(std::size_t axis) const noexcept {
    if (data_.size() != 1 || axis >= 3)
        return std::unexpected(failure(invalid, context_));
    return data_.values()[0].query.value[axis];
}
AeroResult<ScaledEnclosure> P1ActionAdjoint::radius() const noexcept {
    if (data_.size() != 1)
        return std::unexpected(failure(invalid, context_));
    const auto& value = data_.values()[0].radius;
    if (!value)
        return std::unexpected(failure(invalid, context_));
    return *value;
}
std::uint64_t P1ActionAdjoint::used_steps() const noexcept {
    return data_.size() == 1 ? data_.values()[0].steps : 0;
}

AeroResult<P1AssemblyPlan> forecast_shared_p1(const P1ActionSource& source,
                                              std::span<const P1FaceStrengths> faces,
                                              std::span<const P1StrengthId> strengths,
                                              const P1AssemblyLimits& limits,
                                              const std::stop_token& stop) noexcept {
    return SharedP1Engine::forecast(source, faces, strengths, limits, stop);
}
AeroResult<SharedP1Action>
prepare_shared_p1(const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
                  std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy,
                  P1Accuracy accuracy, const P1AssemblyLimits& limits, SharedP1Workspace& workspace,
                  LocalPacketWorkspace& local, P1Workspace& kernel, core::ResourcePlan& memory,
                  const std::stop_token& stop) noexcept {
    return SharedP1Engine::prepare(source, faces, strengths, local_accuracy, accuracy, limits,
                                   workspace, local, kernel, memory, stop);
}
AeroResult<P1AssemblyPlan> forecast_bound_wake_induction(const P1ActionSource& source,
                                                         std::span<const P1FaceStrengths> faces,
                                                         std::span<const P1StrengthId> strengths,
                                                         const P1AssemblyLimits& limits,
                                                         const std::stop_token& stop) noexcept {
    return SharedP1Engine::forecast(source, faces, strengths, limits, stop, true);
}

AeroResult<SharedP1Action> prepare_bound_wake_induction(
    const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
    std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy, P1Accuracy accuracy,
    const P1AssemblyLimits& limits, SharedP1Workspace& workspace, LocalPacketWorkspace& local,
    P1Workspace& kernel, core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return SharedP1Engine::prepare(source, faces, strengths, local_accuracy, accuracy, limits,
                                   workspace, local, kernel, memory, stop, true);
}

AeroResult<SharedP1Action> prepare_bound_wake_induction_reusing_facets(
    const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
    std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy, P1Accuracy accuracy,
    const P1AssemblyLimits& limits, SharedP1Workspace& workspace, LocalPacketWorkspace& local,
    P1Workspace& kernel, BoundWakeFacetStorage& reuse, const std::stop_token& stop) noexcept {
    return SharedP1Engine::prepare_reusing_facets(source, faces, strengths, local_accuracy,
                                                  accuracy, limits, workspace, local, kernel, reuse,
                                                  stop);
}

AeroResult<P1SerialAssemblyPlan> forecast_bound_wake_serial_induction(
    const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
    std::span<const P1StrengthId> strengths, const P1AssemblyLimits& limits,
    const std::stop_token& stop) noexcept {
    return SharedP1Engine::forecast_serial(source, faces, strengths, limits, stop);
}

AeroResult<SharedP1Action>
prepare_bound_wake_values(const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
                          std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy,
                          P1Accuracy accuracy, const P1AssemblyLimits& limits,
                          SharedP1Workspace& workspace, LocalPacketWorkspace& local,
                          P1Workspace& kernel, core::ResourcePlan& memory,
                          BoundWakeFacetStorage* reuse, const std::stop_token& stop) noexcept {
    if (reuse != nullptr)
        return SharedP1Engine::prepare_reusing_facets(source, faces, strengths, local_accuracy,
                                                      accuracy, limits, workspace, local, kernel,
                                                      *reuse, stop, {}, true);
    return SharedP1Engine::prepare(source, faces, strengths, local_accuracy, accuracy, limits,
                                   workspace, local, kernel, memory, stop, true, {}, nullptr, true);
}

AeroResult<std::uint64_t> project_real_strength_transpose(
    const SharedP1Action& action, const P1ActionSource& source, std::span<p1::Interval> output,
    std::uint64_t steps, SharedP1Workspace& workspace, const std::stop_token& stop) noexcept {
    return SharedP1Engine::strength_transpose(action, source, output, steps, workspace, stop);
}

AeroResult<SharedP1Action> prepare_bound_wake_midpoint_induction(
    const P1ActionSource& source, EdgeId self, std::span<const P1FaceStrengths> faces,
    std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy, P1Accuracy accuracy,
    const P1AssemblyLimits& limits, SharedP1Workspace& workspace, LocalPacketWorkspace& local,
    P1Workspace& kernel, core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    return SharedP1Engine::prepare(source, faces, strengths, local_accuracy, accuracy, limits,
                                   workspace, local, kernel, memory, stop, true, self);
}

AeroResult<P1SerialAssemblyPlan> forecast_bound_wake_midpoint_serial_induction(
    const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
    std::span<const P1StrengthId> strengths, const P1AssemblyLimits& limits,
    const std::stop_token& stop) noexcept {
    return SharedP1Engine::forecast_serial(source, faces, strengths, limits, stop, true);
}

AeroResult<SharedP1Action> prepare_bound_wake_midpoint_induction(
    const P1ActionSource& source, EdgeId self, std::span<const P1FaceStrengths> faces,
    std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy, P1Accuracy accuracy,
    const P1AssemblyLimits& limits, SharedP1Workspace& workspace, LocalPacketWorkspace& local,
    P1Workspace& kernel, BoundWakeFacetStorage& storage, const std::stop_token& stop) noexcept {
    return SharedP1Engine::prepare_reusing_facets(source, faces, strengths, local_accuracy,
                                                  accuracy, limits, workspace, local, kernel,
                                                  storage, stop, self);
}

AeroResult<P1ActionValue> apply_shared_p1(const SharedP1Action& action,
                                          const P1ActionSource& source,
                                          std::span<const P1ComplexWords> strengths,
                                          const std::optional<P1ActionDirection>& direction,
                                          std::uint64_t max_steps, SharedP1Workspace& workspace,
                                          const std::stop_token& stop) noexcept {
    return SharedP1Engine::apply(action, source, strengths, direction, max_steps, workspace, stop);
}
AeroResult<P1ActionAdjoint> adjoint_shared_p1(const SharedP1Action& action,
                                              const P1ActionSource& source,
                                              std::span<const P1ComplexWords> strengths,
                                              const std::array<P1ComplexWords, 3>& dual,
                                              std::uint64_t max_steps, SharedP1Workspace& workspace,
                                              core::ResourcePlan& memory,
                                              const std::stop_token& stop) noexcept {
    return SharedP1Engine::adjoint(action, source, strengths, dual, max_steps, workspace, memory,
                                   stop);
}

} // namespace uvlas::aero::detail
