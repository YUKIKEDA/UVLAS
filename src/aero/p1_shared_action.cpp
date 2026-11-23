#include "p1_shared_action.hpp"

#include "geometry_snapshot_storage.hpp"
#include "global_geometry_storage.hpp"
#include "local_geometry_bounds_private.hpp"
#include "p1_interval_private.hpp"
#include "p1_shared_kernel.hpp"

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
struct Facet {
    std::array<std::size_t, 3> point{}, strength{};
    std::array<bool, 3> boundary{true, true, true};
    std::optional<P1Basis> basis;
    bool far = false;
};
struct EdgeOccurrence {
    std::size_t owner = 0;
    std::array<std::size_t, 2> trace{};
    std::size_t facet = 0, edge = 0;
    bool positive = false;
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

template <class T, class Key>
[[nodiscard]] std::size_t locate(std::span<const T> sorted, std::uint64_t id, Key key, Work& work) {
    const auto position =
        std::lower_bound(sorted.begin(), sorted.end(), id, [&](const T& item, std::uint64_t value) {
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
} // namespace shared_p1

struct SharedP1Data {
    std::shared_ptr<const GlobalTopologyData> topology;
    std::optional<OwnedGeometry> geometry;
    QueryRequest query{QueryId{0}, FreePoint{FreeSingle{}}};
    P1ActionStamp stamp;
    shared_p1::MaterialPoint material;
    core::PlannedBuffer<std::uint64_t> points, strengths;
    core::PlannedBuffer<shared_p1::Facet> facets;
    P1AssemblyReport report;
    std::size_t canonical_points = 0;
    core::Precision precision = core::Precision::Binary64;
};

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
                                const P1AssemblyLimits& limits, Work& scan) {
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
                if (i != 0 && strengths[i - 1].value() >= strengths[i].value())
                    throw ArithmeticFailure{invalid};
            }
            for (std::size_t i = 0; i < faces.size(); ++i) {
                scan.step();
                if (i != 0 && faces[i - 1].face.value() >= faces[i].face.value())
                    throw ArithmeticFailure{invalid};
                (void)ordinal(geometry.face_lookup.values(), faces[i].face.value(), scan);
                std::visit(
                    [&](const auto& keys) {
                        for (const auto& key : keys)
                            (void)locate(strengths, key.value(), &P1StrengthId::value, scan);
                    },
                    faces[i].columns);
            }
        };
        finite_outer();
        for (const auto& record : geometry.facets.values()) {
            scan.step();
            const auto face = geometry.faces.values()[static_cast<std::size_t>(record.face)].id;
            const auto row =
                locate(faces, face, [](const P1FaceStrengths& f) { return f.face.value(); }, scan);
            if (std::holds_alternative<QuadOrigin>(record.origin) !=
                std::holds_alternative<std::array<P1StrengthId, 2>>(faces[row].columns))
                throw ArithmeticFailure{invalid};
        }
        const auto n = geometry.facets.size();
        const auto e = product(n, 3);
        const auto p = geometry.points.size();
        const auto s = strengths.size();
        P1AssemblyPlan result;
        result.payload_bytes =
            sum(sum(bytes<SharedP1Data>(1), sum(bytes<std::uint64_t>(p), bytes<std::uint64_t>(s))),
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
        const auto metadata = sum(sum(sum(product(e, 8), sum(product(p, 3), s)), sort),
                                  sum(product(n, sum(30, sum(logarithmic_visits(faces.size()),
                                                             product(3, logarithmic_visits(s))))),
                                      64));
        result.reserved_steps =
            sum(scan.used,
                sum(metadata, product(n, sum(limits.geometry.max_steps, limits.kernel.max_steps))));
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
                                locate(strengths, keys[col].value(), &P1StrengthId::value, work);
                            break;
                        }
                    }
                }
            };
            if (const auto* quad = std::get_if<QuadLRScatter>(&scatter)) {
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
    static AeroResult<P1AssemblyPlan> forecast(const P1ActionSource& source,
                                               std::span<const P1FaceStrengths> faces,
                                               std::span<const P1StrengthId> strengths,
                                               const P1AssemblyLimits& limits,
                                               const std::stop_token& stop) noexcept {
        return guarded<P1AssemblyPlan>(source.query.context(), [&] {
            Work work{0, limits.steps, stop};
            return count(source, faces, strengths, limits, work);
        });
    }

    static AeroResult<SharedP1Action>
    prepare(const P1ActionSource& source, std::span<const P1FaceStrengths> faces,
            std::span<const P1StrengthId> strengths, LocalAccuracy local_accuracy,
            P1Accuracy accuracy, const P1AssemblyLimits& limits, SharedP1Workspace& workspace,
            LocalPacketWorkspace& local, P1Workspace& kernel, core::ResourcePlan& memory,
            const std::stop_token& stop) noexcept {
        const auto context = source.query.context();
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(invalid, context));
        const Lease lease{workspace.entered_};
        return guarded<SharedP1Action>(context, [&] {
            operator_environment();
            Work work{0, limits.steps, stop};
            const auto plan = count(source, faces, strengths, limits, work);
            (void)checked_limit(plan.payload_bytes, memory.remaining_bytes());
            const auto& verified = source_data(source);
            if (!verified.geometry || !verified.geometry->data_)
                throw ArithmeticFailure{invalid};
            const auto& geometry = *verified.geometry->data_;
            const auto facet_count = geometry.facets.size();
            const auto point_slots = geometry.points.size();
            const auto edges = product(facet_count, 3);
            auto storage = allocate<SharedP1Data>(memory, 1, context, stop);
            auto holder = std::make_shared<ActionHolder>(std::move(storage));
            auto& data = holder->payload.values()[0];
            data.points = allocate<std::uint64_t>(memory, point_slots, context, stop);
            data.strengths = allocate<std::uint64_t>(memory, strengths.size(), context, stop);
            data.facets = allocate<Facet>(memory, facet_count, context, stop);
            auto packets =
                allocate<std::optional<LocalGeometryPacket>>(memory, facet_count, context, stop);
            auto point_columns = allocate<std::size_t>(memory, point_slots, context, stop);
            auto occurrences = allocate<EdgeOccurrence>(memory, edges, context, stop);
            auto sorting = allocate<EdgeOccurrence>(memory, edges, context, stop);
            operator_environment();
            data.query = needed(source.query.definition());
            data.stamp = source.stamp;
            data.topology = verified.topology;
            data.geometry.emplace(OwnedGeometry{verified.geometry->data_, context});
            data.precision = accuracy.representation;
            data.report.plan = plan;
            data.report.facets = facet_count;
            data.canonical_points = populate_points(data, geometry, point_columns.values(), work);
            data.material = material(geometry, data.query, point_columns.values(), work);
            for (std::size_t s = 0; s < strengths.size(); ++s) {
                work.step();
                data.strengths.values()[s] = strengths[s].value();
            }
            map_facets(data, geometry, faces, strengths, point_columns.values(),
                       occurrences.values(), work);
            for (std::size_t i = 0; i < facet_count; ++i) {
                auto packet = needed(prepare_local_packet(
                    *verified.geometry, source.query, FacetId{geometry.facets.values()[i].id},
                    local_accuracy, limits.geometry, local, memory, stop));
                charge(work, needed(packet.steps()));
                const bool far = needed(packet.far()) == FarRelation::Far;
                data.facets.values()[i].far = far;
                if (far)
                    ++data.report.far_facets;
                packets.values()[i].emplace(std::move(packet));
            }
            const auto ordered = sort_edges(occurrences.values(), sorting.values(), work);
            data.report.cancelled_analytic_pairs =
                cancel_pairs(ordered, data.facets.values(), work);
            for (std::size_t i = 0; i < facet_count; ++i) {
                auto& entry = data.facets.values()[i];
                const auto& packet = packets.values()[i];
                if (!packet)
                    throw ArithmeticFailure{invalid};
                auto basis = needed(evaluate_p1_partition(*packet, entry.boundary, accuracy,
                                                          limits.kernel, kernel, memory, stop));
                charge(work, needed(basis.steps()));
                entry.basis.emplace(std::move(basis));
            }
            work.step();
            data.report.used_steps = work.used;
            return SharedP1Action{std::shared_ptr<const SharedP1Data>{std::move(holder), &data},
                                  context};
        });
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
                    const auto mu = words(strengths[facet.strength[k]], work);
                    const auto delta_mu = direction
                                              ? words(direction->strengths[facet.strength[k]], work)
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
                    const auto mu = words(strengths[facet.strength[k]], work);
                    auto& mu_gradient = strength_gradient.values()[facet.strength[k]];
                    for (std::size_t a = 0; a < 3; ++a) {
                        const auto value = enclosure(needed(facet.basis->basis(k, a)).get(), work);
                        mu_gradient = add(mu_gradient, multiply(lambda[a], value, work), work);
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
    return P1StrengthId{data_->strengths.values()[index]};
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
