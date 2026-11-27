#include "exact_triangle_private.hpp"
#include "geometry_snapshot_storage.hpp"
#include "local_geometry_input_private.hpp"
#include "query_coordinates_private.hpp"

#include <boost/hash2/sha2.hpp>

#include <algorithm>
#include <bit>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>

namespace uvlas::aero::detail {
namespace {
constexpr auto operation = core::Operation::EvaluateAero;
using Problem = std::optional<GeometryProblem>;
core::Failure complete_details(core::Failure failure, const core::FailureDetails& origin) noexcept {
    auto& d = failure.details;
    if (!d.model)
        d.model = origin.model;
    if (!d.parent)
        d.parent = origin.parent;
    if (!d.trial)
        d.trial = origin.trial;
    if (!d.numerical_generation)
        d.numerical_generation = origin.numerical_generation;
    if (!d.component)
        d.component = origin.component;
    if (!d.problem)
        d.problem = origin.problem;
    if (!d.resources)
        d.resources = origin.resources;
    if (!d.last_accepted)
        d.last_accepted = origin.last_accepted;
    return failure;
}
core::Failure complete(core::Failure failure, const GeometryContext& context) noexcept {
    failure = complete_details(failure, context.diagnostics);
    if (!failure.details.model)
        failure.details.model = context.model;
    if (!failure.details.numerical_generation)
        failure.details.numerical_generation = context.numerical;
    return failure;
}
OwnedGeometryFailure rejected(core::FailureCode code, const GeometryContext& c,
                              Problem p = {}) noexcept {
    return {complete(core::Failure{code, operation}, c), p, c};
}
template <class F>
auto guarded(const GeometryContext& context, F&& action, Problem problem = {}) noexcept
    -> std::invoke_result_t<F> {
    try {
        return std::invoke(std::forward<F>(action));
    } catch (const std::bad_alloc&) {
        return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context, problem));
    } catch (...) {
        return std::unexpected(rejected(core::FailureCode::PrecisionUnresolved, context, problem));
    }
}
bool finite(double x) noexcept {
    return ((std::bit_cast<std::uint64_t>(x) >> 52U) & 2047U) != 2047U;
}
bool finite(ScalarWords words) noexcept {
    if (const auto* p = std::get_if<SingleWord>(&words))
        return finite(p->value);
    if (const auto* p = std::get_if<DoubleWord>(&words))
        return finite(p->high) && finite(p->low);
    return false;
}
void decode_into(double word, unsigned extra, ExactInteger& output) {
    const auto bits = std::bit_cast<std::uint64_t>(word);
    const auto exponent = static_cast<unsigned>((bits >> 52U) & 2047U);
    auto significand = bits & UINT64_C(0xfffffffffffff);
    if (exponent != 0)
        significand |= UINT64_C(0x10000000000000);
    output = significand;
    output <<= (exponent == 0 ? 0U : exponent - 1U) + extra;
    if ((bits >> 63U) != 0)
        output = -output;
}
void scalar_into(ScalarWords words, unsigned extra, ExactInteger& output, ExactInteger& temporary) {
    if (const auto* single = std::get_if<SingleWord>(&words))
        decode_into(single->value, extra, output);
    else if (const auto* pair = std::get_if<DoubleWord>(&words)) {
        decode_into(pair->high, extra, output);
        decode_into(pair->low, extra, temporary);
        output += temporary;
    } else
        throw std::logic_error("Invalid scalar representation");
}
ParameterStorage store(ScalarWords words) {
    if (const auto* p = std::get_if<SingleWord>(&words))
        return p->value;
    if (const auto* p = std::get_if<DoubleWord>(&words))
        return std::array{p->high, p->low};
    throw std::logic_error("Invalid scalar storage");
}
ScalarWords restore(const ParameterStorage& words) {
    if (const auto* p = std::get_if<double>(&words))
        return SingleWord{*p};
    if (const auto* p = std::get_if<std::array<double, 2>>(&words))
        return DoubleWord{(*p)[0], (*p)[1]};
    throw std::logic_error("Invalid parameter storage");
}
// A context is the legacy request stamp. The owned snapshot identity also binds
// the exact time, when present, and is retained by every resolved query.
bool same_context(const GeometryContext& a, const GeometryContext& b) noexcept {
    return a.fingerprint == b.fingerprint && a.model == b.model && a.reference == b.reference &&
           a.current == b.current && a.mapping == b.mapping && a.connectivity == b.connectivity &&
           a.policy == b.policy && a.kernel == b.kernel && a.radius == b.radius &&
           a.precision == b.precision && a.numerical == b.numerical && a.time == b.time &&
           a.diagnostics.parent == b.diagnostics.parent &&
           a.diagnostics.trial == b.diagnostics.trial &&
           a.diagnostics.last_accepted == b.diagnostics.last_accepted;
}
std::optional<core::Index> lookup(std::span<const IdLookup> entries, std::uint64_t id) noexcept {
    const auto found =
        std::lower_bound(entries.begin(), entries.end(), id,
                         [](const IdLookup& row, std::uint64_t value) { return row.id < value; });
    if (found == entries.end() || found->id != id)
        return std::nullopt;
    return found->ordinal;
}
ScalarWords position(const GeometryData& data, core::Index point, Configuration which,
                     std::size_t axis) {
    const auto i = static_cast<std::size_t>(point) * 3 + axis;
    const auto& high = which == Configuration::Reference ? data.reference_high : data.current_high;
    const auto& low = which == Configuration::Reference ? data.reference_low : data.current_low;
    if (low.size() == 0)
        return SingleWord{high.values()[i]};
    return DoubleWord{high.values()[i], low.values()[i]};
}
struct Run {
    const GeometryLimits& limits;
    GeometryContext context;
    std::stop_token stop;
    std::uint64_t steps = 0;
    Run(const GeometryLimits& l, GeometryContext c, std::stop_token s) noexcept
        : limits(l), context(c), stop(std::move(s)) {}
    // This checks a whole phase without charging unvisited work.
    [[nodiscard]] AeroResult<void> check_steps(std::uint64_t count) const noexcept {
        if (stop.stop_requested())
            return std::unexpected(rejected(core::FailureCode::Cancelled, context));
        if (steps > limits.max_steps || count > limits.max_steps - steps)
            return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
        return {};
    }
    AeroResult<void> charge(std::uint64_t count) noexcept {
        if (const auto allowed = check_steps(count); !allowed)
            return allowed;
        steps += count;
        return {};
    }
    AeroResult<void> step(Problem problem = {}) noexcept {
        if (stop.stop_requested())
            return std::unexpected(rejected(core::FailureCode::Cancelled, context, problem));
        if (steps == limits.max_steps)
            return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context, problem));
        ++steps;
        return {};
    }
    AeroResult<void> admission() noexcept {
        if (stop.stop_requested())
            return std::unexpected(rejected(core::FailureCode::Cancelled, context));
        const auto memory = limits.predicate.external_memory.selection();
        if (std::holds_alternative<numerics::linear::RequireProvenBound>(memory))
            return std::unexpected(rejected(core::FailureCode::Unsupported, context));
        const auto* e = std::get_if<numerics::linear::EstimateAdmission>(&memory);
        if (e == nullptr)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
        if (e->bytes < predicate_memory_plan().peak_payload_estimate_bytes())
            return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
        return {};
    }
};
bool fits(std::uint64_t actual, std::uint64_t maximum) noexcept {
    return actual <= maximum &&
           actual <= static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
}
std::uint8_t word_count(PositionView view) noexcept {
    return std::holds_alternative<SinglePositions>(view) ? 1 : 2;
}
bool shape(PositionView view, std::uint64_t count) noexcept {
    if (count > std::numeric_limits<std::size_t>::max() / 3)
        return false;
    const auto size = static_cast<std::size_t>(count) * 3;
    if (const auto* p = std::get_if<SinglePositions>(&view))
        return p->words.size() == size;
    if (const auto* p = std::get_if<DoublePositions>(&view))
        return p->high.size() == size && p->low.size() == size;
    return false;
}
ScalarWords input_position(PositionView view, std::size_t index) {
    if (const auto* p = std::get_if<SinglePositions>(&view))
        return SingleWord{p->words[index]};
    if (const auto* p = std::get_if<DoublePositions>(&view))
        return DoubleWord{p->high[index], p->low[index]};
    throw std::logic_error("Invalid position view");
}
AeroResult<std::uint64_t> footprint(const GeometryCounts& n, std::uint8_t wr, std::uint8_t wc,
                                    std::uint64_t time_words,
                                    const GeometryContext& context) noexcept {
    std::uint64_t total = 0;
    const auto add = [&](std::uint64_t count, std::uint64_t bytes) {
        if (bytes != 0 && count > (std::numeric_limits<std::uint64_t>::max() - total) / bytes)
            return false;
        total += count * bytes;
        return true;
    };
    if (n.quads > (std::numeric_limits<std::uint64_t>::max() - n.triangles) / 2)
        return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
    const auto facets = n.triangles + 2 * n.quads;
    if (!fits(facets, std::numeric_limits<std::uint64_t>::max()) ||
        !add(n.points, 24ULL * (wr + wc)) ||
        !add(n.points, sizeof(PointStorage) + sizeof(IdLookup)) ||
        !add(n.edges, sizeof(EdgeStorage) + sizeof(IdLookup)) ||
        !add(n.faces, sizeof(FaceStorage) + sizeof(IdLookup)) ||
        !add(n.bands, sizeof(std::uint64_t) + sizeof(IdLookup)) ||
        !add(n.triangles, sizeof(TriangleStorage)) || !add(n.quads, sizeof(QuadStorage)) ||
        !add(facets, sizeof(FacetStorage) + sizeof(IdLookup)) ||
        !add(time_words, sizeof(std::uint32_t)))
        return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
    // SHA-256's canonical message has a 64-bit bit count, not an arbitrary wraparound.
    if (total >
        std::numeric_limits<std::uint64_t>::max() / 8 - 1024 - time_words * sizeof(std::uint32_t))
        return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
    return total;
}
// Two finite binary64 encodings summed on the 2^-1074 grid require at
// most 2099 bits. Leave that headroom for the exact rounding comparisons.
constexpr std::size_t time_integer_words = (std::numeric_limits<ExactInteger>::digits - 2099) / 32;

AeroResult<void> validate_time(pchip::RationalView time, Run& run) {
    const auto invalid = [&] {
        return std::unexpected(rejected(core::FailureCode::InvalidInput, run.context));
    };
    if (time.denominator.empty() || time.denominator.back() == 0)
        return invalid();
    const auto raw = std::bit_cast<std::uint64_t>(run.context.time);
    const auto bits = raw & UINT64_C(0x7fffffffffffffff);
    if (time.sign == 0) {
        if (!time.numerator.empty() || time.denominator.size() != 1 || time.denominator[0] != 1 ||
            bits != 0)
            return invalid();
        return {};
    }
    if (time.sign != ((raw >> 63U) == 0 ? 1 : -1) || time.numerator.empty() ||
        time.numerator.back() == 0)
        return invalid();

    const auto decode = [&](std::span<const std::uint32_t> words,
                            ExactInteger& value) -> AeroResult<void> {
        std::uint64_t occupied = 0;
        for (auto it = words.rbegin(); it != words.rend(); ++it) {
            if (const auto step = run.charge(++occupied); !step)
                return step;
            value <<= 32;
            value += *it;
        }
        return {};
    };
    ExactInteger numerator, denominator;
    if (const auto decoded = decode(time.numerator, numerator); !decoded)
        return decoded;
    if (const auto decoded = decode(time.denominator, denominator); !decoded)
        return decoded;
    const auto integer_words = [](const ExactInteger& value) {
        return static_cast<std::uint64_t>(boost::multiprecision::msb(value)) / 32 + 1;
    };
    auto a = numerator, b = denominator;
    while (b != 0) {
        const auto left = integer_words(a), right = integer_words(b);
        // Charge a conservative schoolbook quotient bound in 32-bit limb visits.
        if (const auto step = run.charge(left * right + left + right); !step)
            return step;
        ExactInteger remainder = a % b;
        a = std::move(b);
        b = std::move(remainder);
    }
    if (a != 1)
        return invalid();

    // The binary64 grid midpoints fit in 66 limbs. Charge two full products,
    // the grid decodes/additions, and shifting/comparing the exact numerator.
    constexpr std::uint64_t binary_grid_words = (2099 + 31) / 32;
    const auto widths = time.numerator.size() + time.denominator.size();
    const auto midpoint_work = 2 * binary_grid_words * (time.denominator.size() + 3) + 3 * widths;
    if (const auto step = run.charge(midpoint_work); !step)
        return step;
    // Compare exact adjacent midpoints, including subnormals and the overflow
    // midpoint, without taking a floating-point midpoint. Ties go to even.
    numerator <<= 1075;
    ExactInteger center, neighbor;
    decode_into(std::bit_cast<double>(bits), 0, center);
    const bool odd = (bits & 1U) != 0;
    if (bits != 0) {
        decode_into(std::bit_cast<double>(bits - 1), 0, neighbor);
        neighbor = (center + neighbor) * denominator;
        if (numerator < neighbor || (numerator == neighbor && odd))
            return invalid();
    }
    if (bits == UINT64_C(0x7fefffffffffffff))
        neighbor = ExactInteger{1} << 2098;
    else
        decode_into(std::bit_cast<double>(bits + 1), 0, neighbor);
    neighbor = (center + neighbor) * denominator;
    if (numerator > neighbor || (numerator == neighbor && odd))
        return invalid();
    return {};
}

template <class T>
AeroResult<void> allocate(core::PlannedBuffer<T>& destination, std::uint64_t count,
                          core::ResourcePlan& plan, const GeometryContext& context) noexcept {
    if (count > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
    auto result = plan.allocate<T>(static_cast<core::Index>(count), operation);
    if (!result)
        return std::unexpected(
            OwnedGeometryFailure{complete(result.error(), context), {}, context});
    destination = std::move(*result);
    return {};
}
template <class F, class P>
AeroResult<void> build_lookup(core::PlannedBuffer<IdLookup>& output, F&& id, Run& run,
                              P&& problem) {
    for (std::size_t i = 0; i < output.size(); ++i)
        output.values()[i] = {id(i), static_cast<core::Index>(i)};
    std::sort(output.values().begin(), output.values().end(),
              [](const IdLookup& a, const IdLookup& b) {
                  return a.id < b.id || (a.id == b.id && a.ordinal < b.ordinal);
              });
    for (std::size_t i = 0; i < output.size(); ++i) {
        if (auto step = run.step(problem(output.values()[i].id)); !step)
            return step;
        if (i > 0 && output.values()[i - 1].id == output.values()[i].id)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, run.context,
                                            problem(output.values()[i].id)));
    }
    return {};
}
} // namespace

enum class GeometryWorkAccounting : std::uint8_t {
    GeometryAndPredicateBudgets,
    TotalAndPredicateSubcap
};
class GeometryEngine {
public:
    static OwnedGeometry retain(const LocalGeometryInput& input) noexcept {
        return OwnedGeometry{input.owner, input.owner->context};
    }
    static AeroResult<LocalGeometryInput> local(const OwnedGeometry&, const ResolvedQuery&, FacetId,
                                                const GeometryLimits&, GeometryWorkspace&,
                                                const PredicateLease&, core::ResourcePlan&,
                                                std::stop_token, GeometryWorkAccounting,
                                                TriangleOutcome* = nullptr) noexcept;
    static AeroResult<OwnedGeometry> prepare(const GeometryInput&, const GeometryLimits&,
                                             GeometryWorkspace&, core::ResourcePlan&,
                                             std::stop_token) noexcept;
    static AeroResult<ResolvedQuery> resolve(const OwnedGeometry&, const QueryRequest&,
                                             const GeometryContext&, const GeometryLimits&,
                                             GeometryWorkspace&, core::ResourcePlan&,
                                             std::stop_token) noexcept;
    static AeroResult<ExactTriangleResult> classify(const OwnedGeometry&, const ResolvedQuery&,
                                                    FacetId, const GeometryLimits&,
                                                    GeometryWorkspace&, core::ResourcePlan&,
                                                    std::stop_token) noexcept;

    static AeroResult<ExactQueryCoordinates> coordinates(const OwnedGeometry&, const ResolvedQuery&,
                                                         const GeometryLimits&, GeometryWorkspace&,
                                                         std::stop_token, Configuration) noexcept;

private:
    static AeroResult<void> copy(const GeometryInput&, GeometryData&, Run&);
    static AeroResult<void> lookups(GeometryData&, Run&);
    static AeroResult<void> connect(GeometryData&, GeometryWorkspace&, const PredicateLease&, Run&);
    static AeroResult<void> check_facets(GeometryData&, GeometryWorkspace&, const PredicateLease&,
                                         Run&);
    static AeroResult<void> fill_facet(GeometryData&, std::size_t, std::uint64_t, std::uint64_t,
                                       const std::array<std::uint64_t, 3>&,
                                       const std::array<std::uint64_t, 3>&,
                                       std::variant<IndependentOrigin, QuadOrigin>, Run&);
    static AeroResult<void> load_point(const GeometryData&, core::Index, Configuration,
                                       std::array<ExactInteger, 3>&, QueryArithmeticScratch&, Run&);
    static AeroResult<void> query_position(const GeometryData&, const QueryStorage&, ExactPoints&,
                                           QueryArithmeticScratch&, Run&,
                                           ExactQueryCoordinates* = nullptr,
                                           Configuration = Configuration::Current);
    static void hash(GeometryData&);
};

AeroResult<OwnedGeometry> GeometryEngine::prepare(const GeometryInput& input,
                                                  const GeometryLimits& limits,
                                                  GeometryWorkspace& workspace,
                                                  core::ResourcePlan& plan,
                                                  std::stop_token stop) noexcept {
    const GeometryContext context = input.context;
    auto lease = workspace.acquire(context.diagnostics);
    if (!lease)
        return std::unexpected(OwnedGeometryFailure{complete(lease.error(), context), {}, context});
    return guarded(context, [&]() -> AeroResult<OwnedGeometry> {
        const GeometryCounts n{input.points.size(), input.edges.size(),     input.faces.size(),
                               input.bands.size(),  input.triangles.size(), input.quads.size()};
        Run run{limits, context, stop};
        if (input.bindings.size() != n.points || !shape(input.reference, n.points) ||
            !shape(input.current, n.points))
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
        if (!fits(n.points, limits.counts.points) || !fits(n.edges, limits.counts.edges) ||
            !fits(n.faces, limits.counts.faces) || !fits(n.bands, limits.counts.bands) ||
            !fits(n.triangles, limits.counts.triangles) || !fits(n.quads, limits.counts.quads))
            return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
        const auto instant = input.exact_time;
        if (instant && (instant->numerator.size() > time_integer_words ||
                        instant->denominator.size() > time_integer_words))
            return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
        const auto time_words =
            instant ? instant->numerator.size() + instant->denominator.size() : 0;
        auto bytes = footprint(n, word_count(input.reference), word_count(input.current),
                               time_words, context);
        if (!bytes)
            return std::unexpected(bytes.error());
        // Admit the complete point-word phase before selecting its canonical
        // failure. No lookup/payload is needed and no future step is charged.
        if (stop.stop_requested())
            return std::unexpected(rejected(core::FailureCode::Cancelled, context));
        if (n.points > std::numeric_limits<std::uint64_t>::max() / 6)
            return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
        if (auto admitted = run.check_steps(6 * n.points); !admitted)
            return std::unexpected(admitted.error());
        std::optional<IdLookup> nonfinite;
        for (const auto view : {input.reference, input.current})
            for (std::size_t i = 0; i < input.points.size() * 3; ++i) {
                if (auto step = run.step(); !step)
                    return std::unexpected(step.error());
                if (!finite(input_position(view, i))) {
                    const IdLookup candidate{input.points[i / 3].value(),
                                             static_cast<core::Index>(i / 3)};
                    if (!nonfinite || candidate.id < nonfinite->id ||
                        (candidate.id == nonfinite->id && candidate.ordinal < nonfinite->ordinal))
                        nonfinite = candidate;
                }
            }
        if (nonfinite)
            return std::unexpected(
                rejected(core::FailureCode::Nonfinite, context, PointId{nonfinite->id}));
        if (!finite(context.time) || !finite(input.core_radius))
            return std::unexpected(rejected(core::FailureCode::Nonfinite, context));
        if (auto admitted = run.admission(); !admitted)
            return std::unexpected(admitted.error());
        if (*bytes > plan.remaining_bytes()) {
            auto failure = rejected(core::FailureCode::ResourceExceeded, context);
            const auto committed = plan.byte_limit() - plan.remaining_bytes();
            if (*bytes <= std::numeric_limits<std::uint64_t>::max() - committed)
                failure.cause.details.resources =
                    core::ResourceUsage{core::ResourceKind::Bytes, committed + *bytes, committed};
            return std::unexpected(failure);
        }
        auto data = std::make_shared<GeometryData>(
            context, n, input.core_radius, word_count(input.reference), word_count(input.current));
        data->payload = *bytes;
        const auto facets = n.triangles + 2 * n.quads;
        const auto position_count = n.points * 3;
        // Keep every allocation on the same cumulative preparation plan.
        const auto allocate_all = [&]() -> AeroResult<void> {
            if (auto r = allocate(data->reference_high, position_count, plan, context); !r)
                return r;
            if (auto r = allocate(data->reference_low,
                                  data->reference_words == 2 ? position_count : 0, plan, context);
                !r)
                return r;
            if (auto r = allocate(data->current_high, position_count, plan, context); !r)
                return r;
            if (auto r = allocate(data->current_low, data->current_words == 2 ? position_count : 0,
                                  plan, context);
                !r)
                return r;
            if (auto r = allocate(data->points, n.points, plan, context); !r)
                return r;
            if (auto r = allocate(data->edges, n.edges, plan, context); !r)
                return r;
            if (auto r = allocate(data->faces, n.faces, plan, context); !r)
                return r;
            if (auto r = allocate(data->bands, n.bands, plan, context); !r)
                return r;
            if (auto r = allocate(data->triangles, n.triangles, plan, context); !r)
                return r;
            if (auto r = allocate(data->quads, n.quads, plan, context); !r)
                return r;
            if (auto r = allocate(data->facets, facets, plan, context); !r)
                return r;
            if (auto r = allocate(data->point_lookup, n.points, plan, context); !r)
                return r;
            if (auto r = allocate(data->edge_lookup, n.edges, plan, context); !r)
                return r;
            if (auto r = allocate(data->face_lookup, n.faces, plan, context); !r)
                return r;
            if (auto r = allocate(data->band_lookup, n.bands, plan, context); !r)
                return r;
            if (auto r = allocate(data->facet_lookup, facets, plan, context); !r)
                return r;
            return allocate(data->time_words, time_words, plan, context);
        };
        if (auto r = allocate_all(); !r)
            return std::unexpected(r.error());
        if (auto r = copy(input, *data, run); !r)
            return std::unexpected(r.error());
        if (instant) {
            data->time_sign = instant->sign;
            data->time_numerator_words = instant->numerator.size();
            auto words = data->time_words.values();
            std::size_t i = 0;
            for (const auto part : {instant->numerator, instant->denominator})
                for (const auto word : part) {
                    if (auto step = run.step(); !step)
                        return std::unexpected(step.error());
                    words[i++] = word;
                }
            const auto split = data->time_numerator_words;
            const pchip::RationalView captured{data->time_sign, words.first(split),
                                               words.subspan(split)};
            // A present time requires a denominator, even for exact zero.
            if (auto checked = validate_time(captured, run); !checked)
                return std::unexpected(checked.error());
        }
        // No external input is read after this snapshot boundary.
        if (auto r = lookups(*data, run); !r)
            return std::unexpected(r.error());
        if (auto r = connect(*data, workspace, *lease, run); !r)
            return std::unexpected(r.error());
        auto& arithmetic = workspace.arithmetic_.values()[0];
        scalar_into(data->radius, 0, arithmetic.coordinate, arithmetic.word);
        if (arithmetic.coordinate < 0)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
        if (auto r = check_facets(*data, workspace, *lease, run); !r)
            return std::unexpected(r.error());
        if (auto step = run.step(); !step)
            return std::unexpected(step.error());
        if (instant) {
            if (const auto hashed = run.charge(4 + time_words); !hashed)
                return std::unexpected(hashed.error());
        }
        hash(*data);
        data->preparation_steps = run.steps;
        return OwnedGeometry{std::move(data), context};
    });
}

AeroResult<void> GeometryEngine::copy(const GeometryInput& input, GeometryData& data, Run& run) {
    const auto copy_positions = [&](PositionView view, core::PlannedBuffer<double>& high,
                                    core::PlannedBuffer<double>& low) {
        for (std::size_t i = 0; i < high.size(); ++i) {
            const auto value = input_position(view, i);
            if (const auto* single = std::get_if<SingleWord>(&value))
                high.values()[i] = single->value;
            else if (const auto* pair = std::get_if<DoubleWord>(&value)) {
                high.values()[i] = pair->high;
                low.values()[i] = pair->low;
            }
        }
    };
    copy_positions(input.reference, data.reference_high, data.reference_low);
    copy_positions(input.current, data.current_high, data.current_low);
    for (std::size_t i = 0; i < data.points.size(); ++i) {
        if (auto step = run.step(input.points[i]); !step)
            return step;
        auto& row = data.points.values()[i];
        const auto& b = input.bindings[i];
        row.id = input.points[i].value();
        row.map = b.map.value();
        row.map_generation = b.generation.value();
        if (std::holds_alternative<CanonicalPoint>(b.ownership))
            row.original = CanonicalPoint{};
        else if (const auto* alias = std::get_if<AliasOf>(&b.ownership))
            row.original = AliasRequestStorage{alias->owner.value()};
        else
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, input.points[i]));
    }
    for (std::size_t i = 0; i < data.edges.size(); ++i) {
        if (auto step = run.step(input.edges[i].id); !step)
            return step;
        auto& row = data.edges.values()[i];
        row.id = input.edges[i].id.value();
        for (std::size_t j = 0; j < 2; ++j)
            row.original[j] = input.edges[i].points[j].value();
    }
    for (std::size_t i = 0; i < data.faces.size(); ++i) {
        if (auto step = run.step(input.faces[i].id); !step)
            return step;
        data.faces.values()[i].id = input.faces[i].id.value();
        data.faces.values()[i].band = input.faces[i].band.value();
    }
    for (std::size_t i = 0; i < data.bands.size(); ++i)
        data.bands.values()[i] = input.bands[i].value();
    for (std::size_t i = 0; i < data.triangles.size(); ++i) {
        if (auto step = run.step(input.triangles[i].id); !step)
            return step;
        auto& row = data.triangles.values()[i];
        const auto& source = input.triangles[i];
        row.id = source.id.value();
        row.parent = source.face.value();
        for (std::size_t j = 0; j < 3; ++j) {
            row.points[j] = source.points[j].value();
            row.edges[j] = source.edges[j].value();
        }
    }
    for (std::size_t i = 0; i < data.quads.size(); ++i) {
        if (auto step = run.step(input.quads[i].facets[0]); !step)
            return step;
        auto& row = data.quads.values()[i];
        const auto& source = input.quads[i];
        row.parent = source.face.value();
        row.diagonal = source.diagonal.value();
        for (std::size_t j = 0; j < 4; ++j) {
            row.points[j] = source.points[j].value();
            row.edges[j] = source.edges[j].value();
        }
        for (std::size_t j = 0; j < 2; ++j)
            row.facets[j] = source.facets[j].value();
    }
    return {};
}

AeroResult<void> GeometryEngine::lookups(GeometryData& data, Run& run) {
    for (std::size_t i = 0; i < data.triangles.size(); ++i)
        data.facets.values()[i].id = data.triangles.values()[i].id;
    for (std::size_t i = 0; i < data.quads.size(); ++i)
        for (std::size_t half = 0; half < 2; ++half)
            data.facets.values()[data.triangles.size() + 2 * i + half].id =
                data.quads.values()[i].facets[half];
    if (auto r = build_lookup(
            data.point_lookup, [&](std::size_t i) { return data.points.values()[i].id; }, run,
            [](std::uint64_t id) -> GeometryProblem { return PointId{id}; });
        !r)
        return r;
    if (auto r = build_lookup(
            data.edge_lookup, [&](std::size_t i) { return data.edges.values()[i].id; }, run,
            [](std::uint64_t id) -> GeometryProblem { return EdgeId{id}; });
        !r)
        return r;
    if (auto r = build_lookup(
            data.face_lookup, [&](std::size_t i) { return data.faces.values()[i].id; }, run,
            [](std::uint64_t id) -> GeometryProblem { return PhysicalFaceId{id}; });
        !r)
        return r;
    if (auto r = build_lookup(
            data.band_lookup, [&](std::size_t i) { return data.bands.values()[i]; }, run,
            [](std::uint64_t id) -> GeometryProblem { return BandId{id}; });
        !r)
        return r;
    return build_lookup(
        data.facet_lookup, [&](std::size_t i) { return data.facets.values()[i].id; }, run,
        [](std::uint64_t id) -> GeometryProblem { return FacetId{id}; });
}

AeroResult<void> GeometryEngine::connect(GeometryData& data, GeometryWorkspace& workspace,
                                         const PredicateLease&, Run& run) {
    auto& arithmetic = workspace.arithmetic_.values()[0];
    for (const auto entry : data.point_lookup.values()) {
        auto& point = data.points.values()[static_cast<std::size_t>(entry.ordinal)];
        if (auto step = run.step(PointId{point.id}); !step)
            return step;
        if (std::holds_alternative<CanonicalPoint>(point.original))
            point.owner = entry.ordinal;
        else if (const auto* alias = std::get_if<AliasRequestStorage>(&point.original)) {
            const auto owner = lookup(data.point_lookup.values(), alias->owner);
            if (!owner || *owner == entry.ordinal ||
                !std::holds_alternative<CanonicalPoint>(
                    data.points.values()[static_cast<std::size_t>(*owner)].original))
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, run.context, PointId{point.id}));
            point.owner = *owner;
            const auto& canonical = data.points.values()[static_cast<std::size_t>(*owner)];
            if (point.map != canonical.map)
                return std::unexpected(
                    rejected(core::FailureCode::GeometryInvalid, run.context, PointId{point.id}));
            if (point.map_generation != canonical.map_generation)
                return std::unexpected(
                    rejected(core::FailureCode::StaleGeneration, run.context, PointId{point.id}));
            for (const auto configuration : {Configuration::Reference, Configuration::Current})
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    if (auto step = run.step(PointId{point.id}); !step)
                        return step;
                    scalar_into(position(data, entry.ordinal, configuration, axis), 0,
                                arithmetic.coordinate, arithmetic.word);
                    scalar_into(position(data, *owner, configuration, axis), 0, arithmetic.sum,
                                arithmetic.word);
                    if (arithmetic.coordinate != arithmetic.sum)
                        return std::unexpected(rejected(core::FailureCode::GeometryInvalid,
                                                        run.context, PointId{point.id}));
                }
        } else
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, PointId{point.id}));
    }
    for (const auto entry : data.edge_lookup.values()) {
        auto& edge = data.edges.values()[static_cast<std::size_t>(entry.ordinal)];
        if (auto step = run.step(EdgeId{edge.id}); !step)
            return step;
        for (std::size_t i = 0; i < 2; ++i) {
            const auto p = lookup(data.point_lookup.values(), edge.original[i]);
            if (!p)
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, run.context, EdgeId{edge.id}));
            edge.points[i] = data.points.values()[static_cast<std::size_t>(*p)].owner;
        }
        if (edge.points[0] == edge.points[1])
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, EdgeId{edge.id}));
        if (data.points.values()[static_cast<std::size_t>(edge.points[0])].id >
            data.points.values()[static_cast<std::size_t>(edge.points[1])].id)
            std::swap(edge.points[0], edge.points[1]);
    }
    for (const auto entry : data.face_lookup.values()) {
        const auto& face = data.faces.values()[static_cast<std::size_t>(entry.ordinal)];
        if (auto step = run.step(PhysicalFaceId{face.id}); !step)
            return step;
        if (!lookup(data.band_lookup.values(), face.band))
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, PhysicalFaceId{face.id}));
    }
    // Inspect sorted IDs while writing each source's original expanded ordinal.
    // Independent triangles precede quad children; half order follows child IDs.
    for (const auto entry : data.facet_lookup.values()) {
        const auto ordinal = static_cast<std::size_t>(entry.ordinal);
        if (ordinal >= data.triangles.size())
            continue;
        const auto& t = data.triangles.values()[ordinal];
        if (auto result = fill_facet(data, ordinal, t.id, t.parent, t.points, t.edges,
                                     IndependentOrigin{entry.ordinal}, run);
            !result)
            return result;
    }
    for (const auto entry : data.facet_lookup.values()) {
        const auto ordinal = static_cast<std::size_t>(entry.ordinal);
        if (ordinal < data.triangles.size())
            continue;
        const auto expanded = ordinal - data.triangles.size();
        const auto i = expanded / 2, half = expanded % 2;
        const auto& q = data.quads.values()[i];
        // Fixed positive geometry is independent of validation order.
        if (half == 0) {
            if (auto r = fill_facet(data, ordinal, q.facets[0], q.parent,
                                    {q.points[0], q.points[2], q.points[1]},
                                    {q.diagonal, q.edges[1], q.edges[0]},
                                    QuadOrigin{static_cast<core::Index>(i), 0}, run);
                !r)
                return r;
        } else {
            if (auto r = fill_facet(data, ordinal, q.facets[1], q.parent,
                                    {q.points[0], q.points[3], q.points[2]},
                                    {q.edges[3], q.edges[2], q.diagonal},
                                    QuadOrigin{static_cast<core::Index>(i), 1}, run);
                !r)
                return r;
        }
    }
    return {};
}

AeroResult<void> GeometryEngine::fill_facet(GeometryData& data, std::size_t ordinal,
                                            std::uint64_t id, std::uint64_t face,
                                            const std::array<std::uint64_t, 3>& points,
                                            const std::array<std::uint64_t, 3>& edges,
                                            std::variant<IndependentOrigin, QuadOrigin> origin,
                                            Run& run) {
    if (auto step = run.step(FacetId{id}); !step)
        return step;
    auto& row = data.facets.values()[ordinal];
    row.id = id;
    row.origin = origin;
    const auto parent = lookup(data.face_lookup.values(), face);
    if (!parent)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, run.context, FacetId{id}));
    row.face = *parent;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto p = lookup(data.point_lookup.values(), points[i]);
        const auto e = lookup(data.edge_lookup.values(), edges[i]);
        if (!p || !e)
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, FacetId{id}));
        row.points[i] = data.points.values()[static_cast<std::size_t>(*p)].owner;
        row.edges[i] = *e;
    }
    if (row.points[0] == row.points[1] || row.points[1] == row.points[2] ||
        row.points[2] == row.points[0])
        return std::unexpected(
            rejected(core::FailureCode::GeometryInvalid, run.context, FacetId{id}));
    std::array<int, 3> boundary{};
    for (std::size_t i = 0; i < 3; ++i) {
        const auto& edge = data.edges.values()[static_cast<std::size_t>(row.edges[i])];
        const auto a = row.points[i], b = row.points[(i + 1) % 3];
        if (edge.points[0] == a && edge.points[1] == b)
            row.incidence[i] = 1;
        else if (edge.points[0] == b && edge.points[1] == a)
            row.incidence[i] = -1;
        else
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, FacetId{id}));
        for (std::size_t vertex = 0; vertex < 3; ++vertex) {
            if (row.points[vertex] == edge.points[0])
                boundary[vertex] -= row.incidence[i];
            if (row.points[vertex] == edge.points[1])
                boundary[vertex] += row.incidence[i];
        }
    }
    if (std::any_of(boundary.begin(), boundary.end(),
                    [](int coefficient) { return coefficient != 0; }))
        return std::unexpected(rejected(core::FailureCode::InvalidInput, run.context, FacetId{id}));
    return {};
}

AeroResult<void> GeometryEngine::load_point(const GeometryData& data, core::Index point,
                                            Configuration configuration,
                                            std::array<ExactInteger, 3>& destination,
                                            QueryArithmeticScratch& scratch, Run& run) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (auto step = run.step(PointId{data.points.values()[static_cast<std::size_t>(point)].id});
            !step)
            return step;
        scalar_into(position(data, point, configuration, axis), 1074, destination[axis],
                    scratch.word);
    }
    return {};
}
AeroResult<void> GeometryEngine::check_facets(GeometryData& data, GeometryWorkspace& workspace,
                                              const PredicateLease& lease, Run& run) {
    auto points =
        ExactPredicateAccess::coordinates(*workspace.predicate_, lease, run.context.diagnostics);
    if (!points)
        return std::unexpected(
            OwnedGeometryFailure{complete(points.error(), run.context), {}, run.context});
    for (const auto entry : data.facet_lookup.values()) {
        const auto& facet = data.facets.values()[static_cast<std::size_t>(entry.ordinal)];
        for (const auto configuration : {Configuration::Reference, Configuration::Current}) {
            for (std::size_t vertex = 0; vertex < 3; ++vertex)
                if (auto r =
                        load_point(data, facet.points[vertex], configuration, points->get()[vertex],
                                   workspace.arithmetic_.values()[0], run);
                    !r)
                    return r;
            if (auto step = run.step(FacetId{facet.id}); !step)
                return step;
            const auto valid = ExactPredicateAccess::nondegenerate(*workspace.predicate_, lease,
                                                                   run.context.diagnostics);
            if (!valid)
                return std::unexpected(OwnedGeometryFailure{complete(valid.error(), run.context),
                                                            FacetId{facet.id}, run.context});
            if (!*valid)
                return std::unexpected(
                    rejected(core::FailureCode::GeometryInvalid, run.context, FacetId{facet.id}));
        }
    }
    return {};
}

void GeometryEngine::hash(GeometryData& data) {
    boost::hash2::sha2_256 digest;
    constexpr std::array<unsigned char, 8> domain{'U', 'V', 'L', 'A', 'S', 'G', '0', '1'};
    digest.update(domain.data(), domain.size());
    const auto word = [&](std::uint64_t value) {
        std::array<unsigned char, 8> bytes{};
        for (auto& byte : bytes) {
            byte = static_cast<unsigned char>(value & 255U);
            value >>= 8U;
        }
        digest.update(bytes.data(), bytes.size());
    };
    const auto scalar = [&](ScalarWords value) {
        if (const auto* single = std::get_if<SingleWord>(&value)) {
            word(1);
            word(std::bit_cast<std::uint64_t>(single->value));
        } else if (const auto* pair = std::get_if<DoubleWord>(&value)) {
            word(2);
            word(std::bit_cast<std::uint64_t>(pair->high));
            word(std::bit_cast<std::uint64_t>(pair->low));
        } else
            throw std::logic_error("Invalid scalar hash input");
    };
    const auto optional_id = [&](const auto& id) {
        word(id.has_value() ? 1 : 0);
        if (id)
            word(id->value());
    };
    const auto& c = data.context;
    digest.update(c.fingerprint.bytes.data(), c.fingerprint.bytes.size());
    word(c.model.value());
    word(c.reference.value());
    word(c.current.value());
    word(c.mapping.value());
    word(c.connectivity.value());
    word(c.policy.value());
    word(c.kernel.value());
    word(c.radius.value());
    word(c.precision.value());
    word(c.numerical.value());
    word(std::bit_cast<std::uint64_t>(c.time));
    optional_id(c.diagnostics.parent);
    optional_id(c.diagnostics.trial);
    optional_id(c.diagnostics.last_accepted);
    // Resource caps and diagnostic component/problem/resources are not geometric content.
    for (const auto n : {data.counts.points, data.counts.edges, data.counts.faces,
                         data.counts.bands, data.counts.triangles, data.counts.quads})
        word(n);
    scalar(data.radius);
    word(data.reference_words);
    word(data.current_words);
    for (const auto* plane :
         {&data.reference_high, &data.reference_low, &data.current_high, &data.current_low})
        for (const double value : plane->values())
            word(std::bit_cast<std::uint64_t>(value));
    for (const auto& p : data.points.values()) {
        word(p.id);
        word(p.map);
        word(p.map_generation);
        if (const auto* a = std::get_if<AliasRequestStorage>(&p.original)) {
            word(1);
            word(a->owner);
        } else
            word(0);
    }
    for (const auto& e : data.edges.values()) {
        word(e.id);
        for (const auto id : e.original)
            word(id);
    }
    for (const auto& f : data.faces.values()) {
        word(f.id);
        word(f.band);
    }
    for (const auto id : data.bands.values())
        word(id);
    for (const auto& t : data.triangles.values()) {
        word(t.id);
        word(t.parent);
        for (const auto id : t.points)
            word(id);
        for (const auto id : t.edges)
            word(id);
    }
    for (const auto& q : data.quads.values()) {
        word(q.parent);
        for (const auto id : q.points)
            word(id);
        for (const auto id : q.edges)
            word(id);
        word(q.diagonal);
        for (const auto id : q.facets)
            word(id);
    }
    // Preserve the legacy identity when no explicit exact instant was supplied.
    if (data.time_words.size() != 0) {
        word(UINT64_C(0x55564c415354494d));
        const auto encoded_sign = data.time_sign + 1;
        word(static_cast<std::uint64_t>(encoded_sign));
        word(data.time_numerator_words);
        word(data.time_words.size());
        for (const auto digit : data.time_words.values())
            word(digit);
    }
    const auto result = digest.result();
    for (std::size_t i = 0; i < data.identity.size(); ++i)
        data.identity[i] = static_cast<std::byte>(result[i]);
}

AeroResult<std::optional<pchip::RationalView>> OwnedGeometry::exact_time() const& noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
    if (data_->time_words.size() == 0)
        return std::optional<pchip::RationalView>{};
    const auto words = data_->time_words.values();
    const auto split = data_->time_numerator_words;
    return std::optional{
        pchip::RationalView{data_->time_sign, words.first(split), words.subspan(split)}};
}

AeroResult<std::array<std::byte, 32>> OwnedGeometry::identity() const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
    return data_->identity;
}
AeroResult<GeometryCounts> OwnedGeometry::counts() const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
    return data_->counts;
}
AeroResult<std::uint64_t> OwnedGeometry::payload_bytes() const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
    return data_->payload;
}
AeroResult<std::uint64_t> OwnedGeometry::used_preparation_steps() const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
    return data_->preparation_steps;
}
AeroResult<PointId> OwnedGeometry::owner(PointId id) const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto index = lookup(data_->point_lookup.values(), id.value());
    if (!index)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto owner = data_->points.values()[static_cast<std::size_t>(*index)].owner;
    return PointId{data_->points.values()[static_cast<std::size_t>(owner)].id};
}
AeroResult<std::array<ScalarWords, 3>>
OwnedGeometry::point_words(PointId id, Configuration which) const noexcept {
    return guarded(context_, [&]() -> AeroResult<std::array<ScalarWords, 3>> {
        if (!data_ || (which != Configuration::Reference && which != Configuration::Current))
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
        const auto point = lookup(data_->point_lookup.values(), id.value());
        if (!point)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
        return std::array{position(*data_, *point, which, 0), position(*data_, *point, which, 1),
                          position(*data_, *point, which, 2)};
    });
}
AeroResult<ScalarWords> OwnedGeometry::core_radius() const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
    return data_->radius;
}
AeroResult<PhysicalFaceId> OwnedGeometry::facet_face(FacetId id) const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto index = lookup(data_->facet_lookup.values(), id.value());
    if (!index)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto face = data_->facets.values()[static_cast<std::size_t>(*index)].face;
    return PhysicalFaceId{data_->faces.values()[static_cast<std::size_t>(face)].id};
}
AeroResult<BandId> OwnedGeometry::face_band(PhysicalFaceId id) const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto index = lookup(data_->face_lookup.values(), id.value());
    if (!index)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    return BandId{data_->faces.values()[static_cast<std::size_t>(*index)].band};
}
AeroResult<std::array<PointId, 3>> OwnedGeometry::facet_points(FacetId id) const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto index = lookup(data_->facet_lookup.values(), id.value());
    if (!index)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto& p = data_->facets.values()[static_cast<std::size_t>(*index)].points;
    return std::array{PointId{data_->points.values()[static_cast<std::size_t>(p[0])].id},
                      PointId{data_->points.values()[static_cast<std::size_t>(p[1])].id},
                      PointId{data_->points.values()[static_cast<std::size_t>(p[2])].id}};
}
AeroResult<std::array<EdgeId, 3>> OwnedGeometry::facet_edges(FacetId id) const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto index = lookup(data_->facet_lookup.values(), id.value());
    if (!index)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto& e = data_->facets.values()[static_cast<std::size_t>(*index)].edges;
    return std::array{EdgeId{data_->edges.values()[static_cast<std::size_t>(e[0])].id},
                      EdgeId{data_->edges.values()[static_cast<std::size_t>(e[1])].id},
                      EdgeId{data_->edges.values()[static_cast<std::size_t>(e[2])].id}};
}
AeroResult<std::array<std::int8_t, 3>> OwnedGeometry::incidence(FacetId id) const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto index = lookup(data_->facet_lookup.values(), id.value());
    if (!index)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    return data_->facets.values()[static_cast<std::size_t>(*index)].incidence;
}
AeroResult<FacetScatter> OwnedGeometry::scatter(FacetId id) const noexcept {
    if (!data_)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto index = lookup(data_->facet_lookup.values(), id.value());
    if (!index)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    const auto& origin = data_->facets.values()[static_cast<std::size_t>(*index)].origin;
    if (std::holds_alternative<IndependentOrigin>(origin))
        return FacetScatter{IndependentP1Scatter{{1, 0, 0, 0, 1, 0, 0, 0, 1}}};
    const auto* quad = std::get_if<QuadOrigin>(&origin);
    if (quad == nullptr || quad->half > 1)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, id));
    return FacetScatter{quad->half == 0 ? QuadLRScatter{{1, 0, 0, 1, 1, 0}}
                                        : QuadLRScatter{{1, 0, 0, 1, 0, 1}}};
}

AeroResult<QueryRequest> ResolvedQuery::definition() const noexcept {
    return guarded(context_, [&]() -> AeroResult<QueryRequest> {
        if (!active())
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
        const auto& q = query_.values()[0];
        if (const auto* single = std::get_if<FreeSingleStorage>(&q.definition))
            return QueryRequest{QueryId{q.id}, FreePoint{FreeSingle{*single}}};
        else if (const auto* pair = std::get_if<FreeDoubleStorage>(&q.definition))
            return QueryRequest{QueryId{q.id},
                                FreePoint{FreeDouble{{(*pair)[0], (*pair)[1], (*pair)[2]},
                                                     {(*pair)[3], (*pair)[4], (*pair)[5]}}}};
        else if (const auto* vertex = std::get_if<VertexStorage>(&q.definition))
            return QueryRequest{QueryId{q.id}, SourceVertex{PointId{vertex->point}}};
        else if (const auto* edge = std::get_if<EdgeQueryStorage>(&q.definition))
            return QueryRequest{QueryId{q.id},
                                SourceEdgePoint{EdgeId{edge->edge}, restore(edge->t)}};
        else if (const auto* face = std::get_if<FaceQueryStorage>(&q.definition))
            return QueryRequest{QueryId{q.id}, SourceFacePoint{FacetId{face->facet},
                                                               restore(face->u), restore(face->v)}};
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, QueryId{q.id}));
    });
}
AeroResult<CanonicalSupport> ResolvedQuery::support() const noexcept {
    if (!active())
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context_));
    const auto& q = query_.values()[0];
    if (std::holds_alternative<NoMaterialSupport>(q.support))
        return CanonicalSupport{NoMaterialSupport{}};
    if (const auto* p = std::get_if<PointSupportStorage>(&q.support))
        return CanonicalSupport{PointId{p->id}};
    if (const auto* p = std::get_if<EdgeSupportStorage>(&q.support))
        return CanonicalSupport{EdgeId{p->id}};
    if (const auto* p = std::get_if<FacetSupportStorage>(&q.support))
        return CanonicalSupport{FacetId{p->id}};
    return std::unexpected(rejected(core::FailureCode::InvalidInput, context_, QueryId{q.id}));
}

AeroResult<ResolvedQuery>
GeometryEngine::resolve(const OwnedGeometry& geometry, const QueryRequest& request,
                        const GeometryContext& expected, const GeometryLimits& limits,
                        GeometryWorkspace& workspace, core::ResourcePlan& plan,
                        std::stop_token stop) noexcept {
    const GeometryContext context = expected;
    auto lease = workspace.acquire(context.diagnostics);
    if (!lease)
        return std::unexpected(
            OwnedGeometryFailure{complete(lease.error(), context), request.id, context});
    return guarded(context, [&]() -> AeroResult<ResolvedQuery> {
        if (!geometry.data_)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context, request.id));
        const auto& data = *geometry.data_;
        Run run{limits, context, stop};
        if (!finite(context.time))
            return std::unexpected(rejected(core::FailureCode::Nonfinite, context, request.id));
        QueryStorage query;
        query.id = request.id.value();
        auto& arithmetic = workspace.arithmetic_.values()[0];
        // The complete scalar input is checked before any ID/context/range decision.
        if (const auto* free = std::get_if<FreePoint>(&request.definition)) {
            if (const auto* single = std::get_if<FreeSingle>(&free->words)) {
                for (const auto value : single->words)
                    if (!finite(value))
                        return std::unexpected(
                            rejected(core::FailureCode::Nonfinite, context, request.id));
                query.definition.emplace<FreeSingleStorage>(single->words);
            } else if (const auto* pair = std::get_if<FreeDouble>(&free->words)) {
                for (std::size_t i = 0; i < 3; ++i)
                    if (!finite(pair->high[i]) || !finite(pair->low[i]))
                        return std::unexpected(
                            rejected(core::FailureCode::Nonfinite, context, request.id));
                query.definition.emplace<FreeDoubleStorage>(std::array{pair->high[0], pair->high[1],
                                                                       pair->high[2], pair->low[0],
                                                                       pair->low[1], pair->low[2]});
            } else
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, request.id));
        } else if (const auto* vertex = std::get_if<SourceVertex>(&request.definition)) {
            const auto index = lookup(data.point_lookup.values(), vertex->point.value());
            if (!index || data.points.values()[static_cast<std::size_t>(*index)].owner != *index)
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, vertex->point));
            query.definition = VertexStorage{vertex->point.value()};
            query.support = PointSupportStorage{vertex->point.value()};
        } else if (const auto* edge = std::get_if<SourceEdgePoint>(&request.definition)) {
            if (!finite(edge->t))
                return std::unexpected(rejected(core::FailureCode::Nonfinite, context, request.id));
            if (!lookup(data.edge_lookup.values(), edge->edge.value()))
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, edge->edge));
            query.definition = EdgeQueryStorage{edge->edge.value(), store(edge->t)};
        } else if (const auto* face = std::get_if<SourceFacePoint>(&request.definition)) {
            if (!finite(face->u) || !finite(face->v))
                return std::unexpected(rejected(core::FailureCode::Nonfinite, context, request.id));
            if (!lookup(data.facet_lookup.values(), face->facet.value()))
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, face->facet));
            query.definition =
                FaceQueryStorage{face->facet.value(), store(face->u), store(face->v)};
        } else
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context, request.id));
        if (!same_context(data.context, context))
            return std::unexpected(
                rejected(core::FailureCode::StaleGeneration, context, request.id));
        const auto& n = data.counts;
        if (!fits(n.points, limits.counts.points) || !fits(n.edges, limits.counts.edges) ||
            !fits(n.faces, limits.counts.faces) || !fits(n.bands, limits.counts.bands) ||
            !fits(n.triangles, limits.counts.triangles) || !fits(n.quads, limits.counts.quads))
            return std::unexpected(
                rejected(core::FailureCode::ResourceExceeded, context, request.id));
        if (auto r = run.admission(); !r)
            return std::unexpected(r.error());
        if (auto step = run.step(request.id); !step)
            return std::unexpected(step.error());
        arithmetic.weights[0] = 1;
        arithmetic.weights[0] <<= 1074;
        if (const auto* edge = std::get_if<EdgeQueryStorage>(&query.definition)) {
            scalar_into(restore(edge->t), 0, arithmetic.weights[1], arithmetic.word);
            arithmetic.weights[0] -= arithmetic.weights[1];
            if (arithmetic.weights[0] < 0 || arithmetic.weights[1] < 0)
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, request.id));
            const auto index = lookup(data.edge_lookup.values(), edge->edge);
            if (!index)
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, request.id));
            const auto& source = data.edges.values()[static_cast<std::size_t>(*index)];
            if (arithmetic.weights[1] == 0)
                query.support = PointSupportStorage{
                    data.points.values()[static_cast<std::size_t>(source.points[0])].id};
            else if (arithmetic.weights[0] == 0)
                query.support = PointSupportStorage{
                    data.points.values()[static_cast<std::size_t>(source.points[1])].id};
            else
                query.support = EdgeSupportStorage{edge->edge};
        } else if (const auto* face = std::get_if<FaceQueryStorage>(&query.definition)) {
            scalar_into(restore(face->u), 0, arithmetic.weights[1], arithmetic.word);
            scalar_into(restore(face->v), 0, arithmetic.weights[2], arithmetic.word);
            arithmetic.weights[0] -= arithmetic.weights[1];
            arithmetic.weights[0] -= arithmetic.weights[2];
            if (arithmetic.weights[0] < 0 || arithmetic.weights[1] < 0 || arithmetic.weights[2] < 0)
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, request.id));
            const auto index = lookup(data.facet_lookup.values(), face->facet);
            if (!index)
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, request.id));
            const auto& source = data.facets.values()[static_cast<std::size_t>(*index)];
            std::size_t positive = 0, nonzero = 0, zero = 0;
            for (std::size_t i = 0; i < 3; ++i) {
                if (arithmetic.weights[i] > 0) {
                    ++positive;
                    nonzero = i;
                } else
                    zero = i;
            }
            if (positive == 1)
                query.support = PointSupportStorage{
                    data.points.values()[static_cast<std::size_t>(source.points[nonzero])].id};
            else if (positive == 2)
                query.support = EdgeSupportStorage{
                    data.edges.values()[static_cast<std::size_t>(source.edges[(zero + 1) % 3])].id};
            else if (positive == 3)
                query.support = FacetSupportStorage{face->facet};
            else
                return std::unexpected(
                    rejected(core::FailureCode::InvalidInput, context, request.id));
        }
        auto storage = plan.allocate<QueryStorage>(1, operation);
        if (!storage)
            return std::unexpected(
                OwnedGeometryFailure{complete(storage.error(), context), request.id, context});
        storage->values()[0] = query;
        return ResolvedQuery{geometry.data_, std::move(*storage), context, run.steps};
    });
}

AeroResult<void> GeometryEngine::query_position(const GeometryData& data, const QueryStorage& query,
                                                ExactPoints& points,
                                                QueryArithmeticScratch& scratch, Run& run,
                                                ExactQueryCoordinates* result,
                                                Configuration configuration) {
    if (const auto* single = std::get_if<FreeSingleStorage>(&query.definition)) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (auto r = run.step(QueryId{query.id}); !r)
                return r;
            decode_into((*single)[axis], 1074, points[3][axis]);
        }
        return {};
    }
    if (const auto* pair = std::get_if<FreeDoubleStorage>(&query.definition)) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (auto r = run.step(QueryId{query.id}); !r)
                return r;
            decode_into((*pair)[axis], 1074, points[3][axis]);
            decode_into((*pair)[axis + 3], 1074, scratch.word);
            points[3][axis] += scratch.word;
        }
        return {};
    }
    if (const auto* vertex = std::get_if<VertexStorage>(&query.definition)) {
        const auto index = lookup(data.point_lookup.values(), vertex->point);
        if (!index)
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, QueryId{query.id}));
        if (result != nullptr) {
            result->point_count = 1;
            const auto owner = data.points.values()[static_cast<std::size_t>(*index)].owner;
            result->point_ids[0] = data.points.values()[static_cast<std::size_t>(owner)].id;
            result->weights[0] = ExactInteger{1} << 1074;
        }
        return load_point(data, *index, configuration, points[3], scratch, run);
    }
    std::array<core::Index, 3> vertices{};
    std::size_t count = 0;
    scratch.weights[0] = 1;
    scratch.weights[0] <<= 1074;
    if (const auto* edge = std::get_if<EdgeQueryStorage>(&query.definition)) {
        const auto index = lookup(data.edge_lookup.values(), edge->edge);
        if (!index)
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, QueryId{query.id}));
        const auto& source = data.edges.values()[static_cast<std::size_t>(*index)];
        vertices[0] = source.points[0];
        vertices[1] = source.points[1];
        count = 2;
        scalar_into(restore(edge->t), 0, scratch.weights[1], scratch.word);
        scratch.weights[0] -= scratch.weights[1];
    } else if (const auto* face = std::get_if<FaceQueryStorage>(&query.definition)) {
        const auto index = lookup(data.facet_lookup.values(), face->facet);
        if (!index)
            return std::unexpected(
                rejected(core::FailureCode::InvalidInput, run.context, QueryId{query.id}));
        vertices = data.facets.values()[static_cast<std::size_t>(*index)].points;
        count = 3;
        scalar_into(restore(face->u), 0, scratch.weights[1], scratch.word);
        scalar_into(restore(face->v), 0, scratch.weights[2], scratch.word);
        scratch.weights[0] -= scratch.weights[1];
        scratch.weights[0] -= scratch.weights[2];
    } else
        return std::unexpected(
            rejected(core::FailureCode::InvalidInput, run.context, QueryId{query.id}));
    if (result != nullptr) {
        result->point_count = count;
        for (std::size_t i = 0; i < count; ++i) {
            if (auto charged = run.step(QueryId{query.id}); !charged)
                return charged;
            const auto owner = data.points.values()[static_cast<std::size_t>(vertices[i])].owner;
            result->point_ids[i] = data.points.values()[static_cast<std::size_t>(owner)].id;
            result->weights[i] = scratch.weights[i];
        }
    }
    for (std::size_t axis = 0; axis < 3; ++axis) {
        scratch.sum = 0;
        for (std::size_t i = 0; i < count; ++i) {
            if (auto r = run.step(QueryId{query.id}); !r)
                return r;
            // Decode at delta, multiply complete exact word sums, then add at
            // delta^2. No rounded derived point or independent third weight exists.
            scalar_into(position(data, vertices[i], configuration, axis), 0, scratch.coordinate,
                        scratch.word);
            scratch.product = scratch.weights[i] * scratch.coordinate;
            scratch.sum += scratch.product;
        }
        points[3][axis] = scratch.sum;
    }
    return {};
}

AeroResult<LocalGeometryInput>
GeometryEngine::local(const OwnedGeometry& geometry, const ResolvedQuery& query, FacetId target,
                      const GeometryLimits& limits, GeometryWorkspace& workspace,
                      const PredicateLease& lease, core::ResourcePlan& plan, std::stop_token stop,
                      GeometryWorkAccounting accounting, TriangleOutcome* transient) noexcept {
    const auto context = geometry.context_;
    return guarded(context, [&]() -> AeroResult<LocalGeometryInput> {
        if (!geometry.data_ || !query.active())
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
        const auto& data = *geometry.data_;
        const auto& q = query.query_.values()[0];
        const auto index = lookup(data.facet_lookup.values(), target.value());
        if (!index)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context, target));
        if (data.identity != query.geometry_->identity || !same_context(context, query.context_))
            return std::unexpected(
                rejected(core::FailureCode::StaleGeneration, context, QueryId{q.id}));
        Run run{limits, context, stop};
        const auto& n = data.counts;
        if (!fits(n.points, limits.counts.points) || !fits(n.edges, limits.counts.edges) ||
            !fits(n.faces, limits.counts.faces) || !fits(n.bands, limits.counts.bands) ||
            !fits(n.triangles, limits.counts.triangles) || !fits(n.quads, limits.counts.quads))
            return std::unexpected(
                rejected(core::FailureCode::ResourceExceeded, context, QueryId{q.id}));
        if (auto r = run.admission(); !r)
            return std::unexpected(r.error());
        core::PlannedBuffer<TriangleOutcome> owned;
        if (transient == nullptr) {
            auto added = plan.allocate<TriangleOutcome>(1, operation);
            if (!added)
                return std::unexpected(
                    OwnedGeometryFailure{complete(added.error(), context), QueryId{q.id}, context});
            owned = std::move(*added);
        }
        auto slots =
            ExactPredicateAccess::coordinates(*workspace.predicate_, lease, context.diagnostics);
        if (!slots)
            return std::unexpected(
                OwnedGeometryFailure{complete(slots.error(), context), QueryId{q.id}, context});
        if (auto r = query_position(data, q, slots->get(), workspace.arithmetic_.values()[0], run);
            !r)
            return std::unexpected(r.error());
        const auto& facet = data.facets.values()[static_cast<std::size_t>(*index)];
        for (std::size_t vertex = 0; vertex < 3; ++vertex)
            if (auto r = load_point(data, facet.points[vertex], Configuration::Current,
                                    slots->get()[vertex], workspace.arithmetic_.values()[0], run);
                !r)
                return std::unexpected(r.error());
        auto predicate_limits = limits.predicate;
        if (accounting == GeometryWorkAccounting::TotalAndPredicateSubcap) {
            if (run.steps > limits.max_steps)
                return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
            predicate_limits.max_steps =
                std::min(predicate_limits.max_steps, limits.max_steps - run.steps);
        }
        std::optional<ExactTriangleResult> result;
        std::uint64_t predicate_steps = 0;
        if (transient != nullptr) {
            const auto used =
                ExactPredicateAccess::finish_into(*workspace.predicate_, lease, *transient,
                                                  predicate_limits, stop, context.diagnostics, 0);
            if (!used)
                return std::unexpected(
                    OwnedGeometryFailure{complete(used.error(), context), target, context});
            predicate_steps = *used;
        } else {
            auto classified =
                ExactPredicateAccess::finish(*workspace.predicate_, lease, std::move(owned),
                                             predicate_limits, stop, context.diagnostics, 0);
            if (!classified)
                return std::unexpected(
                    OwnedGeometryFailure{complete(classified.error(), context), target, context});
            predicate_steps = classified->steps();
            result.emplace(std::move(*classified));
        }
        std::uint64_t work = 0;
        if (accounting == GeometryWorkAccounting::TotalAndPredicateSubcap) {
            if (predicate_steps > limits.max_steps - run.steps)
                return std::unexpected(rejected(core::FailureCode::ResourceExceeded, context));
            work = run.steps + predicate_steps;
        }
        if (result)
            return LocalGeometryInput{geometry.data_,     q,   *index, slots->get(),
                                      std::move(*result), work};
        if (transient == nullptr)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
        return LocalGeometryInput{geometry.data_,        q,   *index, slots->get(),
                                  std::cref(*transient), work};
    });
}

AeroResult<ExactTriangleResult>
GeometryEngine::classify(const OwnedGeometry& geometry, const ResolvedQuery& query, FacetId target,
                         const GeometryLimits& limits, GeometryWorkspace& workspace,
                         core::ResourcePlan& plan, std::stop_token stop) noexcept {
    const auto context = geometry.context_;
    auto lease = workspace.acquire(context.diagnostics);
    if (!lease)
        return std::unexpected(OwnedGeometryFailure{complete(lease.error(), context), {}, context});
    auto result = local(geometry, query, target, limits, workspace, *lease, plan, std::move(stop),
                        GeometryWorkAccounting::GeometryAndPredicateBudgets);
    if (!result)
        return std::unexpected(result.error());
    auto* owned = std::get_if<ExactTriangleResult>(&result->classification);
    if (owned == nullptr)
        return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
    return std::move(*owned);
}

OwnedGeometry GeometryPacketAccess::retain(const LocalGeometryInput& input) noexcept {
    return GeometryEngine::retain(input);
}

AeroResult<LocalGeometryInput> GeometryPacketAccess::load(
    const OwnedGeometry& geometry, const ResolvedQuery& query, FacetId target,
    const GeometryLimits& limits, GeometryWorkspace& workspace, const PredicateLease& lease,
    core::ResourcePlan& plan, std::stop_token stop, TriangleOutcome* transient) noexcept {
    return GeometryEngine::local(geometry, query, target, limits, workspace, lease, plan,
                                 std::move(stop), GeometryWorkAccounting::TotalAndPredicateSubcap,
                                 transient);
}
AeroResult<OwnedGeometry> prepare_owned_geometry(const GeometryInput& input,
                                                 const GeometryLimits& limits,
                                                 GeometryWorkspace& workspace,
                                                 core::ResourcePlan& plan,
                                                 std::stop_token stop) noexcept {
    return GeometryEngine::prepare(input, limits, workspace, plan, std::move(stop));
}
AeroResult<ResolvedQuery> resolve_query(const OwnedGeometry& geometry, const QueryRequest& query,
                                        const GeometryContext& expected,
                                        const GeometryLimits& limits, GeometryWorkspace& workspace,
                                        core::ResourcePlan& plan, std::stop_token stop) noexcept {
    return GeometryEngine::resolve(geometry, query, expected, limits, workspace, plan,
                                   std::move(stop));
}
AeroResult<ExactTriangleResult>
classify_resolved_query(const OwnedGeometry& geometry, const ResolvedQuery& query, FacetId target,
                        const GeometryLimits& limits, GeometryWorkspace& workspace,
                        core::ResourcePlan& plan, std::stop_token stop) noexcept {
    return GeometryEngine::classify(geometry, query, target, limits, workspace, plan,
                                    std::move(stop));
}

AeroResult<ExactQueryCoordinates>
GeometryEngine::coordinates(const OwnedGeometry& geometry, const ResolvedQuery& query,
                            const GeometryLimits& limits, GeometryWorkspace& workspace,
                            std::stop_token stop, Configuration configuration) noexcept {
    const auto context = geometry.context_;
    return guarded(context, [&]() -> AeroResult<ExactQueryCoordinates> {
        if (configuration != Configuration::Reference && configuration != Configuration::Current)
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
        if (!geometry.data_ || !query.active())
            return std::unexpected(rejected(core::FailureCode::InvalidInput, context));
        const auto& data = *geometry.data_;
        const auto& q = query.query_.values().front();
        if (data.identity != query.geometry_->identity || !same_context(context, query.context_))
            return std::unexpected(
                rejected(core::FailureCode::StaleGeneration, context, QueryId{q.id}));
        Run run{limits, context, stop};
        if (auto admitted = run.admission(); !admitted)
            return std::unexpected(admitted.error());
        const auto& n = data.counts;
        if (!fits(n.points, limits.counts.points) || !fits(n.edges, limits.counts.edges) ||
            !fits(n.faces, limits.counts.faces) || !fits(n.bands, limits.counts.bands) ||
            !fits(n.triangles, limits.counts.triangles) || !fits(n.quads, limits.counts.quads))
            return std::unexpected(
                rejected(core::FailureCode::ResourceExceeded, context, QueryId{q.id}));
        auto lease = workspace.acquire(context.diagnostics);
        if (!lease)
            return std::unexpected(
                OwnedGeometryFailure{complete(lease.error(), context), QueryId{q.id}, context});
        auto slots =
            ExactPredicateAccess::coordinates(*workspace.predicate_, *lease, context.diagnostics);
        if (!slots)
            return std::unexpected(
                OwnedGeometryFailure{complete(slots.error(), context), QueryId{q.id}, context});
        ExactQueryCoordinates coordinates{};
        if (auto loaded = query_position(data, q, slots->get(), workspace.arithmetic_.values()[0],
                                         run, &coordinates, configuration);
            !loaded)
            return std::unexpected(loaded.error());
        if (auto published = run.step(QueryId{q.id}); !published)
            return std::unexpected(published.error());
        coordinates.numerator = slots->get()[3];
        coordinates.used_steps = run.steps;
        return coordinates;
    });
}

AeroResult<ExactQueryCoordinates>
exact_query_coordinates(const OwnedGeometry& geometry, const ResolvedQuery& query,
                        const GeometryLimits& limits, GeometryWorkspace& workspace,
                        std::stop_token stop, Configuration configuration) noexcept {
    return GeometryEngine::coordinates(geometry, query, limits, workspace, std::move(stop),
                                       configuration);
}

GeometryWorkspace::GeometryWorkspace(Key, std::unique_ptr<PredicateWorkspace> predicate,
                                     core::PlannedBuffer<QueryArithmeticScratch> scratch) noexcept
    : predicate_(std::move(predicate)), arithmetic_(std::move(scratch)) {}
GeometryWorkspace::~GeometryWorkspace() = default;
core::Result<std::unique_ptr<GeometryWorkspace>>
GeometryWorkspace::create(core::ResourcePlan& plan, core::FailureDetails details) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::PrecisionUnresolved, details,
        [&]() -> core::Result<std::unique_ptr<GeometryWorkspace>> {
            auto predicate = PredicateWorkspace::create(plan, details);
            if (!predicate)
                return std::unexpected(predicate.error());
            auto arithmetic = plan.allocate<QueryArithmeticScratch>(1, operation);
            if (!arithmetic)
                return std::unexpected(complete_details(arithmetic.error(), details));
            return std::make_unique<GeometryWorkspace>(Key{}, std::move(*predicate),
                                                       std::move(*arithmetic));
        });
}
core::Result<PredicateLease> GeometryWorkspace::acquire(core::FailureDetails details) noexcept {
    return predicate_->acquire(details);
}
std::uint64_t GeometryWorkspace::payload_bytes() noexcept {
    return PredicateWorkspace::payload_bytes() + sizeof(QueryArithmeticScratch);
}
ResolvedQuery::ResolvedQuery(std::shared_ptr<const GeometryData> data,
                             core::PlannedBuffer<QueryStorage> query, GeometryContext stamp,
                             std::uint64_t steps) noexcept
    : geometry_(std::move(data)), query_(std::move(query)), context_(stamp), steps_(steps) {}
ResolvedQuery::ResolvedQuery(ResolvedQuery&&) noexcept = default;
ResolvedQuery& ResolvedQuery::operator=(ResolvedQuery&&) noexcept = default;
ResolvedQuery::~ResolvedQuery() = default;
bool ResolvedQuery::active() const noexcept {
    return geometry_ != nullptr && query_.size() == 1;
}
std::uint64_t ResolvedQuery::payload_bytes() noexcept {
    return sizeof(QueryStorage);
}
} // namespace uvlas::aero::detail
