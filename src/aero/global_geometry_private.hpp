#pragma once
#include "global_geometry_storage.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>

namespace uvlas::aero::detail {
class GlobalEngine {
public:
    static GlobalResult<GlobalTopology> build(const OwnedGeometry&, const GlobalTopologyInput&,
                                              GlobalLimits, GlobalWorkspace&, core::ResourcePlan&,
                                              std::stop_token) noexcept;
    static GlobalResult<GlobalContactState> check(const GlobalTopology&, const OwnedGeometry&,
                                                  GlobalLimits, GlobalWorkspace&,
                                                  core::ResourcePlan&, std::stop_token) noexcept;
    static GlobalResult<AliasKinematicsState> kinematics(const OwnedGeometry&,
                                                         const AliasKinematicsInput&, GlobalLimits,
                                                         GlobalWorkspace&, core::ResourcePlan&,
                                                         std::stop_token) noexcept;
    static const std::shared_ptr<const GeometryData>& source(const OwnedGeometry& g) noexcept {
        return g.data_;
    }
    static OwnedGeometry retain(const OwnedGeometry& g) noexcept {
        return OwnedGeometry{g.data_, g.context_};
    }
};
namespace global_detail {
constexpr auto operation = core::Operation::EvaluateAero;
constexpr std::uint64_t integer_limbs = 1024;
constexpr std::uint64_t simple_cost = integer_limbs;
// Three cleared arrays, two exports, schoolbook products/carries and result import.
constexpr std::uint64_t multiply_cost = integer_limbs * integer_limbs + 12 * integer_limbs;
constexpr std::uint64_t decode_cost = 8 * simple_cost;
// Each coordinate comparison also adds at most two low words to their highs.
constexpr std::uint64_t node_cost = 18 * simple_cost + 24 * decode_cost + 5;
// Static loop expansion, including paths that cannot all occur together:
// normals/plane values/direction 36 products; six inside tests (3 edges,2 products);
// nine 2D crossings (12 products); four cut points (6 products, doubled conservatively);
// endpoint order/equality 16*2 products; nine possible shared edges *24*11;
// nine possible shared vertices *24*3. Each product site permits eight simple
// adds/sign/order/loop operations. Coordinate loads include all nine edge matches.
constexpr std::uint64_t pair_products =
    36 + 6 * 3 * 2 + 9 * 12 + 2 * 4 * 6 + 16 * 2 + 9 * 24 * 11 + 9 * 24 * 3 +
    (15 * 14 / 2) * 6; // Exact candidate dedup: at most 105 pairs, 3 coordinates each.
constexpr std::uint64_t narrow_cost = pair_products * (multiply_cost + 8 * simple_cost) +
                                      (18 + 9 * 6) * (2 * decode_cost + simple_cost);
struct Abort {
    core::FailureCode code = core::FailureCode::InvalidInput;
    std::optional<GeometryProblem> problem;
    std::optional<FacetId> other;
    std::optional<core::DofId> dof;
};
[[noreturn]] inline void reject(core::FailureCode code, std::optional<GeometryProblem> p = {},
                                std::optional<FacetId> other = {}) {
    throw Abort{code, p, other, {}};
}
[[noreturn]] inline void reject(core::FailureCode code, core::DofId dof) {
    throw Abort{code, {}, {}, dof};
}
inline std::uint64_t plus(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        reject(core::FailureCode::ResourceExceeded);
    return a + b;
}
inline std::uint64_t times(std::uint64_t a, std::uint64_t b) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        reject(core::FailureCode::ResourceExceeded);
    return a * b;
}
inline core::Failure with_details(core::Failure failure,
                                  const core::FailureDetails& origin) noexcept {
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
inline GlobalGeometryFailure failure(core::Failure cause, const GeometryContext& c,
                                     std::optional<GeometryProblem> problem = {},
                                     std::optional<FacetId> other = {}) noexcept {
    auto& d = cause.details;
    const auto& p = c.diagnostics;
    if (!d.model)
        d.model = p.model ? p.model : std::optional{c.model};
    if (!d.parent)
        d.parent = p.parent;
    if (!d.trial)
        d.trial = p.trial;
    if (!d.numerical_generation)
        d.numerical_generation =
            p.numerical_generation ? p.numerical_generation : std::optional{c.numerical};
    if (!d.component)
        d.component = p.component;
    if (!d.problem)
        d.problem = p.problem;
    if (!d.resources)
        d.resources = p.resources;
    if (!d.last_accepted)
        d.last_accepted = p.last_accepted;
    return {{cause, problem, c}, other};
}
inline GlobalGeometryFailure failure(core::FailureCode code, const GeometryContext& c,
                                     std::optional<GeometryProblem> p = {},
                                     std::optional<FacetId> other = {}) noexcept {
    return failure(core::Failure{code, operation}, c, p, other);
}
template <class F>
auto guarded(const GeometryContext& c, F&& action) noexcept -> decltype(action()) {
    try {
        return action();
    } catch (const Abort& e) {
        auto result = failure(e.code, c, e.problem, e.other);
        if (e.dof)
            result.origin.cause.details.problem = *e.dof;
        return std::unexpected(result);
    } catch (const GlobalGeometryFailure& e) {
        return std::unexpected(e);
    } catch (const std::bad_alloc&) {
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, c));
    } catch (...) {
        return std::unexpected(failure(core::FailureCode::PrecisionUnsupported, c));
    }
}
struct Run {
    GlobalLimits limits;
    GeometryContext context;
    GlobalScratch& scratch;
    std::stop_token stop;
    std::uint64_t reserved = 0, used = 0;
    void admit(std::uint64_t n) {
        if (stop.stop_requested())
            reject(core::FailureCode::Cancelled);
        const auto total = plus(reserved, n);
        if (total > limits.max_steps) {
            core::FailureDetails details;
            details.resources =
                core::ResourceUsage{core::ResourceKind::Iterations, total, reserved};
            throw failure(core::Failure{core::FailureCode::ResourceExceeded, operation, details},
                          context);
        }
        reserved = total;
    }
    void step(std::uint64_t n = 1) {
        if (stop.stop_requested())
            reject(core::FailureCode::Cancelled);
        if (n > reserved - used) {
            core::FailureDetails details;
            details.resources =
                core::ResourceUsage{core::ResourceKind::Iterations, plus(used, n), used};
            throw failure(core::Failure{core::FailureCode::ResourceExceeded, operation, details},
                          context);
        }
        used += n;
    }
};
struct Lease {
    std::atomic_flag& flag;
    explicit Lease(std::atomic_flag& f) : flag(f) {
        if (flag.test_and_set(std::memory_order_acquire))
            reject(core::FailureCode::InvalidInput);
    }
    ~Lease() {
        flag.clear(std::memory_order_release);
    }
};
template <class T>
void allocate(core::PlannedBuffer<T>& target, std::uint64_t n, core::ResourcePlan& plan,
              const GeometryContext& c) {
    if (n > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        reject(core::FailureCode::ResourceExceeded);
    auto result = plan.allocate<T>(static_cast<core::Index>(n), operation);
    if (!result)
        throw failure(result.error(), c);
    target = std::move(*result);
}
template <class T>
std::uint64_t bytes(std::uint64_t n) {
    return times(n, sizeof(T));
}
inline void byte_admit(std::uint64_t amount, core::ResourcePlan& plan,
                       const GeometryContext& context) {
    if (amount > plan.remaining_bytes()) {
        core::FailureDetails details;
        const auto committed = plan.byte_limit() - plan.remaining_bytes();
        if (amount <= std::numeric_limits<std::uint64_t>::max() - committed)
            details.resources =
                core::ResourceUsage{core::ResourceKind::Bytes, committed + amount, committed};
        throw failure(core::Failure{core::FailureCode::ResourceExceeded, operation, details},
                      context);
    }
}
inline std::uint64_t metadata_admit(numerics::linear::BackendMemoryLimit limit,
                                    std::uint64_t estimate) {
    const auto selection = limit.selection();
    if (std::holds_alternative<numerics::linear::RequireProvenBound>(selection) && estimate != 0)
        reject(core::FailureCode::Unsupported);
    if (estimate == 0)
        return 0;
    if (estimate > std::get<numerics::linear::EstimateAdmission>(selection).bytes)
        reject(core::FailureCode::ResourceExceeded);
    return estimate;
}
inline bool finite(double x) {
    return ((std::bit_cast<std::uint64_t>(x) >> 52U) & 2047U) != 2047U;
}
inline ExactInteger word(double x, unsigned extra, Run& r) {
    r.step(decode_cost);
    if (!finite(x))
        reject(core::FailureCode::Nonfinite);
    const auto bits = std::bit_cast<std::uint64_t>(x);
    const auto exponent = static_cast<unsigned>((bits >> 52U) & 2047U);
    auto significand = bits & UINT64_C(0xfffffffffffff);
    if (exponent != 0)
        significand |= UINT64_C(0x10000000000000);
    ExactInteger value = significand;
    value <<= (exponent == 0 ? 0U : exponent - 1U) + extra;
    return (bits >> 63U) != 0 ? -value : value;
}
inline ExactInteger scalar(ScalarWords w, unsigned extra, Run& r) {
    if (const auto* v = std::get_if<SingleWord>(&w))
        return word(v->value, extra, r);
    const auto& v = std::get<DoubleWord>(w);
    auto a = word(v.high, extra, r);
    auto b = word(v.low, extra, r);
    r.step(simple_cost);
    return a + b;
}
inline ParameterStorage store(ScalarWords w) {
    if (const auto* v = std::get_if<SingleWord>(&w))
        return v->value;
    const auto& v = std::get<DoubleWord>(w);
    return std::array{v.high, v.low};
}
inline ScalarWords restore(const ParameterStorage& w) {
    if (const auto* v = std::get_if<double>(&w))
        return SingleWord{*v};
    const auto& v = std::get<std::array<double, 2>>(w);
    return DoubleWord{v[0], v[1]};
}
inline ExactInteger coordinate(const GeometryData& g, core::Index p, std::size_t k,
                               Configuration config, Run& r) {
    const auto i = static_cast<std::size_t>(p) * 3 + k;
    const auto high =
        config == Configuration::Reference ? g.reference_high.values() : g.current_high.values();
    const auto low =
        config == Configuration::Reference ? g.reference_low.values() : g.current_low.values();
    auto x = word(high[i], 1074, r);
    if (!low.empty()) {
        auto y = word(low[i], 1074, r);
        r.step(simple_cost);
        x += y;
    }
    return x;
}
inline ExactInteger add(const ExactInteger& a, const ExactInteger& b, Run& r) {
    r.step(simple_cost);
    return a + b;
}
inline ExactInteger sub(const ExactInteger& a, const ExactInteger& b, Run& r) {
    r.step(simple_cost);
    return a - b;
}
inline int compare(const ExactInteger& a, const ExactInteger& b, Run& r) {
    r.step(simple_cost);
    return a < b ? -1 : a > b ? 1 : 0;
}
ExactInteger multiply(const ExactInteger&, const ExactInteger&, Run&);
Vector subtract(const Vector&, const Vector&, Run&);
Vector cross(const Vector&, const Vector&, Run&);
ExactInteger dot(const Vector&, const Vector&, Run&);
bool intersects_outside_shared(const GeometryData&, const GlobalTopologyData&, std::uint64_t,
                               std::uint64_t, Run&);
void build_tree(GlobalTopologyData&, Run&);
void refit_tree(const GeometryData&, std::span<Node>, Run&);
struct TraversalResult {
    std::uint64_t candidates = 0, steps = 0;
};
TraversalResult traverse(const GeometryData&, const GlobalTopologyData&, std::span<const Node>,
                         std::span<NodePair>, bool, Run&);
std::uint64_t tree_shape(std::uint64_t);
std::uint64_t pair_stack(std::uint64_t);
inline std::size_t lookup(std::span<const IdLookup> table, std::uint64_t id, Run& r,
                          std::optional<GeometryProblem> problem = {}) {
    std::size_t lo = 0, hi = table.size();
    while (lo < hi) {
        r.step();
        const auto mid = lo + (hi - lo) / 2;
        if (table[mid].id < id)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo == table.size() || table[lo].id != id)
        reject(core::FailureCode::InvalidInput, problem);
    return static_cast<std::size_t>(table[lo].ordinal);
}
} // namespace global_detail
} // namespace uvlas::aero::detail
