#pragma once

#include <boost/multiprecision/cpp_int.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stop_token>
#include <utility>
#include <uvlas/core/resource_plan.hpp>
#include <uvlas/numerics/backend_memory.hpp>
#include <variant>

namespace uvlas::aero::detail {
/// Internal preparation integer; not a public mathematical API.
using ExactInteger = boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<32768, 32768, boost::multiprecision::signed_magnitude,
                                           boost::multiprecision::checked, void>,
    boost::multiprecision::et_off>;
/// Packed x0.xyz, x1.xyz, x2.xyz, query.xyz, in SI metres.
struct Binary64Coordinates {
    std::span<const double> values;
};
/// The same twelve coordinates as exact sums of corresponding finite words.
struct TwoWordCoordinates {
    std::span<const double> high, low;
};
/// Distinct representations; an empty low plane does not silently change precision.
using CoordinateView = std::variant<Binary64Coordinates, TwoWordCoordinates>;
/// Side of the positive oriented triangle plane.
enum class PlaneSide : std::uint8_t { Below, On, Above };
/// Feature containing the orthogonal projection; Outside is the finite-face exterior.
enum class ProjectionFeature : std::uint8_t {
    Outside,
    Face,
    Edge01,
    Edge12,
    Edge20,
    Vertex0,
    Vertex1,
    Vertex2
};
/// Relation of distance squared to sixteen times the longest edge squared.
enum class FarRelation : std::uint8_t { Near, Boundary, Far };
/// Exact mathematical degeneracy, without a physical geometry-policy decision.
struct DegenerateTriangle {};
/// An exact nonnegative squared distance in SI m²: numerator/denominator * 2^unit_exponent.
struct ExactSquaredDistance {
    ExactInteger numerator, denominator;
    int unit_exponent;
    ExactSquaredDistance(ExactInteger n, ExactInteger d, int exponent)
        : numerator(std::move(n)), denominator(std::move(d)), unit_exponent(exponent) {}
};
/// Mathematical relation only; no physical incidence or global geometry qualification.
struct NondegenerateTriangle {
    PlaneSide plane;
    ProjectionFeature projection;
    ExactSquaredDistance distance_squared;
    ExactInteger longest_edge_squared; ///< Integer in the same 2^-4296 m² units.
    std::uint8_t minimizing_features;  ///< Bits 0..2 vertices, 3..5 finite edges, 6 face.
    FarRelation far;
    NondegenerateTriangle(PlaneSide side, ProjectionFeature support, ExactSquaredDistance distance,
                          ExactInteger longest, std::uint8_t features, FarRelation relation)
        : plane(side),
          projection(support),
          distance_squared(std::move(distance)),
          longest_edge_squared(std::move(longest)),
          minimizing_features(features),
          far(relation) {}
};
using TriangleOutcome = std::variant<DegenerateTriangle, NondegenerateTriangle>;
/// Explicit work and external forecast limits; zero is a real zero budget.
struct PredicateLimits {
    numerics::linear::BackendMemoryLimit external_memory;
    std::uint64_t max_steps;
    PredicateLimits(numerics::linear::BackendMemoryLimit memory, std::uint64_t steps) noexcept
        : external_memory(memory), max_steps(steps) {}
};
/// Aero-only factory plan for external Boost scratch; not a payload/RSS hard cap.
class PredicateMemoryPlan {
public:
    [[nodiscard]] std::uint64_t peak_payload_estimate_bytes() const noexcept {
        return estimate_;
    }
    [[nodiscard]] std::optional<std::uint64_t>
    proven_peak_payload_upper_bound_bytes() const noexcept {
        return std::nullopt;
    }

private:
    friend PredicateMemoryPlan predicate_memory_plan() noexcept;
    explicit PredicateMemoryPlan(std::uint64_t estimate) noexcept : estimate_(estimate) {}
    std::uint64_t estimate_;
};
struct PredicateScratch;
class PredicateLease;
/// Owns actual planned arithmetic scratch and an exclusive lease state.
class PredicateWorkspace {
    class ConstructionKey {
        friend class PredicateWorkspace;
        ConstructionKey() = default;
    };

public:
    /// Factory-only construction key prevents callers from bypassing planned storage.
    explicit PredicateWorkspace(ConstructionKey,
                                core::PlannedBuffer<PredicateScratch> scratch) noexcept;
    /// Allocate real scratch payload through plan; metadata allocation is separate overhead.
    [[nodiscard]] static core::Result<std::unique_ptr<PredicateWorkspace>>
    create(core::ResourcePlan& plan, core::FailureDetails details = {}) noexcept;
    PredicateWorkspace(const PredicateWorkspace&) = delete;
    PredicateWorkspace& operator=(const PredicateWorkspace&) = delete;
    PredicateWorkspace(PredicateWorkspace&&) = delete;
    PredicateWorkspace& operator=(PredicateWorkspace&&) = delete;
    ~PredicateWorkspace();
    /// Acquire before reading mutable state; nested/concurrent use returns InvalidInput.
    [[nodiscard]] core::Result<PredicateLease> acquire(core::FailureDetails details = {}) noexcept;
    /// Exact declared payload charged for one workspace, excluding owner metadata.
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;

private:
    friend class PredicateLease;
    friend core::Result<class ExactTriangleResult>
    classify_triangle(CoordinateView, PredicateWorkspace&, core::ResourcePlan&, PredicateLimits,
                      std::stop_token, core::FailureDetails) noexcept;
    core::PlannedBuffer<PredicateScratch> scratch_;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};
/// A genuine exclusive borrow of a prepared predicate workspace.
class PredicateLease {
public:
    PredicateLease(PredicateLease&& other) noexcept;
    PredicateLease& operator=(PredicateLease&& other) noexcept;
    PredicateLease(const PredicateLease&) = delete;
    PredicateLease& operator=(const PredicateLease&) = delete;
    ~PredicateLease();

private:
    friend class PredicateWorkspace;
    explicit PredicateLease(PredicateWorkspace& owner) noexcept : owner_(&owner) {}
    PredicateWorkspace* owner_;
};
/// Owns its outcome independently of the inputs, workspace and resource-plan lifetime.
class ExactTriangleResult {
public:
    ExactTriangleResult(ExactTriangleResult&&) noexcept = default;
    ExactTriangleResult& operator=(ExactTriangleResult&&) noexcept = default;
    ExactTriangleResult(const ExactTriangleResult&) = delete;
    ExactTriangleResult& operator=(const ExactTriangleResult&) = delete;
    ~ExactTriangleResult() = default;
    /// Borrow a read-only outcome; null after moving from this owner.
    [[nodiscard]] const TriangleOutcome* outcome() const& noexcept;
    [[nodiscard]] const TriangleOutcome* outcome() const&& = delete;
    /// Copy the explicit admitted limit; it is not a measured external allocation.
    [[nodiscard]] numerics::linear::BackendMemoryLimit admitted_limit() const noexcept {
        return limit_;
    }
    /// Copy the factory-generated external scratch forecast, with no proven bound.
    [[nodiscard]] PredicateMemoryPlan memory_plan() const noexcept;
    /// Number of named preparation stages charged before performing their arithmetic.
    [[nodiscard]] std::uint64_t steps() const noexcept {
        return steps_;
    }
    /// Reporting identities copied at the entry boundary; not a membership proof.
    [[nodiscard]] core::FailureDetails diagnostic_context() const noexcept {
        return details_;
    }
    /// Exact payload charged for the independently owned result.
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept {
        return sizeof(TriangleOutcome);
    }

private:
    friend core::Result<ExactTriangleResult> classify_triangle(CoordinateView, PredicateWorkspace&,
                                                               core::ResourcePlan&, PredicateLimits,
                                                               std::stop_token,
                                                               core::FailureDetails) noexcept;
    ExactTriangleResult(core::PlannedBuffer<TriangleOutcome> outcome,
                        numerics::linear::BackendMemoryLimit limit, std::uint64_t steps,
                        core::FailureDetails details) noexcept
        : outcome_(std::move(outcome)), limit_(limit), steps_(steps), details_(details) {}
    core::PlannedBuffer<TriangleOutcome> outcome_;
    numerics::linear::BackendMemoryLimit limit_;
    std::uint64_t steps_;
    core::FailureDetails details_;
};
/// Fixed Boost multiplication-scratch forecast; proof is deliberately absent.
[[nodiscard]] PredicateMemoryPlan predicate_memory_plan() noexcept;
/// Classify the actual twelve coordinate values, without a near-contact tolerance.
/// A degenerate triangle is a mathematical outcome, not a physically accepted geometry.
/// The result's distance is exact and never silently projected to binary64.
/// Boost's external scratch is forecast-admitted preparation work, not a no-allocation kernel.
[[nodiscard]] core::Result<ExactTriangleResult>
classify_triangle(CoordinateView coordinates, PredicateWorkspace& workspace,
                  core::ResourcePlan& result_plan, PredicateLimits limits,
                  std::stop_token cancellation = {}, core::FailureDetails details = {}) noexcept;
} // namespace uvlas::aero::detail
