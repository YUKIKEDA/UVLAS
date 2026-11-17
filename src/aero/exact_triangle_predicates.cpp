#include "exact_triangle_private.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <stdexcept>
#include <utility>

namespace uvlas::aero::detail {
namespace {
constexpr auto operation = core::Operation::EvaluateAero;
constexpr int coordinate_unit = -2148;
constexpr int squared_unit = 2 * coordinate_unit;
using Vector = std::array<ExactInteger, 3>;
[[nodiscard]] core::Failure
complete_allocation_failure(core::Failure failure, const core::FailureDetails& context) noexcept {
    auto& inner = failure.details;
    if (!inner.model)
        inner.model = context.model;
    if (!inner.parent)
        inner.parent = context.parent;
    if (!inner.trial)
        inner.trial = context.trial;
    if (!inner.numerical_generation)
        inner.numerical_generation = context.numerical_generation;
    if (!inner.component)
        inner.component = context.component;
    if (!inner.problem)
        inner.problem = context.problem;
    if (!inner.resources)
        inner.resources = context.resources;
    if (!inner.last_accepted)
        inner.last_accepted = context.last_accepted;
    return failure;
}
struct Candidate {
    ExactInteger numerator, denominator;
    Candidate(ExactInteger n, ExactInteger d)
        : numerator(std::move(n)), denominator(std::move(d)) {}
};
[[nodiscard]] Vector subtract(const Vector& a, const Vector& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
[[nodiscard]] ExactInteger dot(const Vector& a, const Vector& b) {
    ExactInteger value = a[0] * b[0];
    value += a[1] * b[1];
    value += a[2] * b[2];
    return value;
}
[[nodiscard]] Vector cross(const Vector& a, const Vector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
[[nodiscard]] bool finite_word(double value) noexcept {
    return ((std::bit_cast<std::uint64_t>(value) >> 52U) & 2047U) != 2047U;
}
[[nodiscard]] ExactInteger decode(double value) {
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto exponent = static_cast<unsigned>((bits >> 52U) & 2047U);
    auto significand = bits & UINT64_C(0x000fffffffffffff);
    if (exponent != 0)
        significand |= UINT64_C(0x0010000000000000);
    // Shift a positive integer first: checked cpp_int forbids shifting negatives.
    ExactInteger integer{significand};
    integer <<= exponent == 0 ? 1074U : exponent + 1073U;
    return (bits >> 63U) != 0 ? -integer : integer;
}
[[nodiscard]] ProjectionFeature projection_feature(const ExactInteger& n0, const ExactInteger& n1,
                                                   const ExactInteger& n2) {
    if (n0 < 0 || n1 < 0 || n2 < 0)
        return ProjectionFeature::Outside;
    if (n1 == 0 && n2 == 0)
        return ProjectionFeature::Vertex0;
    if (n0 == 0 && n2 == 0)
        return ProjectionFeature::Vertex1;
    if (n0 == 0 && n1 == 0)
        return ProjectionFeature::Vertex2;
    if (n2 == 0)
        return ProjectionFeature::Edge01;
    if (n0 == 0)
        return ProjectionFeature::Edge12;
    if (n1 == 0)
        return ProjectionFeature::Edge20;
    return ProjectionFeature::Face;
}
} // namespace
struct PredicateScratch {
    std::array<Vector, 4> points;
    Vector e, f, r, normal;
    std::array<std::optional<Candidate>, 7> candidates;
};
PredicateWorkspace::PredicateWorkspace(ConstructionKey,
                                       core::PlannedBuffer<PredicateScratch> scratch) noexcept
    : scratch_(std::move(scratch)) {}
PredicateWorkspace::~PredicateWorkspace() = default;
std::uint64_t PredicateWorkspace::payload_bytes() noexcept {
    return sizeof(PredicateScratch);
}
core::Result<std::unique_ptr<PredicateWorkspace>>
PredicateWorkspace::create(core::ResourcePlan& plan, core::FailureDetails details) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::PrecisionUnresolved, details,
        [&]() -> core::Result<std::unique_ptr<PredicateWorkspace>> {
            auto storage = plan.allocate<PredicateScratch>(1, operation);
            if (!storage)
                return std::unexpected(complete_allocation_failure(storage.error(), details));
            return std::make_unique<PredicateWorkspace>(ConstructionKey{}, std::move(*storage));
        });
}
core::Result<PredicateLease> PredicateWorkspace::acquire(core::FailureDetails details) noexcept {
    if (busy_.test_and_set(std::memory_order_acquire))
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation, details});
    return PredicateLease{*this};
}
PredicateLease::PredicateLease(PredicateLease&& other) noexcept
    : owner_(std::exchange(other.owner_, nullptr)) {}
PredicateLease& PredicateLease::operator=(PredicateLease&& other) noexcept {
    if (this != &other) {
        if (owner_ != nullptr)
            owner_->busy_.clear(std::memory_order_release);
        owner_ = std::exchange(other.owner_, nullptr);
    }
    return *this;
}
PredicateLease::~PredicateLease() {
    if (owner_ != nullptr)
        owner_->busy_.clear(std::memory_order_release);
}
const TriangleOutcome* ExactTriangleResult::outcome() const& noexcept {
    return outcome_.size() == 1 ? outcome_.values().data() : nullptr;
}
PredicateMemoryPlan predicate_memory_plan() noexcept {
    return PredicateMemoryPlan{28672};
}
PredicateMemoryPlan ExactTriangleResult::memory_plan() const noexcept {
    return predicate_memory_plan();
}
namespace {
bool coordinate_bound(const ExactInteger& value) {
    if (value == 0)
        return true;
    const ExactInteger magnitude = value < 0 ? -value : value;
    return boost::multiprecision::msb(magnitude) < 3173;
}
} // namespace
core::Result<std::reference_wrapper<ExactPoints>>
ExactPredicateAccess::coordinates(PredicateWorkspace& workspace, const PredicateLease& lease,
                                  core::FailureDetails details) noexcept {
    if (lease.owner_ != &workspace)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation, details});
    return std::ref(workspace.scratch_.values()[0].points);
}
core::Result<bool> ExactPredicateAccess::nondegenerate(PredicateWorkspace& workspace,
                                                       const PredicateLease& lease,
                                                       core::FailureDetails details) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::PrecisionUnresolved, details, [&]() -> core::Result<bool> {
            auto points = coordinates(workspace, lease, details);
            if (!points)
                return std::unexpected(points.error());
            for (std::size_t i = 0; i < 3; ++i)
                for (const auto& value : points->get()[i])
                    if (!coordinate_bound(value))
                        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved,
                                                             operation, details});
            auto& scratch = workspace.scratch_.values()[0];
            scratch.e = subtract(scratch.points[1], scratch.points[0]);
            scratch.f = subtract(scratch.points[2], scratch.points[0]);
            scratch.normal = cross(scratch.e, scratch.f);
            return dot(scratch.normal, scratch.normal) != 0;
        });
}
core::Result<ExactTriangleResult>
ExactPredicateAccess::finish(PredicateWorkspace& workspace, const PredicateLease& lease,
                             core::PlannedBuffer<TriangleOutcome> owned, PredicateLimits limits,
                             std::stop_token cancellation, core::FailureDetails details,
                             std::uint64_t steps) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::PrecisionUnresolved, details,
        [&]() -> core::Result<ExactTriangleResult> {
            auto points = coordinates(workspace, lease, details);
            if (!points)
                return std::unexpected(points.error());
            // Only fixed-q factory values enter here. Check the proved magnitude
            // bound again before any product; raw arbitrary integers cannot bypass it.
            for (const auto& point : points->get())
                for (const auto& value : point)
                    if (!coordinate_bound(value))
                        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved,
                                                             operation, details});
            auto& scratch = workspace.scratch_.values()[0];
            const auto checkpoint = [&] {
                if (cancellation.stop_requested())
                    return core::Result<void>{std::unexpected(
                        core::Failure{core::FailureCode::Cancelled, operation, details})};
                if (steps == limits.max_steps)
                    return core::Result<void>{std::unexpected(
                        core::Failure{core::FailureCode::ResourceExceeded, operation, details})};
                ++steps;
                return core::Result<void>{};
            };
            if (auto check = checkpoint(); !check)
                return std::unexpected(check.error());
            scratch.e = subtract(scratch.points[1], scratch.points[0]);
            scratch.f = subtract(scratch.points[2], scratch.points[0]);
            scratch.r = subtract(scratch.points[3], scratch.points[0]);
            scratch.normal = cross(scratch.e, scratch.f);
            const auto jacobian_squared = dot(scratch.normal, scratch.normal);
            if (jacobian_squared == 0)
                return ExactTriangleResult{std::move(owned), limits.external_memory, steps,
                                           details};
            if (auto check = checkpoint(); !check)
                return std::unexpected(check.error());
            const auto ee = dot(scratch.e, scratch.e), ff = dot(scratch.f, scratch.f);
            const auto ef = dot(scratch.e, scratch.f), er = dot(scratch.e, scratch.r);
            const auto fr = dot(scratch.f, scratch.r);
            const ExactInteger n1 = ff * er - ef * fr;
            const ExactInteger n2 = ee * fr - ef * er;
            const ExactInteger n0 = jacobian_squared - n1 - n2;
            const auto projection = projection_feature(n0, n1, n2);
            const auto height = dot(scratch.normal, scratch.r);
            const auto side = height < 0   ? PlaneSide::Below
                              : height > 0 ? PlaneSide::Above
                                           : PlaneSide::On;
            for (std::size_t i = 0; i < 3; ++i) {
                if (auto check = checkpoint(); !check)
                    return std::unexpected(check.error());
                const auto difference = subtract(scratch.points[3], scratch.points[i]);
                scratch.candidates[i].emplace(dot(difference, difference), ExactInteger{1});
            }
            ExactInteger longest = 0;
            for (std::size_t i = 0; i < 3; ++i) {
                if (auto check = checkpoint(); !check)
                    return std::unexpected(check.error());
                const auto edge = subtract(scratch.points[(i + 1) % 3], scratch.points[i]);
                const auto difference = subtract(scratch.points[3], scratch.points[i]);
                const auto length_squared = dot(edge, edge), along = dot(edge, difference);
                if (length_squared > longest)
                    longest = length_squared;
                auto& candidate = scratch.candidates[3 + i];
                if (along <= 0) {
                    const auto& vertex = scratch.candidates[i];
                    if (!vertex)
                        throw std::logic_error("Missing vertex candidate");
                    candidate.emplace(vertex->numerator, ExactInteger{1});
                } else if (along >= length_squared) {
                    const auto& vertex = scratch.candidates[(i + 1) % 3];
                    if (!vertex)
                        throw std::logic_error("Missing vertex candidate");
                    candidate.emplace(vertex->numerator, ExactInteger{1});
                } else {
                    candidate.emplace(dot(difference, difference) * length_squared - along * along,
                                      length_squared);
                }
            }
            if (auto check = checkpoint(); !check)
                return std::unexpected(check.error());
            auto& face = scratch.candidates[6];
            face.reset();
            if (projection != ProjectionFeature::Outside)
                face.emplace(height * height, jacobian_squared);
            std::size_t selected = 0;
            std::uint8_t mask = 0;
            for (std::size_t i = 0; i < 7; ++i) {
                const auto& a = scratch.candidates[i];
                if (!a)
                    continue;
                if (auto check = checkpoint(); !check)
                    return std::unexpected(check.error());
                const auto& b = scratch.candidates[selected];
                if (!b)
                    throw std::logic_error("Missing selected candidate");
                const ExactInteger difference =
                    a->numerator * b->denominator - b->numerator * a->denominator;
                if (difference < 0 || mask == 0) {
                    selected = i;
                    mask = static_cast<std::uint8_t>(1U << i);
                } else if (difference == 0) {
                    mask = static_cast<std::uint8_t>(mask | static_cast<std::uint8_t>(1U << i));
                }
            }
            if (auto check = checkpoint(); !check)
                return std::unexpected(check.error());
            const auto& best = scratch.candidates[selected];
            if (!best)
                throw std::logic_error("Missing minimum distance");
            const ExactInteger comparison = best->numerator - 16 * longest * best->denominator;
            const auto far = comparison > 0    ? FarRelation::Far
                             : comparison == 0 ? FarRelation::Boundary
                                               : FarRelation::Near;
            if (best->numerator < 0 || best->denominator <= 0)
                throw std::logic_error("Invalid exact distance candidate");
            auto& outcome = owned.values()[0];
            outcome.emplace<NondegenerateTriangle>(
                side, projection,
                ExactSquaredDistance{best->numerator,
                                     best->numerator == 0 ? ExactInteger{1} : best->denominator,
                                     squared_unit},
                longest, mask, far);
            return ExactTriangleResult{std::move(owned), limits.external_memory, steps, details};
        });
}
core::Result<ExactTriangleResult>
classify_triangle(CoordinateView coordinates, PredicateWorkspace& workspace,
                  core::ResourcePlan& result_plan, PredicateLimits limits,
                  std::stop_token cancellation, core::FailureDetails details) noexcept {
    auto lease = workspace.acquire(details);
    if (!lease)
        return std::unexpected(lease.error());
    return core::guard_entry(
        operation, core::FailureCode::PrecisionUnresolved, details,
        [&]() -> core::Result<ExactTriangleResult> {
            const auto fail = [&](core::FailureCode code) -> core::Result<ExactTriangleResult> {
                return std::unexpected(core::Failure{code, operation, details});
            };
            std::span<const double> high, low;
            if (const auto* single = std::get_if<Binary64Coordinates>(&coordinates)) {
                high = single->values;
            } else if (const auto* pair = std::get_if<TwoWordCoordinates>(&coordinates)) {
                high = pair->high;
                low = pair->low;
                if (low.size() != 12)
                    return fail(core::FailureCode::InvalidInput);
            } else {
                return fail(core::FailureCode::InvalidInput);
            }
            if (high.size() != 12)
                return fail(core::FailureCode::InvalidInput);
            // Snapshot all words before allocations or exact arithmetic; no later external read.
            std::array<double, 12> input_high{}, input_low{};
            for (std::size_t i = 0; i < 12; ++i) {
                if (!finite_word(high[i]) || (!low.empty() && !finite_word(low[i])))
                    return fail(core::FailureCode::Nonfinite);
                input_high[i] = high[i];
                if (!low.empty())
                    input_low[i] = low[i];
            }
            if (cancellation.stop_requested())
                return fail(core::FailureCode::Cancelled);
            const auto memory = limits.external_memory.selection();
            if (std::holds_alternative<numerics::linear::RequireProvenBound>(memory))
                return fail(core::FailureCode::Unsupported);
            const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&memory);
            if (estimate == nullptr)
                return fail(core::FailureCode::InvalidInput);
            if (estimate->bytes < predicate_memory_plan().peak_payload_estimate_bytes())
                return fail(core::FailureCode::ResourceExceeded);
            auto owned = result_plan.allocate<TriangleOutcome>(1, operation);
            if (!owned)
                return std::unexpected(complete_allocation_failure(owned.error(), details));
            auto& scratch = workspace.scratch_.values()[0];
            std::uint64_t steps = 0;
            const auto checkpoint = [&] {
                if (cancellation.stop_requested())
                    return core::Result<void>{std::unexpected(
                        core::Failure{core::FailureCode::Cancelled, operation, details})};
                if (steps == limits.max_steps)
                    return core::Result<void>{std::unexpected(
                        core::Failure{core::FailureCode::ResourceExceeded, operation, details})};
                ++steps;
                return core::Result<void>{};
            };
            for (std::size_t i = 0; i < 12; ++i) {
                if (auto check = checkpoint(); !check)
                    return std::unexpected(check.error());
                scratch.points[i / 3][i % 3] = decode(input_high[i]) + decode(input_low[i]);
            }
            return ExactPredicateAccess::finish(workspace, *lease, std::move(*owned), limits,
                                                cancellation, details, steps);
        });
}
} // namespace uvlas::aero::detail
