#include "reference_dependency_walk.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
using Code = core::FailureCode;
using Reason = DependencyWalkFailure;

[[nodiscard]] auto fail(Code code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}

} // namespace

void DependencyWalkDiagnostic::reset() noexcept {
    reason.reset();
    edge.reset();
    required.reset();
    limit.reset();
    cycle = {};
}

DependencyWalk::DependencyWalk(DependencyWalkLimits limits,
                               core::PlannedBuffer<core::Index> positions,
                               core::PlannedBuffer<core::Index> frames,
                               core::PlannedBuffer<core::Index> order, std::uint64_t work) noexcept
    : limits_(limits),
      positions_(std::move(positions)),
      frames_(std::move(frames)),
      order_(std::move(order)),
      work_(work) {}

DependencyWalk::DependencyWalk(DependencyWalk&& other) noexcept
    : limits_(std::exchange(other.limits_, std::nullopt)),
      positions_(std::move(other.positions_)),
      frames_(std::move(other.frames_)),
      order_(std::move(other.order_)),
      active_(other.active_),
      finished_(other.finished_),
      edges_(other.edges_),
      work_(other.work_) {}

core::Result<DependencyWalk> DependencyWalk::create(DependencyWalkLimits limits,
                                                    core::ResourcePlan& plan,
                                                    const std::stop_token& stop) noexcept {
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    const auto frames = std::min(limits.depth, limits.sources);
    constexpr auto index_max = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    if (frames > index_max || limits.sources > (index_max - frames) / 2) {
        return fail(Code::ResourceExceeded);
    }
    // The two source arrays and the explicit stack are value-initialized under the work cap.
    const auto cells = 2 * limits.sources + frames;
    const auto bytes =
        core::checked_byte_count(static_cast<core::Index>(cells), sizeof(core::Index), operation);
    if (!bytes) {
        return std::unexpected(bytes.error());
    }
    if (cells > limits.work || *bytes > plan.remaining_bytes()) {
        return fail(Code::ResourceExceeded);
    }
    const auto allocate =
        [&](std::uint64_t count) -> core::Result<core::PlannedBuffer<core::Index>> {
        auto storage = plan.allocate<core::Index>(static_cast<core::Index>(count), operation);
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        return storage;
    };
    auto positions = allocate(limits.sources);
    if (!positions) {
        return std::unexpected(positions.error());
    }
    auto stack = allocate(frames);
    if (!stack) {
        return std::unexpected(stack.error());
    }
    auto order = allocate(limits.sources);
    if (!order) {
        return std::unexpected(order.error());
    }
    return DependencyWalk{limits, std::move(*positions), std::move(*stack), std::move(*order),
                          cells};
}

core::Result<DependencyStep> DependencyWalk::enter(const DependencyEdge& edge,
                                                   core::ResourcePlan& plan,
                                                   const std::stop_token& stop,
                                                   DependencyWalkDiagnostic* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    const auto reject = [&](Code code, Reason reason, std::optional<std::uint64_t> cap = {}) {
        if (diagnostic) {
            diagnostic->reason = reason;
            diagnostic->edge = edge;
            diagnostic->limit = cap;
            if (cap && *cap != std::numeric_limits<std::uint64_t>::max()) {
                diagnostic->required = *cap + 1;
            }
        }
        return fail(code);
    };
    if (!limits_) {
        return reject(Code::InvalidInput, Reason::InvalidEdge);
    }
    const auto expected =
        active_ == 0 ? std::nullopt : std::optional{frames_.values()[active_ - 1]};
    if (edge.consumer != expected || edge.dependency < 0 ||
        std::cmp_greater_equal(edge.dependency, positions_.size())) {
        return reject(Code::InvalidInput, Reason::InvalidEdge);
    }
    if (edges_ == limits_->edges) {
        return reject(Code::ResourceExceeded, Reason::EdgeLimit, limits_->edges);
    }
    if (work_ == limits_->work) {
        return reject(Code::ResourceExceeded, Reason::Work, limits_->work);
    }

    auto& position = positions_.values()[static_cast<std::size_t>(edge.dependency)];
    if (position > 0) {
        const auto first = static_cast<std::size_t>(position - 1);
        const auto count = active_ - first + 1;
        const auto bytes = core::checked_byte_count(static_cast<core::Index>(count),
                                                    sizeof(core::Index), operation);
        // The closing edge needs one work unit. The failed transition commits no work,
        // but diagnostic copying still has to fit after that unit under the same cap.
        const auto diagnostic_work = limits_->work - work_ - 1;
        if (diagnostic && bytes && *bytes <= limits_->diagnostic_bytes &&
            count <= diagnostic_work) {
            auto cycle = plan.allocate<core::Index>(static_cast<core::Index>(count), operation);
            if (cycle) {
                const auto chain = frames_.values().first(active_).subspan(first);
                // Copy the stack suffix without changing traversal state.
                for (std::size_t i = 0; i < chain.size(); ++i) {
                    if (stop.stop_requested()) {
                        return fail(Code::Cancelled);
                    }
                    cycle->values()[i] = chain[i];
                }
                cycle->values().back() = edge.dependency;
                diagnostic->cycle = std::move(*cycle);
            }
        }
        if (stop.stop_requested()) {
            if (diagnostic) {
                diagnostic->reset();
            }
            return fail(Code::Cancelled);
        }
        return reject(Code::InvalidInput, Reason::Cycle);
    }
    if (position == 0 && active_ == frames_.size()) {
        return reject(Code::ResourceExceeded, Reason::DepthLimit, limits_->depth);
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    ++edges_;
    ++work_;
    if (position < 0) {
        return DependencyStep::Shared;
    }
    frames_.values()[active_++] = edge.dependency;
    position = static_cast<core::Index>(active_);
    return DependencyStep::Visit;
}

core::Result<void> DependencyWalk::finish(core::Index source,
                                          const std::stop_token& stop) noexcept {
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    if (!limits_ || active_ == 0 || frames_.values()[active_ - 1] != source) {
        return fail(Code::InvalidInput);
    }
    if (work_ == limits_->work) {
        return fail(Code::ResourceExceeded);
    }
    // Every active slot was admitted by enter, and a finished source is never pushed again.
    positions_.values()[static_cast<std::size_t>(source)] = -1;
    order_.values()[finished_++] = source;
    --active_;
    ++work_;
    return {};
}

std::optional<DependencyWalkView> DependencyWalk::view() const& noexcept {
    if (!limits_) {
        return std::nullopt;
    }
    return DependencyWalkView{frames_.values().first(active_), order_.values().first(finished_),
                              edges_, work_};
}

} // namespace uvlas::io::detail
