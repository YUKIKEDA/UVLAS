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
                               core::PlannedBuffer<core::Index> heights,
                               core::PlannedBuffer<core::Index> order, std::uint64_t work) noexcept
    : limits_(limits),
      positions_(std::move(positions)),
      frames_(std::move(frames)),
      heights_(std::move(heights)),
      order_(std::move(order)),
      work_(work) {}

DependencyWalk::DependencyWalk(DependencyWalk&& other) noexcept
    : limits_(std::exchange(other.limits_, std::nullopt)),
      positions_(std::move(other.positions_)),
      frames_(std::move(other.frames_)),
      heights_(std::move(other.heights_)),
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
    if (frames > index_max / 2 || limits.sources > index_max / 2 - frames) {
        return fail(Code::ResourceExceeded);
    }
    // Two source arrays and two stack arrays are value-initialized under the work cap.
    const auto cells = 2 * (limits.sources + frames);
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
    auto heights = allocate(frames);
    if (!heights) {
        return std::unexpected(heights.error());
    }
    auto order = allocate(limits.sources);
    if (!order) {
        return std::unexpected(order.error());
    }
    return DependencyWalk{
        limits, std::move(*positions), std::move(*stack), std::move(*heights), std::move(*order),
        cells};
}

core::Result<void> DependencyWalk::charge(std::uint64_t amount) noexcept {
    if (!limits_ || amount > limits_->work - work_) {
        return fail(Code::ResourceExceeded);
    }
    work_ += amount;
    return {};
}

core::Result<void> DependencyWalk::grow_sources(std::uint64_t admitted, core::ResourcePlan& plan,
                                                const std::stop_token& stop) noexcept {
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    if (!limits_ || admitted < positions_.size()) {
        return fail(Code::InvalidInput);
    }
    if (admitted == positions_.size()) {
        return {};
    }
    const auto frames = std::min(limits_->depth, admitted);
    constexpr auto max_index = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    if (frames > max_index / 4 || admitted > max_index / 4 - frames) {
        return fail(Code::ResourceExceeded);
    }
    const auto changes_frames = frames != frames_.size();
    const auto copied = positions_.size() + finished_ + (changes_frames ? 2 * active_ : 0);
    const auto initialized = 2 * (admitted + (changes_frames ? frames : 0));
    const auto work = initialized + copied;
    const auto bytes = core::checked_byte_count(static_cast<core::Index>(initialized),
                                                sizeof(core::Index), operation);
    if (!bytes) {
        return std::unexpected(bytes.error());
    }
    if (work > limits_->work - work_ || *bytes > plan.remaining_bytes()) {
        return fail(Code::ResourceExceeded);
    }
    const auto grow = [&](const auto& from, std::uint64_t size,
                          std::size_t amount) -> core::Result<core::PlannedBuffer<core::Index>> {
        auto target = plan.allocate<core::Index>(static_cast<core::Index>(size), operation);
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        if (target) {
            std::ranges::copy(from.values().first(amount), target->values().begin());
        }
        return target;
    };
    auto positions = grow(positions_, admitted, positions_.size());
    if (!positions) {
        return std::unexpected(positions.error());
    }
    auto order = grow(order_, admitted, finished_);
    if (!order) {
        return std::unexpected(order.error());
    }
    if (changes_frames) {
        auto stack = grow(frames_, frames, active_);
        if (!stack) {
            return std::unexpected(stack.error());
        }
        auto heights = grow(heights_, frames, active_);
        if (!heights) {
            return std::unexpected(heights.error());
        }
        frames_ = std::move(*stack);
        heights_ = std::move(*heights);
    }
    positions_ = std::move(*positions);
    order_ = std::move(*order);
    limits_->sources = admitted;
    work_ += work;
    return {};
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
    const auto reject = [&](Code code, Reason reason, std::optional<std::uint64_t> cap = {},
                            std::optional<std::uint64_t> required = {}) {
        if (diagnostic) {
            diagnostic->reason = reason;
            diagnostic->edge = edge;
            diagnostic->limit = cap;
            diagnostic->required = required;
            if (!required && cap && *cap != std::numeric_limits<std::uint64_t>::max()) {
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
    const auto height = position < 0 ? -position : core::Index{1};
    // The active prefix and completed height fit individually in the frame capacity,
    // which create() bounds by half the maximum Index.
    const auto required_depth = active_ + static_cast<std::uint64_t>(height);
    if (required_depth > limits_->depth) {
        return reject(Code::ResourceExceeded, Reason::DepthLimit, limits_->depth, required_depth);
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    ++edges_;
    ++work_;
    if (position < 0) {
        remember_height(height);
        return DependencyStep::Shared;
    }
    frames_.values()[active_] = edge.dependency;
    heights_.values()[active_++] = height;
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
    const auto height = heights_.values()[active_ - 1];
    positions_.values()[static_cast<std::size_t>(source)] = -height;
    order_.values()[finished_++] = source;
    --active_;
    remember_height(height);
    ++work_;
    return {};
}

void DependencyWalk::remember_height(core::Index dependency_height) noexcept {
    if (active_ != 0) {
        auto& height = heights_.values()[active_ - 1];
        height = std::max(height, dependency_height + 1);
    }
}

std::optional<DependencyWalkView> DependencyWalk::view() const& noexcept {
    if (!limits_) {
        return std::nullopt;
    }
    return DependencyWalkView{frames_.values().first(active_), order_.values().first(finished_),
                              edges_, work_};
}

} // namespace uvlas::io::detail
