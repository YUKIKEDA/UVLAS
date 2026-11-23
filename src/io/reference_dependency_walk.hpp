#pragma once

#include "reference_declarations.hpp"

namespace uvlas::io::detail {

struct DependencyWalkLimits {
    std::uint64_t sources;
    std::uint64_t edges;
    std::uint64_t depth;
    std::uint64_t work;
    std::uint64_t diagnostic_bytes;
};

struct DependencyEdge {
    std::optional<core::Index> consumer; // Absent for an edge from the manifest.
    core::Index dependency;
    SourceId declaring_source;
    ReferenceRole role;
    std::uint64_t declaration_node;
};

enum class DependencyStep : std::uint8_t { Visit, Shared };
enum class DependencyWalkFailure : std::uint8_t { InvalidEdge, Cycle, EdgeLimit, DepthLimit, Work };

struct DependencyWalkDiagnostic {
    std::optional<DependencyWalkFailure> reason;
    std::optional<DependencyEdge> edge;
    std::optional<std::uint64_t> required, limit;
    core::PlannedBuffer<core::Index> cycle; // Active suffix followed by the closing target.

    void reset() noexcept;
};

struct DependencyWalkView {
    std::span<const core::Index> active;
    std::span<const core::Index> dependency_first;
    std::uint64_t edges{};
    std::uint64_t work{};
};

// Intermediate state for ordered dependency discovery. The caller feeds verified edges in
// declaration order and finishes the active source only after processing all its dependencies.
// Views borrow this owner. Observe them before the next successful step.
// The caller must hold exclusive use of the walker during a transition.
class DependencyWalk {
public:
    DependencyWalk(DependencyWalk&&) noexcept;
    DependencyWalk& operator=(DependencyWalk&&) = delete;

    [[nodiscard]] static core::Result<DependencyWalk> create(DependencyWalkLimits limits,
                                                             core::ResourcePlan& plan,
                                                             const std::stop_token& stop) noexcept;

    [[nodiscard]] core::Result<DependencyStep> enter(const DependencyEdge& edge,
                                                     core::ResourcePlan& plan,
                                                     const std::stop_token& stop,
                                                     DependencyWalkDiagnostic* diagnostic) noexcept;

    [[nodiscard]] core::Result<void> finish(core::Index source,
                                            const std::stop_token& stop) noexcept;

    [[nodiscard]] std::optional<DependencyWalkView> view() const& noexcept;
    [[nodiscard]] std::optional<DependencyWalkView> view() const&& = delete;

private:
    DependencyWalk(DependencyWalkLimits limits, core::PlannedBuffer<core::Index> positions,
                   core::PlannedBuffer<core::Index> frames, core::PlannedBuffer<core::Index> order,
                   std::uint64_t work) noexcept;

    std::optional<DependencyWalkLimits> limits_;
    // Zero means unseen, positive is the one-based stack position, -1 means finished.
    core::PlannedBuffer<core::Index> positions_;
    core::PlannedBuffer<core::Index> frames_;
    core::PlannedBuffer<core::Index> order_;
    std::size_t active_ = 0;
    std::size_t finished_ = 0;
    std::uint64_t edges_ = 0;
    std::uint64_t work_ = 0;
};

} // namespace uvlas::io::detail
