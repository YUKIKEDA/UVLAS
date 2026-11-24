#pragma once

#include "gmsh_mesh.hpp"

namespace uvlas::io::detail {

enum class LeafReferenceFailure : std::uint8_t { UnsupportedRole, Format, Graph };

struct LeafReferenceDiagnostic {
    std::optional<LeafReferenceFailure> reason;
    std::optional<RootSnapshotLocation> location;
    std::optional<std::size_t> format_offset;
    DependencyWalkDiagnostic walk;
    core::PlannedBuffer<char> declaration_pointer;

    [[nodiscard]] std::string_view pointer() const& noexcept;
    [[nodiscard]] std::string_view pointer() const&& = delete;
    void reset() noexcept;
};

class ReferencedCase {
public:
    ReferencedCase(ReferencedCase&&) noexcept = default;
    ReferencedCase& operator=(ReferencedCase&&) = delete;

    [[nodiscard]] const AcquiredRootReferences& roots() const& noexcept {
        return roots_;
    }
    [[nodiscard]] const AcquiredRootReferences& roots() const&& = delete;
    [[nodiscard]] std::span<const core::Index> dependency_first() const& noexcept {
        return dependency_first_.values();
    }
    [[nodiscard]] std::span<const core::Index> dependency_first() const&& = delete;

private:
    friend core::Result<ReferencedCase>
    close_leaf_reference_inventory(AcquiredRootReferences&&, DependencyWalkLimits, GmshLimits,
                                   core::ResourcePlan&, const std::stop_token&,
                                   LeafReferenceDiagnostic*) noexcept;

    ReferencedCase(AcquiredRootReferences&&, core::PlannedBuffer<core::Index>&&) noexcept;
    AcquiredRootReferences roots_;
    core::PlannedBuffer<core::Index> dependency_first_;
};

// Qualify an already authenticated inventory. A source is finished only when its role
// proves that its accepted format declares no children. Preserve repeated manifest edges.
// Consume the inventory only on success. The result owns bytes, paths, edges and order.
[[nodiscard]] core::Result<ReferencedCase>
close_leaf_reference_inventory(AcquiredRootReferences&&, DependencyWalkLimits, GmshLimits,
                               core::ResourcePlan&, const std::stop_token&,
                               LeafReferenceDiagnostic*) noexcept;

} // namespace uvlas::io::detail
