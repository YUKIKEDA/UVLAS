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
    resolve_leaf_references(const ReferenceInput&, PreparedRootReferences&&, RootInventoryLimits,
                            DependencyWalkLimits, GmshLimits, core::ResourcePlan&,
                            const std::stop_token&, RootSnapshotDiagnostic*,
                            LeafReferenceDiagnostic*) noexcept;

    ReferencedCase(AcquiredRootReferences&&, core::PlannedBuffer<core::Index>&&) noexcept;
    AcquiredRootReferences roots_;
    core::PlannedBuffer<core::Index> dependency_first_;
};

// Admit the direct graph, then authenticate and inspect each distinct source in order.
// A source is accepted only if its role proves that it declares no children. Preserve
// repeated manifest edges. Move the declarations only on success.
[[nodiscard]] core::Result<ReferencedCase>
resolve_leaf_references(const ReferenceInput&, PreparedRootReferences&&, RootInventoryLimits,
                        DependencyWalkLimits, GmshLimits, core::ResourcePlan&,
                        const std::stop_token&, RootSnapshotDiagnostic*,
                        LeafReferenceDiagnostic*) noexcept;

} // namespace uvlas::io::detail
