#pragma once

#include "reference_declarations.hpp"
#include "reference_path.hpp"

namespace uvlas::io::detail {

struct RootReferenceLimits {
    std::uint64_t declarations;
    std::uint64_t root_edges;
    std::uint64_t unique_sources;
    std::uint64_t normalized_path_bytes;
    std::size_t path_components;
    std::uint64_t work;
    std::uint64_t diagnostic_bytes;
};

enum class RootReferenceFailure : std::uint8_t {
    InvalidOwner,
    Declarations,
    RootEdges,
    UniqueSources,
    Path,
    Work,
    Storage,
};

struct RootReferenceDiagnosticDetail {
    RootReferenceFailure reason;
    std::optional<SchemaIdentity> parent;
    std::optional<SourceId> source;
    std::optional<std::uint64_t> declaration_node;
    std::optional<ReferencePathDiagnostic> path;
    std::optional<std::uint64_t> required;
    std::optional<std::uint64_t> limit;
};

class RootReferenceDiagnostic {
public:
    [[nodiscard]] const std::optional<RootReferenceDiagnosticDetail>& detail() const& noexcept;
    [[nodiscard]] const std::optional<RootReferenceDiagnosticDetail>& detail() const&& = delete;
    [[nodiscard]] std::optional<std::string_view> pointer() const& noexcept;
    [[nodiscard]] std::optional<std::string_view> pointer() const&& = delete;
    void reset() noexcept;

private:
    friend struct RootReferenceContext;
    std::optional<RootReferenceDiagnosticDetail> detail_;
    core::PlannedBuffer<char> pointer_;
};

struct RootReferenceReport {
    std::uint64_t declarations;
    std::uint64_t unique_sources;
    std::uint64_t normalized_path_bytes;
    std::uint64_t work;
    std::uint64_t cumulative_bytes; // New Plan charges, including the temporary ID index.
};

struct RootReferenceDeclarationView {
    DeclaredReferenceView original;
    core::Index source_ordinal{};
};

struct RootReferenceSourceView {
    core::Index first_declaration;
    const NormalizedReferencePath& path;
};

// Owns checked manifest declarations and the ordered mapping to normalized source paths.
// Views borrow the owner until it moves or is destroyed.
class PreparedRootReferences {
public:
    PreparedRootReferences(PreparedRootReferences&&) noexcept = default;
    PreparedRootReferences& operator=(PreparedRootReferences&&) = delete;
    [[nodiscard]] core::Result<SchemaIdentity> identity() const noexcept;
    [[nodiscard]] core::Result<RootReferenceReport> report() const noexcept;
    [[nodiscard]] core::Result<RootReferenceDeclarationView>
        declaration(core::Index) const& noexcept;
    [[nodiscard]] core::Result<RootReferenceDeclarationView>
        declaration(core::Index) const&& = delete;
    [[nodiscard]] core::Result<RootReferenceSourceView> source(core::Index) const& noexcept;
    [[nodiscard]] core::Result<RootReferenceSourceView> source(core::Index) const&& = delete;

private:
    struct Source {
        core::Index first_declaration;
        NormalizedReferencePath path;
    };

    friend core::Result<PreparedRootReferences>
    prepare_root_references(DeclaredReferences&&, RootReferenceLimits, core::ResourcePlan&,
                            const std::stop_token&, RootReferenceDiagnostic*) noexcept;

    PreparedRootReferences(DeclaredReferences declarations, core::PlannedBuffer<core::Index> edges,
                           core::PlannedBuffer<std::optional<Source>> sources,
                           RootReferenceReport report) noexcept;

    DeclaredReferences declarations_;
    core::PlannedBuffer<core::Index> edges_;
    core::PlannedBuffer<std::optional<Source>> sources_;
    RootReferenceReport report_;
};

// Consumes the declarations only on success. The schema already established that all
// occurrences of a SourceId declare the same original path, digest and role.
// Does not open files or infer dependencies. Successful charges stay on the Plan.
[[nodiscard]] core::Result<PreparedRootReferences>
prepare_root_references(DeclaredReferences&& declarations, RootReferenceLimits limits,
                        core::ResourcePlan& plan, const std::stop_token& stop,
                        RootReferenceDiagnostic* diagnostic) noexcept;

} // namespace uvlas::io::detail
