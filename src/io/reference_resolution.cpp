#include "reference_resolution.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
using Code = core::FailureCode;

[[nodiscard]] auto fail(Code code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}

void describe(LeafReferenceDiagnostic& diagnostic, const DeclaredReferenceView& declaration,
              const SchemaIdentity& parent, LeafReferenceFailure reason,
              std::uint64_t pointer_limit, core::ResourcePlan& plan) noexcept {
    diagnostic.reason = reason;
    diagnostic.location = RootSnapshotLocation{parent, declaration.source_id, declaration.role,
                                               declaration.declaration_node};
    const auto pointer = declaration.declaration_pointer;
    if (pointer.size() <= pointer_limit && std::in_range<core::Index>(pointer.size())) {
        auto storage = plan.allocate<char>(static_cast<core::Index>(pointer.size()), operation);
        if (storage) {
            std::ranges::copy(pointer, storage->values().begin());
            diagnostic.declaration_pointer = std::move(*storage);
        }
    }
}

} // namespace

std::string_view LeafReferenceDiagnostic::pointer() const& noexcept {
    const auto bytes = declaration_pointer.values();
    return bytes.empty() ? std::string_view{} : std::string_view{bytes.data(), bytes.size()};
}

void LeafReferenceDiagnostic::reset() noexcept {
    reason.reset();
    location.reset();
    format_offset.reset();
    walk.reset();
    declaration_pointer = {};
}

ReferencedCase::ReferencedCase(AcquiredRootReferences&& roots,
                               core::PlannedBuffer<core::Index>&& order) noexcept
    : roots_(std::move(roots)), dependency_first_(std::move(order)) {}

core::Result<ReferencedCase>
close_leaf_reference_inventory(AcquiredRootReferences&& inventory, DependencyWalkLimits limits,
                               GmshLimits gmsh, core::ResourcePlan& plan,
                               const std::stop_token& stop,
                               LeafReferenceDiagnostic* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    const auto parent = inventory.identity();
    const auto report = inventory.report();
    if (!parent || !report) {
        return fail(Code::InvalidInput);
    }
    if (report->unique_sources > limits.sources || report->unique_sources > limits.work / 2 ||
        !std::in_range<core::Index>(report->declarations) ||
        !std::in_range<core::Index>(report->unique_sources)) {
        return fail(Code::ResourceExceeded);
    }
    // Reserve work for value-initializing and copying the owned completion order.
    limits.work -= 2 * report->unique_sources;
    limits.sources = report->unique_sources;
    auto walk = DependencyWalk::create(limits, plan, stop);
    if (!walk) {
        return std::unexpected(walk.error());
    }
    for (core::Index i = 0; i < static_cast<core::Index>(report->declarations); ++i) {
        const auto entry = inventory.declaration(i);
        if (!entry) {
            return std::unexpected(entry.error());
        }
        const auto& declaration = entry->original;
        const DependencyEdge edge{std::nullopt, entry->source_ordinal, parent->source.source,
                                  declaration.role, declaration.declaration_node};
        const auto reject = [&](core::Failure error,
                                LeafReferenceFailure reason) -> core::Result<ReferencedCase> {
            if (diagnostic && error.code != Code::Cancelled) {
                describe(*diagnostic, declaration, *parent, reason, limits.diagnostic_bytes, plan);
            }
            if (stop.stop_requested() || error.code == Code::Cancelled) {
                if (diagnostic) {
                    diagnostic->reset();
                }
                return fail(Code::Cancelled);
            }
            return std::unexpected(error);
        };
        const auto step = walk->enter(edge, plan, stop, diagnostic ? &diagnostic->walk : nullptr);
        if (!step) {
            return reject(step.error(), LeafReferenceFailure::Graph);
        }
        if (*step == DependencyStep::Shared) {
            continue;
        }
        if (declaration.role != ReferenceRole::GmshMesh) {
            return reject({Code::Unsupported, operation}, LeafReferenceFailure::UnsupportedRole);
        }
        const auto snapshot = inventory.snapshot(entry->source_ordinal);
        if (!snapshot) {
            return reject(snapshot.error(), LeafReferenceFailure::Graph);
        }
        std::size_t offset = 0;
        const auto proof = prove_gmsh_mesh_leaf(snapshot->get(), gmsh, plan, stop, &offset);
        if (!proof) {
            if (diagnostic && proof.error().code != Code::ResourceExceeded) {
                diagnostic->format_offset = offset;
            }
            return reject(proof.error(), LeafReferenceFailure::Format);
        }
        gmsh.work -= proof->work();
        if (const auto finished = walk->finish(edge.dependency, stop); !finished) {
            return reject(finished.error(), LeafReferenceFailure::Graph);
        }
    }
    const auto result = walk->view();
    if (!result || !result->active.empty() || result->edges != report->declarations ||
        result->dependency_first.size() != report->unique_sources) {
        return fail(Code::InvalidInput);
    }
    auto order =
        plan.allocate<core::Index>(static_cast<core::Index>(report->unique_sources), operation);
    if (!order) {
        return std::unexpected(order.error());
    }
    std::ranges::copy(result->dependency_first, order->values().begin());
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    return ReferencedCase{std::move(inventory), std::move(*order)};
}

} // namespace uvlas::io::detail
