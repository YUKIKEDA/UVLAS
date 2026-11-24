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

void describe(LeafReferenceDiagnostic& diagnostic, RootSnapshotLocation location,
              std::string_view pointer, LeafReferenceFailure reason, std::uint64_t pointer_limit,
              core::ResourcePlan& plan) noexcept {
    diagnostic.reason = reason;
    diagnostic.location = location;
    if (pointer.size() <= pointer_limit && std::in_range<core::Index>(pointer.size())) {
        auto storage = plan.allocate<char>(static_cast<core::Index>(pointer.size()), operation);
        if (storage) {
            std::ranges::copy(pointer, storage->values().begin());
            diagnostic.declaration_pointer = std::move(*storage);
        }
    }
}

class LeafInspector final : public RootSnapshotInspector {
public:
    LeafInspector(LeafFormatLimits limits, std::uint64_t pointer_limit, core::ResourcePlan& plan,
                  const std::stop_token& stop, LeafReferenceDiagnostic* diagnostic) noexcept
        : limits_(limits),
          pointer_limit_(pointer_limit),
          plan_(plan),
          stop_(stop),
          diagnostic_(diagnostic) {}

    [[nodiscard]] core::Result<void>
    inspect(const AcquiredRootSnapshot& snapshot) noexcept override {
        switch (snapshot.location.role) {
        case ReferenceRole::GmshMesh:
            return inspect_gmsh(snapshot);
        case ReferenceRole::NastranDeck:
            return inspect_nastran(snapshot);
        default:
            return reject({Code::Unsupported, operation}, snapshot,
                          LeafReferenceFailure::UnsupportedRole);
        }
    }

private:
    [[nodiscard]] core::Result<void> inspect_gmsh(const AcquiredRootSnapshot& snapshot) noexcept {
        std::size_t offset = 0;
        const auto proof = prove_gmsh_mesh_leaf(snapshot, limits_.gmsh, plan_, stop_, &offset);
        if (!proof) {
            return reject_format(proof.error(), snapshot, offset);
        }
        limits_.gmsh.work -= proof->work();
        return {};
    }

    [[nodiscard]] core::Result<void>
    inspect_nastran(const AcquiredRootSnapshot& snapshot) noexcept {
        std::optional<NastranBulkDiagnostic> detail;
        const auto work =
            verify_nastran_bulk_leaf(snapshot, limits_.nastran, plan_, stop_, &detail);
        if (!work) {
            return reject_format(work.error(), snapshot,
                                 detail.transform([](const NastranBulkDiagnostic& value) {
                                     return value.position.byte;
                                 }));
        }
        limits_.nastran.work -= *work;
        return {};
    }

    [[nodiscard]] core::Result<void> reject_format(core::Failure error,
                                                   const AcquiredRootSnapshot& file,
                                                   std::optional<std::size_t> offset) noexcept {
        if (diagnostic_ && error.code != Code::ResourceExceeded) {
            diagnostic_->format_offset = offset;
        }
        return reject(error, file, LeafReferenceFailure::Format);
    }

    [[nodiscard]] core::Result<void> reject(core::Failure error, const AcquiredRootSnapshot& file,
                                            LeafReferenceFailure reason) noexcept {
        if (diagnostic_ && error.code != Code::Cancelled) {
            describe(*diagnostic_, file.location, file.pointer(), reason, pointer_limit_, plan_);
        }
        if (stop_.stop_requested() || error.code == Code::Cancelled) {
            if (diagnostic_) {
                diagnostic_->reset();
            }
            return fail(Code::Cancelled);
        }
        return std::unexpected(error);
    }

    LeafFormatLimits limits_;
    std::uint64_t pointer_limit_;
    core::ResourcePlan& plan_;
    const std::stop_token& stop_;
    LeafReferenceDiagnostic* diagnostic_;
};

[[nodiscard]] core::Result<core::PlannedBuffer<core::Index>>
plan_leaf_order(const PreparedRootReferences& prepared, const SchemaIdentity& parent,
                const RootReferenceReport& report, DependencyWalkLimits limits,
                core::ResourcePlan& plan, const std::stop_token& stop,
                LeafReferenceDiagnostic* diagnostic) noexcept {
    if (report.unique_sources > limits.sources || report.unique_sources > limits.work / 2 ||
        !std::in_range<core::Index>(report.declarations) ||
        !std::in_range<core::Index>(report.unique_sources)) {
        return fail(Code::ResourceExceeded);
    }
    // Account for initialization and copying of the owned completion order.
    limits.work -= 2 * report.unique_sources;
    limits.sources = report.unique_sources;
    auto walk = DependencyWalk::create(limits, plan, stop);
    if (!walk) {
        return std::unexpected(walk.error());
    }
    // Preflight the direct-edge topology only. Format inspection authorizes its publication.
    for (core::Index i = 0; i < static_cast<core::Index>(report.declarations); ++i) {
        const auto entry = prepared.declaration(i);
        if (!entry) {
            return std::unexpected(entry.error());
        }
        const auto& declaration = entry->original;
        const DependencyEdge edge{std::nullopt, entry->source_ordinal, parent.source.source,
                                  declaration.role, declaration.declaration_node};
        const auto step = walk->enter(edge, plan, stop, diagnostic ? &diagnostic->walk : nullptr);
        if (!step) {
            if (diagnostic && step.error().code != Code::Cancelled) {
                describe(
                    *diagnostic,
                    {parent, declaration.source_id, declaration.role, declaration.declaration_node},
                    declaration.declaration_pointer, LeafReferenceFailure::Graph,
                    limits.diagnostic_bytes, plan);
            }
            if (stop.stop_requested()) {
                if (diagnostic) {
                    diagnostic->reset();
                }
                return fail(Code::Cancelled);
            }
            return std::unexpected(step.error());
        }
        if (*step == DependencyStep::Visit) {
            if (const auto finished = walk->finish(edge.dependency, stop); !finished) {
                return std::unexpected(finished.error());
            }
        }
    }
    const auto result = walk->view();
    if (!result || !result->active.empty() || result->edges != report.declarations ||
        result->dependency_first.size() != report.unique_sources) {
        return fail(Code::InvalidInput);
    }
    auto order =
        plan.allocate<core::Index>(static_cast<core::Index>(report.unique_sources), operation);
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    if (order) {
        std::ranges::copy(result->dependency_first, order->values().begin());
    }
    return order;
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
resolve_leaf_references(const ReferenceInput& input, PreparedRootReferences&& prepared,
                        RootInventoryLimits intake, DependencyWalkLimits limits,
                        LeafFormatLimits formats, core::ResourcePlan& plan,
                        const std::stop_token& stop, RootSnapshotDiagnostic* acquisition,
                        LeafReferenceDiagnostic* diagnostic) noexcept {
    if (acquisition) {
        acquisition->reset();
    }
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    const auto parent = prepared.identity();
    const auto report = prepared.report();
    const auto view = input.view();
    if (!parent || !report || !view) {
        return fail(Code::InvalidInput);
    }
    const auto actual = view->schema.identity();
    if (!actual) {
        return std::unexpected(actual.error());
    }
    if (actual->source != parent->source ||
        actual->contract_revision != parent->contract_revision ||
        actual->representation_revision != parent->representation_revision) {
        return fail(Code::StaleGeneration);
    }
    auto order = plan_leaf_order(prepared, *parent, *report, limits, plan, stop, diagnostic);
    if (!order) {
        if (diagnostic && order.error().code != Code::Cancelled) {
            diagnostic->reason = LeafReferenceFailure::Graph;
        }
        return std::unexpected(order.error());
    }
    LeafInspector inspector{formats, limits.diagnostic_bytes, plan, stop, diagnostic};
    auto inventory = acquire_root_references(input, std::move(prepared), intake, plan, stop,
                                             acquisition, &inspector);
    if (!inventory) {
        return std::unexpected(inventory.error());
    }
    return ReferencedCase{std::move(*inventory), std::move(*order)};
}

} // namespace uvlas::io::detail
