#include "reference_root_snapshot.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;

[[nodiscard]] std::string_view text_of(const core::PlannedBuffer<char>& buffer) noexcept {
    return buffer.size() != 0 ? std::string_view{buffer.values().data(), buffer.size()}
                              : std::string_view{};
}

[[nodiscard]] auto failure(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}

[[nodiscard]] bool same_parent(const SchemaIdentity& left, const SchemaIdentity& right) noexcept {
    return left.source == right.source && left.contract_revision == right.contract_revision &&
           left.representation_revision == right.representation_revision;
}

void mark(RootSnapshotDiagnostic* diagnostic, RootSnapshotFailure reason) noexcept {
    if (diagnostic) {
        diagnostic->reason = reason;
    }
}

[[nodiscard]] bool cancelled(const std::stop_token& stop, RootSnapshotDiagnostic* slot) noexcept {
    if (!stop.stop_requested()) {
        return false;
    }
    if (slot) {
        slot->reset();
    }
    return true;
}

#if defined(_WIN32) || defined(__linux__)

[[nodiscard]] core::Result<core::PlannedBuffer<char>>
copy_pointer(std::string_view pointer, core::ResourcePlan& plan) noexcept {
    if (!std::in_range<core::Index>(pointer.size())) {
        return failure(core::FailureCode::ResourceExceeded);
    }
    auto buffer = plan.allocate<char>(static_cast<core::Index>(pointer.size()), operation);
    if (buffer) {
        std::ranges::copy(pointer, buffer->values().begin());
    }
    return buffer;
}

[[nodiscard]] std::unexpected<core::Failure>
describe_failure(core::Failure error, RootSnapshotFailure reason, RootSnapshotLocation location,
                 std::string_view pointer, std::uint64_t diagnostic_bytes, core::ResourcePlan& plan,
                 const std::stop_token& stop, RootSnapshotDiagnostic* slot) noexcept {
    if (cancelled(stop, slot) || error.code == core::FailureCode::Cancelled) {
        if (slot) {
            slot->reset();
        }
        return failure(core::FailureCode::Cancelled);
    }
    if (slot) {
        slot->reason = reason;
        slot->location = location;
        if (pointer.size() <= diagnostic_bytes) {
            if (auto copy = copy_pointer(pointer, plan)) {
                slot->declaration_pointer = std::move(*copy);
            }
        }
    }
    if (cancelled(stop, slot)) {
        return failure(core::FailureCode::Cancelled);
    }
    return std::unexpected(error);
}

#endif

} // namespace

std::string_view RootSnapshotDiagnostic::pointer() const& noexcept {
    return text_of(declaration_pointer);
}

void RootSnapshotDiagnostic::reset() noexcept {
    reason.reset();
    location.reset();
    open.reset();
    read.reset();
    source_quota.reset();
    declaration_pointer = {};
}

std::string_view AcquiredRootSnapshot::pointer() const& noexcept {
    return text_of(declaration_pointer);
}

namespace {

struct BoundRoot {
    ReferenceInputView origin;
    SchemaIdentity parent;
};

[[nodiscard]] core::Result<BoundRoot> bind_root(const ReferenceInput& input,
                                                const PreparedRootReferences& prepared,
                                                RootSnapshotDiagnostic* diagnostic) noexcept {
    const auto origin = input.view();
    const auto declared_parent = prepared.identity();
    if (!origin || !declared_parent) {
        mark(diagnostic, RootSnapshotFailure::InvalidSelection);
        return failure(core::FailureCode::InvalidInput);
    }
    const auto input_parent = origin->schema.identity();
    if (!input_parent) {
        mark(diagnostic, RootSnapshotFailure::InvalidSelection);
        return failure(core::FailureCode::InvalidInput);
    }
    if (!same_parent(*input_parent, *declared_parent)) {
        mark(diagnostic, RootSnapshotFailure::StaleParent);
        return failure(core::FailureCode::StaleGeneration);
    }
    return BoundRoot{*origin, *declared_parent};
}

[[nodiscard]] core::Result<AcquiredRootSnapshot> acquire_bound_root_snapshot(
    const BoundRoot& bound, const PreparedRootReferences& prepared, core::Index source_ordinal,
    std::span<const ReferenceFileId> earlier_sources, RootSnapshotLimits limits,
    core::ResourcePlan& plan, const std::stop_token& stop,
    RootSnapshotDiagnostic* diagnostic) noexcept {
    const auto source = prepared.source(source_ordinal);
    if (!source || bound.origin.root == nullptr) {
        mark(diagnostic, RootSnapshotFailure::InvalidSelection);
        return failure(core::FailureCode::InvalidInput);
    }
    const auto occurrence = prepared.declaration(source->first_declaration);
    if (!occurrence || occurrence->source_ordinal != source_ordinal) {
        mark(diagnostic, RootSnapshotFailure::InvalidSelection);
        return failure(core::FailureCode::InvalidInput);
    }

#if defined(_WIN32) || defined(__linux__)
    const auto& declared = occurrence->original;
    const RootSnapshotLocation location{bound.parent, declared.source_id, declared.role,
                                        declared.declaration_node};
    const auto pointer = declared.declaration_pointer;
    const auto reject = [&](core::Failure error, RootSnapshotFailure reason) {
        return describe_failure(error, reason, location, pointer, limits.diagnostic_bytes, plan,
                                stop, diagnostic);
    };
    const auto acquire = [&]() -> core::Result<ReferenceFileSnapshot> {
        auto reader = bound.origin.root->open_file(source->path, plan, stop,
                                                   diagnostic ? &diagnostic->open : nullptr);
        if (!reader) {
            return reject(reader.error(), RootSnapshotFailure::Open);
        }
        auto snapshot = read_reference_file(
            *reader, declared.declared_digest, limits.file, plan, stop,
            diagnostic ? &diagnostic->read : nullptr,
            {bound.origin.manifest_identity, earlier_sources, limits.maximum_identity_comparisons});
        if (!snapshot) {
            return reject(snapshot.error(), RootSnapshotFailure::Read);
        }
        return snapshot;
    };
    auto file = acquire();
    if (!file) {
        return std::unexpected(file.error());
    }

    auto owned_pointer = copy_pointer(pointer, plan);
    if (cancelled(stop, diagnostic)) {
        return failure(core::FailureCode::Cancelled);
    }
    if (!owned_pointer) {
        return reject(owned_pointer.error(), RootSnapshotFailure::ProvenanceStorage);
    }
    return AcquiredRootSnapshot{location, std::move(*file), std::move(*owned_pointer)};
#else
    static_cast<void>(earlier_sources);
    static_cast<void>(limits);
    static_cast<void>(plan);
    mark(diagnostic, RootSnapshotFailure::Open);
    return failure(core::FailureCode::Unsupported);
#endif
}

} // namespace

core::Result<AcquiredRootSnapshot>
acquire_root_snapshot(const ReferenceInput& input, const PreparedRootReferences& prepared,
                      core::Index source_ordinal, std::span<const ReferenceFileId> earlier_sources,
                      RootSnapshotLimits limits, core::ResourcePlan& plan,
                      const std::stop_token& stop, RootSnapshotDiagnostic* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return failure(core::FailureCode::Cancelled);
    }
    const auto bound = bind_root(input, prepared, diagnostic);
    if (!bound) {
        return std::unexpected(bound.error());
    }
    return acquire_bound_root_snapshot(*bound, prepared, source_ordinal, earlier_sources, limits,
                                       plan, stop, diagnostic);
}

AcquiredRootReferences::AcquiredRootReferences(
    PreparedRootReferences declarations,
    core::PlannedBuffer<std::optional<AcquiredRootSnapshot>> sources) noexcept
    : declarations_(std::move(declarations)), sources_(std::move(sources)) {}

core::Result<SchemaIdentity> AcquiredRootReferences::identity() const noexcept {
    return declarations_.identity();
}

core::Result<RootReferenceReport> AcquiredRootReferences::report() const noexcept {
    return declarations_.report();
}

core::Result<RootReferenceDeclarationView>
AcquiredRootReferences::declaration(core::Index ordinal) const& noexcept {
    return declarations_.declaration(ordinal);
}

core::Result<RootReferenceSourceView>
AcquiredRootReferences::source(core::Index ordinal) const& noexcept {
    return declarations_.source(ordinal);
}

core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
AcquiredRootReferences::snapshot(core::Index ordinal) const& noexcept {
    if (!std::in_range<std::size_t>(ordinal) ||
        static_cast<std::size_t>(ordinal) >= sources_.size()) {
        return failure(core::FailureCode::InvalidInput);
    }
    const auto& source = sources_.values()[static_cast<std::size_t>(ordinal)];
    if (!source) {
        return failure(core::FailureCode::InvalidInput);
    }
    return std::cref(*source);
}

core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
AcquiredRootReferences::find(const SchemaIdentity& parent, SourceId source, ReferenceRole role,
                             const std::stop_token& stop) const& noexcept {
    using Code = core::FailureCode;
    if (stop.stop_requested()) {
        return failure(Code::Cancelled);
    }
    const auto actual = identity();
    if (!actual) {
        return std::unexpected(actual.error());
    }
    if (!same_parent(parent, *actual)) {
        return failure(Code::StaleGeneration);
    }
    for (const auto& entry : sources_.values()) {
        if (stop.stop_requested()) {
            return failure(Code::Cancelled);
        }
        if (entry && entry->location.source == source && entry->location.role == role) {
            return std::cref(*entry);
        }
    }
    return failure(Code::InvalidInput);
}

core::Result<AcquiredRootReferences>
acquire_root_references(const ReferenceInput& input, PreparedRootReferences&& prepared,
                        RootInventoryLimits limits, core::ResourcePlan& plan,
                        const std::stop_token& stop, RootSnapshotDiagnostic* diagnostic,
                        RootSnapshotInspector* inspector) noexcept {
    using Code = core::FailureCode;
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return failure(Code::Cancelled);
    }
    const auto bound = bind_root(input, prepared, diagnostic);
    if (!bound) {
        return std::unexpected(bound.error());
    }
    const auto report = prepared.report();
    if (!report) {
        return std::unexpected(report.error());
    }
    const auto count = report->unique_sources;
    if (count > limits.maximum_sources || !std::in_range<core::Index>(count)) {
        mark(diagnostic, RootSnapshotFailure::SourceLimit);
        if (diagnostic) {
            diagnostic->source_quota = ReferenceFileQuota{count, limits.maximum_sources};
        }
        return failure(Code::ResourceExceeded);
    }
    const auto extent = static_cast<core::Index>(count);
    auto files = plan.allocate<std::optional<AcquiredRootSnapshot>>(extent, operation);
    if (!files) {
        mark(diagnostic, RootSnapshotFailure::InventoryStorage);
        return std::unexpected(files.error());
    }
    auto previous = plan.allocate<ReferenceFileId>(extent, operation);
    if (!previous) {
        mark(diagnostic, RootSnapshotFailure::InventoryStorage);
        return std::unexpected(previous.error());
    }
    const std::uint64_t manifest_checks = bound->origin.manifest_identity.has_value() ? 1 : 0;
    for (std::size_t i = 0; i < files->size(); ++i) {
        if (cancelled(stop, diagnostic)) {
            return failure(Code::Cancelled);
        }
        const auto& per_file = limits.per_source.file;
        const RootSnapshotLimits admitted{
            {std::min(per_file.maximum_file_bytes, limits.total_payload_bytes),
             std::min(per_file.maximum_read_bytes, limits.total_read_bytes),
             std::min(per_file.maximum_read_calls, limits.total_read_calls),
             std::min(per_file.maximum_hash_bytes, limits.total_hash_bytes)},
            std::min(limits.per_source.maximum_identity_comparisons,
                     limits.total_identity_comparisons),
            limits.per_source.diagnostic_bytes};
        auto acquired = acquire_bound_root_snapshot(*bound, prepared, static_cast<core::Index>(i),
                                                    previous->values().first(i), admitted, plan,
                                                    stop, diagnostic);
        if (!acquired) {
            return std::unexpected(acquired.error());
        }
        if (inspector) {
            if (const auto inspected = inspector->inspect(*acquired); !inspected) {
                return std::unexpected(inspected.error());
            }
        }
        const auto file = acquired->file.view();
        if (!file) {
            return failure(Code::InvalidInput);
        }
        const auto bytes = file->file_info.byte_count;
        limits.total_payload_bytes -= bytes;
        limits.total_read_bytes -= bytes + 1;
        limits.total_read_calls -= file->read_calls;
        limits.total_hash_bytes -= bytes;
        limits.total_identity_comparisons -= i + manifest_checks;
        previous->values()[i] = file->file_info.identity;
        files->values()[i].emplace(std::move(*acquired));
    }
    if (cancelled(stop, diagnostic)) {
        return failure(Code::Cancelled);
    }
    return AcquiredRootReferences{std::move(prepared), std::move(*files)};
}

std::string_view ReferenceWalkStartDiagnostic::pointer() const& noexcept {
    return reason == ReferenceWalkStartFailure::Edge ? text_of(edge_pointer) : snapshot.pointer();
}

void ReferenceWalkStartDiagnostic::reset() noexcept {
    reason.reset();
    sources.reset();
    snapshot.reset();
    walk.reset();
    edge_pointer = {};
}

ReferenceWalkStart::ReferenceWalkStart(SchemaIdentity parent, DependencyWalk walk,
                                       std::optional<PendingRootDependencies> pending) noexcept
    : parent_(parent), walk_(std::move(walk)), pending_(std::move(pending)) {}

std::optional<ReferenceWalkStartView> ReferenceWalkStart::view() const& noexcept {
    if (const auto state = walk_.view()) {
        return ReferenceWalkStartView{parent_, *state, pending_ ? &*pending_ : nullptr};
    }
    return std::nullopt;
}

core::Result<ReferenceWalkStart>
start_reference_walk(const ReferenceInput& input, const PreparedRootReferences& prepared,
                     DependencyWalkLimits walk_limits, RootSnapshotLimits snapshot_limits,
                     core::ResourcePlan& plan, const std::stop_token& stop,
                     ReferenceWalkStartDiagnostic* diagnostic) noexcept {
    using Reason = ReferenceWalkStartFailure;
    using Code = core::FailureCode;
    if (diagnostic) {
        diagnostic->reset();
    }
    const auto reject = [&](const core::Failure& cause, Reason reason) {
        if (stop.stop_requested() || cause.code == Code::Cancelled) {
            if (diagnostic) {
                diagnostic->reset();
            }
            return cause.code == Code::Cancelled ? std::unexpected(cause)
                                                 : failure(Code::Cancelled);
        }
        if (diagnostic) {
            diagnostic->reason = reason;
        }
        return std::unexpected(cause);
    };
    if (stop.stop_requested()) {
        return failure(Code::Cancelled);
    }
    const auto bound = bind_root(input, prepared, diagnostic ? &diagnostic->snapshot : nullptr);
    if (!bound) {
        return reject(bound.error(), Reason::Binding);
    }
    const auto report = prepared.report();
    if (!report) {
        return reject(report.error(), Reason::Binding);
    }
    if (report->unique_sources > walk_limits.sources) {
        if (diagnostic) {
            diagnostic->sources = ReferenceFileQuota{report->unique_sources, walk_limits.sources};
        }
        return reject(core::Failure{Code::ResourceExceeded, operation}, Reason::Sources);
    }
    walk_limits.sources = report->unique_sources;
    auto walk = DependencyWalk::create(walk_limits, plan, stop);
    if (!walk) {
        return reject(walk.error(), Reason::WalkStorage);
    }
    if (report->declarations == 0) {
        if (stop.stop_requested()) {
            return failure(Code::Cancelled);
        }
        return ReferenceWalkStart{bound->parent, std::move(*walk), std::nullopt};
    }

    const auto first = prepared.declaration(0);
    if (!first) {
        return reject(first.error(), Reason::Binding);
    }
    const DependencyEdge edge{std::nullopt, first->source_ordinal, bound->parent.source.source,
                              first->original.role, first->original.declaration_node};
    // The snapshot's same-handle identity, bytes and digest precede admission of a verified edge.
    auto snapshot =
        acquire_bound_root_snapshot(*bound, prepared, edge.dependency, {}, snapshot_limits, plan,
                                    stop, diagnostic ? &diagnostic->snapshot : nullptr);
    if (!snapshot) {
        return reject(snapshot.error(), Reason::Snapshot);
    }
    if (const auto entered =
            walk->enter(edge, plan, stop, diagnostic ? &diagnostic->walk : nullptr);
        !entered) {
        if (diagnostic &&
            snapshot->declaration_pointer.size() <= snapshot_limits.diagnostic_bytes) {
            diagnostic->edge_pointer = std::move(snapshot->declaration_pointer);
        }
        return reject(entered.error(), Reason::Edge);
    }
    if (stop.stop_requested()) {
        return failure(Code::Cancelled);
    }
    return ReferenceWalkStart{bound->parent, std::move(*walk),
                              PendingRootDependencies{edge, std::move(*snapshot)}};
}

} // namespace uvlas::io::detail
