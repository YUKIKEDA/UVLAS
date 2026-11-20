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
    declaration_pointer = {};
}

std::string_view AcquiredRootSnapshot::pointer() const& noexcept {
    return text_of(declaration_pointer);
}

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
    if (input_parent->source != declared_parent->source ||
        input_parent->contract_revision != declared_parent->contract_revision ||
        input_parent->representation_revision != declared_parent->representation_revision) {
        mark(diagnostic, RootSnapshotFailure::StaleParent);
        return failure(core::FailureCode::StaleGeneration);
    }

    const auto source = prepared.source(source_ordinal);
    if (!source || origin->root == nullptr) {
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
    const RootSnapshotLocation location{*declared_parent, declared.source_id, declared.role,
                                        declared.declaration_node};
    const auto pointer = declared.declaration_pointer;
    const auto reject = [&](core::Failure error, RootSnapshotFailure reason) {
        return describe_failure(error, reason, location, pointer, limits.diagnostic_bytes, plan,
                                stop, diagnostic);
    };
    const auto acquire = [&]() -> core::Result<ReferenceFileSnapshot> {
        auto reader = origin->root->open_file(source->path, plan, stop,
                                              diagnostic ? &diagnostic->open : nullptr);
        if (!reader) {
            return reject(reader.error(), RootSnapshotFailure::Open);
        }
        auto snapshot = read_reference_file(
            *reader, declared.declared_digest, limits.file, plan, stop,
            diagnostic ? &diagnostic->read : nullptr,
            {origin->manifest_identity, earlier_sources, limits.maximum_identity_comparisons});
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

} // namespace uvlas::io::detail
