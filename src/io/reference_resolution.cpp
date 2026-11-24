#include "reference_resolution.hpp"

#include "utf8.hpp"

#include <algorithm>
#include <limits>
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

class LeafInspector final : public ReferenceDependencyExtractor {
public:
    LeafInspector(LeafFormatLimits limits, std::uint64_t pointer_limit, core::ResourcePlan& plan,
                  const std::stop_token& stop, LeafReferenceDiagnostic* diagnostic) noexcept
        : limits_(limits),
          pointer_limit_(pointer_limit),
          plan_(plan),
          stop_(stop),
          diagnostic_(diagnostic) {}

    [[nodiscard]] core::Result<ExtractedReferenceDependencies>
    extract(const AcquiredRootSnapshot& snapshot) noexcept override {
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
    [[nodiscard]] core::Result<ExtractedReferenceDependencies>
    inspect_gmsh(const AcquiredRootSnapshot& snapshot) noexcept {
        std::size_t offset = 0;
        const auto proof = prove_gmsh_mesh_leaf(snapshot, limits_.gmsh, plan_, stop_, &offset);
        if (!proof) {
            return reject_format(proof.error(), snapshot, offset);
        }
        limits_.gmsh.work -= proof->work();
        return ExtractedReferenceDependencies{};
    }

    [[nodiscard]] core::Result<ExtractedReferenceDependencies>
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
        return ExtractedReferenceDependencies{};
    }

    [[nodiscard]] core::Result<ExtractedReferenceDependencies>
    reject_format(core::Failure error, const AcquiredRootSnapshot& file,
                  std::optional<std::size_t> offset) noexcept {
        if (diagnostic_ && error.code != Code::ResourceExceeded) {
            diagnostic_->format_offset = offset;
        }
        return reject(error, file, LeafReferenceFailure::Format);
    }

    [[nodiscard]] core::Result<ExtractedReferenceDependencies>
    reject(core::Failure error, const AcquiredRootSnapshot& file,
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

constexpr auto index_max = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());

bool same_parent(const SchemaIdentity& a, const SchemaIdentity& b) noexcept {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}

std::string_view text_of(const core::PlannedBuffer<char>& text) noexcept {
    return text.size() == 0 ? std::string_view{}
                            : std::string_view{text.values().data(), text.size()};
}

std::size_t slot_of(std::optional<core::Index> consumer) noexcept {
    return consumer ? static_cast<std::size_t>(*consumer) + 1 : 0;
}

struct Frame {
    ExtractedReferenceDependencies children;
    std::size_t next{};
};

} // namespace

class ReferenceGraphBuilder {
public:
    ReferenceGraphBuilder(const ReferenceInputView& origin, SchemaIdentity parent,
                          const PreparedRootReferences& prepared, ReferenceGraphLimits limits,
                          DependencyWalk walk, ReferenceDependencyExtractor& extractor,
                          core::ResourcePlan& plan, const std::stop_token& stop,
                          RootSnapshotDiagnostic* acquisition,
                          LeafReferenceDiagnostic* diagnostic) noexcept
        : origin_(origin),
          parent_(parent),
          prepared_(prepared),
          limits_(limits),
          walk_(std::move(walk)),
          extractor_(extractor),
          plan_(plan),
          stop_(stop),
          acquisition_(acquisition),
          diagnostic_(diagnostic) {}

    core::Result<void> collect(const RootReferenceReport& report) noexcept {
        if (report.declarations > limits_.walk.edges) {
            const auto declaration =
                prepared_.declaration(static_cast<core::Index>(limits_.walk.edges));
            if (!declaration) {
                return std::unexpected(declaration.error());
            }
            const auto& original = declaration->original;
            const DependencyEdge edge{std::nullopt, declaration->source_ordinal,
                                      parent_.source.source, original.role,
                                      original.declaration_node};
            if (diagnostic_) {
                diagnostic_->walk.reason = DependencyWalkFailure::EdgeLimit;
                diagnostic_->walk.edge = edge;
                diagnostic_->walk.required = limits_.walk.edges + 1;
                diagnostic_->walk.limit = limits_.walk.edges;
            }
            return reject({Code::ResourceExceeded, operation}, edge, original.source_id,
                          original.declaration_pointer);
        }
        if (const auto setup = configure_roots(report); !setup) {
            return setup;
        }
        for (core::Index i = 0; i < static_cast<core::Index>(report.declarations); ++i) {
            const auto declaration = prepared_.declaration(i);
            if (!declaration) {
                return std::unexpected(declaration.error());
            }
            const auto& original = declaration->original;
            const DependencyEdge edge{std::nullopt, declaration->source_ordinal,
                                      parent_.source.source, original.role,
                                      original.declaration_node};
            if (const auto visited = follow(edge, i, original.declaration_pointer, {}); !visited) {
                return visited;
            }
            if (const auto descended = descend(); !descended) {
                return descended;
            }
        }
        return {};
    }

    core::Result<ReferencedCase> publish(PreparedRootReferences&& prepared) noexcept {
        const auto state = walk_.view();
        if (!state || !state->active.empty() || state->edges != edge_count_ ||
            state->dependency_first.size() != static_cast<std::size_t>(vertex_count_)) {
            return fail(Code::InvalidInput);
        }
        if (vertex_count_ == 0) {
            if (stop_.stop_requested()) {
                return fail(Code::Cancelled);
            }
            return ReferencedCase{std::move(prepared), {}, 0, 0, {}, {}, {}};
        }
        const auto n = static_cast<std::uint64_t>(vertex_count_);
        const auto offset_count = n + 2;
        // Initialization, counting/filling, prefix/restore, and the owned completion order.
        for (const auto count : {offset_count, offset_count, offset_count, offset_count,
                                 edge_count_, edge_count_, edge_count_, n, n}) {
            if (const auto work = walk_.charge(count); !work) {
                return std::unexpected(work.error());
            }
        }
        auto offsets =
            plan_.allocate<core::Index>(static_cast<core::Index>(offset_count), operation);
        if (!offsets) {
            return std::unexpected(offsets.error());
        }
        auto ordered = plan_.allocate<Occurrence>(static_cast<core::Index>(edge_count_), operation);
        if (!ordered) {
            return std::unexpected(ordered.error());
        }
        auto order = plan_.allocate<core::Index>(vertex_count_, operation);
        if (!order) {
            return std::unexpected(order.error());
        }
        auto bounds = offsets->values();
        auto edges = edges_.values().first(static_cast<std::size_t>(edge_count_));
        for (const auto& entry : edges) {
            ++bounds[slot_of(entry.edge.consumer) + 1];
        }
        for (std::size_t i = 1; i < bounds.size(); ++i) {
            bounds[i] += bounds[i - 1];
        }
        // Use starts as write cursors, then shift the resulting ends to restore the starts.
        for (auto& entry : edges) {
            auto& next = bounds[slot_of(entry.edge.consumer)];
            ordered->values()[static_cast<std::size_t>(next++)] = std::move(entry);
        }
        std::ranges::move_backward(bounds.first(bounds.size() - 1), bounds.end());
        bounds.front() = 0;
        std::ranges::copy(state->dependency_first, order->values().begin());
        if (stop_.stop_requested()) {
            return fail(Code::Cancelled);
        }
        return ReferencedCase{std::move(prepared), std::move(vertices_), vertex_count_,
                              root_count_,         std::move(*ordered),  std::move(*offsets),
                              std::move(*order)};
    }

private:
    using Vertex = ReferencedCase::Vertex;
    using Occurrence = ReferencedCase::Occurrence;
    core::Result<void> configure_roots(const RootReferenceReport&) noexcept;
    core::Result<core::Index> admit_child(const ChildReferenceDeclaration&,
                                          const ReferenceFileView&, const DependencyEdge&) noexcept;
    core::Result<void> descend() noexcept;
    core::Result<void> follow(const DependencyEdge&, core::Index, std::string_view,
                              std::string_view) noexcept;
    core::Result<void> authenticate(Vertex&, const DependencyEdge&, std::string_view) noexcept;
    core::Result<core::PlannedBuffer<char>> copy_text(std::string_view) noexcept;
    core::Result<std::string_view> original_path(core::Index) const noexcept;
    core::Result<std::reference_wrapper<const NormalizedReferencePath>>
        normalized_path(core::Index) const noexcept;
    std::unexpected<core::Failure> reject(core::Failure, const DependencyEdge&, SourceId,
                                          std::string_view) noexcept;

    template <class T>
    core::Result<void> reserve(core::PlannedBuffer<T>& storage, std::uint64_t count,
                               std::uint64_t used, std::uint64_t cap) noexcept {
        if (stop_.stop_requested()) {
            return fail(Code::Cancelled);
        }
        if (count > cap) {
            return fail(Code::ResourceExceeded);
        }
        if (count <= storage.size()) {
            return {};
        }
        const auto next = count + std::min(static_cast<std::uint64_t>(storage.size()), cap - count);
        if (const auto work = walk_.charge(next + used); !work) {
            return work;
        }
        auto fresh = plan_.allocate<T>(static_cast<core::Index>(next), operation);
        if (stop_.stop_requested()) {
            return fail(Code::Cancelled);
        }
        if (!fresh) {
            return std::unexpected(fresh.error());
        }
        std::ranges::move(storage.values().first(static_cast<std::size_t>(used)),
                          fresh->values().begin());
        storage = std::move(*fresh);
        return {};
    }

    const ReferenceInputView& origin_;
    SchemaIdentity parent_;
    const PreparedRootReferences& prepared_;
    ReferenceGraphLimits limits_;
    DependencyWalk walk_;
    ReferenceDependencyExtractor& extractor_;
    core::ResourcePlan& plan_;
    const std::stop_token& stop_;
    RootSnapshotDiagnostic* acquisition_;
    LeafReferenceDiagnostic* diagnostic_;
    core::PlannedBuffer<Vertex> vertices_;
    core::PlannedBuffer<ReferenceFileId> identities_;
    core::PlannedBuffer<Occurrence> edges_;
    core::PlannedBuffer<Frame> frames_;
    core::Index vertex_count_{};
    core::Index root_count_{};
    std::size_t acquired_{};
    std::uint64_t edge_count_{};
    std::uint64_t path_bytes_{};
};

core::Result<void>
ReferenceGraphBuilder::configure_roots(const RootReferenceReport& report) noexcept {
    const auto count = report.unique_sources;
    const auto cap = std::min(limits_.walk.sources, limits_.inventory.maximum_sources);
    if (count > cap || report.normalized_path_bytes > limits_.normalized_path_bytes) {
        return fail(Code::ResourceExceeded);
    }
    if (const auto stored = reserve(vertices_, count, 0, cap); !stored) {
        return stored;
    }
    if (const auto stored = reserve(identities_, count, 0, cap); !stored) {
        return stored;
    }
    if (const auto work = walk_.charge(count); !work) {
        return work;
    }
    for (core::Index i = 0; i < static_cast<core::Index>(count); ++i) {
        const auto vertex = prepared_.source(i);
        if (!vertex) {
            return std::unexpected(vertex.error());
        }
        const auto declaration = prepared_.declaration(vertex->first_declaration);
        if (!declaration) {
            return std::unexpected(declaration.error());
        }
        const auto& original = declaration->original;
        auto& entry = vertices_.values()[static_cast<std::size_t>(i)];
        entry.id = original.source_id;
        entry.role = original.role;
        entry.digest = original.declared_digest;
    }
    root_count_ = static_cast<core::Index>(count);
    vertex_count_ = root_count_;
    path_bytes_ = report.normalized_path_bytes;
    return reserve(edges_, report.declarations, 0, limits_.walk.edges);
}

core::Result<core::Index>
ReferenceGraphBuilder::admit_child(const ChildReferenceDeclaration& child,
                                   const ReferenceFileView& parent_bytes,
                                   const DependencyEdge& location) noexcept {
    const auto reject_child = [&](core::Failure error) {
        return reject(error, location, child.source, {});
    };
    if (child.original_path.size() > limits_.child_path.input_bytes ||
        child.root_relative_path.size() > limits_.child_path.input_bytes) {
        return reject_child({Code::ResourceExceeded, operation});
    }
    // Bound the original UTF-8 check and the forward/backward normalization passes.
    for (const auto amount : {child.original_path.size(), child.original_path.size(),
                              child.root_relative_path.size(), child.root_relative_path.size(),
                              child.root_relative_path.size(), child.root_relative_path.size(),
                              child.root_relative_path.size(), child.root_relative_path.size()}) {
        if (const auto work = walk_.charge(amount); !work) {
            return reject_child(work.error());
        }
    }
    if (child.declaration_byte >= parent_bytes.bytes.size() || !reference_role_name(child.role) ||
        child.original_path.empty() || child.original_path.find('\0') != std::string_view::npos ||
        !valid_utf8(child.original_path) || !valid_utf8(child.root_relative_path)) {
        return reject_child({Code::InvalidInput, operation});
    }
    std::optional<core::Index> known;
    for (core::Index i = 0; i < vertex_count_; ++i) {
        if (const auto work = walk_.charge(1); !work) {
            return reject_child(work.error());
        }
        if (vertices_.values()[static_cast<std::size_t>(i)].id == child.source) {
            known = i;
            break;
        }
    }
    auto path_limits = limits_.child_path;
    if (!known) {
        const auto remaining = limits_.normalized_path_bytes - path_bytes_;
        if (remaining < path_limits.normalized_bytes) {
            path_limits.normalized_bytes = static_cast<std::size_t>(remaining);
        }
    }
    std::optional<ReferencePathDiagnostic> path_diagnostic;
    auto path = normalize_reference_path(child.root_relative_path, path_limits, plan_, stop_,
                                         diagnostic_ ? &path_diagnostic : nullptr);
    if (!path) {
        if (diagnostic_) {
            diagnostic_->path = path_diagnostic;
        }
        return reject_child(path.error());
    }
    if (known) {
        const auto ordinal = *known;
        const auto& vertex = vertices_.values()[static_cast<std::size_t>(ordinal)];
        const auto original = original_path(ordinal);
        const auto normalized = normalized_path(ordinal);
        if (!original || !normalized) {
            return reject_child({Code::InvalidInput, operation});
        }
        for (const auto amount : {std::min(original->size(), child.original_path.size()),
                                  std::min(normalized->get().text().size(), path->text().size())}) {
            if (const auto work = walk_.charge(amount); !work) {
                return reject_child(work.error());
            }
        }
        if (vertex.role != child.role || vertex.digest != child.digest ||
            *original != child.original_path || normalized->get().text() != path->text()) {
            return reject_child({Code::InvalidInput, operation});
        }
        return ordinal;
    }
    const auto cap = std::min(limits_.walk.sources, limits_.inventory.maximum_sources);
    const auto count = static_cast<std::uint64_t>(vertex_count_);
    if (const auto stored = reserve(vertices_, count + 1, count, cap); !stored) {
        return reject_child(stored.error());
    }
    if (const auto stored = reserve(identities_, count + 1, acquired_, cap); !stored) {
        return reject_child(stored.error());
    }
    if (const auto grown = walk_.grow_sources(vertices_.size(), plan_, stop_); !grown) {
        return reject_child(grown.error());
    }
    auto original = copy_text(child.original_path);
    if (!original) {
        return reject_child(original.error());
    }
    path_bytes_ += path->text().size();
    auto& vertex = vertices_.values()[static_cast<std::size_t>(count)];
    vertex.id = child.source;
    vertex.role = child.role;
    vertex.digest = child.digest;
    vertex.child_path.emplace(std::move(*path));
    vertex.original_child_path = std::move(*original);
    return vertex_count_++;
}

core::Result<void> ReferenceGraphBuilder::descend() noexcept {
    while (true) {
        if (stop_.stop_requested()) {
            return fail(Code::Cancelled);
        }
        const auto state = walk_.view();
        if (!state || state->active.empty()) {
            return {};
        }
        const auto parent = state->active.back();
        auto& frame = frames_.values()[state->active.size() - 1];
        if (frame.next == frame.children.declarations.size()) {
            frame = {};
            if (const auto finished = walk_.finish(parent, stop_); !finished) {
                return finished;
            }
            continue;
        }
        // Copy the descriptor before following it. Descending can grow the frame array.
        const auto child = frame.children.declarations.values()[frame.next++];
        const auto& vertex = vertices_.values()[static_cast<std::size_t>(parent)];
        const auto bytes = vertex.snapshot->file.view();
        if (!bytes) {
            return fail(Code::InvalidInput);
        }
        DependencyEdge edge{parent, -1, vertex.id, child.role, child.declaration_byte};
        const auto dependency = admit_child(child, *bytes, edge);
        if (!dependency) {
            return std::unexpected(dependency.error());
        }
        edge.dependency = *dependency;
        if (const auto visited = follow(edge, -1, {}, child.original_path); !visited) {
            return visited;
        }
    }
}

core::Result<void> ReferenceGraphBuilder::authenticate(Vertex& vertex, const DependencyEdge& edge,
                                                       std::string_view pointer) noexcept {
    const auto path = normalized_path(edge.dependency);
    if (!path) {
        return reject(path.error(), edge, vertex.id, pointer);
    }
    auto snapshot = acquire_reference_snapshot(
        origin_, path->get(), {parent_, vertex.id, vertex.role, edge.declaration_node},
        vertex.digest, pointer, identities_.values().first(acquired_),
        limits_.inventory.next_file(), plan_, stop_, acquisition_);
    if (!snapshot) {
        if (diagnostic_ && edge.consumer && snapshot.error().code != Code::Cancelled) {
            diagnostic_->declaring_source = edge.declaring_source;
        }
        return std::unexpected(snapshot.error());
    }
    const auto file = snapshot->file.view();
    if (!file) {
        return reject({Code::InvalidInput, operation}, edge, vertex.id, pointer);
    }
    limits_.inventory.charge(*file, acquired_ + origin_.manifest_identity.has_value());
    identities_.values()[acquired_++] = file->file_info.identity;
    vertex.snapshot.emplace(std::move(*snapshot));
    return {};
}

core::Result<void> ReferenceGraphBuilder::follow(const DependencyEdge& edge,
                                                 core::Index root_declaration,
                                                 std::string_view pointer,
                                                 std::string_view child_original) noexcept {
    auto& vertex = vertices_.values()[static_cast<std::size_t>(edge.dependency)];
    const auto reject_edge = [&](core::Failure error) {
        return reject(error, edge, vertex.id, pointer);
    };
    // Authenticate an unseen target before the DFS transition. Active and completed targets
    // already have their checked bytes and can be classified without opening another handle.
    if (!vertex.snapshot) {
        if (const auto ready = authenticate(vertex, edge, pointer); !ready) {
            return ready;
        }
    }
    const auto step = walk_.enter(edge, plan_, stop_, diagnostic_ ? &diagnostic_->walk : nullptr);
    if (!step) {
        return reject_edge(step.error());
    }
    if (const auto stored = reserve(edges_, edge_count_ + 1, edge_count_, limits_.walk.edges);
        !stored) {
        return reject_edge(stored.error());
    }
    auto& occurrence = edges_.values()[static_cast<std::size_t>(edge_count_++)];
    occurrence.edge = edge;
    occurrence.root_declaration = root_declaration;
    if (edge.consumer) {
        auto original = copy_text(child_original);
        if (!original) {
            return reject_edge(original.error());
        }
        occurrence.original_child_path = std::move(*original);
    }
    if (*step == DependencyStep::Shared) {
        return {};
    }
    auto extracted = extractor_.extract(*vertex.snapshot);
    if (!extracted) {
        return std::unexpected(extracted.error());
    }
    const auto state = walk_.view();
    if (!state) {
        return reject_edge({Code::InvalidInput, operation});
    }
    const auto depth = state->active.size();
    const auto maximum =
        std::min({limits_.walk.sources, limits_.inventory.maximum_sources, limits_.walk.depth});
    if (const auto stored = reserve(frames_, depth, depth - 1, maximum); !stored) {
        return reject_edge(stored.error());
    }
    frames_.values()[depth - 1] = {std::move(*extracted), 0};
    return {};
}

core::Result<core::PlannedBuffer<char>>
ReferenceGraphBuilder::copy_text(std::string_view text) noexcept {
    if (!std::in_range<core::Index>(text.size())) {
        return fail(Code::ResourceExceeded);
    }
    if (const auto work = walk_.charge(text.size() * std::uint64_t{2}); !work) {
        return std::unexpected(work.error());
    }
    auto result = plan_.allocate<char>(static_cast<core::Index>(text.size()), operation);
    if (stop_.stop_requested()) {
        return fail(Code::Cancelled);
    }
    if (result) {
        std::ranges::copy(text, result->values().begin());
    }
    return result;
}

core::Result<std::string_view>
ReferenceGraphBuilder::original_path(core::Index ordinal) const noexcept {
    if (ordinal < root_count_) {
        const auto source = prepared_.source(ordinal);
        if (!source) {
            return std::unexpected(source.error());
        }
        const auto declaration = prepared_.declaration(source->first_declaration);
        if (!declaration) {
            return std::unexpected(declaration.error());
        }
        return declaration->original.original_path;
    }
    return text_of(vertices_.values()[static_cast<std::size_t>(ordinal)].original_child_path);
}

core::Result<std::reference_wrapper<const NormalizedReferencePath>>
ReferenceGraphBuilder::normalized_path(core::Index ordinal) const noexcept {
    if (ordinal < root_count_) {
        const auto source = prepared_.source(ordinal);
        if (!source) {
            return std::unexpected(source.error());
        }
        return std::cref(source->path);
    }
    const auto& path = vertices_.values()[static_cast<std::size_t>(ordinal)].child_path;
    if (!path) {
        return fail(Code::InvalidInput);
    }
    return std::cref(*path);
}

std::unexpected<core::Failure> ReferenceGraphBuilder::reject(core::Failure error,
                                                             const DependencyEdge& edge,
                                                             SourceId target,
                                                             std::string_view pointer) noexcept {
    if (diagnostic_ && error.code != Code::Cancelled) {
        const auto declaring = edge.consumer ? std::optional{edge.declaring_source} : std::nullopt;
        describe(*diagnostic_, {parent_, target, edge.role, edge.declaration_node, declaring},
                 pointer, LeafReferenceFailure::Graph, limits_.walk.diagnostic_bytes, plan_);
        diagnostic_->declaring_source = edge.declaring_source;
    }
    if (stop_.stop_requested() || error.code == Code::Cancelled) {
        if (diagnostic_) {
            diagnostic_->reset();
        }
        return fail(Code::Cancelled);
    }
    return std::unexpected(error);
}

std::string_view LeafReferenceDiagnostic::pointer() const& noexcept {
    return text_of(declaration_pointer);
}

void LeafReferenceDiagnostic::reset() noexcept {
    reason.reset();
    location.reset();
    format_offset.reset();
    declaring_source.reset();
    path.reset();
    walk.reset();
    declaration_pointer = {};
}

ReferencedCase::ReferencedCase(PreparedRootReferences&& prepared,
                               core::PlannedBuffer<Vertex>&& vertices, core::Index vertex_count,
                               core::Index root_count,
                               core::PlannedBuffer<Occurrence>&& occurrences,
                               core::PlannedBuffer<core::Index>&& offsets,
                               core::PlannedBuffer<core::Index>&& order) noexcept
    : prepared_(std::move(prepared)),
      vertices_(std::move(vertices)),
      vertex_count_(vertex_count),
      root_count_(root_count),
      occurrences_(std::move(occurrences)),
      offsets_(std::move(offsets)),
      dependency_first_(std::move(order)) {}

core::Result<SchemaIdentity> ReferencedCase::identity() const noexcept {
    return prepared_.identity();
}

std::span<const core::Index> ReferencedCase::dependency_first() const& noexcept {
    return dependency_first_.values();
}

core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
ReferencedCase::snapshot(core::Index source) const& noexcept {
    if (source < 0 || std::cmp_greater_equal(source, vertices_.size()) || source >= vertex_count_) {
        return fail(Code::InvalidInput);
    }
    const auto& verified = vertices_.values()[static_cast<std::size_t>(source)].snapshot;
    if (!verified) {
        return fail(Code::InvalidInput);
    }
    return std::cref(*verified);
}

core::Result<std::reference_wrapper<const AcquiredRootSnapshot>>
ReferencedCase::find(const SchemaIdentity& parent, SourceId source, ReferenceRole role,
                     const std::stop_token& stop) const& noexcept {
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    const auto actual = identity();
    if (!actual) {
        return std::unexpected(actual.error());
    }
    if (!same_parent(*actual, parent)) {
        return fail(Code::StaleGeneration);
    }
    for (core::Index i = 0; i < vertex_count_; ++i) {
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        const auto& vertex = vertices_.values()[static_cast<std::size_t>(i)];
        if (vertex.id == source && vertex.role == role) {
            return snapshot(i);
        }
    }
    return fail(Code::InvalidInput);
}

core::Result<std::string_view> ReferencedCase::normalized_path(core::Index source) const& noexcept {
    if (source < 0 || std::cmp_greater_equal(source, vertices_.size()) || source >= vertex_count_) {
        return fail(Code::InvalidInput);
    }
    if (source < root_count_) {
        const auto root = prepared_.source(source);
        if (!root) {
            return std::unexpected(root.error());
        }
        return root->path.text();
    }
    const auto& path = vertices_.values()[static_cast<std::size_t>(source)].child_path;
    return path ? core::Result<std::string_view>{path->text()} : fail(Code::InvalidInput);
}

core::Result<std::span<const ReferencedCase::Occurrence>>
ReferencedCase::occurrences(std::optional<core::Index> consumer) const noexcept {
    if (!identity() || (consumer && (*consumer < 0 || *consumer >= vertex_count_))) {
        return fail(Code::InvalidInput);
    }
    if (offsets_.size() == 0) {
        return std::span<const Occurrence>{};
    }
    const auto slot = slot_of(consumer);
    const auto first = static_cast<std::size_t>(offsets_.values()[slot]);
    const auto end = static_cast<std::size_t>(offsets_.values()[slot + 1]);
    return occurrences_.values().subspan(first, end - first);
}

core::Result<core::Index>
ReferencedCase::dependency_count(std::optional<core::Index> consumer) const noexcept {
    const auto range = occurrences(consumer);
    if (!range) {
        return std::unexpected(range.error());
    }
    return static_cast<core::Index>(range->size());
}

core::Result<ReferenceGraphEdgeView>
ReferencedCase::dependency(std::optional<core::Index> consumer,
                           core::Index ordinal) const& noexcept {
    const auto range = occurrences(consumer);
    if (!range || ordinal < 0 || std::cmp_greater_equal(ordinal, range->size())) {
        return fail(Code::InvalidInput);
    }
    const auto& entry = (*range)[static_cast<std::size_t>(ordinal)];
    if (entry.root_declaration >= 0) {
        const auto declared = prepared_.declaration(entry.root_declaration);
        if (!declared) {
            return std::unexpected(declared.error());
        }
        return ReferenceGraphEdgeView{entry.edge, declared->original.original_path,
                                      declared->original.declaration_pointer};
    }
    return ReferenceGraphEdgeView{entry.edge, text_of(entry.original_child_path), {}};
}

core::Result<ReferencedCase> resolve_reference_graph(
    const ReferenceInput& input, PreparedRootReferences&& prepared, ReferenceGraphLimits limits,
    ReferenceDependencyExtractor& extractor, core::ResourcePlan& plan, const std::stop_token& stop,
    RootSnapshotDiagnostic* acquisition, LeafReferenceDiagnostic* diagnostic) noexcept {
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
    const auto origin = input.view();
    if (!parent || !report || !origin) {
        return fail(Code::InvalidInput);
    }
    const auto actual = origin->schema.identity();
    if (!actual) {
        return std::unexpected(actual.error());
    }
    if (!same_parent(*actual, *parent)) {
        return fail(Code::StaleGeneration);
    }
    const auto graph_failure = [&](core::Failure error) {
        if (stop.stop_requested() || error.code == Code::Cancelled) {
            if (diagnostic) {
                diagnostic->reset();
            }
            if (acquisition) {
                acquisition->reset();
            }
            return fail(Code::Cancelled);
        }
        if (diagnostic && !diagnostic->reason && (!acquisition || !acquisition->reason)) {
            diagnostic->reason = LeafReferenceFailure::Graph;
        }
        return std::unexpected(error);
    };
    if (limits.walk.sources > index_max - 2 || limits.walk.edges > index_max ||
        report->unique_sources > std::min(limits.walk.sources, limits.inventory.maximum_sources) ||
        report->declarations > index_max) {
        return graph_failure({Code::ResourceExceeded, operation});
    }
    auto initial = limits.walk;
    initial.sources = report->unique_sources;
    auto walk = DependencyWalk::create(initial, plan, stop);
    if (!walk) {
        return graph_failure(walk.error());
    }
    ReferenceGraphBuilder builder{*origin,   *parent, prepared, limits,      std::move(*walk),
                                  extractor, plan,    stop,     acquisition, diagnostic};
    if (const auto collected = builder.collect(*report); !collected) {
        return graph_failure(collected.error());
    }
    auto result = builder.publish(std::move(prepared));
    if (!result) {
        return graph_failure(result.error());
    }
    return result;
}

core::Result<ReferencedCase>
resolve_leaf_references(const ReferenceInput& input, PreparedRootReferences&& prepared,
                        RootInventoryLimits intake, DependencyWalkLimits limits,
                        LeafFormatLimits formats, core::ResourcePlan& plan,
                        const std::stop_token& stop, RootSnapshotDiagnostic* acquisition,
                        LeafReferenceDiagnostic* diagnostic) noexcept {
    LeafInspector inspector{formats, limits.diagnostic_bytes, plan, stop, diagnostic};
    const auto report = prepared.report();
    const auto path_bytes = report ? report->normalized_path_bytes : 0;
    return resolve_reference_graph(input, std::move(prepared), {intake, limits, path_bytes, {}},
                                   inspector, plan, stop, acquisition, diagnostic);
}

} // namespace uvlas::io::detail
