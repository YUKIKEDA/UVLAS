#include "global_geometry_private.hpp"
#include "query_coordinates_private.hpp"

#include <utility>

namespace uvlas::aero::detail {
using namespace global_detail;
namespace {
bool same_words(ScalarWords a, ScalarWords b) noexcept {
    if (a.index() != b.index())
        return false;
    const auto bits = [](double x) { return std::bit_cast<std::uint64_t>(x); };
    if (const auto* single = std::get_if<SingleWord>(&a)) {
        const auto* other = std::get_if<SingleWord>(&b);
        return other != nullptr && bits(single->value) == bits(other->value);
    }
    const auto* left = std::get_if<DoubleWord>(&a);
    const auto* right = std::get_if<DoubleWord>(&b);
    return left != nullptr && right != nullptr && bits(left->high) == bits(right->high) &&
           bits(left->low) == bits(right->low);
}
std::uint64_t height(std::uint64_t n) {
    std::uint64_t h = 0;
    for (auto v = n > 0 ? n - 1 : 0; v != 0; v >>= 1U)
        ++h;
    return h;
}
void nondegenerate(const GeometryData& g, Configuration config, Run& r) {
    for (const auto& face : g.facets.values()) {
        Vector a, b, c;
        for (std::size_t k = 0; k < 3; ++k) {
            a[k] = coordinate(g, face.points[0], k, config, r);
            b[k] = coordinate(g, face.points[1], k, config, r);
            c[k] = coordinate(g, face.points[2], k, config, r);
        }
        const auto normal = cross(subtract(b, a, r), subtract(c, a, r), r);
        if (normal[0] == 0 && normal[1] == 0 && normal[2] == 0)
            reject(core::FailureCode::GeometryInvalid, FacetId{face.id});
    }
}
std::uint64_t linear_geometry_work(std::uint64_t f) {
    return times(
        f, plus(times(18, decode_cost), plus(times(24, simple_cost), times(6, multiply_cost))));
}
void same_model(const GeometryData& a, const GeometryData& b, Run& r) {
    auto ca = a.context, cb = b.context;
    ca.current = cb.current;
    // Reference compatibility ignores both trial time representations. The
    // current geometry owner retains its own exact-time identity.
    ca.time = cb.time;
    ca.diagnostics = cb.diagnostics;
    if (ca != cb)
        reject(core::FailureCode::StaleGeneration);
    if (a.points.size() != b.points.size() || a.edges.size() != b.edges.size() ||
        a.facets.size() != b.facets.size() || a.faces.size() != b.faces.size() ||
        a.bands.size() != b.bands.size() || a.quads.size() != b.quads.size() ||
        a.triangles.size() != b.triangles.size() || a.reference_words != b.reference_words)
        reject(core::FailureCode::InvalidInput);
    for (std::size_t i = 0; i < a.reference_high.size(); ++i) {
        r.step();
        if (std::bit_cast<std::uint64_t>(a.reference_high.values()[i]) !=
            std::bit_cast<std::uint64_t>(b.reference_high.values()[i]))
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.reference_low.size(); ++i) {
        r.step();
        if (std::bit_cast<std::uint64_t>(a.reference_low.values()[i]) !=
            std::bit_cast<std::uint64_t>(b.reference_low.values()[i]))
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.points.size(); ++i) {
        r.step();
        const auto& x = a.points.values()[i];
        const auto& y = b.points.values()[i];
        if (x.id != y.id || x.owner != y.owner || x.map != y.map ||
            x.map_generation != y.map_generation || x.original.index() != y.original.index())
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.edges.size(); ++i) {
        r.step();
        const auto& x = a.edges.values()[i];
        const auto& y = b.edges.values()[i];
        if (x.id != y.id || x.original != y.original || x.points != y.points)
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.facets.size(); ++i) {
        r.step();
        const auto& x = a.facets.values()[i];
        const auto& y = b.facets.values()[i];
        if (x.id != y.id || x.face != y.face || x.points != y.points || x.edges != y.edges ||
            x.incidence != y.incidence || x.origin.index() != y.origin.index())
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.faces.size(); ++i) {
        r.step();
        if (a.faces.values()[i].id != b.faces.values()[i].id ||
            a.faces.values()[i].band != b.faces.values()[i].band)
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.bands.size(); ++i) {
        r.step();
        if (a.bands.values()[i] != b.bands.values()[i])
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.quads.size(); ++i) {
        r.step();
        const auto& x = a.quads.values()[i];
        const auto& y = b.quads.values()[i];
        if (x.parent != y.parent || x.points != y.points || x.edges != y.edges ||
            x.diagonal != y.diagonal || x.facets != y.facets)
            reject(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < a.triangles.size(); ++i) {
        r.step();
        const auto& x = a.triangles.values()[i];
        const auto& y = b.triangles.values()[i];
        if (x.id != y.id || x.parent != y.parent || x.points != y.points || x.edges != y.edges)
            reject(core::FailureCode::InvalidInput);
    }
    if (!same_words(a.radius, b.radius))
        reject(core::FailureCode::InvalidInput);
}
template <class T, class Compare>
void sorted(std::span<T> values, Compare less, Run& r) {
    const auto compare = [&](const auto& a, const auto& b) {
        r.step();
        return less(a, b);
    };
    std::make_heap(values.begin(), values.end(), compare);
    std::sort_heap(values.begin(), values.end(), compare);
}
} // namespace

GlobalResult<GlobalTopology> GlobalEngine::build(const OwnedGeometry& source,
                                                 const GlobalTopologyInput& input,
                                                 GlobalLimits limits, GlobalWorkspace& workspace,
                                                 core::ResourcePlan& plan,
                                                 std::stop_token stop) noexcept {
    const auto context = source.context();
    return guarded(context, [&]() -> GlobalResult<GlobalTopology> {
        Lease const lease{workspace.busy_};
        if (!source.data_)
            reject(core::FailureCode::InvalidInput);
        const auto& g = *source.data_;
        const auto f = g.facets.size(), e = g.edges.size(), p = g.points.size(),
                   te = input.trailing_edges.size();
        if (input.roles.size() != g.faces.size() || input.edges.size() != e)
            reject(core::FailureCode::InvalidInput);
        const auto nodes = tree_shape(f), h = height(f);
        const auto model_count = plus(plus(p, e), plus(f, g.faces.size()));
        // Balanced median sort at each depth: <=2n(h+1) sifts, <=2 comparisons
        // each. A centre comparison uses <=8 word decodes, 2 low-word adds,
        // 2 centre sums and 1 integer compare. Union/axis/leaf passes are linear.
        const auto tree_work =
            times(plus(nodes, 1), times(times(plus(h, 1), plus(h, 1)),
                                        plus(times(128, decode_cost), times(96, simple_cost))));
        const auto topology_work =
            plus(tree_work,
                 plus(linear_geometry_work(f),
                      plus(times(plus(model_count, plus(te, 1)), 256),
                           plus(times(times(te, f), 3),
                                times(times(plus(te, 1), plus(height(plus(te, 1)), 1)), 32)))));
        Run run{limits, context, workspace.scratch_.values()[0], stop};
        run.admit(topology_work);
        const auto meta =
            metadata_admit(limits.metadata, sizeof(GlobalTopologyOwner) + 4 * sizeof(void*));
        auto amount = bytes<GlobalTopologyData>(1);
        amount = plus(amount, bytes<std::uint8_t>(g.faces.size()));
        amount = plus(amount, bytes<BoundEdge>(e));
        amount = plus(amount, bytes<Node>(nodes));
        amount = plus(amount, bytes<std::uint64_t>(f));
        amount = plus(amount, bytes<std::uint64_t>(times(te, 2)));
        amount = plus(amount, bytes<TrailingEdgeRow>(te));
        byte_admit(amount, plan, context);
        const auto before = plan.allocated_bytes();
        core::PlannedBuffer<GlobalTopologyData> payload;
        allocate(payload, 1, plan, context);
        auto& data = payload.values()[0];
        data.geometry.emplace(retain(source));
        allocate(data.roles, g.faces.size(), plan, context);
        allocate(data.edges, e, plan, context);
        allocate(data.nodes, nodes, plan, context);
        allocate(data.order, f, plan, context);
        allocate(data.trailing_points, times(te, 2), plan, context);
        allocate(data.trailing, te, plan, context);
        for (const auto& role : input.roles) {
            const auto i = lookup(g.face_lookup.values(), role.face.value(), run, role.face);
            if (data.roles.values()[i] != 0 || static_cast<unsigned>(role.role) > 1)
                reject(core::FailureCode::InvalidInput, role.face);
            data.roles.values()[i] =
                role.role == AeroFaceRole::Bound ? GlobalTopologyData::bound_role : 2;
        }
        for (const auto& binding : input.edges) {
            const auto i = lookup(g.edge_lookup.values(), binding.edge.value(), run, binding.edge);
            auto& edge = data.edges.values()[i];
            if (edge.seen)
                reject(core::FailureCode::InvalidInput, binding.edge);
            edge.seen = true;
            edge.original = binding.edge.value();
            edge.canonical = std::holds_alternative<CanonicalEdge>(binding.ownership);
            edge.owner =
                std::holds_alternative<CanonicalEdge>(binding.ownership)
                    ? i
                    : lookup(g.edge_lookup.values(), std::get<EdgeId>(binding.ownership).value(),
                             run, std::get<EdgeId>(binding.ownership));
            const auto& original = g.edges.values()[i];
            const auto raw0 = lookup(g.point_lookup.values(), original.original[0], run);
            edge.raw_basis = g.points.values()[raw0].owner == original.points[0] ? 1 : -1;
        }
        for (std::size_t i = 0; i < e; ++i) {
            run.step();
            const auto& binding = data.edges.values()[i];
            const auto& owner = data.edges.values()[static_cast<std::size_t>(binding.owner)];
            if (!binding.seen || !owner.seen || !owner.canonical || owner.owner != binding.owner ||
                g.edges.values()[i].points !=
                    g.edges.values()[static_cast<std::size_t>(binding.owner)].points)
                reject(core::FailureCode::InvalidInput, EdgeId{binding.original});
        }
        // A1 incidence is ALREADY canonical. Derive raw->canonical from original
        // endpoints and check it, never multiply A1 incidence by the basis twice.
        for (const auto& facet : g.facets.values()) {
            std::array<int, 3> boundary{};
            for (std::size_t k = 0; k < 3; ++k) {
                run.step();
                const auto edge_index = static_cast<std::size_t>(facet.edges[k]);
                const auto& edge = g.edges.values()[edge_index];
                const auto basis = data.edges.values()[edge_index].raw_basis;
                const auto raw_from = basis == 1 ? edge.points[0] : edge.points[1];
                const auto raw_incidence = facet.points[k] == raw_from ? 1 : -1;
                if (raw_incidence * basis != facet.incidence[k])
                    reject(core::FailureCode::InvalidInput, FacetId{facet.id});
                const auto& canonical =
                    g.edges
                        .values()[static_cast<std::size_t>(data.edges.values()[edge_index].owner)];
                for (std::size_t vertex = 0; vertex < 3; ++vertex) {
                    run.step();
                    if (facet.points[vertex] == canonical.points[0])
                        boundary[vertex] -= facet.incidence[k];
                    if (facet.points[vertex] == canonical.points[1])
                        boundary[vertex] += facet.incidence[k];
                }
            }
            if (std::any_of(boundary.begin(), boundary.end(), [](int entry) { return entry != 0; }))
                reject(core::FailureCode::InvalidInput, FacetId{facet.id});
        }
        for (std::size_t i = 0; i < te; ++i) {
            const auto& request = input.trailing_edges[i];
            const auto raw =
                lookup(g.edge_lookup.values(), request.edge.value(), run, request.edge);
            const auto owner = static_cast<std::size_t>(data.edges.values()[raw].owner);
            const auto& edge = g.edges.values()[owner];
            const auto from =
                g.points
                    .values()[lookup(g.point_lookup.values(), request.points[0].value(), run,
                                     request.points[0])]
                    .owner;
            const auto to = g.points
                                .values()[lookup(g.point_lookup.values(), request.points[1].value(),
                                                 run, request.points[1])]
                                .owner;
            const int sign = from == edge.points[0] && to == edge.points[1]   ? 1
                             : from == edge.points[1] && to == edge.points[0] ? -1
                                                                              : 0;
            if (sign == 0)
                reject(core::FailureCode::InvalidInput, request.edge);
            std::size_t count = 0;
            std::uint64_t face_id = 0;
            int incidence = 0;
            for (const auto& facet : g.facets.values())
                for (std::size_t k = 0; k < 3; ++k) {
                    run.step();
                    if (data.roles.values()[static_cast<std::size_t>(facet.face)] ==
                            GlobalTopologyData::bound_role &&
                        data.edges.values()[static_cast<std::size_t>(facet.edges[k])].owner ==
                            owner) {
                        ++count;
                        face_id = g.faces.values()[static_cast<std::size_t>(facet.face)].id;
                        // An orientation coefficient must retain -1.
                        // NOLINTNEXTLINE(bugprone-signed-char-misuse)
                        incidence = facet.incidence[k];
                    }
                }
            if (count != 1)
                reject(core::FailureCode::InvalidInput, request.edge);
            auto& row = data.trailing.values()[i];
            row.edge = EdgeId{edge.id};
            row.bound_face = PhysicalFaceId{face_id};
            row.wake_sign = static_cast<std::int8_t>(sign);
            row.coefficient = static_cast<std::int8_t>(sign * incidence);
            if (-incidence + sign * row.coefficient != 0)
                reject(core::FailureCode::InvalidInput, request.edge);
            data.trailing_points.values()[2 * i] =
                g.points.values()[static_cast<std::size_t>(edge.points[0])].id;
            data.trailing_points.values()[2 * i + 1] =
                g.points.values()[static_cast<std::size_t>(edge.points[1])].id;
        }
        sorted(
            data.trailing.values(), [](const auto& a, const auto& b) { return a.edge < b.edge; },
            run);
        for (std::size_t i = 1; i < te; ++i) {
            run.step();
            if (data.trailing.values()[i - 1].edge == data.trailing.values()[i].edge)
                reject(core::FailureCode::InvalidInput, data.trailing.values()[i].edge);
        }
        sorted(data.trailing_points.values(), std::less<>{}, run);
        for (const auto id : data.trailing_points.values()) {
            run.step();
            if (data.trailing_point_size == 0 ||
                data.trailing_points
                        .values()[static_cast<std::size_t>(data.trailing_point_size - 1)] != id)
                data.trailing_points
                    .values()[static_cast<std::size_t>(data.trailing_point_size++)] = id;
        }
        nondegenerate(g, Configuration::Reference, run);
        build_tree(data, run);
        data.plan = {plan.allocated_bytes() - before, run.reserved, run.used, 0, meta};
        auto owner = std::make_shared<GlobalTopologyOwner>(std::move(payload));
        // The final allocator can request cancellation or reenter user code.
        // Observe that boundary before publishing an owning qualification.
        run.step(0);
        std::shared_ptr<const GlobalTopologyData> shared{owner, &owner->payload.values()[0]};
        return GlobalTopology{std::move(shared), context};
    });
}
GlobalResult<GlobalContactState>
GlobalEngine::check(const GlobalTopology& topology, const OwnedGeometry& source,
                    GlobalLimits limits, GlobalWorkspace& workspace, core::ResourcePlan& plan,
                    std::stop_token stop) noexcept {
    const auto context = source.context();
    return guarded(context, [&]() -> GlobalResult<GlobalContactState> {
        Lease const lease{workspace.busy_};
        if (!source.data_ || !topology.data_)
            reject(core::FailureCode::InvalidInput);
        const auto& g = *source.data_;
        const auto& owner = topology.data_->geometry;
        if (!owner)
            reject(core::FailureCode::InvalidInput);
        const auto& reference = *GlobalEngine::source(*owner);
        const auto f = g.facets.size(), n = topology.data_->nodes.size();
        Run run{limits, context, workspace.scratch_.values()[0], stop};
        const auto model_count =
            plus(plus(g.points.size(), g.edges.size()),
                 plus(g.facets.size(),
                      plus(g.faces.size(),
                           plus(g.bands.size(), plus(g.quads.size(), g.triangles.size())))));
        const auto v = times(n, plus(n, 1)) / 2;
        const auto traversal_budget = times(2, times(v, node_cost));
        const auto refit = times(plus(n, 1), plus(times(64, decode_cost), times(32, simple_cost)));
        run.admit(plus(traversal_budget, plus(refit, plus(linear_geometry_work(f),
                                                          times(plus(model_count, 1), 32)))));
        same_model(reference, g, run);
        metadata_admit(limits.metadata, 0);
        const auto stack_size = pair_stack(f);
        const auto amount =
            plus(bytes<GlobalContactData>(1), plus(bytes<Node>(n), bytes<NodePair>(stack_size)));
        byte_admit(amount, plan, context);
        const auto before = plan.allocated_bytes();
        core::PlannedBuffer<GlobalContactData> payload;
        core::PlannedBuffer<NodePair> stack;
        core::PlannedBuffer<Node> current_nodes;
        allocate(payload, 1, plan, context);
        allocate(stack, stack_size, plan, context);
        auto& data = payload.values()[0];
        data.geometry.emplace(retain(source));
        data.topology = topology.data_;
        allocate(current_nodes, n, plan, context);
        if (n != 0)
            std::copy(topology.data_->nodes.values().begin(), topology.data_->nodes.values().end(),
                      current_nodes.values().begin());
        nondegenerate(g, Configuration::Current, run);
        refit_tree(g, current_nodes.values(), run);
        const auto candidates =
            traverse(g, *topology.data_, current_nodes.values(), stack.values(), false, run);
        run.admit(times(candidates.candidates, narrow_cost));
        const auto repeated =
            traverse(g, *topology.data_, current_nodes.values(), stack.values(), true, run);
        if (candidates.candidates != repeated.candidates || candidates.steps != repeated.steps)
            reject(core::FailureCode::InvalidInput);
        data.plan = {plan.allocated_bytes() - before, run.reserved, run.used, candidates.candidates,
                     0};
        data.plan.traversal_reserved_steps = traversal_budget;
        data.plan.traversal_executed_steps = plus(candidates.steps, repeated.steps);
        if (data.plan.traversal_executed_steps > traversal_budget)
            reject(core::FailureCode::ResourceExceeded);
        return GlobalContactState{std::move(payload), context};
    });
}

GlobalTopology::GlobalTopology(std::shared_ptr<const GlobalTopologyData> d,
                               GeometryContext c) noexcept
    : data_(std::move(d)), context_(c) {}
GlobalTopology::~GlobalTopology() = default;
bool GlobalTopology::active() const noexcept {
    return static_cast<bool>(data_);
}
std::span<const TrailingEdgeRow> GlobalTopology::trailing_edges() const& noexcept {
    return data_ ? data_->trailing.values() : std::span<const TrailingEdgeRow>{};
}
std::size_t GlobalTopology::trailing_point_count() const noexcept {
    return data_ ? static_cast<std::size_t>(data_->trailing_point_size) : 0;
}
GlobalResult<PointId> GlobalTopology::trailing_point(std::size_t i) const noexcept {
    if (!data_ || i >= data_->trailing_point_size)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    return PointId{data_->trailing_points.values()[i]};
}
GlobalResult<MaterialReleaseKey> GlobalTopology::release_key(PointId p,
                                                             ReleaseEpoch epoch) const noexcept {
    if (!data_)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    const auto& owner = data_->geometry;
    if (!owner)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    auto canonical = owner->owner(p);
    if (!canonical)
        return std::unexpected(GlobalGeometryFailure{canonical.error(), {}});
    const auto values =
        data_->trailing_points.values().first(static_cast<std::size_t>(data_->trailing_point_size));
    if (!std::binary_search(values.begin(), values.end(), canonical->value()))
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_, p));
    return MaterialReleaseKey{*canonical, epoch};
}
GlobalResult<std::reference_wrapper<const OwnedGeometry>>
GlobalTopology::geometry() const& noexcept {
    if (!data_)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    const auto& owner = data_->geometry;
    if (!owner)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    return std::cref(*owner);
}
GlobalResult<GlobalWorkPlan> GlobalTopology::work_plan() const noexcept {
    if (!data_)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    return data_->plan;
}
GlobalContactState::GlobalContactState(core::PlannedBuffer<GlobalContactData> d,
                                       GeometryContext c) noexcept
    : data_(std::move(d)), context_(c) {}
GlobalContactState::~GlobalContactState() = default;
GlobalContactState::GlobalContactState(GlobalContactState&&) noexcept = default;
GlobalContactState& GlobalContactState::operator=(GlobalContactState&&) noexcept = default;
bool GlobalContactState::active() const noexcept {
    return !data_.values().empty();
}
GlobalResult<std::reference_wrapper<const OwnedGeometry>>
GlobalContactState::geometry() const& noexcept {
    if (!active())
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    const auto& owner = data_.values()[0].geometry;
    if (!owner)
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    return std::cref(*owner);
}

core::Result<bool> PreparedTopologyOwner::is_trailing_point(PointId canonical) const noexcept {
    if (!topology_)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::EvaluateAero});
    const auto ids = topology_->trailing_points.values().first(
        static_cast<std::size_t>(topology_->trailing_point_size));
    return std::ranges::binary_search(ids, canonical.value());
}

std::size_t PreparedTopologyOwner::trailing_point_count() const noexcept {
    return topology_ ? static_cast<std::size_t>(topology_->trailing_point_size) : 0;
}

GlobalResult<GlobalContactState>
GlobalContactState::retain_checked(core::ResourcePlan& plan,
                                   const std::stop_token& stop) const noexcept {
    return guarded(context_, [&]() -> GlobalResult<GlobalContactState> {
        if (stop.stop_requested())
            reject(core::FailureCode::Cancelled);
        if (!active())
            reject(core::FailureCode::InvalidInput);
        const auto& source = data_.values()[0];
        if (!source.geometry || !source.topology)
            reject(core::FailureCode::InvalidInput);
        byte_admit(bytes<GlobalContactData>(1), plan, context_);
        core::PlannedBuffer<GlobalContactData> payload;
        allocate(payload, 1, plan, context_);
        auto& result = payload.values()[0];
        result.geometry.emplace(GlobalEngine::retain(*source.geometry));
        result.topology = source.topology;
        result.plan = source.plan;
        if (stop.stop_requested())
            reject(core::FailureCode::Cancelled);
        return GlobalContactState{std::move(payload), context_};
    });
}

core::Result<PreparedTopologyOwner> GlobalContactState::reference_topology() const noexcept {
    const auto state = data_.values();
    if (state.empty() || !state.front().topology)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput,
                                             core::Operation::EvaluateAero, context_.diagnostics});
    return PreparedTopologyOwner{state.front().topology};
}

GlobalResult<UvlmQuad> GlobalContactState::uvlm_quad(FacetId id, std::uint64_t max_steps,
                                                     const std::stop_token& stop) const noexcept {
    return guarded(context_, [&]() -> GlobalResult<UvlmQuad> {
        const auto geometry_owner = geometry();
        if (!geometry_owner)
            return std::unexpected(geometry_owner.error());
        const auto& topology = data_.values()[0].topology;
        const auto& source = GlobalEngine::source(geometry_owner->get());
        if (!topology || !source)
            reject(core::FailureCode::InvalidInput, id);
        std::uint64_t used = 0;
        const auto step = [&] {
            if (stop.stop_requested())
                reject(core::FailureCode::Cancelled, id);
            if (used == max_steps)
                reject(core::FailureCode::ResourceExceeded, id);
            ++used;
        };
        const auto find = [&](std::span<const IdLookup> rows, std::uint64_t key) {
            const auto it = std::lower_bound(rows.begin(), rows.end(), key,
                                             [&](const IdLookup& row, std::uint64_t value) {
                                                 step();
                                                 return row.id < value;
                                             });
            if (it == rows.end() || it->id != key)
                reject(core::FailureCode::InvalidInput, id);
            return static_cast<std::size_t>(it->ordinal);
        };
        step();
        const auto& facet =
            source->facets.values()[find(source->facet_lookup.values(), id.value())];
        const auto* origin = std::get_if<QuadOrigin>(&facet.origin);
        if (origin == nullptr)
            reject(core::FailureCode::InvalidInput, id);
        const auto& quad = source->quads.values()[static_cast<std::size_t>(origin->source)];
        const bool bound = topology->roles.values()[static_cast<std::size_t>(facet.face)] ==
                           GlobalTopologyData::bound_role;
        std::optional<TrailingEdgeRow> attachment;
        if (!bound) {
            const auto rows = topology->trailing.values();
            for (std::size_t side = 0; side < quad.edges.size(); ++side) {
                step();
                const auto raw = find(source->edge_lookup.values(), quad.edges[side]);
                const auto owner = static_cast<std::size_t>(topology->edges.values()[raw].owner);
                const auto& edge = source->edges.values()[owner];
                const auto entry = std::lower_bound(rows.begin(), rows.end(), EdgeId{edge.id},
                                                    [&](const TrailingEdgeRow& row, EdgeId value) {
                                                        step();
                                                        return row.edge < value;
                                                    });
                if (entry == rows.end() || entry->edge.value() != edge.id)
                    continue;
                if (side != 0)
                    reject(core::FailureCode::GeometryInvalid, id);
                const auto from =
                    source->points.values()[find(source->point_lookup.values(), quad.points[0])]
                        .owner;
                const bool positive = from == edge.points[0];
                if ((entry->wake_sign == 1) != positive)
                    reject(core::FailureCode::GeometryInvalid, id);
                attachment = *entry;
            }
        }
        return UvlmQuad{PhysicalFaceId{quad.parent},
                        bound ? AeroFaceRole::Bound : AeroFaceRole::Wake,
                        {FacetId{quad.facets[0]}, FacetId{quad.facets[1]}},
                        {PointId{quad.points[0]}, PointId{quad.points[1]}, PointId{quad.points[2]},
                         PointId{quad.points[3]}},
                        attachment,
                        used};
    });
}
GlobalResult<FaceRoleInput> GlobalContactState::face_role(std::size_t index) const noexcept {
    return guarded(context_, [&]() -> GlobalResult<FaceRoleInput> {
        const auto owner = geometry();
        if (!owner)
            return std::unexpected(owner.error());
        const auto data = data_.values();
        const auto& topology = data.front().topology;
        const auto& source = GlobalEngine::source(owner->get());
        if (!topology || !source || index >= source->faces.size())
            reject(core::FailureCode::InvalidInput);
        const auto id = PhysicalFaceId{source->faces.values()[index].id};
        const auto role = topology->roles.values()[index] == GlobalTopologyData::bound_role
                              ? AeroFaceRole::Bound
                              : AeroFaceRole::Wake;
        return FaceRoleInput{id, role};
    });
}

GlobalResult<PhysicalFaceId>
GlobalContactState::original_triangle_face(std::size_t index) const noexcept {
    return guarded(context_, [&]() -> GlobalResult<PhysicalFaceId> {
        const auto owner = geometry();
        if (!owner)
            return std::unexpected(owner.error());
        const auto& source = GlobalEngine::source(owner->get());
        if (!source || index >= source->triangles.size())
            reject(core::FailureCode::InvalidInput);
        return PhysicalFaceId{source->triangles.values()[index].parent};
    });
}

GlobalResult<UvlmOriginalQuad>
GlobalContactState::uvlm_original_quad(std::size_t index, std::uint64_t maximum,
                                       const std::stop_token& stop) const noexcept {
    return guarded(context_, [&]() -> GlobalResult<UvlmOriginalQuad> {
        const auto owner = geometry();
        if (!owner)
            return std::unexpected(owner.error());
        const auto data = data_.values();
        const auto& topology = data.front().topology;
        const auto& source = GlobalEngine::source(owner->get());
        if (!topology || !source || index >= source->quads.size())
            reject(core::FailureCode::InvalidInput);
        const auto& quad = source->quads.values()[index];
        const auto face = PhysicalFaceId{quad.parent};
        std::uint64_t used = 0;
        const auto step = [&] {
            if (stop.stop_requested())
                reject(core::FailureCode::Cancelled, face);
            if (used == maximum)
                reject(core::FailureCode::ResourceExceeded, face);
            ++used;
        };
        const auto faces = source->face_lookup.values();
        step();
        const auto found = std::lower_bound(faces.begin(), faces.end(), face.value(),
                                            [&](const IdLookup& row, std::uint64_t value) {
                                                step();
                                                return row.id < value;
                                            });
        if (found == faces.end() || found->id != face.value())
            reject(core::FailureCode::InvalidInput, face);
        const auto bound = topology->roles.values()[static_cast<std::size_t>(found->ordinal)] ==
                           GlobalTopologyData::bound_role;
        return UvlmOriginalQuad{face,
                                bound ? AeroFaceRole::Bound : AeroFaceRole::Wake,
                                {PointId{quad.points[0]}, PointId{quad.points[1]},
                                 PointId{quad.points[2]}, PointId{quad.points[3]}},
                                {FacetId{quad.facets[0]}, FacetId{quad.facets[1]}},
                                used};
    });
}
std::size_t GlobalContactState::trailing_edge_count() const noexcept {
    if (!active() || !data_.values()[0].topology)
        return 0;
    return data_.values()[0].topology->trailing.size();
}

GlobalResult<UvlmLoadEdge>
GlobalContactState::uvlm_load_edge(EdgeId id, std::uint64_t max_steps,
                                   const std::stop_token& stop) const noexcept {
    return guarded(context_, [&]() -> GlobalResult<UvlmLoadEdge> {
        const auto geometry_owner = geometry();
        if (!geometry_owner)
            return std::unexpected(geometry_owner.error());
        const auto& source = GlobalEngine::source(geometry_owner->get());
        const auto& topology = data_.values()[0].topology;
        if (!source || !topology)
            reject(core::FailureCode::InvalidInput, id);
        std::uint64_t used = 0;
        const auto step = [&] {
            if (stop.stop_requested())
                reject(core::FailureCode::Cancelled, id);
            if (used == max_steps)
                reject(core::FailureCode::ResourceExceeded, id);
            ++used;
        };
        step();
        const auto lookup = source->edge_lookup.values();
        const auto entry = std::lower_bound(lookup.begin(), lookup.end(), id.value(),
                                            [&](const IdLookup& value, std::uint64_t key) {
                                                step();
                                                return value.id < key;
                                            });
        if (entry == lookup.end() || entry->id != id.value())
            reject(core::FailureCode::InvalidInput, id);
        const auto owner = topology->edges.values()[static_cast<std::size_t>(entry->ordinal)].owner;
        const auto& edge = source->edges.values()[static_cast<std::size_t>(owner)];
        const auto key = EdgeId{edge.id};
        const auto rows = topology->trailing.values();
        const auto te = std::lower_bound(rows.begin(), rows.end(), key,
                                         [&](const TrailingEdgeRow& row, EdgeId edge_id) {
                                             step();
                                             return row.edge < edge_id;
                                         });
        const auto points = source->points.values();
        return UvlmLoadEdge{key,
                            {PointId{points[static_cast<std::size_t>(edge.points[0])].id},
                             PointId{points[static_cast<std::size_t>(edge.points[1])].id}},
                            te != rows.end() && te->edge == key,
                            used};
    });
}

GlobalResult<std::pair<EdgeId, std::uint64_t>>
GlobalContactState::facet_edge_alias(FacetId target, EdgeId canonical, std::uint64_t max_steps,
                                     const std::stop_token& stop) const noexcept {
    return guarded(context_, [&]() -> GlobalResult<std::pair<EdgeId, std::uint64_t>> {
        const auto geometry_owner = geometry();
        if (!geometry_owner)
            return std::unexpected(geometry_owner.error());
        const auto& source = GlobalEngine::source(geometry_owner->get());
        const auto& topology = data_.values()[0].topology;
        if (!source || !topology)
            reject(core::FailureCode::InvalidInput, target);
        std::uint64_t used = 0;
        const auto step = [&] {
            if (stop.stop_requested())
                reject(core::FailureCode::Cancelled, target);
            if (used == max_steps)
                reject(core::FailureCode::ResourceExceeded, target);
            ++used;
        };
        step();
        const auto lookup = source->facet_lookup.values();
        const auto entry = std::lower_bound(lookup.begin(), lookup.end(), target.value(),
                                            [&](const IdLookup& value, std::uint64_t key) {
                                                step();
                                                return value.id < key;
                                            });
        if (entry == lookup.end() || entry->id != target.value())
            reject(core::FailureCode::InvalidInput, target);
        const auto& facet = source->facets.values()[static_cast<std::size_t>(entry->ordinal)];
        const auto edges = source->edges.values();
        std::optional<EdgeId> alias;
        for (const auto index : facet.edges) {
            step();
            const auto raw = static_cast<std::size_t>(index);
            const auto owner = static_cast<std::size_t>(topology->edges.values()[raw].owner);
            if (raw != owner && edges[owner].id == canonical.value()) {
                if (alias)
                    reject(core::FailureCode::GeometryInvalid, target);
                alias = EdgeId{edges[raw].id};
            }
        }
        if (!alias)
            reject(core::FailureCode::GeometryInvalid, target);
        return std::pair{*alias, used};
    });
}

GlobalResult<GlobalWorkPlan> GlobalContactState::work_plan() const noexcept {
    if (!active())
        return std::unexpected(failure(core::FailureCode::InvalidInput, context_));
    return data_.values()[0].plan;
}

ContourContactIndex::ContourContactIndex(GlobalContactState&& contact,
                                         core::PlannedBuffer<Node>&& nodes,
                                         std::uint64_t steps) noexcept
    : contact_(std::move(contact)), nodes_(std::move(nodes)), preparation_steps_(steps) {}
ContourContactIndex::ContourContactIndex(ContourContactIndex&&) noexcept = default;
ContourContactIndex& ContourContactIndex::operator=(ContourContactIndex&&) noexcept = default;
ContourContactIndex::~ContourContactIndex() = default;

GlobalResult<ContourContactIndex> ContourContactIndex::prepare(const GlobalContactState& contact,
                                                               GlobalLimits limits,
                                                               GlobalWorkspace& workspace,
                                                               core::ResourcePlan& plan,
                                                               std::stop_token stop) noexcept {
    return GlobalEngine::prepare_contour(contact, limits, workspace, plan, std::move(stop));
}

std::uint64_t ContourContactIndex::preparation_steps() const noexcept {
    return contact_.active() ? preparation_steps_ : 0;
}

GlobalResult<ContourContactWork>
ContourContactIndex::qualify(const ResolvedQuery& query, std::optional<EdgeId> self,
                             const GeometryLimits& query_limits, GlobalLimits limits,
                             GeometryWorkspace& query_workspace, GlobalWorkspace& workspace,
                             std::stop_token stop) const noexcept {
    return GlobalEngine::qualify_contour(*this, query, self, query_limits, limits, query_workspace,
                                         workspace, std::move(stop));
}

GlobalResult<ContourContactIndex> GlobalEngine::prepare_contour(const GlobalContactState& contact,
                                                                GlobalLimits limits,
                                                                GlobalWorkspace& workspace,
                                                                core::ResourcePlan& plan,
                                                                std::stop_token stop) noexcept {
    const auto context = contact.context_;
    return guarded(context, [&]() -> GlobalResult<ContourContactIndex> {
        const Lease lease{workspace.busy_};
        const auto geometry = contact.geometry();
        if (!geometry)
            return std::unexpected(geometry.error());
        const auto& data = *source(geometry->get());
        const auto& reference_nodes = contact.data_.values().front().topology->nodes;
        const auto count = reference_nodes.size();
        Run run{limits, context, workspace.scratch_.values().front(), stop};
        const auto refit_cost = plus(times(64, decode_cost), times(32, simple_cost));
        run.admit(plus(times(plus(count, 1), refit_cost), times(2, decode_cost)));
        if (scalar(data.radius, 0, run) <= 0)
            reject(core::FailureCode::Unsupported);
        byte_admit(plus(bytes<GlobalContactData>(1), bytes<Node>(count)), plan, context);
        auto retained = contact.retain_checked(plan, stop);
        if (!retained)
            return std::unexpected(retained.error());
        core::PlannedBuffer<Node> nodes;
        allocate(nodes, count, plan, context);
        run.step(count);
        std::ranges::copy(reference_nodes.values(), nodes.values().begin());
        refit_tree(data, nodes.values(), run);
        run.step(0);
        return ContourContactIndex{std::move(*retained), std::move(nodes), run.used};
    });
}

GlobalResult<ContourContactWork>
GlobalEngine::qualify_contour(const ContourContactIndex& index, const ResolvedQuery& query,
                              std::optional<EdgeId> self, const GeometryLimits& query_limits,
                              GlobalLimits limits, GeometryWorkspace& query_workspace,
                              GlobalWorkspace& workspace, std::stop_token stop) noexcept {
    const auto& contact = index.contact_;
    const auto context = contact.context_;
    return guarded(context, [&]() -> GlobalResult<ContourContactWork> {
        const Lease lease{workspace.busy_};
        const auto geometry = contact.geometry();
        if (!geometry)
            return std::unexpected(geometry.error());
        const auto& data = *source(geometry->get());
        const auto& topology = *contact.data_.values().front().topology;
        Run run{limits, context, workspace.scratch_.values().front(), stop};
        run.admit(limits.max_steps);
        const auto expression = query.definition();
        if (!expression)
            return std::unexpected(GlobalGeometryFailure{expression.error(), {}});
        auto capped = query_limits;
        capped.max_steps = std::min(capped.max_steps, run.reserved - run.used);
        const auto coordinates =
            exact_query_coordinates(geometry->get(), query, capped, query_workspace, stop);
        if (!coordinates)
            return std::unexpected(GlobalGeometryFailure{coordinates.error(), {}});
        run.step(coordinates->used_steps);
        if (self) {
            const auto* edge = std::get_if<SourceEdgePoint>(&expression->definition);
            if (edge == nullptr || scalar(edge->t, 0, run) != word(0.5, 0, run))
                reject(core::FailureCode::InvalidInput, expression->id);
            const auto owner = contact.uvlm_load_edge(edge->edge, run.reserved - run.used, stop);
            if (!owner)
                return std::unexpected(owner.error());
            run.step(owner->used_steps);
            if (owner->edge != *self || owner->trailing)
                reject(core::FailureCode::InvalidInput, expression->id);
        } else if (!std::holds_alternative<FreePoint>(expression->definition)) {
            reject(core::FailureCode::Unsupported, expression->id);
        }
        return qualify_contour_point(data, topology, index.nodes_.values(), coordinates->numerator,
                                     expression->id, self, run);
    });
}

GlobalWorkspace::GlobalWorkspace(core::PlannedBuffer<GlobalScratch> s) noexcept
    : scratch_(std::move(s)) {}
GlobalWorkspace::~GlobalWorkspace() = default;
std::uint64_t GlobalWorkspace::payload_bytes() noexcept {
    return sizeof(GlobalScratch);
}
core::Result<std::unique_ptr<GlobalWorkspace>>
GlobalWorkspace::create(core::ResourcePlan& plan, numerics::linear::BackendMemoryLimit limit,
                        core::FailureDetails details) noexcept {
    try {
        metadata_admit(limit, sizeof(GlobalWorkspace) + 2 * sizeof(void*));
        auto s = plan.allocate<GlobalScratch>(1, operation);
        if (!s)
            return std::unexpected(with_details(s.error(), details));
        return std::unique_ptr<GlobalWorkspace>{new GlobalWorkspace{std::move(*s)}};
    } catch (const Abort& e) {
        return std::unexpected(core::Failure{e.code, operation, details});
    } catch (...) {
        return std::unexpected(
            core::Failure{core::FailureCode::ResourceExceeded, operation, details});
    }
}
GlobalResult<GlobalTopology> build_global_topology(const OwnedGeometry& g,
                                                   const GlobalTopologyInput& i, GlobalLimits l,
                                                   GlobalWorkspace& w, core::ResourcePlan& p,
                                                   std::stop_token s) noexcept {
    return GlobalEngine::build(g, i, l, w, p, std::move(s));
}
GlobalResult<GlobalContactState> check_global_geometry(const GlobalTopology& t,
                                                       const OwnedGeometry& g, GlobalLimits l,
                                                       GlobalWorkspace& w, core::ResourcePlan& p,
                                                       std::stop_token s) noexcept {
    return GlobalEngine::check(t, g, l, w, p, std::move(s));
}
GlobalResult<AliasKinematicsState>
check_alias_kinematics(const OwnedGeometry& g, const AliasKinematicsInput& i, GlobalLimits l,
                       GlobalWorkspace& w, core::ResourcePlan& p, std::stop_token s) noexcept {
    return GlobalEngine::kinematics(g, i, l, w, p, std::move(s));
}
} // namespace uvlas::aero::detail
