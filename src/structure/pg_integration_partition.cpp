#include "pg_partition_prepare.hpp"

namespace uvlas::structure {
namespace detail {
const PgPartitionData* PgPartitionFactory::data(const PgIntegrationPartition& p) noexcept {
    return p.data_.size() == 1 ? p.data_.values().data() : nullptr;
}
void PgPartitionFactory::info(PgIntegrationPartition& value, PgPartitionPlan info) noexcept {
    value.info_ = info;
}
PgIntegrationPartition PgPartitionFactory::make(core::PlannedBuffer<PgPartitionData> data,
                                                PgPartitionPlan plan) noexcept {
    return PgIntegrationPartition{std::move(data), plan};
}
} // namespace detail
PgIntegrationPartition::PgIntegrationPartition(core::PlannedBuffer<detail::PgPartitionData> data,
                                               PgPartitionPlan info) noexcept
    : data_(std::move(data)), info_(info) {}
PgIntegrationPartition::PgIntegrationPartition(PgIntegrationPartition&&) noexcept = default;
PgIntegrationPartition&
PgIntegrationPartition::operator=(PgIntegrationPartition&&) noexcept = default;
PgIntegrationPartition::~PgIntegrationPartition() = default;
core::Result<PgPartitionPlan>
PgIntegrationPartition::plan(const PgPartitionRequest& request,
                             const PgPartitionLimits& limits) noexcept {
    return detail::partition_entry([&]() -> core::Result<PgPartitionPlan> {
        return detail::partition_budget(request, limits).report;
    });
}
core::Result<PgIntegrationPartition>
PgIntegrationPartition::prepare(const PgPartitionRequest& request, const PgPartitionLimits& limits,
                                core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    core::StageRecorder recorder{limits.algebra.timing == ExactTimingMode::Enabled,
                                 core::Operation::PrepareModel};
    auto result = detail::partition_entry([&]() -> core::Result<PgIntegrationPartition> {
        const core::StageScope phase{recorder, core::MeasurementStage::Preparation};
        return detail::partition_prepare(request, limits, workspace, plan);
    });
    const auto sample = recorder.snapshot();
    if (!sample)
        return std::unexpected(sample.error());
    if (result && recorder.enabled()) {
        auto info = result->preparation();
        if (!info.timings)
            return std::unexpected(
                core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
        info.timings->whole =
            (*sample)[static_cast<std::size_t>(core::MeasurementStage::Preparation)];
        detail::PgPartitionFactory::info(*result, info);
    }
    return result;
}
core::Result<PgGraphView> PgIntegrationPartition::graph() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr || !d->graph)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->graph->view();
}
core::Result<PgPartitionThicknessView> PgIntegrationPartition::thickness() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return PgPartitionThicknessView{
        d->source_points.values().subspan(static_cast<std::size_t>(d->thickness_offset),
                                          d->thickness_values.size()),
        d->thickness_triangles.values(), d->thickness_values.values()};
}
core::Result<PgPartitionIdentity> PgIntegrationPartition::identity() const noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->identity;
}
const PgPartitionPlan& PgIntegrationPartition::preparation() const& noexcept {
    return info_;
}
core::Result<void>
PgIntegrationPartition::validate_for(const PgPartitionRequest& r) const noexcept {
    return detail::partition_entry([&]() -> core::Result<void> {
        const auto* d = detail::PgPartitionFactory::data(*this);
        if (d == nullptr || !d->graph)
            detail::partition_fail(core::FailureCode::InvalidInput);
        const auto a = detail::partition_take(d->graph->view()),
                   b = detail::partition_take(r.graph.view());
        const auto stale = []() { detail::partition_fail(core::FailureCode::StaleGeneration); };
        if (!detail::partition_same_graph(a, b) || r.domain.size() != d->domain_count ||
            r.thickness.vertices.size() != d->thickness_values.size() ||
            r.thickness.values.size() != d->thickness_values.size() ||
            r.thickness.triangles.size() != d->thickness_triangles.size() ||
            r.layers.size() != d->layers.size() || r.boundaries.size() != d->boundaries.size())
            stale();
        const auto points = d->source_points.values();
        for (std::size_t i = 0; i < r.domain.size(); ++i)
            if (!detail::partition_equal_point(r.domain[i], points[i]))
                stale();
        for (std::size_t i = 0; i < r.thickness.vertices.size(); ++i)
            if (!detail::partition_equal_point(
                    r.thickness.vertices[i],
                    points[static_cast<std::size_t>(d->thickness_offset) + i]) ||
                !(r.thickness.values[i] == d->thickness_values.values()[i]))
                stale();
        if (!std::equal(r.thickness.triangles.begin(), r.thickness.triangles.end(),
                        d->thickness_triangles.values().begin()))
            stale();
        for (std::size_t i = 0; i < r.layers.size(); ++i) {
            const auto& x = r.layers[i];
            const auto& y = d->layers.values()[i];
            if (x.id != y.id || !(x.lower == y.lower) || !(x.upper == y.upper))
                stale();
        }
        for (std::size_t i = 0; i < r.boundaries.size(); ++i) {
            const auto& x = r.boundaries[i];
            const auto& y = d->boundaries.values()[i];
            const auto bytes = d->text.values().subspan(static_cast<std::size_t>(y.text_offset),
                                                        static_cast<std::size_t>(y.text_count));
            if (x.name.size() != bytes.size() ||
                !std::equal(x.name.begin(), x.name.end(), bytes.begin()) ||
                x.vertices.size() != static_cast<std::size_t>(y.point_count))
                stale();
            for (std::size_t j = 0; j < x.vertices.size(); ++j)
                if (!detail::partition_equal_point(
                        x.vertices[j], points[static_cast<std::size_t>(y.first_point) + j]))
                    stale();
        }
        std::size_t count = 0;
        for (const auto& source : d->sources.values())
            if (const auto* x = std::get_if<PgPartitionRefinementSource>(&source.descriptor)) {
                if (count >= r.refinements.size())
                    stale();
                const auto& y = r.refinements[count++];
                if (x->id != y.id || x->version != y.version ||
                    !detail::partition_equal_point(
                        points[static_cast<std::size_t>(source.endpoints[0])], y.first) ||
                    !detail::partition_equal_point(
                        points[static_cast<std::size_t>(source.endpoints[1])], y.second))
                    stale();
            }
        if (count != r.refinements.size())
            stale();
        return {};
    });
}
core::Result<std::span<const PgPartitionPoint>> PgIntegrationPartition::vertices() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->vertices.values().first(static_cast<std::size_t>(d->vertex_count));
}

core::Result<std::span<const PgPartitionCell>>
PgIntegrationPartition::activity_map() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->cells.values();
}

core::Result<std::span<const PgPartitionPiece>> PgIntegrationPartition::pieces() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->pieces.values().first(static_cast<std::size_t>(d->piece_count));
}

core::Result<std::span<const core::Index>>
PgIntegrationPartition::piece_vertex_indices() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->piece_vertices.values().first(static_cast<std::size_t>(d->piece_vertex_count));
}

core::Result<std::span<const PgPartitionEdge>> PgIntegrationPartition::edges() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->edges.values().first(static_cast<std::size_t>(d->edge_count));
}

core::Result<std::span<const core::Index>>
PgIntegrationPartition::edge_source_indices() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->edge_sources.values().first(static_cast<std::size_t>(d->edge_source_count));
}

core::Result<std::span<const PgPartitionSource>> PgIntegrationPartition::sources() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->sources.values();
}

core::Result<std::span<const PgPartitionPoint>>
PgIntegrationPartition::source_points() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->source_points.values();
}

core::Result<std::span<const PgPartitionOwnedLayer>>
PgIntegrationPartition::layers() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->layers.values();
}

core::Result<std::span<const PgPartitionBoundaryRecord>>
PgIntegrationPartition::boundaries() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->boundaries.values();
}

core::Result<std::span<const char>> PgIntegrationPartition::boundary_text() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->text.values();
}

core::Result<std::span<const PgPartitionTriangle>>
PgIntegrationPartition::triangles() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->triangles.values().first(static_cast<std::size_t>(d->triangle_count));
}

core::Result<std::span<const PgPartitionVolume>> PgIntegrationPartition::volumes() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->volumes.values().first(static_cast<std::size_t>(d->volume_count));
}

core::Result<std::span<const PgPartitionHorizontalFace>>
PgIntegrationPartition::horizontal_faces() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->horizontal.values().first(static_cast<std::size_t>(d->horizontal_count));
}

core::Result<std::span<const PgPartitionLateralFace>>
PgIntegrationPartition::lateral_faces() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->lateral.values().first(static_cast<std::size_t>(d->lateral_count));
}

core::Result<std::span<const PgPartitionEdgeTrace>>
PgIntegrationPartition::edge_traces() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->traces.values().first(static_cast<std::size_t>(d->trace_count));
}

core::Result<std::span<const PgPartitionContact>>
PgIntegrationPartition::contacts() const& noexcept {
    const auto* d = detail::PgPartitionFactory::data(*this);
    if (d == nullptr)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return d->contacts.values().first(static_cast<std::size_t>(d->contact_count));
}

} // namespace uvlas::structure
