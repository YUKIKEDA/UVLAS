#include "case_pg_reference.hpp"

#include "structure/pg_partition_common.hpp"

#include <array>
#include <memory_resource>
#include <utility>

namespace uvlas::coupling::detail {
namespace {
namespace pg = structure;
namespace schema = io::schema;
using core::FailureCode;
using pg::detail::partition_add;
using pg::detail::partition_fail;
using pg::detail::partition_index;
using pg::detail::partition_mul;
using pg::detail::partition_take;
constexpr auto operation = core::Operation::PrepareModel;
constexpr io::SchemaDimension metres{{1, 0, 0, 0}};
constexpr io::SchemaDimension dimensionless{{0, 0, 0, 0}};

struct InputReader {
    const io::SiCase& source;
    const CasePgReferenceLimits& limits;
    core::ResourcePlan& plan;
    std::uint64_t work_left;

    void charge(std::uint64_t work = 1) {
        if (limits.cancellation.stop_requested() || limits.graph.cancellation.stop_requested() ||
            limits.coefficients.cancellation.stop_requested() ||
            limits.partition.cancellation.stop_requested() ||
            limits.geometry.cancellation.stop_requested())
            partition_fail(FailureCode::Cancelled);
        if (work > work_left)
            partition_fail(FailureCode::ResourceExceeded);
        work_left -= work;
    }

    template <class T>
    core::PlannedBuffer<T> allocate(std::uint64_t size) {
        charge(size);
        return partition_take(plan.allocate<T>(partition_index(size), operation));
    }

    template <class Field>
    pg::ExactRational rational(Field field) {
        const auto value = field.record();
        const auto numerator = value.numerator.text();
        const auto denominator = value.denominator.text();
        charge(partition_add(1, partition_add(numerator.size(), denominator.size())));
        return partition_take(pg::ExactRational::from_integer_pair(numerator, denominator,
                                                                   limits.graph.algebra, plan));
    }

    template <class Field>
    pg::ExactRational quantity(Field field, io::SchemaDimension dimension) {
        charge();
        const auto value = partition_take(source.quantity(field)).recipe;
        if (value.dimensions.exponents != dimension.exponents || value.pi_power != 0)
            partition_fail(FailureCode::InvalidInput);
        charge(std::visit(
            [](const auto& exact) -> std::uint64_t {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(exact)>,
                                             io::SiDecimalView>)
                    return exact.digits.size();
                else
                    return partition_add(exact.numerator.digits.size(),
                                         exact.denominator.digits.size());
            },
            value.value));
        return partition_take(io::materialize_rational(value, limits.graph.algebra, plan));
    }

    template <class Points>
    void read_points(Points source_points, std::span<pg::PgPartitionPoint> target) {
        for (std::size_t i = 0; i < target.size(); ++i) {
            const auto point = partition_take(source_points.at(partition_index(i)));
            target[i] = {rational(partition_take(point.at(0))),
                         rational(partition_take(point.at(1)))};
        }
    }
};

/// A contiguous buffer for views with mandatory references. The vectors reserve
/// once and the arena has no upstream allocator.
template <class T>
std::uint64_t view_storage_bytes(std::uint64_t count) {
    return partition_add(partition_mul(count, sizeof(T)), alignof(T) - 1);
}

schema::Shell find_shell(schema::Model model, std::uint64_t shell_id, InputReader& reader) {
    for (core::Index i = 0; i < partition_index(model.shells.size()); ++i) {
        reader.charge();
        const auto shell = partition_take(model.shells.at(i)).record();
        if (shell.id.unsigned_value() == shell_id)
            return shell;
    }
    partition_fail(FailureCode::InvalidInput);
}

schema::GeometryPg find_pg(schema::Model model, std::uint64_t geometry_id, InputReader& reader) {
    for (core::Index i = 0; i < partition_index(model.geometry.size()); ++i) {
        reader.charge();
        const auto geometry = partition_take(model.geometry.at(i)).record();
        if (std::visit([](const auto& entry) { return entry.id.unsigned_value(); }, geometry) !=
            geometry_id)
            continue;
        if (const auto* pg_geometry = std::get_if<schema::GeometryPg>(&geometry))
            return *pg_geometry;
        partition_fail(FailureCode::Unsupported);
    }
    partition_fail(FailureCode::InvalidInput);
}

pg::PgGraph prepare_graph(schema::GeometryPg geometry, InputReader& reader) {
    const auto count = geometry.coefficients.size();
    auto values = reader.allocate<pg::ExactRational>(count);
    auto storage = reader.allocate<std::byte>(view_storage_bytes<pg::ExactInputTerm>(count));
    std::pmr::monotonic_buffer_resource arena{storage.values().data(), storage.size(),
                                              std::pmr::null_memory_resource()};
    std::pmr::vector<pg::ExactInputTerm> terms{&arena};
    terms.reserve(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        const auto term = partition_take(geometry.coefficients.at(partition_index(i))).record();
        values.values()[i] = reader.rational(term.coefficient);
        terms.push_back({{partition_index(partition_take(term.powers.at(0)).unsigned_value()),
                          partition_index(partition_take(term.powers.at(1)).unsigned_value())},
                         values.values()[i]});
    }
    const auto height = partition_take(
        pg::ExactChartPolynomial::make(terms, reader.limits.graph.algebra, reader.plan));
    const auto x0 = reader.rational(geometry.x0);
    const auto width = reader.rational(geometry.width);
    const auto skew = reader.rational(geometry.skew);
    const pg::PgGraphRequest request{x0,
                                     width,
                                     skew,
                                     height,
                                     partition_index(geometry.n.unsigned_value()),
                                     pg::PgGraphConvention::Version1};
    return partition_take(pg::PgGraph::prepare(request, reader.limits.graph, reader.plan));
}

struct PartitionInput {
    core::PlannedBuffer<pg::PgPartitionPoint> domain;
    core::PlannedBuffer<pg::PgPartitionPoint> thickness_vertices;
    core::PlannedBuffer<std::array<core::Index, 3>> thickness_triangles;
    core::PlannedBuffer<pg::ExactRational> thickness_values;
    core::PlannedBuffer<pg::PgPartitionPoint> boundary_vertices;
    core::PlannedBuffer<pg::PgPartitionBoundary> boundaries;
};

PartitionInput read_partition(schema::GeometryPg geometry, InputReader& reader) {
    const auto thickness = geometry.thickness.record();
    const auto boundary_count = partition_index(geometry.boundaries.size());
    std::uint64_t boundary_points = 0;
    for (core::Index i = 0; i < boundary_count; ++i) {
        reader.charge();
        const auto boundary = partition_take(geometry.boundaries.at(i)).record();
        boundary_points = partition_add(boundary_points, boundary.vertices.size());
    }
    const auto input_points = partition_add(
        partition_add(geometry.domain.size(), thickness.vertices.size()), boundary_points);
    if (input_points > reader.limits.partition.maximum_input_points)
        partition_fail(FailureCode::ResourceExceeded);
    PartitionInput input{
        reader.allocate<pg::PgPartitionPoint>(geometry.domain.size()),
        reader.allocate<pg::PgPartitionPoint>(thickness.vertices.size()),
        reader.allocate<std::array<core::Index, 3>>(thickness.triangles.size()),
        reader.allocate<pg::ExactRational>(thickness.values.size()),
        reader.allocate<pg::PgPartitionPoint>(boundary_points),
        reader.allocate<pg::PgPartitionBoundary>(geometry.boundaries.size()),
    };
    reader.read_points(geometry.domain, input.domain.values());
    reader.read_points(thickness.vertices, input.thickness_vertices.values());
    for (std::size_t i = 0; i < input.thickness_values.size(); ++i)
        input.thickness_values.values()[i] =
            reader.quantity(partition_take(thickness.values.at(partition_index(i))), metres);
    for (std::size_t i = 0; i < input.thickness_triangles.size(); ++i) {
        reader.charge();
        const auto triangle = partition_take(thickness.triangles.at(partition_index(i)));
        for (core::Index j = 0; j < 3; ++j)
            input.thickness_triangles.values()[i][static_cast<std::size_t>(j)] =
                partition_index(partition_take(triangle.at(j)).unsigned_value());
    }
    std::size_t offset = 0;
    for (std::size_t i = 0; i < input.boundaries.size(); ++i) {
        const auto boundary = partition_take(geometry.boundaries.at(partition_index(i))).record();
        const auto count = static_cast<std::size_t>(boundary.vertices.size());
        auto vertices = input.boundary_vertices.values().subspan(offset, count);
        const auto name = boundary.name.text();
        reader.charge(partition_add(1, name.size()));
        reader.read_points(boundary.vertices, vertices);
        input.boundaries.values()[i] = {name, vertices};
        offset += count;
    }
    return input;
}

pg::PgIntegrationPartition prepare_partition(const pg::PgGraph& graph, schema::GeometryPg geometry,
                                             schema::Shell shell, InputReader& reader,
                                             core::Workspace& workspace) {
    const auto count = shell.layers.size();
    if (count > reader.limits.partition.maximum_layers)
        partition_fail(FailureCode::ResourceExceeded);
    auto endpoints = reader.allocate<std::array<pg::ExactRational, 2>>(count);
    auto storage = reader.allocate<std::byte>(view_storage_bytes<pg::PgPartitionLayer>(count));
    std::pmr::monotonic_buffer_resource arena{storage.values().data(), storage.size(),
                                              std::pmr::null_memory_resource()};
    std::pmr::vector<pg::PgPartitionLayer> layers{&arena};
    layers.reserve(endpoints.size());
    for (std::size_t i = 0; i < endpoints.size(); ++i) {
        const auto layer = partition_take(shell.layers.at(partition_index(i))).record();
        auto& endpoint_pair = endpoints.values()[i];
        endpoint_pair[0] = reader.rational(partition_take(layer.interval.at(0)));
        endpoint_pair[1] = reader.rational(partition_take(layer.interval.at(1)));
        layers.emplace_back(core::LayerId{layer.id.unsigned_value()}, endpoint_pair[0],
                            endpoint_pair[1]);
    }
    const auto input = read_partition(geometry, reader);
    const pg::PgPartitionRequest request{
        graph,
        input.domain.values(),
        {input.thickness_vertices.values(), input.thickness_triangles.values(),
         input.thickness_values.values()},
        layers,
        input.boundaries.values(),
        {},
    };
    return partition_take(pg::PgIntegrationPartition::prepare(request, reader.limits.partition,
                                                              workspace, reader.plan));
}
} // namespace

core::Result<pg::PgReferenceGeometry> prepare_case_pg_reference(const io::SiCase& source,
                                                                std::uint64_t shell_id,
                                                                const CasePgReferenceLimits& limits,
                                                                core::Workspace& workspace,
                                                                core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        operation, FailureCode::InvalidInput, {}, [&]() -> core::Result<pg::PgReferenceGeometry> {
            try {
                InputReader reader{source, limits, plan, limits.maximum_input_work};
                reader.charge();
                const auto model = partition_take(partition_take(source.source_case()).get().root())
                                       .model.record();
                const auto shell = find_shell(model, shell_id, reader);
                const auto geometry = find_pg(model, shell.geometry.unsigned_value(), reader);
                const auto frame = geometry.frame.record();
                const auto length = reader.quantity(geometry.reference_length, metres);
                std::array<pg::ExactRational, 3> origin;
                std::array<pg::ExactRational, 9> rotation;
                for (core::Index i = 0; i < 3; ++i) {
                    origin[static_cast<std::size_t>(i)] =
                        reader.quantity(partition_take(frame.origin.at(i)), metres);
                    const auto row = partition_take(frame.rotation.at(i));
                    for (core::Index j = 0; j < 3; ++j)
                        rotation[static_cast<std::size_t>(3 * i + j)] =
                            reader.quantity(partition_take(row.at(j)), dimensionless);
                }
                auto graph = prepare_graph(geometry, reader);
                auto partition = prepare_partition(graph, geometry, shell, reader, workspace);
                auto coefficients = partition_take(pg::PgCoefficientDictionary::prepare(
                    graph, limits.coefficients, workspace, plan));
                auto result = partition_take(pg::PgReferenceGeometry::prepare(
                    {std::move(graph), std::move(coefficients), std::move(partition)},
                    {rotation, origin, length}, limits.geometry, plan));
                reader.charge(0);
                return result;
            } catch (const pg::detail::PartitionAbort& error) {
                return std::unexpected(error.failure);
            }
        });
}
} // namespace uvlas::coupling::detail
