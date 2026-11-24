#include "external_mesh_sidecar_binding.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = MeshSidecarFailure;
constexpr auto operation = core::Operation::PrepareModel;

struct Rejected {
    core::Failure error;
};

class Binding {
public:
    Binding(const MeshSidecarBindingLimits& configured, core::ResourcePlan& resources,
            const std::stop_token& token, std::optional<MeshSidecarDiagnostic>* diagnostic) noexcept
        : limits(configured), plan(resources), stop(token), diagnostic_(diagnostic) {}

    void step() {
        if (stop.stop_requested()) {
            reject(Reason::Resource, Code::Cancelled);
        }
        if (work == limits.work) {
            reject(Reason::Resource, Code::ResourceExceeded);
        }
        ++work;
    }

    [[noreturn]] void reject(Reason reason, Code code = Code::InvalidInput) const {
        report({code, operation}, reason);
    }

    template <class T>
    [[nodiscard]] T value(core::Result<T> result, Reason reason) const {
        if (!result) {
            report(result.error(), reason);
        }
        return *std::move(result);
    }

    template <class T>
    [[nodiscard]] T require(std::optional<T> result, Reason reason) const {
        if (!result) {
            reject(reason);
        }
        return *std::move(result);
    }

    [[nodiscard]] core::PlannedBuffer<BoundExternalFace> allocate(std::uint64_t count) const {
        if (count > limits.faces || !std::in_range<core::Index>(count)) {
            reject(Reason::Resource, Code::ResourceExceeded);
        }
        return value(plan.allocate<BoundExternalFace>(static_cast<core::Index>(count), operation),
                     Reason::Resource);
    }

    void row(core::Index node, std::optional<GmshKey> physical = {},
             std::optional<std::uint64_t> element = {}) noexcept {
        detail.node = node;
        detail.physical = physical;
        detail.element = element;
        detail.byte.reset();
        detail.target_namespace.reset();
        detail.target.reset();
    }

    void target(SchemaNamespace name_space, std::uint64_t id) noexcept {
        detail.target_namespace = name_space;
        detail.target = id;
    }

    const MeshSidecarBindingLimits& limits;
    core::ResourcePlan& plan;
    const std::stop_token& stop;
    MeshSidecarDiagnostic detail;
    std::uint64_t work = 0;

private:
    [[noreturn]] void report(core::Failure error, Reason reason) const {
        if (diagnostic_ && error.code != Code::Cancelled) {
            *diagnostic_ = detail;
            (*diagnostic_)->reason = reason;
        }
        throw Rejected{error};
    }

    std::optional<MeshSidecarDiagnostic>* diagnostic_;
};

template <class Record>
[[nodiscard]] std::uint64_t record_id(const Record& record) noexcept {
    if constexpr (requires { record.id; }) {
        return record.id.unsigned_value();
    } else {
        return std::visit([](const auto& value) { return value.id.unsigned_value(); }, record);
    }
}

template <class Sequence>
[[nodiscard]] auto find_case(Binding& b, const Sequence& values, SchemaNamespace name_space,
                             std::uint64_t id) {
    b.target(name_space, id);
    for (core::Index i = 0; std::cmp_less(i, values.size()); ++i) {
        b.step();
        const auto item = b.value(values.at(i), Reason::CaseReference).record();
        if (record_id(item) == id) {
            return item;
        }
    }
    b.reject(Reason::CaseReference);
}

[[nodiscard]] bool matches_digest(std::string_view hex, const SourceDigest& digest) noexcept {
    constexpr std::string_view chars = "0123456789abcdef";
    if (hex.size() != 2 * digest.bytes.size()) {
        return false;
    }
    for (std::size_t i = 0; i < digest.bytes.size(); ++i) {
        const auto value = std::to_integer<unsigned>(digest.bytes[i]);
        if (hex[2 * i] != chars[value >> 4U] || hex[2 * i + 1] != chars[value & 15U]) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] const AcquiredRootSnapshot& authenticated(Binding& b, const ReferencedCase& graph,
                                                        const SchemaIdentity& parent,
                                                        const schema::FileReference& file,
                                                        ReferenceRole role) {
    const auto source = SourceId{file.source_id.unsigned_value()};
    const auto reference = b.value(graph.find(parent, source, role, b.stop), Reason::Binding);
    const auto& snapshot = reference.get();
    b.detail.source = snapshot.location;
    const auto actual = b.require(snapshot.file.view(), Reason::Snapshot);
    if (!matches_digest(file.sha256.text(), actual.digest)) {
        b.reject(Reason::Binding);
    }
    return reference;
}

[[nodiscard]] ParsedMeshSidecar parse(Binding& b, const AcquiredRootSnapshot& file,
                                      MeshGeometryId geometry) {
    std::optional<MeshSidecarDiagnostic> diagnostic;
    auto parsed = read_mesh_sidecar(file, b.limits.sidecar, b.plan, b.stop, &diagnostic);
    if (diagnostic) {
        b.detail = *diagnostic;
    }
    auto result = b.value(std::move(parsed), b.detail.reason);
    if (result.geometry() != geometry) {
        b.target(SchemaNamespace::Geometry, geometry.value());
        b.reject(Reason::Binding);
    }
    return result;
}

void check_observations(Binding& b, const schema::Case& root, const ParsedMeshSidecar& sidecar) {
    if (const auto& budget = sidecar.conversion()) {
        for (const auto& observation : budget->observations.values()) {
            b.row(observation.node);
            static_cast<void>(find_case(b, root.observations, SchemaNamespace::Observation,
                                        observation.id.value()));
        }
    }
}

[[nodiscard]] schema::Shell shell(Binding& b, const schema::Model& model, MeshShellId id) {
    const auto found = find_case(b, model.shells, SchemaNamespace::Shell, id.value());
    const auto pg =
        find_case(b, model.geometry, SchemaNamespace::Geometry, found.geometry.unsigned_value());
    if (!std::holds_alternative<schema::GeometryPg>(pg)) {
        b.reject(Reason::CaseReference);
    }
    return found;
}

void gmsh_layer_ids(Binding& b, const schema::Shell& target, const GmshSidecarShell& row) {
    if (const auto& layers = row.material_layers) {
        if (layers->size() != static_cast<std::size_t>(target.layers.size())) {
            b.reject(Reason::CaseReference);
        }
        for (std::size_t i = 0; i < layers->size(); ++i) {
            b.step();
            const auto& declared = layers->values()[i];
            const auto actual =
                b.value(target.layers.at(static_cast<core::Index>(i)), Reason::CaseReference)
                    .record();
            if (declared.layer.value() != actual.id.unsigned_value() ||
                declared.material.value() != actual.material.unsigned_value()) {
                b.row(declared.node, row.physical);
                b.target(SchemaNamespace::Shell, row.shell.value());
                b.reject(Reason::CaseReference);
            }
        }
    }
}

[[nodiscard]] bool member(Binding& b, const GmshMeshView& mesh, const GmshElement& element,
                          GmshKey group) {
    if (element.entity.dimension != group.dimension) {
        return false;
    }
    const auto compare = [&b](const auto& lhs, const auto& rhs) {
        b.step();
        return lhs < rhs;
    };
    const auto entity =
        std::ranges::lower_bound(mesh.entities, element.entity, compare, &GmshEntity::key);
    if (entity == mesh.entities.end() || entity->key != element.entity) {
        return false;
    }
    const auto links = mesh.entity_links.subspan(entity->first_physical, entity->physical_count);
    return std::ranges::any_of(links, [&](const auto& link) {
        b.step();
        return link.key == group;
    });
}

[[nodiscard]] std::size_t element_index(Binding& b, const GmshMeshView& mesh, std::uint64_t tag) {
    const auto found = std::ranges::lower_bound(
        mesh.elements, tag,
        [&b](auto lhs, auto rhs) {
            b.step();
            return lhs < rhs;
        },
        &GmshElement::tag);
    if (found == mesh.elements.end() || found->tag != tag) {
        b.reject(Reason::Selection);
    }
    return static_cast<std::size_t>(found - mesh.elements.begin());
}

[[nodiscard]] bool incident(const GmshElement& element, const BoundExternalFace& face) noexcept {
    const std::span face_nodes{face.nodes.data(), face.count};
    const auto nodes = element.connectivity();
    if (element.kind == GmshElementKind::Point) {
        return std::ranges::find(face_nodes, nodes.front()) != face_nodes.end();
    }
    if (element.kind != GmshElementKind::Line) {
        return face.element == element.tag;
    }
    for (std::size_t i = 0; i < face_nodes.size(); ++i) {
        const auto first = face_nodes[i], second = face_nodes[(i + 1) % face_nodes.size()];
        if ((first == nodes[0] && second == nodes[1]) ||
            (first == nodes[1] && second == nodes[0])) {
            return true;
        }
    }
    return false;
}

void gmsh_boundary_group(Binding& b, const GmshMeshView& mesh, GmshKey group,
                         std::span<const BoundExternalFace> faces) {
    bool nonempty = false;
    for (const auto& element : mesh.elements) {
        b.step();
        if (!member(b, mesh, element, group)) {
            continue;
        }
        nonempty = true;
        const auto touching = std::ranges::any_of(faces, [&](const auto& face) {
            b.step();
            return incident(element, face);
        });
        if (!touching) {
            b.detail.element = element.tag;
            b.reject(Reason::Selection);
        }
    }
    if (!nonempty) {
        b.reject(Reason::Selection);
    }
}

[[nodiscard]] core::PlannedBuffer<BoundExternalFace> gmsh_faces(Binding& b,
                                                                const schema::Model& model,
                                                                const GmshSidecar& input,
                                                                const GmshMeshView& mesh) {
    if (!mesh.entities_declared) {
        b.reject(Reason::Selection);
    }
    std::uint64_t count = 0;
    for (const auto& selection : input.shells.values()) {
        b.step();
        if (selection.faces.size() > std::numeric_limits<std::uint64_t>::max() - count) {
            b.reject(Reason::Resource, Code::ResourceExceeded);
        }
        count += selection.faces.size();
    }
    if (count == 0) {
        b.reject(Reason::Selection);
    }
    auto result = b.allocate(count);
    if (!std::in_range<core::Index>(mesh.elements.size())) {
        b.reject(Reason::Resource, Code::ResourceExceeded);
    }
    auto owners = b.value(b.plan.allocate<std::optional<GmshKey>>(
                              static_cast<core::Index>(mesh.elements.size()), operation),
                          Reason::Resource);
    std::size_t next = 0;
    const auto selections = input.shells.values();
    for (std::size_t row = 0; row < selections.size(); ++row) {
        const auto& selection = selections[row];
        b.row(selection.node, selection.physical);
        const auto target = shell(b, model, selection.shell);
        gmsh_layer_ids(b, target, selection);
        const auto start = next;
        const auto declared_faces = selection.faces.values();
        for (std::size_t declaration = 0; declaration < declared_faces.size(); ++declaration) {
            const auto& face = declared_faces[declaration];
            b.row(face.node, selection.physical, face.element);
            b.target(SchemaNamespace::Shell, selection.shell.value());
            const auto i = element_index(b, mesh, face.element);
            const auto& original = mesh.elements[i];
            if (!member(b, mesh, original, selection.physical) || owners.values()[i] ||
                original.connectivity().size() != face.chart.count) {
                b.reject(Reason::Selection);
            }
            owners.values()[i] = selection.physical;
            result.values()[next++] = {selection.shell,
                                       MeshGeometryId{target.geometry.unsigned_value()},
                                       original.tag,
                                       original.nodes,
                                       face.chart.count,
                                       face.node,
                                       i,
                                       row,
                                       declaration};
        }
        std::ranges::sort(result.values().subspan(start, next - start),
                          [&](const auto& a, const auto& other) {
                              b.step();
                              return a.element < other.element;
                          });
    }
    for (const auto& selection : input.shells.values()) {
        b.row(selection.node, selection.physical);
        for (std::size_t i = 0; i < mesh.elements.size(); ++i) {
            b.step();
            const auto& element = mesh.elements[i];
            if (owners.values()[i] != selection.physical &&
                member(b, mesh, element, selection.physical)) {
                b.detail.element = element.tag;
                b.reject(Reason::Selection);
            }
        }
    }
    const auto supports = input.supports.values();
    for (std::size_t i = 0; i < supports.size(); ++i) {
        const auto& support = supports[i];
        b.row(support.node, support.physical);
        static_cast<void>(
            find_case(b, model.supports, SchemaNamespace::Support, support.support.value()));
        if (std::ranges::any_of(supports.first(i), [&](const auto& previous) {
                b.step();
                return previous.physical == support.physical && previous.support == support.support;
            })) {
            b.reject(Reason::Selection);
        }
        gmsh_boundary_group(b, mesh, support.physical, result.values());
    }
    const auto joints = input.joints.values();
    for (std::size_t i = 0; i < joints.size(); ++i) {
        const auto& joint = joints[i];
        b.row(joint.node, joint.first);
        static_cast<void>(find_case(b, model.joints, SchemaNamespace::Joint, joint.joint.value()));
        if (std::ranges::any_of(joints.first(i), [&](const auto& previous) {
                b.step();
                return previous.first == joint.first && previous.second == joint.second &&
                       previous.joint == joint.joint;
            })) {
            b.reject(Reason::Selection);
        }
        gmsh_boundary_group(b, mesh, joint.first, result.values());
        b.detail.physical = joint.second;
        gmsh_boundary_group(b, mesh, joint.second, result.values());
    }
    return result;
}

[[nodiscard]] core::PlannedBuffer<BoundExternalFace> nastran_faces(Binding& b,
                                                                   const schema::Model& model,
                                                                   const NastranSidecar& input,
                                                                   const NastranShellView& deck) {
    if (input.shells.size() > deck.elements.size()) {
        b.reject(Reason::Selection);
    }
    auto result = b.allocate(input.shells.size());
    std::size_t next = 0;
    for (const auto& face : input.shells.values()) {
        b.row(face.node, {}, face.element.value());
        const auto target = shell(b, model, face.shell);
        const auto found = std::ranges::find_if(deck.elements, [&](const auto& element) {
            b.step();
            return element.id == face.element;
        });
        if (found == deck.elements.end() ||
            std::ranges::any_of(result.values().first(next), [&](const auto& previous) {
                b.step();
                return previous.element == face.element.value();
            })) {
            b.reject(Reason::Selection);
        }
        b.detail.byte = found->id_field.position.byte;
        auto& output = result.values()[next];
        output = {face.shell,
                  MeshGeometryId{target.geometry.unsigned_value()},
                  face.element.value(),
                  {},
                  face.chart.count,
                  face.node,
                  found->card_ordinal,
                  next,
                  0};
        ++next;
        std::visit(
            [&](const auto& corners) {
                if (corners.size() != face.chart.count) {
                    b.reject(Reason::Selection);
                }
                for (std::size_t i = 0; i < corners.size(); ++i) {
                    output.nodes[i] = corners[i].id.value();
                }
            },
            found->corners);
    }
    if (result.size() != deck.elements.size()) {
        const auto missing = std::ranges::find_if(deck.elements, [&](const auto& element) {
            b.step();
            return std::ranges::none_of(result.values(), [&](const auto& seen) {
                b.step();
                return seen.element == element.id.value();
            });
        });
        b.detail.node.reset();
        b.detail.element = missing->id.value();
        b.detail.byte = missing->id_field.position.byte;
        b.reject(Reason::Selection);
    }
    return result;
}

[[nodiscard]] std::size_t nastran_card(Binding& b, const NastranBulkView& bulk, std::uint64_t id,
                                       std::string_view first, std::string_view second = {}) {
    std::optional<std::size_t> ordinal;
    for (std::size_t i = 0; i < bulk.cards.size(); ++i) {
        b.step();
        const auto& card = bulk.cards[i];
        const auto name = bulk.text(card.name);
        if (name != first && name != second) {
            continue;
        }
        const auto fields = bulk.data_fields(card);
        if (fields.empty()) {
            continue;
        }
        const auto parsed = bulk.integer(fields.front());
        if (parsed && *parsed > 0 && static_cast<std::uint64_t>(*parsed) == id) {
            b.detail.byte = fields.front().position.byte;
            if (ordinal) {
                b.reject(Reason::Selection);
            }
            ordinal = i;
        }
    }
    return b.require(ordinal, Reason::Selection);
}

template <SchemaNamespace Space, class Sequence>
void nastran_elements(Binding& b, const NastranBulkView& bulk,
                      std::span<const NastranSidecarElementBinding<Space>> rows,
                      const Sequence& targets, std::string_view name) {
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        b.row(row.node, {}, row.element.value());
        static_cast<void>(find_case(b, targets, Space, row.target.value()));
        static_cast<void>(nastran_card(b, bulk, row.element.value(), name));
        if (std::ranges::any_of(rows.first(i), [&](const auto& previous) {
                b.step();
                return previous.element == row.element;
            })) {
            b.reject(Reason::Selection);
        }
    }
}

void nastran_other(Binding& b, const schema::Model& model, const NastranSidecar& input,
                   const NastranBulkView& bulk) {
    const auto materials = input.materials.values();
    for (std::size_t i = 0; i < materials.size(); ++i) {
        const auto& row = materials[i];
        b.row(row.node, {}, row.mid.value());
        static_cast<void>(
            find_case(b, model.materials, SchemaNamespace::Material, row.material.value()));
        static_cast<void>(nastran_card(b, bulk, row.mid.value(), "MAT1", "MAT8"));
        if (std::ranges::any_of(materials.first(i), [&](const auto& previous) {
                b.step();
                return previous.mid == row.mid;
            })) {
            b.reject(Reason::Selection);
        }
    }
    nastran_elements(b, bulk, input.beams.values(), model.beams, "CBEAM");
    nastran_elements(b, bulk, input.masses.values(), model.masses, "CONM2");
    nastran_elements(b, bulk, input.joints.values(), model.joints, "RBE2");
}

void nastran_supports(Binding& b, const schema::Model& model, const NastranSidecar& input,
                      const NastranSupportView& supports) {
    const auto rows = input.supports.values();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        b.row(row.node);
        static_cast<void>(
            find_case(b, model.supports, SchemaNamespace::Support, row.support.value()));
        if (std::ranges::any_of(rows.first(i), [&](const auto& previous) {
                b.step();
                return previous.set == row.set && previous.grid == row.grid;
            })) {
            b.reject(Reason::Selection);
        }
        const auto in_set = [&](NastranSupportSetId set) {
            if (set == row.set) {
                return true;
            }
            return std::ranges::any_of(supports.union_members, [&](const auto& member) {
                b.step();
                return member.combination == row.set && member.member == set;
            });
        };
        unsigned components = 0;
        for (const auto& term : supports.terms) {
            b.step();
            if (!in_set(term.set)) {
                continue;
            }
            const bool includes = std::visit(
                [&](const auto& grid) {
                    if constexpr (requires { grid.grid; }) {
                        return grid.grid == row.grid;
                    } else {
                        return row.grid >= grid.first && row.grid <= grid.last;
                    }
                },
                term.target);
            if (includes) {
                components |= term.components;
                b.detail.byte = term.component_field.position.byte;
            }
        }
        if (components != row.component_bits) {
            b.reject(Reason::Selection);
        }
    }
}

[[nodiscard]] ExternalMeshSidecarBinding bind_gmsh(Binding& b, const schema::Case& root,
                                                   const ReferencedCase& graph,
                                                   const SchemaIdentity& parent,
                                                   const schema::GeometryGmsh& geometry) {
    const auto sidecar_file =
        authenticated(b, graph, parent, geometry.sidecar.record(), ReferenceRole::GmshSidecar);
    auto sidecar = parse(b, sidecar_file.get(), MeshGeometryId{geometry.id.unsigned_value()});
    if (sidecar.conversion() &&
        std::holds_alternative<schema::ReferenceGeometryFit>(geometry.reference.record())) {
        b.reject(Reason::Binding);
    }
    check_observations(b, root, sidecar);
    const auto mesh_file =
        authenticated(b, graph, parent, geometry.mesh.record(), ReferenceRole::GmshMesh);
    auto mesh = b.value(GmshMesh::prepare(mesh_file.get(), b.limits.gmsh, b.plan, b.stop),
                        Reason::Snapshot);
    b.detail.source = sidecar_file.get().location;
    const auto view = b.require(mesh.view(), Reason::Snapshot);
    const auto& declaration = std::get<GmshSidecar>(sidecar.payload());
    auto faces = gmsh_faces(b, root.model.record(), declaration, view);
    return {std::move(sidecar), std::move(mesh), std::move(faces), b.work};
}

[[nodiscard]] ExternalMeshSidecarBinding bind_nastran(Binding& b, const schema::Case& root,
                                                      const ReferencedCase& graph,
                                                      const SchemaIdentity& parent,
                                                      const schema::GeometryNastran& geometry) {
    const auto sidecar_file =
        authenticated(b, graph, parent, geometry.sidecar.record(), ReferenceRole::NastranSidecar);
    auto sidecar = parse(b, sidecar_file.get(), MeshGeometryId{geometry.id.unsigned_value()});
    check_observations(b, root, sidecar);
    const auto deck_file =
        authenticated(b, graph, parent, geometry.deck.record(), ReferenceRole::NastranDeck);
    auto deck = b.value(NastranBulkDeck::prepare(deck_file.get(), b.limits.nastran, b.plan, b.stop),
                        Reason::Snapshot);
    auto coordinates = b.value(
        NastranCoordinateCards::prepare(std::move(deck), b.limits.coordinates, b.plan, b.stop),
        Reason::Snapshot);
    const auto& declaration = std::get<NastranSidecar>(sidecar.payload());
    const auto model = root.model.record();
    if (declaration.supports.size() != 0) {
        auto supports =
            b.value(NastranSupportCards::prepare(coordinates, b.limits.support, b.plan, b.stop),
                    Reason::Snapshot);
        b.detail.source = sidecar_file.get().location;
        const auto view = b.require(supports.view(), Reason::Snapshot);
        const auto coordinate = b.require(coordinates.view(), Reason::Snapshot);
        for (const auto& row : declaration.supports.values()) {
            b.row(row.node);
            if (std::ranges::none_of(coordinate.grids, [&](const auto& grid) {
                    b.step();
                    return grid.id == row.grid;
                })) {
                b.reject(Reason::Selection);
            }
        }
        nastran_supports(b, model, declaration, view);
    }
    b.detail.source = deck_file.get().location;
    auto shells =
        b.value(NastranShellCards::prepare(std::move(coordinates), b.limits.shells, b.plan, b.stop),
                Reason::Snapshot);
    b.detail.source = sidecar_file.get().location;
    const auto view = b.require(shells.view(), Reason::Snapshot);
    auto faces = nastran_faces(b, model, declaration, view);
    nastran_other(b, model, declaration, view.input.bulk);
    return {std::move(sidecar), std::move(shells), std::move(faces), b.work};
}

} // namespace

core::Result<ExternalMeshSidecarBinding>
bind_external_mesh_sidecar(const SchemaCase& schema, const ReferencedCase& graph,
                           MeshGeometryId geometry, const MeshSidecarBindingLimits& limits,
                           core::ResourcePlan& plan, const std::stop_token& stop,
                           std::optional<MeshSidecarDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<ExternalMeshSidecarBinding> {
            try {
                Binding binding{limits, plan, stop, diagnostic};
                binding.step();
                const auto parent = binding.value(schema.identity(), Reason::Binding);
                const auto root = binding.value(schema.root(), Reason::CaseReference);
                const auto selected = find_case(binding, root.model.record().geometry,
                                                SchemaNamespace::Geometry, geometry.value());
                if (const auto* gmsh = std::get_if<schema::GeometryGmsh>(&selected)) {
                    return bind_gmsh(binding, root, graph, parent, *gmsh);
                }
                if (const auto* nastran = std::get_if<schema::GeometryNastran>(&selected)) {
                    return bind_nastran(binding, root, graph, parent, *nastran);
                }
                binding.reject(Reason::CaseReference);
            } catch (const Rejected& rejected) {
                return std::unexpected(rejected.error);
            }
        });
}

} // namespace uvlas::io::detail
