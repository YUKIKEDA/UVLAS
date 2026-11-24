#include "external_mesh_sidecar_pg.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <memory_resource>
#include <type_traits>
#include <utility>
#include <vector>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;
using Code = core::FailureCode;
using Reason = MeshSidecarFailure;

struct Rejected {
    core::Failure error;
};

class Preparation {
public:
    Preparation(const ExternalSidecarPgLimits& configured, core::ResourcePlan& resources,
                const std::stop_token& token, RootSnapshotLocation location,
                ExternalSidecarPgDiagnostic* output) noexcept
        : limits(configured), plan(resources), stop(token), diagnostic(output) {
        detail.source = location;
    }

    void step() {
        if (stop.stop_requested()) {
            reject(Code::Cancelled, Reason::Resource);
        }
        if (work == limits.work) {
            reject(Code::ResourceExceeded, Reason::Resource);
        }
        ++work;
    }

    [[noreturn]] void reject(Code code, Reason why = Reason::Selection) {
        propagate({code, operation}, why);
    }

    template <class T>
    [[nodiscard]] T value(core::Result<T> result, Reason reason = Reason::Binding) {
        if (!result) {
            propagate(result.error(), reason);
        }
        return *std::move(result);
    }

    template <class T>
    [[nodiscard]] T require(std::optional<T> result, Reason reason = Reason::Snapshot) {
        if (!result) {
            reject(Code::InvalidInput, reason);
        }
        return *std::move(result);
    }

    void check(core::Result<void> result) {
        if (!result) {
            propagate(result.error(), Reason::Binding);
        }
    }

    [[noreturn]] void propagate(core::Failure failure, Reason reason) {
        if (diagnostic) {
            if (failure.code == Code::Cancelled) {
                *diagnostic = {};
            } else {
                detail.reason = reason;
                diagnostic->binding = detail;
            }
        }
        throw Rejected{failure};
    }

    const ExternalSidecarPgLimits& limits;
    core::ResourcePlan& plan;
    const std::stop_token& stop;
    ExternalSidecarPgDiagnostic* diagnostic;
    MeshSidecarDiagnostic detail{};
    std::uint64_t work = 0;
};

[[nodiscard]] schema::Shell find_shell(Preparation& p, const schema::Model& model, MeshShellId id) {
    for (core::Index i = 0; std::cmp_less(i, model.shells.size()); ++i) {
        p.step();
        const auto shell = p.value(model.shells.at(i), Reason::CaseReference).record();
        if (shell.id.unsigned_value() == id.value()) {
            return shell;
        }
    }
    p.reject(Code::InvalidInput, Reason::CaseReference);
}

[[nodiscard]] schema::GeometryPg find_pg(Preparation& p, const schema::Model& model,
                                         MeshGeometryId id) {
    for (core::Index i = 0; std::cmp_less(i, model.geometry.size()); ++i) {
        p.step();
        const auto geometry = p.value(model.geometry.at(i), Reason::CaseReference).record();
        const auto* pg = std::get_if<schema::GeometryPg>(&geometry);
        if (pg && pg->id.unsigned_value() == id.value()) {
            return *pg;
        }
    }
    p.reject(Code::InvalidInput, Reason::CaseReference);
}

[[nodiscard]] ImportedCaseMaterial find_material(Preparation& p, const schema::Model& model,
                                                 MeshMaterialId id) {
    for (core::Index i = 0; std::cmp_less(i, model.materials.size()); ++i) {
        p.step();
        const auto material = p.value(model.materials.at(i), Reason::CaseReference).record();
        if (std::visit([&](const auto& record) { return record.id.unsigned_value() == id.value(); },
                       material)) {
            return material;
        }
    }
    p.reject(Code::InvalidInput, Reason::CaseReference);
}

[[nodiscard]] const MeshSidecarChart& chart(const MeshSidecarPayload& payload,
                                            const BoundExternalFace& face) {
    if (const auto* gmsh = std::get_if<GmshSidecar>(&payload)) {
        return gmsh->shells.values()[face.sidecar_row].faces.values()[face.sidecar_face].chart;
    }
    return std::get<NastranSidecar>(payload).shells.values()[face.sidecar_row].chart;
}

[[nodiscard]] std::size_t selected_faces(Preparation& p, const ExternalMeshSidecarBinding& bound,
                                         MeshShellId shell, MeshGeometryId pg) {
    std::size_t count = 0;
    for (const auto& face : bound.faces()) {
        p.step();
        if (face.shell == shell) {
            p.detail.node = face.sidecar_node;
            p.detail.element = face.element;
            if (face.pg != pg) {
                p.reject(Code::InvalidInput, Reason::CaseReference);
            }
            ++count;
        }
    }
    if (count == 0) {
        p.reject(Code::InvalidInput);
    }
    if (count > p.limits.geometry.maximum_faces || !std::in_range<core::Index>(count)) {
        p.reject(Code::ResourceExceeded, Reason::Resource);
    }
    return count;
}

template <class Evaluate>
[[nodiscard]] auto with_face_selection(Preparation& p, const ExternalMeshSidecarBinding& bound,
                                       MeshShellId shell, std::size_t count, Evaluate&& evaluate) {
    constexpr auto per_face = sizeof(ExternalPgFace) + 4 * sizeof(ExternalPgCorner);
    constexpr auto alignment = alignof(ExternalPgFace) + alignof(ExternalPgCorner);
    const auto bytes =
        p.value(core::checked_byte_count(static_cast<core::Index>(count), per_face, operation),
                Reason::Resource);
    const auto limit = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    if (bytes > limit - alignment || count > std::numeric_limits<std::size_t>::max() / 4U) {
        p.reject(Code::ResourceExceeded, Reason::Resource);
    }
    auto memory =
        p.value(p.plan.allocate<std::byte>(static_cast<core::Index>(bytes + alignment), operation),
                Reason::Resource);
    // Store only the reference-bearing views in the charged slab.
    // Iterator-checking containers may allocate separate bookkeeping.
    static_assert(std::is_trivially_destructible_v<ExternalPgFace>);
    static_assert(std::is_trivially_destructible_v<ExternalPgCorner>);
    std::pmr::monotonic_buffer_resource scratch(memory.values().data(), memory.size(),
                                                std::pmr::null_memory_resource());
    auto* faces = std::pmr::polymorphic_allocator<ExternalPgFace>{&scratch}.allocate(count);
    auto* corners =
        std::pmr::polymorphic_allocator<ExternalPgCorner>{&scratch}.allocate(4U * count);
    std::size_t face_index = 0;
    std::size_t corner_index = 0;
    auto tags = p.value(p.plan.allocate<std::uint64_t>(static_cast<core::Index>(count), operation),
                        Reason::Resource);
    const auto& payload = bound.sidecar().payload();
    for (const auto& face : bound.faces()) {
        p.step();
        if (face.shell != shell) {
            continue;
        }
        const auto& mapping = chart(payload, face);
        auto* ring = corners + corner_index;
        for (std::size_t i = 0; i < face.count; ++i) {
            std::construct_at(ring + i, face.nodes[i], mapping.vertices[i]);
        }
        const auto winding = mapping.orientation == MeshSidecarOrientation::Same
                                 ? ExternalPgWinding::Forward
                                 : ExternalPgWinding::Reversed;
        tags.values()[face_index] = face.element;
        std::construct_at(faces + face_index, face.element, winding,
                          std::span<const ExternalPgCorner>{ring, face.count});
        ++face_index;
        corner_index += face.count;
    }

    return evaluate(std::span<const ExternalPgFace>{faces, count},
                    std::span<const std::uint64_t>{tags.values()});
}

[[nodiscard]] GmshSurfaceTopology certify_topology(Preparation& p, const GmshMesh& mesh,
                                                   std::span<const std::uint64_t> tags) {
    GmshTopologyDiagnostic detail;
    auto topology =
        GmshSurfaceTopology::prepare(mesh, tags, p.limits.topology, p.plan, p.stop, &detail);
    if (!topology && p.diagnostic) {
        p.diagnostic->topology = detail;
    }
    return p.value(std::move(topology), Reason::Selection);
}

[[nodiscard]] ExternalPlanarPg establish_geometry(Preparation& p,
                                                  const ExternalMeshSidecarBinding& bound,
                                                  const structure::PgReferenceGeometry& pg,
                                                  MeshShellId shell, std::size_t count) {
    auto limits = p.limits.geometry;
    limits.arithmetic.cancellation = p.stop;
    auto* diagnostics = p.diagnostic ? &p.diagnostic->geometry : nullptr;
    const auto& payload = bound.sidecar().payload();
    return with_face_selection(p, bound, shell, count, [&](auto faces, auto tags) {
        if (const auto* gmsh = std::get_if<GmshMesh>(&bound.source())) {
            const auto topology = certify_topology(p, *gmsh, tags);
            return p.value(ExternalPlanarPg::prepare(pg, *gmsh, topology,
                                                     std::get<GmshSidecar>(payload).length, faces,
                                                     limits, p.plan, diagnostics),
                           Reason::Selection);
        }
        const auto& nastran = std::get<NastranShellCards>(bound.source());
        return p.value(ExternalPlanarPg::prepare(pg, nastran,
                                                 std::get<NastranSidecar>(payload).units.length,
                                                 faces, limits, p.plan, diagnostics),
                       Reason::Selection);
    });
}

using ConversionGeometry = std::variant<GmshPgBounds, NastranPlanarPgBounds>;
using structure::ExactRational;

struct ConversionBudget {
    ExactRational distance;
    ExactRational normal;
    core::PlannedBuffer<ExternalPgConversionObservation> observations;
};

[[nodiscard]] schema::Case check_case(Preparation& p, const SchemaCase& schema,
                                      const ExternalMeshSidecarBinding& bound) {
    p.step();
    const auto parent = p.value(schema.identity());
    const auto source = bound.sidecar().source();
    if (parent.source != source.parent.source ||
        parent.contract_revision != source.parent.contract_revision ||
        parent.representation_revision != source.parent.representation_revision) {
        p.reject(Code::StaleGeneration, Reason::Binding);
    }
    return p.value(schema.root(), Reason::CaseReference);
}

[[nodiscard]] std::optional<schema::ReferenceGeometryFit>
find_fit(Preparation& p, const schema::Model& model, MeshGeometryId id) {
    for (core::Index i = 0; std::cmp_less(i, model.geometry.size()); ++i) {
        p.step();
        const auto geometry = p.value(model.geometry.at(i), Reason::CaseReference).record();
        if (const auto* gmsh = std::get_if<schema::GeometryGmsh>(&geometry);
            gmsh && gmsh->id.unsigned_value() == id.value()) {
            const auto reference = gmsh->reference.record();
            if (const auto* fit = std::get_if<schema::ReferenceGeometryFit>(&reference)) {
                return *fit;
            }
            return {};
        }
    }
    p.reject(Code::InvalidInput, Reason::CaseReference);
}

template <class ReadObservation>
[[nodiscard]] core::PlannedBuffer<ExternalPgConversionObservation>
read_conversion_observations(Preparation& p, std::uint64_t count, ReadObservation read) {
    if (count == 0) {
        p.reject(Code::InvalidInput, Reason::CaseReference);
    }
    if (count > p.limits.work - p.work || !std::in_range<core::Index>(count)) {
        p.reject(Code::ResourceExceeded, Reason::Resource);
    }
    auto rows = p.value(p.plan.allocate<ExternalPgConversionObservation>(
                            static_cast<core::Index>(count), operation),
                        Reason::Resource);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        p.step();
        rows.values()[i].id = read(i);
    }
    return rows;
}

[[nodiscard]] ConversionBudget
read_conversion_budget(Preparation& p, const SiCase& si, const schema::Case& root,
                       const ParsedMeshSidecar& sidecar,
                       const std::optional<schema::ReferenceGeometryFit>& fit) {
    const auto& exact = p.limits.geometry.arithmetic.algebra;
    const auto& declared = sidecar.conversion();
    if (fit && declared) {
        p.reject(Code::InvalidInput, Reason::CaseReference);
    }
    if (fit) {
        const auto distance = p.value(si.quantity(fit->distance_budget), Reason::CaseReference);
        const auto normal = p.value(si.quantity(fit->normal_budget), Reason::CaseReference);
        auto observations = read_conversion_observations(
            p, fit->conversion_observations.size(), [&](std::size_t i) {
                const auto field =
                    p.value(fit->conversion_observations.at(static_cast<core::Index>(i)),
                            Reason::CaseReference);
                return core::ObservationId{field.unsigned_value()};
            });
        return {p.value(materialize_rational(distance.recipe, exact, p.plan)),
                p.value(materialize_rational(normal.recipe, exact, p.plan)),
                std::move(observations)};
    }
    if (!declared) {
        p.reject(Code::Unsupported, Reason::CaseReference);
    }
    const auto length = static_cast<SiLengthUnit>(root.units.record().length.selected());
    const SiUnits units{length, SiMassUnit::Kilogram, SiTimeUnit::Second, SiAngleUnit::Radian};
    const SiRecipeView scale{SiDecimalView{1, "1", 0}, {}, units, {{1, 0, 0, 0}}, 0};
    const auto metres = p.value(materialize_rational(scale, exact, p.plan));
    auto observations =
        read_conversion_observations(p, declared->observations.size(), [&](std::size_t i) {
            return declared->observations.values()[i].id;
        });
    return {p.value(declared->distance.multiply(metres, exact, p.plan)),
            p.value(declared->normal.clone(exact, p.plan)), std::move(observations)};
}

[[nodiscard]] ConversionGeometry bound_conversion_geometry(Preparation& p,
                                                           const ExternalMeshSidecarBinding& bound,
                                                           const structure::PgReferenceGeometry& pg,
                                                           MeshShellId shell, std::size_t count,
                                                           std::uint64_t maximum_directed_work) {
    auto region = p.limits.geometry;
    region.arithmetic.cancellation = p.stop;
    const ExternalPgBoundLimits limits{region, maximum_directed_work};
    auto* diagnostics = p.diagnostic ? &p.diagnostic->geometry : nullptr;
    const auto& payload = bound.sidecar().payload();
    return with_face_selection(
        p, bound, shell, count, [&](auto faces, auto tags) -> ConversionGeometry {
            if (const auto* gmsh = std::get_if<GmshMesh>(&bound.source())) {
                const auto topology = certify_topology(p, *gmsh, tags);
                return p.value(GmshPgBounds::prepare(pg, *gmsh, topology,
                                                     std::get<GmshSidecar>(payload).length, faces,
                                                     limits, p.plan, diagnostics),
                               Reason::Selection);
            }
            const auto& nastran = std::get<NastranShellCards>(bound.source());
            return p.value(NastranPlanarPgBounds::prepare(
                               pg, nastran, std::get<NastranSidecar>(payload).units.length, faces,
                               limits, p.plan, diagnostics),
                           Reason::Selection);
        });
}

void compare_conversion_bounds(Preparation& p, const ConversionGeometry& geometry,
                               const ConversionBudget& budget) {
    const auto& exact = p.limits.geometry.arithmetic.algebra;
    std::visit(
        [&](const auto& owner) {
            const auto view = p.require(owner.view());
            for (const auto& face : view.faces) {
                p.step();
                if (p.value(face.distance_upper_metres.compare(budget.distance, exact, p.plan)) >
                        0 ||
                    p.value(face.normal_chord_upper.compare(budget.normal, exact, p.plan)) > 0) {
                    p.detail.element = face.element;
                    p.reject(Code::PrecisionUnresolved, Reason::Chart);
                }
            }
        },
        geometry);
}

[[nodiscard]] std::optional<ImportedQuantityMismatch>* mismatch(Preparation& p) noexcept {
    return p.diagnostic == nullptr ? nullptr : &p.diagnostic->quantity;
}

[[nodiscard]] structure::ExactRational length_scale(Preparation& p, SiLengthUnit length,
                                                    const structure::ExactLimits& exact) {
    const SiRecipeView recipe{
        SiDecimalView{1, "1", 0},
        {},
        {length, SiMassUnit::Kilogram, SiTimeUnit::Second, SiAngleUnit::Radian},
        SchemaDimension{{1, 0, 0, 0}},
        0};
    return p.value(materialize_rational(recipe, exact, p.plan));
}

[[nodiscard]] structure::ExactLimits
thickness_input_limits(Preparation& p, std::uint64_t uniform_values, std::uint64_t nodal_faces) {
    auto arithmetic = p.limits.geometry.arithmetic;
    arithmetic.cancellation = p.stop;
    return p.value(
        admit_external_pg_thickness_inputs(arithmetic, uniform_values, nodal_faces, p.plan),
        Reason::Resource);
}

template <class Geometry>
[[nodiscard]] std::optional<ExternalPgThicknessReport>
compare_thickness_regions(Preparation& p, const Geometry& corresponding,
                          const structure::PgReferenceGeometry& reference,
                          std::span<const ExternalPgThicknessClaim> claims) {
    if (claims.empty()) {
        return {};
    }
    auto exact = p.limits.geometry.arithmetic;
    exact.cancellation = p.stop;
    const ExternalPgThicknessLimits limits{exact, p.limits.geometry.maximum_faces,
                                           p.limits.case_thickness_vertices};
    std::optional<ExternalPgThicknessDiagnostic> detail;
    auto verified = compare_external_pg_regional_thickness(corresponding, reference, claims, limits,
                                                           p.plan, &detail);
    if (detail) {
        if (const auto position = detail->field) {
            p.detail.byte = position->byte;
        }
        p.detail.node = detail->sidecar_node;
        p.detail.element = detail->element;
        if (p.diagnostic) {
            p.diagnostic->regional_thickness = detail;
            if (detail->reason == ExternalPgThicknessFailure::Difference) {
                p.diagnostic->quantity = ImportedQuantityMismatch{
                    ImportedQuantityField::FullThickness, detail->field, {}};
            }
        }
    }
    return p.value(std::move(verified), Reason::Binding);
}

template <class Geometry>
[[nodiscard]] std::optional<ExternalPgThicknessReport>
gmsh_quantities(Preparation& p, const ExternalMeshSidecarBinding& bound, const GmshSidecar& source,
                MeshShellId shell, const structure::PgReferenceGeometry& reference,
                const Geometry& corresponding) {
    const auto view = p.require(corresponding.view(), Reason::Binding);
    const auto rows = source.shells.values();
    if (std::ranges::none_of(rows, [&](const auto& row) {
            p.step();
            return row.full_thickness.has_value();
        })) {
        return {};
    }
    if (!std::in_range<core::Index>(rows.size()) ||
        !std::in_range<core::Index>(view.faces.size())) {
        p.reject(Code::ResourceExceeded, Reason::Resource);
    }
    const auto exact = thickness_input_limits(p, std::min(rows.size(), view.faces.size()), 0);
    auto claims = p.value(p.plan.allocate<ExternalPgThicknessClaim>(
                              static_cast<core::Index>(view.faces.size()), operation),
                          Reason::Resource);
    auto metres = p.value(
        p.plan.allocate<structure::ExactRational>(static_cast<core::Index>(rows.size()), operation),
        Reason::Resource);
    structure::ExactRational unit;
    std::size_t declared = 0;
    for (const auto& face : bound.faces()) {
        p.step();
        if (face.shell != shell) {
            continue;
        }
        const auto group = rows[face.sidecar_row].physical;
        const auto definition = std::ranges::find_if(rows, [&](const auto& row) {
            p.step();
            return row.physical == group && row.full_thickness.has_value();
        });
        if (definition == rows.end()) {
            continue;
        }
        if (unit.sign() == 0) {
            unit = length_scale(p, source.length, exact);
        }
        const auto index = static_cast<std::size_t>(definition - rows.begin());
        auto& value = metres.values()[index];
        if (value.sign() == 0) {
            const auto& input = definition->full_thickness;
            if (!input) {
                p.reject(Code::InvalidInput, Reason::Binding);
            }
            value = p.value(input->multiply(unit, exact, p.plan));
        }
        claims.values()[declared++] = {face.element, value, definition->full_thickness_node};
    }
    return compare_thickness_regions(p, corresponding, reference, claims.values().first(declared));
}

[[nodiscard]] const NastranSidecarMaterial& mapped_material(Preparation& p,
                                                            const NastranSidecar& source,
                                                            NastranMaterialId id,
                                                            MeshMaterialId expected) {
    const auto materials = source.materials.values();
    const auto found = std::ranges::find_if(materials, [&](const auto& row) {
        p.step();
        return row.mid == id;
    });
    if (found == materials.end()) {
        p.reject(Code::InvalidInput, Reason::CaseReference);
    }
    p.detail.node = found->node;
    p.detail.target_namespace = SchemaNamespace::Material;
    p.detail.target = expected.value();
    if (found->material != expected) {
        p.reject(Code::InvalidInput, Reason::CaseReference);
    }
    return *found;
}

void nastran_material(Preparation& p, const NastranSidecar& source, const NastranLaminateView& deck,
                      const NastranMaterialReference& use, const schema::Model& model,
                      MeshMaterialId expected, const SiUnits& case_units) {
    static_cast<void>(mapped_material(p, source, use.id, expected));
    const auto target = find_material(p, model, expected);
    const auto& declaration = deck.input.materials[use.declaration_ordinal];
    const auto& limits = p.limits.geometry.arithmetic.algebra;
    p.detail.byte = declaration.id_field.position.byte;
    if (declaration.kind == NastranMaterialKind::Mat1) {
        const auto ordinal = p.require(declaration.mat1_ordinal);
        p.check(check_imported_material(deck.input.mat1[ordinal], target, source.units, case_units,
                                        limits, p.plan, p.stop, mismatch(p)));
        return;
    }
    const auto found = std::ranges::find_if(deck.mat8, [&](const auto& value) {
        p.step();
        return value.declaration_ordinal == use.declaration_ordinal;
    });
    if (found == deck.mat8.end()) {
        p.reject(Code::InvalidInput, Reason::CaseReference);
    }
    p.check(check_imported_material(*found, target, source.units, case_units, limits, p.plan,
                                    p.stop, mismatch(p)));
}

struct QuantityCounts {
    bool thickness{};
    std::uint64_t materials{};
    core::PlannedBuffer<std::optional<NastranPcompLayerSection>> laminate{};
    std::optional<ExternalPgThicknessReport> regional{};
};

struct NastranRegionalFace {
    const NastranShellElement* element{};
    core::Index node{};
    NastranShellPropertyKind kind = NastranShellPropertyKind::Pshell;
    MeshSidecarOrientation orientation = MeshSidecarOrientation::Same;
};

[[nodiscard]] bool has_corner_thickness(const NastranShellElement& element) {
    return std::visit(
        [](const auto& ring) {
            return std::ranges::any_of(
                ring, [](const auto& corner) { return corner.thickness.has_value(); });
        },
        element.corners);
}

[[nodiscard]] core::PlannedBuffer<NastranRegionalFace>
nastran_faces(Preparation& p, const ExternalMeshSidecarBinding& bound, const NastranShellView& deck,
              MeshShellId target, const ExternalPlanarPgView& geometry) {
    const auto count = geometry.faces.size();
    if (!std::in_range<core::Index>(count)) {
        p.reject(Code::ResourceExceeded, Reason::Resource);
    }
    auto selected =
        p.value(p.plan.allocate<NastranRegionalFace>(static_cast<core::Index>(count), operation),
                Reason::Resource);
    std::size_t next = 0;
    for (const auto& face : bound.faces()) {
        p.step();
        if (face.shell != target) {
            continue;
        }
        p.detail.node = face.sidecar_node;
        p.detail.element = face.element;
        const auto element = std::ranges::find_if(deck.elements, [&](const auto& value) {
            p.step();
            return value.card_ordinal == face.source_ordinal;
        });
        if (element == deck.elements.end() || next == count ||
            geometry.faces[next].element != face.element) {
            p.reject(Code::InvalidInput, Reason::Binding);
        }
        const auto kind = deck.properties[element->property.property_ordinal].kind;
        const auto orientation = geometry.faces[next].winding == ExternalPgWinding::Forward
                                     ? MeshSidecarOrientation::Same
                                     : MeshSidecarOrientation::Reversed;
        selected.values()[next++] = {&*element, face.sidecar_node, kind, orientation};
    }
    if (next != count) {
        p.reject(Code::InvalidInput, Reason::Binding);
    }
    return selected;
}

[[nodiscard]] NastranPcompLayerSection
compare_laminate(Preparation& p, const NastranLaminateCards& deck, const NastranSidecar& source,
                 std::size_t property, MeshSidecarOrientation winding,
                 const structure::PgReferenceGeometry& reference) {
    const NastranPcompGeometryLimits limits{p.limits.geometry.arithmetic.algebra,
                                            p.limits.case_thickness_vertices,
                                            p.limits.work - p.work};
    std::optional<NastranPcompGeometryDiagnostic> detail;
    auto match = compare_nastran_pcomp_layer_section(deck, property, source.units.length, winding,
                                                     reference, limits, p.plan, p.stop, &detail);
    if (detail) {
        if (p.diagnostic) {
            p.diagnostic->laminate = detail;
        }
        if (const auto field = detail->field) {
            p.detail.byte = field->byte;
        }
        if (const auto layer = detail->case_layer) {
            p.detail.target_namespace = SchemaNamespace::Layer;
            p.detail.target = layer->value();
        }
    }
    auto result = p.value(std::move(match), Reason::Binding);
    p.work += result.work;
    return result;
}

struct PcompRegion {
    std::size_t property_ordinal{};
    std::size_t first_face{};
    MeshSidecarOrientation orientation = MeshSidecarOrientation::Same;
};

void locate_quantity_face(Preparation& p, const NastranRegionalFace& face) {
    p.detail.element = face.element->id.value();
    p.detail.node = face.node;
    p.detail.byte.reset();
}

void add_material_count(Preparation& p, QuantityCounts& quantities, std::uint64_t count) {
    if (count > std::numeric_limits<std::uint64_t>::max() - quantities.materials) {
        p.reject(Code::ResourceExceeded, Reason::Resource);
    }
    quantities.materials += count;
}

void laminate_claims(Preparation& p, const NastranLaminateCards& cards,
                     const NastranLaminateView& deck, const NastranSidecar& source,
                     const schema::Model& model, const schema::Shell& shell,
                     const structure::PgReferenceGeometry& reference, const SiUnits& case_units,
                     std::span<const NastranRegionalFace> faces, core::Index pcomp_faces,
                     std::span<ExternalPgThicknessClaim> claims, QuantityCounts& out) {
    if (pcomp_faces == 0) {
        return;
    }
    auto regions = p.value(p.plan.allocate<PcompRegion>(pcomp_faces, operation), Reason::Resource);
    auto to_region =
        p.value(p.plan.allocate<std::size_t>(static_cast<core::Index>(faces.size()), operation),
                Reason::Resource);
    std::size_t groups = 0;
    for (std::size_t i = 0; i < faces.size(); ++i) {
        p.step();
        const auto& face = faces[i];
        if (face.kind != NastranShellPropertyKind::Pcomp) {
            continue;
        }
        if (has_corner_thickness(*face.element)) {
            locate_quantity_face(p, face);
            p.reject(Code::Unsupported, Reason::Binding);
        }
        const auto ordinal = face.element->property.property_ordinal;
        const auto previous = regions.values().first(groups);
        const auto group = std::ranges::find_if(previous, [&](const PcompRegion& candidate) {
            p.step();
            return candidate.property_ordinal == ordinal &&
                   candidate.orientation == face.orientation;
        });
        if (group == previous.end()) {
            to_region.values()[i] = groups;
            regions.values()[groups++] = {ordinal, i, face.orientation};
        } else {
            to_region.values()[i] = static_cast<std::size_t>(group - previous.begin());
        }
    }
    out.laminate = p.value(p.plan.allocate<std::optional<NastranPcompLayerSection>>(
                               static_cast<core::Index>(groups), operation),
                           Reason::Resource);
    for (std::size_t i = 0; i < groups; ++i) {
        p.step();
        const auto& region = regions.values()[i];
        locate_quantity_face(p, faces[region.first_face]);
        auto layer_section = compare_laminate(p, cards, source, region.property_ordinal,
                                              region.orientation, reference);
        for (const auto& link : layer_section.layers.values()) {
            p.step();
            const auto layer = p.value(shell.layers.at(static_cast<core::Index>(link.case_ordinal)),
                                       Reason::CaseReference)
                                   .record();
            if (layer.id.unsigned_value() != link.layer.value()) {
                p.reject(Code::StaleGeneration, Reason::CaseReference);
            }
            nastran_material(p, source, deck, link.material, model,
                             MeshMaterialId{layer.material.unsigned_value()}, case_units);
            add_material_count(p, out, 1);
        }
        out.laminate.values()[i].emplace(std::move(layer_section));
    }
    for (std::size_t i = 0; i < faces.size(); ++i) {
        p.step();
        const auto& face = faces[i];
        if (face.kind != NastranShellPropertyKind::Pcomp) {
            continue;
        }
        if (const auto& section = out.laminate.values()[to_region.values()[i]]) {
            auto& claim = claims[i];
            claim = {face.element->id.value(), section->full_thickness_metres, face.node,
                     section->property_field};
            claim.property = section->property;
        } else {
            p.reject(Code::InvalidInput, Reason::Binding);
        }
    }
}

struct PshellRegion {
    const NastranPshell* property{};
    std::uint64_t first_element{};
    core::Index node{};
    structure::ExactRational metres;
    bool corner_qualified{};
};

[[nodiscard]] std::uint64_t check_pshell_materials(Preparation& p, const NastranSidecar& source,
                                                   const NastranLaminateView& deck,
                                                   const PshellRegion& region,
                                                   const schema::Model& model,
                                                   const schema::Shell& shell,
                                                   const SiUnits& case_units) {
    p.detail.element = region.first_element;
    p.detail.node = region.node;
    p.detail.byte.reset();
    const auto& property = *region.property;
    if (std::holds_alternative<NastranBendingDesignation>(property.bending)) {
        p.reject(Code::Unsupported, Reason::Binding);
    }
    const auto* bending = std::get_if<NastranMaterialReference>(&property.bending);
    const std::array<const NastranMaterialReference*, 4> uses{
        property.membrane ? &*property.membrane : nullptr, bending,
        property.transverse_shear ? &*property.transverse_shear : nullptr,
        property.coupling ? &*property.coupling : nullptr};
    std::uint64_t comparisons = 0;
    for (const auto* use : uses) {
        p.step();
        if (!use) {
            continue;
        }
        if (shell.layers.size() != 1) {
            p.reject(Code::Unsupported, Reason::Binding);
        }
        const auto layer = p.value(shell.layers.at(0), Reason::CaseReference).record();
        nastran_material(p, source, deck, *use, model,
                         MeshMaterialId{layer.material.unsigned_value()}, case_units);
        ++comparisons;
    }
    return comparisons;
}

struct PshellCornerBuffer {
    std::array<structure::ExactRational, 4> metres;
    std::array<ExternalPgCornerThicknessClaim, 4> fields{};
};

void require_homogeneous_pshell(Preparation& p, const NastranSidecar& source,
                                const NastranLaminateView& deck, const NastranPshell& property,
                                const schema::Model& model, const schema::Shell& shell,
                                const structure::ExactLimits& exact) {
    p.detail.byte = deck.input.input.properties[property.property_ordinal].id_field.position.byte;
    const auto* bending = std::get_if<NastranMaterialReference>(&property.bending);
    if (!property.membrane || !bending || !property.transverse_shear) {
        if (const auto* designation = std::get_if<NastranBendingDesignation>(&property.bending)) {
            p.detail.byte = designation->field.position.byte;
        }
        p.reject(Code::Unsupported, Reason::Binding);
    }
    const auto& membrane = *property.membrane;
    const auto& shear = *property.transverse_shear;
    if (bending->id != membrane.id || shear.id != membrane.id) {
        p.detail.byte = (bending->id != membrane.id ? bending->field : shear.field).position.byte;
        p.reject(Code::Unsupported, Reason::Binding);
    }
    if (const auto& coupling = property.coupling) {
        p.detail.byte = coupling->field.position.byte;
        p.reject(Code::Unsupported, Reason::Binding);
    }
    if (const auto& shear_ratio = property.shear_ratio) {
        p.detail.byte = shear_ratio->field.position.byte;
        p.reject(Code::Unsupported, Reason::Binding);
    }
    if (const auto& ratio = property.bending_ratio) {
        const auto value =
            p.value(structure::ExactRational::from_decimal(ratio->value.decimal(), exact, p.plan));
        const auto one =
            p.value(structure::ExactRational::from_integer_pair("1", "1", exact, p.plan));
        if (value != one) {
            p.detail.byte = ratio->field.position.byte;
            p.reject(Code::Unsupported, Reason::Binding);
        }
    }
    if (shell.layers.size() != 1 ||
        deck.input.materials[membrane.declaration_ordinal].kind != NastranMaterialKind::Mat1) {
        p.detail.byte = membrane.field.position.byte;
        p.reject(Code::Unsupported, Reason::Binding);
    }
    const auto layer = p.value(shell.layers.at(0), Reason::CaseReference).record();
    const auto expected = MeshMaterialId{layer.material.unsigned_value()};
    static_cast<void>(mapped_material(p, source, membrane.id, expected));
    if (!std::holds_alternative<schema::MaterialIsotropic>(find_material(p, model, expected))) {
        p.reject(Code::Unsupported, Reason::Binding);
    }
}

[[nodiscard]] ExternalPgThicknessClaim corner_claim(Preparation& p, const NastranRegionalFace& face,
                                                    const PshellRegion& group,
                                                    const structure::ExactRational& unit,
                                                    PshellCornerBuffer& storage,
                                                    const structure::ExactLimits& exact) {
    const auto& element = *face.element;
    const auto position = group.property->thickness.transform(
        [](const NastranShellReal& source) { return source.field.position; });
    const bool fraction = element.thickness_mode == NastranShellThicknessMode::Fraction;
    const auto mode_field = element.thickness_mode_field.transform(
        [](const NastranBulkLexeme& field) { return field.position; });
    const auto count = std::visit(
        [&](const auto& ring) {
            for (std::size_t i = 0; i < ring.size(); ++i) {
                p.step();
                const auto& corner = ring[i];
                const auto& field = corner.thickness;
                if (field) {
                    p.detail.byte = field->field.position.byte;
                } else if (const auto original = mode_field) {
                    p.detail.byte = original->byte;
                }
                if ((fraction || !field) && !position) {
                    p.reject(Code::InvalidInput, Reason::Binding);
                }
                auto& value = storage.metres[i];
                if (field) {
                    const auto number = p.value(structure::ExactRational::from_decimal(
                        field->value.decimal(), exact, p.plan));
                    value = p.value(number.multiply(fraction ? group.metres : unit, exact, p.plan));
                } else {
                    value = p.value(group.metres.clone(exact, p.plan));
                }
                storage.fields[i] = {corner.id.value(), &value,
                                     field ? std::optional{field->field.position} : position,
                                     (field && fraction) ? position : std::nullopt};
            }
            return ring.size();
        },
        element.corners);
    auto claim = ExternalPgThicknessClaim{element.id.value(),
                                          std::span{storage.fields}.first(count), face.node};
    claim.field = mode_field;
    return claim;
}

struct PshellClaims {
    core::PlannedBuffer<PshellRegion> regions;
    core::PlannedBuffer<PshellCornerBuffer> corners;
    std::size_t groups{};
};

[[nodiscard]] PshellClaims pshell_claims(Preparation& p, const NastranSidecar& source,
                                         const NastranLaminateView& deck,
                                         std::span<const NastranRegionalFace> faces,
                                         core::Index pshell_faces, const schema::Model& model,
                                         const schema::Shell& shell,
                                         std::span<ExternalPgThicknessClaim> claims) {
    if (pshell_faces == 0) {
        return {};
    }
    std::uint64_t corner_count = 0;
    for (const auto& face : faces) {
        p.step();
        if (face.kind == NastranShellPropertyKind::Pshell && has_corner_thickness(*face.element)) {
            ++corner_count;
        }
    }
    const auto exact =
        thickness_input_limits(p, static_cast<std::uint64_t>(pshell_faces), corner_count);
    auto regions =
        p.value(p.plan.allocate<PshellRegion>(pshell_faces, operation), Reason::Resource);
    auto corners = p.value(
        p.plan.allocate<PshellCornerBuffer>(static_cast<core::Index>(corner_count), operation),
        Reason::Resource);
    std::size_t next_corner = 0;
    structure::ExactRational unit;
    std::size_t groups = 0;
    for (std::size_t i = 0; i < faces.size(); ++i) {
        p.step();
        const auto& face = faces[i];
        if (face.kind != NastranShellPropertyKind::Pshell) {
            continue;
        }
        locate_quantity_face(p, face);
        const auto ordinal = face.element->property.property_ordinal;
        const auto previous = regions.values().first(groups);
        auto group = std::ranges::find_if(previous, [&](const PshellRegion& candidate) {
            p.step();
            return candidate.property->property_ordinal == ordinal;
        });
        if (group == previous.end()) {
            const auto property = std::ranges::find_if(deck.input.pshell, [&](const auto& entry) {
                p.step();
                return entry.property_ordinal == ordinal;
            });
            if (property == deck.input.pshell.end()) {
                p.reject(Code::Unsupported, Reason::Binding);
            }
            auto& new_group = regions.values()[groups++];
            new_group.property = &*property;
            new_group.first_element = face.element->id.value();
            new_group.node = face.node;
            if (const auto& declared = property->thickness) {
                if (unit.sign() == 0) {
                    unit = length_scale(p, source.units.length, exact);
                }
                const auto input = p.value(structure::ExactRational::from_decimal(
                    declared->value.decimal(), exact, p.plan));
                new_group.metres = p.value(input.multiply(unit, exact, p.plan));
            }
            group = regions.values().begin() + static_cast<std::ptrdiff_t>(groups - 1);
        }
        if (has_corner_thickness(*face.element)) {
            if (!group->corner_qualified) {
                require_homogeneous_pshell(p, source, deck, *group->property, model, shell, exact);
                group->corner_qualified = true;
            }
            if (unit.sign() == 0) {
                unit = length_scale(p, source.units.length, exact);
            }
            claims[i] = corner_claim(p, face, *group, unit, corners.values()[next_corner++], exact);
        } else if (const auto& declared = group->property->thickness) {
            claims[i] = {face.element->id.value(), group->metres, face.node,
                         declared->field.position};
        } else {
            continue;
        }
        claims[i].property = face.element->property.id;
    }
    return {std::move(regions), std::move(corners), groups};
}

[[nodiscard]] QuantityCounts nastran_quantities(
    Preparation& p, const ExternalMeshSidecarBinding& bound, const NastranSidecar& source,
    const NastranShellView& original, const NastranLaminateCards& cards, const schema::Model& model,
    const schema::Shell& shell, const SiUnits& case_units,
    const structure::PgReferenceGeometry& reference, const ExternalPlanarPg& corresponding) {
    const auto deck = p.require(cards.view());
    if (!original.input.bulk.same_snapshot(deck.input.input.input.bulk.location,
                                           deck.input.input.input.bulk.digest)) {
        p.reject(Code::StaleGeneration, Reason::Binding);
    }
    const auto target = MeshShellId{shell.id.unsigned_value()};
    const auto shape = p.require(corresponding.view(), Reason::Binding);
    const auto selected = nastran_faces(p, bound, deck.input.input, target, shape);
    const auto faces = selected.values();
    const auto count = static_cast<core::Index>(faces.size());
    core::Index pcomp_faces = 0;
    for (const auto& face : faces) {
        p.step();
        if (face.kind == NastranShellPropertyKind::Pcomp) {
            ++pcomp_faces;
        }
    }
    auto claims =
        p.value(p.plan.allocate<ExternalPgThicknessClaim>(count, operation), Reason::Resource);
    const auto pshell =
        pshell_claims(p, source, deck, faces, count - pcomp_faces, model, shell, claims.values());
    QuantityCounts out;
    laminate_claims(p, cards, deck, source, model, shell, reference, case_units, faces, pcomp_faces,
                    claims.values(), out);
    // Compact the claims in source order. An unclaimed PSHELL region is left out.
    std::size_t used = 0;
    for (const auto& claim : claims.values()) {
        p.step();
        if (claim.property) {
            claims.values()[used++] = claim;
        }
    }
    out.regional =
        compare_thickness_regions(p, corresponding, reference, claims.values().first(used));
    out.thickness = out.regional.has_value();
    for (const auto& group : pshell.regions.values().first(pshell.groups)) {
        p.step();
        const auto compared =
            check_pshell_materials(p, source, deck, group, model, shell, case_units);
        add_material_count(p, out, compared);
    }
    return out;
}

} // namespace

core::Result<QualifiedExternalSidecarPg> qualify_external_sidecar_planar_pg(
    const SchemaCase& schema, const ExternalMeshSidecarBinding& bound, MeshShellId shell,
    const structure::PgReferenceGeometry& pg, const SiUnits& case_units,
    const NastranLaminateCards* quantities, const ExternalSidecarPgLimits& limits,
    core::ResourcePlan& plan, const std::stop_token& stop,
    ExternalSidecarPgDiagnostic* diagnosis) noexcept {
    if (diagnosis) {
        *diagnosis = {};
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<QualifiedExternalSidecarPg> {
            try {
                Preparation call{limits, plan, stop, bound.sidecar().source(), diagnosis};
                const auto root = check_case(call, schema, bound);
                const auto source = bound.sidecar().source();
                call.detail.target_namespace = SchemaNamespace::Shell;
                call.detail.target = shell.value();
                const auto model = root.model.record();
                const auto selected = find_shell(call, model, shell);
                const auto pg_id = MeshGeometryId{selected.geometry.unsigned_value()};
                static_cast<void>(find_pg(call, model, pg_id));
                const auto count = selected_faces(call, bound, shell, pg_id);
                auto geometry = establish_geometry(call, bound, pg, shell, count);
                const auto& payload = bound.sidecar().payload();
                QuantityCounts result;
                if (const auto* gmsh = std::get_if<GmshSidecar>(&payload)) {
                    result.regional = gmsh_quantities(call, bound, *gmsh, shell, pg, geometry);
                    result.thickness = result.regional.has_value();
                } else {
                    if (!quantities) {
                        call.reject(Code::Unsupported, Reason::Binding);
                    }
                    const auto original =
                        call.require(std::get<NastranShellCards>(bound.source()).view());
                    result =
                        nastran_quantities(call, bound, std::get<NastranSidecar>(payload), original,
                                           *quantities, model, selected, case_units, pg, geometry);
                }
                call.step();
                return QualifiedExternalSidecarPg{std::move(geometry),
                                                  shell,
                                                  pg_id,
                                                  source,
                                                  bound.sidecar().digest(),
                                                  result.thickness,
                                                  result.materials,
                                                  call.work,
                                                  std::move(result.laminate),
                                                  std::move(result.regional)};
            } catch (const Rejected& error) {
                return std::unexpected(error.error);
            }
        });
}

core::Result<CertifiedExternalSidecarPgBounds> certify_external_sidecar_pg_bounds(
    const SiCase& si, const ExternalMeshSidecarBinding& bound, MeshShellId shell,
    const structure::PgReferenceGeometry& pg, const ExternalSidecarPgLimits& limits,
    std::uint64_t maximum_directed_work, core::ResourcePlan& plan, const std::stop_token& stop,
    ExternalSidecarPgDiagnostic* diagnosis) noexcept {
    if (diagnosis) {
        *diagnosis = {};
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<CertifiedExternalSidecarPgBounds> {
            try {
                const auto source = bound.sidecar().source();
                Preparation call{limits, plan, stop, source, diagnosis};
                const auto root = check_case(call, call.value(si.source_case()).get(), bound);
                call.detail.target_namespace = SchemaNamespace::Shell;
                call.detail.target = shell.value();
                const auto model = root.model.record();
                const auto selected = find_shell(call, model, shell);
                const auto pg_id = MeshGeometryId{selected.geometry.unsigned_value()};
                const auto count = selected_faces(call, bound, shell, pg_id);
                const auto& sidecar = bound.sidecar();
                const auto fit = std::holds_alternative<GmshMesh>(bound.source())
                                     ? find_fit(call, model, sidecar.geometry())
                                     : std::nullopt;
                auto budget = read_conversion_budget(call, si, root, sidecar, fit);
                auto geometry =
                    bound_conversion_geometry(call, bound, pg, shell, count, maximum_directed_work);
                compare_conversion_bounds(call, geometry, budget);
                call.step();
                return CertifiedExternalSidecarPgBounds{std::move(geometry),
                                                        shell,
                                                        pg_id,
                                                        source,
                                                        bound.sidecar().digest(),
                                                        std::move(budget.distance),
                                                        std::move(budget.normal),
                                                        std::move(budget.observations)};
            } catch (const Rejected& error) {
                return std::unexpected(error.error);
            }
        });
}

core::Result<PreparedGmshSidecarPg>
prepare_gmsh_sidecar_pg(const SiCase& si, const ExternalMeshSidecarBinding& bound,
                        MeshShellId shell, const structure::PgReferenceGeometry& pg,
                        const ExternalSidecarPgLimits& limits, std::uint64_t maximum_directed_work,
                        core::ResourcePlan& plan, const std::stop_token& stop,
                        ExternalSidecarPgDiagnostic* diagnosis) noexcept {
    if (diagnosis) {
        *diagnosis = {};
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<PreparedGmshSidecarPg> {
            try {
                const auto& sidecar = bound.sidecar();
                const auto source = sidecar.source();
                Preparation call{limits, plan, stop, source, diagnosis};
                const auto* gmsh = std::get_if<GmshSidecar>(&sidecar.payload());
                if (!gmsh || !std::holds_alternative<GmshMesh>(bound.source())) {
                    call.reject(Code::Unsupported, Reason::Binding);
                }
                const auto root = check_case(call, call.value(si.source_case()).get(), bound);
                const auto model = root.model.record();
                call.detail.target_namespace = SchemaNamespace::Shell;
                call.detail.target = shell.value();
                const auto selected = find_shell(call, model, shell);
                const auto pg_id = MeshGeometryId{selected.geometry.unsigned_value()};
                static_cast<void>(find_pg(call, model, pg_id));
                const auto count = selected_faces(call, bound, shell, pg_id);
                const auto fit = find_fit(call, model, sidecar.geometry());

                auto membership = call.value(
                    plan.allocate<GmshPgFaceMembership>(static_cast<core::Index>(count), operation),
                    Reason::Resource);
                std::size_t i = 0;
                for (const auto& face : bound.faces()) {
                    call.step();
                    if (face.shell == shell) {
                        membership.values()[i++] = {
                            gmsh->shells.values()[face.sidecar_row].physical, face.element,
                            face.sidecar_node};
                    }
                }
                std::optional<ExternalPgThicknessReport> thickness;
                auto geometry = [&]() -> PreparedGmshSidecarPg::Geometry {
                    if (!fit && !sidecar.conversion()) {
                        auto exact = establish_geometry(call, bound, pg, shell, count);
                        thickness = gmsh_quantities(call, bound, *gmsh, shell, pg, exact);
                        return exact;
                    }
                    auto budget = read_conversion_budget(call, si, root, sidecar, fit);
                    auto measured = bound_conversion_geometry(call, bound, pg, shell, count,
                                                              maximum_directed_work);
                    compare_conversion_bounds(call, measured, budget);
                    auto& bounded = std::get<GmshPgBounds>(measured);
                    thickness = gmsh_quantities(call, bound, *gmsh, shell, pg, bounded);
                    return GmshPgConversionEvidence{std::move(bounded), std::move(budget.distance),
                                                    std::move(budget.normal),
                                                    std::move(budget.observations)};
                }();
                call.step();
                return PreparedGmshSidecarPg{std::move(geometry),
                                             shell,
                                             pg_id,
                                             source,
                                             sidecar.digest(),
                                             std::move(membership),
                                             std::move(thickness)};
            } catch (const Rejected& error) {
                return std::unexpected(error.error);
            }
        });
}

} // namespace uvlas::io::detail
