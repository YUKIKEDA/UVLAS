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

struct ForwardStop {
    std::stop_source& source;

    void operator()() const noexcept {
        source.request_stop();
    }
};

class Preparation {
    std::stop_source cancellation_;
    std::stop_callback<ForwardStop> caller_stop_;
    std::stop_callback<ForwardStop> geometry_stop_;
    std::stop_callback<ForwardStop> reference_stop_;

public:
    Preparation(const ExternalSidecarPgLimits& configured, core::ResourcePlan& resources,
                const std::stop_token& token, RootSnapshotLocation location,
                ExternalSidecarPgDiagnostic* output, const std::stop_token& reference = {})
        : caller_stop_(token, ForwardStop{cancellation_}),
          geometry_stop_(configured.geometry.arithmetic.cancellation, ForwardStop{cancellation_}),
          reference_stop_(reference, ForwardStop{cancellation_}),
          limits(configured),
          plan(resources),
          stop(cancellation_.get_token()),
          diagnostic(output) {
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
    const std::stop_token stop;
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

[[nodiscard]] ConversionBudget read_conversion_budget(Preparation& p, const SiCase& si,
                                                      const schema::Case& root,
                                                      const ExternalMeshSidecarBinding& bound) {
    const auto& exact = p.limits.geometry.arithmetic.algebra;
    const auto& sidecar = bound.sidecar();
    const auto& declared = sidecar.conversion();
    const auto fit = std::holds_alternative<GmshMesh>(bound.source())
                         ? find_fit(p, root.model.record(), sidecar.geometry())
                         : std::nullopt;
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

[[nodiscard]] bool gmsh_quantities(Preparation& p, const GmshSidecar& source,
                                   const schema::GeometryPg& pg, MeshShellId shell,
                                   SiLengthUnit case_length) {
    const structure::ExactRational* full = nullptr;
    std::optional<core::Index> earliest_declaration;
    std::optional<GmshKey> earliest_group;
    bool missing = false;
    for (const auto& row : source.shells.values()) {
        p.step();
        if (row.shell != shell) {
            continue;
        }
        p.detail.node = row.node;
        p.detail.physical = row.physical;
        if (const auto& thickness = row.full_thickness) {
            if (full && *full != *thickness) {
                p.reject(Code::Unsupported, Reason::Binding);
            }
            if (!full) {
                full = &*thickness;
                earliest_declaration = row.full_thickness_node;
                earliest_group = row.physical;
            }
        } else {
            missing = true;
        }
    }
    if (!full) {
        return false;
    }
    if (missing) {
        p.reject(Code::Unsupported, Reason::Binding);
    }
    p.detail.node = earliest_declaration;
    p.detail.physical = earliest_group;
    p.check(check_imported_uniform_thickness(
        *full, source.length, pg.thickness.record(), case_length, p.limits.case_thickness_vertices,
        p.limits.geometry.arithmetic.algebra, p.plan, p.stop, mismatch(p)));
    return true;
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
    std::array<std::optional<NastranPcompGeometryMatch>, 2> laminate{};
};

[[nodiscard]] std::size_t same_property(Preparation& p, const ExternalMeshSidecarBinding& bound,
                                        const NastranShellView& deck, MeshShellId target) {
    std::optional<std::size_t> selected;
    for (const auto& face : bound.faces()) {
        p.step();
        if (face.shell != target) {
            continue;
        }
        p.detail.element = face.element;
        p.detail.node = face.sidecar_node;
        const auto found = std::ranges::find_if(deck.elements, [&](const auto& element) {
            p.step();
            return element.card_ordinal == face.source_ordinal;
        });
        if (found == deck.elements.end()) {
            p.reject(Code::InvalidInput, Reason::Binding);
        }
        const auto ordinal = found->property.property_ordinal;
        if (selected && *selected != ordinal) {
            p.reject(Code::Unsupported, Reason::Binding);
        }
        if (std::visit(
                [](const auto& ring) {
                    return std::ranges::any_of(
                        ring, [](const auto& corner) { return corner.thickness.has_value(); });
                },
                found->corners)) {
            p.reject(Code::Unsupported, Reason::Binding);
        }
        selected = ordinal;
    }
    return p.require(selected, Reason::Binding);
}

[[nodiscard]] NastranPcompGeometryMatch
compare_laminate(Preparation& p, const NastranLaminateCards& deck, const NastranSidecar& source,
                 std::size_t property, MeshSidecarOrientation winding,
                 const structure::PgReferenceGeometry& reference) {
    const NastranPcompGeometryLimits limits{p.limits.geometry.arithmetic.algebra,
                                            p.limits.case_thickness_vertices,
                                            p.limits.work - p.work};
    std::optional<NastranPcompGeometryDiagnostic> detail;
    auto match = compare_nastran_pcomp_geometry(deck, property, source.units.length, winding,
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

[[nodiscard]] QuantityCounts
laminate_quantities(Preparation& p, const NastranLaminateCards& cards,
                    const NastranLaminateView& deck, const NastranSidecar& source,
                    std::size_t property, const schema::Model& model,
                    std::span<const BoundExternalFace> original_faces, const schema::Shell& shell,
                    const structure::PgReferenceGeometry& reference,
                    const ExternalPlanarPg& corresponding, const SiUnits& case_units) {
    QuantityCounts out;
    out.thickness = true;
    const auto view = p.require(corresponding.view(), Reason::Binding);
    for (const auto& face : view.faces) {
        p.step();
        const bool reversed = face.winding == ExternalPgWinding::Reversed;
        auto& accepted = out.laminate[reversed ? 1U : 0U];
        if (accepted) {
            continue;
        }
        p.detail.element = face.element;
        const auto declaration = std::ranges::find_if(original_faces, [&](const auto& original) {
            p.step();
            return original.element == face.element;
        });
        if (declaration == original_faces.end()) {
            p.reject(Code::InvalidInput, Reason::Binding);
        }
        p.detail.node = declaration->sidecar_node;
        const auto orientation =
            reversed ? MeshSidecarOrientation::Reversed : MeshSidecarOrientation::Same;
        auto laminate = compare_laminate(p, cards, source, property, orientation, reference);
        for (const auto& link : laminate.layers.values()) {
            p.step();
            const auto layer = p.value(shell.layers.at(static_cast<core::Index>(link.case_ordinal)),
                                       Reason::CaseReference)
                                   .record();
            if (layer.id.unsigned_value() != link.layer.value()) {
                p.reject(Code::StaleGeneration, Reason::CaseReference);
            }
            nastran_material(p, source, deck, link.material, model,
                             MeshMaterialId{layer.material.unsigned_value()}, case_units);
            ++out.materials;
        }
        accepted.emplace(std::move(laminate));
    }
    return out;
}

[[nodiscard]] QuantityCounts nastran_quantities(
    Preparation& p, const ExternalMeshSidecarBinding& bound, const NastranSidecar& source,
    const NastranShellView& original, const NastranLaminateCards& cards, const schema::Model& model,
    const schema::Shell& shell, const schema::GeometryPg& pg, const SiUnits& case_units,
    const structure::PgReferenceGeometry& reference, const ExternalPlanarPg& corresponding) {
    const auto deck = p.require(cards.view());
    if (!original.input.bulk.same_snapshot(deck.input.input.input.bulk.location,
                                           deck.input.input.input.bulk.digest)) {
        p.reject(Code::StaleGeneration, Reason::Binding);
    }
    const auto target = MeshShellId{shell.id.unsigned_value()};
    const auto ordinal = same_property(p, bound, deck.input.input, target);
    if (deck.input.input.properties[ordinal].kind == NastranShellPropertyKind::Pcomp) {
        return laminate_quantities(p, cards, deck, source, ordinal, model, bound.faces(), shell,
                                   reference, corresponding, case_units);
    }
    const auto property = std::ranges::find_if(deck.input.pshell, [&](const auto& entry) {
        p.step();
        return entry.property_ordinal == ordinal;
    });
    if (property == deck.input.pshell.end()) {
        p.reject(Code::Unsupported, Reason::Binding);
    }
    QuantityCounts out;
    if (const auto& thickness = property->thickness) {
        p.detail.byte = thickness->field.position.byte;
        p.check(check_imported_uniform_thickness(
            {thickness->value.decimal(), source.units.length, thickness->field.position},
            pg.thickness.record(), case_units.length, p.limits.case_thickness_vertices,
            p.limits.geometry.arithmetic.algebra, p.plan, p.stop, mismatch(p)));
        out.thickness = true;
    }
    if (std::holds_alternative<NastranBendingDesignation>(property->bending)) {
        p.reject(Code::Unsupported, Reason::Binding);
    }
    const auto* bending = std::get_if<NastranMaterialReference>(&property->bending);
    const std::array<const NastranMaterialReference*, 4> uses{
        property->membrane ? &*property->membrane : nullptr, bending,
        property->transverse_shear ? &*property->transverse_shear : nullptr,
        property->coupling ? &*property->coupling : nullptr};
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
        ++out.materials;
    }
    return out;
}

} // namespace

core::Result<QualifiedCaseSidecarPg> qualify_external_sidecar_planar_pg(
    const SiCase& si, const ExternalMeshSidecarBinding& bound, MeshShellId shell,
    const NastranLaminateCards* quantities, const CasePgReferenceLimits& reference_limits,
    const ExternalSidecarPgLimits& limits, core::Workspace& workspace, core::ResourcePlan& plan,
    const std::stop_token& stop, ExternalSidecarPgDiagnostic* diagnosis) noexcept {
    if (diagnosis) {
        *diagnosis = {};
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<QualifiedCaseSidecarPg> {
            try {
                Preparation call{limits,    plan,
                                 stop,      bound.sidecar().source(),
                                 diagnosis, reference_limits.cancellation};
                const auto root = check_case(call, call.value(si.source_case()).get(), bound);
                const auto source = bound.sidecar().source();
                call.detail.target_namespace = SchemaNamespace::Shell;
                call.detail.target = shell.value();
                const auto units = root.units.record();
                const SiUnits case_units{static_cast<SiLengthUnit>(units.length.selected()),
                                         static_cast<SiMassUnit>(units.mass.selected()),
                                         static_cast<SiTimeUnit>(units.time.selected()),
                                         static_cast<SiAngleUnit>(units.angle.selected())};
                const auto model = root.model.record();
                const auto selected = find_shell(call, model, shell);
                const auto pg_id = MeshGeometryId{selected.geometry.unsigned_value()};
                const auto pg_record = find_pg(call, model, pg_id);
                const auto count = selected_faces(call, bound, shell, pg_id);
                auto preparation = reference_limits;
                preparation.cancellation = call.stop;
                auto pg = call.value(
                    prepare_case_pg_reference(si, shell.value(), preparation, workspace, plan),
                    Reason::CaseReference);
                auto geometry = establish_geometry(call, bound, pg, shell, count);
                const auto& payload = bound.sidecar().payload();
                QuantityCounts result;
                if (const auto* gmsh = std::get_if<GmshSidecar>(&payload)) {
                    result.thickness =
                        gmsh_quantities(call, *gmsh, pg_record, shell, case_units.length);
                } else {
                    if (!quantities) {
                        call.reject(Code::Unsupported, Reason::Binding);
                    }
                    const auto original =
                        call.require(std::get<NastranShellCards>(bound.source()).view());
                    result = nastran_quantities(call, bound, std::get<NastranSidecar>(payload),
                                                original, *quantities, model, selected, pg_record,
                                                case_units, pg, geometry);
                }
                call.step();
                return QualifiedCaseSidecarPg{
                    std::move(pg),
                    {std::move(geometry), shell, pg_id, source, bound.sidecar().digest(),
                     result.thickness, result.materials, call.work, std::move(result.laminate)}};
            } catch (const Rejected& error) {
                return std::unexpected(error.error);
            }
        });
}

core::Result<CertifiedCaseSidecarPgBounds> certify_external_sidecar_pg_bounds(
    const SiCase& si, const ExternalMeshSidecarBinding& bound, MeshShellId shell,
    const CasePgReferenceLimits& reference_limits, const ExternalSidecarPgLimits& limits,
    std::uint64_t maximum_directed_work, core::Workspace& workspace, core::ResourcePlan& plan,
    const std::stop_token& stop, ExternalSidecarPgDiagnostic* diagnosis) noexcept {
    if (diagnosis) {
        *diagnosis = {};
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<CertifiedCaseSidecarPgBounds> {
            try {
                const auto source = bound.sidecar().source();
                Preparation call{limits, plan,      stop,
                                 source, diagnosis, reference_limits.cancellation};
                const auto root = check_case(call, call.value(si.source_case()).get(), bound);
                call.detail.target_namespace = SchemaNamespace::Shell;
                call.detail.target = shell.value();
                const auto model = root.model.record();
                const auto selected = find_shell(call, model, shell);
                const auto pg_id = MeshGeometryId{selected.geometry.unsigned_value()};
                const auto count = selected_faces(call, bound, shell, pg_id);
                auto budget = read_conversion_budget(call, si, root, bound);
                auto preparation = reference_limits;
                preparation.cancellation = call.stop;
                auto pg = call.value(
                    prepare_case_pg_reference(si, shell.value(), preparation, workspace, plan),
                    Reason::CaseReference);
                auto geometry =
                    bound_conversion_geometry(call, bound, pg, shell, count, maximum_directed_work);
                compare_conversion_bounds(call, geometry, budget);
                call.step();
                return CertifiedCaseSidecarPgBounds{
                    std::move(pg),
                    {std::move(geometry), shell, pg_id, source, bound.sidecar().digest(),
                     std::move(budget.distance), std::move(budget.normal),
                     std::move(budget.observations)}};
            } catch (const Rejected& error) {
                return std::unexpected(error.error);
            }
        });
}

} // namespace uvlas::io::detail
