#include "case_material_surface_model.hpp"

#include "case_material_alias.hpp"
#include "structure/pg_reference_arithmetic.hpp"

#include <algorithm>
#include <tuple>

namespace uvlas::coupling::detail {
namespace {
namespace a = aero::detail;
namespace pg = structure::detail;
namespace s = structure;
using core::FailureCode;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_mul;
using pg::partition_take;

constexpr auto preparation = core::Operation::PrepareModel;

struct AssemblyWork {
    std::uint64_t remaining{};
    std::stop_token stop;

    void visit(std::uint64_t used = 1) {
        if (stop.stop_requested())
            partition_fail(FailureCode::Cancelled);
        if (used > remaining)
            partition_fail(FailureCode::ResourceExceeded);
        remaining -= used;
    }
};

template <class T, class E>
T aero_result(std::expected<T, E> result) {
    if (!result) {
        if constexpr (requires { result.error().cause; })
            throw pg::PartitionAbort{result.error().cause};
        else
            throw pg::PartitionAbort{result.error().origin.cause};
    }
    return std::move(*result);
}

struct SurfaceInput {
    std::uint64_t id{};
    std::span<const a::PhysicalFaceId> faces{};
};

struct FaceInput {
    a::PhysicalFaceId id{0};
    std::size_t surface{};
    std::optional<a::UvlmOriginalQuad> quad;
};

bool same_case(const io::SchemaIdentity& a, const io::SchemaIdentity& b) noexcept {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}

template <class Buffer, class Key, class Projection>
auto lookup(Buffer rows, Key key, Projection project, AssemblyWork& work) {
    return std::ranges::lower_bound(
        rows, key,
        [&](auto a, auto b) {
            work.visit();
            return a < b;
        },
        project);
}

auto corner_key(const CaseMaterialOriginalCorner& row) noexcept {
    return std::pair{row.facet, row.local_vertex};
}

auto site_key(const CaseMaterialVertexSite& site) noexcept {
    return std::tuple{site.chord_column, site.span_knot, site.original_span_cell};
}

} // namespace

struct CaseMaterialSurfaceModel::Surface {
    CaseMaterialChart chart;
    CaseReferenceSurface thin;
    CaseMaterialOriginalLayout layout;
    CaseMaterialChartIdentity id;
    s::PgPreparedPatchIdentity patch{};
    std::size_t patch_ordinal{};
};

struct CaseMaterialSurfaceModel::Route {
    CaseMaterialOriginalCorner corner;
    std::size_t surface{};
};

CaseMaterialSurfaceModel::CaseMaterialSurfaceModel(
    core::PlannedBuffer<std::optional<Surface>> surfaces,
    core::PlannedBuffer<Route> routes) noexcept
    : surfaces_(std::move(surfaces)), routes_(std::move(routes)) {}

CaseMaterialSurfaceModel::CaseMaterialSurfaceModel(CaseMaterialSurfaceModel&&) noexcept = default;
CaseMaterialSurfaceModel&
CaseMaterialSurfaceModel::operator=(CaseMaterialSurfaceModel&&) noexcept = default;
CaseMaterialSurfaceModel::~CaseMaterialSurfaceModel() = default;

core::Result<CaseMaterialSurfaceModel> CaseMaterialSurfaceModel::prepare(
    const io::SiCase& case_input, std::span<const CaseAeroSurfacePanels> assembly,
    const a::GlobalContactState& contact, const CaseStructuralSupport& support,
    std::span<const s::PgConstraintPatchBinding> patches,
    const CaseMaterialSurfaceModelLimits& limits, core::Workspace& workspace,
    core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        preparation, FailureCode::PrecisionUnresolved, {},
        [&]() -> core::Result<CaseMaterialSurfaceModel> {
            try {
                AssemblyWork work{limits.maximum_assembly_steps, limits.alias.cancellation};
                work.visit();
                const auto& source = partition_take(case_input.source_case()).get();
                if (!same_case(partition_take(source.identity()), support.source()))
                    partition_fail(FailureCode::StaleGeneration);
                const auto declared = partition_take(source.root()).model.record().aero_surfaces;
                if (declared.size() != assembly.size() || assembly.empty())
                    partition_fail(FailureCode::InvalidInput);
                if (assembly.size() > limits.maximum_surfaces)
                    partition_fail(FailureCode::ResourceExceeded);
                auto records = partition_take(
                    plan.allocate<SurfaceInput>(partition_index(assembly.size()), preparation));
                std::uint64_t number = 0;
                for (std::size_t i = 0; i < assembly.size(); ++i) {
                    work.visit();
                    records.values()[i] = {assembly[i].surface, assembly[i].physical_panels};
                    number = partition_add(number, assembly[i].physical_panels.size());
                }
                if (number > limits.maximum_bound_panels)
                    partition_fail(FailureCode::ResourceExceeded);
                std::ranges::sort(records.values(), [&](const auto& a, const auto& b) {
                    work.visit();
                    return a.id < b.id;
                });
                auto faces =
                    partition_take(plan.allocate<FaceInput>(partition_index(number), preparation));
                std::size_t offset = 0;
                for (std::size_t i = 0; i < records.size(); ++i) {
                    work.visit();
                    const auto& row = records.values()[i];
                    if (i != 0 && records.values()[i - 1].id == row.id)
                        partition_fail(FailureCode::InvalidInput);
                    for (const auto face : row.faces) {
                        work.visit();
                        faces.values()[offset++] = {face, i, {}};
                    }
                }
                std::ranges::sort(faces.values(), [&](const auto& a, const auto& b) {
                    work.visit();
                    return a.id < b.id;
                });
                for (std::size_t i = 1; i < faces.size(); ++i) {
                    work.visit();
                    if (faces.values()[i - 1].id == faces.values()[i].id)
                        partition_fail(FailureCode::InvalidInput);
                }
                const auto geometry = aero_result(contact.geometry());
                const auto counts = aero_result(geometry.get().counts());
                const auto face_book = faces.values();
                std::uint64_t bound = 0;
                for (std::size_t i = 0; i < counts.faces; ++i) {
                    work.visit();
                    const auto role = aero_result(contact.face_role(i));
                    if (role.role != a::AeroFaceRole::Bound)
                        continue;
                    const auto entry = lookup(face_book, role.face, &FaceInput::id, work);
                    if (entry == face_book.end() || entry->id != role.face)
                        partition_fail(FailureCode::GeometryInvalid);
                    ++bound;
                }
                if (bound != faces.size())
                    partition_fail(FailureCode::GeometryInvalid);
                for (std::size_t i = 0; i < counts.quads; ++i) {
                    const auto quad =
                        aero_result(contact.uvlm_original_quad(i, work.remaining, work.stop));
                    work.visit(quad.used_steps);
                    if (quad.role != a::AeroFaceRole::Bound)
                        continue;
                    const auto entry = lookup(face_book, quad.face, &FaceInput::id, work);
                    if (entry == face_book.end() || entry->id != quad.face || entry->quad)
                        partition_fail(FailureCode::GeometryInvalid);
                    entry->quad = quad;
                }
                for (std::size_t i = 0; i < counts.triangles; ++i) {
                    work.visit();
                    const auto face = aero_result(contact.original_triangle_face(i));
                    const auto entry = lookup(face_book, face, &FaceInput::id, work);
                    if (entry != face_book.end() && entry->id == face)
                        partition_fail(FailureCode::GeometryInvalid);
                }
                for (const auto& entry : face_book) {
                    work.visit();
                    if (!entry.quad)
                        partition_fail(FailureCode::GeometryInvalid);
                }

                auto surfaces = partition_take(plan.allocate<std::optional<Surface>>(
                    partition_index(assembly.size()), preparation));
                const auto declared_book = records.values();
                const auto active = support.support().chart().original().layout();
                auto material = limits.material;
                material.chart.cancellation = work.stop;
                material.support.cancellation = work.stop;
                for (std::uint64_t i = 0; i < declared.size(); ++i) {
                    work.visit();
                    const auto surface = partition_take(declared.at(partition_index(i))).record();
                    const auto id = surface.id.unsigned_value();
                    const auto row = lookup(declared_book, id, &SurfaceInput::id, work);
                    if (row == declared_book.end() || row->id != id)
                        partition_fail(FailureCode::InvalidInput);
                    if (surface.chord.size() < 2 || surface.span.size() < 2 ||
                        row->faces.size() !=
                            partition_mul(surface.chord.size() - 1, surface.span.size() - 1))
                        partition_fail(FailureCode::InvalidInput);
                    const auto attachment = surface.attachment.record();
                    const auto* structural =
                        std::get_if<io::schema::AeroAttachmentStructural>(&attachment);
                    if (structural == nullptr)
                        continue;
                    const auto host = structural->host.unsigned_value();
                    const auto binding = std::ranges::find_if(patches, [&](const auto& candidate) {
                        work.visit();
                        return candidate.id.value() == host;
                    });
                    const auto structural_layout =
                        std::ranges::find_if(active, [&](const auto& candidate) {
                            work.visit();
                            return candidate.id.value() == host;
                        });
                    if (binding == patches.end() || structural_layout == active.end())
                        partition_fail(FailureCode::StaleGeneration);
                    const auto patch = partition_take(binding->patch.info()).identity;
                    if (patch != structural_layout->basis)
                        partition_fail(FailureCode::StaleGeneration);
                    auto chart = partition_take(
                        CaseMaterialChart::prepare(case_input, id, material, workspace, plan));
                    const auto identity = partition_take(chart.identity());
                    auto thin = partition_take(CaseReferenceSurface::prepare(
                        case_input, identity.geometry_id, material.chart, plan));
                    auto layout = partition_take(CaseMaterialOriginalLayout::prepare(
                        chart, thin, contact, row->faces, material.chart, plan));
                    const auto corners = layout.original_corners();
                    if (corners.empty())
                        partition_fail(FailureCode::GeometryInvalid);
                    // A surface without aliases still has to bind this actual
                    // conditioned patch to the Case reference.
                    static_cast<void>(partition_take(CaseMaterialAliasRow::prepare(
                        chart, thin, binding->patch, corners.front().site, limits.alias, plan)));
                    const auto ordinal = static_cast<std::size_t>(row - declared_book.begin());
                    const auto patch_ordinal = static_cast<std::size_t>(binding - patches.begin());
                    surfaces.values()[ordinal].emplace(Surface{std::move(chart), std::move(thin),
                                                               std::move(layout), identity, patch,
                                                               patch_ordinal});
                }
                auto routes = partition_take(
                    plan.allocate<Route>(partition_index(partition_mul(number, 6)), preparation));
                std::size_t route_count = 0;
                for (std::size_t i = 0; i < surfaces.size(); ++i) {
                    work.visit();
                    const auto& owner = surfaces.values()[i];
                    if (!owner)
                        continue;
                    for (const auto& corner : owner->layout.original_corners()) {
                        work.visit();
                        routes.values()[route_count++] = {corner, i};
                    }
                }
                if (route_count == 0)
                    partition_fail(FailureCode::InvalidInput);
                for (const auto& face : face_book) {
                    work.visit();
                    if (surfaces.values()[face.surface])
                        continue;
                    if (!face.quad)
                        partition_fail(FailureCode::GeometryInvalid);
                    for (const auto id : face.quad->facets) {
                        const auto points = aero_result(geometry.get().facet_points(id));
                        for (std::size_t local = 0; local < points.size(); ++local) {
                            work.visit();
                            const auto point = aero_result(geometry.get().owner(points[local]));
                            routes.values()[route_count++] = {
                                {face.id, id, static_cast<std::uint8_t>(local), point, {}},
                                face.surface};
                        }
                    }
                }
                if (route_count != routes.size())
                    partition_fail(FailureCode::GeometryInvalid);
                const auto items = routes.values().first(route_count);
                std::ranges::sort(items, [&](const Route& a, const Route& b) {
                    work.visit();
                    return std::tuple{a.corner.canonical, a.surface, site_key(a.corner.site)} <
                           std::tuple{b.corner.canonical, b.surface, site_key(b.corner.site)};
                });
                const auto original_point = [&](const Route& route) {
                    const auto& item = surfaces.values()[route.surface];
                    if (!item)
                        partition_fail(FailureCode::Unsupported);
                    return CaseMaterialOriginalAliasPoint{
                        item->layout, route.corner.facet, route.corner.local_vertex,
                        item->chart,  item->thin,         patches[item->patch_ordinal].patch};
                };
                std::uint64_t comparisons = 0;
                for (std::size_t first = 0; first < items.size();) {
                    work.visit();
                    auto index = first + 1;
                    while (index < items.size() &&
                           items[index].corner.canonical == items[first].corner.canonical) {
                        work.visit();
                        const auto& here = items[index];
                        const auto& previous = items[index - 1];
                        const bool distinct =
                            here.surface != previous.surface ||
                            site_key(here.corner.site) != site_key(previous.corner.site);
                        if (distinct && (surfaces.values()[previous.surface] ||
                                         surfaces.values()[here.surface])) {
                            if (comparisons == limits.maximum_alias_checks)
                                partition_fail(FailureCode::ResourceExceeded);
                            ++comparisons;
                            partition_take(qualify_case_material_original_alias(
                                original_point(previous), original_point(here), support,
                                limits.alias, workspace, plan));
                        }
                        ++index;
                    }
                    first = index;
                }
                std::ranges::sort(items, [&](const Route& a, const Route& b) {
                    work.visit();
                    return corner_key(a.corner) < corner_key(b.corner);
                });
                work.visit();
                return CaseMaterialSurfaceModel{std::move(surfaces), std::move(routes)};
            } catch (const pg::PartitionAbort& error) {
                return std::unexpected(
                    core::Failure{error.failure.code, preparation, error.failure.details});
            }
        });
}

std::size_t CaseMaterialSurfaceModel::surface_count() const noexcept {
    return surfaces_.size();
}

const CaseMaterialSurfaceModel::Route*
CaseMaterialSurfaceModel::find(a::FacetId facet, std::uint8_t local) const noexcept {
    const auto items = routes_.values();
    const auto key = std::pair{facet, local};
    const auto it = std::ranges::lower_bound(items, key, {},
                                             [](const Route& x) { return corner_key(x.corner); });
    return it != items.end() && corner_key(it->corner) == key ? &*it : nullptr;
}

core::Result<CaseMaterialSurfaceCorner>
CaseMaterialSurfaceModel::corner(a::FacetId facet, std::uint8_t local) const noexcept {
    const auto* route = find(facet, local);
    if (route == nullptr || route->surface >= surfaces_.size())
        return std::unexpected(core::Failure{FailureCode::InvalidInput, preparation});
    const auto& owner = surfaces_.values()[route->surface];
    if (!owner)
        return std::unexpected(core::Failure{FailureCode::Unsupported, preparation});
    return CaseMaterialSurfaceCorner{owner->id, route->corner};
}

core::Result<CaseMaterialOriginalEnclosedForce> CaseMaterialSurfaceModel::pullback(
    const DiscreteCirculationLoads& loads, std::size_t index, const s::PgPreparedPatch& patch,
    const s::PgPreparedReferenceProjection& projection, const s::PgReferenceLimits& prep,
    const s::PgVolumePositionLimits& apply, core::ResourcePlan& plan) const noexcept {
    const auto load = loads.original(index);
    if (!load)
        return std::unexpected(load.error());
    const auto* route = find(load->facet, load->local_vertex);
    constexpr auto operation = core::Operation::EvaluateStep;
    if (!route || route->surface >= surfaces_.size())
        return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
    const auto& owner = surfaces_.values()[route->surface];
    if (!owner)
        return std::unexpected(core::Failure{FailureCode::Unsupported, operation});
    const auto identity = patch.info();
    if (!identity || identity->identity != owner->patch)
        return std::unexpected(core::Failure{FailureCode::StaleGeneration, operation});
    return owner->layout.pullback(loads, index, owner->chart, owner->thin, patch, projection, prep,
                                  apply, plan);
}

} // namespace uvlas::coupling::detail
