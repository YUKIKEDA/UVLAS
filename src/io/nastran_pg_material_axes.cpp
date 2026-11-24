#include "nastran_pg_material_axes.hpp"

#include "../structure/pg_rank_storage.hpp"
#include "../structure/pg_reference_arithmetic.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

namespace geom = structure::detail;
using geom::partition_fail;
using geom::partition_take;
using geom::ReferenceArithmetic;
using Code = core::FailureCode;
using Reason = NastranPgAxisFailure;
using Exact = structure::ExactRational;
using Ray = std::array<Exact, 3>;
constexpr auto operation = core::Operation::PrepareModel;

struct Call {
    ReferenceArithmetic math;
    std::optional<NastranPgAxisDiagnostic>* output{};
    NastranPgAxisDiagnostic detail{};

    [[noreturn]] void fail(Code code, Reason reason) const {
        if (output) {
            auto failure = detail;
            failure.reason = reason;
            *output = failure;
        }
        partition_fail(code);
    }
};

[[nodiscard]] Exact read_decimal(ReferenceArithmetic& math, std::string_view token) {
    math.charge();
    return partition_take(Exact::from_decimal(token, math.leaf_limits, math.plan));
}

[[nodiscard]] Ray subtract(ReferenceArithmetic& math, const Ray& a, const Ray& b) {
    return {math.sub(a[0], b[0]), math.sub(a[1], b[1]), math.sub(a[2], b[2])};
}

[[nodiscard]] Ray cross(ReferenceArithmetic& math, const Ray& a, const Ray& b) {
    return {math.sub(math.mul(a[1], b[2]), math.mul(a[2], b[1])),
            math.sub(math.mul(a[2], b[0]), math.mul(a[0], b[2])),
            math.sub(math.mul(a[0], b[1]), math.mul(a[1], b[0]))};
}

[[nodiscard]] Exact dot(ReferenceArithmetic& math, const Ray& a, const Ray& b) {
    return math.add(math.add(math.mul(a[0], b[0]), math.mul(a[1], b[1])), math.mul(a[2], b[2]));
}

[[nodiscard]] Ray copy_ray(ReferenceArithmetic& math, const Ray& a) {
    return {math.copy(a[0]), math.copy(a[1]), math.copy(a[2])};
}

[[nodiscard]] bool zero_ray(const Ray& ray) noexcept {
    return std::ranges::all_of(ray, [](const Exact& v) { return v.sign() == 0; });
}

[[nodiscard]] bool same_ray(ReferenceArithmetic& math, const Ray& a, const Ray& b) {
    return zero_ray(cross(math, a, b)) && dot(math, a, b).sign() > 0;
}

[[nodiscard]] Ray project(ReferenceArithmetic& math, const Ray& direction, const Ray& normal) {
    const auto square = dot(math, normal, normal);
    const auto axial = dot(math, direction, normal);
    return {math.sub(math.mul(square, direction[0]), math.mul(axial, normal[0])),
            math.sub(math.mul(square, direction[1]), math.mul(axial, normal[1])),
            math.sub(math.mul(square, direction[2]), math.mul(axial, normal[2]))};
}

[[nodiscard]] Ray quarter_rotate(ReferenceArithmetic& math, const Ray& tangent, const Ray& normal,
                                 unsigned quarters) {
    if (quarters == 0U)
        return copy_ray(math, tangent);
    if (quarters == 1U)
        return cross(math, normal, tangent);
    if (quarters == 2U)
        return subtract(math, Ray{}, tangent);
    return cross(math, tangent, normal);
}

[[nodiscard]] unsigned quarter_turns(Call& call, const Exact& degrees, const Exact& ninety) {
    const auto turns = call.math.div(degrees, ninety);
    if (turns.denominator_bits() != 1)
        call.fail(Code::Unsupported, Reason::AngleRepresentation);
    constexpr std::uint64_t maximum_rank_bits = 8191;
    if (turns.numerator_bits() > maximum_rank_bits)
        partition_fail(Code::ResourceExceeded);
    constexpr auto limb_bits = std::numeric_limits<geom::RankWord>::digits;
    call.math.charge(1 + turns.numerator_bits() / limb_bits);
    const auto fraction = geom::rank_fraction<geom::RankWide>(turns);
    const geom::RankWide residue = fraction.numerator % 4;
    const geom::RankWide positive = residue < 0 ? residue + 4 : residue;
    return positive.convert_to<unsigned>();
}

[[nodiscard]] const schema::ReferenceFieldConstant& require_constant(
    Call& call, const std::variant<schema::ReferenceFieldConstant, schema::ReferenceFieldPolynomial,
                                   schema::ReferenceFieldExpression>& source) {
    if (const auto* constant = std::get_if<schema::ReferenceFieldConstant>(&source)) {
        return *constant;
    }
    call.fail(Code::Unsupported, Reason::CaseAxis);
}

[[nodiscard]] Exact read_case_component(Call& call, const SiCase& si,
                                        const schema::ReferenceFieldConstant& values, core::Index i,
                                        bool angle = false) {
    call.math.charge();
    auto recipe = partition_take(si.quantity(partition_take(values.value.at(i)))).recipe;
    if (angle && recipe.pi_power == 1 && recipe.units.angle == SiAngleUnit::Degree) {
        recipe.pi_power = 0;
        recipe.units.angle = SiAngleUnit::Radian;
        return partition_take(materialize_rational(recipe, call.math.leaf_limits, call.math.plan));
    }
    auto value =
        partition_take(materialize_rational(recipe, call.math.leaf_limits, call.math.plan));
    if (angle && value.sign() != 0) {
        call.fail(Code::Unsupported, Reason::AngleRepresentation);
    }
    return value;
}

[[nodiscard]] Ray case_direction(Call& call, const SiCase& si,
                                 const schema::ReferenceFieldConstant& values) {
    return {read_case_component(call, si, values, 0), read_case_component(call, si, values, 1),
            read_case_component(call, si, values, 2)};
}

[[nodiscard]] bool right_handed_basis(ReferenceArithmetic& math, const std::array<Ray, 3>& basis,
                                      const Exact& one) {
    for (std::size_t i = 0; i < basis.size(); ++i) {
        if (dot(math, basis[i], basis[i]) != one ||
            dot(math, basis[i], basis[(i + 1) % basis.size()]).sign() != 0) {
            return false;
        }
    }
    return dot(math, cross(math, basis[0], basis[1]), basis[2]).sign() > 0;
}

void check_case_axis(Call& call, const SiCase& si, const schema::Layer& layer, const Ray& tangent,
                     const Ray& normal, const Exact& degrees, const Exact& one,
                     const Exact& ninety) {
    const auto axis = layer.axis.record();
    if (const auto* projected = std::get_if<schema::AxisProjected>(&axis)) {
        const auto direction = projected->direction.record();
        const auto angle = projected->angle.record();
        const auto ray =
            project(call.math, case_direction(call, si, require_constant(call, direction)), normal);
        if (zero_ray(ray))
            call.fail(Code::InvalidInput, Reason::CaseAxis);
        const auto case_angle =
            read_case_component(call, si, require_constant(call, angle), 0, true);
        // Remove the Case rotation on both sides before comparing the rays.
        const auto difference = call.math.sub(degrees, case_angle);
        if (!same_ray(
                call.math,
                quarter_rotate(call.math, tangent, normal, quarter_turns(call, difference, ninety)),
                ray))
            call.fail(Code::InvalidInput, Reason::Difference);
        return;
    }
    const auto& explicit_axis = std::get<schema::AxisExplicit>(axis);
    const auto columns = explicit_axis.columns.record();
    const auto& values = require_constant(call, columns);
    std::array<Ray, 3> basis;
    for (core::Index i = 0; i < 3; ++i) {
        basis[static_cast<std::size_t>(i)] = {read_case_component(call, si, values, i),
                                              read_case_component(call, si, values, 3 + i),
                                              read_case_component(call, si, values, 6 + i)};
    }
    if (!right_handed_basis(call.math, basis, one))
        call.fail(Code::InvalidInput, Reason::CaseAxis);
    if (!same_ray(call.math, normal, basis[2]) ||
        !same_ray(call.math,
                  quarter_rotate(call.math, tangent, normal, quarter_turns(call, degrees, ninety)),
                  basis[0]))
        call.fail(Code::InvalidInput, Reason::Difference);
}

[[nodiscard]] Ray rectangular_x(Call& call, const NastranCoordinateView& deck,
                                const NastranShellMaterialAxis& axis, const Exact& one) {
    if (axis.id.value() == 0)
        return {call.math.copy(one), {}, {}};
    const auto ordinal = axis.coordinate_ordinal;
    if (!ordinal || *ordinal >= deck.coordinates.size())
        call.fail(Code::StaleGeneration, Reason::Source);
    const auto& frame = deck.coordinates[*ordinal];
    if (frame.id != axis.id)
        call.fail(Code::StaleGeneration, Reason::Source);
    if (frame.kind != NastranCoordinateKind::Rectangular || frame.reference.id.value() != 0)
        call.fail(Code::Unsupported, Reason::SourceAxis);
    const auto tuple = [&](const std::array<NastranCoordinateValue, 3>& input) {
        return Ray{read_decimal(call.math, input[0].decimal()),
                   read_decimal(call.math, input[1].decimal()),
                   read_decimal(call.math, input[2].decimal())};
    };
    const auto origin = tuple(frame.a);
    const auto z = subtract(call.math, tuple(frame.b), origin);
    const auto in_plane = subtract(call.math, tuple(frame.c), origin);
    // Cross products carry the exact positive ray without square-root normalization.
    return cross(call.math, cross(call.math, z, in_plane), z);
}

struct ElementAxis {
    Ray tangent;
    Ray normal;
    Exact degrees;
};

[[nodiscard]] ElementAxis source_axis(Call& call, const NastranShellElement& element,
                                      const NastranCoordinateView& deck,
                                      const ExternalPgOwnedFace& face, const Exact& one) {
    const auto ring = face.ring();
    auto first = subtract(call.math, ring[1].position_metres, ring[0].position_metres);
    const auto other = subtract(call.math, ring[2].position_metres, ring[0].position_metres);
    auto normal = face.winding == ExternalPgWinding::Forward ? cross(call.math, first, other)
                                                             : cross(call.math, other, first);
    if (const auto* axis = std::get_if<NastranShellMaterialAxis>(&element.orientation)) {
        call.detail.element_field = axis->field.position;
        auto tangent = project(call.math, rectangular_x(call, deck, *axis, one), normal);
        if (zero_ray(tangent))
            call.fail(Code::InvalidInput, Reason::SourceAxis);
        return {std::move(tangent), std::move(normal), {}};
    }
    Exact angle;
    if (const auto* theta = std::get_if<NastranShellReal>(&element.orientation)) {
        call.detail.element_field = theta->field.position;
        angle = read_decimal(call.math, theta->value.decimal());
    }
    return {std::move(first), std::move(normal), std::move(angle)};
}

[[nodiscard]] bool same_parent(const SchemaIdentity& a, const SchemaIdentity& b) noexcept {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}

[[nodiscard]] schema::Shell find_shell(Call& call, const SchemaCase& input, MeshShellId target) {
    const auto root = partition_take(input.root());
    const auto shells = root.model.record().shells;
    for (core::Index i = 0; std::cmp_less(i, shells.size()); ++i) {
        call.math.visit(1);
        auto shell = partition_take(shells.at(i)).record();
        if (shell.id.unsigned_value() == target.value())
            return shell;
    }
    call.fail(Code::StaleGeneration, Reason::Source);
}

void check_constant_thickness(Call& call, const QualifiedExternalSidecarPg& admitted,
                              const ExternalPgOwnedFace& face, NastranPropertyId property,
                              const structure::PgReferenceIdentity& pg,
                              const NastranBulkView& deck) {
    const auto& report = admitted.thickness_regions;
    if (!report)
        call.fail(Code::Unsupported, Reason::Thickness);
    if (report->pg != pg || !deck.same_snapshot(report->source, report->digest))
        call.fail(Code::StaleGeneration, Reason::Source);
    const auto matches = report->matches.values();
    const auto found = std::ranges::find_if(matches, [&](const auto& row) {
        call.math.visit(1);
        return row.element == face.element;
    });
    if (found == matches.end() || found->property != property || found->winding != face.winding ||
        !std::holds_alternative<Exact>(found->value_metres))
        call.fail(Code::Unsupported, Reason::Thickness);
}

[[nodiscard]] const NastranPcompLayerSection*
find_laminate(Call& call, const QualifiedExternalSidecarPg& admitted,
              const NastranShellElement& element, ExternalPgWinding winding,
              const structure::PgReferenceIdentity& pg, const NastranBulkView& deck) {
    const auto orientation = winding == ExternalPgWinding::Forward
                                 ? MeshSidecarOrientation::Same
                                 : MeshSidecarOrientation::Reversed;
    for (const auto& section : admitted.laminate.values()) {
        call.math.visit(1);
        if (section && section->property == element.property.id &&
            section->winding == orientation) {
            if (section->pg != pg || !deck.same_snapshot(section->source, section->digest))
                call.fail(Code::StaleGeneration, Reason::Source);
            return &*section;
        }
    }
    return nullptr;
}

} // namespace

core::Result<NastranPgMaterialAxes>
match_nastran_pg_material_axes(const SiCase& si, const ExternalMeshSidecarBinding& bound,
                               const QualifiedExternalSidecarPg& admitted,
                               const NastranPgAxisLimits& limits, core::ResourcePlan& plan,
                               std::optional<NastranPgAxisDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranPgMaterialAxes> {
            try {
                Call call{{limits.arithmetic, plan}, diagnostic};
                call.math.cancel();
                const auto& schema = partition_take(si.source_case()).get();
                const auto identity = partition_take(schema.identity());
                const auto& sidecar = bound.sidecar();
                const auto sidecar_location = sidecar.source();
                const auto* native = std::get_if<NastranShellCards>(&bound.source());
                const auto geometry = admitted.geometry.view();
                if (!native)
                    call.fail(Code::Unsupported, Reason::Source);
                if (!geometry || !same_parent(identity, sidecar_location.parent) ||
                    !same_parent(identity, admitted.sidecar.parent) ||
                    admitted.sidecar.source != sidecar_location.source ||
                    admitted.sidecar.role != sidecar_location.role ||
                    admitted.sidecar.declaration_node != sidecar_location.declaration_node ||
                    admitted.sidecar.declaring_file != sidecar_location.declaring_file ||
                    admitted.sidecar_digest != sidecar.digest())
                    call.fail(Code::StaleGeneration, Reason::Source);
                const auto deck = native->view();
                if (!deck || !deck->input.bulk.same_snapshot(geometry->source, geometry->digest))
                    call.fail(Code::StaleGeneration, Reason::Source);
                const auto shell = find_shell(call, schema, admitted.shell);
                if (shell.geometry.unsigned_value() != admitted.pg.value())
                    call.fail(Code::StaleGeneration, Reason::Source);
                const auto layers = shell.layers.size();
                const auto faces = geometry->faces;
                if (layers == 0 || faces.empty())
                    call.fail(Code::InvalidInput, Reason::Source);
                if (faces.size() > limits.maximum_matches / layers)
                    partition_fail(Code::ResourceExceeded);
                const auto count = geom::partition_index(geom::partition_mul(faces.size(), layers));
                call.math.admit();
                auto output = partition_take(plan.allocate<NastranPgLayerAxis>(count, operation));
                const auto one = read_decimal(call.math, "1");
                const auto ninety = read_decimal(call.math, "90");
                std::size_t next = 0;
                for (const auto& face : faces) {
                    call.math.visit(1);
                    call.detail = {};
                    call.detail.element = face.element;
                    const auto bindings = bound.faces();
                    const auto mapped = std::ranges::find_if(bindings, [&](const auto& row) {
                        call.math.visit(1);
                        return row.element == face.element && row.shell == admitted.shell &&
                               row.pg == admitted.pg;
                    });
                    if (mapped == bindings.end())
                        call.fail(Code::StaleGeneration, Reason::Source);
                    const auto found = std::ranges::find_if(deck->elements, [&](const auto& item) {
                        call.math.visit(1);
                        return item.card_ordinal == mapped->source_ordinal &&
                               item.id.value() == face.element;
                    });
                    if (found == deck->elements.end())
                        call.fail(Code::StaleGeneration, Reason::Source);
                    const auto property = found->property.id;
                    call.detail.property = property;
                    check_constant_thickness(call, admitted, face, property, geometry->patch,
                                             deck->input.bulk);
                    const auto native_axis = source_axis(call, *found, deck->input, face, one);
                    const auto& property_declaration =
                        deck->properties[found->property.property_ordinal];
                    const auto* laminate =
                        property_declaration.kind == NastranShellPropertyKind::Pcomp
                            ? find_laminate(call, admitted, *found, face.winding, geometry->patch,
                                            deck->input.bulk)
                            : nullptr;
                    if ((property_declaration.kind == NastranShellPropertyKind::Pcomp &&
                         !laminate) ||
                        (laminate ? laminate->layers.size() != layers : layers != 1)) {
                        call.fail(Code::StaleGeneration, Reason::Source);
                    }
                    for (std::size_t i = 0; i < layers; ++i) {
                        call.math.visit(1);
                        const auto layer =
                            partition_take(shell.layers.at(static_cast<core::Index>(i))).record();
                        const auto id = core::LayerId{layer.id.unsigned_value()};
                        call.detail.layer = id;
                        auto degrees = call.math.copy(native_axis.degrees);
                        call.detail.ply_field.reset();
                        if (laminate) {
                            const auto& ply = laminate->layers.values()[i];
                            if (ply.layer != id)
                                call.fail(Code::StaleGeneration, Reason::Source);
                            if (const auto& angle = ply.source_angle) {
                                call.detail.ply_field = angle->field.position;
                                degrees = call.math.add(
                                    degrees, read_decimal(call.math, angle->value.decimal()));
                            }
                        }
                        if (face.winding == ExternalPgWinding::Reversed)
                            degrees = call.math.sub(Exact{}, degrees);
                        check_case_axis(call, si, layer, native_axis.tangent, native_axis.normal,
                                        degrees, one, ninety);
                        output.values()[next++] = {face.element,
                                                   property,
                                                   id,
                                                   face.winding,
                                                   call.detail.element_field,
                                                   call.detail.ply_field,
                                                   copy_ray(call.math, native_axis.tangent),
                                                   copy_ray(call.math, native_axis.normal),
                                                   std::move(degrees)};
                    }
                }
                call.math.cancel();
                return NastranPgMaterialAxes{geometry->source,  geometry->digest, sidecar_location,
                                             sidecar.digest(),  geometry->patch,  admitted.shell,
                                             std::move(output), call.math.info()};
            } catch (const geom::PartitionAbort& failure) {
                return std::unexpected(failure.failure);
            }
        });
}

} // namespace uvlas::io::detail
