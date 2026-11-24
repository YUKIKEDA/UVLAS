#include "nastran_pcomp_physical.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Exact = structure::ExactRational;
using Reason = NastranPcompGeometryFailure;
using Code = core::FailureCode;
constexpr auto operation = core::Operation::PrepareModel;

struct Stop {
    core::Failure error;
};

class Comparison {
public:
    Comparison(const NastranPcompGeometryLimits& configured, core::ResourcePlan& resource,
               const std::stop_token& token,
               std::optional<NastranPcompGeometryDiagnostic>* output) noexcept
        : limits(configured), plan(resource), stop(token), diagnostic_(output) {}

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
    [[nodiscard]] T value(core::Result<T> result, Reason reason = Reason::Source) const {
        if (!result) {
            report(result.error(), reason);
        }
        return *std::move(result);
    }

    template <class T>
    [[nodiscard]] T require(std::optional<T> result) const {
        if (!result) {
            reject(Reason::Source);
        }
        return *std::move(result);
    }

    template <class T>
    [[nodiscard]] core::PlannedBuffer<T> memory(std::size_t count) const {
        if (!std::in_range<core::Index>(count)) {
            reject(Reason::Resource, Code::ResourceExceeded);
        }
        return value(plan.allocate<T>(static_cast<core::Index>(count), operation),
                     Reason::Resource);
    }

    [[nodiscard]] Exact number(std::string_view text) const {
        return value(Exact::from_decimal(text, limits.exact, plan));
    }

    const NastranPcompGeometryLimits& limits;
    core::ResourcePlan& plan;
    const std::stop_token& stop;
    NastranPcompGeometryDiagnostic detail{};
    std::uint64_t work{};

private:
    [[noreturn]] void report(core::Failure error, Reason reason) const {
        if (diagnostic_ && error.code != Code::Cancelled) {
            auto position = detail;
            position.reason = reason;
            *diagnostic_ = position;
        }
        throw Stop{error};
    }

    std::optional<NastranPcompGeometryDiagnostic>* diagnostic_;
};

[[nodiscard]] const NastranPcomp& find_property(Comparison& check, const NastranLaminateView& deck,
                                                std::size_t ordinal) {
    const auto found = std::ranges::find_if(deck.pcomp, [&](const auto& property) {
        check.step();
        return property.property_ordinal == ordinal;
    });
    if (found == deck.pcomp.end()) {
        check.reject(Reason::Source);
    }
    return *found;
}

[[nodiscard]] bool symmetry(Comparison& check, const NastranPcomp& property) {
    switch (property.mode) {
    case NastranLaminateMode::FullStack:
        return false;
    case NastranLaminateMode::SymmetricHalf:
        return true;
    default:
        if (const auto& mode = property.mode_field) {
            check.detail.field = mode->position;
        }
        check.reject(Reason::Mode, Code::Unsupported);
    }
}

[[nodiscard]] bool centered_ply(Comparison& check, bool symmetric, std::size_t declared,
                                std::size_t actual) {
    if (!symmetric) {
        if (declared != actual) {
            check.reject(Reason::LayerCount);
        }
        return false;
    }
    if (declared > std::numeric_limits<std::size_t>::max() / 2) {
        check.reject(Reason::Resource, Code::ResourceExceeded);
    }
    const auto mirrored = 2 * declared;
    if (actual != mirrored && actual != mirrored - 1) {
        check.reject(Reason::LayerCount);
    }
    // The deck specifies half of a center ply when present. Its unsplit or
    // separated Case counterpart is explicit in the checked layer intervals.
    return actual != mirrored;
}

[[nodiscard]] Exact length_unit(Comparison& check, SiLengthUnit length) {
    const SiRecipeView unit{{SiDecimalView{1, "1", 0}},
                            {},
                            {length, SiMassUnit::Kilogram, SiTimeUnit::Second, SiAngleUnit::Radian},
                            SchemaDimension{{1, 0, 0, 0}},
                            0};
    return check.value(materialize_rational(unit, check.limits.exact, check.plan));
}

[[nodiscard]] const NastranShellReal&
thickness(Comparison& check, std::span<const NastranPcompPly> stack, const NastranPcompPly& ply) {
    if (const auto& declared = stack[ply.thickness_source_ply].thickness) {
        check.detail.field = declared->field.position;
        return *declared;
    }
    check.reject(Reason::Source);
}

[[nodiscard]] NastranMaterialReference
material(Comparison& check, std::span<const NastranPcompPly> stack, const NastranPcompPly& ply) {
    return check.require(stack[ply.material_source_ply].material);
}

[[nodiscard]] Exact full_thickness(Comparison& check, const NastranPcomp& property,
                                   std::span<const NastranPcompPly> stack, const Exact& unit,
                                   std::span<Exact> plies, bool symmetric) {
    Exact total;
    for (std::size_t i = 0; i < stack.size(); ++i) {
        check.step();
        const auto& source = thickness(check, stack, stack[i]);
        const auto input = check.number(source.value.decimal());
        if (input.sign() <= 0) {
            check.reject(Reason::FullThickness);
        }
        plies[i] = check.value(input.multiply(unit, check.limits.exact, check.plan));
        total = check.value(total.add(plies[i], check.limits.exact, check.plan));
    }
    if (symmetric) {
        total = check.value(total.add(total, check.limits.exact, check.plan));
    }
    if (const auto& offset = property.bottom_offset) {
        check.detail.field = offset->field.position;
        const auto source = check.number(offset->value.decimal());
        const auto metres = check.value(source.multiply(unit, check.limits.exact, check.plan));
        const auto twice = check.value(metres.add(metres, check.limits.exact, check.plan));
        if (check.value(twice.add(total, check.limits.exact, check.plan)).sign() != 0) {
            check.reject(Reason::ReferencePlane);
        }
    }
    return total;
}

void check_active_thickness(Comparison& check, const structure::PgIntegrationPartition& partition,
                            const Exact& expected) {
    const auto source = check.value(partition.thickness());
    const auto pieces = check.value(partition.pieces());
    auto visited = check.memory<bool>(source.values.size());
    std::uint64_t comparisons = 0;
    // An affine triangle constant on a positive-area subset has zero gradient
    // throughout that triangle. This admits cropped active regions without
    // inspecting unrelated, wholly inactive input triangles.
    for (const auto& piece : pieces) {
        check.step();
        const auto& triangle = source.triangles[static_cast<std::size_t>(piece.thickness_triangle)];
        for (const auto vertex : triangle) {
            check.step();
            const auto i = static_cast<std::size_t>(vertex);
            if (visited.values()[i]) {
                continue;
            }
            if (comparisons == check.limits.thickness_vertices) {
                check.reject(Reason::Resource, Code::ResourceExceeded);
            }
            ++comparisons;
            if (source.values[i] != expected) {
                check.detail.thickness_vertex = vertex;
                check.reject(Reason::FullThickness);
            }
            visited.values()[i] = true;
        }
    }
    if (comparisons == 0) {
        check.reject(Reason::Source);
    }
}

[[nodiscard]] core::PlannedBuffer<NastranPcompLayerMatch>
compare_layers(Comparison& check, std::span<const NastranPcompPly> stack,
               std::span<const Exact> plies,
               std::span<const structure::PgPartitionOwnedLayer> layers, const Exact& total,
               bool joined_center, ExternalPgWinding winding) {
    auto matches = check.memory<NastranPcompLayerMatch>(layers.size());
    const auto one = check.number("1");
    const Exact zero;
    auto boundary = check.number("-1");
    Exact elapsed;
    for (std::size_t i = 0; i < layers.size(); ++i) {
        check.step();
        const auto original = winding == ExternalPgWinding::Forward ? i : layers.size() - 1 - i;
        const auto index = original < stack.size() ? original : layers.size() - 1 - original;
        const auto& ply = stack[index];
        const auto& source = thickness(check, stack, ply);
        const auto& layer = layers[i];
        check.detail.case_layer = layer.id;
        if (boundary != layer.lower) {
            check.reject(Reason::Interface);
        }
        const auto& depth = plies[index];
        const bool center = joined_center && original == stack.size() - 1;
        auto realised =
            check.value(depth.add(center ? depth : zero, check.limits.exact, check.plan));
        elapsed = check.value(elapsed.add(realised, check.limits.exact, check.plan));
        const auto twice = check.value(elapsed.add(elapsed, check.limits.exact, check.plan));
        const auto unit_fraction = check.value(twice.divide(total, check.limits.exact, check.plan));
        boundary = check.value(unit_fraction.subtract(one, check.limits.exact, check.plan));
        if (boundary != layer.upper) {
            check.reject(Reason::Interface);
        }
        matches.values()[i] = {
            layer.id,           i, index, material(check, stack, ply), source.field.position,
            std::move(realised)};
    }
    return matches;
}

} // namespace

core::Result<NastranPcompGeometryMatch> compare_nastran_pcomp_geometry(
    const NastranLaminateCards& input, std::size_t property_ordinal, SiLengthUnit length,
    ExternalPgWinding winding, const structure::PgReferenceGeometry& pg,
    const NastranPcompGeometryLimits& limits, core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<NastranPcompGeometryDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranPcompGeometryMatch> {
            try {
                Comparison check{limits, plan, stop, diagnostic};
                check.step();
                const auto deck = check.require(input.view());
                const auto& property = find_property(check, deck, property_ordinal);
                const auto& base = deck.input.input.properties[property.property_ordinal];
                check.detail.field = base.id_field.position;
                const bool symmetric = symmetry(check, property);
                const auto stack = deck.layers(property);
                if (stack.empty()) {
                    check.reject(Reason::Source);
                }
                const auto identity = check.value(pg.identity());
                const auto source = check.value(pg.source());
                const auto layers = check.value(source.partition.layers());
                const auto center = centered_ply(check, symmetric, stack.size(), layers.size());
                auto plies = check.memory<Exact>(stack.size());
                const auto unit = length_unit(check, length);
                auto total =
                    full_thickness(check, property, stack, unit, plies.values(), symmetric);
                check.detail.field = base.id_field.position;
                check_active_thickness(check, source.partition, total);
                auto matches =
                    compare_layers(check, stack, plies.values(), layers, total, center, winding);
                check.step();
                const auto& bulk = deck.input.input.input.bulk;
                return NastranPcompGeometryMatch{
                    bulk.location,    bulk.digest,        identity,  winding,
                    std::move(total), std::move(matches), check.work};
            } catch (const Stop& error) {
                return std::unexpected(error.error);
            }
        });
}

} // namespace uvlas::io::detail
