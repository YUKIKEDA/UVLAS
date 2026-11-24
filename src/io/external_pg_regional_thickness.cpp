#include "external_pg_regional_thickness.hpp"

#include "../structure/pg_reference_domain.hpp"
#include "external_pg_chart_clip.hpp"

#include <algorithm>

namespace uvlas::io::detail {
namespace {
namespace ref = structure::detail;
namespace chart = external_pg_chart;
using ref::partition_fail;
using ref::partition_take;
using Code = core::FailureCode;
using Reason = ExternalPgThicknessFailure;
constexpr auto operation = core::Operation::PrepareModel;

struct Diagnosis {
    std::optional<ExternalPgThicknessDiagnostic>* output;
    ExternalPgThicknessDiagnostic value{};

    void select(const ExternalPgThicknessClaim& claim) noexcept {
        value = {Reason::Selection, claim.element, claim.sidecar_node, claim.field, {}, {}};
    }

    [[noreturn]] void fail(Reason reason, Code code) const {
        if (output) {
            auto location = value;
            location.reason = reason;
            *output = location;
        }
        partition_fail(code);
    }
};

[[nodiscard]] const ExternalPgOwnedFace&
resolve_face(ref::ReferenceArithmetic& math, const ExternalPlanarPgView& source, std::uint64_t id,
             std::span<ExternalPgThicknessMatch> previous, const Diagnosis& diagnosis) {
    const auto selected = std::ranges::find_if(source.faces, [&](const auto& face) {
        math.visit(1);
        return face.element == id;
    });
    const bool duplicate = std::ranges::any_of(previous, [&](const auto& match) {
        math.visit(1);
        return match.element == id;
    });
    if (selected == source.faces.end() || duplicate) {
        diagnosis.fail(Reason::Selection, Code::InvalidInput);
    }
    return *selected;
}

[[nodiscard]] std::array<const chart::Point*, 4>
positive_ring(const ExternalPgOwnedFace& face) noexcept {
    const auto count = face.ring().size();
    std::array<const chart::Point*, 4> ring{};
    for (std::size_t i = 0; i < count; ++i) {
        const auto native = face.winding == ExternalPgWinding::Forward ? i : (count - i) % count;
        ring[i] = &face.corners[native].chart;
    }
    return ring;
}

struct CasePlanes {
    explicit CasePlanes(const structure::PgIntegrationPartition& partition)
        : triangles(partition_take(partition.triangles())),
          vertices(partition_take(partition.vertices())),
          pieces(partition_take(partition.pieces())),
          thickness(partition_take(partition.thickness())) {}

    std::span<const structure::PgPartitionTriangle> triangles;
    std::span<const chart::Point> vertices;
    std::span<const structure::PgPartitionPiece> pieces;
    structure::PgPartitionThicknessView thickness;
};

[[nodiscard]] ExternalPgThicknessMatch
compare_claim(ref::ReferenceArithmetic& math, const CasePlanes& pg, const ExternalPgOwnedFace& face,
              const ExternalPgThicknessClaim& claim, Diagnosis& diagnosis) {
    const auto* value = claim.value_metres;
    if (!value || value->sign() <= 0) {
        diagnosis.fail(Reason::Selection, Code::InvalidInput);
    }
    const auto outward = positive_ring(face);
    const auto domain = std::span{outward}.first(face.ring().size());
    ExternalPgThicknessMatch match{face.element,       face.winding, {}, {}, 0, 0,
                                   claim.sidecar_node, claim.field};
    for (const auto& triangle : pg.triangles) {
        math.visit(1);
        const chart::Triangle basis{&ref::reference_at(pg.vertices, triangle.vertices[0]),
                                    &ref::reference_at(pg.vertices, triangle.vertices[1]),
                                    &ref::reference_at(pg.vertices, triangle.vertices[2])};
        const auto overlap = chart::clip(math, basis, domain);
        const auto area = chart::area2(math, overlap);
        if (area.sign() <= 0) {
            continue;
        }
        const auto original = ref::reference_at(pg.pieces, triangle.piece).thickness_triangle;
        const auto plane = ref::thickness_plane(math, pg.thickness, original);
        // The difference is affine on the entire clipped piece. Matching all
        // its exact vertices proves the full piece, including the chosen trace.
        const auto vertices = overlap.vertices();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            math.visit(1);
            const auto h = plane.evaluate(math, vertices[i].s, vertices[i].t);
            if (math.compare(h, *value) != 0) {
                diagnosis.value.case_triangle = original;
                diagnosis.value.overlap_vertex = i;
                diagnosis.fail(Reason::Difference, Code::InvalidInput);
            }
        }
        match.twice_chart_area = math.add(match.twice_chart_area, area);
        match.pieces = ref::partition_add(match.pieces, 1);
        match.vertices = ref::partition_add(match.vertices, vertices.size());
    }
    if (match.pieces == 0) {
        diagnosis.fail(Reason::Source, Code::GeometryInvalid);
    }
    match.value_metres = math.copy(*value);
    return match;
}

} // namespace

core::Result<ExternalPgThicknessReport> compare_external_pg_regional_thickness(
    const ExternalPlanarPg& certified, const structure::PgReferenceGeometry& reference,
    std::span<const ExternalPgThicknessClaim> claims, const ExternalPgThicknessLimits& limits,
    core::ResourcePlan& plan, std::optional<ExternalPgThicknessDiagnostic>* output) noexcept {
    if (output) {
        output->reset();
    }
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<ExternalPgThicknessReport> {
            try {
                ref::ReferenceArithmetic math{limits.arithmetic, plan};
                math.cancel();
                Diagnosis diagnosis{output};
                const auto source = certified.view();
                if (!source) {
                    diagnosis.fail(Reason::Source, Code::InvalidInput);
                }
                const auto identity = partition_take(reference.identity());
                if (identity != source->patch) {
                    diagnosis.fail(Reason::Source, Code::StaleGeneration);
                }
                if (claims.size() > limits.maximum_claims) {
                    partition_fail(Code::ResourceExceeded);
                }
                auto match = partition_take(plan.allocate<ExternalPgThicknessMatch>(
                    ref::partition_index(claims.size()), operation));
                if (!claims.empty()) {
                    math.admit();
                    const auto input = partition_take(reference.source());
                    const CasePlanes pg{input.partition};
                    for (std::size_t i = 0; i < claims.size(); ++i) {
                        math.visit(1);
                        const auto& claim = claims[i];
                        diagnosis.select(claim);
                        const auto& face = resolve_face(math, *source, claim.element,
                                                        match.values().first(i), diagnosis);
                        match.values()[i] = compare_claim(math, pg, face, claim, diagnosis);
                    }
                }
                math.cancel();
                return ExternalPgThicknessReport{source->source, source->digest, identity,
                                                 std::move(match), math.info()};
            } catch (const ref::PartitionAbort& error) {
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
