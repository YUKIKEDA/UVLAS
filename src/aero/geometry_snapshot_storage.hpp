#pragma once

#include "geometry_snapshot.hpp"

#include <type_traits>

namespace uvlas::aero::detail {
// These are private, budgeted storage records. Zero initializes a builder slot;
// it is never interpreted as an absent physical ID or published before filling.
struct AliasRequestStorage {
    std::uint64_t owner = 0;
};
struct PointStorage {
    std::uint64_t id = 0;
    core::Index owner = 0;
    std::uint64_t map = 0, map_generation = 0;
    std::variant<CanonicalPoint, AliasRequestStorage> original;
};
struct EdgeStorage {
    std::uint64_t id = 0;
    std::array<std::uint64_t, 2> original{};
    std::array<core::Index, 2> points{};
};
struct FaceStorage {
    std::uint64_t id = 0, band = 0;
};
struct TriangleStorage {
    std::uint64_t id = 0, parent = 0;
    std::array<std::uint64_t, 3> points{}, edges{};
};
struct QuadStorage {
    std::uint64_t parent = 0;
    std::array<std::uint64_t, 4> points{}, edges{};
    std::uint64_t diagonal = 0;
    std::array<std::uint64_t, 2> facets{};
};
struct IndependentOrigin {
    core::Index source = 0;
};
struct QuadOrigin {
    core::Index source = 0;
    std::uint8_t half = 0;
};
struct FacetStorage {
    std::uint64_t id = 0;
    core::Index face = 0;
    std::array<core::Index, 3> points{}, edges{};
    std::array<std::int8_t, 3> incidence{};
    std::variant<IndependentOrigin, QuadOrigin> origin;
};
struct IdLookup {
    std::uint64_t id = 0;
    core::Index ordinal = 0;
};
using ParameterStorage = std::variant<double, std::array<double, 2>>;
using FreeSingleStorage = std::array<double, 3>;
using FreeDoubleStorage = std::array<double, 6>;
struct VertexStorage {
    std::uint64_t point = 0;
};
struct EdgeQueryStorage {
    std::uint64_t edge = 0;
    ParameterStorage t;
};
struct FaceQueryStorage {
    std::uint64_t facet = 0;
    ParameterStorage u, v;
};
struct PointSupportStorage {
    std::uint64_t id = 0;
};
struct EdgeSupportStorage {
    std::uint64_t id = 0;
};
struct FacetSupportStorage {
    std::uint64_t id = 0;
};
struct QueryStorage {
    std::uint64_t id = 0;
    std::variant<FreeSingleStorage, FreeDoubleStorage, VertexStorage, EdgeQueryStorage,
                 FaceQueryStorage>
        definition;
    std::variant<NoMaterialSupport, PointSupportStorage, EdgeSupportStorage, FacetSupportStorage>
        support;
};
// These selected alternatives construct/copy/move without a throwing inner
// variant conversion, so no outer strong-guarantee temporary is needed.
static_assert(std::is_nothrow_copy_constructible_v<FreeSingleStorage>);
static_assert(std::is_nothrow_copy_constructible_v<FreeDoubleStorage>);
static_assert(std::is_nothrow_move_constructible_v<EdgeQueryStorage>);
static_assert(std::is_nothrow_copy_constructible_v<EdgeQueryStorage>);
static_assert(std::is_nothrow_move_constructible_v<FaceQueryStorage>);
static_assert(std::is_nothrow_copy_constructible_v<FaceQueryStorage>);
static_assert(std::is_nothrow_copy_assignable_v<QueryStorage>);
struct QueryArithmeticScratch {
    std::array<ExactInteger, 3> weights;
    ExactInteger word, coordinate, product, sum;
};
struct GeometryData {
    GeometryContext context;
    GeometryCounts counts;
    ScalarWords radius;
    std::uint8_t reference_words, current_words;
    std::array<std::byte, 32> identity{};
    std::uint64_t payload = 0;
    core::PlannedBuffer<double> reference_high, reference_low, current_high, current_low;
    core::PlannedBuffer<PointStorage> points;
    core::PlannedBuffer<EdgeStorage> edges;
    core::PlannedBuffer<FaceStorage> faces;
    core::PlannedBuffer<std::uint64_t> bands;
    core::PlannedBuffer<TriangleStorage> triangles;
    core::PlannedBuffer<QuadStorage> quads;
    core::PlannedBuffer<FacetStorage> facets;
    core::PlannedBuffer<IdLookup> point_lookup, edge_lookup, face_lookup, band_lookup, facet_lookup;
    GeometryData(GeometryContext c, GeometryCounts n, ScalarWords r, std::uint8_t rw,
                 std::uint8_t cw) noexcept
        : context(c), counts(n), radius(r), reference_words(rw), current_words(cw) {}
};
static_assert(sizeof(EdgeStorage) == 40 && sizeof(FaceStorage) == 16);
static_assert(sizeof(TriangleStorage) == 64 && sizeof(QuadStorage) == 96 && sizeof(IdLookup) == 16);
static_assert(sizeof(QueryArithmeticScratch) == 7 * sizeof(ExactInteger));
} // namespace uvlas::aero::detail
