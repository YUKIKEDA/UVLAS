#pragma once
#include "exact_triangle_predicates.hpp"
#include "pchip_rational_slots.hpp"

#include <array>
#include <memory>
#include <optional>
#include <span>
#include <uvlas/core/state.hpp>
#include <variant>

namespace uvlas::aero::detail {
namespace geometry_tags {
struct Point;
struct Edge;
struct Facet;
struct Face;
struct Band;
struct Query;
struct Map;
struct Reference;
struct Mapping;
struct Connectivity;
struct Policy;
struct Kernel;
struct Radius;
} // namespace geometry_tags
using PointId = core::Identifier<geometry_tags::Point>;
using EdgeId = core::Identifier<geometry_tags::Edge>;
using FacetId = core::Identifier<geometry_tags::Facet>;
using PhysicalFaceId = core::Identifier<geometry_tags::Face>;
using BandId = core::Identifier<geometry_tags::Band>;
using QueryId = core::Identifier<geometry_tags::Query>;
using ReferenceMapId = core::Identifier<geometry_tags::Map>;
using ReferenceGeometryGeneration = core::Identifier<geometry_tags::Reference>;
using MappingGeneration = core::Identifier<geometry_tags::Mapping>;
using ConnectivityGeneration = core::Identifier<geometry_tags::Connectivity>;
using KernelPolicyId = core::Identifier<geometry_tags::Policy>;
using KernelGeneration = core::Identifier<geometry_tags::Kernel>;
using RadiusGeneration = core::Identifier<geometry_tags::Radius>;

/// Every revision used by a geometry operation, including its diagnostic origin.
struct GeometryContext {
    core::ModelFingerprint fingerprint;
    core::ModelGeneration model;
    ReferenceGeometryGeneration reference;
    core::GeometryGeneration current;
    MappingGeneration mapping;
    ConnectivityGeneration connectivity;
    KernelPolicyId policy;
    KernelGeneration kernel;
    RadiusGeneration radius;
    core::PrecisionGeneration precision;
    core::NumericalGeneration numerical;
    double time;
    core::FailureDetails diagnostics;
    GeometryContext(core::ModelFingerprint fp, core::ModelGeneration m,
                    ReferenceGeometryGeneration r, core::GeometryGeneration c,
                    MappingGeneration map, ConnectivityGeneration connection, KernelPolicyId p,
                    KernelGeneration k, RadiusGeneration rad, core::PrecisionGeneration prec,
                    core::NumericalGeneration num, double seconds,
                    core::FailureDetails details = {}) noexcept
        : fingerprint(fp),
          model(m),
          reference(r),
          current(c),
          mapping(map),
          connectivity(connection),
          policy(p),
          kernel(k),
          radius(rad),
          precision(prec),
          numerical(num),
          time(seconds),
          diagnostics(details) {}
    bool operator==(const GeometryContext&) const = default;
};
using GeometryProblem = std::variant<PointId, EdgeId, FacetId, PhysicalFaceId, BandId, QueryId>;
/// Own the inner failure without recasting an aerodynamic identity to another domain.
struct OwnedGeometryFailure {
    core::Failure cause;
    std::optional<GeometryProblem> problem;
    GeometryContext context;
    OwnedGeometryFailure(core::Failure failure, std::optional<GeometryProblem> entity,
                         GeometryContext origin) noexcept
        : cause(failure), problem(entity), context(origin) {}
};
template <class T>
using AeroResult = std::expected<T, OwnedGeometryFailure>;

struct SingleWord {
    double value;
    explicit SingleWord(double x) noexcept : value(x) {}
};
struct DoubleWord {
    double high, low;
    DoubleWord(double h, double l) noexcept : high(h), low(l) {}
};
using ScalarWords = std::variant<SingleWord, DoubleWord>;
struct SinglePositions {
    std::span<const double> words;
};
struct DoublePositions {
    std::span<const double> high, low;
};
using PositionView = std::variant<SinglePositions, DoublePositions>;
struct CanonicalPoint {};
struct AliasOf {
    PointId owner;
    explicit AliasOf(PointId p) noexcept : owner(p) {}
};
using PointOwnership = std::variant<CanonicalPoint, AliasOf>;
struct PointBinding {
    ReferenceMapId map;
    MappingGeneration generation;
    PointOwnership ownership;
    PointBinding(ReferenceMapId m, MappingGeneration g, PointOwnership o) noexcept
        : map(m), generation(g), ownership(o) {}
};
struct EdgeInput {
    EdgeId id;
    std::array<PointId, 2> points;
    EdgeInput(EdgeId e, std::array<PointId, 2> p) noexcept : id(e), points(p) {}
};
struct FaceInput {
    PhysicalFaceId id;
    BandId band;
    FaceInput(PhysicalFaceId f, BandId b) noexcept : id(f), band(b) {}
};
struct TriangleInput {
    FacetId id;
    PhysicalFaceId face;
    std::array<PointId, 3> points;
    std::array<EdgeId, 3> edges;
    TriangleInput(FacetId f, PhysicalFaceId parent, std::array<PointId, 3> p,
                  std::array<EdgeId, 3> e) noexcept
        : id(f), face(parent), points(p), edges(e) {}
};
struct QuadInput {
    PhysicalFaceId face;
    std::array<PointId, 4> points;
    std::array<EdgeId, 4> edges;
    EdgeId diagonal;
    std::array<FacetId, 2> facets;
    QuadInput(PhysicalFaceId parent, std::array<PointId, 4> p, std::array<EdgeId, 4> e, EdgeId d,
              std::array<FacetId, 2> f) noexcept
        : face(parent), points(p), edges(e), diagonal(d), facets(f) {}
};
struct GeometryInput {
    std::span<const PointId> points;
    std::span<const PointBinding> bindings;
    PositionView reference, current;
    std::span<const EdgeInput> edges;
    std::span<const FaceInput> faces;
    std::span<const BandId> bands;
    std::span<const TriangleInput> triangles;
    std::span<const QuadInput> quads;
    ScalarWords core_radius;
    GeometryContext context;
    std::optional<pchip::RationalView>
        exact_time; ///< Original canonical SI instant, copied on preparation.
    GeometryInput(std::span<const PointId> p, std::span<const PointBinding> binding, PositionView r,
                  PositionView c, std::span<const EdgeInput> e, std::span<const FaceInput> f,
                  std::span<const BandId> b, std::span<const TriangleInput> t,
                  std::span<const QuadInput> q, ScalarWords radius_words, GeometryContext stamp,
                  std::optional<pchip::RationalView> instant = {}) noexcept
        : points(p),
          bindings(binding),
          reference(r),
          current(c),
          edges(e),
          faces(f),
          bands(b),
          triangles(t),
          quads(q),
          core_radius(radius_words),
          context(stamp),
          exact_time(instant) {}
};
struct GeometryCounts {
    std::uint64_t points, edges, faces, bands, triangles, quads;
    GeometryCounts(std::uint64_t p, std::uint64_t e, std::uint64_t f, std::uint64_t b,
                   std::uint64_t t, std::uint64_t q) noexcept
        : points(p), edges(e), faces(f), bands(b), triangles(t), quads(q) {}
};
struct GeometryLimits {
    GeometryCounts counts;
    std::uint64_t max_steps;
    PredicateLimits predicate;
    GeometryLimits(GeometryCounts c, std::uint64_t steps, PredicateLimits p) noexcept
        : counts(c), max_steps(steps), predicate(p) {}
};
struct FreeSingle {
    std::array<double, 3> words;
};
struct FreeDouble {
    std::array<double, 3> high, low;
};
using PointWords = std::variant<FreeSingle, FreeDouble>;
struct FreePoint {
    PointWords words;
    explicit FreePoint(PointWords p) noexcept : words(p) {}
};
struct SourceVertex {
    PointId point;
    explicit SourceVertex(PointId p) noexcept : point(p) {}
};
struct SourceEdgePoint {
    EdgeId edge;
    ScalarWords t;
    SourceEdgePoint(EdgeId e, ScalarWords parameter) noexcept : edge(e), t(parameter) {}
};
struct SourceFacePoint {
    FacetId facet;
    ScalarWords u, v;
    SourceFacePoint(FacetId f, ScalarWords a, ScalarWords b) noexcept : facet(f), u(a), v(b) {}
};
using QueryDefinition = std::variant<FreePoint, SourceVertex, SourceEdgePoint, SourceFacePoint>;
struct QueryRequest {
    QueryId id;
    QueryDefinition definition;
    QueryRequest(QueryId q, QueryDefinition d) noexcept : id(q), definition(d) {}
};
struct NoMaterialSupport {
    bool operator==(const NoMaterialSupport&) const = default;
};
using CanonicalSupport = std::variant<NoMaterialSupport, PointId, EdgeId, FacetId>;
struct GlobalGeometryNotChecked {};
struct MappingDerivativesNotChecked {};
enum class Configuration : std::uint8_t { Reference, Current };
struct IndependentP1Scatter {
    std::array<std::int8_t, 9> coefficients;
};
struct QuadLRScatter {
    std::array<std::int8_t, 6> coefficients;
};
using FacetScatter = std::variant<IndependentP1Scatter, QuadLRScatter>;
struct GeometryData;
struct QueryStorage;
struct QueryArithmeticScratch;
class GeometryEngine;
class GeometryWorkspace;
class SharedP1Engine;
/// Retain immutable point/topology payload; getters never expose mutable arrays.
class OwnedGeometry {
public:
    OwnedGeometry(OwnedGeometry&&) noexcept = default;
    OwnedGeometry& operator=(OwnedGeometry&&) noexcept = default;
    OwnedGeometry(const OwnedGeometry&) = delete;
    OwnedGeometry& operator=(const OwnedGeometry&) = delete;
    ~OwnedGeometry() = default;
    [[nodiscard]] bool active() const noexcept {
        return data_ != nullptr;
    }
    [[nodiscard]] GeometryContext context() const noexcept {
        return context_;
    }
    [[nodiscard]] AeroResult<std::array<std::byte, 32>> identity() const noexcept;
    [[nodiscard]] AeroResult<std::optional<pchip::RationalView>> exact_time() const& noexcept;
    [[nodiscard]] AeroResult<std::optional<pchip::RationalView>> exact_time() const&& = delete;
    [[nodiscard]] AeroResult<GeometryCounts> counts() const noexcept;
    [[nodiscard]] AeroResult<std::uint64_t> payload_bytes() const noexcept;
    [[nodiscard]] AeroResult<PointId> owner(PointId point) const noexcept;
    [[nodiscard]] AeroResult<std::array<ScalarWords, 3>> point_words(PointId,
                                                                     Configuration) const noexcept;
    [[nodiscard]] AeroResult<ScalarWords> core_radius() const noexcept;
    [[nodiscard]] AeroResult<PhysicalFaceId> facet_face(FacetId) const noexcept;
    [[nodiscard]] AeroResult<BandId> face_band(PhysicalFaceId) const noexcept;
    [[nodiscard]] AeroResult<std::array<PointId, 3>> facet_points(FacetId facet) const noexcept;
    [[nodiscard]] AeroResult<std::array<EdgeId, 3>> facet_edges(FacetId facet) const noexcept;
    [[nodiscard]] AeroResult<std::array<std::int8_t, 3>> incidence(FacetId facet) const noexcept;
    /// Columns map to independent P1 basis (3) or quad L/R values (2).
    [[nodiscard]] AeroResult<FacetScatter> scatter(FacetId facet) const noexcept;
    [[nodiscard]] static constexpr GlobalGeometryNotChecked global_status() noexcept {
        return {};
    }
    [[nodiscard]] static constexpr MappingDerivativesNotChecked derivative_status() noexcept {
        return {};
    }

private:
    friend class GeometryEngine;
    friend class P1Engine;
    friend class GlobalEngine;
    friend class SharedP1Engine;
    explicit OwnedGeometry(std::shared_ptr<const GeometryData> data, GeometryContext stamp) noexcept
        : data_(std::move(data)), context_(stamp) {}
    std::shared_ptr<const GeometryData> data_;
    GeometryContext context_;
};
/// Own original query words/expression and its geometry, never a rounded derived point.
class ResolvedQuery {
public:
    ResolvedQuery(ResolvedQuery&&) noexcept;
    ResolvedQuery& operator=(ResolvedQuery&&) noexcept;
    ResolvedQuery(const ResolvedQuery&) = delete;
    ResolvedQuery& operator=(const ResolvedQuery&) = delete;
    ~ResolvedQuery();
    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] AeroResult<QueryRequest> definition() const noexcept;
    [[nodiscard]] AeroResult<CanonicalSupport> support() const noexcept;
    [[nodiscard]] GeometryContext context() const noexcept {
        return context_;
    }
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;

private:
    friend class GeometryEngine;
    friend class SharedP1Engine;
    ResolvedQuery(std::shared_ptr<const GeometryData>, core::PlannedBuffer<QueryStorage>,
                  GeometryContext) noexcept;
    std::shared_ptr<const GeometryData> geometry_;
    core::PlannedBuffer<QueryStorage> query_;
    GeometryContext context_;
};
/// The underlying predicate lease protects both predicate and interpolation scratch.
class GeometryWorkspace {
    class Key {
        friend class GeometryWorkspace;
        Key() = default;
    };

public:
    explicit GeometryWorkspace(Key, std::unique_ptr<PredicateWorkspace>,
                               core::PlannedBuffer<QueryArithmeticScratch>) noexcept;
    [[nodiscard]] static core::Result<std::unique_ptr<GeometryWorkspace>>
    create(core::ResourcePlan&, core::FailureDetails = {}) noexcept;
    GeometryWorkspace(const GeometryWorkspace&) = delete;
    GeometryWorkspace& operator=(const GeometryWorkspace&) = delete;
    GeometryWorkspace(GeometryWorkspace&&) = delete;
    GeometryWorkspace& operator=(GeometryWorkspace&&) = delete;
    ~GeometryWorkspace();
    [[nodiscard]] core::Result<PredicateLease> acquire(core::FailureDetails = {}) noexcept;
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;

private:
    friend class GeometryEngine;
    std::unique_ptr<PredicateWorkspace> predicate_;
    core::PlannedBuffer<QueryArithmeticScratch> arithmetic_;
};
/// Snapshot all finite words and typed connectivity; no global/derivative qualification.
/// Sources use explicit stable facet/diagonal IDs. Output owns all variable payload.
[[nodiscard]] AeroResult<OwnedGeometry>
prepare_owned_geometry(const GeometryInput&, const GeometryLimits&, GeometryWorkspace&,
                       core::ResourcePlan&, std::stop_token = {}) noexcept;
/// Resolve an original expression against this exact immutable content and context.
/// Parameters are exact word sums. Normalized support never changes the query identity.
[[nodiscard]] AeroResult<ResolvedQuery> resolve_query(const OwnedGeometry&, const QueryRequest&,
                                                      const GeometryContext&, const GeometryLimits&,
                                                      GeometryWorkspace&, core::ResourcePlan&,
                                                      std::stop_token = {}) noexcept;
/// Evaluate the original expression at fixed q and classify it mathematically.
/// The result has no epsilon/self/contact qualification or global geometry token.
[[nodiscard]] AeroResult<ExactTriangleResult>
classify_resolved_query(const OwnedGeometry&, const ResolvedQuery&, FacetId, const GeometryLimits&,
                        GeometryWorkspace&, core::ResourcePlan&, std::stop_token = {}) noexcept;
} // namespace uvlas::aero::detail
