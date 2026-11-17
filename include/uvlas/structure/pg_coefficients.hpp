#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <stop_token>
#include <utility>
#include <uvlas/core/instrumentation.hpp>
#include <uvlas/core/resource_plan.hpp>
#include <uvlas/core/state.hpp>
#include <uvlas/structure/pg_graph.hpp>

namespace uvlas::structure {

/// The three independent numerators of R + z P/g^(3/2) + z^2 C/2.
enum class PgCoefficientField : std::uint8_t {
    R, ///< Midsurface numerator.
    P, ///< Director numerator; never the director itself.
    C  ///< Independent quadratic-through-thickness numerator.
};
/// Origin of a scalar column; raw and appended anchor columns remain distinct.
enum class PgCoefficientKind : std::uint8_t {
    BasicR,     ///< Tensor C1 quadratic B-spline.
    ReferenceR, ///< Graph-height completion, including an exact zero column.
    BasicD,     ///< Tensor C0 quadratic B-spline multiplied by g.
    NormalD,    ///< One of the three normal-numerator completions.
    Paired,     ///< One component of an s/t integrated lifting.
    BasicC,     ///< Tensor continuous piecewise-linear thickness field.
    AnchorR,    ///< Separate R-only extension of one raw anchor trace.
    AnchorP     ///< Separate P-only extension of one raw anchor trace.
};
/// Axis of a paired seed; it is not a Cartesian displacement component.
enum class PgPairedAxis : std::uint8_t {
    S, ///< Curvature seed along chart s.
    T  ///< Curvature seed along chart t.
};
/// Unambiguous scalar column index, before Cartesian expansion.
struct PgCoefficientColumn {
    core::Index value = 0; ///< Zero-based scalar ordinal in the returned dictionary.
    bool operator==(const PgCoefficientColumn&) const = default;
};
/// Provenance of one scalar column, including zero or duplicate source columns.
struct PgCoefficientSource {
    PgCoefficientKind kind = PgCoefficientKind::BasicR; ///< Construction family.
    core::Index first = 0;               ///< First tensor basis index, or paired origin.
    core::Index second = 0;              ///< Second tensor index, or Legendre degree.
    std::uint8_t component = 0;          ///< Lifting/normal component j; unrelated to dof alpha.
    PgPairedAxis axis = PgPairedAxis::S; ///< Relevant only for paired columns.
    PgCoefficientColumn original{};      ///< Raw source ordinal; anchor columns retain it.
};
/// One stored nonzero cell numerator; zero fields have no allocated term storage.
struct PgCoefficientRecord {
    PgCoefficientColumn column{}; ///< Scalar column, independent of Cartesian alpha.
    PgCellIndex cell{};           ///< Background half-cell.
    PgCoefficientField field = PgCoefficientField::R; ///< Numerator identity.
    ExactCellPolynomial polynomial;                   ///< Owned local-u/v canonical coefficients.
};
/// Sparse exact coordinates of the 30 mandatory fields and three translations.
/// Required fields are X,Y,f; j-major eight s-global curvatures; normals X/Y/Z;
/// then translations X/Y/Z. Translation columns share the constant scalar map
/// and differ in their Cartesian direction. No rank or independence is implied.
struct PgRequiredMapEntry {
    PgCoefficientColumn source{}; ///< Raw scalar dictionary column.
    std::uint8_t target = 0;      ///< Required-field index in [0,33).
    ExactRational coefficient;    ///< Exact coordinate, never a rounded fit.
};
/// Coordinates of ell(s) times one C anchor trace in the existing C family.
struct PgAnchorCMapEntry {
    PgCoefficientColumn source{}; ///< Original C column whose trace is extended.
    PgCoefficientColumn target{}; ///< Existing C column used for reconstruction.
    ExactRational coefficient;    ///< Exact reconstruction coordinate.
};
/// Whole-call caps. Zero never means unlimited, and no column or small coefficient
/// is dropped to satisfy a cap. External memory covers retained call values and
/// active Boost scratch; ResourcePlan covers cumulative owned payload.
struct PgCoefficientLimits {
    /// Supply every algebra, capacity and cancellation limit explicitly.
    PgCoefficientLimits(ExactLimits exact_limits, std::uint64_t columns, std::uint64_t terms,
                        std::uint64_t map_entries, std::stop_token cancellation_token) noexcept
        : algebra(exact_limits),
          maximum_columns(columns),
          maximum_stored_terms(terms),
          maximum_map_entries(map_entries),
          cancellation(std::move(cancellation_token)) {}
    ExactLimits algebra;           ///< Explicit value/work/raw/degree and external-memory policy.
    std::uint64_t maximum_columns; ///< Output-capacity admission limit; no clipping.
    std::uint64_t maximum_stored_terms; ///< Total retained nonzero-term admission envelope.
    std::uint64_t maximum_map_entries;  ///< Required/translation and C-map capacity limit.
    std::stop_token cancellation;       ///< Checked between bounded stages and before publication.
};
/// Product-owned preparation phases, absent when timing is disabled.
struct PgCoefficientTimings {
    core::StageTiming coefficients; ///< Basis, moments and all raw coefficient construction.
    core::StageTiming anchor;       ///< Trace selection and separate R/P extensions.
    core::StageTiming verification; ///< Required-field reconstruction and exact identities.
};
/// Owning diagnostics. Forecasts are not measured heap or process RSS bounds.
struct PgCoefficientPreparation {
    ExactPreparationInfo arithmetic; ///< Maximum work/raw and whole-call external estimate.
    std::uint64_t planned_owned_payload_upper_bytes = 0; ///< Cumulative allocation envelope.
    std::uint64_t allocated_owned_payload_bytes = 0;     ///< Actual successful ResourcePlan delta.
    std::uint64_t planned_column_capacity = 0; ///< Conservative full-V capacity, not actual rank.
    std::uint64_t observed_columns = 0; ///< Actual raw-plus-anchor columns, including raw zeros.
    std::uint64_t observed_nonzero_terms = 0; ///< Canonical terms in stored numerator records/g.
    std::uint64_t observed_maximum_value_bits = 0; ///< Observed numerator/denominator bit maximum.
    std::optional<PgCoefficientTimings> timings;   ///< Disabled means no clock reads or sample.
};
namespace detail {
struct PgCoefficientData;
struct PgCoefficientFactory;
} // namespace detail

/// Immutable owned raw dictionary and full V. It carries exact algebraic
/// reconstruction, not activity rank, conditioned N, quadrature or physical
/// geometry qualification. Each synchronous preparation holds an exclusive
/// caller workspace lease; outputs never borrow that workspace or input graph.
class PgCoefficientDictionary {
public:
    /// Transfer ownership; borrowed views expire and moved-from queries fail.
    PgCoefficientDictionary(PgCoefficientDictionary&&) noexcept;
    /// Replace owned storage, invalidating prior views.
    PgCoefficientDictionary& operator=(PgCoefficientDictionary&&) noexcept;
    /// Copying cannot bypass explicit preparation budgets.
    PgCoefficientDictionary(const PgCoefficientDictionary&) = delete;
    /// Assignment cannot bypass explicit preparation budgets.
    PgCoefficientDictionary& operator=(const PgCoefficientDictionary&) = delete;
    /// Release owned coefficients and graph snapshot.
    ~PgCoefficientDictionary();
    /// Construct all adopted families, full anchor completion and exact maps.
    /// Invalid tags/moved owners/reentrant workspace use return InvalidInput;
    /// cap, overflow and allocation failures are ResourceExceeded; cancellation
    /// is Cancelled. Unproven hard external-memory bounds return Unsupported.
    /// No input or previous result is modified, even when a stage fails.
    [[nodiscard]] static core::Result<PgCoefficientDictionary>
    prepare(const PgGraph&, const PgCoefficientLimits&, core::Workspace&,
            core::ResourcePlan&) noexcept;
    /// Borrow the owned graph snapshot, or InvalidInput after move.
    [[nodiscard]] core::Result<PgGraphView> graph() const& noexcept;
    [[nodiscard]] core::Result<PgGraphView> graph() const&& = delete;
    /// Verify actual graph contents as well as its identity; mismatch is StaleGeneration.
    [[nodiscard]] core::Result<void> validate_for(const PgGraph&) const noexcept;
    /// Borrow provenance in raw r/d/paired/c order followed by anchor R then P.
    [[nodiscard]] core::Result<std::span<const PgCoefficientSource>> columns() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgCoefficientSource>> columns() const&& = delete;
    /// Number of original columns, including original zero columns; no rank is returned.
    [[nodiscard]] core::Result<core::Index> original_columns() const noexcept;
    /// Borrow nonzero records sorted by (column,cell,field). Missing records mean
    /// exact zero, not an omitted mandatory family or an uncomputed numerator.
    [[nodiscard]] core::Result<std::span<const PgCoefficientRecord>> records() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgCoefficientRecord>> records() const&& = delete;
    /// Borrow one numerator, including owned canonical zero for an absent record.
    /// Invalid indices/field or a moved-from owner return InvalidInput.
    [[nodiscard]] core::Result<const ExactCellPolynomial*>
        coefficient(PgCoefficientColumn, PgCellIndex, PgCoefficientField) const& noexcept;
    [[nodiscard]] core::Result<const ExactCellPolynomial*>
        coefficient(PgCoefficientColumn, PgCellIndex, PgCoefficientField) const&& = delete;
    /// Borrow g on one half-cell. It is stored once per cell, never per column.
    [[nodiscard]] core::Result<const ExactCellPolynomial*> g(PgCellIndex) const& noexcept;
    [[nodiscard]] core::Result<const ExactCellPolynomial*> g(PgCellIndex) const&& = delete;
    /// Borrow exact sparse coordinates of all 33 required/translation fields.
    [[nodiscard]] core::Result<std::span<const PgRequiredMapEntry>> required_map() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgRequiredMapEntry>> required_map() const&& = delete;
    /// Borrow exact C-anchor inclusion coordinates; no new C columns are appended.
    [[nodiscard]] core::Result<std::span<const PgAnchorCMapEntry>> anchor_c_map() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgAnchorCMapEntry>> anchor_c_map() const&& = delete;
    /// Borrow owned diagnostics, including after a move. This metadata remains
    /// readable without asserting valid coefficient ownership or independent rank.
    [[nodiscard]] const PgCoefficientPreparation& preparation() const& noexcept;
    [[nodiscard]] const PgCoefficientPreparation& preparation() const&& = delete;

private:
    friend struct detail::PgCoefficientFactory;
    explicit PgCoefficientDictionary(core::PlannedBuffer<detail::PgCoefficientData>,
                                     PgCoefficientPreparation) noexcept;
    core::PlannedBuffer<detail::PgCoefficientData> data_;
    PgCoefficientPreparation info_;
};
} // namespace uvlas::structure
