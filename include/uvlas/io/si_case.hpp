#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <stop_token>
#include <string_view>
#include <utility>
#include <uvlas/io/case_schema.hpp>
#include <uvlas/structure/exact_polynomial.hpp>
#include <variant>

namespace uvlas::io {
namespace detail {
struct SiFactory;
struct SiCatalogData;
struct SiRequestData;
struct SiData;
struct SiKeyOwner;
} // namespace detail

/// Length units accepted by the input schema, in their declared physical meaning.
enum class SiLengthUnit : std::uint8_t { Metre, Millimetre, Centimetre, Kilometre, Inch, Foot };
/// Mass units; PoundMass denotes lbm, not force.
enum class SiMassUnit : std::uint8_t { Kilogram, Gram, Tonne, PoundMass };
/// Time units used by the source, including event and interpolation coordinates.
enum class SiTimeUnit : std::uint8_t { Second, Millisecond, Minute };
/// Angles remain a distinct dimension; degrees require an enclosed pi factor.
enum class SiAngleUnit : std::uint8_t { Radian, Degree };
/// Original root units, retained alongside each exact conversion recipe.
struct SiUnits {
    SiLengthUnit length; ///< Source unit of length.
    SiMassUnit mass;     ///< Source unit of mass.
    SiTimeUnit time;     ///< Source unit of time.
    SiAngleUnit angle;   ///< Source unit of angles, excluding strains and rotation matrices.
    bool operator==(const SiUnits&) const = default;
};

/// Independent host caps. Zero permits none; Count and unknown values are rejected.
enum class SiQuota : std::uint8_t {
    SourceBytes,         ///< Original document bytes retained for ownership and diagnostics.
    InputNodes,          ///< Nodes in the owned source tree.
    InputMembers,        ///< Source object members, including nested objects.
    InputDepth,          ///< Maximum source nesting depth.
    SourceDigits,        ///< Total normalized exact significand digits retained from the source.
    Quantities,          ///< Logical numeric values; a rational is counted once.
    RequestTokenBytes,   ///< Copied exact request spellings, including exponent characters.
    RequestDigits,       ///< Decimal digits inspected in request tokens, including exponents.
    ExponentMagnitude,   ///< Magnitude of each normalized decimal exponent.
    UnitPower,           ///< Absolute dimension exponent used in unit conversion.
    OwnedBytes,          ///< Immutable payload capacity returned by this call.
    ScratchBytes,        ///< Simultaneously live temporary payload, including fixed integer frames.
    CumulativeBytes,     ///< Successful payload allocations, without refunds.
    Work,                ///< Charged primitives and conservative opaque exact-operation bounds.
    DiagnosticPathBytes, ///< Escaped JSON Pointer bytes retained on failure.
    Count                ///< Sentinel; not a usable quota.
};
/// Explicit SI-call limits, including capped exact algebra and cancellation.
class SiLimits {
public:
    /// Supply every cap and the external exact-memory policy. No implicit tolerances exist.
    SiLimits(std::array<std::uint64_t, static_cast<std::size_t>(SiQuota::Count)> quotas,
             structure::ExactLimits exact, std::stop_token stop = {}) noexcept
        : quotas_(quotas), exact_(exact), stop_(std::move(stop)) {}
    /// Read a named cap; invalid enum values return InvalidInput.
    [[nodiscard]] core::Result<std::uint64_t> get(SiQuota) const noexcept;
    /// Read the caller's exact-rational preparation limits.
    [[nodiscard]] const structure::ExactLimits& exact() const& noexcept {
        return exact_;
    }
    [[nodiscard]] const structure::ExactLimits& exact() const&& = delete;
    /// Cancellation is observed before individual visits and numerical operations.
    [[nodiscard]] std::stop_token stop() const noexcept {
        return stop_;
    }

private:
    friend struct detail::SiFactory;
    std::array<std::uint64_t, static_cast<std::size_t>(SiQuota::Count)> quotas_;
    structure::ExactLimits exact_;
    std::stop_token stop_;
};

/// Factory-issued identity of one original logical numeric field.
/// Equal values at different source fields have different keys; keys cannot be fabricated.
class SiQuantityKey {
public:
    /// Read the qualified source/schema identity, independent of source lifetime.
    [[nodiscard]] SchemaIdentity source() const noexcept {
        return source_;
    }
    /// Compare complete source, field, component and role identities.
    [[nodiscard]] bool operator==(const SiQuantityKey&) const noexcept;

private:
    friend struct detail::SiFactory;
    SiQuantityKey(std::shared_ptr<const detail::SiKeyOwner> owner, SchemaIdentity source,
                  std::uint64_t node, std::uint16_t type, std::uint16_t role) noexcept
        : owner_(std::move(owner)), source_(source), node_(node), type_(type), role_(role) {}
    std::shared_ptr<const detail::SiKeyOwner> owner_;
    SchemaIdentity source_;
    std::uint64_t node_;
    std::uint16_t type_, role_;
};

/// Borrowed original decimal. Numeric value is sign * digits * 10^exponent10.
struct SiDecimalView {
    /// Bind every exact component; omitted values are not inferred.
    SiDecimalView(int value_sign, std::string_view value_digits,
                  std::int64_t value_exponent) noexcept
        : sign(value_sign), digits(value_digits), exponent10(value_exponent) {}
    int sign;                ///< -1, zero or +1. Signed-zero spelling is retained separately.
    std::string_view digits; ///< All normalized significand digits, with no floating conversion.
    std::int64_t exponent10; ///< Exact signed base-10 exponent.
};
/// Borrowed original rational parts; both are exact integers, not separately converted quantities.
struct SiRationalView {
    SiDecimalView numerator;   ///< Signed exact numerator.
    SiDecimalView denominator; ///< Positive exact denominator.
};
/// Exact source and its SI multiplication recipe. Views borrow the owning catalog or SiCase.
struct SiRecipeView {
    std::variant<SiDecimalView, SiRationalView> value; ///< One logical decimal or rational value.
    std::string_view original;  ///< Original decimal token when present; otherwise empty.
    SiUnits units;              ///< Source root units; no source lifetime is required by the owner.
    SchemaDimension dimensions; ///< Length, mass, time and angle exponents.
    std::int16_t pi_power;      ///< Angle exponent for degrees, otherwise zero.
};
/// Metadata for one conversion, preserving its source field and component meaning.
struct SiCatalogEntry {
    SiQuantityKey key;       ///< Stable source-bound key to use in a request row.
    SiRecipeView recipe;     ///< Borrowed exact value and multiplication recipe.
    std::string_view record; ///< Static name of the enclosing typed record.
    std::string_view field;  ///< Static name of its enclosing field; array order is preserved.
};

/// Successful owned-payload accounting; excludes allocator bookkeeping and process RSS.
struct SiMemoryReport {
    std::uint64_t owned_bytes;      ///< Immutable payload including copied exact input.
    std::uint64_t scratch_bytes;    ///< Maximum simultaneous planned scratch payload.
    std::uint64_t cumulative_bytes; ///< Successful allocations charged during this call.
    std::uint64_t work; ///< Charged primitives and conservative opaque exact-operation bounds.
    std::uint64_t planned_owned_bytes;   ///< Immutable capacity admitted before allocation.
    std::uint64_t planned_scratch_bytes; ///< Whole-call simultaneous scratch capacity.
    std::uint64_t planned_work;    ///< Source-derived upper bound admitted for the entire call.
    std::uint64_t admission_bytes; ///< Whole-call cumulative capacity including diagnostic reserve.
};

/// Owned complete inventory of the case's physical and dimensionless logical numbers.
class SiCatalog {
public:
    SiCatalog(SiCatalog&&) noexcept;            ///< Transfer ownership, leaving the source absent.
    SiCatalog& operator=(SiCatalog&&) noexcept; ///< Replace storage and invalidate old views.
    SiCatalog(const SiCatalog&) = delete;
    SiCatalog& operator=(const SiCatalog&) = delete;
    ~SiCatalog(); ///< Release payload; successful plan charges are not refunded.
    [[nodiscard]] core::Result<std::uint64_t> size() const noexcept; ///< Logical numeric count.
    [[nodiscard]] core::Result<SiCatalogEntry>
    entry(std::uint64_t) const& noexcept; ///< Checked row.
    [[nodiscard]] core::Result<SiCatalogEntry> entry(std::uint64_t) const&& = delete;
    [[nodiscard]] core::Result<SiMemoryReport>
    memory_report() const noexcept; ///< Actual accounting.
private:
    friend struct detail::SiFactory;
    explicit SiCatalog(core::PlannedBuffer<detail::SiCatalogData>) noexcept;
    core::PlannedBuffer<detail::SiCatalogData> data_;
};

/// A finite exact decimal supplied to this synchronous construction call.
struct SiDecimalInput {
    std::string_view token; ///< Exact decimal spelling, including its optional exponent.
}; ///< No whitespace or nonfinite spellings.
/// Exact signed decimal integers; the denominator must be positive.
struct SiRationalInput {
    std::string_view numerator;   ///< Signed exact decimal integer.
    std::string_view denominator; ///< Strictly positive exact decimal integer.
};
/// Explicit exact numeric request input; constructors copy its bytes before returning.
using SiExactInput = std::variant<SiDecimalInput, SiRationalInput>;
/// Scale is already in SI, with explicit physical dimensions; it is never converted again.
struct SiScaleInput {
    SiExactInput value;         ///< Strictly positive scale.
    SchemaDimension dimensions; ///< Must agree with the catalog quantity.
};
/// Requirement E <= absolute * scale + relative * abs(true SI value).
struct SiRequirementRow {
    SiQuantityKey key;     ///< Every catalog key occurs exactly once.
    SiScaleInput scale;    ///< Positive exact SI scale.
    SiExactInput absolute; ///< Nonnegative dimensionless absolute coefficient.
    SiExactInput relative; ///< Nonnegative dimensionless relative coefficient.
};
/// Borrowed exact requirements retained by an owning SI result or diagnostic.
struct SiRequirementView {
    SiScaleInput scale;    ///< Original positive SI scale and its dimensions.
    SiExactInput absolute; ///< Original nonnegative dimensionless absolute coefficient.
    SiExactInput relative; ///< Original nonnegative dimensionless relative coefficient.
};

/// Owned, complete and validated requirements, independent of caller rows/catalog lifetime.
class SiRequest {
public:
    SiRequest(SiRequest&&) noexcept;            ///< Transfer ownership, leaving an absent source.
    SiRequest& operator=(SiRequest&&) noexcept; ///< Replace ownership without throwing.
    SiRequest(const SiRequest&) = delete;
    SiRequest& operator=(const SiRequest&) = delete;
    ~SiRequest(); ///< Release payload without refunding plan charges.
    [[nodiscard]] core::Result<SiMemoryReport> memory_report() const noexcept; ///< Call accounting.
private:
    friend struct detail::SiFactory;
    explicit SiRequest(core::PlannedBuffer<detail::SiRequestData>) noexcept;
    core::PlannedBuffer<detail::SiRequestData> data_;
};

/// A normalized center/radius multiplied by 2^exponent2; no implicit projection exists.
struct ScaledBinary64 {
    double center;          ///< Finite normalized center; zero is represented at exponent zero.
    double radius;          ///< Nonnegative outward error in the same scaled coordinates.
    std::int64_t exponent2; ///< Must accompany the center and radius to every consumer.
};
/// Two independently retained words and an outward error, multiplied by 2^exponent2.
struct ScaledTwoWord {
    double high;            ///< Nearest-even high word.
    double low;             ///< Nearest-even exact residual word; never silently discarded.
    double radius;          ///< Outward bound including all input and conversion rounding.
    std::int64_t exponent2; ///< Scale exponent; TwoWord extends precision, not double range.
};
/// The selected center representation is explicit in the type.
using SiEnclosure = std::variant<ScaledBinary64, ScaledTwoWord>;
/// Qualified read-only converted quantity, borrowing the owning SiCase.
struct SiQuantityView {
    SiQuantityKey key;     ///< One uncertainty origin, shared by repeated references to this field.
    SiRecipeView recipe;   ///< Exact original value remains available after conversion.
    SiEnclosure enclosure; ///< Total conversion enclosure; pi is already included once.
    std::uint64_t ledger_origin;   ///< Distinct source fields remain distinct uncertainty origins.
    SiRequirementView requirement; ///< Exact request borrowed from this result, not the caller.
};
/// Structural absence resolved without inventing a physical zero state or a numerical budget.
struct SiResolvedDefaults {
    bool provenance_missing; ///< Supplemental provenance was absent and resolves to empty.
    bool initial_missing;    ///< Original absence, interpreted according to the existing analysis.
    bool restart_missing;    ///< Original absence; no checkpoint or initial state was invented.
};

/// SI-only qualified owning result. Physical preparation and QoI error propagation remain separate.
class SiCase {
public:
    SiCase(SiCase&&) noexcept;            ///< Transfer ownership; source becomes absent.
    SiCase& operator=(SiCase&&) noexcept; ///< Replace payload and invalidate prior borrowed views.
    SiCase(const SiCase&) = delete;
    SiCase& operator=(const SiCase&) = delete;
    ~SiCase(); ///< Release payload without refunding successful plan charges.
    /// Borrow the original-unit typed case owned by this result, not the caller's case.
    [[nodiscard]] core::Result<std::reference_wrapper<const SchemaCase>>
    source_case() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const SchemaCase>>
    source_case() const&& = delete;
    /// Convert only a numeric field obtained from this result's source_case().
    /// Foreign owners are rejected even when source identity and ordinals agree.
    template <detail::FieldType Type>
    [[nodiscard]] core::Result<SiQuantityView> quantity(schema::Field<Type> field) const& noexcept
        requires(detail::field_kinds[static_cast<std::size_t>(Type)] == detail::FieldKind::Number ||
                 detail::field_kinds[static_cast<std::size_t>(Type)] ==
                     detail::FieldKind::ExactQuantity ||
                 Type == detail::FieldType::T62)
    {
        return quantity_impl(field.owner_, field.node_);
    }
    template <detail::FieldType Type>
    core::Result<SiQuantityView> quantity(schema::Field<Type>) const&& = delete;
    [[nodiscard]] core::Result<SiResolvedDefaults>
    defaults() const noexcept; ///< Original presence.
    [[nodiscard]] core::Result<SiMemoryReport> memory_report() const noexcept; ///< Call accounting.
private:
    friend struct detail::SiFactory;
    explicit SiCase(core::PlannedBuffer<detail::SiData>) noexcept;
    core::Result<SiQuantityView> quantity_impl(const detail::SchemaData*,
                                               std::uint64_t) const noexcept;
    core::PlannedBuffer<detail::SiData> data_;
};

/// Precise SI rejection reason, supplementing the unchanged core Failure.
enum class SiReason : std::uint8_t {
    InvalidOwner,
    IdentityMismatch,
    InvalidRequirement,
    MissingRequirement,
    DuplicateRequirement,
    DimensionMismatch,
    InvalidExactInput,
    QuotaExceeded,
    Overflow,
    Allocation,
    Cancelled,
    UnsupportedMemory,
    UnresolvedPrecision,
    IrrationalExactOutput,
    UnexpectedException
};
/// Completeness of the diagnostic's owned exact requirements, independent of its path.
enum class SiRequestState : std::uint8_t { NotAvailable, Complete, Quota, Allocation };
/// Fixed owned diagnostic survives failure-path allocation failure.
struct SiDiagnosticDetail {
    SiReason reason;                      ///< Exact reason, independent of borrowed exception text.
    std::optional<SourceIdentity> source; ///< Original source identity when known.
    std::optional<SiQuantityKey> key;     ///< Logical field responsible, when known.
    std::optional<SiQuota> quota;         ///< Insufficient explicit cap, when applicable.
    std::optional<std::uint64_t> required, limit; ///< Exact count demand and configured bound.
    std::optional<SchemaDimension> dimensions;    ///< Physical dimension of the affected quantity.
    std::optional<SiEnclosure> achieved; ///< Last enclosure, when precision was unresolved.
    std::uint64_t work;                  ///< Work actually charged before failure.
    SiRequestState request_state;        ///< Whether the exact requirements were completely copied.
    SchemaPathState path_state; ///< Complete JSON Pointer or a named failure to construct it.
};
/// Caller-owned move-only diagnostic slot. Reset at entry, empty after success.
class SiDiagnostic {
public:
    SiDiagnostic() noexcept = default;     ///< Construct an empty slot without allocation.
    SiDiagnostic(SiDiagnostic&&) noexcept; ///< Transfer ownership and leave the source empty.
    SiDiagnostic& operator=(SiDiagnostic&&) noexcept; ///< Replace this slot and empty the source.
    SiDiagnostic(const SiDiagnostic&) = delete;
    SiDiagnostic& operator=(const SiDiagnostic&) = delete;
    [[nodiscard]] std::optional<SiDiagnosticDetail> detail() const noexcept {
        return detail_;
    }
    /// Borrow a complete path only; partial paths are never presented as complete.
    [[nodiscard]] std::optional<std::string_view> path() const& noexcept;
    [[nodiscard]] std::optional<std::string_view> path() const&& = delete;
    /// Borrow complete owned requirements, with no caller-row or input lifetime dependency.
    [[nodiscard]] std::optional<SiRequirementView> requirement() const& noexcept;
    [[nodiscard]] std::optional<SiRequirementView> requirement() const&& = delete;
    void reset() noexcept; ///< Discard this slot; no source/result/plan is changed.
private:
    friend struct detail::SiFactory;
    std::optional<SiDiagnosticDetail> detail_;
    core::PlannedBuffer<char> path_;
    core::PlannedBuffer<char> requirement_text_;
    std::optional<SiRequirementView> requirement_;
};

/// Own every logical exact number and its unit/role; no source lifetime is borrowed.
/// Explicit limits and optional diagnostics do not change the input on failure.
[[nodiscard]] core::Result<SiCatalog> make_si_catalog(const SchemaCase&, const SiLimits&,
                                                      core::ResourcePlan&, SiDiagnostic*) noexcept;
/// Copy and validate exactly one scale/absolute/relative row for every catalog key.
/// Scales are already SI; absolute/relative coefficients are dimensionless and nonnegative.
[[nodiscard]] core::Result<SiRequest> make_si_request(const SiCatalog&,
                                                      std::span<const SiRequirementRow>,
                                                      const SiLimits&, core::ResourcePlan&,
                                                      SiDiagnostic*) noexcept;
/// Qualify every original number before returning an owning SI-only result.
/// Never overlap calls sharing a mutable plan/slot. Independent calls may run concurrently.
/// Resource, cancellation and precision failures preserve inputs and all prior results.
[[nodiscard]] core::Result<SiCase> convert_case_to_si(const SchemaCase&, const SiCatalog&,
                                                      const SiRequest&, const SiLimits&,
                                                      core::ResourcePlan&, SiDiagnostic*) noexcept;
/// Materialize only rational recipes under explicit exact/working caps.
/// Nonzero pi powers return PrecisionUnsupported; exact zero always denotes 0/1.
[[nodiscard]] core::Result<structure::ExactRational>
materialize_rational(const SiRecipeView&, const structure::ExactLimits&,
                     core::ResourcePlan&) noexcept;

} // namespace uvlas::io
