#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string_view>
#include <uvlas/core/resource_plan.hpp>
#include <uvlas/io/case_schema_records.hpp>

namespace uvlas::io {
namespace detail {
struct SchemaContext;
struct SiFactory;
}

/// Logical limits are independent; zero permits none of the corresponding item.
/// These ordinals are an in-process API and must not be persisted.
enum class SchemaQuota : std::uint8_t {
    SourceBytes,         ///< Original UTF-8 source bytes.
    InputDepth,          ///< Original container depth, including the root.
    InputNodes,          ///< All source value nodes, excluding member keys.
    InputMembers,        ///< All source object members.
    DecodedStringBytes,  ///< One decoded key or string, in UTF-8 bytes.
    DecimalTokenBytes,   ///< One original decimal token, in ASCII bytes.
    TotalTextBytes,      ///< All source decoded text and decimal tokens.
    TypedBytes,          ///< All additional owned snapshot and typed payload bytes.
    Entities,            ///< Definitions in all typed identifier namespaces.
    References,          ///< Every internal reference occurrence.
    Scalars,             ///< All scalar occurrences, including IDs, enum tags and rational parts.
    PolynomialTerms,     ///< Terms in all polynomial inputs.
    TableSamples,        ///< Time/value rows or polynomial coefficient rows.
    ExpressionNodes,     ///< Nodes summed over all expression DAGs.
    ExpressionEdges,     ///< Operand occurrences summed over all DAGs.
    ExpressionDepth,     ///< Largest expression depth, with leaves at one.
    ExpressionPower,     ///< Largest absolute integer power.
    ExpressionPoints,    ///< Evaluation-point cap retained for later numerical preparation.
    DiagnosticPathBytes, ///< Maximum complete JSON Pointer length.
    ScratchBytes,        ///< All simultaneous temporary payload bytes.
    Work,                ///< Total elementary visits, digit operations and byte copies.
    ExactDigits,         ///< Largest exact temporary significand length.
    ExponentMagnitude,   ///< Largest absolute decimal exponent permitted.
    Count                ///< Number of required limit values; not a quota.
};

/// Explicit limits supplied by the host. A case declaration can only lower them.
class SchemaLimits {
public:
    /// Supply every logical limit; no default or zero-as-unlimited convention exists.
    explicit SchemaLimits(
        std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> values) noexcept
        : values_(values) {}
    /// Return the configured quota, or InvalidInput for Count or an unknown enum value.
    [[nodiscard]] core::Result<std::uint64_t> get(SchemaQuota quota) const noexcept {
        const auto index = static_cast<std::size_t>(quota);
        if (index >= values_.size())
            return std::unexpected(
                core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
        return values_[index];
    }

private:
    friend struct detail::SchemaContext;
    std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> values_;
};

/// Machine-readable schema rejection, independent of exception text and file lifetime.
enum class SchemaReason : std::uint8_t {
    MissingField,        ///< A required field is absent.
    UnknownField,        ///< A field is not allowed by the selected record.
    WrongKind,           ///< The JSON kind is incompatible with the field.
    WrongShape,          ///< Array dimensions or required components disagree.
    UnsupportedTag,      ///< A schema or closed enumeration value is unsupported.
    InvalidValue,        ///< An exact range, order or mutual condition is violated.
    InvalidIdentifier,   ///< Identifier spelling or range is invalid.
    DuplicateIdentifier, ///< A namespace contains two different definitions.
    MissingReference,    ///< A typed reference has no compatible definition.
    ReferenceCycle,      ///< A prohibited dependency contains a cycle.
    DimensionMismatch,   ///< Expression or consumer unit dimensions disagree.
    InvalidExpression,   ///< An operand, output, reachability or domain shape is invalid.
    QuotaExceeded,       ///< An explicit quota is insufficient.
    Overflow,            ///< A checked count, index, exponent or dimension is not representable.
    Allocation,          ///< Storage could not be allocated.
    UnexpectedException, ///< An unexpected exception was translated at the entry.
    InvalidOwner         ///< The supplied document is absent or its identity cannot be read.
};

/// Path status never presents a truncated or missing path as a complete JSON Pointer.
enum class SchemaPathState : std::uint8_t {
    Complete,  ///< The diagnostic owns the complete path.
    NotBuilt,  ///< Failure occurred before a complete path could be built.
    PathQuota, ///< Complete path exceeds the explicit quota.
    Allocation ///< Complete path storage could not be allocated.
};

/// Typed identity namespace. Equal numeric IDs in different domains remain distinct.
enum class SchemaNamespace : std::uint8_t {
    Geometry,
    Material,
    Shell,
    Beam,
    Joint,
    Support,
    Mass,
    AeroSurface,
    Layer,
    Load,
    Observation,
    Source,
    Frame,
    Branch,
    WakeRow,
    WakePoint,
    WakeBand
};

/// Owned fixed diagnostic context; absent positions and demands remain unknown.
struct SchemaDiagnosticDetail {
    SchemaReason reason; ///< Precise schema-level reason supplementing core Failure.
    std::optional<SourceIdentity> source;  ///< Original source identity, when available.
    std::optional<std::uint64_t> node;     ///< Actual document-local node ordinal, when known.
    std::optional<SchemaQuota> quota;      ///< Quota causing rejection, when applicable.
    std::optional<std::uint64_t> required; ///< Exactly computed demand; absent on overflow.
    std::optional<std::uint64_t> limit;    ///< Explicit active bound, when applicable.
    std::optional<detail::FieldType> expected_type; ///< Owned closed type ID; its spelling is
                                                    ///< supplied by the static type formatter.
    std::optional<std::uint16_t>
        expected_field; ///< Owned closed field-name ID; never a borrowed input string.
    std::optional<std::uint16_t>
        expected_variant; ///< Owned closed record-variant ID, when a branch was selected.
    std::optional<InputNodeKind> actual_kind; ///< Observed source kind, when a real node is known.
    std::optional<SchemaNamespace>
        reference_domain; ///< Required target namespace, when resolving a reference.
    std::optional<std::uint64_t>
        reference_id;           ///< Actual user-supplied reference ID, including zero.
    std::uint64_t work;         ///< Elementary work actually performed before failure.
    SchemaPathState path_state; ///< Whether the owned path is complete.
    /// Borrow the expected declared type's static name; no input lifetime is required.
    [[nodiscard]] std::optional<std::string_view> expected_type_name() const noexcept;
    /// Borrow the required field's static name when its selected record is known.
    [[nodiscard]] std::optional<std::string_view> expected_field_name() const noexcept;
    /// Borrow the selected record variant's static name when known.
    [[nodiscard]] std::optional<std::string_view> expected_variant_name() const noexcept;
};

/// Caller-owned move-only diagnostic slot. Reset at entry, empty on success.
/// It owns all data and remains usable after the input, result and work storage die.
class SchemaDiagnostic {
public:
    /// Start empty without allocating.
    SchemaDiagnostic() noexcept = default;
    /// Transfer diagnostic ownership; the source becomes empty.
    SchemaDiagnostic(SchemaDiagnostic&&) noexcept;
    /// Replace this diagnostic with a nonthrowing ownership transfer.
    SchemaDiagnostic& operator=(SchemaDiagnostic&&) noexcept;
    SchemaDiagnostic(const SchemaDiagnostic&) = delete;
    SchemaDiagnostic& operator=(const SchemaDiagnostic&) = delete;
    /// Release owned path; past ResourcePlan charges are not refunded.
    ~SchemaDiagnostic() = default;
    /// Return fixed failure detail by value, or no value following success.
    [[nodiscard]] std::optional<SchemaDiagnosticDetail> detail() const noexcept;
    /// Borrow a complete path only; an unavailable path returns no value.
    [[nodiscard]] std::optional<std::string_view> path() const& noexcept;
    [[nodiscard]] std::optional<std::string_view> path() const&& = delete;
    /// Discard this slot only, without modifying any stored Failure copy.
    void reset() noexcept;

private:
    friend struct detail::SchemaFactory;
    std::optional<SchemaDiagnosticDetail> detail_;
    core::PlannedBuffer<char> path_;
};

/// Source-derived accounting, not allocator/RSS or elapsed-time guarantees.
struct SchemaMemoryReport {
    std::uint64_t owned_bytes;   ///< Additional immutable payload including the deep source copy.
    std::uint64_t scratch_bytes; ///< Maximum simultaneously live scratch payload.
    std::uint64_t cumulative_bytes;      ///< Total successful allocations, including scratch.
    std::uint64_t work;                  ///< Actual elementary operations charged to this call.
    std::uint64_t planned_owned_bytes;   ///< Pre-admitted immutable payload capacity.
    std::uint64_t planned_scratch_bytes; ///< Pre-admitted simultaneous temporary capacity.
    std::uint64_t admission_bytes;       ///< Whole-call cumulative admission, including complete
                                         ///< failure-path reserve.
};

/// Input identity and validator representation, distinct from a physical Model fingerprint.
struct SchemaIdentity {
    SourceIdentity
        source; ///< Exact original source identity; equal values can have different digests.
    std::uint32_t contract_revision;       ///< Closed schema contract revision used for validation.
    std::uint32_t representation_revision; ///< Typed payload representation revision.
};

/// Immutable fully checked case, owning its source, typed records, exact values and references.
/// No numerical SI, formula-domain, geometric or practical qualification is implied.
class SchemaCase {
public:
    /// Transfer storage; the source becomes absent.
    SchemaCase(SchemaCase&&) noexcept;
    /// Replace owned storage and invalidate all earlier borrowed fields.
    SchemaCase& operator=(SchemaCase&&) noexcept;
    SchemaCase(const SchemaCase&) = delete;
    SchemaCase& operator=(const SchemaCase&) = delete;
    /// Release storage without refunding preparation charges.
    ~SchemaCase();
    /// Borrow the named root record. An absent/moved-from owner returns InvalidInput.
    [[nodiscard]] core::Result<schema::Case> root() const& noexcept;
    [[nodiscard]] core::Result<schema::Case> root() const&& = delete;
    /// Return the original source identity; an absent owner returns InvalidInput.
    [[nodiscard]] core::Result<SourceIdentity> source_identity() const noexcept;
    /// Return source and schema/representation revisions; an absent owner is invalid.
    [[nodiscard]] core::Result<SchemaIdentity> identity() const noexcept;
    /// Return actual accounting for the successful call; an absent owner is invalid.
    [[nodiscard]] core::Result<SchemaMemoryReport> memory_report() const noexcept;

private:
    friend struct detail::SchemaFactory;
    explicit SchemaCase(core::PlannedBuffer<detail::SchemaData>) noexcept;
    friend struct detail::SiFactory;
    core::PlannedBuffer<detail::SchemaData> data_;
};

/// Validate the complete closed schema, preserving exact decimals and an owned source copy.
/// Input is read-only; do not overlap or concurrently share the mutable plan/diagnostic.
/// Unknown tags are Unsupported; invalid fields/references are InvalidInput; insufficient
/// quotas, overflow and allocation failure are ResourceExceeded. No exception escapes.
/// The diagnostic may be null, without changing acceptance or the core failure category.
/// Failure preserves the input and all old results; successful plan charges are not refunded.
/// A schema success cannot substitute for numerical preparation or full CLI validation.
[[nodiscard]] core::Result<SchemaCase> validate_case_schema(const InputDocument& input,
                                                            SchemaLimits limits,
                                                            core::ResourcePlan& resources,
                                                            SchemaDiagnostic* diagnostic) noexcept;

} // namespace uvlas::io
