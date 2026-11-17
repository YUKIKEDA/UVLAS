#pragma once

#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <uvlas/core/identifier.hpp>
#include <uvlas/core/resource_plan.hpp>
#include <uvlas/io/input_scalar.hpp>

namespace uvlas::io {
namespace detail {
struct InputSourceTag;
struct InputDocumentData;
struct InputDocumentFactory;
} // namespace detail

/// Caller-assigned source identity; zero is valid. It is neither a physical ID
/// nor proof that two sources have equal bytes.
using SourceId = core::Identifier<detail::InputSourceTag>;
/// SHA-256 of the exact source bytes, independent of path and serialization.
struct SourceDigest {
    std::array<std::byte, 32> bytes; ///< Digest bytes in standard SHA-256 output order.
    /// Compare all digest bytes exactly.
    bool operator==(const SourceDigest&) const = default;
};
/// Owned reporting identity of one immutable source snapshot.
struct SourceIdentity {
    SourceId source;          ///< Caller-assigned source identity.
    std::uint64_t byte_count; ///< Exact original byte length; no added terminator.
    SourceDigest digest;      ///< Hash of those exact bytes.
    /// Compare source identity, length and digest; no model qualification is implied.
    bool operator==(const SourceIdentity&) const = default;
};
/// Select explicit forecast admission or require a genuine external bound.
enum class ParserAllocationMode : std::uint8_t {
    EstimateAdmission, ///< Admit the declared forecast; not an allocator/RSS guarantee.
    RequireProvenBound ///< Require a proven bound; unavailable models return Unsupported.
};
/// External parser allocation policy; zero bytes is a real zero cap.
struct ParserAllocationPolicy {
    ParserAllocationMode mode;   ///< Explicit mutually exclusive admission request.
    std::uint64_t maximum_bytes; ///< Maximum acceptable forecast or proven bound.
};
/// Explicit logical quotas, checked before writes or narrowing conversions.
/// Root depth is one; text totals exclude the separately counted source bytes.
struct InputLimits {
    /// Supply all seven logical limits and the external-memory request explicitly.
    InputLimits(std::uint64_t source_bytes, std::uint64_t depth, std::uint64_t nodes,
                std::uint64_t members, std::uint64_t decoded_string_bytes,
                std::uint64_t decimal_token_bytes, std::uint64_t total_text_bytes,
                ParserAllocationPolicy parser_memory) noexcept;
    std::uint64_t maximum_source_bytes;         ///< Maximum original source length.
    std::uint64_t maximum_depth;                ///< Maximum container depth, root object=1.
    std::uint64_t maximum_nodes;                ///< All values including root; keys excluded.
    std::uint64_t maximum_members;              ///< Sum of members over all objects.
    std::uint64_t maximum_decoded_string_bytes; ///< UTF-8 bytes in one decoded key or String.
    std::uint64_t maximum_decimal_token_bytes;  ///< ASCII bytes in one float-event token.
    std::uint64_t maximum_total_text_bytes;     ///< Sum of owned keys, Strings and decimal tokens.
    ParserAllocationPolicy external_memory;     ///< Separate parser-library allocation policy.
};
/// Fixed source-model family; an unrecognized implementation is Unsupported.
enum class ParserLibraryFamily : std::uint8_t {
    MsvcStl,  ///< Microsoft STL and its separately identified iterator-debug setting.
    Libstdcxx ///< GNU libstdc++ with the identified revision and ABI.
};
/// Reporting metadata of the source-derived external memory model.
struct ParserEnvironment {
    ParserLibraryFamily family;              ///< Matched implementation family.
    std::uint64_t standard_library_revision; ///< Vendor revision inspected by the factory.
    std::uint32_t model_revision;            ///< Version of the source-derived forecast model.
    std::uint32_t parser_major;              ///< Upstream parser major version.
    std::uint32_t parser_minor;              ///< Upstream parser minor version.
    std::uint32_t parser_patch;              ///< Upstream parser patch version.
    std::uint32_t parser_repair;             ///< Identified allocation-exception source repair.
    std::uint32_t pointer_bits;              ///< Pointer width of the matched ABI.
    std::uint32_t size_bits;                 ///< Width of size_t in that ABI.
    std::uint32_t word_bits;                 ///< Parser vector<bool> storage word width.
    std::uint32_t char_bits;                 ///< Bits per source byte; only eight is supported.
    std::uint32_t iterator_debug_level;      ///< Matched STL iterator-debug configuration.
};
/// Roles distinguish the fixed dependency sources used to justify a memory model.
enum class ParserSourceRole : std::uint8_t {
    Parser,           ///< JSON SAX grammar and state stack.
    Lexer,            ///< JSON scanner, token buffers and number conversion boundary.
    Exceptions,       ///< JSON exception creation and message ownership.
    StringConcat,     ///< JSON message concatenation and copy-return implementation.
    Vector,           ///< Standard vector capacity/constructor implementation.
    VectorGrowth,     ///< Separate standard vector growth implementation, when present.
    BitVector,        ///< Standard packed-bool vector representation.
    String,           ///< Standard string growth and storage implementation.
    CowString,        ///< GNU exception message's legacy COW string representation.
    ExceptionRuntime, ///< Runtime message copying and allocation behavior.
    ExceptionHeader,  ///< Runtime exception/fallback interface.
    IteratorStorage   ///< Standard iterator-proxy allocation and storage.
};
/// Owned role and source hash; not a caller-supplied proof capability.
struct ParserSourceEvidence {
    ParserSourceRole role; ///< Function of this inspected dependency source.
    SourceDigest sha256;   ///< SHA-256 of the source used to justify the model.
};
/// Simultaneous payload forecasts, including old/new container capacities.
/// Syntax/range alternatives are combined by maximum, not added as live heaps.
struct ParserMemoryBreakdown {
    std::uint64_t token_vector_bytes; ///< Lexer diagnostic token vector, including growth.
    std::uint64_t token_buffer_bytes; ///< Lexer decoded/token string, including growth.
    std::uint64_t state_stack_bytes;  ///< Packed parser state stack and its proxies.
    std::uint64_t syntax_error_bytes; ///< Simultaneous syntax-message construction forecast.
    std::uint64_t number_range_bytes; ///< Alternative number-range message forecast.
};
/// Factory-generated external forecast and its fixed source/model provenance.
/// Estimates are not observed allocation counts or an allocator/RSS hard cap.
class ParserMemoryPlan {
public:
    /// Copy the complete simultaneous external payload forecast in bytes.
    [[nodiscard]] std::uint64_t estimate_bytes() const noexcept;
    /// Return no value when a genuine external bound has not been established.
    [[nodiscard]] std::optional<std::uint64_t> proven_upper_bytes() const noexcept;
    /// Copy the implementation, ABI and model identification used for admission.
    [[nodiscard]] ParserEnvironment environment() const noexcept;
    /// Copy the named forecast components; none is an observed allocation count.
    [[nodiscard]] ParserMemoryBreakdown breakdown() const noexcept;
    /// Borrow fixed evidence owned by this plan until it moves or is destroyed.
    [[nodiscard]] std::span<const ParserSourceEvidence> sources() const& noexcept;
    [[nodiscard]] std::span<const ParserSourceEvidence> sources() const&& = delete;

private:
    friend struct detail::InputDocumentFactory;
    ParserMemoryPlan(ParserEnvironment, ParserMemoryBreakdown, std::uint64_t,
                     std::array<ParserSourceEvidence, 12>, std::size_t) noexcept;
    ParserEnvironment environment_;
    ParserMemoryBreakdown breakdown_;
    std::uint64_t estimate_;
    std::array<ParserSourceEvidence, 12> sources_;
    std::size_t source_count_;
};

/// Logical quota named by an owned input diagnostic.
enum class InputQuota : std::uint8_t {
    SourceBytes,        ///< Original source length.
    Depth,              ///< Active container depth.
    Nodes,              ///< Total value-node count.
    Members,            ///< Total object-member count.
    DecodedStringBytes, ///< UTF-8 length of one decoded key or string.
    DecimalTokenBytes,  ///< ASCII length of one retained decimal token.
    TotalTextBytes      ///< Sum of all retained decoded and decimal text.
};
/// Known logical quota demand; a missing demand was not measured or computed.
struct InputQuotaDiagnostic {
    InputQuota kind;                       ///< Logical count whose admission failed.
    std::uint64_t limit;                   ///< Explicit requested cap.
    std::optional<std::uint64_t> required; ///< Known demand, absent if not representable.
};
/// Fixed diagnostic reasons, independent of third-party exception lifetimes.
enum class InputDiagnosticReason : std::uint8_t {
    InvalidLimits,          ///< An admission tag or explicit setting is invalid.
    UnsupportedMemoryModel, ///< This library/ABI/configuration has no inspected model.
    UnprovenMemoryBound,    ///< The request requires an unavailable proven bound.
    ExternalMemoryLimit,    ///< The declared forecast exceeds the explicit cap.
    DimensionRange,         ///< A size/count cannot be represented or calculated safely.
    QuotaExceeded,          ///< A named logical quota failed.
    Allocation,             ///< Payload/library allocation failed or payload cap was insufficient.
    EmptyInput,             ///< No source bytes were supplied.
    ByteOrderMark,          ///< The input begins with a prohibited UTF-8 BOM.
    RawNull,                ///< A raw zero byte would violate strict whole-source parsing.
    InvalidUtf8,            ///< JSON string bytes are not valid UTF-8.
    InvalidStringEscape,    ///< A string escape or surrogate sequence is malformed.
    JsonSyntax,             ///< Other strict JSON syntax or trailing-input failure.
    DuplicateKey,           ///< Decoded keys repeat within one object.
    RootType,               ///< Root is not an object.
    SchemaType,             ///< A present schema field is not a String.
    UnsupportedSchema,      ///< Schema text names an unsupported version.
    UnknownRootKey,         ///< A root member is not in the allowed envelope.
    MissingRootKey,         ///< A mandatory root member is absent.
    NumberRange,            ///< Parser cannot produce a finite float event.
    UnexpectedException ///< Other contained standard/nonstandard exception or invariant failure.
};
/// Allocation-free owned detail accompanying core::Failure. A byte offset is
/// zero-based in the original source (EOF equals byte_count); event ordinal is
/// zero-based among callbacks actually invoked. Missing positions stay absent.
struct InputDiagnostic {
    /// Start with a reason/source and no fabricated digest, position or demand.
    InputDiagnostic(InputDiagnosticReason, SourceId) noexcept;
    InputDiagnosticReason reason;                ///< Stable machine-readable cause.
    SourceId source;                             ///< Supplied source identity, including zero.
    std::optional<SourceDigest> digest;          ///< Present only after snapshot hashing.
    std::optional<std::uint64_t> byte_offset;    ///< Known zero-based original-byte position.
    std::optional<std::uint64_t> event_ordinal;  ///< An actual callback's zero-based ordinal.
    std::optional<InputQuotaDiagnostic> quota;   ///< Logical admission detail, when relevant.
    std::optional<ParserMemoryPlan> parser_plan; ///< Produced forecast, when available.
};

/// Reporting ordinal generated by a document. No lookup accepts a caller-made
/// ordinal as a cross-document capability, and physical IDs cannot substitute it.
class InputNodeOrdinal {
public:
    /// Report the zero-based ordinal; it is not a reusable lookup capability.
    [[nodiscard]] core::Index value() const noexcept {
        return value_;
    }
    /// Compare reporting numbers without proving document identity.
    auto operator<=>(const InputNodeOrdinal&) const = default;

private:
    friend struct detail::InputDocumentFactory;
    explicit InputNodeOrdinal(core::Index value) noexcept : value_(value) {}
    core::Index value_;
};
/// All JSON value kinds; missing members and null remain distinct.
enum class InputNodeKind : std::uint8_t {
    Null,            ///< Explicit JSON null, distinct from an absent member.
    Boolean,         ///< JSON true or false.
    String,          ///< Decoded UTF-8 bytes with explicit length.
    SignedInteger,   ///< Exact native signed event, including -0 as value zero.
    UnsignedInteger, ///< Exact native unsigned event without a double round trip.
    DecimalNumber,   ///< Original float-event token, including large-integer fallback.
    Array,           ///< Ordered direct element sequence.
    Object           ///< Ordered unique decoded-key/member sequence.
};
class InputArrayView;
class InputObjectView;
/// Immutable borrow of a parsed node and its originating storage. Views expire
/// when the document owner moves, is assigned or is destroyed.
class InputNodeView {
public:
    /// Report this node's document-local ordinal.
    [[nodiscard]] InputNodeOrdinal ordinal() const noexcept;
    /// Return the retained JSON kind without numeric coercion.
    [[nodiscard]] InputNodeKind kind() const noexcept;
    /// A wrong-kind scalar query returns no value, without converting JSON types.
    [[nodiscard]] std::optional<bool> boolean() const noexcept;
    /// Return only a native signed event; other kinds return no value.
    [[nodiscard]] std::optional<std::int64_t> signed_integer() const noexcept;
    /// Return only a native unsigned event; other kinds return no value.
    [[nodiscard]] std::optional<std::uint64_t> unsigned_integer() const noexcept;
    /// Only native integer nodes can feed the strict Index converter.
    [[nodiscard]] std::optional<NativeInteger> native_integer() const noexcept;
    /// Decoded UTF-8 bytes, including escaped NUL; no normalization is performed.
    [[nodiscard]] std::optional<std::string_view> string() const& noexcept;
    [[nodiscard]] std::optional<std::string_view> string() const&& = delete;
    /// Exact original float-event token; no as_double conversion is provided.
    [[nodiscard]] std::optional<std::string_view> decimal_token() const& noexcept;
    [[nodiscard]] std::optional<std::string_view> decimal_token() const&& = delete;
    /// Borrow an Array sequence, or no value for another kind.
    [[nodiscard]] std::optional<InputArrayView> array() const& noexcept;
    [[nodiscard]] std::optional<InputArrayView> array() const&& = delete;
    /// Borrow an Object sequence, or no value for another kind.
    [[nodiscard]] std::optional<InputObjectView> object() const& noexcept;
    [[nodiscard]] std::optional<InputObjectView> object() const&& = delete;

private:
    friend struct detail::InputDocumentFactory;
    InputNodeView(const detail::InputDocumentData*, core::Index) noexcept;
    const detail::InputDocumentData* owner_;
    core::Index node_;
};
/// Borrowed direct array elements in input order, including container-valued elements.
/// The sequence is not a preorder interval of node ordinals.
class InputArrayView {
public:
    /// Number of direct elements; descendants are not counted again.
    [[nodiscard]] core::Index size() const noexcept;
    /// Return a direct child, or InvalidInput for an out-of-range position.
    [[nodiscard]] core::Result<InputNodeView> at(core::Index) const& noexcept;
    [[nodiscard]] core::Result<InputNodeView> at(core::Index) const&& = delete;

private:
    friend struct detail::InputDocumentFactory;
    InputArrayView(const detail::InputDocumentData*, core::Index) noexcept;
    const detail::InputDocumentData* owner_;
    core::Index node_;
};
/// Borrowed decoded key and value; both belong to the same document.
struct InputMemberView {
    std::string_view key; ///< Decoded key bytes, including escaped NUL when present.
    InputNodeView value;  ///< This member's direct value node.
};
/// Borrowed unique decoded keys and their values, retaining original member order.
class InputObjectView {
public:
    /// Number of members of this object, excluding descendant members.
    [[nodiscard]] core::Index size() const noexcept;
    /// Return an ordered member, or InvalidInput for an out-of-range position.
    [[nodiscard]] core::Result<InputMemberView> at(core::Index) const& noexcept;
    [[nodiscard]] core::Result<InputMemberView> at(core::Index) const&& = delete;
    /// Compare decoded UTF-8 bytes exactly. Missing is distinct from a null node.
    [[nodiscard]] std::optional<InputNodeView> find(std::string_view) const& noexcept;
    [[nodiscard]] std::optional<InputNodeView> find(std::string_view) const&& = delete;

private:
    friend struct detail::InputDocumentFactory;
    InputObjectView(const detail::InputDocumentData*, core::Index) noexcept;
    const detail::InputDocumentData* owner_;
    core::Index node_;
};
/// One immutable borrow of an owned source, tree root and admitted parser model.
/// Envelope success qualifies neither nested schema nor a physical model.
struct InputDocumentView {
    SourceIdentity identity;                 ///< Owned copy of source ID, length and digest.
    std::span<const std::byte> source_bytes; ///< Borrow of exact original bytes.
    InputNodeView root;                    ///< Root object; nested physical schemas are unchecked.
    ParserAllocationPolicy policy;         ///< Owned copy of the admitted request.
    const ParserMemoryPlan& parser_memory; ///< Borrow of the produced external forecast.
};
/// Immutable owned source/tree/root-envelope result, independent of parser,
/// caller bytes and ResourcePlan lifetimes. Implicit copying is disabled.
class InputDocument {
public:
    /// Transfer ownership; the source becomes absent and rejects view().
    InputDocument(InputDocument&&) noexcept;
    /// Release old storage and transfer ownership, invalidating its previous views.
    InputDocument& operator=(InputDocument&&) noexcept;
    /// Copying is forbidden so ownership cannot bypass an explicit preparation budget.
    InputDocument(const InputDocument&) = delete;
    /// Copy assignment is forbidden for the same ownership reason.
    InputDocument& operator=(const InputDocument&) = delete;
    /// Release owned buffers without refunding successful preparation charges.
    ~InputDocument();
    /// Borrow this document; a moved-from owner returns InvalidInput.
    [[nodiscard]] core::Result<InputDocumentView> view() const& noexcept;
    [[nodiscard]] core::Result<InputDocumentView> view() const&& = delete;

private:
    friend struct detail::InputDocumentFactory;
    explicit InputDocument(core::PlannedBuffer<detail::InputDocumentData>) noexcept;
    core::PlannedBuffer<detail::InputDocumentData> data_;
};

/// Own exact source bytes and all JSON kinds, and check only the root envelope.
/// Input is borrowed synchronously and must remain immutable during the call.
/// Do not overlap input bytes with limits, diagnostic or ResourcePlan storage;
/// reentrant/concurrent sharing of a mutable plan or diagnostic is prohibited.
/// Diagnostic is reset on entry and empty on success; failures own all available
/// detail. Use a distinct diagnostic slot and ResourcePlan in concurrent calls.
/// Failures preserve input/old documents and never refund earlier successful
/// ResourcePlan allocations. Unknown memory models/unproved strict requests are
/// Unsupported; quotas/overflow/allocation are ResourceExceeded; parser numerical
/// overflow is Nonfinite; malformed text/envelope is InvalidInput.
[[nodiscard]] core::Result<InputDocument>
read_input_document(std::span<const std::byte> bytes, SourceId source, const InputLimits& limits,
                    core::ResourcePlan& resources, std::optional<InputDiagnostic>& diagnostic,
                    core::FailureDetails details = {}) noexcept;
} // namespace uvlas::io
