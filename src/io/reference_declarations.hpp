#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <stop_token>
#include <string_view>
#include <uvlas/io/case_schema.hpp>

namespace uvlas::io::detail {
struct ReferenceDeclarationData;
struct ReferenceDeclarationFactory;

/// Purpose declared by a validated FileRef; no file parser or byte validation is implied.
enum class ReferenceRole : std::uint8_t {
    NastranDeck,
    NastranSidecar,
    GmshMesh,
    GmshSidecar,
    ReferenceGeometry,
    PreparationHistory,
    InitialHistory,
    ObservationReference
};
/// Return the exact schema spelling, or no value for an invalid enum representation.
[[nodiscard]] std::optional<std::string_view> reference_role_name(ReferenceRole) noexcept;
/// Independent explicit limits for the internal declaration-extraction call.
enum class DeclarationQuota : std::uint8_t {
    SourceNodes,     ///< All nodes in the source case, including unreferenced values.
    SourceDepth,     ///< Maximum source container depth; the root object is depth one.
    Declarations,    ///< Every FileRef occurrence, including repeated source IDs.
    SourcePathBytes, ///< Sum of decoded original paths, without normalization.
    OwnedTextBytes,  ///< Copied IDs, paths, digest spellings and JSON Pointers.
    OwnedBytes,      ///< Returned immutable payload capacity.
    ScratchBytes,    ///< Simultaneous automatic value-payload envelope, not process stack/RSS.
    CumulativeBytes, ///< Successful payload allocations, without refunds.
    Work,            ///< Visits, compared key bytes, ancestry digits and copied bytes.
    DiagnosticBytes, ///< Complete failure JSON Pointer capacity.
    Count            ///< Sentinel; not a usable quota.
};
class DeclarationLimits {
public:
    /// Supply every bound. Zero permits none and never means unlimited.
    explicit DeclarationLimits(
        std::array<std::uint64_t, static_cast<std::size_t>(DeclarationQuota::Count)> caps) noexcept
        : caps_(caps) {}
    /// Reject Count and unknown enum values rather than indexing outside the array.
    [[nodiscard]] core::Result<std::uint64_t> get(DeclarationQuota) const noexcept;

private:
    std::array<std::uint64_t, static_cast<std::size_t>(DeclarationQuota::Count)> caps_;
};
enum class DeclarationReason : std::uint8_t {
    InvalidOwner,
    UnsupportedRole,
    QuotaExceeded,
    Overflow,
    Allocation,
    Cancelled,
    UnexpectedException
};
/// Fixed owned context. Node ordinals always refer to the identified original input document.
struct DeclarationDiagnosticDetail {
    DeclarationReason reason;
    std::optional<SchemaIdentity> parent;
    std::optional<std::uint64_t> node;
    std::optional<DeclarationQuota> quota;
    std::optional<std::uint64_t> required, limit;
    std::uint64_t work;
    SchemaPathState path_state;
};
class DeclarationDiagnostic {
public:
    DeclarationDiagnostic() noexcept = default;
    DeclarationDiagnostic(DeclarationDiagnostic&&) noexcept;
    DeclarationDiagnostic& operator=(DeclarationDiagnostic&&) noexcept;
    DeclarationDiagnostic(const DeclarationDiagnostic&) = delete;
    DeclarationDiagnostic& operator=(const DeclarationDiagnostic&) = delete;
    /// Return owned fixed context, or no value after a successful call.
    [[nodiscard]] std::optional<DeclarationDiagnosticDetail> detail() const noexcept {
        return detail_;
    }
    /// Borrow a complete original-document JSON Pointer; no truncated path is returned.
    [[nodiscard]] std::optional<std::string_view> path() const& noexcept;
    [[nodiscard]] std::optional<std::string_view> path() const&& = delete;
    /// Clear this slot without refunding earlier successful plan charges.
    void reset() noexcept;

private:
    friend struct ReferenceDeclarationFactory;
    std::optional<DeclarationDiagnosticDetail> detail_;
    core::PlannedBuffer<char> path_;
};
/// One occurrence in schema field/array order. Views borrow DeclaredReferences only.
struct DeclaredReferenceView {
    SourceId source_id;           ///< Canonical external source ID; zero is valid.
    ReferenceRole role;           ///< The role required by the referring schema field.
    SourceDigest declared_digest; ///< Declared digest bytes, not a hash computed from a file.
    std::string_view source_id_spelling;  ///< Original decoded decimal ID spelling.
    std::string_view original_path;       ///< Original decoded UTF-8 path, byte-for-byte unchanged.
    std::string_view digest_spelling;     ///< Original lowercase hexadecimal digest.
    std::uint64_t declaration_node;       ///< Original FileRef object node ordinal.
    std::string_view declaration_pointer; ///< Original-document JSON Pointer to that object.
};
struct DeclarationMemoryReport {
    std::uint64_t owned_bytes, scratch_bytes, cumulative_bytes, work;
    std::uint64_t planned_owned_bytes, planned_scratch_bytes, planned_work, admission_bytes;
};
/// Internal metadata only. This is not a resolved, normalized or acquired file snapshot.
class DeclaredReferences {
public:
    DeclaredReferences(DeclaredReferences&&) noexcept;
    DeclaredReferences& operator=(DeclaredReferences&&) noexcept;
    DeclaredReferences(const DeclaredReferences&) = delete;
    DeclaredReferences& operator=(const DeclaredReferences&) = delete;
    ~DeclaredReferences();
    [[nodiscard]] core::Result<SchemaIdentity> identity() const noexcept;
    [[nodiscard]] core::Result<std::uint64_t> size() const noexcept;
    [[nodiscard]] core::Result<DeclaredReferenceView> entry(core::Index) const& noexcept;
    [[nodiscard]] core::Result<DeclaredReferenceView> entry(core::Index) const&& = delete;
    [[nodiscard]] core::Result<DeclarationMemoryReport> memory_report() const noexcept;

private:
    friend struct ReferenceDeclarationFactory;
    explicit DeclaredReferences(core::PlannedBuffer<ReferenceDeclarationData>) noexcept;
    core::PlannedBuffer<ReferenceDeclarationData> data_;
};
/// Copy every validated FileRef in dictionary-field and array order, retaining all occurrences.
/// No SI tolerances, path normalization, file I/O, File/Memory origin or role parser is used.
/// Independent plans/slots may run concurrently. Never share or reenter a mutable plan/slot.
[[nodiscard]] core::Result<DeclaredReferences>
collect_case_file_references(const SchemaCase&, const DeclarationLimits&, core::ResourcePlan&,
                             std::stop_token, DeclarationDiagnostic*) noexcept;
} // namespace uvlas::io::detail
