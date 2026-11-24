#pragma once

#include "nastran_number_fields.hpp"
#include "reference_root_snapshot.hpp"

namespace uvlas::io::detail {

enum class NastranBulkFormat : std::uint8_t { Small, Large, Free, FreeLarge };

struct NastranBulkPosition {
    std::size_t byte = 0;   // Zero based.
    std::size_t line = 1;   // One based, counting comments and blank physical lines.
    std::size_t column = 1; // One based, including implied fixed-width padding.
    bool operator==(const NastranBulkPosition&) const = default;
};

struct NastranBulkLexeme {
    NastranBulkPosition position;
    std::size_t size = 0; // Zero for a blank field. Implied padding points at the end of data.
};

struct NastranBulkCard {
    NastranBulkFormat format = NastranBulkFormat::Small; // Format of the first physical line.
    NastranBulkLexeme name;
    std::size_t first_field = 0;
    std::size_t field_count = 0;
};

enum class NastranBulkFailure : std::uint8_t {
    Snapshot,
    SourceByte,
    LineWidth,
    Header,
    Field,
    Continuation,
    UnsupportedForm,
    Limit
};

struct NastranBulkDiagnostic {
    NastranBulkFailure reason = NastranBulkFailure::Snapshot;
    NastranBulkPosition position;
};

struct NastranBulkLimits {
    std::uint64_t file_bytes;
    std::uint64_t physical_line_bytes;
    std::uint64_t physical_lines;
    std::uint64_t cards;
    std::uint64_t fields;
    std::uint64_t work;
};

struct NastranBulkView {
    RootSnapshotLocation location;
    SourceDigest digest;
    std::string_view bytes;
    std::span<const NastranBulkCard> cards;
    std::span<const NastranBulkLexeme> fields;
    std::uint64_t work;

    [[nodiscard]] std::string_view text(NastranBulkLexeme) const noexcept;
    // The card dictionary selects the field's type. The spelling and source position
    // remain available independently of conversion.
    [[nodiscard]] core::Result<std::int64_t> integer(NastranBulkLexeme field) const noexcept {
        return read_nastran_integer_field(text(field));
    }
    [[nodiscard]] core::Result<NastranRealField> real(NastranBulkLexeme field) const noexcept {
        return read_nastran_real_field(text(field));
    }
    [[nodiscard]] std::span<const NastranBulkLexeme>
    data_fields(const NastranBulkCard&) const noexcept;
};

// Owns the exact authenticated nastran-deck snapshot and lexical card/field records.
// The input owner is moved only on success. No files are opened and no card is
// physically interpreted. Fields borrow this owner until it moves or is destroyed.
class NastranBulkDeck {
public:
    NastranBulkDeck(NastranBulkDeck&&) noexcept = default;
    NastranBulkDeck& operator=(NastranBulkDeck&&) = delete;

    [[nodiscard]] static core::Result<NastranBulkDeck>
    prepare(AcquiredRootSnapshot&&, NastranBulkLimits, core::ResourcePlan&, const std::stop_token&,
            std::optional<NastranBulkDiagnostic>* = nullptr) noexcept;
    [[nodiscard]] std::optional<NastranBulkView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranBulkView> view() const&& = delete;

private:
    NastranBulkDeck(AcquiredRootSnapshot&&, core::PlannedBuffer<NastranBulkCard>&&,
                    core::PlannedBuffer<NastranBulkLexeme>&&, std::uint64_t) noexcept;
    AcquiredRootSnapshot source_;
    core::PlannedBuffer<NastranBulkCard> cards_;
    core::PlannedBuffer<NastranBulkLexeme> fields_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
