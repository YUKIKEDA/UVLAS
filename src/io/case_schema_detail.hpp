#pragma once

#include "input_document_detail.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stop_token>
#include <string_view>
#include <uvlas/io/case_schema.hpp>

namespace uvlas::io::detail {

struct DecimalRange {
    int sign = 0;
    std::uint64_t offset = 0;
    std::uint64_t count = 1;
    std::int64_t exponent = 0;
};
struct SchemaNode {
    FieldType type = FieldType::T0;
    std::uint16_t record_kind = 0;
    std::uint64_t record_row = 0;
    std::optional<std::uint64_t> parent;
    std::optional<std::uint64_t> parent_member;
    std::optional<std::uint64_t> parent_array_index;
    std::optional<std::uint64_t> reference_node;
    std::optional<SchemaNamespace> reference_domain;
    std::uint64_t array_offset = 0;
    std::uint64_t unsigned_value = 0;
    std::int64_t signed_value = 0;
    DecimalRange number;
    DecimalRange rational_numerator;
    DecimalRange rational_denominator;
    SchemaDimension dimensions{{0, 0, 0, 0}};
    std::uint64_t components = 1;
    bool assigned = false;
};
struct FieldDefinition {
    std::string_view name;
    FieldType type;
    bool optional;
};
struct RecordDefinition {
    std::string_view name;
    std::string_view family;
    std::span<const FieldDefinition> fields;
    std::string_view selector;
    std::string_view alternatives;
};
struct TypeDefinition {
    std::string_view spelling;
    FieldKind kind;
    std::optional<FieldType> element;
    std::uint64_t extent;
    std::string_view family;
    SchemaDimension dimensions;
    bool dynamic_dimension;
    std::uint64_t components;
};
struct SchemaIdentifier {
    std::uint8_t domain = 0;
    std::uint64_t value = 0;
    std::uint64_t node = 0;
};
struct SchemaReference {
    std::uint8_t domain = 0;
    std::uint64_t value = 0;
    std::uint64_t node = 0;
};
struct ScanFrame {
    std::uint64_t node = 0;
    FieldType type = FieldType::T0;
    std::uint64_t depth = 0;
    SchemaDimension dimensions{{0, 0, 0, 0}};
    std::uint64_t components = 1;
    bool arithmetic = false;
    std::uint64_t path_bytes = 0;
};
struct SchemaBootstrap {
    ScanFrame frame;
    std::optional<std::uint64_t> previous_sibling;
    std::optional<std::uint64_t> last_child;
};
struct ExactScratch {
    core::PlannedBuffer<char> digits;
    std::uint64_t stride = 0;
};
struct SchemaData;
struct SchemaTokenLimit {
    std::uint64_t node;
    std::uint64_t bytes;
    FieldType type;
};
struct SchemaContext {
    const InputDocumentData* input;
    std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> limits;
    std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> used{};
    core::ResourcePlan& resources;
    SchemaDiagnostic* diagnostic;
    SchemaData* output = nullptr;
    std::optional<std::uint64_t> current;
    std::optional<std::uint16_t> expected_field;
    std::optional<FieldType> expected_type;
    std::optional<std::uint16_t> expected_variant;
    std::optional<SourceIdentity> source;
    std::uint64_t scratch_bytes = 0;
    std::uint64_t owned_bytes = 0;
    std::uint64_t initial_charged = 0;
    ExactScratch exact;
    std::uint64_t path_required = 0;
    std::uint64_t planned_identifiers = 0;
    std::uint64_t planned_references = 0;
    std::uint64_t rational_digits = 1;
    std::optional<std::int64_t> arithmetic_low;
    std::optional<std::int64_t> arithmetic_high;
    bool planning_overflow = false;
    std::optional<SchemaTokenLimit> forecast_token_excess;
    std::uint64_t planned_owned = 0;
    std::uint64_t planned_scratch = 0;
    std::uint64_t planned_cumulative = 0;
    std::uint64_t payload_root = 0;
    FieldType payload_type = FieldType::T0;
    std::uint64_t payload_components = 1;
    SchemaDimension payload_dimensions{{0, 0, 0, 0}};
    bool whole_case = true;
    bool case_registry = false;
    std::stop_token cancellation;
    SchemaContext(const InputDocumentData*, SchemaLimits, core::ResourcePlan&, SchemaDiagnostic*);
    void work(std::uint64_t amount = 1);
    void charge(SchemaQuota, std::uint64_t amount);
    [[noreturn]] void reject(SchemaReason, core::FailureCode = core::FailureCode::InvalidInput,
                             std::optional<SchemaQuota> = {}, std::optional<std::uint64_t> = {});
    std::uint64_t add(std::uint64_t, std::uint64_t);
    std::uint64_t multiply(std::uint64_t, std::uint64_t);
    template <class T>
    core::PlannedBuffer<T> allocate(std::uint64_t count, bool scratch) {
        const auto bytes = multiply(count, sizeof(T));
        if (scratch) {
            scratch_bytes = add(scratch_bytes, bytes);
            if (scratch_bytes > limits[static_cast<std::size_t>(SchemaQuota::ScratchBytes)])
                reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                       SchemaQuota::ScratchBytes, scratch_bytes);
        } else {
            owned_bytes = add(owned_bytes, bytes);
            if (owned_bytes > limits[static_cast<std::size_t>(SchemaQuota::TypedBytes)])
                reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                       SchemaQuota::TypedBytes, owned_bytes);
        }
        if (count > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
            reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
        work(bytes);
        auto value =
            resources.allocate<T>(static_cast<core::Index>(count), core::Operation::PrepareModel);
        if (!value)
            throw value.error();
        return std::move(*value);
    }
};

std::span<const TypeDefinition> schema_types() noexcept;
std::span<const RecordDefinition> schema_records() noexcept;
FieldType schema_root_type() noexcept;
std::uint64_t schema_record_bytes(SchemaContext&);
void schema_prepare_records(SchemaContext&);
void schema_build_records(SchemaContext&);
std::optional<std::uint64_t> schema_find(SchemaContext&, std::uint64_t, std::string_view);
InputNodeKind node_kind(const InputDocumentData&, std::uint64_t) noexcept;
std::string_view node_text(const InputDocumentData&, std::uint64_t) noexcept;
std::uint64_t array_size(const InputDocumentData&, std::uint64_t) noexcept;
DecimalRange schema_decimal(SchemaContext&, std::uint64_t, std::string_view,
                            bool numerator = false);
int schema_compare(SchemaContext&, const DecimalRange&, const DecimalRange&);
int schema_compare_integer(SchemaContext&, const DecimalRange&, std::int64_t);
void schema_exact_semantics(SchemaContext&);
void schema_semantics(SchemaContext&);

} // namespace uvlas::io::detail
