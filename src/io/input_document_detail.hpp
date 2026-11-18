#pragma once

#include <optional>
#include <uvlas/io/input_document.hpp>
#include <variant>

namespace uvlas::io::detail {
enum class ParserFailureClass : std::uint8_t { Syntax, NumberRange, Other };
struct ParserFailureClassification {
    core::FailureCode code;
    InputDiagnosticReason reason;
};
// Used by the real SAX handler and direct fallback controls, not a test-only hook.
[[nodiscard]] ParserFailureClassification classify_parser_failure(ParserFailureClass, int,
                                                                  std::string_view) noexcept;
struct TextRange {
    std::size_t offset = 0;
    std::size_t size = 0;
};
struct StringValue {
    TextRange text;
};
struct DecimalValue {
    TextRange text;
};
struct ChildList {
    core::Index head = -1;
    core::Index tail = -1;
    core::Index size = 0;
};
struct ArrayValue {
    ChildList children;
};
struct ObjectValue {
    ChildList members;
};
using NodeValue = std::variant<std::monostate, bool, StringValue, std::int64_t, std::uint64_t,
                               DecimalValue, ArrayValue, ObjectValue>;
struct InputNode {
    NodeValue value;
    std::uint64_t event = 0;
};
struct ArrayEdge {
    core::Index child = -1;
    core::Index next = -1;
};
struct ObjectMember {
    TextRange key;
    core::Index child = -1;
    core::Index next = -1;
    std::uint64_t event = 0;
};
struct KeySlot {
    core::Index parent = -1;
    core::Index member = -1;
    std::uint64_t hash = 0;
};
struct BuildFrame {
    core::Index node = -1;
    core::Index pending_member = -1;
};
struct InputDocumentData {
    std::optional<SourceIdentity> identity;
    std::optional<ParserMemoryPlan> parser_plan;
    ParserAllocationPolicy policy{ParserAllocationMode::EstimateAdmission, 0};
    core::PlannedBuffer<char> source;
    core::PlannedBuffer<InputNode> nodes;
    core::PlannedBuffer<ArrayEdge> edges;
    core::PlannedBuffer<ObjectMember> members;
    core::PlannedBuffer<char> text;
    core::PlannedBuffer<KeySlot> keys;
    std::size_t node_count = 0, edge_count = 0, member_count = 0, text_count = 0;
};
// Only these factories can bind a node ordinal/view to actual owned storage.
struct InputDocumentFactory {
    static const InputDocumentData* storage(const InputDocument& document) noexcept {
        return document.data_.size() ? &document.data_.values()[0] : nullptr;
    }
    static const InputDocumentData* storage(const InputNodeView& node) noexcept {
        return node.owner_;
    }
    static core::Index position(const InputNodeView& node) noexcept {
        return node.node_;
    }
    static InputNodeOrdinal ordinal(core::Index value) noexcept;
    static InputNodeView node(const InputDocumentData*, core::Index) noexcept;
    static InputArrayView array(const InputDocumentData*, core::Index) noexcept;
    static InputObjectView object(const InputDocumentData*, core::Index) noexcept;
    static InputDocument document(core::PlannedBuffer<InputDocumentData>) noexcept;
    static ParserMemoryPlan plan(ParserEnvironment, ParserMemoryBreakdown, std::uint64_t,
                                 std::array<ParserSourceEvidence, 12>, std::size_t) noexcept;
};
} // namespace uvlas::io::detail
