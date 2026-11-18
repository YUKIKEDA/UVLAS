#include "si_case_detail.hpp"

#include <algorithm>

namespace uvlas::io::detail {
bool same_identity(const SchemaIdentity& a, const SchemaIdentity& b) noexcept {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}
std::uint64_t si_snapshot_bytes(SiContext& c, const InputDocumentData& p) {
    auto bytes = c.add(p.source.size(), p.text.size());
    bytes = c.add(bytes, c.multiply(p.nodes.size(), sizeof(InputNode)));
    bytes = c.add(bytes, c.multiply(p.edges.size(), sizeof(ArrayEdge)));
    bytes = c.add(bytes, c.multiply(p.members.size(), sizeof(ObjectMember)));
    return c.add(bytes, c.multiply(p.keys.size(), sizeof(KeySlot)));
}
template <class T>
void copy_buffer(SiContext& c, const core::PlannedBuffer<T>& from, core::PlannedBuffer<T>& to) {
    to = c.allocate<T>(from.size());
    c.work.step(c.multiply(from.size(), sizeof(T)));
    std::ranges::copy(from.values(), to.values().begin());
}
void si_clone_snapshot(SiContext& c, const InputDocumentData& from, InputDocumentData& to) {
    to.identity = from.identity;
    to.parser_plan = from.parser_plan;
    to.policy = from.policy;
    to.node_count = from.node_count;
    to.edge_count = from.edge_count;
    to.member_count = from.member_count;
    to.text_count = from.text_count;
    copy_buffer(c, from.source, to.source);
    copy_buffer(c, from.nodes, to.nodes);
    copy_buffer(c, from.edges, to.edges);
    copy_buffer(c, from.members, to.members);
    copy_buffer(c, from.text, to.text);
    copy_buffer(c, from.keys, to.keys);
}
namespace {
SchemaLimits clone_limits(SiContext& c) {
    std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> values{};
    // Only storage construction and typed rebind are called, never schema validation.
    values[static_cast<std::size_t>(SchemaQuota::TypedBytes)] =
        c.limit(SiQuota::OwnedBytes) - c.owned;
    values[static_cast<std::size_t>(SchemaQuota::ScratchBytes)] = c.limit(SiQuota::ScratchBytes);
    values[static_cast<std::size_t>(SchemaQuota::Work)] = c.work.limit - c.work.used;
    return SchemaLimits{values};
}
} // namespace
std::uint64_t si_schema_clone_bytes(SiContext& c, const SchemaData& source) {
    auto bytes = c.add(sizeof(SchemaData), sizeof(InputDocumentData));
    bytes = c.add(bytes, si_snapshot_bytes(c, source.snapshot.values()[0]));
    bytes = c.add(bytes, c.multiply(source.nodes.size(), sizeof(SchemaNode)));
    bytes = c.add(bytes, c.multiply(source.array_elements.size(), sizeof(std::uint64_t)));
    bytes = c.add(bytes, source.digits.size());
    bytes = c.add(bytes, c.multiply(source.identifiers.size(), sizeof(SchemaIdentifier)));
    bytes = c.add(bytes, c.multiply(source.references.size(), sizeof(SchemaReference)));
    SchemaData counter;
    counter.record_counts = source.record_counts;
    SchemaContext context{&source.snapshot.values()[0], clone_limits(c), c.resources, nullptr};
    context.output = &counter;
    context.cancellation = c.limits.stop();
    bytes = c.add(bytes, schema_record_bytes(context));
    c.work.step(context.used[static_cast<std::size_t>(SchemaQuota::Work)]);
    return bytes;
}
SchemaCase si_clone_schema(SiContext& c, const SchemaData& source) {
    auto data = c.allocate<SchemaData>(1);
    auto& target = data.values()[0];
    target.snapshot = c.allocate<InputDocumentData>(1);
    si_clone_snapshot(c, source.snapshot.values()[0], target.snapshot.values()[0]);
    copy_buffer(c, source.nodes, target.nodes);
    copy_buffer(c, source.array_elements, target.array_elements);
    copy_buffer(c, source.digits, target.digits);
    copy_buffer(c, source.identifiers, target.identifiers);
    copy_buffer(c, source.references, target.references);
    target.record_counts = source.record_counts;
    target.digit_count = source.digit_count;
    target.identifier_count = source.identifier_count;
    target.reference_count = source.reference_count;
    target.report = source.report;
    SchemaContext context{&target.snapshot.values()[0], clone_limits(c), c.resources, nullptr};
    context.output = &target;
    context.cancellation = c.limits.stop();
    // The producer already proved every shape and reference. Rebuild only owner-bound fields.
    schema_prepare_records(context);
    schema_build_records(context);
    c.owned = c.add(c.owned, context.owned_bytes);
    c.require(SiQuota::OwnedBytes, c.owned);
    c.work.step(context.used[static_cast<std::size_t>(SchemaQuota::Work)]);
    return SchemaFactory::case_result(std::move(data));
}
SiUnits si_units(const SchemaData& data) {
    const auto& root = ready_record(data.records_Case.values()[0]);
    const auto units = root.units.record();
    return {static_cast<SiLengthUnit>(units.length.selected()),
            static_cast<SiMassUnit>(units.mass.selected()),
            static_cast<SiTimeUnit>(units.time.selected()),
            static_cast<SiAngleUnit>(units.angle.selected())};
}
bool si_numeric_node(const SchemaNode& value) noexcept {
    if (!value.assigned)
        return false;
    const auto kind = field_kinds[static_cast<std::size_t>(value.type)];
    return kind == FieldKind::Number || kind == FieldKind::ExactQuantity ||
           value.type == FieldType::T62;
}
SiEntryStorage si_entry(SiContext& c, std::span<const SchemaNode> nodes, std::uint64_t node) {
    std::uint64_t child = node;
    while (nodes[static_cast<std::size_t>(child)].parent) {
        c.work.step();
        const auto& info = nodes[static_cast<std::size_t>(child)];
        const auto parent = *info.parent;
        const auto& owner = nodes[static_cast<std::size_t>(parent)];
        if (info.parent_member) {
            const auto& source =
                c.snapshot->members.values()[static_cast<std::size_t>(*info.parent_member)];
            const std::string_view key{c.snapshot->text.values().data() + source.key.offset,
                                       source.key.size};
            const auto& record = schema_records()[owner.record_kind];
            for (std::size_t f = 0; f < record.fields.size(); ++f) {
                c.work.step();
                if (record.fields[f].name == key) {
                    if (f > std::numeric_limits<std::uint16_t>::max())
                        c.reject(SiReason::InvalidOwner);
                    return {node, parent, owner.record_kind, static_cast<std::uint16_t>(f)};
                }
            }
            c.reject(SiReason::InvalidOwner);
        }
        child = parent;
    }
    c.reject(SiReason::InvalidOwner);
}
SiRecipeView si_recipe(const InputDocumentData& snapshot, std::span<const SchemaNode> nodes,
                       std::span<const char> digits, std::uint64_t node, SiUnits units) {
    const auto& n = nodes[static_cast<std::size_t>(node)];
    const auto view = [digits](const DecimalRange& value) {
        return SiDecimalView{value.sign,
                             {digits.data() + value.offset, static_cast<std::size_t>(value.count)},
                             value.exponent};
    };
    const auto kind = node_kind(snapshot, node);
    const auto value = [&]() -> std::variant<SiDecimalView, SiRationalView> {
        if (kind == InputNodeKind::Object)
            return SiRationalView{view(n.rational_numerator), view(n.rational_denominator)};
        return view(n.number);
    }();
    return {value, node_text(snapshot, node), units, n.dimensions,
            units.angle == SiAngleUnit::Degree ? n.dimensions.exponents[3] : std::int16_t{0}};
}
SiCatalogEntry si_catalog_entry(const SiCatalogData& owner, std::uint64_t row) {
    const auto& entry = owner.entries.values()[static_cast<std::size_t>(row)];
    const auto& node = owner.nodes.values()[static_cast<std::size_t>(entry.node)];
    const auto& record = schema_records()[entry.record];
    return {SiFactory::key(owner.key_owner, *owner.identity, entry.node,
                           static_cast<std::uint16_t>(node.type), entry.record),
            si_recipe(owner.snapshot.values()[0], owner.nodes.values(), owner.digits.values(),
                      entry.node, owner.units),
            record.name, record.fields[entry.field].name};
}
SiRequirementView si_requirement_view(const SiRequirementStorage& row, std::span<const char> text) {
    const auto read = [text](SiTextRange range) {
        return std::string_view{text.data() + range.offset, static_cast<std::size_t>(range.count)};
    };
    const auto input = [&](const SiExactStorage& value) -> SiExactInput {
        if (const auto* decimal = std::get_if<SiExactStorage::Decimal>(&value.value))
            return SiDecimalInput{read(decimal->token)};
        const auto& rational = std::get<SiExactStorage::Rational>(value.value);
        return SiRationalInput{read(rational.numerator), read(rational.denominator)};
    };
    return {{input(row.scale), row.dimensions}, input(row.absolute), input(row.relative)};
}
} // namespace uvlas::io::detail
