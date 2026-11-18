#include "case_schema_payload.hpp"

#include <algorithm>
#include <charconv>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {
constexpr auto operation = core::Operation::PrepareModel;
constexpr std::size_t qi(SchemaQuota quota) noexcept {
    return static_cast<std::size_t>(quota);
}
std::string_view text_range(const InputDocumentData& data, TextRange range) noexcept {
    return {data.text.values().data() + range.offset, range.size};
}
} // namespace
InputNodeKind node_kind(const InputDocumentData& data, std::uint64_t node) noexcept {
    return static_cast<InputNodeKind>(
        data.nodes.values()[static_cast<std::size_t>(node)].value.index());
}
std::string_view node_text(const InputDocumentData& data, std::uint64_t node) noexcept {
    const auto& value = data.nodes.values()[static_cast<std::size_t>(node)].value;
    if (const auto* text = std::get_if<StringValue>(&value))
        return text_range(data, text->text);
    if (const auto* text = std::get_if<DecimalValue>(&value))
        return text_range(data, text->text);
    return {};
}
std::uint64_t array_size(const InputDocumentData& data, std::uint64_t node) noexcept {
    const auto* array =
        std::get_if<ArrayValue>(&data.nodes.values()[static_cast<std::size_t>(node)].value);
    return array ? static_cast<std::uint64_t>(array->children.size) : 0;
}
std::optional<std::uint64_t> schema_find(SchemaContext& context, std::uint64_t node,
                                         std::string_view key) {
    const auto* object = std::get_if<ObjectValue>(
        &context.input->nodes.values()[static_cast<std::size_t>(node)].value);
    if (!object)
        return std::nullopt;
    auto member = object->members.head;
    while (member >= 0) {
        context.work();
        const auto& entry = context.input->members.values()[static_cast<std::size_t>(member)];
        const auto text = text_range(*context.input, entry.key);
        bool equal = text.size() == key.size();
        if (equal)
            for (std::size_t i = 0; i < key.size(); ++i) {
                context.work();
                if (text[i] != key[i]) {
                    equal = false;
                    break;
                }
            }
        if (equal)
            return static_cast<std::uint64_t>(entry.child);
        member = entry.next;
    }
    return std::nullopt;
}
struct DiagnosticParent {
    std::uint64_t parent;
    std::optional<std::uint64_t> member;
    std::optional<std::uint64_t> array_index;
};
bool diagnostic_work(SchemaContext& c, std::uint64_t amount = 1) noexcept {
    const auto index = qi(SchemaQuota::Work);
    if (amount > c.limits[index] - c.used[index])
        return false;
    c.used[index] += amount;
    return true;
}
std::optional<DiagnosticParent> diagnostic_parent(SchemaContext& c, std::uint64_t child) noexcept {
    if (c.output && child < c.output->nodes.size()) {
        const auto& node = c.output->nodes.values()[static_cast<std::size_t>(child)];
        if (node.parent && (node.parent_member || node.parent_array_index))
            return DiagnosticParent{*node.parent, node.parent_member, node.parent_array_index};
    }
    if (!c.input)
        return std::nullopt;
    for (std::uint64_t node = 0; node < c.input->node_count; ++node) {
        if (!diagnostic_work(c))
            return std::nullopt;
        const auto& value = c.input->nodes.values()[static_cast<std::size_t>(node)].value;
        if (const auto* object = std::get_if<ObjectValue>(&value)) {
            auto member = object->members.head;
            while (member >= 0) {
                if (!diagnostic_work(c))
                    return std::nullopt;
                const auto& entry = c.input->members.values()[static_cast<std::size_t>(member)];
                if (static_cast<std::uint64_t>(entry.child) == child)
                    return DiagnosticParent{node, static_cast<std::uint64_t>(member), std::nullopt};
                member = entry.next;
            }
        }
        if (const auto* array = std::get_if<ArrayValue>(&value)) {
            auto edge = array->children.head;
            std::uint64_t index = 0;
            while (edge >= 0) {
                if (!diagnostic_work(c))
                    return std::nullopt;
                const auto& entry = c.input->edges.values()[static_cast<std::size_t>(edge)];
                if (static_cast<std::uint64_t>(entry.child) == child)
                    return DiagnosticParent{node, std::nullopt, index};
                edge = entry.next;
                ++index;
            }
        }
    }
    return std::nullopt;
}
void failure_path(SchemaContext& c, SchemaDiagnosticDetail fixed) noexcept {
    if (!c.diagnostic)
        return;
    const auto store = [&](SchemaPathState state) {
        fixed.path_state = state;
        fixed.work = c.used[qi(SchemaQuota::Work)];
        SchemaFactory::diagnostic(*c.diagnostic, fixed);
    };
    if (!c.input || !c.current || *c.current >= c.input->node_count) {
        store(SchemaPathState::NotBuilt);
        return;
    }
    std::uint64_t length = 0, current = *c.current, ancestors = 0;
    while (current != 0) {
        if (++ancestors > c.input->node_count) {
            store(SchemaPathState::NotBuilt);
            return;
        }
        const auto parent = diagnostic_parent(c, current);
        if (!parent) {
            store(SchemaPathState::NotBuilt);
            return;
        }
        std::uint64_t segment = 1;
        if (parent->member) {
            const auto& member =
                c.input->members.values()[static_cast<std::size_t>(*parent->member)];
            const auto key = text_range(*c.input, member.key);
            for (const char ch : key) {
                if (!diagnostic_work(c)) {
                    store(SchemaPathState::NotBuilt);
                    return;
                }
                const auto increment = (ch == '~' || ch == '/') ? 2U : 1U;
                if (segment > std::numeric_limits<std::uint64_t>::max() - increment) {
                    store(SchemaPathState::PathQuota);
                    return;
                }
                segment += increment;
            }
        } else if (parent->array_index) {
            std::array<char, 32> buffer{};
            const auto converted =
                std::to_chars(buffer.data(), buffer.data() + buffer.size(), *parent->array_index);
            segment += static_cast<std::uint64_t>(converted.ptr - buffer.data());
        } else {
            store(SchemaPathState::NotBuilt);
            return;
        }
        if (length > std::numeric_limits<std::uint64_t>::max() - segment) {
            store(SchemaPathState::PathQuota);
            return;
        }
        length += segment;
        current = parent->parent;
    }
    if (length > c.limits[qi(SchemaQuota::DiagnosticPathBytes)] ||
        length > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max())) {
        store(SchemaPathState::PathQuota);
        return;
    }
    if (!diagnostic_work(c, length)) {
        store(SchemaPathState::NotBuilt);
        return;
    }
    auto allocated = c.resources.allocate<char>(static_cast<core::Index>(length), operation);
    if (!allocated) {
        store(SchemaPathState::Allocation);
        return;
    }
    auto path = std::move(*allocated);
    std::size_t position = path.size();
    current = *c.current;
    ancestors = 0;
    while (current != 0) {
        if (++ancestors > c.input->node_count) {
            store(SchemaPathState::NotBuilt);
            return;
        }
        const auto parent = diagnostic_parent(c, current);
        if (!parent) {
            store(SchemaPathState::NotBuilt);
            return;
        }
        if (parent->member) {
            const auto& member =
                c.input->members.values()[static_cast<std::size_t>(*parent->member)];
            const auto key = text_range(*c.input, member.key);
            for (std::size_t j = key.size(); j > 0; --j) {
                const auto ch = key[j - 1];
                if (ch == '~' || ch == '/') {
                    path.values()[--position] = ch == '~' ? '0' : '1';
                    path.values()[--position] = '~';
                } else
                    path.values()[--position] = ch;
            }
        } else if (parent->array_index) {
            std::array<char, 32> buffer{};
            const auto converted =
                std::to_chars(buffer.data(), buffer.data() + buffer.size(), *parent->array_index);
            for (auto p = converted.ptr; p != buffer.data();)
                path.values()[--position] = *--p;
        } else {
            store(SchemaPathState::NotBuilt);
            return;
        }
        path.values()[--position] = '/';
        current = parent->parent;
    }
    fixed.path_state = SchemaPathState::Complete;
    fixed.work = c.used[qi(SchemaQuota::Work)];
    SchemaFactory::diagnostic(*c.diagnostic, fixed, std::move(path));
}
void apply_case_limits(SchemaContext& c) {
    constexpr std::array names{"maximum_source_bytes",
                               "maximum_depth",
                               "maximum_nodes",
                               "maximum_members",
                               "maximum_decoded_string_bytes",
                               "maximum_decimal_token_bytes",
                               "maximum_total_text_bytes",
                               "maximum_typed_bytes",
                               "maximum_entities",
                               "maximum_references",
                               "maximum_scalars",
                               "maximum_polynomial_terms",
                               "maximum_table_samples",
                               "maximum_expression_nodes",
                               "maximum_expression_edges",
                               "maximum_expression_depth",
                               "maximum_expression_power",
                               "maximum_expression_points",
                               "maximum_diagnostic_path_bytes"};
    const auto resource = schema_find(c, 0, "resources");
    if (!resource)
        return;
    const auto input = schema_find(c, *resource, "input");
    if (!input)
        return;
    for (std::size_t i = 0; i < names.size(); ++i) {
        const auto value = schema_find(c, *input, names[i]);
        if (!value)
            continue;
        const auto& raw = c.input->nodes.values()[static_cast<std::size_t>(*value)].value;
        std::optional<std::uint64_t> limit;
        if (const auto* integer = std::get_if<std::uint64_t>(&raw))
            limit = *integer;
        if (const auto* integer = std::get_if<std::int64_t>(&raw); integer && *integer >= 0)
            limit = static_cast<std::uint64_t>(*integer);
        if (limit)
            c.limits[i] = std::min(c.limits[i], *limit);
    }
}
SchemaContext::SchemaContext(const InputDocumentData* document, SchemaLimits request,
                             core::ResourcePlan& plan, SchemaDiagnostic* slot)
    : input(document),
      limits(request.values_),
      resources(plan),
      diagnostic(slot),
      source(document ? document->identity : std::nullopt),
      initial_charged(plan.allocated_bytes()) {}
std::uint64_t SchemaContext::add(std::uint64_t left, std::uint64_t right) {
    if (right > std::numeric_limits<std::uint64_t>::max() - left)
        reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
    return left + right;
}
std::uint64_t SchemaContext::multiply(std::uint64_t left, std::uint64_t right) {
    if (left && right > std::numeric_limits<std::uint64_t>::max() / left)
        reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
    return left * right;
}
void SchemaContext::work(std::uint64_t amount) {
    if (cancellation.stop_requested())
        throw core::Failure{core::FailureCode::Cancelled, core::Operation::PrepareModel};
    const auto index = qi(SchemaQuota::Work);
    if (amount > limits[index] - used[index])
        reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded, SchemaQuota::Work,
               add(used[index], amount));
    used[index] += amount;
}
void SchemaContext::charge(SchemaQuota quota, std::uint64_t amount) {
    const auto index = qi(quota);
    const auto needed = add(used[index], amount);
    if (needed > limits[index])
        reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded, quota, needed);
    used[index] = needed;
}
[[noreturn]] void SchemaContext::reject(SchemaReason reason, core::FailureCode code,
                                        std::optional<SchemaQuota> quota,
                                        std::optional<std::uint64_t> required) {
    if (diagnostic) {
        std::optional<FieldType> expected = expected_type;
        std::optional<InputNodeKind> actual;
        std::optional<std::uint16_t> selected = expected_variant;
        if (current && input && *current < input->node_count) {
            actual = node_kind(*input, *current);
            if (output && *current < output->nodes.size() &&
                output->nodes.values()[static_cast<std::size_t>(*current)].assigned) {
                expected = output->nodes.values()[static_cast<std::size_t>(*current)].type;
                if (*actual == InputNodeKind::Object)
                    selected =
                        output->nodes.values()[static_cast<std::size_t>(*current)].record_kind;
            }
        }
        const SchemaDiagnosticDetail fixed{
            reason,
            source,
            current,
            quota,
            required,
            quota ? std::optional{limits[qi(*quota)]} : std::nullopt,
            expected,
            expected_field,
            selected,
            actual,
            output && current && *current < output->nodes.size()
                ? output->nodes.values()[static_cast<std::size_t>(*current)].reference_domain
                : std::nullopt,
            output && current && *current < output->nodes.size() &&
                    output->nodes.values()[static_cast<std::size_t>(*current)].reference_domain
                ? std::optional{output->nodes.values()[static_cast<std::size_t>(*current)]
                                    .unsigned_value}
                : std::nullopt,
            used[qi(SchemaQuota::Work)],
            SchemaPathState::NotBuilt};
        failure_path(*this, fixed);
    }
    throw core::Failure{code, operation};
}
namespace {
bool token_member(SchemaContext& context, std::string_view alternatives,
                  std::string_view selected) {
    std::size_t start = 0;
    for (;;) {
        const auto end = alternatives.find(',', start);
        const auto part =
            alternatives.substr(start, end == std::string_view::npos ? end : end - start);
        context.work(std::min(part.size(), selected.size()) + 1);
        if (part == selected)
            return true;
        if (end == std::string_view::npos)
            return false;
        start = end + 1;
    }
}
const RecordDefinition& select_record(SchemaContext& context, FieldType type, std::uint64_t node) {
    const auto& definition = schema_types()[static_cast<std::size_t>(type)];
    const std::string_view family =
        definition.kind == FieldKind::ExactQuantity ? "Rational" : definition.family;
    for (const auto& candidate : schema_records()) {
        context.work(std::min(candidate.family.size(), family.size()) + 1);
        if (candidate.family != family)
            continue;
        if (candidate.selector.empty())
            return candidate;
        auto selector = schema_find(context, node, candidate.selector);
        if (candidate.selector == "target.kind") {
            const auto target = schema_find(context, node, "target");
            selector = target ? schema_find(context, *target, "kind") : std::nullopt;
        }
        if (!selector)
            context.reject(SchemaReason::MissingField);
        if (node_kind(*context.input, *selector) != InputNodeKind::String)
            context.reject(SchemaReason::WrongKind);
        if (token_member(context, candidate.alternatives, node_text(*context.input, *selector)))
            return candidate;
    }
    context.reject(SchemaReason::UnsupportedTag, core::FailureCode::Unsupported);
}
void set_parent(SchemaData& data, std::uint64_t child, std::uint64_t parent,
                std::optional<std::uint64_t> member, std::optional<std::uint64_t> index) {
    auto& target = data.nodes.values()[static_cast<std::size_t>(child)];
    target.parent = parent;
    target.parent_member = member;
    target.parent_array_index = index;
}
SchemaDimension field_dimensions(SchemaContext& c, std::uint64_t node,
                                 const RecordDefinition& record, std::string_view field,
                                 SchemaDimension inherited) {
    const auto tag = [&](std::string_view key) -> std::string_view {
        const auto value = schema_find(c, node, key);
        if (!value || node_kind(*c.input, *value) != InputNodeKind::String)
            return {};
        return node_text(*c.input, *value);
    };
    const auto family = record.family;
    constexpr SchemaDimension none{{0, 0, 0, 0}}, length{{1, 0, 0, 0}}, velocity{{1, 0, -1, 0}},
        angle{{0, 0, 0, 1}}, force{{1, 1, -2, 0}}, moment{{2, 1, -2, 0}};
    if (family == "Frame" && field == "rotation")
        return none;
    if (family == "FrameStatePayload" && field == "rotation")
        return none;
    if (family == "Motion" && field == "translation")
        return length;
    if (family == "Motion" && field == "rotation_vector")
        return angle;
    if (family == "Flow" && field == "uniform")
        return velocity;
    if (family == "Gust" && field == "amplitude")
        return velocity;
    if (family == "Convection" && field == "velocity")
        return velocity;
    if (family == "Load") {
        const auto kind = tag("kind");
        if (field == "series")
            return kind == "couple" ? moment : force;
        if (field == "acceleration")
            return {{1, 0, -2, 0}};
        if (field == "multiplier")
            return none;
        if (field == "expression")
            return {{-1, 1, -2, 0}};
    }
    if (family == "ShellPhysicalState")
        return field == "displacement" ? length : (field == "velocity" ? velocity : inherited);
    if (family == "Perturbation")
        return field == "displacements" ? length : (field == "velocities" ? velocity : none);
    if (family == "CoefficientState")
        return field == "position" ? length : (field == "velocity" ? velocity : inherited);
    if (family == "ConstraintScalePolicy" && field == "position")
        return length;
    if (family == "ConstraintScalePolicy" && field == "orientation")
        return angle;
    if (family == "ReferenceField" && (field == "origin" || field == "lengths"))
        return tag("coordinates") == "chart" ? none : length;
    if (family == "ExpressionDomain") {
        const auto kind = tag("kind");
        if (kind == "reference" && (field == "lower" || field == "upper"))
            return tag("coordinates") == "chart" ? none : length;
    }
    if (family == "Observation" && (field == "budget" || field == "comparison")) {
        const auto kind = tag("kind");
        if (kind == "reaction")
            return force;
        if (kind == "stress")
            return tag("measure") == "green" ? none : SchemaDimension{{-1, 1, -2, 0}};
        const auto quantity = tag("quantity");
        if (quantity == "position" || quantity == "displacement")
            return length;
        if (quantity == "velocity")
            return velocity;
        if (quantity == "acceleration")
            return {{1, 0, -2, 0}};
        if (quantity == "force")
            return force;
        if (quantity == "moment")
            return moment;
        if (quantity == "frequency" || quantity == "growth")
            return {{0, 0, -1, 0}};
        if (quantity == "subspace")
            return angle;
        return none;
    }
    if ((family == "ObservationBudget" || family == "ErrorGroup") && field == "relative")
        return none;
    if (family == "Scale" && field == "values") {
        const auto name = tag("field");
        if (name == "shell_position" || name == "beam_translation" || name == "rigid_translation" ||
            name == "wake_position")
            return length;
        if (name == "beam_rotation" || name == "rigid_rotation")
            return angle;
        if (name == "structural_velocity" || name == "aero_residual" ||
            name == "transport_residual")
            return velocity;
        if (name == "structural_angular_velocity")
            return {{0, 0, -1, 1}};
        if (name == "circulation" || name == "wake_strength")
            return {{2, 0, -1, 0}};
        if (name == "force_balance")
            return force;
        if (name == "moment_balance")
            return moment;
    }
    return inherited;
}
SchemaDimension element_dimensions(SchemaContext& c, std::uint64_t array, std::uint64_t index,
                                   SchemaDimension inherited) {
    const auto& parent = c.output->nodes.values()[static_cast<std::size_t>(array)];
    if (parent.parent && parent.parent_member) {
        const auto& object = c.output->nodes.values()[static_cast<std::size_t>(*parent.parent)];
        if (node_kind(*c.input, *parent.parent) == InputNodeKind::Object) {
            const auto family = schema_records()[object.record_kind].family;
            const auto& member =
                c.input->members.values()[static_cast<std::size_t>(*parent.parent_member)];
            const auto key = text_range(*c.input, member.key);
            if (family == "Beam" && key == "initial_strain")
                return index < 3 ? SchemaDimension{{0, 0, 0, 0}} : SchemaDimension{{-1, 0, 0, 0}};
            if ((family == "MaterialPoint" || family == "InitialPoint") && key == "coordinates") {
                const auto target = schema_find(c, *parent.parent, "target");
                if (target) {
                    const auto kind = schema_find(c, *target, "kind");
                    if (kind) {
                        const auto value = node_text(*c.input, *kind);
                        if (value == "shell")
                            return {{0, 0, 0, 0}};
                        if (value == "rigid" || index == 0)
                            return {{1, 0, 0, 0}};
                        return {{0, 0, 0, 0}};
                    }
                }
            }
            if ((family == "ObservationBudget" || family == "ErrorGroup") &&
                (key == "absolute" || key == "scale") && index >= 3) {
                auto ancestor = parent.parent;
                while (ancestor) {
                    c.work();
                    const auto& a = c.output->nodes.values()[static_cast<std::size_t>(*ancestor)];
                    if (node_kind(*c.input, *ancestor) == InputNodeKind::Object &&
                        schema_records()[a.record_kind].family == "Observation") {
                        const auto kind = schema_find(c, *ancestor, "kind");
                        if (kind && node_text(*c.input, *kind) == "reaction")
                            return {{2, 1, -2, 0}};
                        break;
                    }
                    ancestor = a.parent;
                }
            }
        }
    }
    if (parent.parent && parent.parent_array_index) {
        const auto matrix = *parent.parent;
        const auto& matrix_info = c.output->nodes.values()[static_cast<std::size_t>(matrix)];
        if (matrix_info.parent && matrix_info.parent_member &&
            node_kind(*c.input, *matrix_info.parent) == InputNodeKind::Object) {
            const auto family =
                schema_records()[c.output->nodes
                                     .values()[static_cast<std::size_t>(*matrix_info.parent)]
                                     .record_kind]
                    .family;
            const auto& member =
                c.input->members.values()[static_cast<std::size_t>(*matrix_info.parent_member)];
            const auto key = text_range(*c.input, member.key);
            const auto row = *parent.parent_array_index;
            if ((family == "Section" && key == "stiffness") ||
                (family == "Beam" && key == "viscosity"))
                return {{static_cast<std::int16_t>(1 + (row >= 3 ? 1 : 0) + (index >= 3 ? 1 : 0)),
                         1, static_cast<std::int16_t>(key == "viscosity" ? -1 : -2), 0}};
            if (family == "Support" && (key == "stiffness" || key == "damping"))
                return {{static_cast<std::int16_t>(2 - (row < 3 ? 1 : 0) - (index < 3 ? 1 : 0)), 1,
                         static_cast<std::int16_t>(key == "damping" ? -1 : -2),
                         static_cast<std::int16_t>(-(row >= 3 ? 1 : 0) - (index >= 3 ? 1 : 0))}};
        }
    }
    return inherited;
}
void forecast_decimal(SchemaContext& c, std::string_view token) {
    std::uint64_t coefficient = 0, fraction = 0;
    bool point = false;
    std::size_t at = 0;
    bool nonzero = false;
    if (!token.empty() && token[0] == '-')
        ++at;
    while (at < token.size() && token[at] != 'e' && token[at] != 'E') {
        const auto ch = token[at++];
        c.work();
        if (ch == '.') {
            point = true;
            continue;
        }
        if (ch < '0' || ch > '9')
            return;
        ++coefficient;
        if (point)
            ++fraction;
        nonzero = nonzero || ch != '0';
    }
    std::int64_t exponent = 0;
    if (at < token.size()) {
        ++at;
        bool negative = false;
        if (at < token.size() && (token[at] == '+' || token[at] == '-'))
            negative = token[at++] == '-';
        std::uint64_t magnitude = 0;
        while (at < token.size()) {
            c.work();
            const auto ch = token[at++];
            if (ch < '0' || ch > '9')
                return;
            const auto digit = static_cast<std::uint64_t>(ch - '0');
            if (magnitude >
                (static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) - digit) /
                    10) {
                c.planning_overflow = true;
                return;
            }
            magnitude = magnitude * 10 + digit;
        }
        exponent =
            negative ? -static_cast<std::int64_t>(magnitude) : static_cast<std::int64_t>(magnitude);
    }
    if (!nonzero)
        return;
    if (fraction > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) ||
        coefficient > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
        c.planning_overflow = true;
        return;
    }
    const auto frac = static_cast<std::int64_t>(fraction);
    if (exponent < std::numeric_limits<std::int64_t>::min() + frac) {
        c.planning_overflow = true;
        return;
    }
    exponent -= frac;
    if (exponent >
        std::numeric_limits<std::int64_t>::max() - static_cast<std::int64_t>(coefficient)) {
        c.planning_overflow = true;
        return;
    }
    const auto upper = exponent + static_cast<std::int64_t>(coefficient);
    c.arithmetic_low = c.arithmetic_low ? std::min(*c.arithmetic_low, exponent) : exponent;
    c.arithmetic_high = c.arithmetic_high ? std::max(*c.arithmetic_high, upper) : upper;
}
void count_numeric_forecast(SchemaContext& c, std::uint64_t node, FieldType type) {
    if (node_kind(*c.input, node) == InputNodeKind::DecimalNumber) {
        const auto token = node_text(*c.input, node);
        if (token.size() > c.limits[qi(SchemaQuota::DecimalTokenBytes)]) {
            // Finish the complete key/kind/shape pass before reporting this limit.
            if (!c.forecast_token_excess)
                c.forecast_token_excess =
                    SchemaTokenLimit{node, static_cast<std::uint64_t>(token.size()), type};
            return;
        }
        forecast_decimal(c, token);
        return;
    }
    std::array<char, 32> buffer{};
    const auto& value = c.input->nodes.values()[static_cast<std::size_t>(node)].value;
    const auto converted = node_kind(*c.input, node) == InputNodeKind::SignedInteger
                               ? std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                               std::get<std::int64_t>(value))
                               : std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                               std::get<std::uint64_t>(value));
    forecast_decimal(c, {buffer.data(), static_cast<std::size_t>(converted.ptr - buffer.data())});
}
std::uint64_t escaped_key_size(SchemaContext& c, std::string_view key) {
    std::uint64_t bytes = 0;
    for (const char ch : key) {
        c.work();
        bytes = c.add(bytes, (ch == '~' || ch == '/') ? 2 : 1);
    }
    return bytes;
}
bool arithmetic_field(const RecordDefinition& record, std::string_view field) noexcept {
    if (record.family == "TimeInterval")
        return field == "start" || field == "end" || field == "step";
    if (record.family == "ObservationBudget")
        return field == "absolute" || field == "relative" || field == "groups";
    if (record.family == "ErrorGroup")
        return field == "absolute" || field == "relative";
    if (record.family == "Layer")
        return field == "interval";
    if (record.family == "Geometry")
        return field == "domain";
    if (record.family == "TimeSeries")
        return field == "times" || field == "outside";
    return false;
}
void shape_pass(SchemaContext& context, core::PlannedBuffer<SchemaBootstrap>& stack,
                bool planning) {
    auto& data = *context.output;
    std::uint64_t top = 0;
    if (context.case_registry) {
        // The registry is qualified in the same pass and owner as the payload.
        const auto& root = schema_records()[0];
        for (const auto& field : root.fields) {
            if (field.name == "units" || field.name == "schema")
                continue;
            const auto child = schema_find(context, 0, field.name);
            if (!child) {
                if (!field.optional)
                    context.reject(SchemaReason::MissingField);
                continue;
            }
            if (top >= stack.size())
                context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
            if (!planning)
                set_parent(data, *child, 0, std::nullopt, std::nullopt);
            stack.values()[static_cast<std::size_t>(top++)].frame = {
                *child,
                field.type,
                2,
                SchemaDimension{{0, 0, 0, 0}},
                1,
                false,
                static_cast<std::uint64_t>(field.name.size() + 1)};
        }
    }
    if (top >= stack.size())
        context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
    stack.values()[static_cast<std::size_t>(top++)].frame = {context.payload_root,
                                                             context.payload_type,
                                                             1,
                                                             context.payload_dimensions,
                                                             context.payload_components,
                                                             false,
                                                             context.case_registry ? 6ULL : 0ULL};
    while (top) {
        context.work();
        auto frame = stack.values()[static_cast<std::size_t>(--top)].frame;
        context.current = frame.node;
        context.expected_field.reset();
        context.expected_type = frame.type;
        context.expected_variant.reset();
        SchemaNode counted;
        SchemaNode* selected_node = &counted;
        if (!planning)
            selected_node = &data.nodes.values()[static_cast<std::size_t>(frame.node)];
        auto& typed = *selected_node;
        typed.type = frame.type;
        typed.assigned = true;
        typed.components = frame.components;
        typed.dimensions = frame.dimensions;
        const auto& type = schema_types()[static_cast<std::size_t>(frame.type)];
        if (!type.dynamic_dimension &&
            (type.spelling.starts_with("Q<") || type.spelling.starts_with("ReferenceField<")))
            typed.dimensions = type.dimensions;
        if (type.components)
            typed.components = type.components;
        if (type.spelling.starts_with("PolyTerm<")) {
            const auto comma = type.spelling.rfind(',');
            if (comma != std::string_view::npos && type.spelling.substr(comma + 1) == "3>")
                typed.components = 3;
        }
        const auto actual = node_kind(*context.input, frame.node);
        if (planning) {
            context.path_required = std::max(context.path_required, frame.path_bytes);
            if (type.kind == FieldKind::Identifier)
                context.planned_identifiers = context.add(context.planned_identifiers, 1);
            if (type.kind == FieldKind::Reference)
                context.planned_references = context.add(context.planned_references, 1);
            if (frame.arithmetic &&
                (actual == InputNodeKind::SignedInteger ||
                 actual == InputNodeKind::UnsignedInteger ||
                 actual == InputNodeKind::DecimalNumber) &&
                (type.kind == FieldKind::Number || type.kind == FieldKind::ExactQuantity))
                count_numeric_forecast(context, frame.node, frame.type);
        }
        if ((actual == InputNodeKind::Object || actual == InputNodeKind::Array) &&
            frame.depth > context.limits[qi(SchemaQuota::InputDepth)])
            context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                           SchemaQuota::InputDepth, frame.depth);
        switch (type.kind) {
        case FieldKind::Object:
        case FieldKind::ExactQuantity: {
            if (type.kind == FieldKind::ExactQuantity && actual != InputNodeKind::Object) {
                if (actual != InputNodeKind::SignedInteger &&
                    actual != InputNodeKind::UnsignedInteger &&
                    actual != InputNodeKind::DecimalNumber)
                    context.reject(SchemaReason::WrongKind);
                break;
            }
            if (actual != InputNodeKind::Object)
                context.reject(SchemaReason::WrongKind);
            const auto& selected = select_record(context, frame.type, frame.node);
            typed.record_kind = static_cast<std::uint16_t>(&selected - schema_records().data());
            context.expected_variant = typed.record_kind;
            typed.record_row = data.record_counts[typed.record_kind]++;
            if (planning && frame.arithmetic && selected.family == "Rational") {
                for (auto key : {"numerator", "denominator"}) {
                    const auto value = schema_find(context, frame.node, key);
                    if (value)
                        context.rational_digits = std::max(
                            context.rational_digits,
                            static_cast<std::uint64_t>(node_text(*context.input, *value).size()));
                }
            }
            const auto& object = std::get<ObjectValue>(
                context.input->nodes.values()[static_cast<std::size_t>(frame.node)].value);
            auto member = object.members.head;
            while (member >= 0) {
                context.work();
                const auto& entry =
                    context.input->members.values()[static_cast<std::size_t>(member)];
                const auto key = text_range(*context.input, entry.key);
                if (!planning)
                    set_parent(data, static_cast<std::uint64_t>(entry.child), frame.node,
                               static_cast<std::uint64_t>(member), std::nullopt);
                if (key.size() > context.limits[qi(SchemaQuota::DecodedStringBytes)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::DecodedStringBytes, key.size());
                const auto found =
                    std::ranges::find_if(selected.fields, [&](const FieldDefinition& field) {
                        context.work(std::min(key.size(), field.name.size()) + 1);
                        return field.name == key;
                    });
                if (found == selected.fields.end()) {
                    context.current = static_cast<std::uint64_t>(entry.child);
                    context.reject(SchemaReason::UnknownField);
                }
                member = entry.next;
            }
            for (const auto& field : selected.fields) {
                context.work();
                const auto child = schema_find(context, frame.node, field.name);
                if (!child) {
                    if (!field.optional) {
                        context.expected_field =
                            static_cast<std::uint16_t>(&field - selected.fields.data());
                        context.reject(SchemaReason::MissingField);
                    }
                    continue;
                }
                if (!planning)
                    set_parent(data, *child, frame.node, std::nullopt, std::nullopt);
            }
            // Reverse pushes preserve the original member order when the stack is consumed.
            for (core::Index position = object.members.size; position > 0; --position) {
                auto selected_member = object.members.head;
                for (core::Index at = 1; at < position; ++at) {
                    context.work();
                    selected_member =
                        context.input->members.values()[static_cast<std::size_t>(selected_member)]
                            .next;
                }
                const auto& entry =
                    context.input->members.values()[static_cast<std::size_t>(selected_member)];
                const auto key = text_range(*context.input, entry.key);
                const auto field =
                    std::ranges::find_if(selected.fields, [&](const FieldDefinition& candidate) {
                        context.work(std::min(key.size(), candidate.name.size()) + 1);
                        return candidate.name == key;
                    });
                if (!planning)
                    set_parent(data, static_cast<std::uint64_t>(entry.child), frame.node,
                               static_cast<std::uint64_t>(selected_member), std::nullopt);
                if (top >= stack.size())
                    context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
                stack.values()[static_cast<std::size_t>(top++)].frame = {
                    static_cast<std::uint64_t>(entry.child),
                    field->type,
                    frame.depth + 1,
                    planning ? typed.dimensions
                             : field_dimensions(context, frame.node, selected, field->name,
                                                typed.dimensions),
                    typed.components,
                    frame.arithmetic || arithmetic_field(selected, field->name),
                    context.add(frame.path_bytes, context.add(1, escaped_key_size(context, key)))};
            }
            break;
        }
        case FieldKind::Array: {
            if (actual != InputNodeKind::Array)
                context.reject(SchemaReason::WrongKind);
            const auto count = array_size(*context.input, frame.node);
            auto expected = type.extent;
            if (!expected && !type.spelling.starts_with("List<"))
                expected = typed.components;
            if (expected && count != expected)
                context.reject(SchemaReason::WrongShape);
            if (!type.element)
                context.reject(SchemaReason::WrongShape);
            auto next_child = stack.values()[static_cast<std::size_t>(frame.node)].last_child;
            for (std::uint64_t i = count; i > 0; --i) {
                context.work();
                if (!next_child)
                    context.reject(SchemaReason::InvalidOwner);
                const auto child = *next_child;
                next_child = stack.values()[static_cast<std::size_t>(child)].previous_sibling;
                if (!planning)
                    set_parent(data, child, frame.node, std::nullopt, i - 1);
                if (top >= stack.size())
                    context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
                std::array<char, 32> index_buffer{};
                stack.values()[static_cast<std::size_t>(top++)].frame = {
                    child,
                    *type.element,
                    frame.depth + 1,
                    planning ? typed.dimensions
                             : element_dimensions(context, frame.node, i - 1, typed.dimensions),
                    typed.components,
                    frame.arithmetic,
                    context.add(
                        frame.path_bytes,
                        context.add(
                            1, static_cast<std::uint64_t>(
                                   std::to_chars(index_buffer.data(),
                                                 index_buffer.data() + index_buffer.size(), i - 1)
                                       .ptr -
                                   index_buffer.data())))};
            }
            break;
        }
        case FieldKind::Number:
            if (actual != InputNodeKind::SignedInteger &&
                actual != InputNodeKind::UnsignedInteger && actual != InputNodeKind::DecimalNumber)
                context.reject(SchemaReason::WrongKind);
            break;
        case FieldKind::Integer:
            if (actual != InputNodeKind::SignedInteger && actual != InputNodeKind::UnsignedInteger)
                context.reject(SchemaReason::WrongKind);
            break;
        case FieldKind::Identifier:
        case FieldKind::Reference:
        case FieldKind::Enumeration:
        case FieldKind::Text:
            if (actual != InputNodeKind::String)
                context.reject(SchemaReason::WrongKind);
            break;
        case FieldKind::Boolean:
            if (actual != InputNodeKind::Boolean)
                context.reject(SchemaReason::WrongKind);
            break;
        }
    }
}
void clone_document(SchemaContext& context) {
    auto& target = context.output->snapshot.values()[0];
    target.identity = context.input->identity;
    target.parser_plan = context.input->parser_plan;
    target.policy = context.input->policy;
    target.source = context.allocate<char>(context.input->source.size(), false);
    target.nodes = context.allocate<InputNode>(context.input->node_count, false);
    target.edges = context.allocate<ArrayEdge>(context.input->edge_count, false);
    target.members = context.allocate<ObjectMember>(context.input->member_count, false);
    target.text = context.allocate<char>(context.input->text_count, false);
    auto copy_bytes = context.add(context.input->source.size(), context.input->text_count);
    copy_bytes =
        context.add(copy_bytes, context.multiply(context.input->node_count, sizeof(InputNode)));
    copy_bytes =
        context.add(copy_bytes, context.multiply(context.input->edge_count, sizeof(ArrayEdge)));
    copy_bytes = context.add(copy_bytes,
                             context.multiply(context.input->member_count, sizeof(ObjectMember)));
    context.work(copy_bytes);
    std::ranges::copy(context.input->source.values(), target.source.values().begin());
    std::ranges::copy(context.input->nodes.values().first(context.input->node_count),
                      target.nodes.values().begin());
    std::ranges::copy(context.input->edges.values().first(context.input->edge_count),
                      target.edges.values().begin());
    std::ranges::copy(context.input->members.values().first(context.input->member_count),
                      target.members.values().begin());
    std::ranges::copy(context.input->text.values().first(context.input->text_count),
                      target.text.values().begin());
    target.node_count = context.input->node_count;
    target.edge_count = context.input->edge_count;
    target.member_count = context.input->member_count;
    target.text_count = context.input->text_count;
}
} // namespace
std::string_view schema_text(const SchemaData* data, std::uint64_t node) noexcept {
    return node_text(data->snapshot.values()[0], node);
}
std::uint64_t schema_count(const SchemaData* data, std::uint64_t node) noexcept {
    return data->nodes.values()[static_cast<std::size_t>(node)].unsigned_value;
}
std::int64_t schema_integer(const SchemaData* data, std::uint64_t node) noexcept {
    return data->nodes.values()[static_cast<std::size_t>(node)].signed_value;
}
bool schema_boolean(const SchemaData* data, std::uint64_t node) noexcept {
    const auto* value = std::get_if<bool>(
        &data->snapshot.values()[0].nodes.values()[static_cast<std::size_t>(node)].value);
    if (!value)
        std::terminate(); // Only qualified Boolean fields can reach this getter.
    return *value;
}
core::Result<std::uint64_t> schema_child(const SchemaData* data, std::uint64_t node,
                                         core::Index index) noexcept {
    const auto& array = data->nodes.values()[static_cast<std::size_t>(node)];
    if (index < 0 || static_cast<std::uint64_t>(index) >= array.unsigned_value)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    return data->array_elements
        .values()[static_cast<std::size_t>(array.array_offset + static_cast<std::uint64_t>(index))];
}
SchemaNumberView schema_number(const SchemaData* data, std::uint64_t node) noexcept {
    const auto& typed = data->nodes.values()[static_cast<std::size_t>(node)];
    return {typed.number.sign,
            {data->digits.values().data() + typed.number.offset,
             static_cast<std::size_t>(typed.number.count)},
            typed.number.exponent,
            node_kind(data->snapshot.values()[0], node),
            typed.dimensions};
}
} // namespace uvlas::io::detail

namespace uvlas::io {
std::optional<std::string_view> SchemaDiagnosticDetail::expected_type_name() const noexcept {
    if (!expected_type)
        return std::nullopt;
    const auto index = static_cast<std::size_t>(*expected_type);
    if (index >= detail::schema_types().size())
        return std::nullopt;
    return detail::schema_types()[index].spelling;
}
std::optional<std::string_view> SchemaDiagnosticDetail::expected_variant_name() const noexcept {
    if (!expected_variant || *expected_variant >= detail::schema_records().size())
        return std::nullopt;
    return detail::schema_records()[*expected_variant].name;
}
std::optional<std::string_view> SchemaDiagnosticDetail::expected_field_name() const noexcept {
    if (!expected_variant || !expected_field ||
        *expected_variant >= detail::schema_records().size())
        return std::nullopt;
    const auto fields = detail::schema_records()[*expected_variant].fields;
    if (*expected_field >= fields.size())
        return std::nullopt;
    return fields[*expected_field].name;
}
SchemaDiagnostic::SchemaDiagnostic(SchemaDiagnostic&& other) noexcept
    : detail_(std::exchange(other.detail_, std::nullopt)), path_(std::move(other.path_)) {}
SchemaDiagnostic& SchemaDiagnostic::operator=(SchemaDiagnostic&& other) noexcept {
    if (this != &other) {
        detail_ = std::exchange(other.detail_, std::nullopt);
        path_ = std::move(other.path_);
    }
    return *this;
}
void SchemaDiagnostic::reset() noexcept {
    detail_.reset();
    path_ = {};
}
std::optional<SchemaDiagnosticDetail> SchemaDiagnostic::detail() const noexcept {
    return detail_;
}
std::optional<std::string_view> SchemaDiagnostic::path() const& noexcept {
    if (!detail_ || detail_->path_state != SchemaPathState::Complete)
        return std::nullopt;
    return std::string_view{path_.values().data(), path_.size()};
}
SchemaCase::SchemaCase(core::PlannedBuffer<detail::SchemaData> data) noexcept
    : data_(std::move(data)) {}
SchemaCase::SchemaCase(SchemaCase&&) noexcept = default;
SchemaCase& SchemaCase::operator=(SchemaCase&&) noexcept = default;
SchemaCase::~SchemaCase() = default;
core::Result<schema::Case> SchemaCase::root() const& noexcept {
    if (data_.size() == 0)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return detail::RecordAccess<detail::FieldType::T389>::read(&data_.values()[0], 0);
}
core::Result<SourceIdentity> SchemaCase::source_identity() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    const auto& source = data_.values()[0].snapshot.values()[0].identity;
    if (!source)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return *source;
}
core::Result<SchemaIdentity> SchemaCase::identity() const noexcept {
    const auto source = source_identity();
    if (!source)
        return std::unexpected(source.error());
    return SchemaIdentity{*source, 1, 1};
}
core::Result<SchemaMemoryReport> SchemaCase::memory_report() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    return data_.values()[0].report;
}
core::Result<detail::SchemaPayload>
detail::validate_schema_payload(const InputDocument& input, SchemaPayloadRequest request,
                                SchemaLimits limits, core::ResourcePlan& resources,
                                SchemaDiagnostic* diagnostic) noexcept {
    using namespace detail;
    if (diagnostic)
        diagnostic->reset();
    const auto* source = InputDocumentFactory::storage(input);
    SchemaContext context{source, limits, resources, diagnostic};
    auto result = core::guard_entry(
        core::Operation::PrepareModel, core::FailureCode::InvalidInput, {},
        [&]() -> core::Result<SchemaPayload> {
            try {
                if (!source || !source->identity)
                    context.reject(SchemaReason::InvalidOwner);
                const auto position = InputDocumentFactory::position(request.input);
                if (InputDocumentFactory::storage(request.input) != source || position < 0 ||
                    static_cast<std::uint64_t>(position) >= source->node_count ||
                    static_cast<std::size_t>(request.type) >= schema_types().size() ||
                    request.components == 0)
                    context.reject(SchemaReason::InvalidOwner);
                context.payload_root = static_cast<std::uint64_t>(position);
                context.payload_type = request.type;
                context.payload_components = request.components;
                context.payload_dimensions = request.dimensions;
                context.case_registry = request.case_registry;
                context.whole_case = context.payload_root == 0 &&
                                     request.type == schema_root_type() && !request.case_registry;
                if (source->node_count >
                    context.limits[static_cast<std::size_t>(SchemaQuota::InputNodes)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::InputNodes, source->node_count);
                if (context.whole_case || context.case_registry)
                    apply_case_limits(context);
                if (source->node_count >
                    context.limits[static_cast<std::size_t>(SchemaQuota::InputNodes)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::InputNodes, source->node_count);
                if (request.case_registry &&
                    schema_find(context, 0, "units") != context.payload_root)
                    context.reject(SchemaReason::InvalidOwner);
                context.charge(SchemaQuota::SourceBytes, source->source.size());
                context.charge(SchemaQuota::InputMembers, source->member_count);
                context.charge(SchemaQuota::TotalTextBytes, source->text_count);
                auto stack =
                    context.allocate<SchemaBootstrap>(context.add(source->node_count, 1), true);
                for (std::uint64_t node = 0; node < source->node_count; ++node) {
                    context.work();
                    const auto* array = std::get_if<ArrayValue>(
                        &source->nodes.values()[static_cast<std::size_t>(node)].value);
                    if (!array)
                        continue;
                    auto edge = array->children.head;
                    std::optional<std::uint64_t> previous;
                    while (edge >= 0) {
                        context.work();
                        const auto& item = source->edges.values()[static_cast<std::size_t>(edge)];
                        const auto child = static_cast<std::uint64_t>(item.child);
                        stack.values()[static_cast<std::size_t>(child)].previous_sibling = previous;
                        previous = child;
                        edge = item.next;
                    }
                    stack.values()[static_cast<std::size_t>(node)].last_child = previous;
                }
                SchemaData counted;
                context.output = &counted;
                shape_pass(context, stack, true);
                if (context.forecast_token_excess) {
                    const auto excess = *context.forecast_token_excess;
                    context.current = excess.node;
                    context.expected_type = excess.type;
                    context.expected_variant.reset();
                    context.expected_field.reset();
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::DecimalTokenBytes, excess.bytes);
                }
                const auto counted_records = counted.record_counts;
                std::uint64_t h = context.multiply(context.rational_digits, 2);
                if (context.arithmetic_low && context.arithmetic_high)
                    h = std::max(
                        h, context.add(static_cast<std::uint64_t>(*context.arithmetic_high) -
                                           static_cast<std::uint64_t>(*context.arithmetic_low),
                                       4));
                h = std::max<std::uint64_t>(h, 4);
                if (context.planning_overflow ||
                    h > context.limits[static_cast<std::size_t>(SchemaQuota::ExactDigits)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::ExactDigits,
                                   context.planning_overflow ? std::nullopt : std::optional{h});
                context.exact.stride = context.add(h, 1);
                auto p = context.add(sizeof(SchemaData), sizeof(InputDocumentData));
                p = context.add(p, source->source.size());
                p = context.add(p, context.multiply(source->node_count, sizeof(InputNode)));
                p = context.add(p, context.multiply(source->edge_count, sizeof(ArrayEdge)));
                p = context.add(p, context.multiply(source->member_count, sizeof(ObjectMember)));
                p = context.add(p, source->text_count);
                p = context.add(p, context.multiply(source->node_count, sizeof(SchemaNode)));
                p = context.add(p, context.multiply(source->edge_count, sizeof(std::uint64_t)));
                p = context.add(p, schema_record_bytes(context));
                p = context.add(
                    p, context.add(source->text_count, context.multiply(source->node_count, 21)));
                p = context.add(
                    p, context.multiply(context.planned_identifiers, sizeof(SchemaIdentifier)));
                p = context.add(
                    p, context.multiply(context.planned_references, sizeof(SchemaReference)));
                auto scratch =
                    context.add(context.scratch_bytes,
                                context.multiply(source->node_count,
                                                 sizeof(std::uint8_t) + sizeof(std::uint64_t)));
                scratch = context.add(scratch, context.multiply(context.exact.stride, 6));
                context.planned_owned = p;
                context.planned_scratch = scratch;
                context.planned_cumulative =
                    context.add(context.add(p, scratch), context.path_required);
                if (p > context.limits[static_cast<std::size_t>(SchemaQuota::TypedBytes)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::TypedBytes, p);
                if (scratch > context.limits[static_cast<std::size_t>(SchemaQuota::ScratchBytes)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::ScratchBytes, scratch);
                if (context.path_required >
                    context.limits[static_cast<std::size_t>(SchemaQuota::DiagnosticPathBytes)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::DiagnosticPathBytes, context.path_required);
                const auto remaining_requirement =
                    context.planned_cumulative - context.scratch_bytes;
                if (remaining_requirement > resources.remaining_bytes())
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::TypedBytes, remaining_requirement);
                auto data = context.allocate<SchemaData>(1, false);
                context.output = &data.values()[0];
                auto& owner = *context.output;
                owner.nodes = context.allocate<SchemaNode>(source->node_count, false);
                owner.snapshot = context.allocate<InputDocumentData>(1, false);
                shape_pass(context, stack, false);
                if (owner.record_counts != counted_records)
                    context.reject(SchemaReason::InvalidOwner);
                context.expected_type.reset();
                context.expected_variant.reset();
                context.expected_field.reset();
                schema_prepare_records(context);
                clone_document(context);
                owner.array_elements = context.allocate<std::uint64_t>(source->edge_count, false);
                std::uint64_t element_offset = 0;
                for (std::uint64_t node = 0; node < source->node_count; ++node) {
                    context.work();
                    const auto* array = std::get_if<ArrayValue>(
                        &source->nodes.values()[static_cast<std::size_t>(node)].value);
                    if (!array)
                        continue;
                    auto& typed = owner.nodes.values()[static_cast<std::size_t>(node)];
                    typed.array_offset = element_offset;
                    typed.unsigned_value = static_cast<std::uint64_t>(array->children.size);
                    auto edge = array->children.head;
                    while (edge >= 0) {
                        context.work(sizeof(std::uint64_t));
                        const auto& item = source->edges.values()[static_cast<std::size_t>(edge)];
                        owner.array_elements.values()[static_cast<std::size_t>(element_offset++)] =
                            static_cast<std::uint64_t>(item.child);
                        edge = item.next;
                    }
                }
                if (element_offset != source->edge_count)
                    context.reject(SchemaReason::InvalidOwner);
                const auto digit_capacity =
                    context.add(source->text_count, context.multiply(source->node_count, 21));
                owner.digits = context.allocate<char>(digit_capacity, false);
                owner.identifiers =
                    context.allocate<SchemaIdentifier>(context.planned_identifiers, false);
                owner.references =
                    context.allocate<SchemaReference>(context.planned_references, false);
                schema_exact_semantics(context);
                schema_semantics(context);
                schema_build_records(context);
                if (context.whole_case)
                    for (const auto& node : owner.nodes.values()) {
                        context.work();
                        if (!node.assigned)
                            context.reject(SchemaReason::InvalidOwner);
                    }
                owner.report = {context.owned_bytes,
                                context.scratch_bytes,
                                resources.allocated_bytes() - context.initial_charged,
                                context.used[static_cast<std::size_t>(SchemaQuota::Work)],
                                context.planned_owned,
                                context.planned_scratch,
                                context.planned_cumulative};
                if (diagnostic)
                    diagnostic->reset();
                return SchemaPayload{std::move(data), context.payload_root, context.payload_type};
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            }
        });
    if (!result && diagnostic && !diagnostic->detail()) {
        const auto reason = result.error().code == core::FailureCode::ResourceExceeded
                                ? SchemaReason::Allocation
                                : SchemaReason::UnexpectedException;
        SchemaFactory::diagnostic(
            *diagnostic,
            {reason, context.source, context.current, std::nullopt, std::nullopt, std::nullopt,
             std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt,
             context.used[static_cast<std::size_t>(SchemaQuota::Work)], SchemaPathState::NotBuilt});
    }
    return result;
}
core::Result<SchemaCase> validate_case_schema(const InputDocument& input, SchemaLimits limits,
                                              core::ResourcePlan& resources,
                                              SchemaDiagnostic* diagnostic) noexcept {
    const auto view = input.view();
    if (!view) {
        if (diagnostic) {
            diagnostic->reset();
            detail::SchemaFactory::diagnostic(
                *diagnostic,
                {SchemaReason::InvalidOwner, std::nullopt, std::nullopt, std::nullopt, std::nullopt,
                 std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt,
                 std::nullopt, 0, SchemaPathState::NotBuilt});
        }
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
    }
    auto payload = detail::validate_schema_payload(
        input, {view->root, detail::schema_root_type(), 1, SchemaDimension{{0, 0, 0, 0}}, false},
        limits, resources, diagnostic);
    if (!payload)
        return std::unexpected(payload.error());
    return detail::SchemaFactory::case_result(std::move(payload->data));
}
} // namespace uvlas::io
