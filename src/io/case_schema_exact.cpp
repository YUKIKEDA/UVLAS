#include "case_schema_storage.hpp"

#include <algorithm>
#include <charconv>
#include <limits>

namespace uvlas::io::detail {
namespace {
constexpr auto qi(SchemaQuota quota) noexcept {
    return static_cast<std::size_t>(quota);
}
std::string_view digits(const SchemaContext& context, const DecimalRange& value) noexcept {
    return {context.output->digits.values().data() + value.offset,
            static_cast<std::size_t>(value.count)};
}
std::uint64_t magnitude(std::int64_t value) noexcept {
    return value < 0 ? std::uint64_t{0} - static_cast<std::uint64_t>(value)
                     : static_cast<std::uint64_t>(value);
}
std::int64_t signed_add(SchemaContext& context, std::int64_t left, std::int64_t right) {
    if ((right > 0 && left > std::numeric_limits<std::int64_t>::max() - right) ||
        (right < 0 && left < std::numeric_limits<std::int64_t>::min() - right))
        context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
    return left + right;
}
std::int64_t signed_subtract(SchemaContext& context, std::int64_t left, std::int64_t right) {
    if ((right > 0 && left < std::numeric_limits<std::int64_t>::min() + right) ||
        (right < 0 && left > std::numeric_limits<std::int64_t>::max() + right))
        context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
    return left - right;
}
int compare_parts(SchemaContext& context, int left_sign, std::string_view left,
                  std::int64_t left_exp, int right_sign, std::string_view right,
                  std::int64_t right_exp) {
    if (left_sign != right_sign)
        return left_sign > right_sign ? 1 : -1;
    if (left_sign == 0)
        return 0;
    if (left.size() > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()) ||
        right.size() > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
    const auto left_order = signed_add(context, left_exp, static_cast<std::int64_t>(left.size()));
    const auto right_order =
        signed_add(context, right_exp, static_cast<std::int64_t>(right.size()));
    if (left_order != right_order)
        return left_sign * (left_order > right_order ? 1 : -1);
    for (std::size_t i = 0; i < std::max(left.size(), right.size()); ++i) {
        context.work();
        const auto a = i < left.size() ? left[i] : '0';
        const auto b = i < right.size() ? right[i] : '0';
        if (a != b)
            return left_sign * (a > b ? 1 : -1);
    }
    return 0;
}
bool decimal_kind(InputNodeKind kind) noexcept {
    return kind == InputNodeKind::SignedInteger || kind == InputNodeKind::UnsignedInteger ||
           kind == InputNodeKind::DecimalNumber;
}
std::uint64_t read_digits(SchemaContext& context, std::string_view text) {
    if (text.empty())
        context.reject(SchemaReason::InvalidIdentifier);
    std::uint64_t value = 0;
    for (const char c : text) {
        context.work();
        if (c < '0' || c > '9')
            context.reject(SchemaReason::InvalidIdentifier);
        const auto d = static_cast<std::uint64_t>(c - '0');
        if (value > (std::numeric_limits<std::uint64_t>::max() - d) / 10)
            context.reject(SchemaReason::InvalidIdentifier);
        value = value * 10 + d;
    }
    return value;
}
std::uint64_t enum_index(SchemaContext& context, std::string_view spelling,
                         std::string_view selected) {
    auto allowed = spelling.substr(4, spelling.size() - 5);
    if (allowed == "同じ集合")
        allowed = "geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,"
                  "algebra,newton,eigen_range,search_width,load_mapping";
    std::uint64_t index = 0;
    std::size_t start = 0;
    for (;;) {
        const auto end = allowed.find(',', start);
        const auto part = allowed.substr(start, end == std::string_view::npos ? end : end - start);
        context.work(std::min(part.size(), selected.size()) + 1);
        if (part == selected)
            return index;
        if (end == std::string_view::npos)
            break;
        start = end + 1;
        ++index;
    }
    context.reject(SchemaReason::UnsupportedTag, core::FailureCode::Unsupported);
}
std::uint8_t id_domain(std::string_view type) noexcept {
    constexpr std::array names{"Geometry",    "Material", "Shell",       "Beam",   "Joint",
                               "Support",     "Mass",     "AeroSurface", "Layer",  "Load",
                               "Observation", "Source",   "Frame",       "Branch", "WakeRow",
                               "WakePoint",   "WakeBand"};
    const auto begin = type.find('<');
    const auto end = type.find('>');
    const auto name =
        begin == std::string_view::npos ? type : type.substr(begin + 1, end - begin - 1);
    for (std::size_t i = 0; i < names.size(); ++i)
        if (name == names[i])
            return static_cast<std::uint8_t>(i);
    return 255;
}
} // namespace
DecimalRange schema_decimal(SchemaContext& context, std::uint64_t node, std::string_view token,
                            bool numerator) {
    context.current = node;
    std::size_t at = 0;
    int sign = 1;
    if (!token.empty() && token[0] == '-') {
        sign = -1;
        ++at;
    }
    if (at == token.size())
        context.reject(SchemaReason::InvalidValue);
    std::uint64_t fractional = 0;
    bool point = false;
    bool seen = false;
    auto& data = *context.output;
    const auto begin = data.digit_count;
    while (at < token.size() && token[at] != 'e' && token[at] != 'E') {
        const auto c = token[at++];
        context.work();
        if (c == '.' && !point && !numerator) {
            point = true;
            continue;
        }
        if (c < '0' || c > '9')
            context.reject(SchemaReason::InvalidValue);
        if (data.digit_count >= data.digits.size())
            context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                           SchemaQuota::TypedBytes);
        data.digits.values()[static_cast<std::size_t>(data.digit_count++)] = c;
        seen = true;
        if (point)
            fractional = context.add(fractional, 1);
    }
    if (!seen)
        context.reject(SchemaReason::InvalidValue);
    std::int64_t exponent = 0;
    if (at < token.size()) {
        if (numerator)
            context.reject(SchemaReason::InvalidValue);
        ++at;
        bool negative = false;
        if (at < token.size() && (token[at] == '+' || token[at] == '-'))
            negative = token[at++] == '-';
        if (at == token.size())
            context.reject(SchemaReason::InvalidValue);
        std::uint64_t absolute = 0;
        while (at < token.size()) {
            context.work();
            const auto c = token[at++];
            if (c < '0' || c > '9')
                context.reject(SchemaReason::InvalidValue);
            const auto digit = static_cast<std::uint64_t>(c - '0');
            if (absolute > (std::numeric_limits<std::uint64_t>::max() - digit) / 10)
                context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded,
                               SchemaQuota::ExponentMagnitude);
            absolute = absolute * 10 + digit;
            if (absolute > context.limits[qi(SchemaQuota::ExponentMagnitude)] ||
                absolute > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
                context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                               SchemaQuota::ExponentMagnitude, absolute);
        }
        exponent =
            negative ? -static_cast<std::int64_t>(absolute) : static_cast<std::int64_t>(absolute);
    }
    if (fractional > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
        context.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
    exponent = signed_add(context, exponent, -static_cast<std::int64_t>(fractional));
    auto first = begin;
    auto end = data.digit_count;
    while (first < end && data.digits.values()[static_cast<std::size_t>(first)] == '0') {
        context.work();
        ++first;
    }
    if (first == end)
        return {0, begin, 1, 0};
    while (end > first + 1 && data.digits.values()[static_cast<std::size_t>(end - 1)] == '0') {
        context.work();
        --end;
        exponent = signed_add(context, exponent, 1);
    }
    const auto count = end - first;
    if (count > context.limits[qi(SchemaQuota::ExactDigits)])
        context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                       SchemaQuota::ExactDigits, count);
    if (magnitude(exponent) > context.limits[qi(SchemaQuota::ExponentMagnitude)])
        context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                       SchemaQuota::ExponentMagnitude, magnitude(exponent));
    return {sign, first, count, exponent};
}
int schema_compare(SchemaContext& context, const DecimalRange& left, const DecimalRange& right) {
    return compare_parts(context, left.sign, digits(context, left), left.exponent, right.sign,
                         digits(context, right), right.exponent);
}
int schema_compare_integer(SchemaContext& context, const DecimalRange& left, std::int64_t right) {
    std::array<char, 32> buffer{};
    const auto converted =
        std::to_chars(buffer.data(), buffer.data() + buffer.size(), magnitude(right));
    return compare_parts(context, left.sign, digits(context, left), left.exponent,
                         right == 0 ? 0 : (right < 0 ? -1 : 1),
                         {buffer.data(), static_cast<std::size_t>(converted.ptr - buffer.data())},
                         0);
}
void schema_exact_semantics(SchemaContext& context) {
    auto& data = *context.output;
    for (std::uint64_t node = 0; node < data.nodes.size(); ++node) {
        context.work();
        context.current = node;
        auto& typed = data.nodes.values()[static_cast<std::size_t>(node)];
        if (!typed.assigned) {
            if (context.whole_case)
                context.reject(SchemaReason::WrongKind);
            continue;
        }
        const auto& type = schema_types()[static_cast<std::size_t>(typed.type)];
        const auto actual = node_kind(*context.input, node);
        if (actual == InputNodeKind::String) {
            const auto text = node_text(*context.input, node);
            if (text.size() > context.limits[qi(SchemaQuota::DecodedStringBytes)])
                context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                               SchemaQuota::DecodedStringBytes, text.size());
        }
        switch (type.kind) {
        case FieldKind::Number:
        case FieldKind::ExactQuantity: {
            if (!decimal_kind(actual))
                break;
            context.charge(SchemaQuota::Scalars, 1);
            if (actual == InputNodeKind::DecimalNumber) {
                const auto token = node_text(*context.input, node);
                if (token.size() > context.limits[qi(SchemaQuota::DecimalTokenBytes)])
                    context.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                   SchemaQuota::DecimalTokenBytes, token.size());
                typed.number = schema_decimal(context, node, token);
            } else {
                std::array<char, 32> buffer{};
                const auto& value =
                    context.input->nodes.values()[static_cast<std::size_t>(node)].value;
                auto converted = actual == InputNodeKind::SignedInteger
                                     ? std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                                     std::get<std::int64_t>(value))
                                     : std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                                     std::get<std::uint64_t>(value));
                typed.number = schema_decimal(
                    context, node,
                    {buffer.data(), static_cast<std::size_t>(converted.ptr - buffer.data())});
            }
            break;
        }
        case FieldKind::Integer: {
            context.charge(SchemaQuota::Scalars, 1);
            const auto& value = context.input->nodes.values()[static_cast<std::size_t>(node)].value;
            if (actual == InputNodeKind::SignedInteger) {
                const auto integer = std::get<std::int64_t>(value);
                typed.signed_value = integer;
                if (type.spelling != "signed-native-integer" && integer < 0)
                    context.reject(SchemaReason::InvalidValue);
                typed.unsigned_value = static_cast<std::uint64_t>(integer);
            } else {
                typed.unsigned_value = std::get<std::uint64_t>(value);
                if ((type.spelling == "Index" || type.spelling == "signed-native-integer") &&
                    typed.unsigned_value >
                        static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
                    context.reject(SchemaReason::InvalidValue);
                if (typed.unsigned_value <=
                    static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
                    typed.signed_value = static_cast<std::int64_t>(typed.unsigned_value);
            }
            if (type.spelling == "PositiveCount" && typed.unsigned_value == 0)
                context.reject(SchemaReason::InvalidValue);
            break;
        }
        case FieldKind::Identifier:
        case FieldKind::Reference: {
            context.charge(SchemaQuota::Scalars, 1);
            typed.unsigned_value = read_digits(context, node_text(*context.input, node));
            auto domain = id_domain(type.spelling);
            if (domain == 255 && typed.parent) {
                const auto kind = schema_find(context, *typed.parent, "kind");
                if (kind) {
                    const auto tag = node_text(*context.input, *kind);
                    domain = tag == "shell" ? 2 : (tag == "beam" ? 3 : 6);
                }
            }
            if (domain == 255)
                context.reject(SchemaReason::InvalidIdentifier);
            typed.reference_domain = static_cast<SchemaNamespace>(domain);
            // Registration is deferred until every exact scalar has passed phase three.
            break;
        }
        case FieldKind::Enumeration:
            context.charge(SchemaQuota::Scalars, 1);
            typed.unsigned_value =
                enum_index(context, type.spelling, node_text(*context.input, node));
            break;
        case FieldKind::Boolean:
            context.charge(SchemaQuota::Scalars, 1);
            break;
        case FieldKind::Text: {
            const auto text = node_text(*context.input, node);
            if (type.spelling == "Digest") {
                if (text.size() != 64)
                    context.reject(SchemaReason::InvalidValue);
                for (const char c : text) {
                    context.work();
                    if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')))
                        context.reject(SchemaReason::InvalidValue);
                }
            }
            break;
        }
        default:
            break;
        }
    }
    // Rational values are parsed only after the complete tree kind/shape phase.
    for (std::uint64_t node = 0; node < data.nodes.size(); ++node) {
        context.work();
        context.current = node;
        auto& typed = data.nodes.values()[static_cast<std::size_t>(node)];
        if (!typed.assigned || node_kind(*context.input, node) != InputNodeKind::Object)
            continue;
        const auto& selected = schema_records()[typed.record_kind];
        if (selected.family != "Rational")
            continue;
        const auto numerator = schema_find(context, node, "numerator");
        const auto denominator = schema_find(context, node, "denominator");
        if (!numerator || !denominator)
            context.reject(SchemaReason::MissingField);
        context.charge(SchemaQuota::Scalars, 2);
        typed.rational_numerator =
            schema_decimal(context, *numerator, node_text(*context.input, *numerator), true);
        typed.rational_denominator =
            schema_decimal(context, *denominator, node_text(*context.input, *denominator), true);
        if (typed.rational_denominator.sign <= 0)
            context.reject(SchemaReason::InvalidValue);
    }
}
} // namespace uvlas::io::detail

namespace uvlas::io::detail {
namespace {
class SemanticCheck;
int compare_exact(SemanticCheck&, std::uint64_t, std::uint64_t);
class SemanticCheck {
public:
    explicit SemanticCheck(SchemaContext& context)
        : c(context), d(*context.input), s(*context.output) {}
    SchemaContext& c;
    const InputDocumentData& d;
    SchemaData& s;
    std::uint64_t field(std::uint64_t node, std::string_view name) {
        const auto value = schema_find(c, node, name);
        if (!value)
            c.reject(SchemaReason::MissingField);
        return *value;
    }
    std::optional<std::uint64_t> maybe(std::uint64_t node, std::string_view name) {
        return schema_find(c, node, name);
    }
    [[nodiscard]] std::string_view text(std::uint64_t node) const noexcept {
        return node_text(d, node);
    }
    std::string_view tag(std::uint64_t node, std::string_view key = "kind") {
        return text(field(node, key));
    }
    [[nodiscard]] std::uint64_t count(std::uint64_t node) const noexcept {
        return s.nodes.values()[static_cast<std::size_t>(node)].unsigned_value;
    }
    [[nodiscard]] std::uint64_t size(std::uint64_t node) const noexcept {
        return s.nodes.values()[static_cast<std::size_t>(node)].unsigned_value;
    }
    std::uint64_t at(std::uint64_t node, std::uint64_t index) {
        c.work();
        const auto& array = s.nodes.values()[static_cast<std::size_t>(node)];
        return s.array_elements.values()[static_cast<std::size_t>(array.array_offset + index)];
    }
    [[nodiscard]] const DecimalRange& number(std::uint64_t node) const noexcept {
        return s.nodes.values()[static_cast<std::size_t>(node)].number;
    }
    [[nodiscard]] const RecordDefinition& record(std::uint64_t node) const noexcept {
        return schema_records()[s.nodes.values()[static_cast<std::size_t>(node)].record_kind];
    }
    void require(bool valid, SchemaReason reason = SchemaReason::InvalidValue) {
        if (!valid)
            c.reject(reason);
    }
    void nonempty(std::uint64_t node) {
        require(size(node) > 0, SchemaReason::WrongShape);
    }
    void length(std::uint64_t node, std::uint64_t expected) {
        require(size(node) == expected, SchemaReason::WrongShape);
    }
    void numeric_bound(std::uint64_t node, std::int64_t bound, bool strict) {
        c.current = node;
        const auto order = schema_compare_integer(c, number(node), bound);
        require(strict ? order > 0 : order >= 0);
    }
    void each_bound(std::uint64_t node, std::int64_t bound, bool strict) {
        for (std::uint64_t j = 0; j < size(node); ++j)
            numeric_bound(at(node, j), bound, strict);
    }
    void order(std::uint64_t node, bool strict) {
        for (std::uint64_t i = 1; i < size(node); ++i) {
            const auto comparison = schema_compare(c, number(at(node, i - 1)), number(at(node, i)));
            require(strict ? comparison < 0 : comparison <= 0);
        }
    }
    void interval(std::uint64_t node, bool strict) {
        length(node, 2);
        const auto comparison = schema_compare(c, number(at(node, 0)), number(at(node, 1)));
        require(strict ? comparison < 0 : comparison <= 0);
    }
    void text_not_empty(std::uint64_t node, bool path = false) {
        const auto value = text(node);
        require(!value.empty());
        if (path) {
            for (const char ch : value) {
                c.work();
                require(ch != '\0');
            }
        }
    }
    bool same_scalar(std::uint64_t first, std::uint64_t second) {
        const auto left = node_kind(d, first), right = node_kind(d, second);
        const auto& lt = s.nodes.values()[static_cast<std::size_t>(first)];
        const auto& rt = s.nodes.values()[static_cast<std::size_t>(second)];
        const auto lkind = schema_types()[static_cast<std::size_t>(lt.type)].kind;
        const auto rkind = schema_types()[static_cast<std::size_t>(rt.type)].kind;
        if ((lkind == FieldKind::Number || lkind == FieldKind::ExactQuantity) &&
            (rkind == FieldKind::Number || rkind == FieldKind::ExactQuantity) &&
            decimal_kind(left) && decimal_kind(right))
            return schema_compare(c, number(first), number(second)) == 0;
        if (lkind == FieldKind::Integer && rkind == FieldKind::Integer) {
            if (schema_types()[static_cast<std::size_t>(lt.type)].spelling ==
                    "signed-native-integer" ||
                schema_types()[static_cast<std::size_t>(rt.type)].spelling ==
                    "signed-native-integer")
                return lt.signed_value == rt.signed_value;
            return lt.unsigned_value == rt.unsigned_value;
        }
        if ((lkind == FieldKind::Identifier || lkind == FieldKind::Reference) &&
            (rkind == FieldKind::Identifier || rkind == FieldKind::Reference))
            return lt.reference_domain == rt.reference_domain &&
                   lt.unsigned_value == rt.unsigned_value;
        if (left != right)
            return false;
        if (left == InputNodeKind::String) {
            c.work(std::min(text(first).size(), text(second).size()) + 1);
            return text(first) == text(second);
        }
        if (left == InputNodeKind::Boolean)
            return std::get<bool>(d.nodes.values()[static_cast<std::size_t>(first)].value) ==
                   std::get<bool>(d.nodes.values()[static_cast<std::size_t>(second)].value);
        return false;
    }
    bool same_value(std::uint64_t left, std::uint64_t right) {
        c.work();
        const auto lk = node_kind(d, left), rk = node_kind(d, right);
        const bool left_rational = lk == InputNodeKind::Object && record(left).family == "Rational";
        const bool right_rational =
            rk == InputNodeKind::Object && record(right).family == "Rational";
        const auto left_type =
            schema_types()[static_cast<std::size_t>(
                               s.nodes.values()[static_cast<std::size_t>(left)].type)]
                .kind;
        const auto right_type =
            schema_types()[static_cast<std::size_t>(
                               s.nodes.values()[static_cast<std::size_t>(right)].type)]
                .kind;
        const bool left_number =
            (left_type == FieldKind::Number || left_type == FieldKind::ExactQuantity) &&
            decimal_kind(lk);
        const bool right_number =
            (right_type == FieldKind::Number || right_type == FieldKind::ExactQuantity) &&
            decimal_kind(rk);
        if ((left_rational || left_number) && (right_rational || right_number) &&
            (left_rational || right_rational))
            return compare_exact(*this, left, right) == 0;
        if (lk == InputNodeKind::Array && rk == InputNodeKind::Array) {
            if (size(left) != size(right))
                return false;
            for (std::uint64_t i = 0; i < size(left); ++i)
                if (!same_value(at(left, i), at(right, i)))
                    return false;
            return true;
        }
        if (lk == InputNodeKind::Object && rk == InputNodeKind::Object) {
            const auto& a =
                std::get<ObjectValue>(d.nodes.values()[static_cast<std::size_t>(left)].value);
            const auto& b =
                std::get<ObjectValue>(d.nodes.values()[static_cast<std::size_t>(right)].value);
            if (a.members.size != b.members.size)
                return false;
            auto member = a.members.head;
            while (member >= 0) {
                c.work();
                const auto& entry = d.members.values()[static_cast<std::size_t>(member)];
                const auto key =
                    std::string_view{d.text.values().data() + entry.key.offset, entry.key.size};
                const auto other = maybe(right, key);
                if (!other || !same_value(static_cast<std::uint64_t>(entry.child), *other))
                    return false;
                member = entry.next;
            }
            return true;
        }
        return same_scalar(left, right);
    }
    void unique(std::uint64_t array) {
        for (std::uint64_t j = 0; j < size(array); ++j)
            for (std::uint64_t k = 0; k < j; ++k)
                require(!same_value(at(array, j), at(array, k)));
    }
    void powers(std::uint64_t terms,
                std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max()) {
        for (std::uint64_t j = 0; j < size(terms); ++j) {
            const auto item = at(terms, j);
            const auto values = field(item, "powers");
            std::uint64_t total = 0;
            for (std::uint64_t k = 0; k < size(values); ++k) {
                const auto value = count(at(values, k));
                require(value <= maximum);
                total = c.add(total, value);
            }
            require(total <= maximum);
            for (std::uint64_t k = 0; k < j; ++k)
                require(!same_value(values, field(at(terms, k), "powers")));
        }
    }
    std::uint64_t referent(std::uint64_t reference) {
        const auto target = s.nodes.values()[static_cast<std::size_t>(reference)].reference_node;
        if (!target)
            c.reject(SchemaReason::MissingReference);
        return *target;
    }
    std::uint64_t model_list(std::string_view name) {
        return field(field(0, "model"), name);
    }
    void complete_states(std::uint64_t states, std::uint64_t model, std::string_view key,
                         bool only_rigid = false) {
        std::uint64_t expected = 0;
        for (std::uint64_t j = 0; j < size(model); ++j) {
            const auto entity = at(model, j);
            if (only_rigid && tag(entity) != "rigid")
                continue;
            ++expected;
            const auto id = count(field(entity, "id"));
            std::uint64_t found = 0;
            for (std::uint64_t k = 0; k < size(states); ++k)
                if (count(field(at(states, k), key)) == id)
                    ++found;
            require(found == 1);
        }
        require(size(states) == expected);
    }
};

void resolve_references(SemanticCheck& check) {
    auto& c = check.c;
    auto& owner = check.s;
    for (std::uint64_t node = 0; node < owner.nodes.size(); ++node) {
        c.work();
        c.current = node;
        auto& typed = owner.nodes.values()[static_cast<std::size_t>(node)];
        if (!typed.assigned)
            continue;
        const auto kind = schema_types()[static_cast<std::size_t>(typed.type)].kind;
        if (kind != FieldKind::Identifier && kind != FieldKind::Reference)
            continue;
        if (!typed.reference_domain)
            c.reject(SchemaReason::InvalidIdentifier);
        const auto domain = static_cast<std::uint8_t>(*typed.reference_domain);
        if (kind == FieldKind::Reference) {
            c.charge(SchemaQuota::References, 1);
            owner.references.values()[static_cast<std::size_t>(owner.reference_count++)] = {
                domain, typed.unsigned_value, node};
            continue;
        }
        bool repeated_source = false;
        if (*typed.reference_domain == SchemaNamespace::Source) {
            for (std::uint64_t k = 0; k < owner.identifier_count; ++k) {
                c.work();
                const auto& prior = owner.identifiers.values()[static_cast<std::size_t>(k)];
                if (prior.domain != domain || prior.value != typed.unsigned_value)
                    continue;
                const auto first =
                    owner.nodes.values()[static_cast<std::size_t>(prior.node)].parent;
                if (!first || !typed.parent)
                    c.reject(SchemaReason::DuplicateIdentifier);
                for (auto key : {"path", "sha256", "role"})
                    if (!check.same_scalar(check.field(*first, key),
                                           check.field(*typed.parent, key)))
                        c.reject(SchemaReason::DuplicateIdentifier);
                typed.reference_node = *first;
                repeated_source = true;
                break;
            }
        }
        if (!repeated_source) {
            c.charge(SchemaQuota::Entities, 1);
            owner.identifiers.values()[static_cast<std::size_t>(owner.identifier_count++)] = {
                domain, typed.unsigned_value, node};
        }
    }
    auto entries =
        owner.identifiers.values().first(static_cast<std::size_t>(owner.identifier_count));
    const auto less = [&](const SchemaIdentifier& a, const SchemaIdentifier& b) {
        c.work(4);
        if (a.domain != b.domain)
            return a.domain < b.domain;
        if (a.value != b.value)
            return a.value < b.value;
        return a.node < b.node;
    };
    const auto exchange = [&](std::size_t a, std::size_t b) {
        c.work(3 * sizeof(SchemaIdentifier));
        std::swap(entries[a], entries[b]);
    };
    const auto sift = [&](std::size_t root, std::size_t count) {
        // root<count/2 proves both 2*root+1 and the optional right child fit.
        while (root < count / 2) {
            c.work();
            auto child = root * 2 + 1;
            if (child + 1 < count && less(entries[child], entries[child + 1]))
                ++child;
            if (!less(entries[root], entries[child]))
                break;
            exchange(root, child);
            root = child;
        }
    };
    for (std::size_t j = entries.size() / 2; j > 0; --j)
        sift(j - 1, entries.size());
    for (std::size_t count = entries.size(); count > 1; --count) {
        exchange(0, count - 1);
        sift(0, count - 1);
    }
    for (std::size_t j = 1; j < entries.size(); ++j) {
        const auto& previous = entries[j - 1];
        const auto& current = entries[j];
        if (previous.domain != current.domain || previous.value != current.value)
            continue;
        c.current = current.node;
        if (current.domain != static_cast<std::uint8_t>(SchemaNamespace::Source))
            c.reject(SchemaReason::DuplicateIdentifier);
        const auto pa = owner.nodes.values()[static_cast<std::size_t>(previous.node)].parent;
        const auto pb = owner.nodes.values()[static_cast<std::size_t>(current.node)].parent;
        if (!pa || !pb)
            c.reject(SchemaReason::DuplicateIdentifier);
        for (auto key : {"path", "sha256", "role"})
            if (!check.same_scalar(check.field(*pa, key), check.field(*pb, key)))
                c.reject(SchemaReason::DuplicateIdentifier);
    }
    for (std::uint64_t j = 0; j < owner.reference_count; ++j) {
        c.work();
        const auto reference = owner.references.values()[static_cast<std::size_t>(j)];
        c.current = reference.node;
        std::size_t low = 0, high = entries.size();
        while (low < high) {
            c.work();
            const auto middle = low + (high - low) / 2;
            const auto& item = entries[middle];
            if (item.domain < reference.domain ||
                (item.domain == reference.domain && item.value < reference.value))
                low = middle + 1;
            else
                high = middle;
        }
        if (low == entries.size() || entries[low].domain != reference.domain ||
            entries[low].value != reference.value)
            c.reject(SchemaReason::MissingReference);
        const auto definition =
            owner.nodes.values()[static_cast<std::size_t>(entries[low].node)].parent;
        if (!definition)
            c.reject(SchemaReason::MissingReference);
        owner.nodes.values()[static_cast<std::size_t>(reference.node)].reference_node = *definition;
        const auto parent = owner.nodes.values()[static_cast<std::size_t>(reference.node)].parent;
        if (parent && check.record(*parent).family == "Target" && check.tag(*parent) == "rigid")
            check.require(check.tag(*definition) == "rigid");
    }
}

void expression_checks(SemanticCheck& check, core::PlannedBuffer<std::uint8_t>& visited,
                       core::PlannedBuffer<std::uint64_t>& depths) {
    auto& c = check.c;
    for (std::uint64_t expression = 0; expression < check.s.nodes.size(); ++expression) {
        c.work();
        if (!check.s.nodes.values()[static_cast<std::size_t>(expression)].assigned ||
            node_kind(check.d, expression) != InputNodeKind::Object ||
            check.record(expression).family != "Expression")
            continue;
        c.current = expression;
        const auto nodes = check.field(expression, "nodes");
        const auto outputs = check.field(expression, "outputs");
        check.nonempty(nodes);
        const auto domain = check.field(expression, "domain");
        const auto domain_kind = check.tag(domain);
        const auto parent = check.s.nodes.values()[static_cast<std::size_t>(expression)].parent;
        if (parent) {
            const auto family = check.record(*parent).family;
            if (family == "ReferenceField")
                check.require(domain_kind == "reference", SchemaReason::InvalidExpression);
            if (family == "TimeSeries")
                check.require(domain_kind == "time", SchemaReason::InvalidExpression);
            if (family == "Load")
                check.require(domain_kind == "current_position_time",
                              SchemaReason::InvalidExpression);
        }
        for (std::uint64_t j = 0; j < check.size(nodes); ++j) {
            visited.values()[static_cast<std::size_t>(j)] = 0;
            depths.values()[static_cast<std::size_t>(j)] = 1;
        }
        for (std::uint64_t j = 0; j < check.size(nodes); ++j) {
            const auto node = check.at(nodes, j);
            c.current = node;
            c.charge(SchemaQuota::ExpressionNodes, 1);
            const auto op = check.tag(node, "op");
            SchemaDimension dimension{{0, 0, 0, 0}};
            std::array<std::uint64_t, 2> operands{};
            std::size_t arity = 0;
            if (op == "constant") {
                const auto dimensions = check.field(node, "dimensions");
                for (std::size_t k = 0; k < 4; ++k) {
                    const auto value =
                        check.s.nodes.values()[static_cast<std::size_t>(check.at(dimensions, k))]
                            .signed_value;
                    check.require(value >= std::numeric_limits<std::int16_t>::min() &&
                                  value <= std::numeric_limits<std::int16_t>::max());
                    dimension.exponents[k] = static_cast<std::int16_t>(value);
                }
            } else if (op == "variable") {
                const auto axis = check.count(check.field(node, "axis"));
                if (domain_kind == "time") {
                    check.require(axis == 0);
                    dimension.exponents[2] = 1;
                } else if (domain_kind == "reference") {
                    check.require(axis < 3);
                    if (check.tag(domain, "coordinates") != "chart")
                        dimension.exponents[0] = 1;
                } else {
                    check.require(axis < 4);
                    dimension.exponents[axis == 3 ? 2 : 0] = 1;
                }
            } else {
                if (op == "add" || op == "subtract" || op == "multiply" || op == "divide") {
                    operands = {check.count(check.field(node, "first")),
                                check.count(check.field(node, "second"))};
                    arity = 2;
                } else {
                    operands[0] = check.count(check.field(node, "operand"));
                    arity = 1;
                }
                for (std::size_t k = 0; k < arity; ++k) {
                    check.require(operands[k] < j, SchemaReason::InvalidExpression);
                    c.charge(SchemaQuota::ExpressionEdges, 1);
                    depths.values()[static_cast<std::size_t>(j)] =
                        std::max(depths.values()[static_cast<std::size_t>(j)],
                                 c.add(depths.values()[static_cast<std::size_t>(operands[k])], 1));
                }
                dimension =
                    check.s.nodes.values()[static_cast<std::size_t>(check.at(nodes, operands[0]))]
                        .dimensions;
                SchemaDimension second{{0, 0, 0, 0}};
                if (arity == 2)
                    second = check.s.nodes
                                 .values()[static_cast<std::size_t>(check.at(nodes, operands[1]))]
                                 .dimensions;
                if (op == "add" || op == "subtract")
                    check.require(dimension == second, SchemaReason::DimensionMismatch);
                if (op == "multiply" || op == "divide")
                    for (std::size_t k = 0; k < 4; ++k) {
                        c.work();
                        const auto value =
                            static_cast<std::int32_t>(dimension.exponents[k]) +
                            (op == "multiply" ? static_cast<std::int32_t>(second.exponents[k])
                                              : -static_cast<std::int32_t>(second.exponents[k]));
                        if (value < std::numeric_limits<std::int16_t>::min() ||
                            value > std::numeric_limits<std::int16_t>::max())
                            c.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
                        dimension.exponents[k] = static_cast<std::int16_t>(value);
                    }
                if (op == "integer_power") {
                    const auto exponent =
                        check.s.nodes
                            .values()[static_cast<std::size_t>(check.field(node, "exponent"))]
                            .signed_value;
                    check.require(exponent >= std::numeric_limits<std::int32_t>::min() &&
                                  exponent <= std::numeric_limits<std::int32_t>::max());
                    const auto absolute = magnitude(exponent);
                    if (absolute > c.limits[qi(SchemaQuota::ExpressionPower)])
                        c.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                                 SchemaQuota::ExpressionPower, absolute);
                    for (auto& power : dimension.exponents) {
                        c.work();
                        const auto value = static_cast<std::int64_t>(power) * exponent;
                        if (value < std::numeric_limits<std::int16_t>::min() ||
                            value > std::numeric_limits<std::int16_t>::max())
                            c.reject(SchemaReason::Overflow, core::FailureCode::ResourceExceeded);
                        power = static_cast<std::int16_t>(value);
                    }
                }
                if (op == "sqrt")
                    for (auto& power : dimension.exponents) {
                        c.work();
                        check.require(power % 2 == 0, SchemaReason::DimensionMismatch);
                        power = static_cast<std::int16_t>(power / 2);
                    }
                if (op == "sin" || op == "cos") {
                    check.require(dimension == SchemaDimension{{0, 0, 0, 1}},
                                  SchemaReason::DimensionMismatch);
                    dimension = {{0, 0, 0, 0}};
                }
                if (op == "exp" || op == "log")
                    check.require(dimension == SchemaDimension{{0, 0, 0, 0}},
                                  SchemaReason::DimensionMismatch);
            }
            const auto depth = depths.values()[static_cast<std::size_t>(j)];
            if (depth > c.limits[qi(SchemaQuota::ExpressionDepth)])
                c.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                         SchemaQuota::ExpressionDepth, depth);
            check.s.nodes.values()[static_cast<std::size_t>(node)].dimensions = dimension;
            if (op == "constant")
                check.s.nodes.values()[static_cast<std::size_t>(check.field(node, "value"))]
                    .dimensions = dimension;
        }
        for (std::uint64_t k = 0; k < check.size(outputs); ++k) {
            const auto index = check.count(check.at(outputs, k));
            check.require(index < check.size(nodes), SchemaReason::InvalidExpression);
            visited.values()[static_cast<std::size_t>(index)] = 1;
            check.require(
                check.s.nodes.values()[static_cast<std::size_t>(check.at(nodes, index))]
                        .dimensions ==
                    check.s.nodes.values()[static_cast<std::size_t>(expression)].dimensions,
                SchemaReason::DimensionMismatch);
        }
        for (std::uint64_t j = check.size(nodes); j > 0; --j) {
            const auto index = j - 1;
            check.require(visited.values()[static_cast<std::size_t>(index)] != 0,
                          SchemaReason::InvalidExpression);
            const auto node = check.at(nodes, index);
            for (auto key : {"first", "second", "operand"}) {
                const auto operand = check.maybe(node, key);
                if (operand)
                    visited.values()[static_cast<std::size_t>(check.count(*operand))] = 1;
            }
        }
    }
}

void generic_record_checks(SemanticCheck& check, std::uint64_t node) {
    auto& c = check.c;
    c.current = node;
    const auto& record = check.record(node);
    const auto family = record.family;
    const auto present = [&](std::string_view key) { return check.maybe(node, key); };
    const auto nonnegative = [&](std::string_view key) {
        const auto value = present(key);
        if (value)
            check.numeric_bound(*value, 0, false);
    };
    const auto positive = [&](std::string_view key) {
        const auto value = present(key);
        if (value) {
            if (node_kind(check.d, *value) == InputNodeKind::Array)
                check.each_bound(*value, 0, true);
            else
                check.numeric_bound(*value, 0, true);
        }
    };
    if (family == "Rational")
        return;
    if (family == "FileReference") {
        check.text_not_empty(check.field(node, "path"), true);
        const auto& spelling =
            schema_types()[static_cast<std::size_t>(
                               check.s.nodes.values()[static_cast<std::size_t>(node)].type)]
                .spelling;
        if (spelling.starts_with("FileRef<"))
            check.require(check.text(check.field(node, "role")) ==
                          spelling.substr(8, spelling.size() - 9));
        return;
    }
    if (family == "PolyTerm" || family == "PgTerm" || family == "MapTerm" || family == "VectorTerm")
        c.charge(SchemaQuota::PolynomialTerms, 1);
    if (family == "Geometry") {
        const auto kind = check.tag(node);
        if (kind == "pg" || kind == "surface_graph") {
            positive("reference_length");
            check.require(
                check.s.nodes.values()[static_cast<std::size_t>(check.field(node, "width"))]
                    .rational_numerator.sign > 0);
            const auto degree = check.field(node, "degree");
            for (std::uint64_t k = 0; k < 2; ++k)
                check.require(check.count(check.at(degree, k)) <= 4);
            const auto terms = check.field(node, "coefficients");
            check.powers(terms, 4);
            for (std::uint64_t j = 0; j < check.size(terms); ++j) {
                const auto powers = check.field(check.at(terms, j), "powers");
                for (std::uint64_t k = 0; k < 2; ++k)
                    check.require(check.count(check.at(powers, k)) <=
                                  check.count(check.at(degree, k)));
            }
            const auto domain = check.field(node, "domain");
            check.require(check.size(domain) >= 3);
            check.require(
                !check.same_value(check.at(domain, 0), check.at(domain, check.size(domain) - 1)));
            if (kind == "pg") {
                const auto n = check.count(check.field(node, "n"));
                check.require(n == 1 || n == 2 || n == 4);
            }
            const auto boundaries = check.field(node, "boundaries");
            for (std::uint64_t j = 0; j < check.size(boundaries); ++j)
                for (std::uint64_t k = 0; k < j; ++k)
                    check.require(!check.same_scalar(check.field(check.at(boundaries, j), "name"),
                                                     check.field(check.at(boundaries, k), "name")));
        }
    }
    if (family == "Thickness") {
        const auto vertices = check.field(node, "vertices"), values = check.field(node, "values"),
                   triangles = check.field(node, "triangles");
        check.require(check.size(vertices) >= 3);
        check.length(values, check.size(vertices));
        check.each_bound(values, 0, true);
        for (std::uint64_t j = 0; j < check.size(triangles); ++j) {
            const auto triangle = check.at(triangles, j);
            check.unique(triangle);
            for (std::uint64_t k = 0; k < 3; ++k)
                check.require(check.count(check.at(triangle, k)) < check.size(vertices));
        }
    }
    if (family == "Boundary") {
        check.text_not_empty(check.field(node, "name"), true);
        check.require(check.size(check.field(node, "vertices")) >= 2);
    }
    if (family == "ReferenceGeometry" && check.tag(node) == "fit") {
        check.require(check.count(check.field(node, "degree")) <= 4);
        positive("distance_budget");
        positive("normal_budget");
        check.nonempty(check.field(node, "conversion_observations"));
    }
    if (family == "ReferenceField") {
        const auto kind = check.tag(node);
        if (kind == "polynomial") {
            check.each_bound(check.field(node, "lengths"), 0, true);
            check.each_bound(check.field(node, "error"), 0, false);
            const auto components = check.field(node, "components");
            for (std::uint64_t j = 0; j < check.size(components); ++j)
                check.powers(check.at(components, j));
        }
        if (kind == "expression")
            check.each_bound(check.field(node, "error"), 0, false);
    }
    if (family == "TimeSeries") {
        const auto kind = check.tag(node);
        if (kind == "polynomial") {
            positive("scale");
            check.interval(check.field(node, "domain"), true);
            const auto rows = check.field(node, "coefficients");
            check.nonempty(rows);
            c.charge(SchemaQuota::TableSamples, check.size(rows));
        }
        if (kind == "pchip") {
            const auto times = check.field(node, "times");
            check.require(check.size(times) >= 2);
            check.order(times, true);
            check.length(check.field(node, "values"), check.size(times));
            c.charge(SchemaQuota::TableSamples, check.size(times));
        }
    }
    if (family == "Outside")
        positive("period");
    if (family == "Motion" || family == "Flow") {
        const auto events = present("events");
        if (events)
            check.order(*events, true);
    }
    if (family == "Material" || family == "Section") {
        for (auto key : {"young", "density", "radius", "line_density", "torsion"})
            positive(key);
        for (auto key : {"young", "shear", "sides", "effective_shear"}) {
            const auto value = present(key);
            if (value && node_kind(check.d, *value) == InputNodeKind::Array)
                check.each_bound(*value, 0, true);
        }
        const auto poisson = present("poisson");
        if (poisson && node_kind(check.d, *poisson) != InputNodeKind::Array) {
            check.numeric_bound(*poisson, -1, true);
            check.require(compare_parts(c, check.number(*poisson).sign,
                                        digits(c, check.number(*poisson)),
                                        check.number(*poisson).exponent, 1, "5", -1) < 0);
        }
    }
    if (family == "Applicability" && check.tag(node) == "declared") {
        const auto lower = check.field(node, "strain_lower"),
                   upper = check.field(node, "strain_upper");
        for (std::uint64_t k = 0; k < 6; ++k)
            check.require(schema_compare(c, check.number(check.at(lower, k)),
                                         check.number(check.at(upper, k))) <= 0);
    }
    if (family == "Shell")
        check.nonempty(check.field(node, "layers"));
    if (family == "AnalysisSpeed")
        nonnegative("value");
    if (family == "Analysis" && check.tag(node) == "flutter") {
        check.nonempty(check.field(node, "branches"));
        check.nonempty(check.field(node, "intervals"));
        positive("step");
    }
    if (family == "Region" && check.tag(node) == "point")
        check.require(check.count(check.field(node, "boundary")) == 0);
    if (family == "CoefficientState") {
        const auto position = check.field(node, "position"),
                   velocity = check.field(node, "velocity");
        check.length(check.field(velocity, "high"), check.size(check.field(position, "high")));
    }
    if (family == "BeamCoordinate") {
        positive("length");
        positive("arc_scale");
    }
    if (family == "Joint" && check.tag(node) == "adhesive")
        check.nonempty(check.field(node, "mapping"));
    if (family == "TraceMap") {
        check.require(check.size(check.field(node, "domain")) >= 3);
        for (auto key : {"first", "second"}) {
            const auto components = check.field(node, key);
            for (std::uint64_t j = 0; j < check.size(components); ++j)
                check.powers(check.at(components, j));
        }
    }
    if (family == "Support") {
        const auto components = present("components");
        if (components) {
            check.nonempty(*components);
            check.unique(*components);
            const auto maximum = check.tag(node) == "trace_constraint" ? 2U : 5U;
            for (std::uint64_t k = 0; k < check.size(*components); ++k)
                check.require(check.count(check.at(*components, k)) <= maximum);
        }
    }
    if (family == "Support" && check.tag(node) == "constraint") {
        const auto target = check.field(check.field(node, "target"), "target");
        if (check.tag(target) == "shell") {
            const auto components = check.field(node, "components");
            for (std::uint64_t k = 0; k < check.size(components); ++k)
                check.require(check.count(check.at(components, k)) < 3);
        }
    }
    if (family == "Load" && check.tag(node) == "couple")
        check.require(check.tag(check.field(check.field(node, "target"), "target")) != "shell");
    if (family == "Mass")
        positive("mass");
    if (family == "Ledger" && check.tag(node) == "replace") {
        check.nonempty(check.field(node, "replaced"));
        check.unique(check.field(node, "replaced"));
    }
    if (family == "AeroSurface") {
        for (auto key : {"chord", "span"}) {
            const auto values = check.field(node, key);
            check.require(check.size(values) >= 2);
            check.order(values, true);
        }
    }
    if (family == "Flow")
        positive("density");
    if (family == "FlowReference") {
        for (auto key : {"density", "area", "chord", "span"})
            positive(key);
        nonnegative("speed");
    }
    if (family == "Wake")
        positive("core_radius");
    if (family == "Retention") {
        positive("duration");
        positive("distance");
        positive("error_budget");
    }
    if (family == "Continuation") {
        const auto minimum = present("minimum_increment"), initial = present("initial_increment"),
                   maximum = present("maximum_increment");
        if (minimum && initial && maximum) {
            check.numeric_bound(*minimum, 0, true);
            check.require(schema_compare(c, check.number(*minimum), check.number(*initial)) <= 0);
            check.require(schema_compare(c, check.number(*initial), check.number(*maximum)) <= 0);
        }
    }
    if (family == "ArcLength")
        positive("speed_scale");
    if (family == "TimeInterval") {
        const auto begin = check.field(node, "start"), end = check.field(node, "end");
        check.require(schema_compare(c, check.number(begin), check.number(end)) < 0);
        positive("step");
    }
    if (family == "SpeedInterval") {
        nonnegative("lower");
        check.require(schema_compare(c, check.number(check.field(node, "lower")),
                                     check.number(check.field(node, "upper"))) < 0);
    }
    if (family == "Search") {
        for (auto key : {"speed_absolute", "speed_relative", "growth_budget"})
            nonnegative(key);
    }
    if (family == "Association") {
        check.nonempty(check.field(node, "physical_observations"));
        for (auto key : {"mac_minimum", "subspace_angle_maximum", "eigenvalue_absolute",
                         "eigenvalue_relative"})
            nonnegative(key);
        check.require(
            schema_compare_integer(c, check.number(check.field(node, "mac_minimum")), 1) <= 0);
    }
    if (family == "Coverage") {
        check.nonempty(check.field(node, "branches"));
        check.nonempty(check.field(node, "mode_observations"));
    }

    if (family == "Planes" && check.tag(node) == "two_word")
        check.length(check.field(node, "low"), check.size(check.field(node, "high")));

    if (family == "Restart") {
        check.text_not_empty(check.field(node, "directory"), true);
        const auto budgets = check.field(node, "budget");
        check.require(check.tag(node, "compatibility") == "same_environment_bitwise"
                          ? check.size(budgets) == 0
                          : check.size(budgets) > 0);
    }
    if (family == "FullShellField") {
        if (check.tag(node) == "half_cells") {
            const auto cells = check.field(node, "cells");
            for (std::uint64_t j = 0; j < check.size(cells); ++j)
                for (std::uint64_t k = 0; k < j; ++k)
                    check.require(!check.same_value(check.field(check.at(cells, j), "cell"),
                                                    check.field(check.at(cells, k), "cell")));
        }
    }
    if (family == "VectorPolynomial") {
        std::uint64_t maximum = 21;
        const auto parent = check.s.nodes.values()[static_cast<std::size_t>(node)].parent;
        if (parent)
            for (auto name : {"r", "p", "c"}) {
                const auto candidate = check.maybe(*parent, name);
                if (candidate && *candidate == node)
                    maximum = name == std::string_view{"r"}
                                  ? 13
                                  : (name == std::string_view{"p"} ? 21 : 2);
            }
        check.powers(check.field(node, "terms"), maximum);
    }
    if (family == "PointCheck")
        check.each_bound(check.field(node, "tolerance"), 0, false);
    if (family == "InitialPoint" || family == "MaterialPoint") {
        const auto target = check.field(node, "target"),
                   coordinates = check.field(node, "coordinates");
        const auto kind = check.tag(target);
        if (kind == "beam") {
            for (std::uint64_t k = 1; k < 3; ++k)
                check.require(
                    schema_compare_integer(c, check.number(check.at(coordinates, k)), 0) == 0);
        }
        if (family == "InitialPoint" && kind == "shell") {
            for (std::uint64_t k = 0; k < 2; ++k) {
                check.numeric_bound(check.at(coordinates, k), 0, false);
                check.require(
                    schema_compare_integer(c, check.number(check.at(coordinates, k)), 1) <= 0);
            }
            check.numeric_bound(check.at(coordinates, 2), -1, false);
            check.require(schema_compare_integer(c, check.number(check.at(coordinates, 2)), 1) <=
                          0);
        }
    }
    if (family == "Observation") {
        if (check.tag(node) == "stress")
            check.require(check.tag(node, "measure") != "green" ||
                          check.tag(node, "contribution") == "total");
        for (auto key : {"surfaces", "targets", "points"}) {
            const auto value = present(key);
            if (value) {
                check.nonempty(*value);
                check.unique(*value);
            }
        }
    }
    if (family == "Comparison") {
        const auto window = check.field(node, "window");
        check.interval(window, check.tag(node, "statistic") != "instantaneous");
        if (check.tag(node, "statistic") == "instantaneous")
            check.require(schema_compare(c, check.number(check.at(window, 0)),
                                         check.number(check.at(window, 1))) == 0);
    }
    if (family == "ComparisonReference" && check.tag(node) == "file") {
        const auto columns = check.field(node, "columns");
        check.nonempty(columns);
        check.length(check.field(node, "uncertainty"), check.size(columns));
        check.each_bound(check.field(node, "uncertainty"), 0, false);
    }
    if (family == "ObservationBudget") {
        const auto absolute = check.field(node, "absolute"),
                   relative = check.field(node, "relative"), scale = check.field(node, "scale");
        check.length(relative, check.size(absolute));
        check.length(scale, check.size(absolute));
        check.each_bound(absolute, 0, false);
        check.each_bound(relative, 0, false);
        check.each_bound(scale, 0, true);
    }
    if (family == "ErrorGroup") {
        const auto status = check.tag(node, "status");
        if (status == "not_applicable")
            check.text_not_empty(check.field(node, "reason"));
        if (status == "required") {
            const auto covers = check.field(node, "covers");
            check.nonempty(covers);
            check.unique(covers);
            check.each_bound(check.field(node, "absolute"), 0, false);
            check.each_bound(check.field(node, "relative"), 0, false);
        }
    }
    if (family == "PatchPolicy" && check.tag(node) == "active") {
        check.require(check.count(check.field(node, "maximum_degree")) <= 4);
        const auto partition = check.count(check.field(node, "maximum_partition"));
        check.require(partition == 1 || partition == 2 || partition == 4);
        positive("geometry_error");
        positive("trace_error");
    }
    if (family == "StructuralQuadrature" && check.tag(node) == "active") {
        check.require(check.count(check.field(node, "independent_surface_degree")) >
                      check.count(check.field(node, "surface_degree")));
        check.require(check.count(check.field(node, "independent_thickness_points")) >
                      check.count(check.field(node, "thickness_points")));
    }
    if (family == "StructuralQuadrature" || family == "AeroQuadrature" ||
        family == "LinearAccuracy" || family == "PrecisionPolicy")
        for (auto key :
             {"absolute", "relative", "backward", "projection_absolute", "projection_relative"})
            nonnegative(key);
    if (family == "Scale") {
        if (check.tag(node) != "physical_layout_policy")
            check.each_bound(check.field(node, "values"), 0, true);
        if (check.tag(node) == "explicit_layout") {
            const auto dimensions = check.field(node, "dimensions");
            check.length(dimensions, check.size(check.field(node, "values")));
            for (std::uint64_t j = 0; j < check.size(dimensions); ++j)
                for (std::uint64_t k = 0; k < 4; ++k) {
                    const auto power = check.s.nodes
                                           .values()[static_cast<std::size_t>(
                                               check.at(check.at(dimensions, j), k))]
                                           .signed_value;
                    check.require(power >= std::numeric_limits<std::int16_t>::min() &&
                                  power <= std::numeric_limits<std::int16_t>::max());
                }
        }
        const auto targets = present("targets");
        if (targets)
            check.unique(*targets);
    }
    if (family == "ComponentScale3" && check.tag(node) == "values")
        check.each_bound(check.field(node, "values"), 0, true);
    if (family == "ConstraintScalePolicy" && check.tag(node) == "physical") {
        positive("multiplier_work_scale");
        check.each_bound(check.field(node, "reaction_force"), 0, true);
        check.each_bound(check.field(node, "reaction_moment"), 0, true);
    }
    if (family == "LinearPolicy") {
        const auto refinements = present("maximum_refinements");
        if (refinements)
            check.require(check.count(*refinements) <= 5);
    }
    if (family == "NewtonPolicy") {
        for (auto key : {"absolute", "relative", "correction"})
            nonnegative(key);
        positive("minimum_step");
        check.require(
            schema_compare_integer(c, check.number(check.field(node, "minimum_step")), 1) <= 0);
    }
    if (family == "EigenPolicy") {
        positive("eigenvalue_scale");
        positive("backend_tolerance");
    }
    if (family == "ArnoldiSeed") {
        const auto nonzero = [&](std::uint64_t z) {
            return schema_compare_integer(c, check.number(check.field(z, "real")), 0) != 0 ||
                   schema_compare_integer(c, check.number(check.field(z, "imaginary")), 0) != 0;
        };
        if (check.tag(node) == "uniform_nonzero")
            check.require(nonzero(check.field(node, "value")));
        else {
            const auto values = check.field(node, "values");
            bool any = false;
            for (std::uint64_t j = 0; j < check.size(values); ++j)
                any = nonzero(check.at(values, j)) || any;
            check.require(any);
        }
    }
    if (family == "EigenRequirements") {
        for (auto key : {"residual", "normalization", "orthogonality", "association_absolute",
                         "association_relative"})
            nonnegative(key);
        check.require(
            schema_compare_integer(c, check.number(check.field(node, "orthogonality")), 1) < 0);
    }
    if (family == "Backward")
        nonnegative("budget");
    if (family == "EigenCluster") {
        const auto members = check.field(node, "members");
        check.nonempty(members);
        check.unique(members);
        check.require(check.count(check.field(node, "expected_dimension")) == check.size(members));
    }
    if (family == "Resources")
        positive("maximum_wall_seconds");
    if (family == "Output" || family == "CheckpointOutput") {
        const auto directory = present("directory");
        if (directory)
            check.text_not_empty(*directory, true);
    }
    if (family == "OutputCadence")
        check.unique(check.field(node, "observations"));
    if (family == "ExpressionDomain") {
        const auto kind = check.tag(node);
        if (kind == "reference") {
            const auto lower = check.field(node, "lower"), upper = check.field(node, "upper");
            for (std::uint64_t j = 0; j < 3; ++j)
                check.require(schema_compare(c, check.number(check.at(lower, j)),
                                             check.number(check.at(upper, j))) <= 0);
        } else {
            const auto lower = check.field(node, kind == "time" ? "lower" : "time_lower"),
                       upper = check.field(node, kind == "time" ? "upper" : "time_upper");
            check.require(schema_compare(c, check.number(lower), check.number(upper)) < 0);
            const auto events = check.field(node, "events");
            check.order(events, true);
            for (std::uint64_t j = 0; j < check.size(events); ++j) {
                check.require(
                    schema_compare(c, check.number(check.at(events, j)), check.number(lower)) >= 0);
                check.require(
                    schema_compare(c, check.number(check.at(events, j)), check.number(upper)) <= 0);
            }
            if (kind == "current_position_time") {
                const auto a = check.field(node, "position_lower"),
                           b = check.field(node, "position_upper");
                for (std::uint64_t j = 0; j < 3; ++j)
                    check.require(schema_compare(c, check.number(check.at(a, j)),
                                                 check.number(check.at(b, j))) < 0);
            }
        }
    }
}
} // namespace
namespace {
struct ExactParts {
    int sign;
    std::string_view digits;
    std::int64_t exponent;
};
ExactParts owned_parts(SchemaContext& c, const DecimalRange& value) {
    return {value.sign, digits(c, value), value.exponent};
}
void plan_exact_scratch(SemanticCheck& check) {
    auto& c = check.c;
    c.exact.digits = c.allocate<char>(c.multiply(c.exact.stride, 6), true);
}
std::span<char> exact_buffer(SchemaContext& c, std::size_t slot) {
    return c.exact.digits.values().subspan(slot * static_cast<std::size_t>(c.exact.stride),
                                           static_cast<std::size_t>(c.exact.stride));
}
void zero_buffer(SchemaContext& c, std::size_t slot) {
    auto buffer = exact_buffer(c, slot);
    c.work(buffer.size());
    std::ranges::fill(buffer, '0');
}
std::string_view normalized_buffer(SchemaContext& c, std::size_t slot) {
    auto buffer = exact_buffer(c, slot);
    std::size_t begin = 0;
    while (begin + 1 < buffer.size() && buffer[begin] == '0') {
        c.work();
        ++begin;
    }
    return {buffer.data() + begin, buffer.size() - begin};
}
ExactParts product(SchemaContext& c, ExactParts a, ExactParts b, std::size_t slot) {
    const auto exponent = signed_add(c, a.exponent, b.exponent);
    const auto width = c.add(a.digits.size(), b.digits.size());
    if (width >= c.exact.stride)
        c.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                 SchemaQuota::ExactDigits, width);
    zero_buffer(c, slot);
    auto out = exact_buffer(c, slot);
    const auto offset = out.size() - static_cast<std::size_t>(width);
    for (std::size_t i = a.digits.size(); i > 0; --i) {
        int carry = 0;
        for (std::size_t j = b.digits.size(); j > 0; --j) {
            c.work();
            const auto k = offset + i + j - 1;
            const int value =
                (a.digits[i - 1] - '0') * (b.digits[j - 1] - '0') + (out[k] - '0') + carry;
            out[k] = static_cast<char>('0' + value % 10);
            carry = value / 10;
        }
        out[offset + i - 1] = static_cast<char>('0' + carry);
    }
    return {a.sign * b.sign, normalized_buffer(c, slot), exponent};
}
std::pair<ExactParts, ExactParts> rational_parts(SemanticCheck& check, std::uint64_t node) {
    const auto& typed = check.s.nodes.values()[static_cast<std::size_t>(node)];
    if (node_kind(check.d, node) == InputNodeKind::Object)
        return {owned_parts(check.c, typed.rational_numerator),
                owned_parts(check.c, typed.rational_denominator)};
    return {owned_parts(check.c, typed.number), ExactParts{1, "1", 0}};
}
int compare_exact(SemanticCheck& check, std::uint64_t left, std::uint64_t right) {
    auto& c = check.c;
    const auto a = rational_parts(check, left), b = rational_parts(check, right);
    if (a.first.sign != b.first.sign)
        return a.first.sign > b.first.sign ? 1 : -1;
    if (a.first.sign == 0)
        return 0;
    const auto l = product(c, a.first, b.second, 0), r = product(c, b.first, a.second, 1);
    return compare_parts(c, l.sign, l.digits, l.exponent, r.sign, r.digits, r.exponent);
}
int compare_exact_integer(SemanticCheck& check, std::uint64_t node, std::int64_t integer) {
    auto& c = check.c;
    const auto parts = rational_parts(check, node);
    std::array<char, 32> buffer{};
    const auto conversion =
        std::to_chars(buffer.data(), buffer.data() + buffer.size(), magnitude(integer));
    const auto right = product(
        c,
        ExactParts{integer == 0 ? 0 : (integer < 0 ? -1 : 1),
                   {buffer.data(), static_cast<std::size_t>(conversion.ptr - buffer.data())},
                   0},
        parts.second, 1);
    return compare_parts(c, parts.first.sign, parts.first.digits, parts.first.exponent, right.sign,
                         right.digits, right.exponent);
}
int compare_buffers(SchemaContext& c, std::size_t left, std::size_t right) {
    const auto a = exact_buffer(c, left), b = exact_buffer(c, right);
    for (std::size_t j = 0; j < a.size(); ++j) {
        c.work();
        if (a[j] != b[j])
            return a[j] > b[j] ? 1 : -1;
    }
    return 0;
}
void load_scaled(SemanticCheck& check, std::uint64_t node, std::int64_t base, std::size_t slot) {
    auto& c = check.c;
    zero_buffer(c, slot);
    const auto& value = check.number(node);
    if (value.sign == 0)
        return;
    const auto difference = signed_subtract(c, value.exponent, base);
    check.require(difference >= 0);
    const auto required = c.add(value.count, static_cast<std::uint64_t>(difference));
    if (required >= c.exact.stride)
        c.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                 SchemaQuota::ExactDigits, required);
    const auto input = digits(c, value);
    auto out = exact_buffer(c, slot);
    const auto begin = out.size() - static_cast<std::size_t>(required);
    for (std::size_t j = 0; j < input.size(); ++j) {
        c.work();
        out[begin + j] = input[j];
    }
}
void add_buffers(SchemaContext& c, std::size_t target, std::size_t value) {
    auto a = exact_buffer(c, target);
    const auto b = exact_buffer(c, value);
    int carry = 0;
    for (std::size_t j = a.size(); j > 0; --j) {
        c.work();
        const auto sum = a[j - 1] - '0' + b[j - 1] - '0' + carry;
        a[j - 1] = static_cast<char>('0' + sum % 10);
        carry = sum / 10;
    }
    if (carry)
        c.reject(SchemaReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                 SchemaQuota::ExactDigits, c.add(c.exact.stride, 1));
}
void subtract_buffers(SchemaContext& c, std::size_t target, std::size_t value) {
    auto a = exact_buffer(c, target);
    const auto b = exact_buffer(c, value);
    int borrow = 0;
    for (std::size_t j = a.size(); j > 0; --j) {
        c.work();
        int difference = a[j - 1] - b[j - 1] - borrow;
        if (difference < 0) {
            difference += 10;
            borrow = 1;
        } else
            borrow = 0;
        a[j - 1] = static_cast<char>('0' + difference);
    }
    if (borrow)
        c.reject(SchemaReason::InvalidValue);
}
std::size_t positive_difference(SemanticCheck& check, std::uint64_t lower, std::uint64_t upper,
                                std::int64_t base) {
    load_scaled(check, lower, base, 0);
    load_scaled(check, upper, base, 1);
    const auto ls = check.number(lower).sign, us = check.number(upper).sign;
    if (us >= 0 && ls <= 0) {
        add_buffers(check.c, 1, 0);
        return 1;
    }
    if (us >= 0) {
        subtract_buffers(check.c, 1, 0);
        return 1;
    }
    subtract_buffers(check.c, 0, 1);
    return 0;
}
void time_divisibility(SemanticCheck& check, std::uint64_t node) {
    auto& c = check.c;
    const auto start = check.field(node, "start"), end = check.field(node, "end"),
               step = check.field(node, "step");
    const auto base = std::min(
        {check.number(start).exponent, check.number(end).exponent, check.number(step).exponent});
    const auto numerator = positive_difference(check, start, end, base);
    load_scaled(check, step, base, 2);
    zero_buffer(c, 3);
    const auto source = normalized_buffer(c, numerator);
    for (const char digit : source) {
        auto remainder = exact_buffer(c, 3);
        check.require(remainder[0] == '0');
        for (std::size_t j = 0; j + 1 < remainder.size(); ++j) {
            c.work();
            remainder[j] = remainder[j + 1];
        }
        remainder.back() = digit;
        while (compare_buffers(c, 3, 2) >= 0)
            subtract_buffers(c, 3, 2);
    }
    for (const char digit : exact_buffer(c, 3)) {
        c.work();
        check.require(digit == '0');
    }
}
void layer_checks(SemanticCheck& check, std::uint64_t shell) {
    const auto layers = check.field(shell, "layers");
    check.nonempty(layers);
    std::optional<std::uint64_t> previous;
    for (std::uint64_t j = 0; j < check.size(layers); ++j) {
        const auto layer = check.at(layers, j);
        const auto interval = check.field(layer, "interval");
        const auto lower = check.at(interval, 0), upper = check.at(interval, 1);
        check.require(compare_exact(check, lower, upper) < 0);
        check.require(compare_exact_integer(check, lower, -1) >= 0);
        check.require(compare_exact_integer(check, upper, 1) <= 0);
        if (previous)
            check.require(compare_exact(check, *previous, lower) == 0);
        else
            check.require(compare_exact_integer(check, lower, -1) == 0);
        previous = upper;
    }
    if (!previous)
        check.c.reject(SchemaReason::InvalidValue);
    check.require(compare_exact_integer(check, *previous, 1) == 0);
}
void periodic_checks(SemanticCheck& check, std::uint64_t series) {
    const auto outside = check.field(series, "outside");
    if (check.tag(outside) != "periodic")
        return;
    const auto times = check.field(series, "times");
    const auto lower = check.at(times, 0), upper = check.at(times, check.size(times) - 1),
               period = check.field(outside, "period");
    const auto base = std::min({check.number(lower).exponent, check.number(upper).exponent,
                                check.number(period).exponent});
    const auto difference = positive_difference(check, lower, upper, base);
    load_scaled(check, period, base, 2);
    check.require(compare_buffers(check.c, difference, 2) == 0);
}
} // namespace
namespace {
void budget_sum(SemanticCheck& check, std::uint64_t groups, std::uint64_t parent,
                std::string_view field, std::uint64_t component) {
    auto& c = check.c;
    const auto parent_value = check.at(check.field(parent, field), component);
    auto base = check.number(parent_value).exponent;
    for (std::uint64_t j = 0; j < check.size(groups); ++j) {
        const auto group = check.at(groups, j);
        if (check.tag(group, "status") == "required")
            base = std::min(base,
                            check.number(check.at(check.field(group, field), component)).exponent);
    }
    zero_buffer(c, 0);
    for (std::uint64_t j = 0; j < check.size(groups); ++j) {
        const auto group = check.at(groups, j);
        if (check.tag(group, "status") != "required")
            continue;
        load_scaled(check, check.at(check.field(group, field), component), base, 1);
        add_buffers(c, 0, 1);
    }
    load_scaled(check, parent_value, base, 2);
    check.require(compare_buffers(c, 0, 2) <= 0);
}
void error_groups(SemanticCheck& check, std::uint64_t groups, std::optional<std::uint64_t> budget,
                  std::uint64_t components) {
    check.length(groups, 12);
    std::array<std::optional<std::uint64_t>, 12> indexed{};
    for (std::uint64_t j = 0; j < 12; ++j) {
        const auto group = check.at(groups, j);
        const auto kind = check.count(check.field(group, "kind"));
        check.require(kind < 12);
        check.require(!indexed[static_cast<std::size_t>(kind)]);
        indexed[static_cast<std::size_t>(kind)] = group;
    }
    for (std::size_t j = 0; j < 12; ++j) {
        const auto entry = indexed[j];
        if (!entry)
            check.c.reject(SchemaReason::InvalidValue);
        const auto group = *entry;
        const auto status = check.tag(group, "status");
        if (status == "covered_by") {
            const auto target = check.count(check.field(group, "representative"));
            check.require(target < 12 && target != j);
            const auto other = indexed[static_cast<std::size_t>(target)];
            if (!other)
                check.c.reject(SchemaReason::InvalidValue);
            check.require(check.tag(*other, "status") == "required");
        }
        if (status == "required") {
            check.length(check.field(group, "absolute"), components);
            check.length(check.field(group, "relative"), components);
            const auto covers = check.field(group, "covers");
            std::array<bool, 12> declared{};
            for (std::uint64_t k = 0; k < check.size(covers); ++k) {
                const auto kind = check.count(check.at(covers, k));
                check.require(kind < 12 && !declared[static_cast<std::size_t>(kind)]);
                declared[static_cast<std::size_t>(kind)] = true;
            }
            for (std::size_t k = 0; k < 12; ++k) {
                const auto other = indexed[k];
                if (!other)
                    check.c.reject(SchemaReason::InvalidValue);
                const bool belongs =
                    k == j || (check.tag(*other, "status") == "covered_by" &&
                               check.count(check.field(*other, "representative")) == j);
                check.require(declared[k] == belongs);
            }
        }
    }
    if (budget)
        for (std::uint64_t component = 0; component < components; ++component)
            for (auto field : {"absolute", "relative"})
                budget_sum(check, groups, *budget, field, component);
}
std::uint64_t observation_components(SemanticCheck& check, std::uint64_t observation) {
    const auto kind = check.tag(observation);
    if (kind == "reaction")
        return 6;
    if (kind == "mode")
        return 1;
    if (kind == "stress")
        return check.tag(observation, "measure") == "first_pk" ? 9 : 6;
    return 3;
}
void observation_checks(SemanticCheck& check, std::uint64_t observation) {
    auto& c = check.c;
    c.current = observation;
    const auto n = observation_components(check, observation);
    const auto budget = check.field(observation, "budget");
    for (auto key : {"absolute", "relative", "scale"})
        check.length(check.field(budget, key), n);
    error_groups(check, check.field(budget, "groups"), budget, n);
    const auto kind = check.tag(observation);
    const auto basis = check.tag(check.field(observation, "basis"));
    if (basis == "reference_material" || basis == "current_material")
        check.require(kind == "stress");
    if (kind == "stress") {
        const auto point = check.field(observation, "point");
        const auto target = check.field(point, "target");
        check.require(check.tag(target) == "shell");
        const auto shell = check.referent(check.field(target, "id"));
        const auto layer = check.referent(check.field(observation, "layer"));
        const auto layers = check.field(shell, "layers");
        bool owns = false;
        for (std::uint64_t j = 0; j < check.size(layers); ++j)
            owns = check.at(layers, j) == layer || owns;
        check.require(owns);
    }
}
bool scale_targets(SemanticCheck& check, std::uint64_t scale, std::string_view kind,
                   std::uint64_t id) {
    const auto targets = check.field(scale, "targets");
    if (check.size(targets) == 0)
        return true;
    for (std::uint64_t j = 0; j < check.size(targets); ++j) {
        const auto target = check.at(targets, j);
        if (check.tag(target) == kind && check.count(check.field(target, "id")) == id)
            return true;
    }
    return false;
}
void scale_checks(SemanticCheck& check, std::uint64_t scales) {
    std::optional<std::uint64_t> physical;
    bool explicit_mixed = false;
    for (std::uint64_t j = 0; j < check.size(scales); ++j) {
        const auto scale = check.at(scales, j);
        const auto kind = check.tag(scale), field = check.tag(scale, "field");
        if (kind == "physical_layout_policy") {
            check.require(!physical);
            physical = scale;
            continue;
        }
        const bool mixed =
            field == "constraint" || field == "newton_correction" || field == "arclength_state";
        explicit_mixed = mixed || explicit_mixed;
        if (kind == "uniform_components") {
            check.require(!mixed);
            const bool scalar =
                field == "circulation" || field == "wake_strength" || field == "aero_residual";
            check.length(check.field(scale, "values"), scalar ? 1 : 3);
        }
        const auto targets = check.field(scale, "targets");
        for (std::uint64_t k = 0; k < check.size(targets); ++k) {
            const auto target = check.at(targets, k);
            const auto target_kind = check.tag(target);
            if (field == "shell_position")
                check.require(target_kind == "shell");
            if (field == "beam_translation" || field == "beam_rotation")
                check.require(target_kind == "beam");
            if (field == "rigid_translation" || field == "rigid_rotation")
                check.require(target_kind == "rigid");
            if (field == "structural_angular_velocity" || field == "moment_balance")
                check.require(target_kind != "shell");
            if (field == "circulation" || field == "wake_position" || field == "wake_strength" ||
                field == "aero_residual" || field == "transport_residual")
                check.c.reject(SchemaReason::InvalidValue);
        }
        for (std::uint64_t k = 0; k < j; ++k) {
            const auto prior = check.at(scales, k);
            if (check.tag(prior) == "physical_layout_policy" || check.tag(prior, "field") != field)
                continue;
            const auto previous_targets = check.field(prior, "targets");
            if (check.size(targets) == 0 || check.size(previous_targets) == 0)
                check.c.reject(SchemaReason::InvalidValue);
            for (std::uint64_t t = 0; t < check.size(targets); ++t) {
                const auto target = check.at(targets, t);
                check.require(!scale_targets(check, prior, check.tag(target),
                                             check.count(check.field(target, "id"))));
            }
        }
    }
    check.require(!physical || !explicit_mixed);
    const auto analysis = check.field(0, "analysis");
    const auto analysis_kind = check.tag(analysis);
    const auto require_scale = [&](std::string_view field, std::string_view target_kind,
                                   std::uint64_t target_id) {
        bool found = false;
        for (std::uint64_t j = 0; j < check.size(scales); ++j) {
            const auto candidate = check.at(scales, j);
            if (check.tag(candidate) == "physical_layout_policy")
                continue;
            if (check.tag(candidate, "field") == field &&
                scale_targets(check, candidate, target_kind, target_id))
                found = true;
        }
        check.require(found);
    };
    for (auto collection : {"shells", "beams", "masses"}) {
        const auto entities = check.model_list(collection);
        for (std::uint64_t j = 0; j < check.size(entities); ++j) {
            const auto entity = check.at(entities, j);
            if (collection == std::string_view{"masses"} && check.tag(entity) != "rigid")
                continue;
            const auto id = check.count(check.field(entity, "id"));
            const auto kind = collection == std::string_view{"shells"}
                                  ? "shell"
                                  : (collection == std::string_view{"beams"} ? "beam" : "rigid");
            require_scale(
                kind == std::string_view{"shell"}
                    ? "shell_position"
                    : (kind == std::string_view{"beam"} ? "beam_translation" : "rigid_translation"),
                kind, id);
            require_scale("force_balance", kind, id);
            if (kind != std::string_view{"shell"}) {
                require_scale(kind == std::string_view{"beam"} ? "beam_rotation" : "rigid_rotation",
                              kind, id);
                require_scale("moment_balance", kind, id);
            }
            if (analysis_kind == "transient" || analysis_kind == "flutter") {
                require_scale("structural_velocity", kind, id);
                if (kind != std::string_view{"shell"})
                    require_scale("structural_angular_velocity", kind, id);
            }
        }
    }
    if (check.size(check.model_list("aero_surfaces")))
        for (auto field : {"circulation", "wake_position", "wake_strength", "aero_residual",
                           "transport_residual"})
            require_scale(field, "", 0);
    bool position = false, orientation = false;
    const auto joints = check.model_list("joints");
    for (std::uint64_t j = 0; j < check.size(joints); ++j) {
        position = true;
        if (check.tag(check.at(joints, j)) != "adhesive")
            orientation = true;
    }
    const auto supports = check.model_list("supports");
    for (std::uint64_t j = 0; j < check.size(supports); ++j) {
        const auto support = check.at(supports, j);
        const auto kind = check.tag(support);
        if (kind == "trace_constraint")
            position = true;
        if (kind == "constraint") {
            const auto components = check.field(support, "components");
            for (std::uint64_t k = 0; k < check.size(components); ++k) {
                const auto component = check.count(check.at(components, k));
                if (component < 3)
                    position = true;
                else
                    orientation = true;
            }
        }
    }
    if (physical) {
        const auto constraint = check.field(*physical, "constraint");
        const auto kind = check.tag(constraint);
        if (position || orientation) {
            check.require(kind == "physical");
            check.require(check.tag(check.field(constraint, "position")) ==
                          (position ? "values" : "none"));
            check.require(check.tag(check.field(constraint, "orientation")) ==
                          (orientation ? "values" : "none"));
        } else
            check.require(kind == "none");
    } else {
        bool correction = false, constraint = !(position || orientation);
        for (std::uint64_t j = 0; j < check.size(scales); ++j) {
            const auto scale = check.at(scales, j);
            if (check.tag(scale) != "explicit_layout")
                continue;
            const auto field = check.tag(scale, "field");
            correction = correction || field == "newton_correction";
            constraint = constraint || field == "constraint";
        }
        check.require(correction && constraint);
    }
}
void root_checks(SemanticCheck& check) {
    auto& c = check.c;
    c.current = 0;
    const auto analysis = check.field(0, "analysis"), numerics = check.field(0, "numerics");
    const auto kind = check.tag(analysis), physics = check.tag(analysis, "physics");
    const auto surfaces = check.model_list("aero_surfaces");
    const bool aero = check.size(surfaces) > 0;
    const bool structure = check.size(check.model_list("shells")) > 0 ||
                           check.size(check.model_list("beams")) > 0 ||
                           check.size(check.model_list("masses")) > 0;
    if (physics == "structural") {
        check.require(!aero);
        check.require(structure);
    }
    if (physics == "aerodynamic" || physics == "coupled")
        check.require(aero);
    if (physics == "coupled") {
        check.require(structure);
        bool attached = false;
        for (std::uint64_t j = 0; j < check.size(surfaces); ++j)
            attached =
                check.tag(check.field(check.at(surfaces, j), "attachment")) == "structural" ||
                attached;
        check.require(attached);
    }
    check.require(check.tag(check.field(0, "flow")) == (aero ? "known" : "none"));
    check.require(check.tag(check.field(analysis, "wake")) == (aero ? "active" : "none"));
    check.require(check.tag(check.field(numerics, "wake_rows")) == (aero ? "limited" : "none"));
    const auto patch = check.field(numerics, "patch");
    check.require(check.tag(patch) == (check.size(check.model_list("shells")) ? "active" : "none"));
    const auto quadrature = check.field(numerics, "quadrature");
    const bool structural_integration =
        check.size(check.model_list("shells")) > 0 || check.size(check.model_list("beams")) > 0;
    check.require(check.tag(check.field(quadrature, "structure")) ==
                  (structural_integration ? "active" : "none"));
    check.require(check.tag(check.field(quadrature, "aero")) ==
                  (aero ? "wp1_induction_v1" : "none"));
    const auto initial = check.maybe(0, "initial"), restart = check.maybe(0, "restart");
    if (kind == "transient")
        check.require(initial.has_value() != restart.has_value());
    else
        check.require(!initial);
    if (kind == "static")
        check.require(check.tag(check.field(analysis, "speed")) ==
                      (physics == "structural" ? "none" : "value"));
    if (kind == "flutter") {
        const auto eigen = check.field(numerics, "eigen");
        check.require(check.tag(eigen) != "none");
        // Flutter growth comes from the coupled step map, not the structural K/M pencil.
        c.current = eigen;
        check.require(check.tag(eigen, "target") == "step_map");
        const auto requirements = check.field(eigen, "requirements");
        check.require(std::get<bool>(
            check.d.nodes
                .values()[static_cast<std::size_t>(check.field(requirements, "physical_coverage"))]
                .value));
        const auto coverage = check.field(check.field(analysis, "search"), "coverage");
        check.require(std::get<bool>(
            check.d.nodes
                .values()[static_cast<std::size_t>(check.field(coverage, "require_physical"))]
                .value));
    }
    const auto observations = check.field(0, "observations");
    for (std::uint64_t j = 0; j < check.size(observations); ++j)
        observation_checks(check, check.at(observations, j));
    error_groups(check, check.field(numerics, "error_groups"), std::nullopt, 1);
    scale_checks(check, check.field(numerics, "scales"));
    if (kind == "transient") {
        const auto intervals = check.field(analysis, "intervals");
        check.nonempty(intervals);
        std::optional<std::uint64_t> previous;
        for (std::uint64_t j = 0; j < check.size(intervals); ++j) {
            const auto interval = check.at(intervals, j);
            if (previous)
                check.require(schema_compare(c, check.number(check.field(*previous, "end")),
                                             check.number(check.field(interval, "start"))) == 0);
            previous = interval;
        }
    }
    if (initial && check.tag(*initial) == "equilibrium")
        check.require(check.count(check.field(*initial, "branch")) ==
                      check.count(check.field(check.field(analysis, "branch"), "id")));
}
} // namespace
namespace {
std::optional<std::uint64_t> ancestor_record(SemanticCheck& check, std::uint64_t node,
                                             std::string_view family) {
    auto ancestor = check.s.nodes.values()[static_cast<std::size_t>(node)].parent;
    while (ancestor) {
        check.c.work();
        const auto& typed = check.s.nodes.values()[static_cast<std::size_t>(*ancestor)];
        if (typed.assigned && node_kind(check.d, *ancestor) == InputNodeKind::Object &&
            check.record(*ancestor).family == family)
            return ancestor;
        ancestor = typed.parent;
    }
    return std::nullopt;
}
void boundary_index(SemanticCheck& check, std::uint64_t geometry, std::uint64_t index) {
    const auto boundaries = check.maybe(geometry, "boundaries");
    if (boundaries)
        check.require(index < check.size(*boundaries));
}
void cross_record_checks(SemanticCheck& check, std::uint64_t node) {
    auto& c = check.c;
    c.current = node;
    const auto family = check.record(node).family;
    if (family == "AeroSurface") {
        const auto geometry = check.referent(check.field(node, "geometry"));
        boundary_index(check, geometry, check.count(check.field(node, "leading")));
        boundary_index(check, geometry, check.count(check.field(node, "trailing")));
    }
    if (family == "Continuity") {
        const auto own = ancestor_record(check, node, "Geometry");
        if (!own)
            c.reject(SchemaReason::InvalidValue);
        const auto other = check.referent(check.field(node, "other"));
        check.require(*own != other);
        boundary_index(check, *own, check.count(check.field(node, "own_boundary")));
        boundary_index(check, other, check.count(check.field(node, "other_boundary")));
    }
    if (family == "AeroSharing") {
        const auto own = ancestor_record(check, node, "AeroSurface");
        if (!own)
            c.reject(SchemaReason::InvalidValue);
        const auto other = check.referent(check.field(node, "other"));
        boundary_index(check, check.referent(check.field(*own, "geometry")),
                       check.count(check.field(node, "own_boundary")));
        boundary_index(check, check.referent(check.field(other, "geometry")),
                       check.count(check.field(node, "other_boundary")));
    }
    if (family == "Region") {
        const auto geometry = check.referent(check.field(node, "geometry"));
        const auto kind = check.tag(node);
        const auto boundary = check.count(check.field(node, "boundary"));
        if (kind == "face" || kind == "edge")
            boundary_index(check, geometry, boundary);
        const auto target = check.field(node, "target");
        if (check.tag(target) == "shell") {
            const auto shell = check.referent(check.field(target, "id"));
            check.require(check.count(check.field(shell, "geometry")) ==
                          check.count(check.field(node, "geometry")));
        }
    }
    if (family == "InitialState" || (family == "Perturbation" && check.tag(node) == "physical")) {
        check.complete_states(check.field(node, "shells"), check.model_list("shells"), "shell");
        check.complete_states(check.field(node, "beams"), check.model_list("beams"), "beam");
        check.complete_states(check.field(node, "rigids"), check.model_list("masses"), "rigid",
                              true);
    }
    if (family == "Perturbation" && check.tag(node) == "physical")
        for (auto field : {"displacements", "velocities"}) {
            const auto points = check.field(node, field);
            for (std::uint64_t j = 0; j < check.size(points); ++j)
                for (std::uint64_t k = 0; k < j; ++k)
                    check.require(!check.same_value(check.field(check.at(points, j), "point"),
                                                    check.field(check.at(points, k), "point")));
        }
    if (family == "InitialWake")
        check.complete_states(check.field(node, "surfaces"), check.model_list("aero_surfaces"),
                              "surface");
    if (family == "SurfaceWake") {
        const auto surface = check.referent(check.field(node, "surface"));
        const auto span = check.field(surface, "span");
        check.require(check.size(span) >= 2);
        const auto n = check.size(span) - 1;
        const auto rows = check.count(check.field(node, "rows"));
        const auto points = c.multiply(rows, c.add(n, 1));
        check.length(check.field(node, "points"), points);
        check.length(check.field(node, "coefficients"), c.multiply(rows, n));
        check.length(check.field(node, "row_ids"), rows);
        check.length(check.field(node, "birth_times"), rows);
        check.length(check.field(node, "point_ids"), points);
        check.length(check.field(node, "band_ids"), n);
        const auto analysis = check.field(0, "analysis");
        if (check.tag(analysis) == "transient") {
            const auto intervals = check.field(analysis, "intervals");
            check.nonempty(intervals);
            const auto start = check.field(check.at(intervals, 0), "start");
            const auto births = check.field(node, "birth_times");
            for (std::uint64_t j = 0; j < check.size(births); ++j)
                check.require(
                    schema_compare(c, check.number(check.at(births, j)), check.number(start)) <= 0);
        }
    }
    if (family == "ShellPhysicalState") {
        const auto shell = check.referent(check.field(node, "shell"));
        const auto geometry = check.referent(check.field(shell, "geometry"));
        if (check.tag(geometry) == "pg") {
            const auto n = check.count(check.field(geometry, "n"));
            check.require(n == 1 || n == 2 || n == 4);
            const auto maximum = c.multiply(n, 2);
            for (auto key : {"displacement", "velocity"}) {
                const auto value = check.field(node, key);
                if (check.tag(value) != "half_cells")
                    continue;
                const auto cells = check.field(value, "cells");
                for (std::uint64_t j = 0; j < check.size(cells); ++j) {
                    const auto cell = check.field(check.at(cells, j), "cell");
                    for (std::uint64_t k = 0; k < 2; ++k)
                        check.require(check.count(check.at(cell, k)) < maximum);
                }
            }
        }
    }
    if (family == "ArcLength") {
        const auto declared = check.field(node, "state_scales"),
                   actual = check.field(check.field(0, "numerics"), "scales");
        for (std::uint64_t j = 0; j < check.size(declared); ++j) {
            const auto item = check.at(declared, j);
            bool found = false;
            for (std::uint64_t k = 0; k < check.size(actual); ++k)
                found = check.same_value(item, check.at(actual, k)) || found;
            check.require(found);
        }
    }
    if (family == "Motion")
        for (auto key : {"translation", "rotation_vector"}) {
            const auto series = check.field(node, key);
            if (check.tag(series) != "pchip")
                continue;
            const auto times = check.field(series, "times"), events = check.field(node, "events");
            for (std::uint64_t j = 0; j < check.size(times); ++j) {
                bool found = false;
                for (std::uint64_t k = 0; k < check.size(events); ++k)
                    found = schema_compare(c, check.number(check.at(times, j)),
                                           check.number(check.at(events, k))) == 0 ||
                            found;
                check.require(found);
            }
        }
}
void mass_cycles(SemanticCheck& check, core::PlannedBuffer<std::uint8_t>& colors,
                 core::PlannedBuffer<std::uint64_t>& stack) {
    for (auto& color : colors.values()) {
        check.c.work();
        color = 0;
    }
    const auto masses = check.model_list("masses");
    for (std::uint64_t j = 0; j < check.size(masses); ++j) {
        const auto root = check.at(masses, j);
        if (colors.values()[static_cast<std::size_t>(root)] != 0)
            continue;
        std::uint64_t top = 1;
        stack.values()[0] = root;
        while (top) {
            check.c.work();
            const auto node = stack.values()[static_cast<std::size_t>(top - 1)];
            check.c.current = node;
            colors.values()[static_cast<std::size_t>(node)] = 1;
            const auto ledger = check.field(node, "ledger");
            bool descended = false;
            if (check.tag(ledger) == "replace") {
                const auto replaced = check.field(ledger, "replaced");
                for (std::uint64_t k = 0; k < check.size(replaced); ++k) {
                    const auto target = check.referent(check.at(replaced, k));
                    const auto color = colors.values()[static_cast<std::size_t>(target)];
                    check.require(color != 1, SchemaReason::ReferenceCycle);
                    if (color == 0) {
                        check.require(top < stack.size());
                        stack.values()[static_cast<std::size_t>(top++)] = target;
                        descended = true;
                        break;
                    }
                }
            }
            if (!descended) {
                colors.values()[static_cast<std::size_t>(node)] = 2;
                --top;
            }
        }
    }
}
void event_and_initial_checks(SemanticCheck& check) {
    auto& c = check.c;
    const auto analysis = check.field(0, "analysis");
    if (check.tag(analysis) == "transient") {
        const auto intervals = check.field(analysis, "intervals");
        check.nonempty(intervals);
        for (std::uint64_t j = 0; j < check.size(intervals); ++j) {
            const auto interval = check.at(intervals, j);
            check.require(check.tag(interval, "restart") != (j ? "initial" : ""));
            const auto start = check.field(interval, "start"), end = check.field(interval, "end");
            for (std::uint64_t node = 0; node < check.s.nodes.size(); ++node) {
                c.work();
                if (node_kind(check.d, node) != InputNodeKind::Object)
                    continue;
                const auto family = check.record(node).family;
                if (family != "Flow" && family != "Motion" && family != "ExpressionDomain")
                    continue;
                const auto events = check.maybe(node, "events");
                if (!events)
                    continue;
                for (std::uint64_t k = 0; k < check.size(*events); ++k) {
                    const auto event = check.at(*events, k);
                    check.require(
                        !(schema_compare(c, check.number(start), check.number(event)) < 0 &&
                          schema_compare(c, check.number(event), check.number(end)) < 0));
                }
            }
        }
        const auto initial = check.maybe(0, "initial");
        if (initial && check.tag(*initial) == "cold_start") {
            std::uint64_t circulation = 0;
            const auto surfaces = check.model_list("aero_surfaces");
            for (std::uint64_t j = 0; j < check.size(surfaces); ++j) {
                const auto surface = check.at(surfaces, j);
                const auto chord = check.size(check.field(surface, "chord")),
                           span = check.size(check.field(surface, "span"));
                check.require(chord >= 2 && span >= 2);
                circulation = c.add(circulation, c.multiply(chord - 1, span - 1));
            }
            check.length(check.field(*initial, "circulation"), circulation);
        }
    }
    if (check.tag(analysis) == "static") {
        std::optional<std::uint64_t> instant;
        const auto observations = check.field(0, "observations");
        for (std::uint64_t j = 0; j < check.size(observations); ++j) {
            const auto comparison = check.field(check.at(observations, j), "comparison"),
                       window = check.field(comparison, "window");
            const auto first = check.at(window, 0), last = check.at(window, 1);
            check.require(schema_compare(c, check.number(first), check.number(last)) == 0);
            if (instant)
                check.require(schema_compare(c, check.number(first), check.number(*instant)) == 0);
            else
                instant = first;
        }
    }
}
} // namespace
void schema_semantics(SchemaContext& context) {
    SemanticCheck check{context};
    plan_exact_scratch(check);
    for (std::uint64_t node = 0; node < check.s.nodes.size(); ++node) {
        context.work();
        if (check.s.nodes.values()[static_cast<std::size_t>(node)].assigned &&
            node_kind(check.d, node) == InputNodeKind::Object) {
            generic_record_checks(check, node);
            const auto family = check.record(node).family;
            if (family == "Shell")
                layer_checks(check, node);
            if (family == "TimeInterval")
                time_divisibility(check, node);
            if (family == "TimeSeries" && check.tag(node) == "pchip")
                periodic_checks(check, node);
        }
    }
    resolve_references(check);
    for (std::uint64_t node = 0; node < check.s.nodes.size(); ++node) {
        context.work();
        if (check.s.nodes.values()[static_cast<std::size_t>(node)].assigned &&
            node_kind(check.d, node) == InputNodeKind::Object)
            cross_record_checks(check, node);
    }
    auto visited = context.allocate<std::uint8_t>(check.s.nodes.size(), true);
    auto depths = context.allocate<std::uint64_t>(check.s.nodes.size(), true);
    if (context.whole_case || context.case_registry)
        mass_cycles(check, visited, depths);
    expression_checks(check, visited, depths);
    if (context.whole_case) {
        root_checks(check);
        event_and_initial_checks(check);
    } else if (context.case_registry) {
        const auto numerics = check.field(0, "numerics");
        error_groups(check, check.field(numerics, "error_groups"), std::nullopt, 1);
        scale_checks(check, check.field(numerics, "scales"));
        const auto observations = check.field(0, "observations");
        for (std::uint64_t j = 0; j < check.size(observations); ++j)
            observation_checks(check, check.at(observations, j));
        if (node_kind(check.d, context.payload_root) == InputNodeKind::Object &&
            check.record(context.payload_root).family == "Observation")
            observation_checks(check, context.payload_root);
    }
}
} // namespace uvlas::io::detail
