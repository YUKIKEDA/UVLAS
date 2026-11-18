"""Generate a test-only traversal of every declared, public typed case member.

This enumerates the public interface, never expected numerical answers. Scalar
encoding and recursion are independent of the product's private node storage.
"""
import argparse
import hashlib
from pathlib import Path
import re

PREFIX = r'''#pragma once
#include <bit>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <variant>
#include <vector>
#include <uvlas/io/case_schema.hpp>

namespace uvlas::test::schema_snapshot {
using Bytes = std::vector<std::uint8_t>;
inline void integer(Bytes& out, std::uint64_t value) {
    for (unsigned shift = 0; shift < 64; shift += 8)
        out.push_back(static_cast<std::uint8_t>(value >> shift));
}
inline void text(Bytes& out, std::string_view value) {
    integer(out, value.size());
    for (const char c : value) out.push_back(static_cast<std::uint8_t>(c));
}
inline void append(Bytes& out, const io::SchemaNumberView& value) {
    integer(out, static_cast<std::uint64_t>(value.sign));
    text(out, value.digits);
    integer(out, static_cast<std::uint64_t>(value.exponent10));
    integer(out, static_cast<std::uint64_t>(value.original_kind));
    for (const auto exponent : value.dimensions.exponents)
        integer(out, static_cast<std::uint64_t>(exponent));
}
'''
TEMPLATES = r'''
template <io::detail::FieldType Type> void append(Bytes&, const io::schema::Field<Type>&);
template <class... T> void append(Bytes& out, const std::variant<T...>& value) {
    integer(out, value.index());
    std::visit([&](const auto& selected) { append(out, selected); }, value);
}
template <class T> void append(Bytes& out, const std::optional<T>& value) {
    integer(out, value.has_value() ? 1U : 0U);
    if (value) append(out, *value);
}
template <io::detail::FieldType Type> void append(Bytes& out, const io::schema::Field<Type>& value) {
    // The static field type distinguishes units, ID namespaces, and selected schemas.
    integer(out, static_cast<std::uint64_t>(Type));
    constexpr auto kind = io::detail::field_kinds[static_cast<std::size_t>(Type)];
    if constexpr (kind == io::detail::FieldKind::Number) append(out, value.number());
    else if constexpr (kind == io::detail::FieldKind::Integer) {
        if constexpr (io::detail::field_is_signed[static_cast<std::size_t>(Type)])
            integer(out, static_cast<std::uint64_t>(value.signed_value()));
        else integer(out, value.unsigned_value());
    } else if constexpr (kind == io::detail::FieldKind::Identifier || kind == io::detail::FieldKind::Reference)
        integer(out, value.unsigned_value());
    else if constexpr (kind == io::detail::FieldKind::Text)
        text(out, value.text());
    else if constexpr (kind == io::detail::FieldKind::Enumeration) {
        text(out, value.text());
        integer(out, static_cast<std::uint64_t>(value.selected()));
    }
    else if constexpr (kind == io::detail::FieldKind::Boolean)
        integer(out, value.boolean() ? 1U : 0U);
    else if constexpr (kind == io::detail::FieldKind::Array) {
        integer(out, value.size());
        for (std::uint64_t i = 0; i < value.size(); ++i) {
            const auto element = value.at(static_cast<core::Index>(i));
            if (!element) throw std::runtime_error("Checked schema array became invalid");
            append(out, *element);
        }
    } else append(out, value.record());
}
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--header', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    records = []
    for name, body in re.findall(r'struct (\w+) \{\s*\n(.*?)\n\};',
                                 args.header.read_text(encoding='utf-8'), re.S):
        fields = re.findall(r'(?:std::optional<)?Field<detail::FieldType::T\d+>>?\s+(\w+);', body)
        if fields:
            records.append((name, fields))
    assert len(records) == 191, len(records)
    assert len({name for name, _ in records}) == len(records)
    assert sum(len(fields) for _, fields in records) > 500
    declarations = '\n'.join(f'inline void append(Bytes&, const io::schema::{name}&);'
                             for name, _ in records)
    definitions = []
    for name, fields in records:
        lines = [f'inline void append(Bytes& out, const io::schema::{name}& value) {{',
                 f'    text(out, "{name}");']
        for field in fields:
            lines.append(f'    text(out, "{field}"); append(out, value.{field});')
        lines.append('}')
        definitions.append('\n'.join(lines))
    provenance = ('// Generated public traversal; no product value is used as an oracle.\n'
                  '// Interface SHA-256: ' + hashlib.sha256(args.header.read_bytes()).hexdigest() + '\n'
                  '// Generator SHA-256: ' + hashlib.sha256(Path(__file__).read_bytes()).hexdigest() + '\n')
    suffix = '''
inline Bytes serialize(const io::schema::Case& value) {
    Bytes out;
    append(out, value);
    return out;
}
} // namespace uvlas::test::schema_snapshot
'''
    args.output.write_bytes((provenance + PREFIX + declarations + '\n' + TEMPLATES + '\n'
                            + '\n'.join(definitions) + suffix).encode('utf-8'))
    print(f'{len(records)} named records / {sum(len(fields) for _, fields in records)} members')


if __name__ == '__main__':
    main()
