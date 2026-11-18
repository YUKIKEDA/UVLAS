"""Independent record/member inventory from the adopted Markdown dictionaries."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def split_fields(text):
    result, current, depth = [], [], 0
    for char in text:
        if char in '<{(':
            depth += 1
        elif char in '>})':
            depth -= 1
        if char == ',' and depth == 0:
            result.append(''.join(current).strip())
            current = []
        else:
            current.append(char)
    result.append(''.join(current).strip())
    assert depth == 0
    return result


def record_name(text):
    plain, depth = [], 0
    for char in text:
        if char == '<':
            depth += 1
        elif char == '>':
            depth -= 1
        elif depth == 0:
            plain.append(char)
    return ''.join(part[:1].upper() + part[1:]
                   for part in re.split(r'[\s/_]+', ''.join(plain).strip('`')) if part)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--header', type=Path)
    args = parser.parse_args()
    expected = {}
    sources = {}
    for name in ('case-schema-fields.md', 'case-initial-physical-state-contract.md'):
        path = args.root / 'docs/specs/foundation' / name
        sources[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        text = path.read_text(encoding='utf-8')
        for line in text.splitlines():
            if not line.startswith('| '):
                continue
            parts = line.split('|')
            if len(parts) != 5:
                continue
            label, fields = parts[1].strip(), parts[2].strip()
            if label == '型' or ':' not in fields:
                continue
            entries = []
            for field in split_fields(fields):
                if ':' not in field:
                    raise ValueError((label, field))
                key, spelling = field.split(':', 1)
                optional = key.endswith('?')
                key = key.removesuffix('?').strip()
                assert re.fullmatch('[a-z][a-z0-9_]*', key), (label, key)
                entries.append({'name': key, 'optional': optional, 'spelling': spelling.strip()})
            expected[record_name(label)] = entries
    # The full-state supplement supersedes the old point-only value path.
    expected.pop('PointValue', None)
    contract = args.root / 'docs/specs/foundation/case-schema-contract.md'
    sources[contract.name] = hashlib.sha256(contract.read_bytes()).hexdigest()
    contract_text = contract.read_text(encoding='utf-8')
    assert '`{numerator: Text, denominator: Text}`' in contract_text
    assert '`{source_id: Id<Source>, path: Text, sha256: Digest, role: Tag{指定role}}`' in contract_text
    expected['Rational'] = [{'name': key, 'optional': False, 'spelling': 'Text'}
                            for key in ('numerator', 'denominator')]
    expected['FileReference'] = [
        {'name': key, 'optional': False, 'spelling': spelling}
        for key, spelling in [('source_id', 'Id<Source>'), ('path', 'Text'),
                              ('sha256', 'Digest'), ('role', 'Tag{指定role}')]]
    assert len(expected) == 191, len(expected)
    assert sum(map(len, expected.values())) == 806
    result = {'sources': sources, 'records': expected}
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.write_bytes((json.dumps(result, indent=2) + '\n').encode())
    if args.header:
        if args.header.exists():
            raise FileExistsError(args.header)
        lines = ['#pragma once', '#include <array>', '#include <string_view>',
                 'namespace uvlas::test::schema_field_inventory {',
                 'struct Field { std::string_view record, name, spelling; bool optional; };',
                 'inline constexpr std::array<Field, 806> fields{{']
        for name, fields in expected.items():
            for field in fields:
                values = ', '.join(json.dumps(value, ensure_ascii=False)
                                   for value in (name, field['name'], field['spelling']))
                lines.append('    {' + values + ', ' + str(field['optional']).lower() + '},')
        lines += ['}};', '} // namespace uvlas::test::schema_field_inventory', '']
        args.header.write_bytes('\n'.join(lines).encode('utf-8'))
    print(json.dumps({'records': len(expected), 'members': sum(map(len, expected.values()))}))


if __name__ == '__main__':
    main()
