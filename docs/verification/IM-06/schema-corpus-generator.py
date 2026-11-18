"""Freeze schema examples as token-preserving C++ literals; no product imports."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path


class Number:
    def __init__(self, token: str):
        self.token = token


def encode(value: object) -> str:
    if isinstance(value, Number):
        return value.token
    if isinstance(value, dict):
        return '{'+','.join(json.dumps(key,ensure_ascii=False)+':'+encode(item)
                            for key,item in value.items())+'}'
    if isinstance(value, list):
        return '['+','.join(encode(item) for item in value)+']'
    return json.dumps(value,ensure_ascii=False,separators=(',',':'))


def literal(value: object) -> str:
    text = encode(value)
    assert ')UVLAS_SCHEMA"' not in text
    # Adjacent tokens avoid MSVC's individual string-literal length limit.
    return '\n        '.join('R"UVLAS_SCHEMA('+text[i:i+7000]+')UVLAS_SCHEMA"'
                              for i in range(0,len(text),7000))


SCHEMA_NEGATIVE = {
    'point_only_underdetermined','missing_shell_state','duplicate_shell_state',
    'missing_velocity_field','missing_director_numerator','nonpositive_denominator',
    'degree_over_limit','duplicate_cell','negative_point_tolerance','missing_side_t',
    'unknown_point_target','missing_beam_endpoint','missing_rigid_state',
}
PREPARE_NEGATIVE = {
    'nonmember_c_squared','missing_active_cell','discontinuous_cell_field',
    'point_value_contradiction','outward_side_at_s1','support_contradiction',
    'non_so3_tip','rigid_attachment_contradiction','body_world_velocity_contradiction',
}


def patched(base: dict, patches: list[dict]) -> dict:
    result = copy.deepcopy(base)
    for patch in patches:
        parts = [part.replace('~1','/').replace('~0','~') for part in patch['path'].split('/')[1:]]
        parent = result
        for key in parts[:-1]:
            parent = parent[int(key)] if isinstance(parent,list) else parent[key]
        key = parts[-1]
        if isinstance(parent,list):
            if patch['op']=='remove':
                del parent[int(key)]
            elif patch['op']=='replace':
                parent[int(key)] = copy.deepcopy(patch['value'])
            elif patch['op']=='add':
                parent.insert(len(parent) if key=='-' else int(key),copy.deepcopy(patch['value']))
            else:
                raise AssertionError('unknown patch operation')
        elif patch['op']=='remove':
            del parent[key]
        elif patch['op'] in ('add','replace'):
            if patch['op']=='replace':
                assert key in parent
            parent[key] = copy.deepcopy(patch['value'])
        else:
            raise AssertionError('unknown patch operation')
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--schema',type=Path,required=True)
    parser.add_argument('--initial',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('output must be new')
    load = lambda path:json.loads(path.read_bytes(),parse_int=Number,parse_float=Number)
    schema,initial = load(args.schema),load(args.initial)
    assert len(schema['fragments'])==144 and len(schema['schema_examples'])==5
    assert len(initial['cases'])==4 and len(initial['negative_cases'])==22
    assert {item['id'] for item in initial['negative_cases']}==SCHEMA_NEGATIVE|PREPARE_NEGATIVE
    # These checks reject the usual float-round-trip corruption without relying on product code.
    witness = json.loads('{"tiny":1e-400,"low":1.000000000000000000000000000001,"signed":-0}',
                         parse_int=Number,parse_float=Number)
    assert encode(witness)=='{"tiny":1e-400,"low":1.000000000000000000000000000001,"signed":-0}'
    assert json.loads(encode(witness),parse_int=Number,parse_float=Number)['tiny'].token=='1e-400'
    header = ['#pragma once','#include <array>','#include <cstdint>','#include <string_view>',
        'namespace uvlas::test::case_schema_corpus {',
        'enum class ExpectedStage : std::uint8_t { Schema, Prepare };',
        'struct Fragment { std::string_view name, json; bool superseded; };',
        'struct Root { std::string_view name, json; };',
        'struct Physical { std::string_view name, context_json, perturbation_json; };',
        'struct Negative { std::string_view name, base_name, context_json, perturbation_json; ExpectedStage stage; };']
    for name,path in [('schema_source',args.schema),('initial_source',args.initial),('generator',Path(__file__))]:
        header.append(f'inline constexpr std::string_view {name}_sha256 = "{hashlib.sha256(path.read_bytes()).hexdigest()}";')
    header.append('inline constexpr std::array<Fragment,144> fragments{{')
    for name,value in schema['fragments'].items():
        superseded = name in ('Perturbation.physical','PointValue')
        header.append('{'+json.dumps(name)+','+literal(value)+','+str(superseded).lower()+'},')
    header.extend(['}};','inline constexpr std::array<Root,5> roots{{'])
    for name,value in schema['schema_examples'].items():
        header.append('{'+json.dumps(name)+','+literal(value)+'},')
    header.extend(['}};','inline constexpr std::array<Physical,4> physical{{'])
    for name,value in initial['cases'].items():
        header.append('{'+json.dumps(name)+','+literal(value['context'])+','+literal(value['perturbation'])+'},')
    header.extend(['}};','inline constexpr std::array<Negative,22> negatives{{'])
    for item in initial['negative_cases']:
        stage = 'Schema' if item['id'] in SCHEMA_NEGATIVE else 'Prepare'
        assert item['expected_stage']==stage.lower() and item['expected_failure']=='InvalidInput'
        value = patched(initial['cases'][item['base_case']],item['patches'])
        header.append('{'+json.dumps(item['id'])+','+json.dumps(item['base_case'])+','
                      +literal(value['context'])+','+literal(value['perturbation'])+',ExpectedStage::'+stage+'},')
    header.extend(['}};','} // namespace uvlas::test::case_schema_corpus',''])
    args.output.write_bytes('\n'.join(header).encode('utf-8'))
    print('144 fragments / 5 roots / 4 physical states / 22 patched negatives; exact tokens retained')


if __name__=='__main__':
    main()
