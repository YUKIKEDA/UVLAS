# Build, test and gate commands. Recipes marked "gate" are what CI runs.
# Windows recipes run inside the MSVC environment (tools/dev/msvc.ps1, ADR 0002).

set windows-shell := ["pwsh.exe", "-NoLogo", "-NoProfile", "-Command"]

msvc := "pwsh -NoProfile -File tools/dev/msvc.ps1"
clang_format := "uvx clang-format@21.1.2"
clang_tidy := "uvx clang-tidy@21.1.1"
python312 := "uv run --python 3.12 --no-project python"

default:
    @just --list

# gate (Windows): configure, build and test Debug and Release
[windows]
gate:
    {{msvc}} cmake --preset windows-msvc
    {{msvc}} cmake --build --preset windows-debug
    {{msvc}} ctest --preset windows-debug
    {{msvc}} cmake --build --preset windows-release
    {{msvc}} ctest --preset windows-release

# gate (Linux): configure, build and test Debug, Release and the sanitizer build
[linux]
gate:
    cmake --preset linux-debug && cmake --build --preset linux-debug && ctest --preset linux-debug
    cmake --preset linux-release && cmake --build --preset linux-release && ctest --preset linux-release
    cmake --preset linux-sanitize && cmake --build --preset linux-sanitize && ctest --preset linux-sanitize

# gate: formatting of product sources and tests
[windows]
format-check:
    {{clang_format}} --dry-run --Werror (git ls-files "src/*.cpp" "include/*.hpp" "tests/*.cpp" "qa/*.cpp")

[linux]
format-check:
    {{clang_format}} --dry-run --Werror $(git ls-files 'src/*.cpp' 'include/*.hpp' 'tests/*.cpp' 'qa/*.cpp')

# not a gate: rewrite product sources and tests in place
[windows]
format:
    {{clang_format}} -i (git ls-files "src/*.cpp" "include/*.hpp" "tests/*.cpp" "qa/*.cpp")

[linux]
format:
    {{clang_format}} -i $(git ls-files 'src/*.cpp' 'include/*.hpp' 'tests/*.cpp' 'qa/*.cpp')

# gate: Python tools: ruff (lint and format) and pyright in strict mode
lint-python:
    uvx ruff@0.17.4 check .
    uvx ruff@0.17.4 format --check .
    uvx pyright@1.1.408

# not a gate: apply ruff fixes and formatting to the Python tools
format-python:
    uvx ruff@0.17.4 check --fix .
    uvx ruff@0.17.4 format .

# gate: clang-tidy on product sources and tests, with the compile database of the configured preset
[windows]
lint:
    {{clang_tidy}} -p out/build/windows-msvc --quiet (git ls-files "src/*.cpp" "tests/*.cpp" "qa/*.cpp")

[linux]
lint:
    {{clang_tidy}} -p out/build/linux-debug --quiet $(git ls-files 'src/*.cpp' 'tests/*.cpp' 'qa/*.cpp')

# gate: every ID a test or a record names is defined; the traceability page is current
trace-check:
    {{python312}} -X utf8 tools/qa/trace.py check

# not a gate: where an ID is defined, tested and recorded, for example `just trace N4-06`
trace id:
    {{python312}} -X utf8 tools/qa/trace.py show {{id}}

# not a gate: rewrite docs/quality/traceability.md after tests or documents change
trace-matrix:
    {{python312}} -X utf8 tools/qa/trace.py matrix

# gate: .claude/skills is a copy of .agents/skills
skills-check:
    python tools/dev/sync_skills.py --check

# not a gate: copy .agents/skills to .claude/skills
skills-sync:
    python tools/dev/sync_skills.py

# gate: links, script compilation and recorded hashes of the R&D tree
rd-check:
    {{python312}} -X utf8 tools/rd/check_tree.py

# not a gate (slow): the organized R&D tree equals the transform of the baseline tag
rd-verify:
    {{python312}} -X utf8 tools/rd/reorganize.py verify --rev HEAD
