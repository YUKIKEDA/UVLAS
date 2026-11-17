"""Extract actual allocation preambles for a limited GCC11/ASan regression."""
import hashlib
import json
from pathlib import Path

root = Path.cwd().resolve()
area = root / "out/issue104-ci-resource"
old = (area / "old-resource.cpp").read_bytes()
new = (root / "tests/io/input_document_resource_test.cpp").read_bytes()
route_start = new.index(b'TEST_CASE("Document fault allocator covers')
route_end = new.index(b'TEST_CASE("Document allocation failures', route_start)
routes = new[route_start:route_end]
headers = b"""#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <memory>
#include <new>
"""
controls = b"""
unsigned check_count = 0;
void check(bool value, unsigned line) {
    ++check_count;
    if (!value) {
        std::fprintf(stderr, "failed control at line %u\\n", line);
        std::exit(42);
    }
}
#define TEST_CASE(...) void check_routes()
#define CHECK(...) check((__VA_ARGS__), __LINE__)
#define REQUIRE(...) check((__VA_ARGS__), __LINE__)
#define CAPTURE(...) static_cast<void>(__VA_ARGS__)
"""
main = b"""
int main() {
    // Match the library's nothrow temporary-buffer allocation + sized delete.
    void* const memory = ::operator new(32, std::nothrow);
    ::operator delete(memory, std::size_t{32});
    check_routes();
    std::printf("allocation control checks=%u\\n", check_count);
}
"""
report = {"scope": "Actual fault-TU preambles with one identical standalone wrapper; no product or Catch2 build on GCC11.",
          "cases": {}}
for label, source in (("old", old), ("new", new)):
    start = source.index(b"namespace {")
    end = source.index(b"namespace {", start + 1)
    preamble = source[start:end]
    (area / f"{label}-preamble.bytes").write_bytes(preamble)
    unit = headers + preamble + controls + routes + main
    (area / f"{label}-asan.cpp").write_bytes(unit)
    report["cases"][label] = {"full_tu_sha256": hashlib.sha256(source).hexdigest(),
                              "preamble_sha256": hashlib.sha256(preamble).hexdigest(),
                              "preamble_bytes": len(preamble),
                              "standalone_sha256": hashlib.sha256(unit).hexdigest()}
(area / "asan-inputs.json").write_bytes((json.dumps(report, indent=2) + "\n").encode())
