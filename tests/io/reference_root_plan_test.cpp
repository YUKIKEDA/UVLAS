#include "reference_declarations_fixture.hpp"
#include "reference_root_plan.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <limits>

namespace {
using namespace uvlas;
namespace fixture = test::declaration_fixture;
namespace rd = io::detail;
using fixture::take;

constexpr rd::RootReferenceLimits wide{
    fixture::budget, fixture::budget, fixture::budget, fixture::budget, 64,
    fixture::budget, fixture::budget};

rd::DeclaredReferences declared(const fixture::JsonToken& json) {
    const auto schema = fixture::checked(json);
    core::ResourcePlan plan{fixture::budget};
    return take(rd::collect_case_file_references(schema, fixture::limits(), plan, {}, nullptr));
}

fixture::JsonToken repeated() {
    auto json = fixture::all_roles();
    auto duplicate = json.field("observations").elements[0];
    duplicate.set("id", fixture::json("\"2\""));
    duplicate.field("comparison")
        .field("reference")
        .field("source")
        .set("source_id", fixture::json("\"0900\""));
    json.field("observations").elements.push_back(std::move(duplicate));
    return json;
}

fixture::JsonToken& mesh(fixture::JsonToken& json) {
    return json.field("model").field("geometry").elements[0].field("mesh");
}

fixture::JsonToken& observation(fixture::JsonToken& json, std::size_t index) {
    return json.field("observations")
        .elements[index]
        .field("comparison")
        .field("reference")
        .field("source");
}

void expect_failure(const core::Result<rd::PreparedRootReferences>& result,
                    core::FailureCode code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    CHECK(result.error().operation == core::Operation::PrepareModel);
}

void check_rejection(rd::DeclaredReferences& owner, rd::RootReferenceLimits caps,
                     core::ResourcePlan& plan, core::FailureCode code,
                     rd::RootReferenceDiagnostic* diagnostic, const std::stop_token& stop = {}) {
    // The rvalue factory transfers ownership only after admission succeeds.
    expect_failure(rd::prepare_root_references(std::move(owner), caps, plan, stop, diagnostic),
                   code);
}

void expect_reason(const rd::RootReferenceDiagnostic& diagnostic,
                   rd::RootReferenceFailure expected) {
    if (const auto& detail = diagnostic.detail()) {
        CHECK(detail->reason == expected);
    } else {
        FAIL("Missing reason");
    }
}

} // namespace

TEST_CASE("Root reference preparation owns paths and every declaration in dictionary order",
          "[IO4-01][N4-02][integration][property]") {
    constexpr std::array<std::uint64_t, 9> ids{300, 301, 302, 200, 201, 700, 900, 900, 800};
    constexpr std::array<core::Index, 9> ordinals{0, 1, 2, 3, 4, 5, 6, 6, 7};
    constexpr std::array<std::string_view, 9> roles{
        "gmsh-msh4.1-ascii",     "gmsh-sidecar",          "reference-geometry",
        "nastran-deck",          "nastran-sidecar",       "preparation-history",
        "observation-reference", "observation-reference", "initial-history"};
    constexpr std::array<std::string_view, 9> pointers{
        "/model/geometry/0/mesh",
        "/model/geometry/0/sidecar",
        "/model/geometry/0/reference/file",
        "/model/geometry/1/deck",
        "/model/geometry/1/sidecar",
        "/analysis/branch/seed/history",
        "/observations/0/comparison/reference/source",
        "/observations/1/comparison/reference/source",
        "/initial/file"};
    constexpr std::array<std::string_view, 8> paths{
        "folder/mesh.msh",    "sidecar.json", "参照/geometry.json", "model.bdf",
        "model-sidecar.json", "seed.history", "observation.dat",    "initial.history"};

    auto json = repeated();
    fixture::reverse_keys(json);
    auto original = declared(json);
    const auto parent = take(original.identity());
    auto prepared = [&] {
        core::ResourcePlan plan{fixture::budget};
        return take(rd::prepare_root_references(std::move(original), wide, plan, {}, nullptr));
    }();
    // NOLINTNEXTLINE(bugprone-use-after-move): successful preflight transfers source ownership.
    CHECK_FALSE(original.size());
    const auto moved = std::move(prepared);
    // NOLINTNEXTLINE(bugprone-use-after-move): the moved-from owner is intentionally tested.
    CHECK_FALSE(prepared.report());
    CHECK(fixture::same_identity(parent, take(moved.identity())));
    const auto report = take(moved.report());
    CHECK(report.declarations == ids.size());
    CHECK(report.unique_sources == paths.size());
    // 9 declaration slots, 9 possible sources, 19 buckets, 9 declarations, and 11
    // lookups because 700 probes the buckets occupied by 301 and 302 first.
    // The eight unique original paths total 128 bytes, each charged at 16 units.
    constexpr std::uint64_t expected_work = 9 + 9 + 19 + 9 + 11 + 16 * 128;
    CHECK(report.work == expected_work);

    for (std::size_t i = 0; i < ids.size(); ++i) {
        const auto entry = take(moved.declaration(static_cast<core::Index>(i)));
        CHECK(entry.original.source_id.value() == ids[i]);
        CHECK(rd::reference_role_name(entry.original.role) == roles[i]);
        CHECK(entry.original.declaration_pointer == pointers[i]);
        CHECK(entry.source_ordinal == ordinals[i]);
    }
    CHECK(take(moved.declaration(0)).original.original_path == R"(folder\./mesh.msh)");
    CHECK(take(moved.declaration(7)).original.source_id_spelling == "0900");
    std::uint64_t total = 0;
    for (std::size_t i = 0; i < paths.size(); ++i) {
        const auto expected_first = i == 7 ? 8 : i;
        const auto source = take(moved.source(static_cast<core::Index>(i)));
        CHECK(source.first_declaration == static_cast<core::Index>(expected_first));
        CHECK(source.path.text() == paths[i]);
        total += paths[i].size();
    }
    CHECK(report.normalized_path_bytes == total);
    CHECK_FALSE(moved.declaration(-1));
    CHECK_FALSE(moved.declaration(9));
    CHECK_FALSE(moved.source(-1));
    CHECK_FALSE(moved.source(8));

    auto empty = declared(fixture::root());
    core::ResourcePlan no_payload{0};
    const auto zero =
        take(rd::prepare_root_references(std::move(empty), {}, no_payload, {}, nullptr));
    CHECK(take(zero.report()).declarations == 0);
    CHECK(take(zero.report()).unique_sources == 0);
    CHECK(take(zero.report()).work == 0);
    CHECK(no_payload.allocated_bytes() == 0);
}

TEST_CASE("Equal normalized path stays distinct by SourceId and bad paths keep first provenance",
          "[IO4-01][N4-08][integration][negative]") {
    auto json = fixture::all_roles();
    fixture::reference(mesh(json), "0", "./observation.dat", '3');
    fixture::reference(observation(json, 0), "18446744073709551615", "observation.dat", '8');
    auto original = declared(json);
    core::ResourcePlan plan{fixture::budget};
    const auto paths =
        take(rd::prepare_root_references(std::move(original), wide, plan, {}, nullptr));
    CHECK(take(paths.declaration(0)).original.source_id.value() == 0);
    CHECK(take(paths.declaration(6)).original.source_id.value() ==
          std::numeric_limits<std::uint64_t>::max());
    CHECK(take(paths.declaration(0)).source_ordinal != take(paths.declaration(6)).source_ordinal);
    const auto mesh_source = take(paths.source(0));
    const auto observation_source = take(paths.source(6));
    CHECK(mesh_source.path.text() == observation_source.path.text());

    fixture::reference(mesh(json), "0", "../outside", '3');
    fixture::reference(observation(json, 0), "900", "/absolute", '8');
    auto bad = declared(json);
    rd::RootReferenceDiagnostic diagnostic;
    constexpr std::string_view first_pointer = "/model/geometry/0/mesh";
    auto exact_pointer = wide;
    exact_pointer.diagnostic_bytes = first_pointer.size();
    check_rejection(bad, exact_pointer, plan, core::FailureCode::InvalidInput, &diagnostic);
    expect_reason(diagnostic, rd::RootReferenceFailure::Path);
    CHECK(diagnostic.pointer() == first_pointer);
    if (const auto& detail = diagnostic.detail(); detail && detail->path) {
        CHECK(detail->source == io::SourceId{0});
        CHECK(detail->declaration_node == take(bad.entry(0)).declaration_node);
        CHECK(detail->path->reason == rd::ReferencePathReason::EscapeRoot);
        CHECK(detail->path->input_offset == 0);
    } else {
        FAIL("Missing path detail");
    }
    check_rejection(bad, wide, plan, core::FailureCode::InvalidInput, nullptr);

    auto short_pointer = exact_pointer;
    --short_pointer.diagnostic_bytes;
    check_rejection(bad, short_pointer, plan, core::FailureCode::InvalidInput, &diagnostic);
    expect_reason(diagnostic, rd::RootReferenceFailure::Path);
    CHECK_FALSE(diagnostic.pointer());
    CHECK(take(paths.report()).unique_sources == 8);
}

TEST_CASE("Case admission rejects same SourceId with a different original declaration",
          "[IO4-01][integration][negative]") {
    for (const auto field :
         {std::string_view{"path"}, std::string_view{"sha256"}, std::string_view{"role"}}) {
        CAPTURE(field);
        auto json = repeated();
        auto& duplicate = observation(json, 1);
        if (field == "role") {
            fixture::reference(mesh(json), "900", "observation.dat", '8');
        } else if (field == "path") {
            duplicate.set("path", fixture::json("\"./observation.dat\""));
        } else {
            duplicate.set("sha256", fixture::json("\"" + std::string(64, '9') + "\""));
        }
        const auto parsed = fixture::parsed(json);
        core::ResourcePlan plan{fixture::budget};
        std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> caps{};
        caps.fill(fixture::budget);
        io::SchemaDiagnostic diagnostic;
        const auto result =
            io::validate_case_schema(parsed, io::SchemaLimits{caps}, plan, &diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::InvalidInput);
        if (const auto detail = diagnostic.detail()) {
            CHECK(detail->reason == io::SchemaReason::DuplicateIdentifier);
        } else {
            FAIL("Missing identity detail");
        }
    }
}

TEST_CASE("Root reference limits admit exact resources and reject one-unit deficits",
          "[IO4-01][N4-02][N4-08][negative]") {
    const auto json = repeated();
    auto baseline = declared(json);
    core::ResourcePlan probe{fixture::budget};
    const auto stable =
        take(rd::prepare_root_references(std::move(baseline), wide, probe, {}, nullptr));
    const auto result = take(stable.report());
    CHECK(result.cumulative_bytes == probe.allocated_bytes());

    const rd::RootReferenceLimits exact{9, 9, 8, result.normalized_path_bytes, 2, wide.work, 0};
    auto first = declared(json);
    core::ResourcePlan initial{fixture::budget};
    const auto prepared =
        take(rd::prepare_root_references(std::move(first), exact, initial, {}, nullptr));
    const auto totals = take(prepared.report());
    auto measured = exact;
    measured.work = totals.work;
    auto on_bound = declared(json);
    core::ResourcePlan exact_plan{totals.cumulative_bytes};
    const auto exact_result =
        rd::prepare_root_references(std::move(on_bound), measured, exact_plan, {}, nullptr);
    REQUIRE(exact_result);
    CHECK(exact_plan.remaining_bytes() == 0);

    using R = rd::RootReferenceFailure;
    constexpr std::array causes{R::Declarations, R::RootEdges, R::UniqueSources,
                                R::Path,         R::Path,      R::Work};
    auto against = declared(json);
    for (std::size_t i = 0; i < causes.size(); ++i) {
        auto reduced = measured;
        // Changing capacity also changes the index size, so allow full probe work in that case.
        switch (i) {
        case 0:
            --reduced.declarations;
            break;
        case 1:
            --reduced.root_edges;
            break;
        case 2:
            --reduced.unique_sources;
            reduced.work = wide.work;
            break;
        case 3:
            --reduced.normalized_path_bytes;
            break;
        case 4:
            --reduced.path_components;
            break;
        default:
            --reduced.work;
            break;
        }
        core::ResourcePlan plan{fixture::budget};
        rd::RootReferenceDiagnostic diagnostic;
        check_rejection(against, reduced, plan, core::FailureCode::ResourceExceeded, &diagnostic);
        expect_reason(diagnostic, causes[i]);
        CHECK_FALSE(diagnostic.pointer());
        CHECK(take(against.size()) == 9);
    }
    core::ResourcePlan short_plan{totals.cumulative_bytes - 1};
    check_rejection(against, measured, short_plan, core::FailureCode::ResourceExceeded, nullptr);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() requires a mutable source.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    core::ResourcePlan cancelled_plan{fixture::budget};
    check_rejection(against, measured, cancelled_plan, core::FailureCode::Cancelled, nullptr,
                    cancel.get_token());
    CHECK(cancelled_plan.allocated_bytes() == 0);
    CHECK(take(against.size()) == 9);
    CHECK(take(stable.report()).work == result.work);
}
