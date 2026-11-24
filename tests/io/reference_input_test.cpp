#include "../structure/pg_partition_test_support.hpp"
#include "coupling/case_external_pg.hpp"
#include "coupling/case_nastran_reference_inputs.hpp"
#include "external_mesh_sidecar_pg.hpp"
#include "gmsh_mesh.hpp"
#include "nastran_coordinate_cards.hpp"
#include "nastran_pg_material_axes.hpp"
#include "reference_declarations_fixture.hpp"
#include "reference_input.hpp"
#include "reference_resolution.hpp"
#include "reference_root_snapshot.hpp"
#include "source_digest.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <compare>
#include <fstream>
#include <type_traits>

namespace {
namespace io = uvlas::io;
namespace core = uvlas::core;
namespace rd = uvlas::io::detail;
namespace fixture = uvlas::test::declaration_fixture;
namespace fs = std::filesystem;

constexpr io::SourceId source{19};
constexpr rd::ManifestLimits manifest_limits{
    128, {fixture::budget, fixture::budget, fixture::budget, fixture::budget}};
constexpr std::string_view name = "\xc3\xa9.json";
constexpr char8_t native_name[] = u8"\u00e9.json";

static_assert(!std::is_constructible_v<rd::ReferenceInput, io::SchemaCase, rd::ReadRoot,
                                       rd::ReferenceFileId>);
static_assert(!std::is_copy_constructible_v<rd::ReferenceInput>);

struct Tree {
    fs::path root;

    Tree() {
        const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
        for (int attempt = 0; attempt < 16; ++attempt) {
            root = fs::temp_directory_path() / ("uvlas-reference-input-" + std::to_string(tick) +
                                                "-" + std::to_string(attempt));
            if (fs::create_directory(root)) {
                return;
            }
        }
        FAIL("Could not create a private test directory");
    }
    ~Tree() {
        std::error_code ignore;
        fs::remove_all(root, ignore);
    }
    Tree(const Tree&) = delete;
    Tree& operator=(const Tree&) = delete;
};

void write_file(const fs::path& target, std::string_view bytes) {
    std::ofstream out(target, std::ios::binary);
    out << bytes;
    out.close();
    REQUIRE(out.good());
}

std::span<const std::byte> bytes_of(std::string_view text) noexcept {
    return std::as_bytes(std::span{text});
}

rd::ReferenceInputLimits limits() {
    const io::InputLimits input{
        10'000'000, 256,
        1'000'000,  1'000'000,
        10'000'000, 10'000'000,
        10'000'000, {io::ParserAllocationMode::EstimateAdmission, fixture::budget}};
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> caps{};
    caps.fill(fixture::budget);
    return {input, io::SchemaLimits{caps}};
}

constexpr std::string_view static_case() {
    for (const auto& example : uvlas::test::case_schema_corpus::roots) {
        if (example.name == "structural_static") {
            return example.json;
        }
    }
    return {};
}

consteval io::SourceDigest digest_literal(std::string_view hex) {
    constexpr std::string_view digits = "0123456789abcdef";
    io::SourceDigest digest{};
    for (std::size_t i = 0; i < digest.bytes.size(); ++i) {
        digest.bytes[i] =
            static_cast<std::byte>(16 * digits.find(hex[2 * i]) + digits.find(hex[2 * i + 1]));
    }
    return digest;
}

// Python 3.12 hashlib SHA-256 over both concatenated structural_static corpus literals.
constexpr auto static_case_digest =
    digest_literal("ecfba5a3cdf129b195b7e1504c74c2981852f6ee57617ac320ea14420ffd9a38");
constexpr auto abc_digest =
    digest_literal("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");

void check_failure(const core::Result<rd::ReferenceInput>& result, core::FailureCode code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().operation == core::Operation::PrepareModel);
    CHECK(result.error().code == code);
}

#ifdef __linux__
std::size_t held_descriptors(const fs::path& path) {
    std::size_t result = 0;
    for (const auto& entry : fs::directory_iterator("/proc/self/fd")) {
        std::error_code ignore;
        if (fs::read_symlink(entry.path(), ignore) == path) {
            ++result;
        }
    }
    return result;
}
#endif

} // namespace

TEST_CASE("File and Memory reference inputs retain one checked case and their actual origins",
          "[IO4-01][N4-02][integration][property]") {
    static_assert(static_case().size() == 8611);
    const Tree tree;
    const auto sub = tree.root / "child";
    REQUIRE(fs::create_directory(sub));
    write_file(tree.root / native_name, "not the manifest");
    write_file(sub / native_name, static_case());
    write_file(tree.root / "copy.bin", "abd");
    write_file(sub / "copy.bin", "abc");
    fs::create_hard_link(sub / native_name, sub / "another-name.json");

    {
        auto file = [&] {
            core::ResourcePlan plan{fixture::budget};
            rd::ReferenceInputDiagnostic diagnostic;
            auto result = rd::ReferenceInput::from_file(sub, name, source, manifest_limits,
                                                        limits(), plan, {}, &diagnostic);
            CHECK_FALSE(diagnostic.stage);
            return fixture::take(std::move(result));
        }();
        auto owned = std::move(file);
        // NOLINTNEXTLINE(bugprone-use-after-move): moved-from rejection is part of the contract.
        CHECK_FALSE(file.view());

        const auto view = fixture::take(owned.view());
        REQUIRE(view.root != nullptr);
        const auto identity = fixture::take(view.schema.identity());
        CHECK(identity.source.source == source);
        CHECK(identity.source.byte_count == static_case().size());
        CHECK(identity.source.digest == static_case_digest);

        auto memory = [&] {
            std::string buffer{static_case()};
            core::ResourcePlan plan{fixture::budget};
            auto result = fixture::take(rd::ReferenceInput::from_memory(
                bytes_of(buffer), source, std::nullopt, limits(), plan, {}, nullptr));
            buffer.assign(buffer.size(), '!');
            return result;
        }();
        const auto memory_view = fixture::take(memory.view());
        CHECK(memory_view.root == nullptr);
        CHECK_FALSE(memory_view.manifest_identity);
        CHECK(fixture::same_identity(identity, fixture::take(memory_view.schema.identity())));

#ifdef __linux__
        CHECK(held_descriptors(sub) == 1);
        CHECK(held_descriptors(sub / native_name) == 1);
#else
        std::error_code manifest_rename, root_rename;
        fs::rename(sub / native_name, sub / "renamed.json", manifest_rename);
        fs::rename(sub, tree.root / "renamed-child", root_rename);
        REQUIRE(manifest_rename);
        REQUIRE(root_rename);
#endif

        core::ResourcePlan plan{fixture::budget};
        if (const auto manifest_identity = view.manifest_identity) {
            const auto alias_path = fixture::take(
                rd::normalize_reference_path("another-name.json", {64, 64, 1}, plan, {}, nullptr));
            auto alias = fixture::take(view.root->open_file(alias_path, plan, {}));
            CHECK(fixture::take(alias.info()).identity == *manifest_identity);
        } else {
            FAIL("File input did not retain the manifest identity");
        }

        const auto copy_path =
            fixture::take(rd::normalize_reference_path("copy.bin", {64, 64, 1}, plan, {}, nullptr));
        auto copy = fixture::take(view.root->open_file(copy_path, plan, {}));
        const auto snapshot = fixture::take(
            rd::read_reference_file(copy, abc_digest, {3, 4, 2, 3}, plan, {}, nullptr));
        if (const auto payload = snapshot.view()) {
            CHECK(payload->digest == abc_digest);
            CHECK(std::ranges::equal(payload->bytes, bytes_of("abc")));
        } else {
            FAIL("Read did not retain the external payload");
        }
    }
#ifdef __linux__
    CHECK(held_descriptors(sub) == 0);
    CHECK(held_descriptors(sub / native_name) == 0);
#else
    fs::rename(sub / native_name, sub / "renamed.json");
    fs::rename(sub, tree.root / "renamed-child");
#endif
}

TEST_CASE("Memory reference input requires a root for the eight declared roles",
          "[IO4-01][N4-02][integration][negative]") {
    const Tree tree;
    const auto supplied = tree.root / "supplied";
    REQUIRE(fs::create_directory(supplied));
    write_file(supplied / "copy.bin", "abc");
    const auto wire = fixture::all_roles().wire();
    const auto bytes = bytes_of(wire);
    core::ResourcePlan plan{fixture::budget};
    rd::ReferenceInputDiagnostic diagnostic;
    check_failure(rd::ReferenceInput::from_memory(bytes, source, std::nullopt, limits(), plan, {},
                                                  &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Root);
    check_failure(
        rd::ReferenceInput::from_memory(bytes, source, std::nullopt, limits(), plan, {}, nullptr),
        core::FailureCode::InvalidInput);

    auto held = fixture::take(rd::ReadRoot::open(supplied, {}));
#ifdef __linux__
    const auto relocated = tree.root / "relocated";
    fs::rename(supplied, relocated);
    REQUIRE(fs::create_directory(supplied));
    write_file(supplied / "copy.bin", "abd");
#endif
    const auto input = fixture::take(rd::ReferenceInput::from_memory(
        bytes, source, std::move(held), limits(), plan, {}, &diagnostic));
    CHECK_FALSE(diagnostic.stage);
    const auto view = fixture::take(input.view());
    REQUIRE(view.root != nullptr);
    CHECK_FALSE(view.manifest_identity);
    const auto copy_path =
        fixture::take(rd::normalize_reference_path("copy.bin", {64, 64, 1}, plan, {}, nullptr));
    auto copy = fixture::take(view.root->open_file(copy_path, plan, {}));
    const auto snapshot =
        fixture::take(rd::read_reference_file(copy, abc_digest, {3, 4, 2, 3}, plan, {}, nullptr));
    if (const auto payload = snapshot.view()) {
        CHECK(std::ranges::equal(payload->bytes, bytes_of("abc")));
    } else {
        FAIL("Read did not retain the external payload");
    }

    const auto declarations = fixture::take(
        rd::collect_case_file_references(view.schema, fixture::limits(), plan, {}, nullptr));
    CHECK(fixture::take(declarations.size()) == 8);
    CHECK(fixture::take(declarations.identity()).source.source == source);

    // NOLINTNEXTLINE(bugprone-use-after-move): a moved-from root cannot qualify an input.
    auto closed = std::optional{std::move(held)};
    check_failure(rd::ReferenceInput::from_memory(bytes, source, std::move(closed), limits(), plan,
                                                  {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Root);
}

TEST_CASE("Reference acquisition rejects invalid names, bytes, quota and cancellation",
          "[IO4-01][N4-08][integration][negative]") {
    const Tree tree;
    core::ResourcePlan plan{fixture::budget};
    rd::ReferenceInputDiagnostic diagnostic;
    for (const std::string_view bad :
         {"", ".", "..", "dir/file", "dir\\file", "/absolute", "a:b", "p/../file", "\xc0\x80",
          "\xe0\x9f\xbf", "\xed\xa0\x80", "\xf4\x90\x80\x80"}) {
        CAPTURE(bad);
        check_failure(rd::ReferenceInput::from_file(tree.root, bad, source, manifest_limits,
                                                    limits(), plan, {}, &diagnostic),
                      core::FailureCode::InvalidInput);
        CHECK(diagnostic.stage == rd::ReferenceInputStage::ManifestName);
    }
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                                plan, {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::ManifestOpen);

    write_file(tree.root / native_name, "{");
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                                plan, {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Document);
    CHECK(diagnostic.document.has_value());

    auto broken = fixture::root();
    broken.field("model").erase("materials");
    const auto wire = broken.wire();
    check_failure(rd::ReferenceInput::from_memory(bytes_of(wire), source, std::nullopt, limits(),
                                                  plan, {}, &diagnostic),
                  core::FailureCode::InvalidInput);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::Schema);
    CHECK(diagnostic.schema.detail().has_value());
    check_failure(rd::ReferenceInput::from_memory(bytes_of(wire), source, std::nullopt, limits(),
                                                  plan, {}, nullptr),
                  core::FailureCode::InvalidInput);

    write_file(tree.root / native_name, static_case());
    auto insufficient = manifest_limits;
    insufficient.file.maximum_file_bytes = static_case().size() - 1;
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, insufficient, limits(),
                                                plan, {}, &diagnostic),
                  core::FailureCode::ResourceExceeded);
    CHECK(diagnostic.stage == rd::ReferenceInputStage::ManifestRead);
    if (const auto file = diagnostic.file) {
        CHECK(file->reason == rd::ReferenceFileFailure::FileBytes);
    } else {
        FAIL("File quota rejection did not retain its reason");
    }

    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() is non-const.
    std::stop_source stopped;
    REQUIRE(stopped.request_stop());
    const auto before = plan.allocated_bytes();
    check_failure(rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                                plan, stopped.get_token(), &diagnostic),
                  core::FailureCode::Cancelled);
    CHECK_FALSE(diagnostic.stage);
    check_failure(rd::ReferenceInput::from_memory(bytes_of(static_case()), source, std::nullopt,
                                                  limits(), plan, stopped.get_token(), &diagnostic),
                  core::FailureCode::Cancelled);
    CHECK_FALSE(diagnostic.stage);
    CHECK(plan.allocated_bytes() == before);
}

TEST_CASE("Both acquisition routes account for their complete allocation on one plan",
          "[IO4-01][N4-02][integration][negative]") {
    const Tree tree;
    write_file(tree.root / native_name, static_case());
    const auto check_budget = [](const auto& acquire) {
        core::ResourcePlan probe{fixture::budget};
        const auto input = fixture::take(acquire(probe));
        const auto view = fixture::take(input.view());
        const auto schema = fixture::take(view.schema.memory_report());
        REQUIRE(probe.allocated_bytes() >= schema.cumulative_bytes);
        REQUIRE(schema.admission_bytes >= schema.cumulative_bytes);
        // Schema admits a failure reserve as well as charging successful allocations.
        // The acquisition and parsing charges must still be present at schema admission.
        const auto required =
            probe.allocated_bytes() - schema.cumulative_bytes + schema.admission_bytes;
        REQUIRE(required > 0);

        core::ResourcePlan exact{required};
        REQUIRE(acquire(exact));
        CHECK(exact.allocated_bytes() == probe.allocated_bytes());

        core::ResourcePlan one_short{required - 1};
        check_failure(acquire(one_short), core::FailureCode::ResourceExceeded);
    };
    check_budget([&](core::ResourcePlan& plan) {
        return rd::ReferenceInput::from_file(tree.root, name, source, manifest_limits, limits(),
                                             plan, {}, nullptr);
    });
    check_budget([](core::ResourcePlan& plan) {
        return rd::ReferenceInput::from_memory(bytes_of(static_case()), source, std::nullopt,
                                               limits(), plan, {}, nullptr);
    });
}

namespace {

constexpr std::string_view abc_hex =
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
constexpr std::string_view valid_manifest_hex =
    "ecfba5a3cdf129b195b7e1504c74c2981852f6ee57617ac320ea14420ffd9a38";
constexpr rd::RootSnapshotLimits acquisition_limits{{3, 4, 2, 3}, 3, 256};
constexpr rd::RootReferenceLimits root_limits{
    fixture::budget, fixture::budget, fixture::budget, fixture::budget, 64, fixture::budget, 0};

void declare_path(fixture::JsonToken& file, std::string_view path, std::string_view hash) {
    file.set("path", fixture::string(path));
    file.set("sha256", fixture::string(hash));
}

fixture::JsonToken manifest_with_reference(std::string_view first, std::string_view digest) {
    auto model = fixture::all_roles();
    auto& geometry = model.field("model").field("geometry").elements;
    declare_path(geometry[0].field("mesh"), first, digest);
    declare_path(geometry[0].field("sidecar"), "alias.bin", abc_hex);
    declare_path(geometry[0].field("reference").field("file"), "copy.bin", abc_hex);
    declare_path(geometry[1].field("sidecar"), "wrong.bin", abc_hex);
    return model;
}

rd::PreparedRootReferences prepare(const rd::ReferenceInput& input) {
    core::ResourcePlan plan{fixture::budget};
    const auto view = fixture::take(input.view());
    auto declared = fixture::take(
        rd::collect_case_file_references(view.schema, fixture::limits(), plan, {}, nullptr));
    return fixture::take(
        rd::prepare_root_references(std::move(declared), root_limits, plan, {}, nullptr));
}

rd::ReferenceInput from_memory(const fs::path& directory, std::string_view wire,
                               io::SourceId source_id = source) {
    core::ResourcePlan plan{fixture::budget};
    auto root = fixture::take(rd::ReadRoot::open(directory, {}));
    return fixture::take(rd::ReferenceInput::from_memory(bytes_of(wire), source_id, std::move(root),
                                                         limits(), plan, {}, nullptr));
}

template <class T>
void check_acquisition_failure(const core::Result<T>& outcome, core::FailureCode expected) {
    REQUIRE_FALSE(outcome);
    CHECK(outcome.error().code == expected);
    CHECK(outcome.error().operation == core::Operation::PrepareModel);
}

void check_read_reason(const rd::RootSnapshotDiagnostic& diagnostic,
                       rd::ReferenceFileFailure reason) {
    CHECK(diagnostic.reason == rd::RootSnapshotFailure::Read);
    if (const auto& observed = diagnostic.read) {
        CHECK(observed->reason == reason);
        CHECK(observed->before.has_value());
    } else {
        FAIL("Missing same-handle diagnostic");
    }
}

constexpr std::string_view inventory_mesh = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 11 99
2 1 0 3
99
11
42
0 1 0
0 0 0
1.000 0 0
$EndNodes
$Elements
1 1 91 91
2 1 2 1
91 11 42 99
$EndElements
)";
constexpr std::string_view inventory_mesh_hash =
    "fb2e2468f3fdb30ac0980510dba4adf79db7f8578b3ca88984465297f00c7d79";
constexpr std::array inventory_ids{300U, 301U, 302U, 200U, 201U, 700U, 900U, 800U};
constexpr std::array inventory_roles{rd::ReferenceRole::GmshMesh,
                                     rd::ReferenceRole::GmshSidecar,
                                     rd::ReferenceRole::ReferenceGeometry,
                                     rd::ReferenceRole::NastranDeck,
                                     rd::ReferenceRole::NastranSidecar,
                                     rd::ReferenceRole::PreparationHistory,
                                     rd::ReferenceRole::ObservationReference,
                                     rd::ReferenceRole::InitialHistory};
constexpr std::uint64_t inventory_payload = inventory_mesh.size() + std::uint64_t{7} * 3;
constexpr rd::RootInventoryLimits memory_inventory_bounds{
    {{inventory_mesh.size(), inventory_mesh.size() + 1, 2, inventory_mesh.size()}, 8, 256},
    inventory_ids.size(),
    inventory_payload,
    inventory_payload + inventory_ids.size(),
    2 * inventory_ids.size(),
    inventory_payload,
    28};

std::string write_root_inventory(const fs::path& directory, std::string_view deck_body = "abc",
                                 std::string_view deck_hash = abc_hex) {
    auto json = fixture::all_roles();
    auto& geometry = json.field("model").field("geometry").elements;
    auto& observation = json.field("observations")
                            .elements[0]
                            .field("comparison")
                            .field("reference")
                            .field("source");
    const auto write = [&](fixture::JsonToken& declaration, std::string_view filename,
                           std::string_view body = "abc", std::string_view hash = abc_hex) {
        write_file(directory / filename, body);
        declare_path(declaration, filename, hash);
    };
    write(geometry[0].field("mesh"), "skin.msh", inventory_mesh, inventory_mesh_hash);
    write(geometry[0].field("sidecar"), "sidecar.json");
    write(geometry[0].field("reference").field("file"), "geometry.json");
    write(geometry[1].field("deck"), "model.bdf", deck_body, deck_hash);
    write(geometry[1].field("sidecar"), "model-sidecar.json");
    write(json.field("analysis").field("branch").field("seed").field("history"), "seed.history");
    write(observation, "observation.dat");
    write(json.field("initial").field("file"), "initial.history");
    auto repeat = json.field("observations").elements[0];
    repeat.set("id", fixture::json("\"2\""));
    repeat.field("comparison")
        .field("reference")
        .field("source")
        .set("source_id", fixture::json("\"0900\""));
    json.field("observations").elements.push_back(std::move(repeat));
    return json.wire();
}

} // namespace

TEST_CASE("The root inventory owns every direct declaration and feeds the verified Gmsh mesh",
          "[IO4-01][IO4-04][N4-02][integration][verification]") {
    auto [inventory, parent] = [] {
        const Tree tree;
        const auto text = write_root_inventory(tree.root);
        write_file(tree.root / "case.json", text);
        core::ResourcePlan source_memory{fixture::budget};
        const auto input = fixture::take(rd::ReferenceInput::from_file(
            tree.root, "case.json", source, manifest_limits, limits(), source_memory, {}, nullptr));
        auto declarations = prepare(input);
        const auto expected_parent = fixture::take(declarations.identity());
        auto bounds = memory_inventory_bounds;
        bounds.total_identity_comparisons += inventory_ids.size();
        core::ResourcePlan inventory_memory{fixture::budget};
        rd::RootSnapshotDiagnostic diagnostic;
        auto files = fixture::take(rd::acquire_root_references(
            input, std::move(declarations), bounds, inventory_memory, {}, &diagnostic));
        CHECK_FALSE(diagnostic.reason);
        return std::pair{std::move(files), expected_parent};
    }(); // Manifest, all native handles, plans, and the tree are gone before consumption.
    CHECK(fixture::same_identity(parent, fixture::take(inventory.identity())));
    const auto count = fixture::take(inventory.report());
    CHECK(count.declarations == 9);
    CHECK(count.unique_sources == inventory_ids.size());
    constexpr std::array ordinals{0, 1, 2, 3, 4, 5, 6, 6, 7};
    for (std::size_t i = 0; i < ordinals.size(); ++i) {
        const auto index = static_cast<std::size_t>(ordinals[i]);
        const auto declaration = fixture::take(inventory.declaration(static_cast<core::Index>(i)));
        CHECK(declaration.source_ordinal == ordinals[i]);
        CHECK(declaration.original.source_id == io::SourceId{inventory_ids[index]});
        CHECK(declaration.original.role == inventory_roles[index]);
    }
    CHECK(fixture::take(inventory.declaration(6)).original.source_id_spelling == "900");
    CHECK(fixture::take(inventory.declaration(7)).original.source_id_spelling == "0900");
    for (std::size_t i = 0; i < inventory_ids.size(); ++i) {
        const auto reference = fixture::take(
            inventory.find(parent, io::SourceId{inventory_ids[i]}, inventory_roles[i], {}));
        const auto payload = reference.get().file.view();
        if (!payload) {
            FAIL("Missing authenticated file");
            return;
        }
        CHECK(std::ranges::equal(payload->bytes, bytes_of(i == 0 ? inventory_mesh : "abc")));
        CHECK(payload->digest == (i == 0 ? digest_literal(inventory_mesh_hash) : abc_digest));
        CHECK(payload->read_calls == 2);
        CHECK_FALSE(reference.get().pointer().empty());
    }
    const auto mesh =
        fixture::take(inventory.find(parent, io::SourceId{300}, rd::ReferenceRole::GmshMesh, {}));
    core::ResourcePlan mesh_memory{fixture::budget};
    const auto parsed = fixture::take(rd::GmshMesh::prepare(
        mesh.get(), {inventory_mesh.size(), 0, 0, 0, 3, 1, 4096}, mesh_memory, {}));
    const auto view = parsed.view();
    if (!view) {
        FAIL("Missing owned mesh");
        return;
    }
    CHECK(view->bytes == inventory_mesh);
    CHECK(view->location.source == io::SourceId{300});
    REQUIRE(view->nodes.size() == 3);
    CHECK(view->nodes[0].tag == 11);
    CHECK(view->nodes[1].tag == 42);
    CHECK(view->nodes[2].tag == 99);
    CHECK(view->text(view->nodes[1].coordinates[0]) == "1.000");
    REQUIRE(view->elements.size() == 1);
    CHECK(view->elements[0].tag == 91);
    constexpr std::array<std::uint64_t, 3> nodes{11, 42, 99};
    CHECK(std::ranges::equal(view->elements[0].connectivity(), nodes));

    auto wrong_parent = parent;
    ++wrong_parent.representation_revision;
    check_acquisition_failure(
        inventory.find(wrong_parent, io::SourceId{300}, rd::ReferenceRole::GmshMesh, {}),
        core::FailureCode::StaleGeneration);
    check_acquisition_failure(
        inventory.find(parent, io::SourceId{300}, rd::ReferenceRole::GmshSidecar, {}),
        core::FailureCode::InvalidInput);
    check_acquisition_failure(
        inventory.find(parent, io::SourceId{999}, rd::ReferenceRole::GmshMesh, {}),
        core::FailureCode::InvalidInput);
}

TEST_CASE("Direct inventory checks cumulative limits and retains the input on any refusal",
          "[IO4-01][N4-02][N4-08][integration][negative]") {
    using enum core::FailureCode;
    using Reason = rd::ReferenceFileFailure;
    const Tree tree;
    const auto text = write_root_inventory(tree.root);
    const auto input = from_memory(tree.root, text);
    auto prepared = prepare(input);
    auto bounds = memory_inventory_bounds;
    auto budget = fixture::budget;
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 request_stop() is non-const.
    std::stop_source cancelled;
    auto expected = ResourceExceeded;
    auto stage = rd::RootSnapshotFailure::Read;
    std::optional<Reason> read;
    std::optional<io::SourceId> location = io::SourceId{800};
    SECTION("source count is decided before reads") {
        --bounds.maximum_sources;
        stage = rd::RootSnapshotFailure::SourceLimit;
        location.reset();
    }
    SECTION("record arrays must be admitted before reads") {
        budget = 0;
        stage = rd::RootSnapshotFailure::InventoryStorage;
        location.reset();
    }
    SECTION("payload bytes include previous files") {
        --bounds.total_payload_bytes;
        read = Reason::FileBytes;
    }
    SECTION("read bytes include an EOF probe per file") {
        --bounds.total_read_bytes;
        read = Reason::ReadBytes;
    }
    SECTION("read call charges include every previous file") {
        --bounds.total_read_calls;
        read = Reason::ReadCalls;
    }
    SECTION("hash bytes include all authenticated bytes") {
        --bounds.total_hash_bytes;
        read = Reason::HashBytes;
    }
    SECTION("identity comparisons are cumulative") {
        --bounds.total_identity_comparisons;
        read = Reason::IdentityComparisons;
    }
    SECTION("the last distinct ID cannot alias an earlier file") {
        REQUIRE(fs::remove(tree.root / "initial.history"));
        fs::create_hard_link(tree.root / "sidecar.json", tree.root / "initial.history");
        expected = InvalidInput;
        read = Reason::DuplicateFile;
    }
    SECTION("the last digest cannot differ even after earlier files verified") {
        write_file(tree.root / "initial.history", "abd");
        expected = InvalidInput;
        read = Reason::DigestMismatch;
    }
    SECTION("a cancelled collection performs no reads or allocation") {
        REQUIRE(cancelled.request_stop());
        expected = Cancelled;
        location.reset();
    }
    core::ResourcePlan memory{budget};
    rd::RootSnapshotDiagnostic diagnostic;
    const auto result = rd::acquire_root_references(input, std::move(prepared), bounds, memory,
                                                    cancelled.get_token(), &diagnostic);
    check_acquisition_failure(result, expected);
    // NOLINTNEXTLINE(bugprone-use-after-move): the rvalue is consumed only on success.
    CHECK(fixture::take(prepared.report()).declarations == 9);
    CHECK(input.view());
    if (expected == Cancelled) {
        CHECK_FALSE(diagnostic.reason);
        CHECK(memory.allocated_bytes() == 0);
    } else {
        CHECK(diagnostic.reason == stage);
    }
    if (location) {
        if (const auto observed = diagnostic.location)
            CHECK(observed->source == location);
        else
            FAIL("Missing original declaration location");
        CHECK(diagnostic.pointer() == "/initial/file");
    }
    if (read) {
        if (const auto observed = diagnostic.read)
            CHECK(observed->reason == read);
        else
            FAIL("Missing same-handle failure");
    }
}

TEST_CASE("File acquisition rejects its manifest and manifest alias before payload admission",
          "[IO4-01][N4-08][integration][negative]") {
    constexpr rd::RootSnapshotLimits no_payload{{0, 0, 0, 0}, 1, 256};
    for (const std::string_view reference : {"manifest.json", "linked.json"}) {
        CAPTURE(reference);
        const Tree tree;
        const auto wire = manifest_with_reference(reference, abc_hex).wire();
        write_file(tree.root / "manifest.json", wire);
        fs::create_hard_link(tree.root / "manifest.json", tree.root / "linked.json");
        core::ResourcePlan creation{fixture::budget};
        const auto input = fixture::take(rd::ReferenceInput::from_file(
            tree.root, "manifest.json", source, manifest_limits, limits(), creation, {}, nullptr));
        const auto prepared = prepare(input);
        core::ResourcePlan plan{fixture::budget};
        rd::RootSnapshotDiagnostic diagnostic;
        const auto acquired =
            rd::acquire_root_snapshot(input, prepared, 0, {}, no_payload, plan, {}, &diagnostic);
        check_acquisition_failure(acquired, core::FailureCode::InvalidInput);
        check_read_reason(diagnostic, rd::ReferenceFileFailure::Cycle);
        CHECK(diagnostic.pointer() == "/model/geometry/0/mesh");
        CHECK(fixture::take(prepared.declaration(0)).original.original_path == reference);
        if (const auto& observed = diagnostic.read) {
            CHECK(observed->bytes_read == 0);
        }
    }
}

TEST_CASE("Memory acquisition owns the selected bytes and rejects only repeated file identities",
          "[IO4-01][N4-02][integration][property]") {
    const auto retained = [] {
        const Tree tree;
        write_file(tree.root / "payload.bin", "abc");
        fs::create_hard_link(tree.root / "payload.bin", tree.root / "alias.bin");
        write_file(tree.root / "copy.bin", "abc");
        write_file(tree.root / "wrong.bin", "abd");
        const auto wire = manifest_with_reference("./payload.bin", abc_hex).wire();
        const auto input = from_memory(tree.root, wire);
        const auto prepared = prepare(input);
        core::ResourcePlan plan{fixture::budget};
        auto first = fixture::take(rd::acquire_root_snapshot(
            input, prepared, 0, {}, acquisition_limits, plan, {}, nullptr));
        CHECK(first.location.source == io::SourceId{300});
        CHECK(first.location.role == rd::ReferenceRole::GmshMesh);
        CHECK(first.location.declaration_node ==
              fixture::take(prepared.declaration(0)).original.declaration_node);
        const auto original = first.file.view();
        if (!original) {
            FAIL("Missing payload");
        } else {
            const std::array seen{original->file_info.identity};
            rd::RootSnapshotDiagnostic diagnostic;
            const auto duplicate = rd::acquire_root_snapshot(
                input, prepared, 1, seen, acquisition_limits, plan, {}, &diagnostic);
            check_acquisition_failure(duplicate, core::FailureCode::InvalidInput);
            check_read_reason(diagnostic, rd::ReferenceFileFailure::DuplicateFile);
            CHECK(diagnostic.pointer() == "/model/geometry/0/sidecar");
            const auto independent = fixture::take(rd::acquire_root_snapshot(
                input, prepared, 2, seen, acquisition_limits, plan, {}, &diagnostic));
            CHECK_FALSE(diagnostic.reason);
            if (const auto other = independent.file.view()) {
                CHECK(other->file_info.identity != original->file_info.identity);
                CHECK(std::ranges::equal(other->bytes, original->bytes));
            } else {
                FAIL("Missing separate equal-content file");
            }

            const auto missing = rd::acquire_root_snapshot(
                input, prepared, 3, seen, acquisition_limits, plan, {}, &diagnostic);
            check_acquisition_failure(missing, core::FailureCode::InvalidInput);
            CHECK(diagnostic.reason == rd::RootSnapshotFailure::Open);
            CHECK(diagnostic.pointer() == "/model/geometry/1/deck");

            auto no_pointer = acquisition_limits;
            no_pointer.diagnostic_bytes = 0;
            const auto wrong = rd::acquire_root_snapshot(input, prepared, 4, seen, no_pointer, plan,
                                                         {}, &diagnostic);
            check_acquisition_failure(wrong, core::FailureCode::InvalidInput);
            check_read_reason(diagnostic, rd::ReferenceFileFailure::DigestMismatch);
            CHECK(diagnostic.pointer().empty());
            const auto without_slot =
                rd::acquire_root_snapshot(input, prepared, 4, seen, no_pointer, plan, {}, nullptr);
            check_acquisition_failure(without_slot, core::FailureCode::InvalidInput);

            const auto spent = plan.allocated_bytes();
            const auto different = from_memory(tree.root, wire, io::SourceId{20});
            const auto stale = rd::acquire_root_snapshot(different, prepared, 0, {},
                                                         acquisition_limits, plan, {}, &diagnostic);
            check_acquisition_failure(stale, core::FailureCode::StaleGeneration);
            CHECK(diagnostic.reason == rd::RootSnapshotFailure::StaleParent);
            const auto invalid = rd::acquire_root_snapshot(
                input, prepared, 8, {}, acquisition_limits, plan, {}, &diagnostic);
            check_acquisition_failure(invalid, core::FailureCode::InvalidInput);
            CHECK(diagnostic.reason == rd::RootSnapshotFailure::InvalidSelection);
            // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop requires mutable source.
            std::stop_source stop;
            REQUIRE(stop.request_stop());
            const auto cancelled = rd::acquire_root_snapshot(
                input, prepared, 0, {}, acquisition_limits, plan, stop.get_token(), &diagnostic);
            check_acquisition_failure(cancelled, core::FailureCode::Cancelled);
            CHECK_FALSE(diagnostic.reason);
            CHECK(plan.allocated_bytes() == spent);
        }
        // Windows also checks that the file handle is already released while the result lives.
        REQUIRE(fs::remove(tree.root / "payload.bin"));
        return first;
    }();
    if (const auto observed = retained.file.view()) {
        CHECK(std::ranges::equal(observed->bytes, bytes_of("abc")));
        CHECK(observed->digest == abc_digest);
        CHECK(observed->file_info.byte_count == 3);
    } else {
        FAIL("Missing detached file snapshot");
    }
    CHECK(retained.pointer() == "/model/geometry/0/mesh");
}

TEST_CASE("Memory root accepts a hard link to a separate manifest as ordinary file bytes",
          "[IO4-01][N4-02][integration][property]") {
    const Tree tree;
    write_file(tree.root / "other-manifest.json", static_case());
    fs::create_hard_link(tree.root / "other-manifest.json", tree.root / "manifest-link.bin");
    const auto wire = manifest_with_reference("manifest-link.bin", valid_manifest_hex).wire();
    const auto input = from_memory(tree.root, wire);
    const auto prepared = prepare(input);
    constexpr rd::RootSnapshotLimits exact{{8611, 8612, 2, 8611}, 0, 0};
    core::ResourcePlan plan{fixture::budget};
    const auto snapshot =
        fixture::take(rd::acquire_root_snapshot(input, prepared, 0, {}, exact, plan, {}, nullptr));
    if (const auto observed = snapshot.file.view()) {
        CHECK(std::ranges::equal(observed->bytes, bytes_of(static_case())));
        CHECK(observed->digest == static_case_digest);
    } else {
        FAIL("Expected separately declared manifest bytes");
    }

    auto too_short = exact;
    --too_short.file.maximum_read_bytes;
    rd::RootSnapshotDiagnostic diagnostic;
    const auto failed =
        rd::acquire_root_snapshot(input, prepared, 0, {}, too_short, plan, {}, &diagnostic);
    check_acquisition_failure(failed, core::FailureCode::ResourceExceeded);
    check_read_reason(diagnostic, rd::ReferenceFileFailure::ReadBytes);
}

TEST_CASE("The first authenticated root enters the walk while later roots stay pending",
          "[IO4-01][N4-02][N4-08][integration][property]") {
    const auto retained = [] {
        constexpr rd::DependencyWalkLimits walk_limits{8, 1, 1, 19, 0};
        const Tree tree;
        const auto wire = manifest_with_reference("first-input.msh", abc_hex).wire();
        write_file(tree.root / "manifest.json", wire);
        write_file(tree.root / "first-input.msh", "abc");
        core::ResourcePlan creation{fixture::budget};
        const auto input = fixture::take(rd::ReferenceInput::from_file(
            tree.root, "manifest.json", source, manifest_limits, limits(), creation, {}, nullptr));
        const auto prepared = prepare(input);
        CHECK(fixture::take(prepared.report()).unique_sources == 8);
        const auto first = fixture::take(prepared.declaration(0));
        core::ResourcePlan plan{fixture::budget};
        rd::ReferenceWalkStartDiagnostic diagnostic;
        auto started = fixture::take(rd::start_reference_walk(
            input, prepared, walk_limits, acquisition_limits, plan, {}, &diagnostic));
        CHECK_FALSE(diagnostic.reason);
        if (const auto state = started.view()) {
            CHECK(fixture::same_identity(state->parent, fixture::take(prepared.identity())));
            if (const auto* pending = state->pending_root) {
                CHECK(pending->edge.declaration_node == first.original.declaration_node);
                CHECK(pending->snapshot.location.declaration_node ==
                      pending->edge.declaration_node);
            } else {
                FAIL("Expected pending dependencies");
            }
        } else {
            FAIL("Expected an owned traversal");
        }

        const auto charges = plan.allocated_bytes();
        REQUIRE(charges > 18 * sizeof(core::Index));
        core::ResourcePlan exact{charges};
        REQUIRE(rd::start_reference_walk(input, prepared, walk_limits, acquisition_limits, exact,
                                         {}, nullptr));
        CHECK(exact.remaining_bytes() == 0);
        core::ResourcePlan short_plan{charges - 1};
        check_acquisition_failure(rd::start_reference_walk(input, prepared, walk_limits,
                                                           acquisition_limits, short_plan, {},
                                                           &diagnostic),
                                  core::FailureCode::ResourceExceeded);
        CHECK(diagnostic.reason == rd::ReferenceWalkStartFailure::Snapshot);
        CHECK(diagnostic.snapshot.reason == rd::RootSnapshotFailure::ProvenanceStorage);
        REQUIRE(fs::remove(tree.root / "first-input.msh"));
        return started;
    }();
    if (const auto state = retained.view()) {
        constexpr std::array<core::Index, 1> active{0};
        CHECK(state->parent.source.source == source);
        CHECK(std::ranges::equal(state->walk.active, active));
        CHECK(state->walk.dependency_first.empty());
        CHECK(state->walk.edges == 1);
        CHECK(state->walk.work == 19);
        if (const auto* pending = state->pending_root) {
            CHECK_FALSE(pending->edge.consumer);
            CHECK(pending->edge.declaring_source == source);
            CHECK(pending->edge.dependency == 0);
            CHECK(pending->edge.role == rd::ReferenceRole::GmshMesh);
            CHECK(pending->snapshot.location.source == io::SourceId{300});
            CHECK(pending->snapshot.location.role == rd::ReferenceRole::GmshMesh);
            CHECK(pending->snapshot.pointer() == "/model/geometry/0/mesh");
            if (const auto payload = pending->snapshot.file.view()) {
                CHECK(std::ranges::equal(payload->bytes, bytes_of("abc")));
                CHECK(payload->digest == abc_digest);
            } else {
                FAIL("Missing authenticated bytes");
            }
        } else {
            FAIL("Expected pending dependencies after dropping the input");
        }
    } else {
        FAIL("Expected a traversal after dropping the input");
    }
}

TEST_CASE("An empty root walk still requires the original checked parent",
          "[IO4-01][N4-02][N4-08][integration][negative]") {
    const auto make_input = [](std::string_view wire, io::SourceId id) {
        core::ResourcePlan plan{fixture::budget};
        return fixture::take(rd::ReferenceInput::from_memory(bytes_of(wire), id, std::nullopt,
                                                             limits(), plan, {}, nullptr));
    };
    const auto input = make_input(static_case(), source);
    const auto other = make_input(static_case(), io::SourceId{20});
    const auto different_bytes = make_input(std::string{static_case()} + " ", source);
    const auto prepared = prepare(input);
    core::ResourcePlan plan{0};
    rd::ReferenceWalkStartDiagnostic diagnostic;
    check_acquisition_failure(
        rd::start_reference_walk(other, prepared, {}, {}, plan, {}, &diagnostic),
        core::FailureCode::StaleGeneration);
    CHECK(diagnostic.reason == rd::ReferenceWalkStartFailure::Binding);
    CHECK(diagnostic.snapshot.reason == rd::RootSnapshotFailure::StaleParent);
    check_acquisition_failure(
        rd::start_reference_walk(different_bytes, prepared, {}, {}, plan, {}, nullptr),
        core::FailureCode::StaleGeneration);
    const auto empty =
        fixture::take(rd::start_reference_walk(input, prepared, {}, {}, plan, {}, &diagnostic));
    CHECK_FALSE(diagnostic.reason);
    CHECK_FALSE(diagnostic.snapshot.reason);
    if (const auto state = empty.view()) {
        CHECK(fixture::same_identity(state->parent, fixture::take(prepared.identity())));
        CHECK(state->pending_root == nullptr);
        CHECK(state->walk.active.empty());
        CHECK(state->walk.dependency_first.empty());
        CHECK(state->walk.edges == 0);
        CHECK(state->walk.work == 0);
    } else {
        FAIL("Missing empty traversal");
    }
    CHECK(plan.allocated_bytes() == 0);
}

TEST_CASE("Only format-proven root leaves enter a closed reference graph",
          "[IO4-01][IO4-04][N4-02][integration][verification][negative]") {
    const auto empty = [] {
        core::ResourcePlan input_plan{fixture::budget};
        const auto input = fixture::take(rd::ReferenceInput::from_memory(
            bytes_of(static_case()), source, std::nullopt, limits(), input_plan, {}, nullptr));
        auto roots = prepare(input);
        core::ResourcePlan no_storage{0};
        auto closed = fixture::take(rd::resolve_leaf_references(input, std::move(roots), {}, {}, {},
                                                                no_storage, {}, nullptr, nullptr));
        CHECK(no_storage.allocated_bytes() == 0);
        return closed;
    }();
    CHECK(fixture::take(empty.identity()).source.source == source);
    CHECK(fixture::take(empty.dependency_count(std::nullopt)) == 0);
    CHECK(empty.dependency_first().empty());
    CHECK_FALSE(empty.snapshot(0));

    const Tree tree;
    const auto wire = write_root_inventory(tree.root);
    write_file(tree.root / "case.json", wire);
    core::ResourcePlan input_plan{fixture::budget};
    const auto input = fixture::take(rd::ReferenceInput::from_file(
        tree.root, "case.json", source, manifest_limits, limits(), input_plan, {}, nullptr));
    auto roots = prepare(input);
    auto acquisition = memory_inventory_bounds;
    acquisition.total_identity_comparisons += inventory_ids.size();
    core::ResourcePlan inventory_plan{fixture::budget};
    constexpr rd::GmshLimits gmsh{1024, 16, 16, 32, 16, 16, 4096};
    constexpr rd::DependencyWalkLimits closure{8, 9, 1, 128, 128};
    core::ResourcePlan graph_plan{fixture::budget};
    rd::LeafReferenceDiagnostic diagnostic;
    rd::RootSnapshotDiagnostic file_diagnostic;
    const auto resolve = [&](rd::PreparedRootReferences& pending, rd::DependencyWalkLimits bounds) {
        return rd::resolve_leaf_references(input, std::move(pending), acquisition, bounds,
                                           {gmsh, {}}, graph_plan, {}, &file_diagnostic,
                                           &diagnostic);
    };
    check_acquisition_failure(resolve(roots, closure), core::FailureCode::Unsupported);
    CHECK_FALSE(file_diagnostic.reason);
    CHECK(diagnostic.reason == rd::LeafReferenceFailure::UnsupportedRole);
    CHECK(diagnostic.pointer() == "/model/geometry/0/sidecar");
    if (const auto location = diagnostic.location) {
        CHECK(location->source == io::SourceId{301});
        CHECK(location->role == rd::ReferenceRole::GmshSidecar);
    } else {
        FAIL("Missing unsupported declaration");
    }
    auto first_edge = closure;
    first_edge.edges = 1;
    check_acquisition_failure(resolve(roots, first_edge), core::FailureCode::ResourceExceeded);
    CHECK(diagnostic.reason == rd::LeafReferenceFailure::Graph);
    CHECK(diagnostic.walk.reason == rd::DependencyWalkFailure::EdgeLimit);
    if (const auto edge = diagnostic.walk.edge) {
        CHECK(edge->declaring_source == source);
        CHECK(edge->dependency == 1);
    } else {
        FAIL("Missing original manifest edge");
    }
    CHECK(diagnostic.pointer() == "/model/geometry/0/sidecar");
    auto inventory = fixture::take(rd::acquire_root_references(input, std::move(roots), acquisition,
                                                               inventory_plan, {}, nullptr));
    CHECK(fixture::take(inventory.source(0)).path.text() == "skin.msh");
    const auto mesh = fixture::take(inventory.snapshot(0));
    CHECK(mesh.get().location.source == io::SourceId{300});
    CHECK_FALSE(inventory.snapshot(-1));
}

TEST_CASE("Root walk applies preflight, snapshot, and edge failures in order",
          "[IO4-01][N4-02][N4-08][integration][negative]") {
    const Tree tree;
    const auto wire = manifest_with_reference("first-input.msh", abc_hex).wire();
    const auto payload_path = tree.root / "first-input.msh";
    const auto input = from_memory(tree.root, wire);
    const auto prepared = prepare(input);
    rd::ReferenceWalkStartDiagnostic diagnostic;
    core::ResourcePlan initial{fixture::budget};
    const auto start = [&](rd::DependencyWalkLimits bounds, core::ResourcePlan& plan,
                           const std::stop_token& stop = {}) {
        return rd::start_reference_walk(input, prepared, bounds, acquisition_limits, plan, stop,
                                        &diagnostic);
    };
    using Code = core::FailureCode;
    using Stage = rd::ReferenceWalkStartFailure;
    check_acquisition_failure(start({7, 1, 1, 19, 0}, initial), Code::ResourceExceeded);
    CHECK(diagnostic.reason == Stage::Sources);
    if (const auto& quota = diagnostic.sources) {
        CHECK(quota->required == 8);
        CHECK(quota->limit == 7);
    } else {
        FAIL("Missing source capacity");
    }
    check_acquisition_failure(start({8, 1, 1, 17, 0}, initial), Code::ResourceExceeded);
    CHECK(diagnostic.reason == Stage::WalkStorage);
    CHECK(initial.allocated_bytes() == 0);
    core::ResourcePlan insufficient{18 * sizeof(core::Index) - 1};
    check_acquisition_failure(start({8, 1, 1, 19, 0}, insufficient), Code::ResourceExceeded);
    CHECK(diagnostic.reason == Stage::WalkStorage);
    CHECK(insufficient.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop requires mutable source.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    check_acquisition_failure(start({8, 1, 1, 19, 0}, initial, cancel.get_token()),
                              Code::Cancelled);
    CHECK_FALSE(diagnostic.reason);
    CHECK_FALSE(diagnostic.snapshot.reason);
    CHECK(initial.allocated_bytes() == 0);

    core::ResourcePlan missing_plan{fixture::budget};
    check_acquisition_failure(start({8, 0, 0, 16, 0}, missing_plan), Code::InvalidInput);
    CHECK(diagnostic.reason == Stage::Snapshot);
    CHECK(diagnostic.snapshot.reason == rd::RootSnapshotFailure::Open);
    CHECK_FALSE(diagnostic.walk.reason);
    CHECK(missing_plan.allocated_bytes() >= 16 * sizeof(core::Index));

    write_file(payload_path, "abd");
    core::ResourcePlan wrong_digest{fixture::budget};
    check_acquisition_failure(start({8, 0, 0, 16, 0}, wrong_digest), Code::InvalidInput);
    CHECK(diagnostic.reason == Stage::Snapshot);
    check_read_reason(diagnostic.snapshot, rd::ReferenceFileFailure::DigestMismatch);
    CHECK(diagnostic.pointer() == "/model/geometry/0/mesh");
    CHECK_FALSE(diagnostic.walk.reason);
    CHECK(wrong_digest.allocated_bytes() > 16 * sizeof(core::Index));

    write_file(payload_path, "abc");
    using WalkReason = rd::DependencyWalkFailure;
    struct Rejection {
        rd::DependencyWalkLimits caps;
        WalkReason reason;
        std::uint64_t limit;
    };
    constexpr std::array cases{
        Rejection{{8, 0, 1, 19, 0}, WalkReason::EdgeLimit, 0},
        Rejection{{8, 1, 1, 18, 0}, WalkReason::Work, 18},
        Rejection{{8, 1, 0, 17, 0}, WalkReason::DepthLimit, 0},
    };
    for (const auto& candidate : cases) {
        CAPTURE(candidate.reason);
        core::ResourcePlan plan{fixture::budget};
        check_acquisition_failure(start(candidate.caps, plan), Code::ResourceExceeded);
        CHECK(diagnostic.reason == Stage::Edge);
        CHECK_FALSE(diagnostic.snapshot.reason);
        CHECK(diagnostic.walk.reason == candidate.reason);
        CHECK(diagnostic.walk.required == candidate.limit + 1);
        CHECK(diagnostic.walk.limit == candidate.limit);
        CHECK(diagnostic.pointer() == "/model/geometry/0/mesh");
        CHECK(plan.allocated_bytes() > 16 * sizeof(core::Index));
    }

    auto tiny_diag = acquisition_limits;
    tiny_diag.diagnostic_bytes = 0;
    core::ResourcePlan compact_plan{fixture::budget};
    check_acquisition_failure(rd::start_reference_walk(input, prepared, cases.front().caps,
                                                       tiny_diag, compact_plan, {}, &diagnostic),
                              Code::ResourceExceeded);
    CHECK(diagnostic.reason == Stage::Edge);
    CHECK(diagnostic.walk.reason == WalkReason::EdgeLimit);
    if (const auto& declaration = diagnostic.walk.edge) {
        CHECK(declaration->declaring_source == source);
        CHECK(declaration->role == rd::ReferenceRole::GmshMesh);
        CHECK(declaration->declaration_node ==
              fixture::take(prepared.declaration(0)).original.declaration_node);
    } else {
        FAIL("Missing declaration with no pointer budget");
    }
    CHECK(diagnostic.pointer().empty());

    core::ResourcePlan final_plan{fixture::budget};
    REQUIRE(start({8, 1, 1, 19, 0}, final_plan));
    CHECK_FALSE(diagnostic.reason);
    CHECK_FALSE(diagnostic.walk.reason);
    CHECK(diagnostic.pointer().empty());
}

namespace {

// The topology is the independent IO4-01 design case: manifest → A,C, A → B,C,
// C → B. Existing fixture roles provide the manifest roots, including a repeated edge.
// This extractor supplies typed declarations and asserts the exact authenticated payload
// it receives. It specifies no product format.
class DesignDependencies final : public rd::ReferenceDependencyExtractor {
public:
    explicit DesignDependencies(core::ResourcePlan& plan) noexcept : plan_(plan) {}

    core::Result<rd::ExtractedReferenceDependencies>
    extract(const rd::AcquiredRootSnapshot& snapshot) noexcept override {
        const auto file = snapshot.file.view();
        if (!file || visits == order.size()) {
            return std::unexpected(
                core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
        }
        const auto id = snapshot.location.source.value();
        order[visits++] = id;
        const auto row = id == 300 ? 2 : id == 302 ? 1 : 0;
        auto entries =
            plan_.allocate<rd::ChildReferenceDeclaration>(row, core::Operation::PrepareModel);
        if (!entries) {
            return std::unexpected(entries.error());
        }
        const rd::ChildReferenceDeclaration b{
            io::SourceId{901}, rd::ReferenceRole::ReferenceGeometry, abc_digest, 0, child_path,
            child_path};
        if (id == 300) {
            sees_exact_bytes = std::ranges::equal(file->bytes, bytes_of(inventory_mesh));
            entries->values()[0] = b;
            entries->values()[1] = {io::SourceId{302}, rd::ReferenceRole::ReferenceGeometry,
                                    abc_digest,        1,
                                    "geometry.json",   "geometry.json"};
        } else if (id == 302) {
            entries->values()[0] = b;
            if (cycle) {
                entries->values()[0] = {io::SourceId{300},
                                        rd::ReferenceRole::GmshMesh,
                                        digest_literal(inventory_mesh_hash),
                                        2,
                                        "skin.msh",
                                        "skin.msh"};
            }
        }
        return rd::ExtractedReferenceDependencies{std::move(*entries), {}};
    }

    std::array<std::uint64_t, 9> order{};
    std::size_t visits{};
    std::string_view child_path{"child/./side"};
    bool cycle{};
    bool sees_exact_bytes{};

private:
    core::ResourcePlan& plan_;
};

rd::ReferenceGraphLimits design_limits() {
    auto inventory = memory_inventory_bounds;
    inventory.maximum_sources += 1;
    inventory.total_payload_bytes += 3;
    inventory.total_hash_bytes += 3;
    inventory.total_read_bytes += 4;
    inventory.total_read_calls += 2;
    inventory.total_identity_comparisons += 8;
    return {inventory, {9, 12, 3, 4096, 256}, 512, {64, 64, 4}};
}

} // namespace

TEST_CASE("An authenticated dependency graph owns sources, provenance and adjacency after cleanup",
          "[IO4-01][N4-02][integration][verification]") {
    const auto graph = [] {
        const Tree tree;
        const auto wire = write_root_inventory(tree.root);
        REQUIRE(fs::create_directory(tree.root / "child"));
        write_file(tree.root / "child" / "side", "abc");
        const auto input = from_memory(tree.root, wire);
        auto prepared = prepare(input);
        core::ResourcePlan plan{fixture::budget};
        DesignDependencies extractor{plan};
        auto result = rd::resolve_reference_graph(input, std::move(prepared), design_limits(),
                                                  extractor, plan, {}, nullptr, nullptr);
        CHECK(extractor.sees_exact_bytes);
        CHECK(extractor.visits == 9);
        CHECK(extractor.order ==
              std::array<std::uint64_t, 9>{300, 901, 302, 301, 200, 201, 700, 900, 800});
        return fixture::take(std::move(result));
    }();

    constexpr std::array<core::Index, 9> expected_order{8, 2, 0, 1, 3, 4, 5, 6, 7};
    CHECK(std::ranges::equal(graph.dependency_first(), expected_order));
    CHECK(fixture::take(graph.dependency_count(std::nullopt)) == 9);
    CHECK(fixture::take(graph.dependency_count(0)) == 2);
    CHECK(fixture::take(graph.dependency_count(2)) == 1);
    CHECK(fixture::take(graph.dependency_count(8)) == 0);
    CHECK_FALSE(graph.dependency(8, 0));
    CHECK_FALSE(graph.dependency(-1, 0));
    CHECK_FALSE(graph.snapshot(9));
    CHECK(fixture::take(graph.dependency(std::nullopt, 7)).topology.dependency == 6);
    CHECK(fixture::take(graph.dependency(std::nullopt, 8)).topology.dependency == 7);
    const auto manifest_edge = fixture::take(graph.dependency({}, 0));
    CHECK(manifest_edge.declaration_pointer == "/model/geometry/0/mesh");
    CHECK(manifest_edge.topology.declaring_source == source);
    CHECK(manifest_edge.original_path == "skin.msh");
    const auto from_a = fixture::take(graph.dependency(0, 0));
    const auto from_c = fixture::take(graph.dependency(2, 0));
    CHECK(from_a.topology.dependency == 8);
    CHECK(from_a.topology.consumer == 0);
    CHECK(from_a.topology.declaring_source == io::SourceId{300});
    CHECK(from_a.topology.declaration_node == 0);
    CHECK(from_a.original_path == "child/./side");
    CHECK(from_a.declaration_pointer.empty());
    CHECK(from_c.topology.dependency == 8);
    CHECK(from_c.topology.declaring_source == io::SourceId{302});
    CHECK(fixture::take(graph.dependency(0, 1)).topology.dependency == 2);

    const auto parent = fixture::take(graph.identity());
    const auto child = fixture::take(
        graph.find(parent, io::SourceId{901}, rd::ReferenceRole::ReferenceGeometry, {}));
    if (const auto child_file = child.get().file.view()) {
        CHECK(std::ranges::equal(child_file->bytes, bytes_of("abc")));
        CHECK(child_file->digest == abc_digest);
    } else {
        FAIL("Missing child payload");
    }
    CHECK(fixture::take(graph.normalized_path(8)) == "child/side");
    CHECK(fixture::take(graph.normalized_path(2)) == "geometry.json");
    CHECK_FALSE(graph.find(parent, io::SourceId{901}, rd::ReferenceRole::GmshSidecar, {}));
    auto stale = parent;
    ++stale.representation_revision;
    check_acquisition_failure(
        graph.find(stale, io::SourceId{901}, rd::ReferenceRole::ReferenceGeometry, {}),
        core::FailureCode::StaleGeneration);
}

TEST_CASE("Reference DFS rejects the first child failure and publishes no partial owner",
          "[IO4-01][N4-02][integration][verification][negative]") {
    using Code = core::FailureCode;
    const Tree tree;
    const auto wire = write_root_inventory(tree.root);
    REQUIRE(fs::create_directory(tree.root / "child"));
    const auto child = tree.root / "child" / "side";
    // Without DFS the missing later manifest file would hide the child failure.
    REQUIRE(fs::remove(tree.root / "sidecar.json"));
    write_file(child, "abd");
    const auto input = from_memory(tree.root, wire);
    auto prepared = prepare(input);
    rd::RootSnapshotDiagnostic acquisition;
    rd::LeafReferenceDiagnostic diagnostic;
    core::ResourcePlan plan{fixture::budget};
    DesignDependencies extractor{plan};
    const auto run = [&](rd::ReferenceGraphLimits caps) {
        extractor.visits = 0;
        return rd::resolve_reference_graph(input, std::move(prepared), caps, extractor, plan, {},
                                           &acquisition, &diagnostic);
    };

    check_acquisition_failure(run(design_limits()), Code::InvalidInput);
    check_read_reason(acquisition, rd::ReferenceFileFailure::DigestMismatch);
    CHECK(diagnostic.declaring_source == io::SourceId{300});
    if (const auto location = acquisition.location) {
        CHECK(location->source == io::SourceId{901});
        CHECK(location->declaration_node == 0);
        CHECK(location->declaring_file == io::SourceId{300});
    } else {
        FAIL("Missing first child location");
    }
    CHECK(extractor.visits == 1);
    CHECK(fixture::take(prepared.report()).unique_sources == 8);

    REQUIRE(fs::remove(child));
    fs::create_hard_link(tree.root / "skin.msh", child);
    check_acquisition_failure(run(design_limits()), Code::InvalidInput);
    check_read_reason(acquisition, rd::ReferenceFileFailure::DuplicateFile);
    CHECK(extractor.visits == 1);
    REQUIRE(fs::remove(child));
    write_file(child, "abc");

    auto too_few_sources = design_limits();
    too_few_sources.inventory.maximum_sources = 8;
    check_acquisition_failure(run(too_few_sources), Code::ResourceExceeded);
    CHECK(diagnostic.reason == rd::LeafReferenceFailure::Graph);
    CHECK_FALSE(acquisition.reason);
    CHECK(extractor.visits == 1);

    auto shallow = design_limits();
    shallow.walk.depth = 1;
    check_acquisition_failure(run(shallow), Code::ResourceExceeded);
    CHECK(diagnostic.walk.reason == rd::DependencyWalkFailure::DepthLimit);
    CHECK(extractor.visits == 1);

    extractor.child_path = "../outside";
    check_acquisition_failure(run(design_limits()), Code::InvalidInput);
    if (const auto path = diagnostic.path) {
        CHECK(path->reason == rd::ReferencePathReason::EscapeRoot);
    } else {
        FAIL("Missing child path reason");
    }
    CHECK(diagnostic.declaring_source == io::SourceId{300});
    CHECK_FALSE(acquisition.reason);
    extractor.child_path = "child/./side";
    extractor.cycle = true;
    const auto cyclic = run(design_limits());
    check_acquisition_failure(cyclic, Code::InvalidInput);
    CHECK(diagnostic.walk.reason == rd::DependencyWalkFailure::Cycle);
    constexpr std::array<core::Index, 3> cycle{0, 2, 0};
    CHECK(std::ranges::equal(diagnostic.walk.cycle.values(), cycle));
    CHECK(diagnostic.declaring_source == io::SourceId{302});
    if (const auto position = diagnostic.location) {
        CHECK(position->declaration_node == 2);
        CHECK(position->declaring_file == io::SourceId{302});
    } else {
        FAIL("Missing child cycle position");
    }
    CHECK_FALSE(acquisition.reason);
    CHECK(fixture::take(prepared.report()).unique_sources == 8);
}

namespace {
constexpr std::string_view owned_deck_source = "GRID,10,,1.25,0.0,0.0\n";
constexpr std::string_view owned_deck_digest =
    "e1b714a50b5bebe1ff4a43cdaaf7749bfed4b197393ad567165ea106584ae09a";
constexpr rd::NastranBulkLimits owned_deck_limits{128, 80, 4, 4, 32, 4096};

// Other roles are test fixtures. The admitted deck runs the production leaf verifier.
class DeckFixtureLeaves final : public rd::ReferenceDependencyExtractor {
public:
    explicit DeckFixtureLeaves(core::ResourcePlan& plan) noexcept : plan_(plan) {}

    core::Result<rd::ExtractedReferenceDependencies>
    extract(const rd::AcquiredRootSnapshot& snapshot) noexcept override {
        if (snapshot.location.role == rd::ReferenceRole::NastranDeck) {
            const auto accepted =
                rd::verify_nastran_bulk_leaf(snapshot, owned_deck_limits, plan_, {});
            if (!accepted) {
                return std::unexpected(accepted.error());
            }
        }
        return rd::ExtractedReferenceDependencies{};
    }

private:
    core::ResourcePlan& plan_;
};
} // namespace

TEST_CASE("A referenced Nastran deck survives native file removal and release of the graph",
          "[IO4-01][IO4-03][integration][negative][verification]") {
    constexpr io::SourceId deck_id{200};
    const auto graph_and_path = [] {
        const Tree tree;
        const auto wire = write_root_inventory(tree.root, owned_deck_source, owned_deck_digest);
        const auto input = from_memory(tree.root, wire);
        auto prepared = prepare(input);
        auto caps = design_limits();
        constexpr auto extra = owned_deck_source.size() - 3;
        caps.inventory.total_payload_bytes += extra;
        caps.inventory.total_read_bytes += extra;
        caps.inventory.total_hash_bytes += extra;
        core::ResourcePlan plan{fixture::budget};
        DeckFixtureLeaves extractor{plan};
        auto graph = fixture::take(rd::resolve_reference_graph(
            input, std::move(prepared), caps, extractor, plan, {}, nullptr, nullptr));
        return std::pair{std::move(graph), tree.root / "model.bdf"};
    };

    auto deck = [&] {
        const auto [graph, old_path] = graph_and_path();
        REQUIRE_FALSE(fs::exists(old_path));
        const auto identity = fixture::take(graph.identity());
        const auto original =
            fixture::take(graph.find(identity, deck_id, rd::ReferenceRole::NastranDeck, {}));
        const auto authenticated = original.get().file.view();
        if (!authenticated) {
            FAIL("Missing authenticated deck bytes");
        }
        CHECK(authenticated->digest == digest_literal(owned_deck_digest));

        auto other = identity;
        ++other.representation_revision;
        core::ResourcePlan plan{fixture::budget};
        const auto stale =
            rd::prepare_referenced_nastran_deck(graph, other, deck_id, owned_deck_limits, plan, {});
        check_acquisition_failure(stale, core::FailureCode::StaleGeneration);
        const auto wrong_role = rd::prepare_referenced_nastran_deck(
            graph, identity, io::SourceId{201}, owned_deck_limits, plan, {});
        check_acquisition_failure(wrong_role, core::FailureCode::InvalidInput);
        CHECK(plan.allocated_bytes() == 0);

        const auto acquired = rd::prepare_referenced_nastran_deck(graph, identity, deck_id,
                                                                  owned_deck_limits, plan, {});
        REQUIRE(acquired);
        if (const auto view = acquired->view()) {
            CHECK(view->bytes == owned_deck_source);
            CHECK(view->digest == authenticated->digest);
            CHECK(view->location.declaration_node == original.get().location.declaration_node);
            CHECK(view->declaration_pointer == "/model/geometry/1/deck");
            CHECK(view->declaration_pointer == original.get().pointer());
            CHECK(view->declaration_pointer.data() != original.get().pointer().data());
            CHECK(bytes_of(view->bytes).data() != authenticated->bytes.data());
        } else {
            FAIL("Missing owned Nastran bytes");
        }

        const auto required = plan.allocated_bytes();
        REQUIRE(required > owned_deck_source.size());
        core::ResourcePlan short_plan{required - 1};
        std::optional<rd::NastranBulkDiagnostic> diagnostic;
        const auto short_copy = rd::prepare_referenced_nastran_deck(
            graph, identity, deck_id, owned_deck_limits, short_plan, {}, &diagnostic);
        check_acquisition_failure(short_copy, core::FailureCode::ResourceExceeded);
        if (diagnostic) {
            CHECK(diagnostic->reason == rd::NastranBulkFailure::Limit);
        } else {
            FAIL("Missing retained-byte admission diagnostic");
        }
        CHECK(short_plan.allocated_bytes() == 0);

        core::ResourcePlan exact_plan{required};
        auto result = fixture::take(rd::prepare_referenced_nastran_deck(
            graph, identity, deck_id, owned_deck_limits, exact_plan, {}, &diagnostic));
        CHECK(exact_plan.remaining_bytes() == 0);
        CHECK_FALSE(diagnostic);
        const auto still_authenticated = original.get().file.view();
        if (still_authenticated) {
            CHECK(std::ranges::equal(still_authenticated->bytes, bytes_of(owned_deck_source)));
        } else {
            FAIL("Import changed its parent snapshot");
        }
        return result;
    }();

    core::ResourcePlan grid_plan{fixture::budget};
    const auto coordinates = fixture::take(
        rd::NastranCoordinateCards::prepare(std::move(deck), {4, 4, 4, 4096}, grid_plan, {}));
    if (const auto values = coordinates.view()) {
        REQUIRE(values->grids.size() == 1);
        CHECK(values->grids.front().id == rd::NastranGridId{10});
        CHECK(values->grids.front().position[0].decimal() == "1.25");
        CHECK(values->bulk.bytes == owned_deck_source);
        CHECK(values->bulk.declaration_pointer == "/model/geometry/1/deck");
    } else {
        FAIL("The card owner did not survive graph release");
    }
}

TEST_CASE("Authenticated mesh sidecars bind original faces to existing Case PGs",
          "[FD-03][MAT-D05][IO4-01][IO4-03][IO4-04][R-05][integration][verification][negative]") {
    constexpr std::string_view mesh = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 1 2 0
11 0 0 0 1000 0 0 1 7 0
17 0 0 0 1000 1000 0 1 7 0
18 0 0 0 1000 1000 0 1 8 0
$EndEntities
$Nodes
1 4 10 99
2 17 0 4
10
20
40
99
0 0 0
1000 0 0
1000 1000 0
0 1000 0
$EndNodes
$Elements
3 4 501 999
1 11 1 1
501 10 20
2 17 2 2
702 10 40 99
701 10 20 40
2 18 2 1
999 10 20 99
$EndElements
)MSH";
    constexpr std::string_view bulk =
        "GRID,10,,0.0,0.0,0.0\nGRID,20,,1.0,0.0,0.0\nGRID,40,,1.0,1.0,0.0\n"
        "GRID,99,,0.0,1.0,0.0\nPSHELL,90,6,0.25\nMAT1,6,1.0,,0.0,1.0\n"
        "CTRIA3,701,90,10,20,40\nCTRIA3,702,90,10,40,99\n"
        "SPC,42,10,1,0.0\nSPC1,51,36,10,40\nSPC1,61,5,10\nSPCADD,61,42,51\n";
    const auto zero = nlohmann::json{{"numerator", "0"}, {"denominator", "1"}};
    const auto one = nlohmann::json{{"numerator", "1"}, {"denominator", "1"}};
    const auto chart = [&](bool upper) {
        nlohmann::json result{{"orientation", "same"}};
        result["vertices"] = nlohmann::json::array(
            {nlohmann::json::array({zero, zero}),
             upper ? nlohmann::json::array({one, one}) : nlohmann::json::array({one, zero}),
             upper ? nlohmann::json::array({zero, one}) : nlohmann::json::array({one, one})});
        return result;
    };

    bool is_gmsh = false;
    bool omitted = false;
    bool duplicate = false;
    bool lexical = false;
    bool wrong_case = false;
    bool foreign_element = false;
    bool wrong_components = false;
    bool budget_repeat = false;
    bool wrong_dimension = false;
    bool wrong_chart = false;
    bool version = false;
    bool physical = false;
    bool thickness_shift = false;
    bool material_shift = false;
    bool chart_shift = false;
    bool frame_shift = false;
    bool conversion = false;
    bool distance_short = false;
    bool normal_short = false;
    bool pcomp = false;
    enum class PcompRegion : std::uint8_t { None, Paired, WrongTotal, WrongInterface, Mixed };
    PcompRegion pcomp_region = PcompRegion::None;
    bool reversed = false;
    bool wrong_layer_material = false;
    bool partial_group = false;
    bool shared_group_thickness = false;
    bool contradictory_group_thickness = false;
    bool partial_pshell = false;
    bool wrong_pshell_thickness = false;
    enum class Corners : std::uint8_t {
        None,
        Absolute,
        Fraction,
        WrongFraction,
        MissingBase,
        MixedMaterial
    };
    Corners corners = Corners::None;
    enum class Axes : std::uint8_t {
        None,
        Quarter,
        CancelledAngles,
        Reversed,
        Mcid,
        Difference,
        Degenerate,
        ResidualAngle,
        Thirty,
        Sixty,
        WrongThirtySign,
        FortyFive,
        WrongFortyFiveQuadrant,
        VariableThickness
    };
    Axes axes = Axes::None;
    bool shell_set = false;
    bool late_axis_mismatch = false;
    bool native_cards = false;
    bool thermal_link = false;
    SECTION("Gmsh group expands and preserves the original rings") {
        is_gmsh = true;
    }
    SECTION("One Gmsh group can be partitioned, but the uncovered member is rejected") {
        is_gmsh = true;
        omitted = true;
    }
    SECTION("Gmsh duplicate ownership is rejected in the declaring row") {
        is_gmsh = true;
        duplicate = true;
    }
    SECTION("A face in another Gmsh physical group cannot bind") {
        is_gmsh = true;
        foreign_element = true;
    }
    SECTION("The adopted schema cannot assert a dependency-free leaf with an unknown field") {
        is_gmsh = true;
        lexical = true;
    }
    SECTION("An unknown revision of the correct sidecar family is unsupported") {
        is_gmsh = true;
        lexical = true;
        version = true;
    }
    SECTION("An exact chart denominator cannot be zero") {
        is_gmsh = true;
        lexical = true;
        wrong_chart = true;
    }
    SECTION("A physical tag cannot be resolved in a different dimension") {
        is_gmsh = true;
        wrong_dimension = true;
    }
    SECTION("Case fit and sidecar conversion cannot both supply a budget") {
        is_gmsh = true;
        budget_repeat = true;
    }
    SECTION("Nastran resolves deck EIDs, MID and the expanded SPCADD component set") {}
    SECTION("Matching an incomplete Nastran support component set fails") {
        wrong_components = true;
    }
    SECTION("Gmsh admitted planes compare 250 mm with one quarter metre") {
        is_gmsh = true;
        physical = true;
    }
    SECTION("Gmsh exact thickness detects a sub-ulp difference") {
        is_gmsh = true;
        physical = true;
        thickness_shift = true;
    }
    SECTION("Gmsh mismatched chart fails before its different thickness") {
        is_gmsh = true;
        physical = true;
        thickness_shift = true;
        chart_shift = true;
    }
    SECTION("A translated Case frame cannot qualify the unchanged Gmsh source") {
        is_gmsh = true;
        physical = true;
        frame_shift = true;
    }
    SECTION("Nastran admitted planes compare PSHELL thickness and explicit MAT1 quantities") {
        physical = true;
    }
    SECTION("Nastran MAT1 mismatch reports the original deck field") {
        physical = true;
        material_shift = true;
    }
    SECTION("A Gmsh fit budget admits a curved PG above the original planar faces") {
        is_gmsh = true;
        physical = true;
        conversion = true;
    }
    SECTION("A Gmsh fit budget below the independent interior difference rejects the face") {
        is_gmsh = true;
        physical = true;
        conversion = true;
        distance_short = true;
    }
    SECTION("A Nastran sidecar budget admits the same curved PG in Case units") {
        physical = true;
        conversion = true;
    }
    SECTION("A Nastran normal budget below the true slope difference rejects the face") {
        physical = true;
        conversion = true;
        normal_short = true;
    }
    SECTION("PCOMP compares distinct MAT1 and MAT8 plies in Case order") {
        physical = true;
        pcomp = true;
    }
    SECTION("Opposite Nastran rings reverse the explicit materials and thickness intervals") {
        physical = true;
        pcomp = true;
        reversed = true;
    }
    SECTION("A reversed stack still detects a changed MAT8 constitutive value") {
        physical = true;
        pcomp = true;
        reversed = true;
        material_shift = true;
    }
    SECTION("A source MID mapped to the wrong Case layer is rejected") {
        physical = true;
        pcomp = true;
        reversed = true;
        wrong_layer_material = true;
    }
    SECTION("A Gmsh group's declared thickness excludes a neighboring unclaimed group") {
        is_gmsh = true;
        physical = true;
        partial_group = true;
    }
    SECTION("A second row of the same Gmsh group inherits the declared thickness") {
        is_gmsh = true;
        physical = true;
        shared_group_thickness = true;
    }
    SECTION("Two conflicting thickness declarations for one Gmsh group fail") {
        is_gmsh = true;
        contradictory_group_thickness = true;
    }
    SECTION("A PSHELL with no explicit T leaves its own adjacent region unclaimed") {
        physical = true;
        partial_pshell = true;
    }
    SECTION("A distinct PSHELL T is checked on its own positive-area element") {
        physical = true;
        partial_pshell = true;
        wrong_pshell_thickness = true;
    }
    SECTION("Absolute triangle corner thickness and omitted corners bind to the Case") {
        physical = true;
        corners = Corners::Absolute;
    }
    SECTION("Native quadrilateral thickness factors own the value and the scale position") {
        physical = true;
        corners = Corners::Fraction;
    }
    SECTION("A wrong native corner fraction is distinguished from a uniform property") {
        physical = true;
        corners = Corners::WrongFraction;
    }
    SECTION("Fractional corner thickness requires the referenced PSHELL thickness") {
        physical = true;
        corners = Corners::MissingBase;
    }
    SECTION("Separate membrane and bending materials have no Case total thickness match") {
        physical = true;
        corners = Corners::MixedMaterial;
    }
    SECTION("Element orientation and ply angle produce the same Case quarter turn") {
        physical = true;
        axes = Axes::Quarter;
    }
    SECTION("General source and Case rotations with the same exact angle cancel symbolically") {
        physical = true;
        axes = Axes::CancelledAngles;
    }
    SECTION("Reversed native winding reverses the represented angle about the Case normal") {
        physical = true;
        axes = Axes::Reversed;
    }
    SECTION("A translated rectangular MCID projects onto the Case explicit triad") {
        physical = true;
        axes = Axes::Mcid;
    }
    SECTION("Valid distinct material axes are rejected even with identical thickness") {
        physical = true;
        axes = Axes::Difference;
    }
    SECTION("A Case direction along the layer normal has a zero projection") {
        physical = true;
        axes = Axes::Degenerate;
    }
    SECTION("A rational direction cannot realize a fifteen-degree difference here") {
        physical = true;
        axes = Axes::ResidualAngle;
    }
    SECTION("Tilted rational directions realize an exact negative thirty degrees") {
        physical = true;
        axes = Axes::Thirty;
    }
    SECTION("Tilted rational tangents also realize an exact negative sixty degrees") {
        physical = true;
        axes = Axes::Sixty;
    }
    SECTION("Thirty degrees with the opposite sine sign cannot pass") {
        physical = true;
        axes = Axes::WrongThirtySign;
    }
    SECTION("Rational planar diagonals make an exact forty-five degrees") {
        physical = true;
        axes = Axes::FortyFive;
    }
    SECTION("The same squared cosine in the wrong quadrant cannot pass") {
        physical = true;
        axes = Axes::WrongFortyFiveQuadrant;
    }
    SECTION("Varying membrane thickness cannot certify one full-depth layer normal") {
        physical = true;
        axes = Axes::VariableThickness;
        corners = Corners::Fraction;
    }

    SECTION("One deck prepares every declared Shell in first occurrence order") {
        physical = true;
        shell_set = true;
    }
    SECTION("A later Shell mismatch leaves both inputs available for the same preparation") {
        physical = true;
        shell_set = true;
        late_axis_mismatch = true;
    }

    SECTION("One authenticated input owns its Shell, beam, mass and rigid reference evidence") {
        shell_set = true;
        physical = true;
        native_cards = true;
    }
    SECTION("The late native RBE2 error keeps its source field and Case declaration") {
        shell_set = true;
        physical = true;
        native_cards = true;
        thermal_link = true;
    }

    const bool corner_claim = corners != Corners::None;
    const bool axis_pcomp = axes != Axes::None && axes != Axes::VariableThickness;
    const bool two_axis_layers = axes == Axes::Quarter || axes == Axes::Difference;
    const bool tilted_axes =
        axes == Axes::Thirty || axes == Axes::Sixty || axes == Axes::WrongThirtySign;
    const bool forty_five_axes = axes == Axes::FortyFive || axes == Axes::WrongFortyFiveQuadrant;
    const bool quad = axis_pcomp || (corner_claim && corners != Corners::Absolute);
    SECTION("Distinct PCOMP properties with the same physical layers cover a Case shell") {
        physical = true;
        pcomp = true;
        pcomp_region = PcompRegion::Paired;
    }
    SECTION("PCOMP full thickness disagreement is attached to its own face and property") {
        physical = true;
        pcomp = true;
        pcomp_region = PcompRegion::WrongTotal;
    }
    SECTION("Equal local PCOMP totals do not excuse an incompatible layer interface") {
        physical = true;
        pcomp = true;
        pcomp_region = PcompRegion::WrongInterface;
    }
    SECTION("PCOMP leaves the adjacent undeclared PSHELL thickness unconstrained") {
        physical = true;
        pcomp = true;
        pcomp_region = PcompRegion::Mixed;
    }
    const bool separate_pcomp =
        pcomp_region != PcompRegion::None && pcomp_region != PcompRegion::Mixed;
    SECTION("Another checked Case cannot consume the same graph") {
        wrong_case = true;
    }

    const auto empty = nlohmann::json::array();
    nlohmann::json sidecar = {
        {"schema", is_gmsh ? "uvlas.gmsh-sidecar/1" : "uvlas.nastran-sidecar/1"},
        {"geometry", "2"},
        {"shells", empty},
        {"supports", empty},
        {"joints", empty}};
    if (is_gmsh) {
        sidecar["length"] = "mm";
        sidecar["shells"].push_back(
            {{"physical", {{"dimension", 2}, {"tag", 7}}},
             {"shell", "1"},
             {"faces", nlohmann::json::array({{{"element", "702"}, {"chart", chart(true)}}})},
             {"material_layers", nlohmann::json::array({{{"layer", "1"}, {"material", "1"}}})},
             {"full_thickness", 250}});
        if (!omitted) {
            auto row = sidecar["shells"][0];
            // Reverse MSH and declaration order for the same selected group.
            row["faces"] = nlohmann::json::array({{{"element", duplicate         ? "702"
                                                               : foreign_element ? "999"
                                                                                 : "701"},
                                                   {"chart", chart(false)}}});
            if (partial_group) {
                row["physical"]["tag"] = 9;
            }
            if (partial_group || shared_group_thickness) {
                row.erase("full_thickness");
            }
            if (contradictory_group_thickness) {
                row["full_thickness"] = 251;
            }
            sidecar["shells"].push_back(std::move(row));
        }
        if (chart_shift) {
            sidecar["shells"][0]["faces"][0]["chart"]["vertices"][1][1] = {{"numerator", "1"},
                                                                           {"denominator", "2"}};
        }
        sidecar["supports"].push_back(
            {{"physical", {{"dimension", wrong_dimension ? 0 : 1}, {"tag", 7}}}, {"support", "1"}});
        if (version) {
            sidecar["schema"] = "uvlas.gmsh-sidecar/2";
        } else if (wrong_chart) {
            sidecar["shells"][0]["faces"][0]["chart"]["vertices"][1][0]["denominator"] = "0";
        } else if (lexical) {
            sidecar["include"] = "unread.json";
        }
        if (budget_repeat) {
            sidecar["conversion"] = {{"distance_budget", 0.1},
                                     {"normal_budget", 0.1},
                                     {"conversion_observations", {"1"}}};
        }
    } else {
        sidecar["generator"] = "test-deck";
        sidecar["units"] = {{"length", "m"}, {"mass", "kg"}, {"time", "s"}, {"angle", "deg"}};
        sidecar["shells"].push_back({{"element", "702"}, {"shell", "1"}, {"chart", chart(true)}});
        sidecar["shells"].push_back({{"element", "701"}, {"shell", "1"}, {"chart", chart(false)}});
        if (quad) {
            sidecar["shells"].erase(0);
            auto& vertices = sidecar["shells"][0]["chart"]["vertices"];
            vertices.push_back(nlohmann::json::array({zero, one}));
            if (axes == Axes::Reversed) {
                vertices[1].swap(vertices[3]);
                sidecar["shells"][0]["chart"]["orientation"] = "reversed";
            }
        }
        if (shell_set) {
            auto square = chart(false);
            square["vertices"].push_back(nlohmann::json::array({zero, one}));
            sidecar["shells"] = nlohmann::json::array({
                {{"element", "900"}, {"shell", "2"}, {"chart", chart(true)}},
                {{"element", "701"}, {"shell", "1"}, {"chart", std::move(square)}},
                {{"element", "901"}, {"shell", "2"}, {"chart", chart(false)}},
            });
        }
        sidecar["materials"] = nlohmann::json::array({{{"mid", "6"}, {"material", "1"}}});
        if (corners == Corners::MixedMaterial) {
            sidecar["materials"].push_back({{"mid", "9"}, {"material", "1"}});
        }
        if (pcomp && pcomp_region != PcompRegion::Mixed) {
            sidecar["materials"].push_back(
                {{"mid", "9"}, {"material", wrong_layer_material ? "1" : "2"}});
        }
        if (reversed) {
            for (auto& row : sidecar["shells"]) {
                auto& mapping = row["chart"];
                mapping["orientation"] = "reversed";
                mapping["vertices"][1].swap(mapping["vertices"][2]);
            }
        }
        sidecar["beams"] = empty;
        sidecar["masses"] = empty;
        if (native_cards) {
            sidecar["beams"].push_back({{"element", "905"}, {"beam", "1"}});
            sidecar["masses"].push_back({{"element", "906"}, {"mass", "1"}});
            sidecar["joints"].push_back({{"element", "907"}, {"joint", "1"}});
        }
        sidecar["supports"] = nlohmann::json::array(
            {{{"set", "61"},
              {"grid", "10"},
              {"support", "1"},
              {"components", wrong_components ? nlohmann::json::array({1, 3})
                                              : nlohmann::json::array({1, 3, 6})}}});
        if (conversion) {
            sidecar["conversion"] = {{"distance_budget", 26},
                                     {"normal_budget", normal_short ? 0.01 : 1.0},
                                     {"conversion_observations", nlohmann::json::array({"1"})}};
        }
    }

    const auto hex = [](std::string_view body) {
        constexpr std::string_view digits = "0123456789abcdef";
        const auto digest = rd::source_digest(bytes_of(body));
        std::string text;
        for (const auto byte : digest.bytes) {
            const auto value = std::to_integer<unsigned>(byte);
            text.push_back(digits[value >> 4U]);
            text.push_back(digits[value & 15U]);
        }
        return text;
    };
    auto sidecar_text = sidecar.dump();
    if (thickness_shift) {
        constexpr std::string_view shared = R"("full_thickness":250)";
        constexpr std::string_view different = R"("full_thickness":250.0000000000000000001)";
        std::size_t first_position = 0;
        for (int i = 0; i != 2; ++i) {
            const auto offset = sidecar_text.find(shared, first_position);
            REQUIRE(offset != std::string::npos);
            sidecar_text.replace(offset, shared.size(), different);
            first_position = offset + different.size();
        }
    }
    std::string source_text{is_gmsh ? mesh : bulk};
    if (partial_group) {
        const auto replace_group = [&](std::string_view old_text, std::string_view new_text) {
            const auto offset = source_text.find(old_text);
            REQUIRE(offset != std::string::npos);
            source_text.replace(offset, old_text.size(), new_text);
        };
        replace_group("$Entities\n0 1 2 0\n", "$Entities\n0 1 3 0\n");
        replace_group("$EndEntities\n", "19 0 0 0 1000 1000 0 1 9 0\n$EndEntities\n");
        replace_group("$Elements\n3 4 501 999\n", "$Elements\n4 4 501 999\n");
        replace_group("2 17 2 2\n702 10 40 99\n701 10 20 40\n",
                      "2 17 2 1\n702 10 40 99\n2 19 2 1\n701 10 20 40\n");
    }
    if (shell_set) {
        constexpr std::string_view original = "PSHELL,90,6,0.25\n";
        const auto property = source_text.find(original);
        REQUIRE(property != std::string::npos);
        source_text.replace(property, original.size(), "PCOMP,90,-.125,,,,,,\n,6,.25\n");
        constexpr std::string_view triangles = "CTRIA3,701,90,10,20,40\nCTRIA3,702,90,10,40,99\n";
        const auto face = source_text.find(triangles);
        REQUIRE(face != std::string::npos);
        source_text.replace(face, triangles.size(),
                            "GRID,110,,0.0,0.0,1.0\nGRID,120,,1.0,0.0,1.0\n"
                            "GRID,140,,1.0,1.0,1.0\nGRID,199,,0.0,1.0,1.0\n"
                            "CTRIA3,900,90,110,140,199,0\nCQUAD4,701,90,10,20,40,99,0\n"
                            "CTRIA3,901,90,110,120,140,0\n");
    }
    if (native_cards) {
        source_text += "CBEAM,905,95,10,20,99\nPBEAM,95,6,2.0,3.0,4.0,,5.0\n"
                       "CONM2,906,10,0,2.0\n";
        source_text += thermal_link ? "RBE2,907,10,123,20,0.001,293.15\n" : "RBE2,907,10,123,20\n";
    }
    if (axis_pcomp) {
        constexpr std::string_view pshell = "PSHELL,90,6,0.25\n";
        const auto property = source_text.find(pshell);
        REQUIRE(property != std::string::npos);
        if (tilted_axes) {
            const std::string_view layer =
                axes == Axes::Thirty  ? "PCOMP,90,-.125,,,,,,\n,6,.25,-30.0\n"
                : axes == Axes::Sixty ? "PCOMP,90,-.125,,,,,,\n,6,.25,-60.0\n"
                                      : "PCOMP,90,-.125,,,,,,\n,6,.25,30.0\n";
            source_text.replace(property, pshell.size(), layer);
            source_text += "CORD2R,12,,1.0,2.0,3.0,1.0,2.0,8.0\n,2.0,3.0,3.0\n";
            const auto world_grid = [&](std::string_view old_grid, std::string_view new_grid) {
                const auto index = source_text.find(old_grid);
                REQUIRE(index != std::string::npos);
                source_text.replace(index, old_grid.size(), new_grid);
            };
            world_grid("GRID,20,,1.0,0.0,0.0\n", "GRID,20,,1.0,0.0,1.0\n");
            world_grid("GRID,99,,0.0,1.0,0.0\n", "GRID,99,,0.0,1.0,-1.0\n");
        } else if (forty_five_axes) {
            source_text.replace(property, pshell.size(), "PCOMP,90,-.125,,,,,,\n,6,.25,30.0\n");
        } else if (two_axis_layers) {
            source_text.replace(property, pshell.size(),
                                "PCOMP,90,-.125,,,,,,\n,6,.125,60.0,,6,.125,-30.0\n");
        } else if (axes == Axes::CancelledAngles) {
            source_text.replace(property, pshell.size(), "PCOMP,90,-.125,,,,,,\n,6,.25,10.0\n");
        } else if (axes == Axes::Mcid) {
            source_text.replace(property, pshell.size(),
                                "PCOMP,90,-.125,,,,,,\n,6,.25\n"
                                "CORD2R,12,,1.0,2.0,3.0,1.0,2.0,8.0\n,1.0,5.0,3.0\n");
        } else {
            source_text.replace(property, pshell.size(), "PCOMP,90,-.125,,,,,,\n,6,.25,60.0\n");
        }
        constexpr std::string_view triangles = "CTRIA3,701,90,10,20,40\nCTRIA3,702,90,10,40,99\n";
        const auto position = source_text.find(triangles);
        REQUIRE(position != std::string::npos);
        const std::string_view face =
            axes == Axes::Reversed                 ? "CQUAD4,701,90,10,99,40,20,30.0,0.0\n"
            : (axes == Axes::Mcid || tilted_axes)  ? "CQUAD4,701,90,10,20,40,99,12,0.0\n"
            : axes == Axes::CancelledAngles        ? "CQUAD4,701,90,10,20,40,99,20.0,0.0\n"
            : axes == Axes::FortyFive              ? "CQUAD4,701,90,10,20,40,99,15.0,0.0\n"
            : axes == Axes::WrongFortyFiveQuadrant ? "CQUAD4,701,90,10,20,40,99,105.0,0.0\n"
                                                   : "CQUAD4,701,90,10,20,40,99,30.0,0.0\n";
        source_text.replace(position, triangles.size(), face);
    }
    if (corner_claim) {
        constexpr std::string_view original = "PSHELL,90,6,0.25\n";
        const auto position = source_text.find(original);
        REQUIRE(position != std::string::npos);
        const std::string_view property = corners == Corners::MissingBase ? "PSHELL,90,6,,6,,6\n"
                                          : corners == Corners::MixedMaterial
                                              ? "PSHELL,90,6,0.25,9,,6\nMAT1,9,1.0,,0.0,1.0\n"
                                              : "PSHELL,90,6,0.25,6,,6\n";
        source_text.replace(position, original.size(), property);
        constexpr std::string_view triangles = "CTRIA3,701,90,10,20,40\nCTRIA3,702,90,10,40,99\n";
        const auto face_position = source_text.find(triangles);
        REQUIRE(face_position != std::string::npos);
        const std::string_view faces = !quad ? "CTRIA3,701,90,10,20,40,,,\n,,0,,0.5,0.5\n"
                                               "CTRIA3,702,90,10,40,99,,,\n,,0,,0.5,\n"
                                       : corners == Corners::WrongFraction
                                           ? "CQUAD4,701,90,10,20,40,99,0.0,0.0\n,,1,,2.1,2.0\n"
                                           : "CQUAD4,701,90,10,20,40,99,0.0,0.0\n,,1,,2.0,2.0\n";
        source_text.replace(face_position, triangles.size(), faces);
    }
    if (partial_pshell) {
        constexpr std::string_view original = "CTRIA3,701,90";
        const auto offset = source_text.find(original);
        REQUIRE(offset != std::string::npos);
        source_text.replace(offset, original.size(), "CTRIA3,701,92");
        source_text += wrong_pshell_thickness ? "PSHELL,92,6,0.5\n" : "PSHELL,92,6\n";
    }
    if (pcomp) {
        constexpr std::string_view pshell = "PSHELL,90,6,0.25\n";
        const auto offset = source_text.find(pshell);
        REQUIRE(offset != std::string::npos);
        if (pcomp_region == PcompRegion::Mixed) {
            source_text.replace(offset, pshell.size(), "PCOMP,90,-.125,,,,,,\n,6,.25\n");
        } else {
            source_text.replace(offset, pshell.size(),
                                "PCOMP,90,-.125,,,,,,\n,6,.0625,,,9,.1875\n"
                                "MAT8,9,130.0,9.0,.28\n");
        }
        if (pcomp_region != PcompRegion::None) {
            constexpr std::string_view first = "CTRIA3,701,90";
            const auto position = source_text.find(first);
            REQUIRE(position != std::string::npos);
            source_text.replace(position, first.size(), "CTRIA3,701,92");
            switch (pcomp_region) {
            case PcompRegion::Paired:
                source_text += "PCOMP,92,-.125,,,,,,\n,6,.0625,,,9,.1875\n";
                break;
            case PcompRegion::WrongTotal:
                source_text += "PCOMP,92,-.25,,,,,,\n,6,.125,,,9,.375\n";
                break;
            case PcompRegion::WrongInterface:
                source_text += "PCOMP,92,-.125,,,,,,\n,6,.125,,,9,.125\n";
                break;
            case PcompRegion::Mixed:
                source_text += "PSHELL,92,6\n";
                break;
            case PcompRegion::None:
                break;
            }
        }
    }
    if (reversed) {
        constexpr std::array originals{std::string_view{"701,90,10,20,40"},
                                       std::string_view{"702,90,10,40,99"}};
        constexpr std::array opposites{std::string_view{"701,90,10,40,20"},
                                       std::string_view{"702,90,10,99,40"}};
        for (std::size_t i = 0; i < originals.size(); ++i) {
            const auto offset = source_text.find(originals[i]);
            REQUIRE(offset != std::string::npos);
            source_text.replace(offset, originals[i].size(), opposites[i]);
        }
    }
    if (material_shift) {
        const std::string_view original = pcomp ? "MAT8,9,130.0" : "MAT1,6,1.0";
        const auto offset = source_text.find(original);
        REQUIRE(offset != std::string::npos);
        source_text.replace(offset, original.size(),
                            pcomp ? "MAT8,9,130.0000000000000000001"
                                  : "MAT1,6,1.0000000000000000001");
    }
    auto config = fixture::root("coupled_transient");
    config.field("model")
        .field("geometry")
        .elements[0]
        .field("thickness")
        .set("values", fixture::json("[0.25,0.25,0.25,0.25]"));
    if (tilted_axes) {
        auto& surface = config.field("model").field("geometry").elements[0];
        surface.set("degree", fixture::json("[1,1]"));
        surface.set(
            "coefficients",
            fixture::json(
                R"([{"powers":[0,1],"coefficient":{"numerator":"-1","denominator":"1"}},{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}])"));
    }
    if (axes != Axes::None) {
        config.field("units").set("angle", fixture::json(R"("deg")"));
        auto& layers = config.field("model").field("shells").elements[0].field("layers").elements;
        auto& layer = layers.front();
        if (axes == Axes::Mcid) {
            layer.set(
                "axis",
                fixture::json(
                    R"({"kind":"explicit","columns":{"kind":"constant","value":[0,-1,0,1,0,0,0,0,1]}})"));
        } else {
            const std::string_view direction =
                axes == Axes::Degenerate                               ? "[0,0,1]"
                : axes == Axes::Sixty                                  ? "[1,0,1]"
                : tilted_axes                                          ? "[2,1,1]"
                : forty_five_axes                                      ? "[1,1,0]"
                : (two_axis_layers || axes == Axes::VariableThickness) ? "[0,1,2]"
                                                                       : "[1,0,2]";
            const std::string_view angle = axes == Axes::CancelledAngles ? "[30]"
                                           : axes == Axes::ResidualAngle ? "[15]"
                                                                         : "[0]";
            layer.set(
                "axis",
                fixture::json(
                    std::string{R"({"kind":"projected","direction":{"kind":"constant","value":)"} +
                    std::string{direction} + R"(},"angle":{"kind":"constant","value":)" +
                    std::string{angle} + "}}"));
        }
        if (two_axis_layers) {
            auto top = layer;
            top.set("id", fixture::json(R"("2")"));
            auto interface = fixture::json(R"({"numerator":"0","denominator":"1"})");
            layer.field("interval").elements[1] = interface;
            top.field("interval").elements[0] = std::move(interface);
            if (axes == Axes::Quarter) {
                top.field("axis").field("direction").set("value", fixture::json("[1,0,1]"));
            }
            layers.push_back(std::move(top));
        }
    }
    if (frame_shift) {
        config.field("model")
            .field("geometry")
            .elements.front()
            .field("frame")
            .set("origin", fixture::json("[0,0,0.05]"));
    }
    if (pcomp && pcomp_region != PcompRegion::Mixed) {
        auto orthotropic = fixture::fragment("Material.orthotropic3d");
        orthotropic.set("id", fixture::json(R"("2")"));
        orthotropic.set("young", fixture::json("[130,9,1]"));
        orthotropic.set("poisson", fixture::json("[0.28,0,0]"));
        auto& model = config.field("model");
        model.field("materials").elements.push_back(std::move(orthotropic));
        auto& layers = model.field("shells").elements[0].field("layers").elements;
        auto top = layers.front();
        top.set("id", fixture::json(R"("2")"));
        const std::string sign = reversed ? "1" : "-1";
        auto interface =
            fixture::json(std::string{R"({"numerator":")"} + sign + R"(","denominator":"2"})");
        auto& bottom = layers.front();
        bottom.field("interval").elements[1] = interface;
        top.field("interval").elements[0] = std::move(interface);
        if (reversed) {
            bottom.set("material", fixture::json(R"("2")"));
        } else {
            top.set("material", fixture::json(R"("2")"));
        }
        layers.push_back(std::move(top));
    }
    if (corner_claim) {
        config.field("model")
            .field("geometry")
            .elements[0]
            .field("thickness")
            .set("values", fixture::json("[0.25,0.5,0.5,0.25]"));
    }
    if (partial_group || shared_group_thickness || partial_pshell ||
        pcomp_region == PcompRegion::Mixed) {
        config.field("model")
            .field("geometry")
            .elements[0]
            .field("thickness")
            .set("values", fixture::json("[0.25,0.5,0.25,0.25]"));
    }
    auto external = fixture::fragment(is_gmsh ? "Geometry.gmsh" : "Geometry.nastran");
    if (is_gmsh) {
        external.set("reference", fixture::fragment("ReferenceGeometry.fit"));
    }
    if (conversion) {
        config.field("units").set("length", fixture::json(R"("mm")"));
        auto& reference = config.field("model").field("geometry").elements[0];
        reference.set("reference_length", fixture::json("1000"));
        reference.field("thickness").set("values", fixture::json("[250,250,250,250]"));
        reference.set("degree", fixture::json("[2,0]"));
        reference.set(
            "coefficients",
            fixture::json(
                R"([{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"10"}},{"powers":[2,0],"coefficient":{"numerator":"-1","denominator":"10"}}])"));
        if (is_gmsh) {
            auto& fit = external.field("reference");
            fit.set("degree", fixture::json("2"));
            fit.set("distance_budget", fixture::json(distance_short ? "24" : "26"));
        }
    }
    declare_path(external.field(is_gmsh ? "mesh" : "deck"), "external.dat", hex(source_text));
    declare_path(external.field("sidecar"), "sidecar.json", hex(sidecar_text));
    config.field("model").field("geometry").elements.push_back(std::move(external));
    if (shell_set) {
        auto& model = config.field("model");
        auto second_pg = model.field("geometry").elements.front();
        second_pg.set("id", fixture::json(R"("3")"));
        second_pg.field("frame").set("origin", fixture::json("[0,0,1]"));
        model.field("geometry").elements.push_back(std::move(second_pg));

        auto& shells = model.field("shells").elements;
        auto& first_layer = shells.front().field("layers").elements.front();
        first_layer.set(
            "axis",
            fixture::json(
                R"({"kind":"projected","direction":{"kind":"constant","value":[1,0,0]},"angle":{"kind":"constant","value":[0]}})"));
        auto second_shell = shells.front();
        second_shell.set("id", fixture::json(R"("2")"));
        second_shell.set("geometry", fixture::json(R"("3")"));
        second_shell.field("layers").elements.front().set("id", fixture::json(R"("2")"));
        if (late_axis_mismatch) {
            first_layer.field("axis").field("direction").set("value", fixture::json("[0,1,0]"));
        }
        shells.push_back(std::move(second_shell));
    }
    if (native_cards) {
        auto& model = config.field("model");
        model.field("beams").elements.push_back(fixture::fragment("Beam"));
        auto mass = fixture::fragment("Mass.particle");
        mass.set("mass", fixture::json("2"));
        model.field("masses").elements.push_back(std::move(mass));
        model.field("joints").elements.push_back(fixture::fragment("Joint.rigid_attachment"));
        auto& scales = config.field("numerics").field("scales").elements;
        scales.back()
            .field("constraint")
            .set("orientation", fixture::fragment("ComponentScale3.values"));
        for (const auto* scale_name : {"beam_translation", "beam_rotation", "moment_balance",
                                       "structural_angular_velocity"}) {
            scales.push_back(
                fixture::json(std::string{R"({"kind":"uniform_components","field":")"} +
                              scale_name + R"(","targets":[],"values":[1,1,1]})"));
        }
    }

    const rd::MeshSidecarBindingLimits binding{
        {limits().document,
         {256, 256, 256, 8, 8, 4096, 8192,
          uvlas::structure::ExactEstimateAdmission{fixture::budget},
          uvlas::structure::ExactTimingMode::Disabled},
         32,
         8,
         32'000},
        {4096, 16, 16, 64, 16, 16, 32'000},
        {4096, 80, 32, 32, 256, 32'000},
        {32, 32, 8, 32'000},
        {32, 32, 32'000},
        {32, 32, 32'000},
        8,
        32'000};
    auto graph = [&] {
        const Tree tree;
        write_file(tree.root / "case.json", config.wire());
        write_file(tree.root / "external.dat", source_text);
        write_file(tree.root / "sidecar.json", sidecar_text);
        core::ResourcePlan input_plan{fixture::budget};
        const auto input = fixture::take(
            rd::ReferenceInput::from_file(tree.root, "case.json", io::SourceId{12}, manifest_limits,
                                          limits(), input_plan, {}, nullptr));
        auto declarations = prepare(input);
        core::ResourcePlan graph_plan{fixture::budget};
        rd::LeafReferenceDiagnostic diagnostic;
        auto result = rd::resolve_leaf_references(
            input, std::move(declarations),
            {{{4096, 4097, 2, 4096}, 8, 256}, 2, 8192, 8194, 4, 8192, 8}, {2, 2, 1, 128, 128},
            {binding.gmsh, binding.nastran, binding.sidecar}, graph_plan, {}, nullptr, &diagnostic);
        if (lexical) {
            REQUIRE_FALSE(result);
            CHECK(result.error().code ==
                  (version ? core::FailureCode::Unsupported : core::FailureCode::InvalidInput));
            CHECK(diagnostic.reason == rd::LeafReferenceFailure::Sidecar);
            CHECK(diagnostic.pointer() == "/model/geometry/1/sidecar");
            REQUIRE(diagnostic.sidecar);
            CHECK(diagnostic.sidecar->reason == (version       ? rd::MeshSidecarFailure::Schema
                                                 : wrong_chart ? rd::MeshSidecarFailure::Scalar
                                                               : rd::MeshSidecarFailure::Key));
        }
        return result;
    }(); // The graph now owns the exact source bytes. All input files have been removed.

    if (lexical) {
        return;
    }
    REQUIRE(graph);
    auto case_copy = config;
    if (wrong_case) {
        fixture::reverse_keys(case_copy);
    }
    const auto schema = fixture::checked(case_copy);
    core::ResourcePlan plan{fixture::budget};
    std::optional<rd::MeshSidecarDiagnostic> diagnostic;
    const auto result = rd::bind_external_mesh_sidecar(schema, *graph, rd::MeshGeometryId{2},
                                                       binding, plan, {}, &diagnostic);
    if (omitted || duplicate || foreign_element || wrong_components || budget_repeat ||
        wrong_case || wrong_dimension || contradictory_group_thickness) {
        REQUIRE_FALSE(result);
        CHECK(result.error().operation == core::Operation::PrepareModel);
        CHECK(result.error().code ==
              (wrong_case ? core::FailureCode::StaleGeneration : core::FailureCode::InvalidInput));
        if (const auto& detail = diagnostic) {
            if (contradictory_group_thickness) {
                CHECK(detail->reason == rd::MeshSidecarFailure::Binding);
                CHECK(detail->physical == rd::GmshKey{2, 7});
                CHECK(detail->node);
            }
            if (omitted || duplicate || foreign_element || wrong_components || wrong_dimension) {
                CHECK(detail->reason == rd::MeshSidecarFailure::Selection);
                CHECK(detail->node);
            }
            if (omitted || duplicate || foreign_element) {
                CHECK(detail->element == (foreign_element ? 999U : omitted ? 701U : 702U));
            }
        } else {
            FAIL("Expected the sidecar binding diagnostic");
        }
        return;
    }

    REQUIRE(result);
    REQUIRE(result->faces().size() == (shell_set ? 3U : (quad ? 1U : 2U)));
    CHECK_FALSE(diagnostic);
    CHECK(result->sidecar().source().source == io::SourceId{2});
    CHECK(result->sidecar().geometry() == rd::MeshGeometryId{2});
    const auto faces = result->faces();
    if (!shell_set) {
        CHECK(faces[0].element == (quad ? 701U : 702U));
        CHECK(faces[0].pg == rd::MeshGeometryId{1});
        CHECK(faces[0].shell == rd::MeshShellId{1});
        CHECK(faces[0].sidecar_row == 0);
        if (quad) {
            const auto native_ring = axes == Axes::Reversed
                                         ? std::array<std::uint64_t, 4>{10, 99, 40, 20}
                                         : std::array<std::uint64_t, 4>{10, 20, 40, 99};
            CHECK(faces[0].nodes == native_ring);
        } else {
            CHECK(faces[0].nodes == (reversed ? std::array<std::uint64_t, 4>{10, 99, 40, 0}
                                              : std::array<std::uint64_t, 4>{10, 40, 99, 0}));
            CHECK(faces[1].element == 701);
            CHECK(faces[1].nodes == (reversed ? std::array<std::uint64_t, 4>{10, 40, 20, 0}
                                              : std::array<std::uint64_t, 4>{10, 20, 40, 0}));
            CHECK(faces[1].sidecar_node != faces[0].sidecar_node);
            CHECK(faces[1].sidecar_row == 1);
        }
    }
    if (physical) {
        using namespace uvlas::structure;
        const PgReferenceLimits pg_limits{
            partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
        const uvlas::coupling::detail::CasePgReferenceLimits preparation{
            {pg_limits.algebra, 4, {}},
            {pg_limits.algebra, 2000, 20'000'000, 100'000, {}},
            partition_test::limits(),
            pg_limits,
            100'000,
            {}};
        core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
        const rd::ExternalSidecarPgLimits qualification{
            {pg_limits, 8}, {8, 32, 32, 8, 32'000}, 8, 32'000};
        rd::ExternalSidecarPgDiagnostic failure;
        const auto to_si = [&](const io::SchemaCase& typed_case) {
            std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
            si_caps.fill(fixture::budget);
            const io::SiLimits si_limits{si_caps, pg_limits.algebra};
            const auto catalog =
                fixture::take(io::make_si_catalog(typed_case, si_limits, plan, nullptr));
            std::vector<io::SiRequirementRow> requirements;
            for (std::uint64_t i = 0; i < fixture::take(catalog.size()); ++i) {
                const auto entry = fixture::take(catalog.entry(i));
                requirements.push_back({entry.key,
                                        {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                                        io::SiDecimalInput{"1e-25"},
                                        io::SiDecimalInput{"1e-25"}});
            }
            const auto request =
                fixture::take(io::make_si_request(catalog, requirements, si_limits, plan, nullptr));
            return fixture::take(
                io::convert_case_to_si(typed_case, catalog, request, si_limits, plan, nullptr));
        };
        const auto si = to_si(schema);
        if (conversion) {
            const auto prepared = uvlas::coupling::detail::prepare_case_external_pg_bounds(
                si, *result, rd::MeshShellId{1}, preparation, qualification, 200'000, workspace,
                plan, &failure);
            if (distance_short || normal_short) {
                REQUIRE_FALSE(prepared);
                CHECK(prepared.error().code == core::FailureCode::PrecisionUnresolved);
                if (const auto& detail = failure.binding) {
                    CHECK(detail->reason == rd::MeshSidecarFailure::Chart);
                    CHECK(detail->element.has_value());
                } else {
                    FAIL("Missing conversion budget diagnostic");
                }
                return;
            }
            REQUIRE(prepared);
            const auto* admitted = &prepared->match();
            const auto expected_budget = fixture::take(
                ExactRational::from_integer_pair("13", "500", pg_limits.algebra, plan));
            CHECK(admitted->distance_budget_metres == expected_budget);
            REQUIRE(admitted->observations.size() == 1);
            CHECK(admitted->observations.values()[0].id == core::ObservationId{1});
            std::visit(
                [&](const auto& original) {
                    const auto view = original.view();
                    REQUIRE(view);
                    CHECK(view->patch == fixture::take(prepared->reference().identity()));
                    REQUIRE(view->faces.size() == 2);
                    const auto interior_distance = fixture::take(
                        ExactRational::from_integer_pair("1", "40", pg_limits.algebra, plan));
                    for (const auto& face : view->faces) {
                        CHECK(fixture::take(face.distance_upper_metres.compare(
                                  interior_distance, pg_limits.algebra, plan)) !=
                              std::strong_ordering::less);
                    }
                },
                admitted->geometry);
            return;
        }
        std::optional<rd::NastranLaminateCards> quantities;
        if (!is_gmsh) {
            auto from_snapshot = fixture::take(
                rd::prepare_referenced_nastran_deck(*graph, fixture::take(schema.identity()),
                                                    io::SourceId{1}, binding.nastran, plan, {}));
            auto coordinates = fixture::take(rd::NastranCoordinateCards::prepare(
                std::move(from_snapshot), binding.coordinates, plan, {}));
            auto cards = fixture::take(
                rd::NastranShellCards::prepare(std::move(coordinates), binding.shells, plan, {}));
            auto isotropic = fixture::take(
                rd::NastranIsotropicShellCards::prepare(std::move(cards), {32, 32'000}, plan, {}));
            quantities.emplace(fixture::take(
                rd::NastranLaminateCards::prepare(std::move(isotropic), {32, 32'000}, plan, {})));
        }
        if (native_cards) {
            using uvlas::coupling::detail::CaseNastranPgSetLimits;
            using uvlas::coupling::detail::CaseNastranReferenceDiagnostic;
            using uvlas::coupling::detail::CaseNastranReferenceInputLimits;
            using uvlas::coupling::detail::CaseNastranReferenceInputs;
            const auto* native_quantities = quantities ? &*quantities : nullptr;
            REQUIRE(native_quantities != nullptr);
            const CaseNastranReferenceInputLimits native_limits{{2, 5, 2},
                                                                {8, 8, 32'000},
                                                                {8, 8, 32'000},
                                                                {8, 32'000},
                                                                {8, 8, 32'000},
                                                                {pg_limits, 300'000, 8},
                                                                {pg_limits, 300'000},
                                                                {pg_limits, 300'000, 8},
                                                                3};
            CaseNastranReferenceDiagnostic provenance;
            const auto native = CaseNastranReferenceInputs::prepare(
                si, *result, *native_quantities, preparation, qualification, native_limits,
                workspace, plan, &provenance);
            if (thermal_link) {
                REQUIRE_FALSE(native);
                CHECK(native.error().code == core::FailureCode::Unsupported);
                if (const auto& binding_error = provenance.selection) {
                    CHECK(binding_error->target_space == io::SchemaNamespace::Joint);
                    CHECK(binding_error->case_id == 1);
                    CHECK(binding_error->element == rd::NastranElementId{907});
                } else {
                    FAIL("Missing Case binding for rejected rigid card");
                }
                const auto* card = std::get_if<rd::NastranAttachmentDiagnostic>(&provenance.card);
                REQUIRE(card != nullptr);
                CHECK(card->reason == rd::NastranAttachmentFailure::Thermal);
                CHECK(card->link == rd::NastranElementId{907});
                if (const auto field = card->position) {
                    CHECK(field->byte == source_text.find("0.001,293.15"));
                } else {
                    FAIL("Missing original ALPHA field");
                }
                const auto same_inputs = uvlas::coupling::detail::CaseNastranPlanarPgSet::prepare(
                    si, *result, *native_quantities, preparation, qualification,
                    native_limits.shells, workspace, plan);
                REQUIRE(same_inputs);
                CHECK(same_inputs->size() == 2);
                return;
            }
            REQUIRE(native);
            CHECK_FALSE(provenance.selection);
            CHECK(native->sidecar_source().source == io::SourceId{2});
            CHECK(native->sidecar_digest() == result->sidecar().digest());
            CHECK(native->deck_source().source == io::SourceId{1});
            if (const auto& shells = native->shells()) {
                REQUIRE(shells->size() == 2);
                CHECK(fixture::take(shells->at(0)).get().match().planar.shell ==
                      rd::MeshShellId{2});
            } else {
                FAIL("Missing owned Shell set");
            }
            const auto exact = [&](std::string_view number) {
                return fixture::take(ExactRational::from_decimal(number, pg_limits.algebra, plan));
            };
            REQUIRE(native->beams().size() == 1);
            if (const auto& beam = native->beams().front()) {
                CHECK(beam->binding.element == rd::NastranElementId{905});
                CHECK(beam->binding.target.value() == 1);
                CHECK(beam->physical.source.location.source == io::SourceId{1});
                CHECK(beam->physical.source.digest == native->deck_digest());
                CHECK(beam->physical.length_metres.lower == exact("1"));
                REQUIRE(beam->physical.stations.size() == 1);
                CHECK(beam->physical.stations.values().front().section.area_m2 == exact("2"));
            } else {
                FAIL("Missing owned native beam");
            }
            REQUIRE(native->masses().size() == 1);
            if (const auto& mass = native->masses().front()) {
                CHECK(mass->binding.element == rd::NastranElementId{906});
                CHECK(mass->binding.target.value() == 1);
                CHECK(mass->physical.mass_kg == exact("2"));
            } else {
                FAIL("Missing owned native mass");
            }
            REQUIRE(native->joints().size() == 1);
            if (const auto& joint = native->joints().front()) {
                CHECK(joint->binding.element == rd::NastranElementId{907});
                CHECK(joint->binding.target.value() == 1);
                if (const auto joint_view = joint->physical.view()) {
                    REQUIRE(joint_view->dependents.size() == 1);
                    CHECK(joint_view->dependents.front().point_metres[0].lower == exact("1"));
                } else {
                    FAIL("Missing owned RBE2 reference");
                }
            } else {
                FAIL("Missing joint binding");
            }
            auto short_limit = native_limits;
            --short_limit.maximum_bound_cards;
            const auto before = plan.allocated_bytes();
            const auto shortage =
                CaseNastranReferenceInputs::prepare(si, *result, *native_quantities, preparation,
                                                    qualification, short_limit, workspace, plan);
            REQUIRE_FALSE(shortage);
            CHECK(shortage.error().code == core::FailureCode::ResourceExceeded);
            CHECK(plan.allocated_bytes() == before);
            return;
        }
        if (shell_set) {
            using uvlas::coupling::detail::CaseNastranPgSetDiagnostic;
            using uvlas::coupling::detail::CaseNastranPgSetLimits;
            using uvlas::coupling::detail::CaseNastranPlanarPgSet;
            using uvlas::coupling::detail::prepare_case_nastran_planar_pg;
            const auto* native_quantities = quantities ? &*quantities : nullptr;
            REQUIRE(native_quantities != nullptr);
            const CaseNastranPgSetLimits batch_limits{2, 5, 2};
            CaseNastranPgSetDiagnostic batch_diagnostic;
            const auto prepared = CaseNastranPlanarPgSet::prepare(
                si, *result, *native_quantities, preparation, qualification, batch_limits,
                workspace, plan, &batch_diagnostic);
            CHECK(faces[0].shell == rd::MeshShellId{2});
            CHECK(faces[1].shell == rd::MeshShellId{1});
            CHECK(faces[2].shell == rd::MeshShellId{2});
            if (late_axis_mismatch) {
                REQUIRE_FALSE(prepared);
                CHECK(prepared.error().code == core::FailureCode::InvalidInput);
                CHECK(batch_diagnostic.shell == rd::MeshShellId{1});
                CHECK(batch_diagnostic.first_declaration == faces[1].sidecar_node);
                if (const auto& mismatch = batch_diagnostic.member.material_axes) {
                    CHECK(mismatch->reason == rd::NastranPgAxisFailure::Difference);
                    CHECK(mismatch->element == 701);
                    CHECK(mismatch->layer == core::LayerId{1});
                } else {
                    FAIL("Missing rejected Shell provenance");
                }
                const auto first_again = prepare_case_nastran_planar_pg(
                    si, *result, rd::MeshShellId{2}, *native_quantities, preparation, qualification,
                    2, workspace, plan);
                REQUIRE(first_again);
                CHECK(first_again->match().material_axes.matches.size() == 2);
                const auto again = CaseNastranPlanarPgSet::prepare(
                    si, *result, *native_quantities, preparation, qualification, batch_limits,
                    workspace, plan, &batch_diagnostic);
                REQUIRE_FALSE(again);
                CHECK(again.error().code == prepared.error().code);
                CHECK(batch_diagnostic.shell == rd::MeshShellId{1});
                return;
            }
            REQUIRE(prepared);
            REQUIRE(prepared->size() == 2);
            CHECK_FALSE(batch_diagnostic.shell);
            CHECK_FALSE(batch_diagnostic.member.material_axes);
            const auto& first = fixture::take(prepared->at(0)).get();
            const auto& second = fixture::take(prepared->at(1)).get();
            CHECK(first.match().planar.shell == rd::MeshShellId{2});
            CHECK(second.match().planar.shell == rd::MeshShellId{1});
            CHECK(first.match().planar.pg == rd::MeshGeometryId{3});
            CHECK(second.match().planar.pg == rd::MeshGeometryId{1});
            CHECK(first.match().material_axes.deck_digest ==
                  second.match().material_axes.deck_digest);
            CHECK(first.match().material_axes.sidecar_digest ==
                  second.match().material_axes.sidecar_digest);
            CHECK(first.match().material_axes.pg == fixture::take(first.reference().identity()));
            CHECK(second.match().material_axes.pg == fixture::take(second.reference().identity()));
            const auto first_axes = first.match().material_axes.matches.values();
            REQUIRE(first_axes.size() == 2);
            CHECK(first_axes[0].element == 900);
            CHECK(first_axes[1].element == 901);
            CHECK(first_axes[0].layer == core::LayerId{2});
            const auto second_axes = second.match().material_axes.matches.values();
            REQUIRE(second_axes.size() == 1);
            CHECK(second_axes.front().element == 701);
            CHECK(second_axes.front().layer == core::LayerId{1});
            CHECK_FALSE(prepared->at(-1));
            CHECK_FALSE(prepared->at(2));

            auto short_limit = batch_limits;
            --short_limit.maximum_shells;
            const auto shortages = CaseNastranPlanarPgSet::prepare(
                si, *result, *native_quantities, preparation, qualification, short_limit, workspace,
                plan, &batch_diagnostic);
            REQUIRE_FALSE(shortages);
            CHECK(shortages.error().code == core::FailureCode::ResourceExceeded);
            CHECK(batch_diagnostic.shell == rd::MeshShellId{1});

            short_limit = batch_limits;
            --short_limit.maximum_selection_work;
            const auto short_work =
                CaseNastranPlanarPgSet::prepare(si, *result, *native_quantities, preparation,
                                                qualification, short_limit, workspace, plan);
            REQUIRE_FALSE(short_work);
            CHECK(short_work.error().code == core::FailureCode::ResourceExceeded);

            auto foreign = case_copy;
            fixture::reverse_keys(foreign);
            const auto foreign_schema = fixture::checked(foreign);
            const auto foreign_si = to_si(foreign_schema);
            const auto before_failure = plan.allocated_bytes();
            const auto stale = CaseNastranPlanarPgSet::prepare(
                foreign_si, *result, *native_quantities, preparation, qualification, batch_limits,
                workspace, plan);
            REQUIRE_FALSE(stale);
            CHECK(stale.error().code == core::FailureCode::StaleGeneration);
            CHECK(plan.allocated_bytes() == before_failure);

            // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() is non-const.
            std::stop_source cancelled;
            REQUIRE(cancelled.request_stop());
            auto cancelled_preparation = preparation;
            cancelled_preparation.cancellation = cancelled.get_token();
            const auto stopped = CaseNastranPlanarPgSet::prepare(
                si, *result, *native_quantities, cancelled_preparation, qualification, batch_limits,
                workspace, plan);
            REQUIRE_FALSE(stopped);
            CHECK(stopped.error().code == core::FailureCode::Cancelled);
            CHECK(plan.allocated_bytes() == before_failure);
            return;
        }
        if (axes != Axes::None) {
            using uvlas::coupling::detail::prepare_case_nastran_planar_pg;
            const auto* native_quantities = quantities ? &*quantities : nullptr;
            REQUIRE(native_quantities != nullptr);
            uvlas::coupling::detail::CaseNastranPgDiagnostic diagnostic_result;
            const auto limit = two_axis_layers ? 2U : 1U;
            const auto matched = prepare_case_nastran_planar_pg(
                si, *result, rd::MeshShellId{1}, *native_quantities, preparation, qualification,
                limit, workspace, plan, &diagnostic_result);
            const auto& axis_diagnostic = diagnostic_result.material_axes;
            CHECK_FALSE(diagnostic_result.planar.binding);
            CHECK_FALSE(diagnostic_result.planar.geometry);
            const bool unsupported = axes == Axes::VariableThickness;
            const bool disagrees = axes == Axes::Difference || axes == Axes::ResidualAngle ||
                                   axes == Axes::WrongThirtySign ||
                                   axes == Axes::WrongFortyFiveQuadrant;
            if (disagrees || axes == Axes::Degenerate || unsupported) {
                REQUIRE_FALSE(matched);
                CHECK(matched.error().code == (unsupported ? core::FailureCode::Unsupported
                                                           : core::FailureCode::InvalidInput));
                if (const auto& failure_detail = axis_diagnostic) {
                    CHECK(failure_detail->element == 701);
                    CHECK(failure_detail->property == rd::NastranPropertyId{90});
                    CHECK(failure_detail->reason ==
                          (disagrees                  ? rd::NastranPgAxisFailure::Difference
                           : axes == Axes::Degenerate ? rd::NastranPgAxisFailure::CaseAxis
                                                      : rd::NastranPgAxisFailure::Thickness));
                    if (axes != Axes::VariableThickness) {
                        CHECK(failure_detail->layer ==
                              core::LayerId{axes == Axes::Difference ? 2U : 1U});
                        CHECK(failure_detail->element_field);
                        CHECK(failure_detail->ply_field);
                    }
                } else {
                    FAIL("Missing material-axis provenance");
                }
                if (axes == Axes::Difference) {
                    const auto geometry_after =
                        uvlas::coupling::detail::prepare_case_external_planar_pg(
                            si, *result, rd::MeshShellId{1}, native_quantities, preparation,
                            qualification, workspace, plan);
                    REQUIRE(geometry_after);
                    CHECK(geometry_after->match().thickness_compared);
                    REQUIRE(geometry_after->match().laminate.size() == 1);
                    if (const auto& section = geometry_after->match().laminate.values().front()) {
                        REQUIRE(section->layers.size() == 2);
                        CHECK(section->property == rd::NastranPropertyId{90});
                    } else {
                        FAIL("Missing laminate after rejected axis match");
                    }
                    const auto repeated = prepare_case_nastran_planar_pg(
                        si, *result, rd::MeshShellId{1}, *native_quantities, preparation,
                        qualification, limit, workspace, plan, &diagnostic_result);
                    REQUIRE_FALSE(repeated);
                    CHECK(repeated.error().code == matched.error().code);
                    if (const auto& detail = diagnostic_result.material_axes) {
                        CHECK(detail->layer == core::LayerId{2});
                        CHECK(detail->reason == rd::NastranPgAxisFailure::Difference);
                    } else {
                        FAIL("Axis provenance changed after retry");
                    }
                }
            } else {
                REQUIRE(matched);
                CHECK_FALSE(axis_diagnostic);
                const auto& evidence = matched->match();
                CHECK(evidence.material_axes.pg == fixture::take(matched->reference().identity()));
                CHECK(evidence.material_axes.sidecar.source == io::SourceId{2});
                CHECK(evidence.material_axes.deck.source == io::SourceId{1});
                CHECK(evidence.planar.sidecar_digest == evidence.material_axes.sidecar_digest);
                CHECK(evidence.planar.thickness_compared);
                REQUIRE(evidence.planar.laminate.size() == 1);
                if (const auto& section = evidence.planar.laminate.values().front()) {
                    CHECK(section->pg == evidence.material_axes.pg);
                    CHECK(section->layers.size() == limit);
                } else {
                    FAIL("Missing laminate in the combined result");
                }
                const auto records = evidence.material_axes.matches.values();
                REQUIRE(records.size() == (two_axis_layers ? 2U : 1U));
                const auto& row = records.front();
                CHECK(row.element == 701);
                CHECK(row.layer == core::LayerId{1});
                CHECK(row.element_field);
                CHECK(row.ply_field.has_value() == (axes != Axes::Mcid));
                CHECK(row.normal[2].sign() > 0);
                const std::string_view angle = axes == Axes::CancelledAngles ? "30"
                                               : axes == Axes::Thirty        ? "-30"
                                               : axes == Axes::Sixty         ? "-60"
                                               : axes == Axes::FortyFive     ? "45"
                                               : axes == Axes::Reversed      ? "-90"
                                               : axes == Axes::Mcid          ? "0"
                                                                             : "90";
                CHECK(row.degrees ==
                      fixture::take(ExactRational::from_decimal(angle, pg_limits.algebra, plan)));
                if (axes == Axes::Mcid) {
                    CHECK(row.tangent[0].sign() == 0);
                    CHECK(row.tangent[1].sign() > 0);
                } else if (tilted_axes) {
                    CHECK(row.tangent[0].sign() > 0);
                    CHECK(row.tangent[0] == row.tangent[1]);
                    CHECK(row.tangent[2].sign() == 0);
                }
                if (axes == Axes::Quarter) {
                    const auto& upper = records[1];
                    CHECK(upper.layer == core::LayerId{2});
                    CHECK(upper.degrees.sign() == 0);
                    CHECK(upper.ply_field != row.ply_field);
                    const auto shortage = prepare_case_nastran_planar_pg(
                        si, *result, rd::MeshShellId{1}, *native_quantities, preparation,
                        qualification, 1, workspace, plan);
                    REQUIRE_FALSE(shortage);
                    CHECK(shortage.error().code == core::FailureCode::ResourceExceeded);
                    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() is non-const.
                    std::stop_source stopped;
                    REQUIRE(stopped.request_stop());
                    auto cancelled_preparation = preparation;
                    cancelled_preparation.cancellation = stopped.get_token();
                    const auto cancelled = prepare_case_nastran_planar_pg(
                        si, *result, rd::MeshShellId{1}, *native_quantities, cancelled_preparation,
                        qualification, limit, workspace, plan);
                    REQUIRE_FALSE(cancelled);
                    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
                }
            }
            return;
        }
        const auto qualified = uvlas::coupling::detail::prepare_case_external_planar_pg(
            si, *result, rd::MeshShellId{1}, quantities ? &*quantities : nullptr, preparation,
            qualification, workspace, plan, &failure);
        if (corners == Corners::MissingBase || corners == Corners::MixedMaterial ||
            corners == Corners::WrongFraction) {
            REQUIRE_FALSE(qualified);
            CHECK(qualified.error().code == (corners == Corners::MixedMaterial
                                                 ? core::FailureCode::Unsupported
                                                 : core::FailureCode::InvalidInput));
            if (corners == Corners::WrongFraction) {
                if (const auto& detail = failure.regional_thickness) {
                    CHECK(detail->reason == rd::ExternalPgThicknessFailure::Difference);
                    CHECK(detail->node == 20);
                    CHECK(detail->field);
                    CHECK(detail->scale_field);
                } else {
                    FAIL("Missing native corner mismatch");
                }
            } else if (const auto& detail = failure.binding) {
                CHECK(detail->element == 701);
                CHECK(detail->byte);
            } else {
                FAIL("Missing property or corner diagnostic");
            }
        } else if (pcomp_region == PcompRegion::WrongTotal ||
                   pcomp_region == PcompRegion::WrongInterface) {
            REQUIRE_FALSE(qualified);
            CHECK(qualified.error().code == core::FailureCode::InvalidInput);
            if (const auto& detail = failure.binding) {
                CHECK(detail->element == 701);
                CHECK(detail->byte);
                CHECK(detail->node == faces[1].sidecar_node);
            } else {
                FAIL("Missing second property location");
            }
            if (pcomp_region == PcompRegion::WrongTotal) {
                if (const auto& regional = failure.regional_thickness) {
                    CHECK(regional->reason == rd::ExternalPgThicknessFailure::Difference);
                    CHECK(regional->element == 701);
                    CHECK(regional->property == rd::NastranPropertyId{92});
                    CHECK(regional->case_triangle);
                } else {
                    FAIL("Missing local PCOMP thickness mismatch");
                }
            } else if (const auto& laminate = failure.laminate) {
                CHECK(laminate->reason == rd::NastranPcompGeometryFailure::Interface);
                CHECK(laminate->case_layer);
                CHECK(laminate->field);
            } else {
                FAIL("Missing second PCOMP's layer mismatch");
            }
        } else if (wrong_layer_material) {
            REQUIRE_FALSE(qualified);
            CHECK(qualified.error().code == core::FailureCode::InvalidInput);
            if (const auto& binding_error = failure.binding) {
                CHECK(binding_error->target_namespace == io::SchemaNamespace::Material);
                CHECK(binding_error->target == 2);
                CHECK(binding_error->node);
            } else {
                FAIL("Missing mapped material disagreement");
            }
        } else if (material_shift || thickness_shift || shared_group_thickness ||
                   wrong_pshell_thickness || frame_shift) {
            REQUIRE_FALSE(qualified);
            CHECK(qualified.error().code == ((chart_shift || frame_shift)
                                                 ? core::FailureCode::GeometryInvalid
                                                 : core::FailureCode::InvalidInput));
            if (chart_shift || frame_shift) {
                CHECK(failure.geometry.has_value());
                CHECK_FALSE(failure.quantity.has_value());
            } else if (const auto& mismatch = failure.quantity) {
                CHECK(mismatch->quantity == (material_shift
                                                 ? rd::ImportedQuantityField::Young1
                                                 : rd::ImportedQuantityField::FullThickness));
                if (material_shift) {
                    CHECK(mismatch->card_position.has_value());
                } else if (const auto& local = failure.regional_thickness) {
                    CHECK(local->reason == rd::ExternalPgThicknessFailure::Difference);
                    CHECK(local->element ==
                          ((shared_group_thickness || wrong_pshell_thickness) ? 701U : 702U));
                    CHECK(local->field.has_value() == wrong_pshell_thickness);
                    CHECK(local->sidecar_node);
                    CHECK(local->case_triangle);
                    CHECK(local->overlap_vertex);
                } else {
                    FAIL("Missing regional thickness provenance");
                }
            } else {
                FAIL("Missing exact quantity mismatch");
            }
        } else {
            REQUIRE(qualified);
            const auto& matched_case = qualified->match();
            CHECK(matched_case.thickness_compared);
            if (const auto& regional = matched_case.thickness_regions) {
                const auto claims = regional->matches.values();
                REQUIRE(claims.size() == ((partial_group || partial_pshell || quad ||
                                           pcomp_region == PcompRegion::Mixed)
                                              ? 1U
                                              : 2U));
                CHECK(claims.front().element == (quad ? 701U : 702U));
                const auto base = fixture::take(
                    ExactRational::from_integer_pair("1", "4", pg_limits.algebra, plan));
                if (corner_claim) {
                    const auto* samples =
                        std::get_if<rd::ExternalPgNodalThickness>(&claims.front().value_metres);
                    REQUIRE(samples);
                    CHECK(samples->count == (quad ? 4U : 3U));
                    CHECK(samples->values[0].node == 10);
                    CHECK(samples->values[0].value_metres == base);
                    CHECK(samples->values[0].field);
                    CHECK_FALSE(samples->values[0].scale_field);
                    CHECK(samples->values[1].node == (quad ? 20U : 40U));
                    CHECK(samples->values[1].value_metres ==
                          fixture::take(
                              ExactRational::from_integer_pair("1", "2", pg_limits.algebra, plan)));
                    CHECK(samples->values[1].field);
                    CHECK(samples->values[1].scale_field.has_value() == quad);
                } else {
                    CHECK(std::get<ExactRational>(claims.front().value_metres) == base);
                }
                if (claims.size() == 2 && is_gmsh) {
                    CHECK(claims[1].sidecar_node == claims[0].sidecar_node);
                }
                if (claims.size() == 2 && !is_gmsh && !corner_claim) {
                    CHECK((claims[1].field == claims[0].field) == !separate_pcomp);
                }
                if (pcomp) {
                    CHECK(claims.front().property == rd::NastranPropertyId{90});
                    if (separate_pcomp) {
                        REQUIRE(claims.size() == 2);
                        CHECK(claims[1].element == 701);
                        CHECK(claims[1].property == rd::NastranPropertyId{92});
                    }
                }
            } else {
                FAIL("Missing whole-face regional thickness evidence");
            }
            const std::uint64_t expected_materials = is_gmsh                     ? 0
                                                     : separate_pcomp            ? 4
                                                     : corner_claim              ? 3
                                                     : (pcomp || partial_pshell) ? 2
                                                                                 : 1;
            CHECK(matched_case.materials_compared == expected_materials);
            if (pcomp) {
                const auto sections = matched_case.laminate.values();
                REQUIRE(sections.size() == (separate_pcomp ? 2U : 1U));
                for (std::size_t i = 0; i < sections.size(); ++i) {
                    if (const auto& matched = sections[i]) {
                        CHECK(matched->winding == (reversed ? rd::MeshSidecarOrientation::Reversed
                                                            : rd::MeshSidecarOrientation::Same));
                        CHECK(matched->property == rd::NastranPropertyId{i == 0 ? 90U : 92U});
                        const auto plies = matched->layers.values();
                        if (pcomp_region == PcompRegion::Mixed) {
                            REQUIRE(plies.size() == 1);
                            CHECK(plies.front().material.id == rd::NastranMaterialId{6});
                        } else {
                            REQUIRE(plies.size() == 2);
                            CHECK(plies.front().material.id ==
                                  rd::NastranMaterialId{reversed ? 9U : 6U});
                            CHECK(plies.front().ply_thickness_metres ==
                                  fixture::take(ExactRational::from_integer_pair(
                                      reversed ? "3" : "1", "16", pg_limits.algebra, plan)));
                            CHECK(plies.back().material.id ==
                                  rd::NastranMaterialId{reversed ? 6U : 9U});
                        }
                        CHECK(matched->pg == fixture::take(qualified->reference().identity()));
                    } else {
                        FAIL("Missing exact PCOMP layer comparison");
                    }
                }
            }
            CHECK(matched_case.pg == rd::MeshGeometryId{1});
            if (const auto admitted = matched_case.geometry.view()) {
                CHECK(admitted->patch == fixture::take(qualified->reference().identity()));
                CHECK(admitted->faces.size() == faces.size());
                CHECK(admitted->faces.front().element == (quad ? 701U : 702U));
            } else {
                FAIL("Missing owned geometric correspondence");
            }
        }
    }
    auto boundary = binding;
    boundary.faces = static_cast<std::uint64_t>(faces.size());
    boundary.work = result->work();
    core::ResourcePlan exact{fixture::budget};
    CHECK(rd::bind_external_mesh_sidecar(schema, *graph, rd::MeshGeometryId{2}, boundary, exact, {},
                                         nullptr));
    auto short_faces = boundary;
    --short_faces.faces;
    core::ResourcePlan short_storage{fixture::budget};
    const auto shortage = rd::bind_external_mesh_sidecar(schema, *graph, rd::MeshGeometryId{2},
                                                         short_faces, short_storage, {}, nullptr);
    REQUIRE_FALSE(shortage);
    CHECK(shortage.error().code == core::FailureCode::ResourceExceeded);
    --boundary.work;
    core::ResourcePlan short_work{fixture::budget};
    const auto short_result = rd::bind_external_mesh_sidecar(schema, *graph, rd::MeshGeometryId{2},
                                                             boundary, short_work, {}, nullptr);
    REQUIRE_FALSE(short_result);
    CHECK(short_result.error().code == core::FailureCode::ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() is non-const.
    std::stop_source cancelled;
    REQUIRE(cancelled.request_stop());
    const auto stopped = rd::bind_external_mesh_sidecar(
        schema, *graph, rd::MeshGeometryId{2}, binding, plan, cancelled.get_token(), nullptr);
    REQUIRE_FALSE(stopped);
    CHECK(stopped.error().code == core::FailureCode::Cancelled);
}
