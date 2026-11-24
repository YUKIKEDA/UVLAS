#include "gmsh_mesh.hpp"
#include "reference_declarations_fixture.hpp"
#include "reference_input.hpp"
#include "reference_root_snapshot.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
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

std::string write_root_inventory(const fs::path& directory) {
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
    write(geometry[1].field("deck"), "model.bdf");
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
