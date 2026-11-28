#include "case_schema_test_adapter.hpp"
#include "checkpoint_directory.hpp"
#include "checkpoint_directory_name.hpp"
#include "checkpoint_test_fixture.hpp"
#include "reference_file_test_fixture.hpp"

#ifdef _WIN32
#include "windows_checkpoint_query.hpp"
#include "windows_reference_test_fixture.hpp"
#endif

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace {

namespace cp = uvlas::io::detail;
namespace core = uvlas::core;
namespace fs = std::filesystem;
namespace fixture = uvlas::test::checkpoint_fixture;
namespace reference = uvlas::test::reference_file_fixture;

class TestDirectory {
public:
    TestDirectory() {
        static std::atomic_uint_fast64_t serial{0};
        const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
        path_ =
            fs::current_path() / ("checkpoint-directory-test-" + std::to_string(now) + "-" +
                                  std::to_string(serial.fetch_add(1, std::memory_order_relaxed)));
        if (!fs::create_directory(path_))
            throw std::logic_error("Checkpoint test directory is already present");
    }
    TestDirectory(const TestDirectory&) = delete;
    TestDirectory& operator=(const TestDirectory&) = delete;
    ~TestDirectory() {
        std::error_code ignored;
        (void)fs::remove_all(path_, ignored);
    }
    [[nodiscard]] const fs::path& path() const noexcept {
        return path_;
    }

private:
    fs::path path_;
};

template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::logic_error("Expected a checkpoint result");
    return std::move(*value);
}

const cp::CheckpointManifestLimits limits{
    {8192,
     8,
     256,
     128,
     64,
     32,
     4096,
     {uvlas::io::ParserAllocationMode::EstimateAdmission, 4'000'000}},
    {4, 1024},
    8,
    2048};
constexpr cp::CheckpointDirectoryLimits files{{8192, 8193, 16, 8192}, 16};
constexpr cp::CheckpointPublicationLimits publication{8, 8192, 16'384};
constexpr cp::CheckpointDiscoveryLimits discovery{16, 8};
constexpr uvlas::io::SourceId input{554};
enum class Newest : std::uint8_t { Corrupt, Incompatible, ExtraChunk };

#ifdef _WIN32
TEST_CASE("Windows checkpoint directory query classifies native resource exhaustion",
          "[IO4-06][negative]") {
    CHECK_FALSE(cp::windows_checkpoint_query_failure(ERROR_NO_MORE_FILES));
    for (const auto error : {ERROR_NOT_ENOUGH_MEMORY, ERROR_OUTOFMEMORY, ERROR_NO_SYSTEM_RESOURCES,
                             ERROR_NOT_ENOUGH_QUOTA, ERROR_COMMITMENT_LIMIT}) {
        CHECK(cp::windows_checkpoint_query_failure(error) == core::FailureCode::ResourceExceeded);
    }
    CHECK(cp::windows_checkpoint_query_failure(ERROR_INVALID_DATA) ==
          core::FailureCode::CheckpointCorrupt);
}
#endif

TEST_CASE("Latest checkpoint lookup stops if a candidate exhausts its enumeration limit",
          "[IO4-06][negative]") {
    const TestDirectory root;
    core::ResourcePlan plan{1'000'000};
    const cp::CheckpointFingerprints context{};
    constexpr std::array<std::uint64_t, 1> values{42};
    const auto chunk = take(cp::CheckpointChunk::encode({}, std::span<const std::uint64_t>{values},
                                                        limits.chunk, plan, {}));
    const std::array named{cp::NamedCheckpointChunk{"node_id", chunk}};
    REQUIRE(cp::save_checkpoint_directory(root.path(), 1, named, context, limits, publication, plan,
                                          {}));
    REQUIRE(cp::save_checkpoint_directory(root.path(), 2, named, context, limits, publication, plan,
                                          {}));
    const auto latest = root.path() / "checkpoint-2";
    REQUIRE(fs::copy_file(latest / "manifest.json", latest / "notes.txt"));

    auto entry_cap = files;
    entry_cap.maximum_entries = 2;
    REQUIRE(cp::read_checkpoint_directory(root.path() / "checkpoint-1", input, context, limits,
                                          entry_cap, plan, {}));
    const auto result = cp::read_latest_checkpoint(root.path(), input, context, limits, entry_cap,
                                                   discovery, plan, {});
    REQUIRE_FALSE(result);
    CHECK(result.error().code == core::FailureCode::ResourceExceeded);
    CHECK(result.error().operation == core::Operation::Recover);
}

TEST_CASE("Checkpoint directories round-trip typed data and discover the latest valid number",
          "[IO4-06][N4-02]") {
    const auto newest = GENERATE(Newest::Corrupt, Newest::Incompatible, Newest::ExtraChunk);
    CAPTURE(newest);
    const TestDirectory root;
    core::ResourcePlan plan{6'000'000};
    auto identity = cp::CheckpointFingerprints{};
    identity.version.bytes[0] = std::byte{0xa7};
    const auto chunk = [&](std::uint64_t value) {
        const std::array numbers{value};
        return take(cp::CheckpointChunk::encode({}, std::span<const std::uint64_t>{numbers},
                                                limits.chunk, plan, {}));
    };
    const auto a = chunk(42);
    const auto b = chunk(47);
    const std::array first{cp::NamedCheckpointChunk{"node_id", a}};
    const std::array second{cp::NamedCheckpointChunk{"node_id", b}};
    const auto save = [&](std::uint64_t number, std::span<const cp::NamedCheckpointChunk> chunks,
                          const cp::CheckpointFingerprints& context) {
        std::optional<cp::CheckpointPublicationDiagnostic> diagnostic;
        const auto result = cp::save_checkpoint_directory(
            root.path(), number, chunks, context, limits, publication, plan, {}, &diagnostic);
        if (diagnostic) {
            UNSCOPED_INFO("native save stage " << static_cast<unsigned>(diagnostic->stage)
                                               << ", error "
                                               << diagnostic->native_error.value_or(0));
        }
        return result;
    };
    REQUIRE(save(4, first, identity));
    REQUIRE(save(12, second, identity));
    CHECK_FALSE(save(12, first, identity));

    const auto good = root.path() / "checkpoint-12";
    if (newest == Newest::Incompatible) {
        auto other = identity;
        other.model.bytes[0] = std::byte{1};
        REQUIRE(save(14, first, other));
        const auto incompatible = cp::read_checkpoint_directory(
            root.path() / "checkpoint-14", input, identity, limits, files, plan, {});
        REQUIRE_FALSE(incompatible);
        CHECK(incompatible.error().code == core::FailureCode::CheckpointIncompatible);
    } else if (newest == Newest::Corrupt) {
        const auto broken = root.path() / "checkpoint-13";
        REQUIRE(fs::create_directory(broken));
        REQUIRE(fs::copy_file(good / "manifest.json", broken / "manifest.json"));
        std::ofstream damage(broken / "chunk-node_id.bin", std::ios::binary | std::ios::trunc);
        REQUIRE(damage.is_open());
        const std::array<char, 8> wrong{};
        damage.write(wrong.data(), static_cast<std::streamsize>(wrong.size()));
        damage.close();
        REQUIRE(damage.good());
        const auto corrupt =
            cp::read_checkpoint_directory(broken, input, identity, limits, files, plan, {});
        REQUIRE_FALSE(corrupt);
        CHECK(corrupt.error().code == core::FailureCode::CheckpointCorrupt);
    } else {
        REQUIRE(save(14, first, identity));
        const auto extra = root.path() / "checkpoint-14";
        REQUIRE(fs::copy_file(good / "chunk-node_id.bin", extra / "chunk-extra.bin"));
        const auto corrupted =
            cp::read_checkpoint_directory(extra, input, identity, limits, files, plan, {});
        REQUIRE_FALSE(corrupted);
        CHECK(corrupted.error().code == core::FailureCode::CheckpointCorrupt);
    }
    // Even a complete set is invisible under a staging name.
    fs::copy(good, root.path() / ".uvlas-checkpoint-fixture");
    const auto staging = cp::read_checkpoint_directory(root.path() / ".uvlas-checkpoint-fixture",
                                                       input, identity, limits, files, plan, {});
    REQUIRE_FALSE(staging);
    CHECK(staging.error().code == core::FailureCode::CheckpointCorrupt);

    const auto latest = take(cp::read_latest_checkpoint(root.path(), input, identity, limits, files,
                                                        discovery, plan, {}));
    CHECK(latest.number == 12);
    CHECK(take(latest.manifest.at(0)).name == "node_id");
    const auto& accepted = take(latest.payloads.at(0)).get();
    auto data = take(accepted.decode(plan, {}));
    CHECK(std::get<core::PlannedBuffer<std::uint64_t>>(data).values().front() == 47);
    const auto first_checkpoint = take(cp::read_checkpoint_directory(
        root.path() / "checkpoint-4", input, identity, limits, files, plan, {}));
    CHECK(first_checkpoint.number == 4);
    auto previous = take(take(first_checkpoint.payloads.at(0)).get().decode(plan, {}));
    CHECK(std::get<core::PlannedBuffer<std::uint64_t>>(previous).values().front() == 42);

    auto entry_cap = files;
    entry_cap.maximum_entries = 2;
    REQUIRE(cp::read_checkpoint_directory(good, input, identity, limits, entry_cap, plan, {}));
    REQUIRE(fs::copy_file(good / "manifest.json", good / "notes.txt"));
    const auto excess_entries =
        cp::read_checkpoint_directory(good, input, identity, limits, entry_cap, plan, {});
    REQUIRE_FALSE(excess_entries);
    CHECK(excess_entries.error().code == core::FailureCode::ResourceExceeded);
    REQUIRE(cp::read_checkpoint_directory(good, input, identity, limits, files, plan, {}));

    const auto denied =
        cp::read_latest_checkpoint(root.path(), input, identity, limits, files, {16, 1}, plan, {});
    REQUIRE_FALSE(denied);
    CHECK(denied.error().code == core::FailureCode::ResourceExceeded);

    if (newest == Newest::ExtraChunk) {
        for (int i = 0; i < 40; ++i)
            REQUIRE(fs::copy_file(good / "manifest.json",
                                  good / ("note-" + std::to_string(i) + ".txt")));
        entry_cap.maximum_entries = 43;
        REQUIRE(cp::read_checkpoint_directory(good, input, identity, limits, entry_cap, plan, {}));
        --entry_cap.maximum_entries;
        const auto bounded =
            cp::read_checkpoint_directory(good, input, identity, limits, entry_cap, plan, {});
        REQUIRE_FALSE(bounded);
        CHECK(bounded.error().code == core::FailureCode::ResourceExceeded);
    }
}

constexpr cp::CaseCheckpointDirectoryLimits case_files{{256, 128, 8}, files};

[[nodiscard]] std::string case_wire(std::string_view path, std::string_view digest) {
    using namespace uvlas::test::schema_adapter;
    auto wire_case = root("coupled_flutter");
    if (!digest.empty()) {
        auto restart = fragment("Restart");
        restart.set("directory", string(path));
        restart.set("manifest_sha256", string(digest));
        wire_case.set("restart", std::move(restart));
    }
    return wire_case.wire();
}

void write_independent_checkpoint(const fs::path& directory) {
    const auto write_bytes = [&](const char* name, const auto& bytes) {
        const auto data = reinterpret_cast<const char*>(bytes.data());
        reference::write_file(directory / name, {data, bytes.size()});
    };
    REQUIRE(fs::create_directories(directory));
    reference::write_file(directory / "manifest.json", fixture::manifest);
    write_bytes("chunk-velocity.bin", fixture::one);
    write_bytes("chunk-node_id.bin", fixture::forty_two);
}

[[nodiscard]] cp::ReferenceInput acquired_case(const fs::path& parent, std::string_view filename,
                                               std::string_view wire, bool from_file) {
    core::ResourcePlan plan{128ULL * 1024 * 1024};
    if (from_file) {
        reference::write_file(parent / filename, wire);
        return take(cp::ReferenceInput::from_file(parent, filename, input,
                                                  reference::manifest_limits, reference::limits(),
                                                  plan, {}, nullptr));
    }
    return take(cp::ReferenceInput::from_memory(reference::bytes_of(wire), input,
                                                take(cp::ReadRoot::open(parent, {})),
                                                reference::limits(), plan, {}, nullptr));
}

class RestoreWorkingDirectory {
public:
    RestoreWorkingDirectory() : previous_(fs::current_path()) {}
    RestoreWorkingDirectory(const RestoreWorkingDirectory&) = delete;
    RestoreWorkingDirectory& operator=(const RestoreWorkingDirectory&) = delete;
    ~RestoreWorkingDirectory() {
        std::error_code ignored;
        fs::current_path(previous_, ignored);
    }

private:
    fs::path previous_;
};

TEST_CASE("Case restart selects the owned relative directory from file or memory",
          "[FD-08][IO4-06][IO4-07][integration][negative][verification]") {
    const bool from_file = GENERATE(false, true);
    CAPTURE(from_file);
    const TestDirectory dir;
    const RestoreWorkingDirectory working_directory;
    const auto initial_root = dir.path() / "case-root";
    const auto decoy = dir.path() / "different-cwd";
    const auto owner_case =
        case_wire(R"(runs\discard/.././checkpoint-12)", fixture::manifest_sha256);
    write_independent_checkpoint(initial_root / "runs" / "checkpoint-12");
    write_independent_checkpoint(initial_root / "runs" / "checkpoint-13");
    REQUIRE(fs::create_directories(decoy / "case-root" / "runs" / "checkpoint-12"));
    REQUIRE(fs::create_directories(decoy / "runs" / "checkpoint-12"));
    reference::write_file(decoy / "case-root" / "runs" / "checkpoint-12" / "manifest.json", "{}");
    reference::write_file(decoy / "runs" / "checkpoint-12" / "manifest.json", "{}");
    core::ResourcePlan plan{12'000'000};

    fs::current_path(dir.path());
    const auto declared = acquired_case("case-root", "original.json", owner_case, from_file);
    fs::current_path(decoy);
    // Linux retains the original inode across rename. Windows roots instead deny delete sharing.
    auto active = initial_root;
#ifdef __linux__
    const auto relocated = dir.path() / "relocated-root";
    fs::rename(initial_root, relocated);
    REQUIRE(fs::create_directories(initial_root / "runs" / "checkpoint-12"));
    reference::write_file(initial_root / "runs" / "checkpoint-12" / "manifest.json", "{}");
    active = relocated;
#endif
    const auto checkpoint = active / "runs" / "checkpoint-12";
    cp::CaseCheckpointDirectoryDiagnostic diagnostic;
    const auto accepted = take(cp::read_case_checkpoint_directory(
        declared, input, fixture::fingerprints, limits, case_files, plan, {}, &diagnostic));
    CHECK(accepted.number == 12);
    CHECK(diagnostic.pointer().empty());
    const auto original_view = accepted.manifest_file.view();
    REQUIRE(original_view);
    CHECK(std::ranges::equal(original_view->bytes, reference::bytes_of(fixture::manifest)));

    // Whitespace preserves the parsed manifest but must fail the source digest before JSON and
    // enumeration limits get a chance to reject the envelope.
    reference::write_file(checkpoint / "manifest.json", std::string{fixture::manifest} + ' ');
    REQUIRE(cp::read_checkpoint_directory(checkpoint, input, fixture::fingerprints, limits, files,
                                          plan, {}));
    auto no_enumeration = case_files;
    no_enumeration.container.maximum_entries = 0;
    const auto rejected = cp::read_case_checkpoint_directory(
        declared, input, fixture::fingerprints, limits, no_enumeration, plan, {}, &diagnostic);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::CheckpointCorrupt);
    CHECK(rejected.error().operation == core::Operation::Recover);
    CHECK(diagnostic.pointer() == "/restart/manifest_sha256");

    const auto updated = acquired_case(
        active, "updated.json",
        case_wire("runs/checkpoint-12", fixture::manifest_sha256_with_space), from_file);
    auto no_parse = limits;
    no_parse.json.maximum_nodes = 1;
    const auto pre_parse = cp::read_case_checkpoint_directory(
        declared, input, fixture::fingerprints, no_parse, case_files, plan, {}, &diagnostic);
    REQUIRE_FALSE(pre_parse);
    CHECK(pre_parse.error().code == core::FailureCode::CheckpointCorrupt);
    CHECK(diagnostic.pointer() == "/restart/manifest_sha256");
    const auto parser_limit = cp::read_case_checkpoint_directory(
        updated, input, fixture::fingerprints, no_parse, case_files, plan, {}, &diagnostic);
    REQUIRE_FALSE(parser_limit);
    CHECK(parser_limit.error().code == core::FailureCode::ResourceExceeded);
    CHECK(diagnostic.stage == cp::CaseCheckpointDirectoryStage::Envelope);

    const auto bounded = cp::read_case_checkpoint_directory(
        updated, input, fixture::fingerprints, limits, no_enumeration, plan, {}, &diagnostic);
    REQUIRE_FALSE(bounded);
    CHECK(bounded.error().code == core::FailureCode::ResourceExceeded);
    const auto matching = take(cp::read_case_checkpoint_directory(
        updated, input, fixture::fingerprints, limits, case_files, plan, {}));
    CHECK(matching.number == 12);

    auto context = fixture::fingerprints;
    context.material.bytes[0] = std::byte{1};
    const auto incompatible =
        cp::read_case_checkpoint_directory(updated, input, context, limits, case_files, plan, {});
    REQUIRE_FALSE(incompatible);
    CHECK(incompatible.error().code == core::FailureCode::CheckpointIncompatible);

    REQUIRE(fs::remove(checkpoint / "chunk-node_id.bin"));
    const auto missing = cp::read_case_checkpoint_directory(updated, input, fixture::fingerprints,
                                                            limits, case_files, plan, {});
    REQUIRE_FALSE(missing);
    CHECK(missing.error().code == core::FailureCode::CheckpointCorrupt);
    // No retry chooses the intact newer checkpoint. The last verified result owns its bytes.
    const auto previous = take(take(accepted.payloads.at(0)).get().decode(plan, {}));
    CHECK(std::get<core::PlannedBuffer<std::uint64_t>>(previous).values().front() == 42);
    const auto retained_view = accepted.manifest_file.view();
    REQUIRE(retained_view);
    CHECK(std::ranges::equal(retained_view->bytes, reference::bytes_of(fixture::manifest)));
}

TEST_CASE("Case restart rejects invalid relative paths and enforces bounded acquisition",
          "[FD-08][IO4-06][integration][negative][verification]") {
    using Code = core::FailureCode;
    using Stage = cp::CaseCheckpointDirectoryStage;

    const TestDirectory dir;
    const auto parent = dir.path() / "root";
    REQUIRE(fs::create_directory(parent));
    write_independent_checkpoint(parent / "runs" / "checkpoint-12");
    reference::write_file(parent / "regular", "");
    reference::write_file(parent / "checkpoint-16", "");

    core::ResourcePlan plan{128ULL * 1024 * 1024};
    cp::CaseCheckpointDirectoryDiagnostic diagnostic;
    const auto acquire = [&](std::string_view path) {
        return acquired_case(parent, "", case_wire(path, fixture::manifest_sha256), false);
    };
    const auto expect = [&](const cp::ReferenceInput& declared,
                            cp::CaseCheckpointDirectoryLimits caps, Code code, Stage stage) {
        const auto result = cp::read_case_checkpoint_directory(
            declared, input, fixture::fingerprints, limits, caps, plan, {}, &diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == code);
        CHECK(result.error().operation == core::Operation::Recover);
        CHECK(diagnostic.stage == stage);
    };

    for (std::string_view path :
         {"/checkpoint-12", R"(X:checkpoint-12)", R"(X:\checkpoint-12)",
          R"(\\network\share\checkpoint-12)", "run:stream/checkpoint-12", "../checkpoint-12",
          "runs/../../checkpoint-12", "runs/checkpoint-12/", "runs/..", "runs/checkpoint-012",
          "runs/checkpoint-18446744073709551616", "runs/latest"}) {
        CAPTURE(path);
        expect(acquire(path), case_files, Code::InvalidInput, Stage::Path);
        CHECK(diagnostic.pointer() == "/restart/directory");
        CHECK(diagnostic.path);
    }
    for (std::string_view path : {"missing/checkpoint-12", "regular/checkpoint-12", "checkpoint-16",
                                  "runs/checkpoint-99"}) {
        CAPTURE(path);
        expect(acquire(path), case_files, Code::CheckpointCorrupt, Stage::Open);
        CHECK(diagnostic.pointer() == "/restart/directory");
        CHECK(diagnostic.open);
    }
#ifdef _WIN32
    using uvlas::test::windows_reference_fixture::Junction;
    const Junction link{parent / "link", parent / "runs"};
    const Junction leaf{parent / "checkpoint-12", parent / "runs" / "checkpoint-12"};
#elif defined(__linux__)
    fs::create_directory_symlink(parent / "runs", parent / "link");
    fs::create_directory_symlink(parent / "runs" / "checkpoint-12", parent / "checkpoint-12");
#endif
    for (std::string_view path : {"link/checkpoint-12", "checkpoint-12"}) {
        CAPTURE(path);
        expect(acquire(path), case_files, Code::CheckpointCorrupt, Stage::Open);
        REQUIRE(diagnostic.open);
#ifdef _WIN32
        CHECK(diagnostic.open->reason == cp::WindowsReferenceOpenReason::Reparse);
        CHECK(diagnostic.open->normalized_offset == 0);
#elif defined(__linux__)
        CHECK(*diagnostic.open == cp::LinuxReferenceOpenReason::Link);
#endif
    }

    const auto declared = acquire("runs//./checkpoint-12");
    auto caps = case_files;
    caps.path.input_bytes = 20;
    expect(declared, caps, Code::ResourceExceeded, Stage::Path);
    REQUIRE(diagnostic.path);
    CHECK(diagnostic.path->reason == cp::ReferencePathReason::InputBytes);
    caps = case_files;
    caps.path.normalized_bytes = 17;
    expect(declared, caps, Code::ResourceExceeded, Stage::Path);
    REQUIRE(diagnostic.path);
    CHECK(diagnostic.path->reason == cp::ReferencePathReason::NormalizedBytes);
    caps = case_files;
    caps.path.normalized_components = 1;
    expect(declared, caps, Code::ResourceExceeded, Stage::Path);
    REQUIRE(diagnostic.path);
    CHECK(diagnostic.path->reason == cp::ReferencePathReason::NormalizedDepth);
    REQUIRE(cp::read_case_checkpoint_directory(declared, input, fixture::fingerprints, limits,
                                               {{21, 18, 2}, files}, plan, {}));

    core::ResourcePlan small_plan{18}; // Just the normalized bytes, no native scratch.
    const auto small = cp::read_case_checkpoint_directory(
        declared, input, fixture::fingerprints, limits, case_files, small_plan, {}, &diagnostic);
    REQUIRE_FALSE(small);
    CHECK(small.error().code == Code::ResourceExceeded);
    CHECK(small.error().operation == core::Operation::Recover);
    CHECK(diagnostic.stage == Stage::Open);

    const auto absent = acquired_case(parent, "", case_wire("", ""), false);
    expect(absent, case_files, Code::InvalidInput, Stage::Restart);
    const auto no_root_wire = case_wire("runs/checkpoint-12", fixture::manifest_sha256);
    const auto no_root =
        take(cp::ReferenceInput::from_memory(reference::bytes_of(no_root_wire), input, std::nullopt,
                                             reference::limits(), plan, {}, nullptr));
    expect(no_root, case_files, Code::InvalidInput, Stage::Path);
    CHECK(diagnostic.pointer() == "/restart/directory");

    std::stop_source cancellation;
    REQUIRE(cancellation.request_stop());
    core::ResourcePlan no_io{0};
    const auto stopped = cp::read_case_checkpoint_directory(declared, input, fixture::fingerprints,
                                                            limits, case_files, no_io,
                                                            cancellation.get_token(), &diagnostic);
    REQUIRE_FALSE(stopped);
    CHECK(stopped.error().code == Code::Cancelled);
    CHECK(stopped.error().operation == core::Operation::Recover);
    CHECK(diagnostic.pointer().empty());
    CHECK_FALSE(diagnostic.path);
    CHECK_FALSE(diagnostic.open);
}

#ifdef __linux__
TEST_CASE("Linux checkpoint save and recovery reject an ancestor symlink", "[IO4-06][negative]") {
    const TestDirectory root;
    const auto actual = root.path() / "actual";
    const auto parent = actual / "runs";
    const auto alias = root.path() / "alias";
    REQUIRE(fs::create_directories(parent));
    fs::create_directory_symlink(actual, alias);
    const auto via_alias = alias / "runs";
    core::ResourcePlan plan{1'000'000};
    const cp::CheckpointFingerprints context{};
    const std::array<std::uint64_t, 1> values{42};
    const auto chunk = take(cp::CheckpointChunk::encode({}, std::span<const std::uint64_t>{values},
                                                        limits.chunk, plan, {}));
    const std::array named{cp::NamedCheckpointChunk{"node_id", chunk}};
    CHECK_FALSE(
        cp::save_checkpoint_directory(via_alias, 1, named, context, limits, publication, plan, {}));
    CHECK(fs::is_empty(parent));
    REQUIRE(
        cp::save_checkpoint_directory(parent, 1, named, context, limits, publication, plan, {}));
    CHECK_FALSE(cp::read_checkpoint_directory(via_alias / "checkpoint-1", input, context, limits,
                                              files, plan, {}));
    REQUIRE(cp::read_checkpoint_directory(parent / "checkpoint-1", input, context, limits, files,
                                          plan, {}));
}
#endif

} // namespace
