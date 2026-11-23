#include "checkpoint_directory.hpp"
#include "checkpoint_directory_name.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <array>
#include <atomic>
#include <chrono>
#include <fstream>
#include <stdexcept>

namespace {

namespace cp = uvlas::io::detail;
namespace core = uvlas::core;
namespace fs = std::filesystem;

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
enum class Newest { Corrupt, Incompatible, ExtraChunk };

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
