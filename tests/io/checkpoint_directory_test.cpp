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
constexpr cp::ReferenceFileLimits files{8192, 8193, 16, 8192};
constexpr cp::CheckpointPublicationLimits publication{8, 8192, 16'384};
constexpr cp::CheckpointDiscoveryLimits discovery{16, 8};
constexpr uvlas::io::SourceId input{554};

TEST_CASE("Checkpoint directories round-trip typed data and discover the latest valid number",
          "[IO4-06][N4-02]") {
    const bool incompatible_newest = GENERATE(false, true);
    CAPTURE(incompatible_newest);
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
    if (incompatible_newest) {
        auto other = identity;
        other.model.bytes[0] = std::byte{1};
        REQUIRE(save(14, first, other));
        const auto incompatible = cp::read_checkpoint_directory(
            root.path() / "checkpoint-14", input, identity, limits, files, plan, {});
        REQUIRE_FALSE(incompatible);
        CHECK(incompatible.error().code == core::FailureCode::CheckpointIncompatible);
    } else {
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

    const auto denied =
        cp::read_latest_checkpoint(root.path(), input, identity, limits, files, {16, 1}, plan, {});
    REQUIRE_FALSE(denied);
    CHECK(denied.error().code == core::FailureCode::ResourceExceeded);
}

} // namespace
