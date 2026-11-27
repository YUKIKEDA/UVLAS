#include "checkpoint_directory.hpp"

#include <array>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

namespace cp = uvlas::io::detail;
namespace core = uvlas::core;
namespace fs = std::filesystem;

namespace {
const cp::CheckpointManifestLimits manifest{
    {8192, 8, 256, 128, 64, 32, 4096, {uvlas::io::ParserAllocationMode::EstimateAdmission, 4'000'000}},
    {4, 1024}, 8, 2048};
constexpr cp::CheckpointDirectoryLimits file_limits{{8192, 8193, 16, 8192}, 16};
constexpr cp::CheckpointPublicationLimits publication{8, 8192, 16'384};
constexpr cp::CheckpointDiscoveryLimits discovery{16, 8};
constexpr uvlas::io::SourceId source{554};

void require(bool condition, const char* failure) {
    if (!condition) {
        throw std::runtime_error(failure);
    }
}

template <class T> T take(core::Result<T> result, const char* failure) {
    require(result.has_value(), failure);
    return std::move(*result);
}

std::uint64_t first_value(const cp::LoadedCheckpointDirectory& result, core::ResourcePlan& plan) {
    const auto payload = take(result.payloads.at(0), "checkpoint first payload");
    auto data = take(payload.get().decode(plan, {}), "checkpoint decode");
    return std::get<core::PlannedBuffer<std::uint64_t>>(data).values().front();
}

void run(const fs::path& parent) {
    require(fs::create_directory(parent), "create fresh ext4 test directory");
    core::ResourcePlan plan{10'000'000};
    auto identity = cp::CheckpointFingerprints{};
    identity.version.bytes[0] = std::byte{0xa7};

    const std::array<std::uint64_t, 1> first_data{42};
    const std::array<std::uint64_t, 2> second_data{100, 200};
    const std::array<core::Index, 1> shape{2};
    const auto first = take(cp::CheckpointChunk::encode({}, std::span{first_data}, manifest.chunk, plan, {}),
                            "first chunk encoding");
    const auto second = take(cp::CheckpointChunk::encode(shape, std::span{second_data}, manifest.chunk, plan, {}),
                             "second chunk encoding");
    const std::array chunks{cp::NamedCheckpointChunk{"first", first},
                            cp::NamedCheckpointChunk{"second", second}};
    const auto save = [&](std::uint64_t number, std::span<const cp::NamedCheckpointChunk> values,
                          core::ResourcePlan& resources) {
        return cp::save_checkpoint_directory(parent, number, values, identity, manifest, publication,
                                             resources, {});
    };
    require(save(12, std::span{chunks}.first(1), plan).has_value(), "save checkpoint 12");
    const auto child = fork();
    require(child >= 0, "fork");
    if (child == 0) {
        const rlimit no_core{0, 0};
        const rlimit file_limit{12, 12};
        if (setrlimit(RLIMIT_CORE, &no_core) != 0 || setrlimit(RLIMIT_FSIZE, &file_limit) != 0 ||
            std::signal(SIGXFSZ, SIG_DFL) == SIG_ERR) {
            _exit(91);
        }
        core::ResourcePlan child_plan{10'000'000};
        _exit(save(13, chunks, child_plan) ? 0 : 92);
    }
    int status = 0;
    require(waitpid(child, &status, 0) == child && WIFSIGNALED(status) && WTERMSIG(status) == SIGXFSZ,
            "writer did not stop during the second chunk");
    require(!fs::exists(parent / "checkpoint-13"), "interrupted publication visible");
    fs::path staging;
    for (const auto& entry : fs::directory_iterator(parent)) {
        if (entry.path().filename().string().starts_with(".uvlas-checkpoint-")) {
            require(staging.empty(), "more than one staging directory");
            staging = entry.path();
        }
    }
    require(!staging.empty() && fs::file_size(staging / "chunk-first.bin") == 8 &&
                fs::file_size(staging / "chunk-second.bin") == 12 &&
                !fs::exists(staging / "manifest.json"),
            "wrong interrupted file sizes or manifest");
    const auto rejected = cp::read_checkpoint_directory(staging, source, identity, manifest, file_limits, plan, {});
    require(!rejected && rejected.error().code == core::FailureCode::CheckpointCorrupt, "staging accepted");
    auto latest = take(cp::read_latest_checkpoint(parent, source, identity, manifest, file_limits, discovery, plan, {}),
                       "read latest after interruption");
    require(latest.number == 12 && first_value(latest, plan) == 42, "latest after interruption changed");
    std::cout << "interrupted: SIGXFSZ, chunks=8/12 bytes, no manifest, explicit=CheckpointCorrupt, latest=12, value=42\n";

    require(save(13, chunks, plan).has_value(), "publish retry");
    require(!save(13, std::span{chunks}.first(1), plan), "retry replaced an existing checkpoint");
    const auto after_retry = take(cp::read_latest_checkpoint(parent, source, identity, manifest, file_limits, discovery, plan, {}),
                  "read latest after retry");
    require(after_retry.number == 13 && first_value(after_retry, plan) == 42, "latest after retry wrong");
    auto payload = take(after_retry.payloads.at(1), "checkpoint second payload");
    auto data = take(payload.get().decode(plan, {}), "decode second payload");
    const auto values = std::get<core::PlannedBuffer<std::uint64_t>>(data).values();
    require(values.size() == 2 && values[0] == 100 && values[1] == 200, "published second values wrong");
    std::cout << "published retry: latest=13, values=42/100/200, duplicate number rejected\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "expected one fresh directory argument");
        run(argv[1]);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
