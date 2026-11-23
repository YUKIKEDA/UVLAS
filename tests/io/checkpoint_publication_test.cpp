#include "checkpoint_directory_name.hpp"
#include "checkpoint_publication.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <barrier>
#include <chrono>
#include <fstream>
#include <limits>
#include <thread>
#ifdef __linux__
#include <csignal>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace {

namespace fs = std::filesystem;
using namespace uvlas::io::detail;
using uvlas::core::FailureCode;
using uvlas::core::ResourcePlan;

constexpr CheckpointPublicationLimits limits{3, 1024, 4096};
constexpr std::string_view manifest = R"({"schema":"uvlas.checkpoint/1","byte_order":"little"})";

struct Directory {
    fs::path path;
    Directory() {
        const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
        for (int attempt = 0; attempt < 32; ++attempt) {
            const auto candidate =
                fs::temp_directory_path() / ("uvlas-io4-06-publication-" + std::to_string(tick) +
                                             "-" + std::to_string(attempt));
            if (fs::create_directory(candidate)) {
                path = candidate;
                return;
            }
        }
        FAIL("Could not reserve a private output directory");
    }
    Directory(const Directory&) = delete;
    Directory& operator=(const Directory&) = delete;
    ~Directory() {
        if (!path.empty()) {
            std::error_code ignored;
            fs::remove_all(path, ignored);
        }
    }
};

auto bytes_of(std::string_view text) noexcept {
    return std::as_bytes(std::span{text});
}

consteval uvlas::io::SourceDigest digest_of(std::string_view text) {
    constexpr std::string_view digits = "0123456789abcdef";
    uvlas::io::SourceDigest result{};
    for (std::size_t i = 0; i < result.bytes.size(); ++i) {
        result.bytes[i] =
            static_cast<std::byte>(16 * digits.find(text[2 * i]) + digits.find(text[2 * i + 1]));
    }
    return result;
}

constexpr std::string_view payload{"\0\1\2\3\4\5\6\7", 8};
constexpr std::array<uvlas::core::Index, 1> words{1};
constexpr std::array<uvlas::core::Index, 1> zero{0};
constexpr auto payload_digest =
    digest_of("8a851ff82ee7048ad09ec3847f1ddf44944104d2cbd17ef4e3db22c6785a0d45");
constexpr auto empty_digest =
    digest_of("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");

const std::array chunks{
    CheckpointPublicationChunk{
        "first", {{CheckpointDtype::U64, words, 8, payload_digest}, bytes_of(payload)}},
    CheckpointPublicationChunk{"empty", {{CheckpointDtype::I64, zero, 0, empty_digest}, {}}}};

[[nodiscard]] std::string read(const fs::path& file) {
    std::ifstream stream(file, std::ios::binary);
    REQUIRE(stream.is_open());
    return {std::istreambuf_iterator<char>{stream}, {}};
}

[[nodiscard]] std::size_t temp_count(const fs::path& path) {
    std::size_t count = 0;
    for (const auto& entry : fs::directory_iterator{path}) {
        if (entry.path().filename().string().starts_with(".uvlas-checkpoint-")) {
            ++count;
        }
    }
    return count;
}

TEST_CASE("Native publication makes exactly the verified files visible", "[IO4-06]") {
    const Directory dir;
    const auto stale = dir.path / ".uvlas-checkpoint-interrupted";
    REQUIRE(fs::create_directory(stale));
    std::ofstream(stale / "chunk-old.bin", std::ios::binary) << "incomplete";
    ResourcePlan plan{2 * sizeof(std::string_view) + manifest.size()};
    std::optional<CheckpointPublicationDiagnostic> diagnostic;
    const auto saved = publish_checkpoint_directory(dir.path, 42, chunks, bytes_of(manifest),
                                                    limits, plan, {}, &diagnostic);
    if (diagnostic) {
        INFO("native stage " << static_cast<int>(diagnostic->stage) << ", error "
                             << diagnostic->native_error.value_or(0));
        REQUIRE(saved);
    }
    REQUIRE(saved);
    CHECK_FALSE(diagnostic);
    CHECK(read(dir.path / "checkpoint-42" / "manifest.json") == manifest);
    CHECK(read(dir.path / "checkpoint-42" / "chunk-first.bin") == payload);
    CHECK(read(dir.path / "checkpoint-42" / "chunk-empty.bin").empty());
    CHECK(std::distance(fs::directory_iterator{dir.path / "checkpoint-42"},
                        fs::directory_iterator{}) == 3);
    CHECK(read(stale / "chunk-old.bin") == "incomplete");
    CHECK_FALSE(fs::exists(stale / "manifest.json"));
    CHECK(temp_count(dir.path) == 1);
}

TEST_CASE("Failed checkpoint saves preserve prior paths and remove only their staging files",
          "[IO4-06]") {
    const Directory dir;
    std::optional<CheckpointPublicationDiagnostic> diagnostic;
    ResourcePlan plan{4096};
    REQUIRE(publish_checkpoint_directory(dir.path, 41, chunks, bytes_of(manifest), limits, plan, {},
                                         &diagnostic));
    REQUIRE(fs::create_directory(dir.path / "checkpoint-42"));
    const auto rejected = publish_checkpoint_directory(dir.path, 42, chunks, bytes_of(manifest),
                                                       limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(rejected);
    REQUIRE(diagnostic);
    if (diagnostic) {
        CHECK(diagnostic->stage == CheckpointPublicationStage::DestinationExists);
    }
    CHECK(fs::is_empty(dir.path / "checkpoint-42"));

    auto mismatch = chunks;
    mismatch[1].chunk.header.digest = payload_digest;
    const auto corrupt = publish_checkpoint_directory(dir.path, 43, mismatch, bytes_of(manifest),
                                                      limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(corrupt);
    CHECK(corrupt.error().code == FailureCode::CheckpointCorrupt);
    REQUIRE(diagnostic);
    if (diagnostic) {
        CHECK(diagnostic->stage == CheckpointPublicationStage::Verify);
        CHECK(diagnostic->chunk == 1);
    }
    CHECK_FALSE(fs::exists(dir.path / "checkpoint-43"));

    // MSVC request_stop is non-const.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source cancel;
    static_cast<void>(cancel.request_stop());
    const auto cancelled = publish_checkpoint_directory(
        dir.path, 44, chunks, bytes_of(manifest), limits, plan, cancel.get_token(), &diagnostic);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == FailureCode::Cancelled);
    CHECK_FALSE(fs::exists(dir.path / "checkpoint-44"));

    auto duplicate = chunks;
    duplicate[1].logical_name = duplicate[0].logical_name;
    const auto duplicated = publish_checkpoint_directory(
        dir.path, 45, duplicate, bytes_of(manifest), limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(duplicated);
    CHECK(duplicated.error().code == FailureCode::InvalidInput);
    auto invalid_name = chunks;
    invalid_name[0].logical_name = "../escape";
    const auto traversal = publish_checkpoint_directory(
        dir.path, 48, invalid_name, bytes_of(manifest), limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(traversal);
    CHECK(traversal.error().code == FailureCode::InvalidInput);
    ResourcePlan insufficient{2 * sizeof(std::string_view) + manifest.size() - 1};
    const auto limited = publish_checkpoint_directory(dir.path, 46, chunks, bytes_of(manifest),
                                                      limits, insufficient, {}, &diagnostic);
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == FailureCode::ResourceExceeded);
    CHECK_FALSE(fs::exists(dir.path / "checkpoint-46"));
    CHECK(read(dir.path / "checkpoint-41" / "manifest.json") == manifest);
    CHECK(read(dir.path / "checkpoint-41" / "chunk-first.bin") == payload);
    CHECK(temp_count(dir.path) == 0);
}

TEST_CASE("Simultaneous checkpoint writers cannot replace the winning directory", "[IO4-06]") {
    const Directory dir;
    constexpr std::string_view first = R"({"schema":"uvlas.checkpoint/1","writer":"first"})";
    constexpr std::string_view second = R"({"schema":"uvlas.checkpoint/1","writer":"second"})";
    std::array<bool, 2> success{};
    std::array<std::optional<CheckpointPublicationDiagnostic>, 2> diagnostics{};
    std::barrier begin{2};
    const auto write = [&](std::size_t ordinal, std::string_view text) {
        ResourcePlan plan{4096};
        begin.arrive_and_wait();
        success[ordinal] = publish_checkpoint_directory(dir.path, 99, chunks, bytes_of(text),
                                                        limits, plan, {}, &diagnostics[ordinal])
                               .has_value();
    };
    {
        const std::jthread a{write, 0, first};
        const std::jthread b{write, 1, second};
    }
    REQUIRE(success[0] != success[1]);
    const auto winner = success[0] ? 0U : 1U;
    const auto loser = 1U - winner;
    CHECK_FALSE(diagnostics[winner]);
    const auto& lost = diagnostics[loser];
    REQUIRE(lost);
    if (lost) {
        CHECK(lost->stage == CheckpointPublicationStage::DestinationExists);
    }
    CHECK(read(dir.path / "checkpoint-99" / "manifest.json") == (success[0] ? first : second));
    CHECK(read(dir.path / "checkpoint-99" / "chunk-first.bin") == payload);
    CHECK(temp_count(dir.path) == 0);
}

TEST_CASE("Only canonical checkpoint sequence names can select a published directory", "[IO4-06]") {
    constexpr auto largest = std::numeric_limits<std::uint64_t>::max();
    CHECK(std::string_view{make_checkpoint_directory_name(0).data()} == "checkpoint-0");
    CHECK(std::string_view{make_checkpoint_directory_name(largest).data()} ==
          "checkpoint-18446744073709551615");
    CHECK(parse_checkpoint_directory_name("checkpoint-0") == 0);
    CHECK(parse_checkpoint_directory_name("checkpoint-18446744073709551615") == largest);
    for (const auto name : {"checkpoint-01", "checkpoint--1", "checkpoint-+1", "checkpoint-1/",
                            "checkpoint-18446744073709551616", ".uvlas-checkpoint-42-0"}) {
        CHECK_FALSE(parse_checkpoint_directory_name(name));
    }
}

#ifdef __linux__
TEST_CASE("Linux rejects checkpoint output through an ancestor symlink", "[IO4-06]") {
    const Directory dir;
    const auto actual = dir.path / "actual";
    const auto runs = actual / "runs";
    const auto alias = dir.path / "alias";
    REQUIRE(fs::create_directories(runs));
    fs::create_directory_symlink(actual, alias);
    REQUIRE(fs::is_directory(alias / "runs"));

    ResourcePlan plan{4096};
    std::optional<CheckpointPublicationDiagnostic> diagnostic;
    const auto denied = publish_checkpoint_directory(alias / "runs", 7, chunks, bytes_of(manifest),
                                                     limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(denied);
    CHECK(denied.error().code == FailureCode::InvalidInput);
    REQUIRE(diagnostic);
    if (diagnostic) {
        CHECK(diagnostic->stage == CheckpointPublicationStage::OpenParent);
    }
    CHECK(fs::is_empty(runs));

    REQUIRE(publish_checkpoint_directory(runs, 7, chunks, bytes_of(manifest), limits, plan, {}));
    CHECK(read(runs / "checkpoint-7" / "manifest.json") == manifest);
    CHECK(read(runs / "checkpoint-7" / "chunk-first.bin") == payload);
}

TEST_CASE("A killed writer cannot make a partial directory a checkpoint", "[IO4-06]") {
    const Directory dir;
    ResourcePlan plan{4096};
    REQUIRE(
        publish_checkpoint_directory(dir.path, 12, chunks, bytes_of(manifest), limits, plan, {}));
    constexpr std::string_view second_payload{"\0\1\2\3\4\5\6\7\0\1\2\3\4\5\6\7", 16};
    constexpr std::array<uvlas::core::Index, 1> two_words{2};
    constexpr auto second_digest =
        digest_of("d2a1c23b0d2514a61539937da65e56acf4fbd42f22631050e07284a8703dccb2");
    auto interrupted_chunks = chunks;
    interrupted_chunks[1] = {
        "second", {{CheckpointDtype::U64, two_words, 16, second_digest}, bytes_of(second_payload)}};
    const auto child = fork();
    REQUIRE(child >= 0);
    if (child == 0) {
        const rlimit no_core{0, 0};
        const rlimit file_limit{12, 12};
        if (setrlimit(RLIMIT_CORE, &no_core) != 0 || setrlimit(RLIMIT_FSIZE, &file_limit) != 0 ||
            std::signal(SIGXFSZ, SIG_DFL) == SIG_ERR) {
            _exit(91);
        }
        ResourcePlan child_plan{4096};
        const auto result = publish_checkpoint_directory(
            dir.path, 13, interrupted_chunks, bytes_of(manifest), limits, child_plan, {});
        _exit(result ? 0 : 92);
    }
    int status = 0;
    REQUIRE(waitpid(child, &status, 0) == child);
    REQUIRE(WIFSIGNALED(status));
    CHECK(WTERMSIG(status) == SIGXFSZ);
    CHECK_FALSE(fs::exists(dir.path / "checkpoint-13"));
    CHECK(read(dir.path / "checkpoint-12" / "manifest.json") == manifest);
    CHECK(temp_count(dir.path) == 1);
    for (const auto& entry : fs::directory_iterator{dir.path}) {
        const auto base = entry.path().filename().string();
        if (base.starts_with(".uvlas-checkpoint-")) {
            CHECK_FALSE(parse_checkpoint_directory_name(base));
            CHECK(read(entry.path() / "chunk-first.bin") == payload);
            CHECK(read(entry.path() / "chunk-second.bin") == second_payload.substr(0, 12));
            CHECK_FALSE(fs::exists(entry.path() / "manifest.json"));
        }
    }
    REQUIRE(publish_checkpoint_directory(dir.path, 13, interrupted_chunks, bytes_of(manifest),
                                         limits, plan, {}));
    CHECK(read(dir.path / "checkpoint-13" / "chunk-first.bin") == payload);
    CHECK(read(dir.path / "checkpoint-13" / "chunk-second.bin") == second_payload);
}
#endif
} // namespace
