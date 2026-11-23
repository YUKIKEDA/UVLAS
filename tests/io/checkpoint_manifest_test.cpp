#include "checkpoint_manifest.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <uvlas/core/failure.hpp>

namespace {

namespace cp = uvlas::io::detail;
namespace core = uvlas::core;
using uvlas::io::SourceId;

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
constexpr auto fingerprints = []() noexcept {
    cp::CheckpointFingerprints result{};
    result.model.bytes[0] = std::byte{0xb2};
    result.version.bytes[31] = std::byte{0x09};
    return result;
}();

constexpr std::array<std::byte, 8> forty_two{std::byte{42}};
constexpr std::array<std::byte, 8> one{std::byte{0},    std::byte{0},   std::byte{0},
                                       std::byte{0},    std::byte{0},   std::byte{0},
                                       std::byte{0xf0}, std::byte{0x3f}};
constexpr std::string_view fixture = R"({
  "schema": "uvlas.checkpoint/1",
  "byte_order": "little",
  "fingerprints": {
    "model": "b200000000000000000000000000000000000000000000000000000000000000",
    "basis": "0000000000000000000000000000000000000000000000000000000000000000",
    "material": "0000000000000000000000000000000000000000000000000000000000000000",
    "constraint": "0000000000000000000000000000000000000000000000000000000000000000",
    "algorithm": "0000000000000000000000000000000000000000000000000000000000000000",
    "precision": "0000000000000000000000000000000000000000000000000000000000000000",
    "version": "0000000000000000000000000000000000000000000000000000000000000009"
  },
  "chunks": [
    {"name":"velocity","dtype":"f64","shape":[1],"bytes":8,
     "sha256":"6c3c396ed6b5c36dcae172271f462051b1266b851e92df3deea8ac65478fd712"},
    {"name":"node_id","dtype":"u64","shape":[],"bytes":8,
     "sha256":"ed049108bc18f2c64369e8d0ea42850bdd1a7d1dd340cfde716315579702a76c"}
  ]
})";

template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::logic_error("Expected a checkpoint result");
    return std::move(*value);
}

[[nodiscard]] std::span<const std::byte> wire(std::string_view text) noexcept {
    return std::as_bytes(std::span{text.data(), text.size()});
}

[[nodiscard]] std::string replaced(std::string_view old, std::string_view replacement) {
    auto text = std::string{fixture};
    const auto at = text.find(old);
    REQUIRE(at != std::string::npos);
    text.replace(at, old.size(), replacement);
    return text;
}

TEST_CASE("Checkpoint manifests validate an independent payload set and retain its words",
          "[IO4-06][N4-02]") {
    core::ResourcePlan plan{1'000'000};
    auto input = std::string{fixture};
    auto manifest = take(
        cp::CheckpointManifest::read(wire(input), SourceId{1}, fingerprints, limits, plan, {}));
    input.clear();
    REQUIRE(manifest.size() == 2);
    CHECK(take(manifest.at(0)).name == "node_id");
    const auto velocity_entry = take(manifest.at(1));
    CHECK(velocity_entry.name == "velocity");
    CHECK(velocity_entry.header.dtype == cp::CheckpointDtype::F64);
    REQUIRE(velocity_entry.header.shape.size() == 1);
    CHECK(velocity_entry.header.shape.front() == 1);

    std::array velocity = one;
    std::array id = forty_two;
    const std::array payloads{cp::CheckpointPayload{"velocity", velocity},
                              cp::CheckpointPayload{"node_id", id}};
    auto chunks =
        take(cp::authenticate_checkpoint_payloads(manifest, payloads, limits.chunk, plan, {}));
    velocity.fill(std::byte{});
    id.fill(std::byte{});
    REQUIRE(chunks.size() == 2);
    auto recovered_id = take(take(chunks.at(0)).get().decode(plan, {}));
    auto recovered_vel = take(take(chunks.at(1)).get().decode(plan, {}));
    CHECK(std::get<core::PlannedBuffer<std::uint64_t>>(recovered_id).values().front() == 42);
    CHECK(std::get<core::PlannedBuffer<double>>(recovered_vel).values().front() == 1.0);

    const auto& original_owner = manifest;
    auto moved = std::move(manifest);
    CHECK(original_owner.size() == 0);
    CHECK_FALSE(original_owner.at(0));
    CHECK(take(moved.at(0)).header.shape.empty());
}

TEST_CASE("A checkpoint manifest rejects incompatible products and corrupt declarations",
          "[IO4-06][N4-02]") {
    auto expected = fingerprints;
    auto input = std::string{fixture};
    auto caps = limits;
    auto code = core::FailureCode::CheckpointCorrupt;
    SECTION("schema version") {
        input = replaced("uvlas.checkpoint/1", "uvlas.checkpoint/2");
        code = core::FailureCode::CheckpointIncompatible;
    }
    SECTION("byte order") {
        input = replaced("\"little\"", "\"big\"");
        code = core::FailureCode::CheckpointIncompatible;
    }
    SECTION("current context fingerprint changed") {
        expected.algorithm.bytes[16] = std::byte{1};
        code = core::FailureCode::CheckpointIncompatible;
    }
    SECTION("undefined word type") {
        input = replaced("\"u64\"", "\"i32\"");
        code = core::FailureCode::CheckpointIncompatible;
    }
    SECTION("path traversal") {
        input = replaced("node_id", "../node_id");
    }
    SECTION("embedded nul") {
        input = replaced("node_id", "node_id\\u0000");
    }
    SECTION("noncanonical case") {
        input = replaced("node_id", "Node_id");
    }
    SECTION("duplicate chunk") {
        input = replaced("node_id", "velocity");
    }
    SECTION("duplicate JSON key") {
        input = replaced("\"bytes\":8", "\"bytes\":8,\"bytes\":8");
    }
    SECTION("unknown field") {
        input = replaced("\"bytes\":8", "\"surplus\":8");
    }
    SECTION("float is not a length") {
        input = replaced("\"bytes\":8", "\"bytes\":8.0");
    }
    SECTION("shape differs from length") {
        input = replaced("\"shape\":[1]", "\"shape\":[2]");
    }
    SECTION("invalid shape extent") {
        input = replaced("\"shape\":[1]", "\"shape\":[-1]");
    }
    SECTION("malformed digest") {
        input = replaced("6c3c396e", "xc3c396e");
    }
    SECTION("truncated document") {
        input.pop_back();
    }
    SECTION("payload quota") {
        caps.maximum_total_payload_bytes = 15;
        code = core::FailureCode::ResourceExceeded;
    }
    SECTION("chunk quota") {
        caps.maximum_chunks = 1;
        code = core::FailureCode::ResourceExceeded;
    }
    SECTION("rank quota") {
        caps.chunk.maximum_rank = 0;
        code = core::FailureCode::ResourceExceeded;
    }
    core::ResourcePlan plan{1'000'000};
    const auto result =
        cp::CheckpointManifest::read(wire(input), SourceId{2}, expected, caps, plan, {});
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    CHECK(result.error().operation == core::Operation::Recover);
}

TEST_CASE("Checkpoint recovery validates every named payload before publishing a set",
          "[IO4-06][N4-02]") {
    core::ResourcePlan plan{1'000'000};
    const auto manifest = take(
        cp::CheckpointManifest::read(wire(fixture), SourceId{3}, fingerprints, limits, plan, {}));
    std::array payloads{cp::CheckpointPayload{"velocity", one},
                        cp::CheckpointPayload{"node_id", forty_two},
                        cp::CheckpointPayload{"extra", one}};
    auto input = std::span<const cp::CheckpointPayload>{payloads}.first(2);
    auto code = core::FailureCode::CheckpointCorrupt;
    std::stop_source stop;
    SECTION("missing entry") {
        input = input.first(1);
    }
    SECTION("extra entry") {
        input = payloads;
    }
    SECTION("duplicate entry") {
        payloads[1].name = "velocity";
    }
    SECTION("unexpected entry") {
        payloads[0].name = "pressure";
    }
    SECTION("bad last chunk digest") {
        payloads[0].bytes = forty_two;
    }
    SECTION("wrong file size") {
        payloads[0].bytes = std::span{one}.first(4);
    }
    SECTION("interrupted read") {
        REQUIRE(stop.request_stop());
        code = core::FailureCode::Cancelled;
    }
    const auto restored =
        cp::authenticate_checkpoint_payloads(manifest, input, limits.chunk, plan, stop.get_token());
    REQUIRE_FALSE(restored);
    CHECK(restored.error().code == code);
}

TEST_CASE("Checkpoint manifest writing is canonical and honors output limits", "[IO4-06][N4-02]") {
    core::ResourcePlan plan{1'000'000};
    constexpr std::array<std::uint64_t, 1> id{42};
    constexpr std::array speed{1.0};
    constexpr std::array<core::Index, 1> shape{1};
    const auto a = take(cp::CheckpointChunk::encode({}, std::span<const std::uint64_t>{id},
                                                    limits.chunk, plan, {}));
    const auto b = take(
        cp::CheckpointChunk::encode(shape, std::span<const double>{speed}, limits.chunk, plan, {}));
    const std::array unsorted{cp::NamedCheckpointChunk{"velocity", b},
                              cp::NamedCheckpointChunk{"node_id", a}};
    const auto ordered = std::array{unsorted[1], unsorted[0]};
    auto encoded = take(cp::write_checkpoint_manifest(unsorted, fingerprints, limits, plan, {}));
    const auto another =
        take(cp::write_checkpoint_manifest(ordered, fingerprints, limits, plan, {}));
    CHECK(std::ranges::equal(encoded.values(), another.values()));
    const auto manifest = take(cp::CheckpointManifest::read(
        std::as_bytes(encoded.values()), SourceId{4}, fingerprints, limits, plan, {}));
    CHECK(take(manifest.at(0)).name == "node_id");
    CHECK(take(manifest.at(1)).header.byte_count == 8);

    auto inputs = unsorted;
    auto caps = limits;
    auto code = core::FailureCode::InvalidInput;
    SECTION("name grammar") {
        inputs[1].name = "bad\\name";
    }
    SECTION("duplicate names") {
        inputs[1].name = "velocity";
    }
    SECTION("output size limit") {
        caps.json.maximum_source_bytes = encoded.size() - 1;
        code = core::FailureCode::ResourceExceeded;
    }
    const auto result = cp::write_checkpoint_manifest(inputs, fingerprints, caps, plan, {});
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    CHECK(result.error().operation == core::Operation::SaveCheckpoint);
}

} // namespace
