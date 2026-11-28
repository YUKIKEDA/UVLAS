#pragma once

#include "checkpoint_manifest.hpp"

#include <array>

namespace uvlas::test::checkpoint_fixture {
namespace cp = io::detail;

// Independent little-endian payloads and digest fixture. See the native Case
// checkpoint verification record for the hashlib replay.
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
constexpr std::string_view manifest = R"({
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
constexpr std::string_view manifest_sha256 =
    "9aae026dedadccda54af409c1c16abc143fcfa1e67293e11bdd297d56a24c208";
constexpr std::string_view manifest_sha256_with_space =
    "9d885403109f994cd4d7dbee6eb9c95c267a0ec1989cafd96d3e2ce230190a02";
} // namespace uvlas::test::checkpoint_fixture
