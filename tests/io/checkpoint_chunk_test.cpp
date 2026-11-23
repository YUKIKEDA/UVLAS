#include "checkpoint_chunk.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

using namespace uvlas::io::detail;
using uvlas::core::FailureCode;
using uvlas::core::Index;
using uvlas::core::ResourcePlan;

constexpr CheckpointChunkLimits limits{4, std::uint64_t{128} * 1024};
constexpr std::array<Index, 1> vector_shape{3};

[[nodiscard]] CheckpointChunkView chunk_view(const CheckpointChunk& chunk) {
    const auto view = chunk.view();
    if (!view) {
        throw std::logic_error("Expected an owned chunk");
    }
    return *view;
}

[[nodiscard]] std::vector<std::byte> from_hex(std::string_view text) {
    const auto digit = [](char value) {
        return static_cast<unsigned>(value <= '9' ? value - '0' : value - 'a' + 10);
    };
    std::vector<std::byte> bytes;
    for (std::size_t i = 0; i < text.size(); i += 2) {
        bytes.push_back(static_cast<std::byte>((digit(text[i]) << 4) | digit(text[i + 1])));
    }
    return bytes;
}

template <class T>
void check_integer_payload(std::span<const T> values, CheckpointDtype dtype,
                           std::string_view golden_hex) {
    ResourcePlan plan{88};
    auto chunk = CheckpointChunk::encode(vector_shape, values, limits, plan, {});
    REQUIRE(chunk);
    const auto view = chunk_view(*chunk);
    CHECK(view.header.dtype == dtype);
    CHECK(std::ranges::equal(view.bytes, from_hex(golden_hex)));
    const auto wire = from_hex(golden_hex);
    const auto received = CheckpointChunk::authenticate({view.header, wire}, limits, plan, {});
    REQUIRE(received);
    auto recovered = received->decode(plan, {});
    REQUIRE(recovered);
    CHECK(std::ranges::equal(std::get<uvlas::core::PlannedBuffer<T>>(*recovered).values(), values));
}

TEST_CASE("Typed chunks preserve the little-endian representation and SHA-256", "[IO4-06][N4-02]") {
    constexpr std::array<std::uint64_t, 3> words{0x8000000000000000ULL, 0x3ff0000000000000ULL,
                                                 0x3ff0000000000001ULL};
    std::array<double, 3> values{std::bit_cast<double>(words[0]), std::bit_cast<double>(words[1]),
                                 std::bit_cast<double>(words[2])};
    const auto expected = from_hex("0000000000000080000000000000f03f010000000000f03f");
    const auto digest =
        from_hex("25034d7dda9aa7a57fca39b65bf810be1185a6e645002f08341c7edd69066253");
    ResourcePlan plan{88};
    auto chunk =
        CheckpointChunk::encode(vector_shape, std::span<const double>{values}, limits, plan, {});
    REQUIRE(chunk);
    values.fill(42.0);
    const auto view = chunk_view(*chunk);
    CHECK(view.header.dtype == CheckpointDtype::F64);
    CHECK(view.header.byte_count == 24);
    CHECK(std::ranges::equal(view.header.shape, vector_shape));
    CHECK(std::ranges::equal(view.bytes, expected));
    CHECK(std::ranges::equal(view.header.digest.bytes, digest));
    auto verified = CheckpointChunk::authenticate(view, limits, plan, {});
    REQUIRE(verified);
    auto moved = std::move(*verified);
    CHECK_FALSE(verified->view());
    const auto unused = verified->decode(plan, {});
    REQUIRE_FALSE(unused);
    CHECK(unused.error().code == FailureCode::InvalidInput);
    auto decoded = moved.decode(plan, {});
    REQUIRE(decoded);
    const auto restored = std::get<uvlas::core::PlannedBuffer<double>>(*decoded).values();
    for (std::size_t i = 0; i < words.size(); ++i) {
        CHECK(std::bit_cast<std::uint64_t>(restored[i]) == words[i]);
    }
    CHECK(plan.allocated_bytes() == 88);

    constexpr std::array<std::int64_t, 3> signed_values{std::numeric_limits<std::int64_t>::min(),
                                                        -1, 0};
    constexpr std::array<std::uint64_t, 3> unsigned_values{0, 0x20000000000001ULL,
                                                           0xffffffffffffffffULL};
    check_integer_payload(std::span<const std::int64_t>{signed_values}, CheckpointDtype::I64,
                          "0000000000000080ffffffffffffffff0000000000000000");
    check_integer_payload(std::span<const std::uint64_t>{unsigned_values}, CheckpointDtype::U64,
                          "00000000000000000100000000002000ffffffffffffffff");
}

TEST_CASE("Recovery rejects corrupt declarations and bytes before allocating a chunk",
          "[IO4-06][N4-02]") {
    constexpr std::array<Index, 1> shape{1};
    constexpr std::array<std::uint64_t, 1> values{42};
    ResourcePlan plan{16};
    auto chunk =
        CheckpointChunk::encode(shape, std::span<const std::uint64_t>{values}, limits, plan, {});
    REQUIRE(chunk);
    const auto valid = chunk_view(*chunk);
    auto candidate = valid;
    auto expected_error = FailureCode::CheckpointCorrupt;
    std::stop_source stop;
    std::array<std::byte, 9> changed{};
    std::ranges::copy(valid.bytes, changed.begin());
    constexpr std::array<Index, 2> invalid_shape{0, -1};
    SECTION("declared byte count differs from shape") {
        candidate.header.byte_count = 16;
    }
    SECTION("truncated payload") {
        candidate.bytes = valid.bytes.first(7);
    }
    SECTION("trailing bytes") {
        candidate.bytes = changed;
    }
    SECTION("changed byte") {
        changed[0] ^= std::byte{1};
        candidate.bytes = std::span{changed}.first(8);
    }
    SECTION("changed digest") {
        candidate.header.digest.bytes[31] ^= std::byte{1};
    }
    SECTION("negative extent is corrupt even after a zero extent") {
        candidate.header.shape = invalid_shape;
    }
    SECTION("unsupported dtype") {
        candidate.header.dtype = static_cast<CheckpointDtype>(255);
        expected_error = FailureCode::CheckpointIncompatible;
    }
    SECTION("cancellation precedes reading the payload") {
        REQUIRE(stop.request_stop());
        expected_error = FailureCode::Cancelled;
    }
    ResourcePlan ingress{16};
    const auto result = CheckpointChunk::authenticate(candidate, limits, ingress, stop.get_token());
    REQUIRE_FALSE(result);
    CHECK(result.error().code == expected_error);
    CHECK(ingress.allocated_bytes() == 0);
}

TEST_CASE("Chunk limits admit scalar and empty shapes and reject oversized products",
          "[IO4-06][N4-02]") {
    constexpr std::array<std::int64_t, 1> scalar{9};
    constexpr std::array<Index, 2> empty{std::numeric_limits<Index>::max(), 0};
    ResourcePlan plan{24};
    const auto one =
        CheckpointChunk::encode({}, std::span<const std::int64_t>{scalar}, {0, 8}, plan, {});
    REQUIRE(one);
    CHECK(chunk_view(*one).header.shape.empty());
    const auto none = CheckpointChunk::encode(empty, std::span<const double>{}, {2, 0}, plan, {});
    REQUIRE(none);
    CHECK(chunk_view(*none).bytes.empty());
    auto zero = none->decode(plan, {});
    REQUIRE(zero);
    CHECK(std::get<uvlas::core::PlannedBuffer<double>>(*zero).values().empty());
    CHECK(plan.allocated_bytes() == 24);

    auto shape = empty;
    auto requested_limits = limits;
    std::stop_source cancelled;
    auto expected_error = FailureCode::ResourceExceeded;
    SECTION("element count overflows") {
        shape[1] = 2;
    }
    SECTION("byte count overflows") {
        shape[1] = 1;
    }
    SECTION("rank exceeds limit even for a zero payload") {
        requested_limits.maximum_rank = 1;
    }
    SECTION("byte limit is below the payload") {
        shape = {1, 1};
        requested_limits.maximum_payload_bytes = 7;
    }
    SECTION("aggregate plan is one byte too small") {
        shape = {1, 1};
    }
    SECTION("shape and typed value count disagree") {
        shape = {1, 2};
        expected_error = FailureCode::InvalidInput;
    }
    SECTION("cancellation precedes allocation") {
        REQUIRE(cancelled.request_stop());
        expected_error = FailureCode::Cancelled;
    }
    ResourcePlan small{23};
    const auto rejected = CheckpointChunk::encode(shape, std::span<const std::int64_t>{scalar},
                                                  requested_limits, small, cancelled.get_token());
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == expected_error);
    CHECK(small.allocated_bytes() == 0);
}

TEST_CASE("Chunk hashes and ownership survive a partial final block", "[IO4-06][N4-02]") {
    constexpr Index count = 8193;
    constexpr auto payload_bytes = static_cast<std::uint64_t>(count) * 8;
    std::array<Index, 1> shape{count};
    std::vector<std::byte> source(payload_bytes);
    const auto golden =
        from_hex("9e0e77a31a045874b27217e09ec7690038c97f4c752c724a8d3ea1670c13ca6d");
    CheckpointChunkHeader header{CheckpointDtype::U64, shape, payload_bytes, {}};
    std::ranges::copy(golden, header.digest.bytes.begin());
    ResourcePlan plan{2 * payload_bytes + 8};
    auto verified = CheckpointChunk::authenticate({header, source}, limits, plan, {});
    REQUIRE(verified);
    shape[0] = 0;
    source.back() = std::byte{1};
    CHECK(chunk_view(*verified).header.shape[0] == count);
    CHECK(chunk_view(*verified).bytes.back() == std::byte{0});
    ResourcePlan insufficient{payload_bytes - 1};
    const auto rejected = verified->decode(insufficient, {});
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == FailureCode::ResourceExceeded);
    std::stop_source stopped;
    REQUIRE(stopped.request_stop());
    const auto cancelled = verified->decode(plan, stopped.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == FailureCode::Cancelled);
    CHECK(plan.allocated_bytes() == payload_bytes + 8);
    auto decoded = verified->decode(plan, {});
    REQUIRE(decoded);
    const auto decoded_values =
        std::get<uvlas::core::PlannedBuffer<std::uint64_t>>(*decoded).values();
    CHECK(decoded_values.size() == static_cast<std::size_t>(count));
    CHECK(std::ranges::all_of(decoded_values, [](auto value) { return value == 0; }));
}

} // namespace
