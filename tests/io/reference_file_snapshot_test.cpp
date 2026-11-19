#include "../support/allocation_audit.hpp"
#include "reference_file_snapshot.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <string_view>
#include <utility>

namespace {

using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using namespace uvlas::test_support;
using Reason = ReferenceFileFailure;
using Diagnostic = std::optional<ReferenceFileDiagnostic>;
constexpr auto operation = core::Operation::PrepareModel;

consteval SourceDigest digest_literal(std::string_view text) {
    constexpr std::string_view digits = "0123456789abcdef";
    SourceDigest digest{};
    for (std::size_t i = 0; i < digest.bytes.size(); ++i) {
        digest.bytes[i] =
            static_cast<std::byte>(16 * digits.find(text[2 * i]) + digits.find(text[2 * i + 1]));
    }
    return digest;
}

// Independent hashlib SHA-256 over the original byte sequences.
constexpr auto binary_hash =
    digest_literal("e9489f37fb3051e9efa1dc916004d7274e7b63975e3209708947267f2393a9be");
constexpr auto empty_hash =
    digest_literal("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
constexpr auto large_hash =
    digest_literal("ac72112c832fa4683b15ebff51a8f5f2ca08226c0d59bdb9ac739c2cdc28a05c");
constexpr std::array binary{std::byte{0}, std::byte{255}, std::byte{13}, std::byte{10}};
constexpr ReferenceFileLimits limits{4, 5, 5, 4};
constexpr ReferenceFileInfo binary_info{{42, {std::byte{3}}}, binary.size(), {1, 2, 3, 4}, true};

struct Reader final : ReferenceFileReader {
    Reader(std::span<const std::byte> payload, ReferenceFileInfo metadata) noexcept
        : input(payload), first(metadata), last(metadata) {}

    [[nodiscard]] core::Result<ReferenceFileInfo> info() noexcept override {
        ++info_calls;
        if (info_calls == failing_info) {
            return std::unexpected(
                core::Failure{core::FailureCode::Unsupported, core::Operation::Recover});
        }
        return info_calls == 1 ? first : last;
    }

    [[nodiscard]] core::Result<std::size_t> read(std::span<std::byte> out) noexcept override {
        ++read_calls;
        if (cancel_source != nullptr && read_calls == cancel_read) {
            cancel_source->request_stop();
        }
        if (read_calls == failing_read) {
            return std::unexpected(
                core::Failure{core::FailureCode::InvalidInput, core::Operation::Recover});
        }
        if (invalid_count) {
            return out.size() + 1;
        }
        const auto count = std::min({out.size(), input.size() - at, packet});
        std::ranges::copy(input.subspan(at, count), out.begin());
        at += count;
        return count;
    }

    const std::span<const std::byte> input;
    ReferenceFileInfo first;
    ReferenceFileInfo last;
    std::size_t at = 0;
    std::size_t packet = std::numeric_limits<std::size_t>::max();
    std::size_t read_calls = 0;
    std::size_t info_calls = 0;
    std::size_t failing_read = 0;
    std::size_t failing_info = 0;
    std::size_t cancel_read = 0;
    std::stop_source* cancel_source = nullptr;
    bool invalid_count = false;
};

void require_failure(const core::Result<ReferenceFileSnapshot>& outcome, const Diagnostic& diag,
                     Reason reason, core::FailureCode code, std::uint64_t bytes_read) {
    REQUIRE_FALSE(outcome);
    CHECK(outcome.error().code == code);
    CHECK(outcome.error().operation == operation);
    if (diag) {
        CHECK(diag->reason == reason);
        CHECK(diag->bytes_read == bytes_read);
    } else {
        FAIL("Missing reference reading diagnostic");
    }
}

} // namespace

TEST_CASE("Reference reading owns the exact bytes across short reads and owner destruction",
          "[IO4-01][N4-02][verification][property]") {
    const auto get = [] {
        Reader reader(binary, binary_info);
        reader.packet = 1;
        core::ResourcePlan plan{binary.size()};
        Diagnostic diag;
        AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
        auto result = read_reference_file(reader, binary_hash, limits, plan, {}, &diag);
        const auto allocation = audit.finish();
        REQUIRE(result);
        CHECK(allocation.cpp_allocation_attempts == 1);
        CHECK_FALSE(diag);
        CHECK(reader.info_calls == 2);
        CHECK(reader.read_calls == 5);
        CHECK(plan.allocated_bytes() == 4);
        return std::move(*result);
    };
    auto snapshot = get();
    const auto moved = std::move(snapshot);
    // NOLINTNEXTLINE(bugprone-use-after-move): the moved-from owner has an explicitly empty view.
    CHECK_FALSE(snapshot.view());
    if (const auto view = moved.view()) {
        CHECK(std::ranges::equal(view->bytes, binary));
        CHECK(view->digest == binary_hash);
        CHECK(view->file_info == binary_info);
    } else {
        FAIL("Expected the verified snapshot after releasing its reader and plan");
    }
}

TEST_CASE("Reference reading handles zero length and crosses the fixed read chunk boundary",
          "[IO4-01][N4-02][verification]") {
    const std::string large(65537, 'A');
    const std::array payloads{std::span<const std::byte>{}, std::as_bytes(std::span{large})};
    constexpr std::array digests{empty_hash, large_hash};
    constexpr std::array<std::size_t, 2> calls{1, 3};
    for (std::size_t i = 0; i < payloads.size(); ++i) {
        const auto payload = payloads[i];
        auto info = binary_info;
        info.byte_count = payload.size();
        Reader reader(payload, info);
        core::ResourcePlan plan{payload.size()};
        const auto result = read_reference_file(
            reader, digests[i], {payload.size(), payload.size() + 1, calls[i], payload.size()},
            plan, {}, nullptr);
        REQUIRE(result);
        if (const auto view = result->view()) {
            CHECK(std::ranges::equal(view->bytes, payload));
            CHECK(view->digest == digests[i]);
        } else {
            FAIL("Missing verified view");
        }
        CHECK(reader.read_calls == calls[i]);
        CHECK(plan.allocated_bytes() == payload.size());
    }
}

TEST_CASE("Reference reading rejects evidence that does not describe the bytes on one handle",
          "[IO4-01][negative][verification]") {
    constexpr std::array trailing{binary[0], binary[1], binary[2], binary[3], std::byte{42}};
    for (const auto reason :
         {Reason::Metadata, Reason::NotRegular, Reason::Read, Reason::InvalidReadCount,
          Reason::PrematureEof, Reason::TrailingByte, Reason::Modified, Reason::DigestMismatch}) {
        CAPTURE(static_cast<int>(reason));
        const auto payload = reason == Reason::PrematureEof   ? std::span{binary}.first(3)
                             : reason == Reason::TrailingByte ? std::span<const std::byte>{trailing}
                                                              : std::span<const std::byte>{binary};
        Reader reader(payload, binary_info);
        reader.packet = 1;
        auto expected_offset = std::uint64_t{4};
        if (reason == Reason::Metadata) {
            reader.failing_info = 2;
        } else if (reason == Reason::NotRegular) {
            reader.first.regular_file = false;
            expected_offset = 0;
        } else if (reason == Reason::Read) {
            reader.failing_read = 2;
            expected_offset = 1;
        } else if (reason == Reason::InvalidReadCount) {
            reader.invalid_count = true;
            expected_offset = 0;
        } else if (reason == Reason::PrematureEof) {
            expected_offset = 3;
        } else if (reason == Reason::Modified) {
            ++reader.last.change_indicators[2];
        }
        const auto expected_digest = reason == Reason::DigestMismatch ? empty_hash : binary_hash;
        const auto code = reason == Reason::Metadata ? core::FailureCode::Unsupported
                                                     : core::FailureCode::InvalidInput;
        core::ResourcePlan plan{4};
        Diagnostic diag;
        const auto result = read_reference_file(reader, expected_digest, limits, plan, {}, &diag);
        require_failure(result, diag, reason, code, expected_offset);
    }

    for (const auto& changed : {ReferenceFileInfo{{43, {std::byte{3}}}, 4, {1, 2, 3, 4}, true},
                                ReferenceFileInfo{binary_info.identity, 3, {1, 2, 3, 4}, true}}) {
        Reader reader(binary, binary_info);
        reader.last = changed;
        core::ResourcePlan plan{4};
        Diagnostic diag;
        const auto result = read_reference_file(reader, binary_hash, limits, plan, {}, &diag);
        require_failure(result, diag, Reason::Modified, core::FailureCode::InvalidInput, 4);
        if (diag) {
            CHECK(diag->before == binary_info);
            CHECK(diag->after == changed);
            CHECK_FALSE(diag->quota);
        }
    }
}

TEST_CASE("Reference reading admits capacity before allocation and meters extra short reads",
          "[IO4-01][N4-02][N4-09][negative][property]") {
    struct Case {
        ReferenceFileLimits cap;
        std::uint64_t plan_bytes;
        Reason reason;
        std::uint64_t required;
        std::uint64_t limit;
    };
    constexpr std::array cases{
        Case{{3, 5, 5, 4}, 4, Reason::FileBytes, 4, 3},
        Case{{4, 4, 5, 4}, 4, Reason::ReadBytes, 5, 4},
        Case{{4, 5, 1, 4}, 4, Reason::ReadCalls, 2, 1},
        Case{{4, 5, 5, 3}, 4, Reason::HashBytes, 4, 3},
        Case{{4, 5, 5, 4}, 3, Reason::PlanBytes, 4, 3},
    };
    for (const auto& test : cases) {
        Reader reader(binary, binary_info);
        core::ResourcePlan plan{test.plan_bytes};
        Diagnostic diag;
        AllocationAudit audit{cpp_backend, AuditThreadDomain::AllThreads};
        const auto result = read_reference_file(reader, binary_hash, test.cap, plan, {}, &diag);
        const auto allocation = audit.finish();
        require_failure(result, diag, test.reason, core::FailureCode::ResourceExceeded, 0);
        if (diag) {
            CHECK(diag->before == binary_info);
            CHECK_FALSE(diag->after);
            if (const auto quota = diag->quota) {
                CHECK(quota->required == test.required);
                CHECK(quota->limit == test.limit);
            } else {
                FAIL("Missing rejected resource bound");
            }
        }
        CHECK(allocation.status == AllocationAuditStatus::Clean);
        CHECK(allocation.cpp_allocation_attempts == 0);
        CHECK(plan.allocated_bytes() == 0);
        CHECK(reader.read_calls == 0);
    }
    Reader short_reads(binary, binary_info);
    short_reads.packet = 1;
    core::ResourcePlan plan{4};
    Diagnostic diag;
    const auto capped =
        read_reference_file(short_reads, binary_hash, {4, 5, 4, 4}, plan, {}, &diag);
    require_failure(capped, diag, Reason::ReadCalls, core::FailureCode::ResourceExceeded, 4);
    CHECK(short_reads.read_calls == 4);
    CHECK(plan.allocated_bytes() == 4);
}

TEST_CASE("Cancellation before and during reference reading cannot publish a payload",
          "[IO4-01][N4-02][negative]") {
    Reader reader(binary, binary_info);
    std::stop_source source;
    reader.cancel_source = &source;
    reader.cancel_read = 2; // Cancel as the adapter returns the terminal EOF probe.
    core::ResourcePlan plan{4};
    Diagnostic diag = ReferenceFileDiagnostic{Reason::DigestMismatch, 0};
    const auto result =
        read_reference_file(reader, binary_hash, limits, plan, source.get_token(), &diag);
    REQUIRE_FALSE(result);
    CHECK(result.error().code == core::FailureCode::Cancelled);
    CHECK_FALSE(diag);
    CHECK(reader.info_calls == 1);
    CHECK(plan.allocated_bytes() == 4);

    Reader skipped(binary, binary_info);
    core::ResourcePlan next{4};
    const auto pre =
        read_reference_file(skipped, binary_hash, limits, next, source.get_token(), nullptr);
    REQUIRE_FALSE(pre);
    CHECK(pre.error().code == core::FailureCode::Cancelled);
    CHECK(skipped.info_calls == 0);
    CHECK(next.allocated_bytes() == 0);
}
