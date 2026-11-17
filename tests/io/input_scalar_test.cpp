#include <catch2/catch_test_macros.hpp>

#include <array>
#include <limits>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <uvlas/io/input_scalar.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using Diagnostic = std::optional<ScalarConversionDiagnostic>;
using Reason = ScalarConversionReason;
using Ids = std::tuple<core::DofId, core::PanelId, core::WakeVertexId, core::WakeFaceId,
                       core::LayerId, core::ObservationId, core::ScopeId, core::ComponentId>;
struct SourceTag;
using SourceId = core::Identifier<SourceTag>;

template <class T>
concept ReadableId = requires(std::string_view text, Diagnostic& diagnostic) {
    read_physical_id<T>(text, diagnostic);
};
template <class Left, class... Right>
constexpr bool distinct_domains(std::tuple<Right...>*) {
    return ((std::same_as<Left, Right> || !std::is_convertible_v<Left, Right>) && ...);
}
template <class... T>
constexpr bool type_contract(std::tuple<T...>* types) {
    return ((ReadableId<T> && distinct_domains<T>(types) &&
             !std::is_convertible_v<T, core::Index> &&
             !std::is_constructible_v<NativeInteger, T>) &&
            ...);
}
static_assert(type_contract(static_cast<Ids*>(nullptr)));
static_assert(!ReadableId<core::Index> && !ReadableId<std::uint64_t> && !ReadableId<SourceId>);
static_assert(!ReadableId<core::ModelGeneration> && !ReadableId<core::NumericalGeneration>);
static_assert(!ReadableId<core::TrialId> && !ReadableId<core::AcceptedStateId>);
static_assert(!std::is_default_constructible_v<NativeInteger>);
static_assert(std::is_constructible_v<NativeInteger, std::int64_t>);
static_assert(std::is_constructible_v<NativeInteger, std::uint64_t>);
static_assert(!std::is_convertible_v<std::int64_t, NativeInteger>);
static_assert(!std::is_convertible_v<std::uint64_t, NativeInteger>);
static_assert(!std::is_constructible_v<NativeInteger, bool>);
static_assert(!std::is_constructible_v<NativeInteger, float>);
static_assert(!std::is_constructible_v<NativeInteger, double>);
static_assert(!std::is_constructible_v<NativeInteger, std::string_view>);
static_assert(!std::is_constructible_v<NativeInteger, std::int32_t>);
static_assert(!std::is_constructible_v<NativeInteger, std::uint32_t>);
static_assert(std::is_nothrow_copy_constructible_v<ScalarConversionDiagnostic>);
static_assert(noexcept(read_physical_id<core::PanelId>(std::string_view{},
                                                       std::declval<Diagnostic&>())));
static_assert(noexcept(read_index(NativeInteger{std::int64_t{0}}, std::declval<Diagnostic&>())));

[[nodiscard]] core::FailureDetails context() {
    return {core::ModelGeneration{1},
            core::AcceptedStateId{2},
            core::TrialId{3},
            core::NumericalGeneration{4},
            core::ComponentId{5},
            core::ProblemId{core::LayerId{6}},
            core::ResourceUsage{core::ResourceKind::Bytes, 8, 7},
            core::AcceptedStateId{9}};
}
void check_diagnostic(const Diagnostic& diagnostic, Reason reason,
                      std::optional<std::size_t> position) {
    if (!diagnostic.has_value()) {
        FAIL("Missing owned conversion diagnostic");
        return;
    }
    CHECK(diagnostic->reason == reason);
    CHECK(diagnostic->decoded_byte_position == position);
}
template <class T>
void check_failure(const core::Result<T>& result, const Diagnostic& diagnostic, Reason reason,
                   std::optional<std::size_t> position) {
    if (result.has_value()) {
        FAIL("Unexpected scalar conversion success");
        return;
    }
    CHECK(result.error().code == core::FailureCode::InvalidInput);
    CHECK(result.error().operation == core::Operation::PrepareModel);
    CHECK(result.error().details == context());
    check_diagnostic(diagnostic, reason, position);
}
struct ValidId {
    std::string_view text;
    std::uint64_t expected;
};
// Exact decimal literals are independent of the product accumulator; tolerance is zero.
constexpr std::array valid_ids{
    ValidId{"0", 0},
    ValidId{"1", 1},
    ValidId{"18446744073709551615", 18446744073709551615ULL},
    ValidId{"9007199254740991", 9007199254740991ULL},
    ValidId{"9007199254740992", 9007199254740992ULL},
    ValidId{"9007199254740993", 9007199254740993ULL},
    ValidId{"00", 0},
    ValidId{"0001", 1},
    ValidId{"018446744073709551615", 18446744073709551615ULL},
    ValidId{"0000000000000000000000000000000000000000000000000000000000000000", 0}};
struct InvalidId {
    std::string_view text;
    Reason reason;
    std::optional<std::size_t> position;
};
constexpr std::array invalid_ids{
    InvalidId{"", Reason::EmptyId, std::nullopt},
    InvalidId{"18446744073709551616", Reason::IdOverflow, 19},
    InvalidId{"018446744073709551616", Reason::IdOverflow, 20},
    InvalidId{"99999999999999999999", Reason::IdOverflow, 19},
    InvalidId{"18446744073709551616X", Reason::IdOverflow, 19},
    InvalidId{"18446744073709551615X", Reason::InvalidIdCharacter, 20},
    InvalidId{"+1", Reason::InvalidIdCharacter, 0},
    InvalidId{"-1", Reason::InvalidIdCharacter, 0},
    InvalidId{"-0", Reason::InvalidIdCharacter, 0},
    InvalidId{" 1", Reason::InvalidIdCharacter, 0},
    InvalidId{"\t1", Reason::InvalidIdCharacter, 0},
    InvalidId{"\r1", Reason::InvalidIdCharacter, 0},
    InvalidId{"\n1", Reason::InvalidIdCharacter, 0},
    InvalidId{"1 ", Reason::InvalidIdCharacter, 1},
    InvalidId{"1\t", Reason::InvalidIdCharacter, 1},
    InvalidId{"1\r", Reason::InvalidIdCharacter, 1},
    InvalidId{"1\n", Reason::InvalidIdCharacter, 1},
    InvalidId{"1.0", Reason::InvalidIdCharacter, 1},
    InvalidId{"1e0", Reason::InvalidIdCharacter, 1},
    InvalidId{"1E0", Reason::InvalidIdCharacter, 1},
    InvalidId{"1_0", Reason::InvalidIdCharacter, 1},
    InvalidId{"1,0", Reason::InvalidIdCharacter, 1},
    InvalidId{"\xef\xbc\x91", Reason::InvalidIdCharacter, 0},
    InvalidId{"\xd9\xa1", Reason::InvalidIdCharacter, 0},
    InvalidId{"1\xef\xbc\x91", Reason::InvalidIdCharacter, 1},
    InvalidId{"1\xd9\xa1", Reason::InvalidIdCharacter, 1},
    InvalidId{std::string_view{"\0", 1}, Reason::InvalidIdCharacter, 0},
    InvalidId{std::string_view{"1\0"
                               "2",
                               3},
              Reason::InvalidIdCharacter, 1},
    InvalidId{std::string_view{"1\0", 2}, Reason::InvalidIdCharacter, 1}};

template <class Id>
void check_id_domain() {
    for (const auto& fixture : valid_ids) {
        Diagnostic diagnostic{ScalarConversionDiagnostic{Reason::IndexOverflow, std::nullopt}};
        const auto result = read_physical_id<Id>(fixture.text, diagnostic, context());
        if (!result.has_value()) {
            FAIL("Rejected fixed valid ID");
            return;
        }
        CHECK(result->value() == fixture.expected);
        CHECK_FALSE(diagnostic.has_value());
    }
    for (const auto& fixture : invalid_ids) {
        Diagnostic diagnostic;
        check_failure(read_physical_id<Id>(fixture.text, diagnostic, context()), diagnostic,
                      fixture.reason, fixture.position);
    }
}
} // namespace

TEST_CASE("Decoded identifiers preserve each physical domain and every fixed integer",
          "[IO4-01][N4-02][N4-08][verification][negative]") {
    check_id_domain<core::DofId>();
    check_id_domain<core::PanelId>();
    check_id_domain<core::WakeVertexId>();
    check_id_domain<core::WakeFaceId>();
    check_id_domain<core::LayerId>();
    check_id_domain<core::ObservationId>();
    check_id_domain<core::ScopeId>();
    check_id_domain<core::ComponentId>();
}

TEST_CASE("Native signed and unsigned indices have fixed strict boundaries",
          "[IO4-01][N4-02][N4-08][verification][negative]") {
    struct ValidIndex {
        NativeInteger input;
        core::Index expected;
        ValidIndex(NativeInteger value, core::Index answer) : input(value), expected(answer) {}
    };
    const std::array valid{
        ValidIndex{NativeInteger{std::int64_t{0}}, 0},
        ValidIndex{NativeInteger{std::int64_t{9223372036854775807LL}}, 9223372036854775807LL},
        ValidIndex{NativeInteger{std::uint64_t{0}}, 0},
        ValidIndex{NativeInteger{std::uint64_t{9223372036854775807ULL}}, 9223372036854775807LL}};
    for (const auto& fixture : valid) {
        Diagnostic diagnostic{ScalarConversionDiagnostic{Reason::EmptyId, std::nullopt}};
        const auto result = read_index(fixture.input, diagnostic, context());
        if (!result.has_value()) {
            FAIL("Rejected fixed valid index");
            return;
        }
        CHECK(result.value() == fixture.expected);
        CHECK_FALSE(diagnostic.has_value());
    }
    struct InvalidIndex {
        NativeInteger input;
        Reason reason;
        InvalidIndex(NativeInteger value, Reason cause) : input(value), reason(cause) {}
    };
    const std::array invalid{
        InvalidIndex{NativeInteger{std::numeric_limits<std::int64_t>::min()},
                     Reason::NegativeIndex},
        InvalidIndex{NativeInteger{std::int64_t{-1}}, Reason::NegativeIndex},
        InvalidIndex{NativeInteger{std::uint64_t{9223372036854775808ULL}}, Reason::IndexOverflow},
        InvalidIndex{NativeInteger{std::uint64_t{18446744073709551615ULL}}, Reason::IndexOverflow}};
    for (const auto& fixture : invalid) {
        Diagnostic diagnostic;
        check_failure(read_index(fixture.input, diagnostic, context()), diagnostic, fixture.reason,
                      std::nullopt);
    }
    // A native integer has no lexical minus-zero bit; the document retains that spelling.
    Diagnostic diagnostic;
    CHECK(read_index(NativeInteger{std::int64_t{0}}, diagnostic).value() == 0);
    CHECK(std::holds_alternative<std::int64_t>(NativeInteger{std::int64_t{0}}.selection()));
    CHECK(std::holds_alternative<std::uint64_t>(NativeInteger{std::uint64_t{0}}.selection()));
}

TEST_CASE("Scalar failures own diagnostics and preserve inputs and prior values",
          "[IO4-01][N4-03][N4-08][property][negative]") {
    Diagnostic diagnostic;
    std::string text = "0001";
    const auto old_text = text;
    const auto old_result = read_physical_id<core::PanelId>(text, diagnostic, context());
    REQUIRE(old_result.has_value());
    CHECK(text == old_text);
    text = "18446744073709551616";
    const auto invalid_text = text;
    const auto failed = read_physical_id<core::PanelId>(text, diagnostic, context());
    check_failure(failed, diagnostic, Reason::IdOverflow, 19);
    CHECK(text == invalid_text);
    const auto owned_diagnostic = diagnostic;
    const auto owned_failure = failed;
    text.clear();
    text.shrink_to_fit();
    check_failure(owned_failure, owned_diagnostic, Reason::IdOverflow, 19);
    CHECK(old_result->value() == 1);
    const auto next = read_physical_id<core::PanelId>("2", diagnostic);
    REQUIRE(next.has_value());
    CHECK_FALSE(diagnostic.has_value());
    check_diagnostic(owned_diagnostic, Reason::IdOverflow, 19);
    CHECK(old_result->value() == 1);
}

TEST_CASE("Independent scalar calls share only immutable input",
          "[IO4-01][N4-03][N4-08][property]") {
    constexpr std::string_view text = "9007199254740993";
    std::array<std::optional<core::Result<core::PanelId>>, 2> results;
    std::array<Diagnostic, 2> diagnostics;
    {
        const std::jthread first{
            [&] { results[0] = read_physical_id<core::PanelId>(text, diagnostics[0]); }};
        const std::jthread second{
            [&] { results[1] = read_physical_id<core::PanelId>(text, diagnostics[1]); }};
    }
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];
        if (!result.has_value() || !result->has_value()) {
            FAIL("Independent scalar call failed");
            return;
        }
        CHECK(result->value().value() == 9007199254740993ULL);
        CHECK_FALSE(diagnostics[i].has_value());
    }
}
