#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <uvlas/io/input_document.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace std::literals;
constexpr std::uint64_t budget = 64ULL * 1024 * 1024;
constexpr std::string_view minimal =
    R"({"schema":"uvlas.case/1","units":null,"model":null,"flow":null,"analysis":null,"observations":null,"numerics":null,"resources":null,"output":null})";
InputLimits limits() {
    return {100000,
            64,
            10000,
            10000,
            10000,
            10000,
            100000,
            {ParserAllocationMode::EstimateAdmission, std::numeric_limits<std::uint64_t>::max()}};
}
template <class T>
T take(core::Result<T> value) {
    if (!value) {
        FAIL_CHECK("Unexpected document failure " << static_cast<unsigned>(value.error().code));
        throw std::runtime_error("Unexpected document failure");
    }
    return std::move(*value);
}
core::Result<InputDocument> read(std::string_view source, const InputLimits& request,
                                 core::ResourcePlan& plan,
                                 std::optional<InputDiagnostic>& diagnostic,
                                 SourceId id = SourceId{0}, core::FailureDetails details = {}) {
    return read_input_document(std::as_bytes(std::span{source.data(), source.size()}), id, request,
                               plan, diagnostic, details);
}
InputDocument document(std::string_view source) {
    core::ResourcePlan plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    auto value = take(read(source, limits(), plan, diagnostic));
    CHECK_FALSE(diagnostic);
    return value;
}
std::string envelope(std::string_view model) {
    return std::string{R"({"schema":"uvlas.case/1","units":null,"model":)"} + std::string{model} +
           R"(,"flow":null,"analysis":null,"observations":null,"numerics":null,"resources":null,"output":null})";
}
InputNodeView model(const InputDocument& value) {
    const auto view = take(value.view());
    const auto root = view.root.object();
    if (!root) {
        FAIL_CHECK("Missing root object");
        throw std::runtime_error("Missing root object");
    }
    const auto found = root->find("model");
    if (!found) {
        FAIL_CHECK("Missing model node");
        throw std::runtime_error("Missing model node");
    }
    return *found;
}
void rejected(std::string_view text, core::FailureCode code,
              std::optional<InputDiagnosticReason> reason = std::nullopt) {
    const std::string preserved{text};
    core::ResourcePlan plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    const auto value = read(text, limits(), plan, diagnostic);
    REQUIRE_FALSE(value);
    CHECK(value.error().code == code);
    CHECK(value.error().operation == core::Operation::PrepareModel);
    if (!diagnostic) {
        FAIL_CHECK("Missing failure diagnostic");
        return;
    }
    if (reason)
        CHECK(diagnostic->reason == *reason);
    CHECK(text == preserved);
}
std::string hex(SourceDigest digest) {
    constexpr std::string_view alphabet = "0123456789abcdef";
    std::string output;
    for (const auto byte : digest.bytes) {
        const auto value = std::to_integer<unsigned>(byte);
        output.push_back(alphabet[value >> 4U]);
        output.push_back(alphabet[value & 15U]);
    }
    return output;
}
} // namespace

TEST_CASE("Owned documents preserve all value kinds exact tokens and nested direct order",
          "[verification][IO4-01][N4-03][N4-08]") {
    const auto value = document(envelope(
        R"([null,true,false,-9223372036854775808,-0,0,9223372036854775807,18446744073709551615,18446744073709551616,-9223372036854775809,1.0,1e0,1e-400,4.9406564584124654e-324,"NaN","a\u0000b","日本語","\uD83D\uDE00",{"a":[1,{"b":2}],"c":3},[],{}])"));
    const auto node = model(value);
    const auto array = node.array();
    if (!array) {
        FAIL_CHECK("Missing expected array");
        return;
    }
    REQUIRE(array->size() == 21);
    const auto null = take(array->at(0));
    CHECK(null.kind() == InputNodeKind::Null);
    CHECK_FALSE(null.boolean());
    CHECK_FALSE(null.string());
    CHECK_FALSE(null.native_integer());
    CHECK(take(array->at(1)).boolean() == true);
    CHECK(take(array->at(2)).boolean() == false);
    CHECK(take(array->at(3)).signed_integer() == std::numeric_limits<std::int64_t>::min());
    CHECK(take(array->at(4)).signed_integer() == 0);
    CHECK(take(array->at(5)).unsigned_integer() == 0);
    CHECK(take(array->at(6)).unsigned_integer() == UINT64_C(9223372036854775807));
    CHECK(take(array->at(7)).unsigned_integer() == std::numeric_limits<std::uint64_t>::max());
    constexpr std::array<std::string_view, 6> raw{
        "18446744073709551616",   "-9223372036854775809", "1.0", "1e0", "1e-400",
        "4.9406564584124654e-324"};
    for (std::size_t i = 0; i < raw.size(); ++i) {
        const auto decimal = take(array->at(static_cast<core::Index>(i + 8)));
        CHECK(decimal.kind() == InputNodeKind::DecimalNumber);
        CHECK(decimal.decimal_token() == raw[i]);
        CHECK_FALSE(decimal.native_integer());
    }
    const auto word = take(array->at(14)), nul = take(array->at(15)),
               japanese = take(array->at(16)), emoji = take(array->at(17));
    CHECK(word.string() == "NaN");
    CHECK(nul.string() == "a\0b"sv);
    CHECK(japanese.string() == "日本語");
    CHECK(emoji.string() == "\xf0\x9f\x98\x80");
    const auto outer = take(array->at(18));
    const auto object = outer.object();
    if (!object) {
        FAIL_CHECK("Missing expected object");
        return;
    }
    REQUIRE(object->size() == 2);
    const auto a = take(object->at(0)), c = take(object->at(1));
    CHECK(a.key == "a");
    CHECK(c.key == "c");
    CHECK(c.value.unsigned_integer() == 3);
    const auto child_array = a.value.array();
    if (!child_array) {
        FAIL_CHECK("Missing expected child array");
        return;
    }
    CHECK(child_array->size() == 2);
    const auto nested = take(child_array->at(1));
    const auto nested_object = nested.object();
    if (!nested_object) {
        FAIL_CHECK("Missing expected nested object");
        return;
    }
    CHECK(take(nested_object->at(0)).value.unsigned_integer() == 2);
    CHECK(c.value.ordinal().value() > nested.ordinal().value());
    CHECK_FALSE(object->find("missing"));
    CHECK_FALSE(array->at(-1));
    CHECK_FALSE(array->at(21));
    const auto empty_array = take(array->at(19)), empty_object = take(array->at(20));
    const auto empty_a = empty_array.array();
    const auto empty_o = empty_object.object();
    if (!empty_a || !empty_o) {
        FAIL_CHECK("Missing empty containers");
        return;
    }
    CHECK(empty_a->size() == 0);
    CHECK(empty_o->size() == 0);
}

TEST_CASE("Strict JSON and UTF8 reject every malformed literal without accepting a prefix",
          "[negative][verification][IO4-01][N4-08]") {
    using Code = core::FailureCode;
    rejected("", Code::InvalidInput, InputDiagnosticReason::EmptyInput);
    for (const std::string_view invalid :
         {" ", "null", "true", "1", "[]", "{}", "[1,]", "[/*x*/1]", "[//x\n1]", "{", "["})
        rejected(invalid, Code::InvalidInput);
    for (const std::string_view invalid :
         {R"("\uD800")", R"("\uDC00")", R"("\uD800\u0000")", R"("\x00")", R"("\u12G4")", "NaN",
          "Infinity", "+1", "01", "1.", ".1"})
        rejected(envelope(invalid), Code::InvalidInput);
    rejected(envelope("1e400"), Code::Nonfinite, InputDiagnosticReason::NumberRange);
    rejected(envelope("-1e400"), Code::Nonfinite, InputDiagnosticReason::NumberRange);
    rejected(std::string{minimal} + " false", Code::InvalidInput,
             InputDiagnosticReason::JsonSyntax);
    rejected(std::string{"\xef\xbb\xbf"} + std::string{minimal}, Code::InvalidInput,
             InputDiagnosticReason::ByteOrderMark);
    for (const auto& bytes : {std::string{"\xc0\xaf", 2}, std::string{"\xed\xa0\x80", 3},
                              std::string{"\xf4\x90\x80\x80", 4}, std::string{"\xe3\x81", 2}})
        rejected(envelope("\"" + bytes + "\""), Code::InvalidInput,
                 InputDiagnosticReason::InvalidUtf8);
    for (const auto& bytes :
         {std::string{"\xff\xfe{\0}\0", 6}, std::string{"\xfe\xff\0{\0}", 6},
          std::string{"\0\0\xfe\xff", 4}, std::string{"\xff\xfe\0\0", 4}, std::string{"{\0}\0", 4}})
        rejected(bytes, Code::InvalidInput);
    rejected(std::string{minimal} + std::string{"\0false", 6}, Code::InvalidInput,
             InputDiagnosticReason::RawNull);
    rejected(envelope(std::string{"\"a\0b\"", 5}), Code::InvalidInput,
             InputDiagnosticReason::RawNull);
    const auto valid = document(" \n\t\r" + std::string{minimal} + "\r\n ");
    CHECK(take(valid.view()).identity.byte_count == minimal.size() + 7);
}

TEST_CASE("Parser failure reasons ignore diagnostic words copied from source text",
          "[verification][negative][IO4-01][N4-08]") {
    constexpr std::array<std::string_view, 7> prefixes{
        "hello",
        "ill-formed UTF-8",
        "surrogate",
        "4 hex digits",
        "after backslash",
        "; last read: 'ill-formed UTF-8 surrogate 4 hex digits after backslash",
        "; last read: 'surrogate'; last read: 'ill-formed UTF-8"};
    for (const auto prefix : prefixes) {
        CAPTURE(prefix);
        // The fixed lexer cause, not a word inside its last-read quotation,
        // determines the reason. Check both a model value and strict EOF scan.
        const auto bad_escape = "\"" + std::string{prefix} + R"(\q")";
        rejected(envelope(bad_escape), core::FailureCode::InvalidInput,
                 InputDiagnosticReason::InvalidStringEscape);
        rejected(std::string{minimal} + " " + bad_escape, core::FailureCode::InvalidInput,
                 InputDiagnosticReason::InvalidStringEscape);
        rejected(envelope("\"" + std::string{prefix} + "\n\""), core::FailureCode::InvalidInput,
                 InputDiagnosticReason::JsonSyntax);
        rejected(envelope("\"" + std::string{prefix} + "\xc0\xaf\""),
                 core::FailureCode::InvalidInput, InputDiagnosticReason::InvalidUtf8);
    }
    for (const auto suffix : {R"(\uZZZZ")", R"(\uD800")", R"(\uDC00")"})
        rejected(envelope("\"ill-formed UTF-8" + std::string{suffix}),
                 core::FailureCode::InvalidInput, InputDiagnosticReason::InvalidStringEscape);
    // A valid extra token takes the parser's unexpected-token path, which does
    // not append the token's user text to the generated diagnostic.
    rejected(std::string{minimal} + R"( "ill-formed UTF-8 surrogate 4 hex digits")",
             core::FailureCode::InvalidInput, InputDiagnosticReason::JsonSyntax);
}

TEST_CASE("Native integer nodes preserve the entire binary64 boundary before strict conversion",
          "[verification][negative][IO4-01][N4-03][N4-08]") {
    const auto value = document(envelope("[9007199254740991,9007199254740992,9007199254740993,"
                                         "9223372036854775807,18446744073709551615,-0]"));
    const auto root = model(value);
    const auto array = root.array();
    if (!array) {
        FAIL_CHECK("Missing integer series");
        return;
    }
    constexpr std::array<std::uint64_t, 3> expected{
        UINT64_C(9007199254740991), UINT64_C(9007199254740992), UINT64_C(9007199254740993)};
    std::optional<ScalarConversionDiagnostic> diagnostic;
    for (core::Index i = 0; i < 3; ++i) {
        const auto node = take(array->at(i));
        CHECK(node.kind() == InputNodeKind::UnsignedInteger);
        CHECK(node.unsigned_integer() == expected[static_cast<std::size_t>(i)]);
        const auto native = node.native_integer();
        if (!native) {
            FAIL_CHECK("Integer node lost its native kind");
            return;
        }
        CHECK(take(read_index(*native, diagnostic)) ==
              static_cast<core::Index>(expected[static_cast<std::size_t>(i)]));
        CHECK_FALSE(diagnostic);
    }
    const auto largest = take(array->at(3)), overflow = take(array->at(4)),
               negative_zero = take(array->at(5));
    const auto native_max = largest.native_integer(), native_overflow = overflow.native_integer(),
               native_zero = negative_zero.native_integer();
    if (!native_max || !native_overflow || !native_zero) {
        FAIL_CHECK("Native boundary kind missing");
        return;
    }
    CHECK(take(read_index(*native_max, diagnostic)) == std::numeric_limits<core::Index>::max());
    CHECK_FALSE(read_index(*native_overflow, diagnostic));
    CHECK(take(read_index(*native_zero, diagnostic)) == 0);
    const auto id_document = document(envelope(R"("9007199254740993")"));
    const auto id_node = model(id_document);
    const auto id_text = id_node.string();
    if (!id_text) {
        FAIL_CHECK("Missing physical ID string");
        return;
    }
    CHECK(take(read_physical_id<core::DofId>(*id_text, diagnostic)).value() ==
          UINT64_C(9007199254740993));
}

TEST_CASE("Decoded key uniqueness is local while object and array order remain observable",
          "[verification][negative][IO4-01][N4-03][N4-08]") {
    rejected(R"({"schema":"uvlas.case/1","\u0073chema":null})", core::FailureCode::InvalidInput,
             InputDiagnosticReason::DuplicateKey);
    rejected(envelope(R"({"k":1,"\u006b":NaN})"), core::FailureCode::InvalidInput,
             InputDiagnosticReason::DuplicateKey);
    const auto value = document(envelope(R"({"z":{"k":1},"a":{"k":2},"k":3,"K":4,"a\u0000b":5})"));
    const auto node = model(value);
    const auto object = node.object();
    if (!object) {
        FAIL_CHECK("Missing expected object");
        return;
    }
    constexpr std::array<std::string_view, 5> order{"z", "a", "k", "K", "a\0b"sv};
    for (core::Index i = 0; i < 5; ++i)
        CHECK(take(object->at(i)).key == order[static_cast<std::size_t>(i)]);
    const auto nul_key = object->find("a\0b"sv);
    if (!nul_key) {
        FAIL_CHECK("Missing decoded NUL key");
        return;
    }
    CHECK(nul_key->unsigned_integer() == 5);
    const auto uppercase = object->find("K");
    if (!uppercase) {
        FAIL_CHECK("Missing case-sensitive key");
        return;
    }
    CHECK(uppercase->unsigned_integer() == 4);
}

TEST_CASE("Root envelope validates schema then unknown keys then every mandatory key",
          "[verification][negative][IO4-01][N4-08]") {
    constexpr std::array<std::string_view, 9> required{"schema",   "units",     "model",
                                                       "flow",     "analysis",  "observations",
                                                       "numerics", "resources", "output"};
    for (const auto omitted : required) {
        std::string source = "{";
        for (const auto key : required) {
            if (key == omitted)
                continue;
            if (source.size() != 1)
                source += ',';
            source += '"';
            source += key;
            source += "\":";
            source += key == "schema" ? "\"uvlas.case/1\"" : "null";
        }
        source += '}';
        rejected(source, core::FailureCode::InvalidInput, InputDiagnosticReason::MissingRootKey);
    }
    for (const std::string_view schema : {"null", "true", "0", "[]", "{}"})
        rejected("{\"schema\":" + std::string{schema} + ",\"unknown\":0}",
                 core::FailureCode::InvalidInput, InputDiagnosticReason::SchemaType);
    rejected(R"({"schema":"uvlas.case/2","unknown":0})", core::FailureCode::Unsupported,
             InputDiagnosticReason::UnsupportedSchema);
    rejected(R"({"schema":"uvlas.case/1","unknown":0})", core::FailureCode::InvalidInput,
             InputDiagnosticReason::UnknownRootKey);
    struct UnknownCase {
        UnknownCase(std::string_view bytes, std::uint64_t event) : input(bytes), ordinal(event) {}
        std::string_view input;
        std::uint64_t ordinal;
    };
    const std::array unknown_order{UnknownCase{R"({"schema":"uvlas.case/1","z":0,"a":0})", 5},
                                   UnknownCase{R"({"schema":"uvlas.case/1","a":0,"z":0})", 3},
                                   UnknownCase{R"({"schema":"uvlas.case/1","\u007a":0,"a":0})", 5},
                                   UnknownCase{R"({"schema":"uvlas.case/1","a":0,"\u007a":0})", 3}};
    for (const auto& test : unknown_order) {
        core::ResourcePlan plan{budget};
        std::optional<InputDiagnostic> diagnostic;
        const auto failed = read(test.input, limits(), plan, diagnostic);
        REQUIRE_FALSE(failed);
        if (!diagnostic) {
            FAIL_CHECK("Missing unknown-key diagnostic");
            return;
        }
        CHECK(diagnostic->reason == InputDiagnosticReason::UnknownRootKey);
        CHECK(diagnostic->event_ordinal == test.ordinal); // Always decoded key a's event.
        CHECK_FALSE(diagnostic->byte_offset);
    }
    constexpr std::array<std::string_view, 3> optional{"initial", "restart", "provenance"};
    for (unsigned mask = 0; mask < 8; ++mask) {
        std::string source{minimal.substr(0, minimal.size() - 1)};
        for (unsigned i = 0; i < 3; ++i)
            if ((mask & (1U << i)) != 0)
                source += ",\"" + std::string{optional[i]} + "\":null";
        source += '}';
        const auto value = document(source);
        const auto view = take(value.view());
        const auto root = view.root.object();
        if (!root) {
            FAIL_CHECK("Missing expected root object");
            return;
        }
        CHECK(root->size() == 9 + static_cast<core::Index>(std::popcount(mask)));
    }
    // Nested physical schemas intentionally remain unqualified at this boundary.
    const auto unqualified = document(envelope(R"({"arbitrary_nested_key":true})"));
    CHECK(model(unqualified).kind() == InputNodeKind::Object);
}

TEST_CASE("Each document quota admits equality and rejects one less with owned quota detail",
          "[verification][negative][IO4-01][N4-02][N4-08]") {
    // Independently counted literal: 146 bytes, root+9 values, 9 keys, depth1,
    // maximum decoded key/value12, key bytes63 + schema bytes12 =75 text bytes.
    CHECK(minimal.size() == 146);
    const InputLimits exact{146, 1, 10, 9, 12, 0, 75, limits().external_memory};
    core::ResourcePlan good{budget};
    std::optional<InputDiagnostic> diagnostic;
    CHECK(read(minimal, exact, good, diagnostic));
    struct QuotaCase {
        InputQuota kind;
        std::uint64_t InputLimits::* field;
    };
    const std::array cases{
        QuotaCase{InputQuota::SourceBytes, &InputLimits::maximum_source_bytes},
        QuotaCase{InputQuota::Depth, &InputLimits::maximum_depth},
        QuotaCase{InputQuota::Nodes, &InputLimits::maximum_nodes},
        QuotaCase{InputQuota::Members, &InputLimits::maximum_members},
        QuotaCase{InputQuota::DecodedStringBytes, &InputLimits::maximum_decoded_string_bytes},
        QuotaCase{InputQuota::TotalTextBytes, &InputLimits::maximum_total_text_bytes}};
    for (const auto test : cases) {
        auto low = exact;
        --(low.*test.field);
        core::ResourcePlan plan{budget};
        const auto result = read(minimal, low, plan, diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::ResourceExceeded);
        if (!diagnostic || !diagnostic->quota) {
            FAIL_CHECK("Missing quota diagnostic");
            return;
        }
        CHECK(diagnostic->quota->kind == test.kind);
        CHECK(diagnostic->quota->limit == low.*test.field);
    }
    const auto numeric = envelope("1.25");
    auto decimal = exact;
    decimal.maximum_decimal_token_bytes = 4;
    decimal.maximum_total_text_bytes = 79;
    core::ResourcePlan equal{budget};
    CHECK(read(numeric, decimal, equal, diagnostic));
    decimal.maximum_decimal_token_bytes = 3;
    core::ResourcePlan low{budget};
    const auto result = read(numeric, decimal, low, diagnostic);
    REQUIRE_FALSE(result);
    if (!diagnostic || !diagnostic->quota) {
        FAIL_CHECK("Missing decimal quota diagnostic");
        return;
    }
    CHECK(diagnostic->quota->kind == InputQuota::DecimalTokenBytes);
    for (const auto test : cases) {
        auto zero = limits();
        zero.*test.field = 0;
        core::ResourcePlan plan{budget};
        CHECK_FALSE(read(minimal, zero, plan, diagnostic));
    }
}

TEST_CASE("Nested quotas count actual containers members and UTF8 bytes without normalization",
          "[verification][negative][IO4-01][N4-02][N4-08]") {
    const auto accepted = [](std::string_view input, const InputLimits& request) {
        core::ResourcePlan plan{budget};
        std::optional<InputDiagnostic> diagnostic;
        CHECK(read(input, request, plan, diagnostic));
        CHECK_FALSE(diagnostic);
    };
    const auto denied = [](std::string_view input, const InputLimits& request, InputQuota kind) {
        core::ResourcePlan plan{budget};
        std::optional<InputDiagnostic> diagnostic;
        const auto failed = read(input, request, plan, diagnostic);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().code == core::FailureCode::ResourceExceeded);
        if (!diagnostic || !diagnostic->quota) {
            FAIL_CHECK("Missing nested quota detail");
            return;
        }
        CHECK(diagnostic->quota->kind == kind);
    };
    const auto nested = envelope("[null,[null]]");
    auto request = limits();
    request.maximum_depth = 3;
    request.maximum_nodes = 13;
    accepted(nested, request);
    request.maximum_depth = 2;
    denied(nested, request, InputQuota::Depth);
    request.maximum_depth = 3;
    request.maximum_nodes = 12;
    denied(nested, request, InputQuota::Nodes);
    const auto members = envelope(R"({"x":{"y":null}})");
    request = limits();
    request.maximum_members = 11;
    accepted(members, request);
    request.maximum_members = 10;
    denied(members, request, InputQuota::Members);
    const auto unicode = envelope(R"("日本語日本語")");
    request = limits();
    request.maximum_decoded_string_bytes = 18;
    request.maximum_total_text_bytes = 93;
    accepted(unicode, request);
    request.maximum_decoded_string_bytes = 17;
    denied(unicode, request, InputQuota::DecodedStringBytes);
    request.maximum_decoded_string_bytes = 18;
    request.maximum_total_text_bytes = 92;
    denied(unicode, request, InputQuota::TotalTextBytes);
    const auto distinct = document(envelope(R"({"é":1,"e\u0301":2,"":"","empty":""})"));
    const auto node = model(distinct);
    const auto object = node.object();
    if (!object) {
        FAIL_CHECK("Missing normalization control object");
        return;
    }
    CHECK(object->size() == 4);
    const auto composed = object->find("é"), decomposed = object->find("e\xcc\x81"),
               empty = object->find("");
    if (!composed || !decomposed || !empty) {
        FAIL_CHECK("Decoded key was normalized or lost");
        return;
    }
    CHECK(composed->unsigned_integer() == 1);
    CHECK(decomposed->unsigned_integer() == 2);
    CHECK(empty->string() == "");
}

TEST_CASE("Parser forecast payload boundaries and unknown requests fail before parsing",
          "[verification][negative][IO4-01][N4-02][N4-08]") {
    const auto kept = document(minimal);
    const auto view = take(kept.view());
    const auto& memory = view.parser_memory;
    CHECK_FALSE(memory.proven_upper_bytes());
    CHECK(memory.sources().size() >= 9);
    const auto b = memory.breakdown();
    CHECK(memory.estimate_bytes() == b.token_vector_bytes + b.token_buffer_bytes +
                                         b.state_stack_bytes +
                                         std::max(b.syntax_error_bytes, b.number_range_bytes));
    const auto platform = memory.environment();
    CHECK(platform.parser_repair == 1);
    CHECK(platform.model_revision == 1);
    auto exact = limits();
    exact.external_memory.maximum_bytes = memory.estimate_bytes();
    core::ResourcePlan plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    CHECK(read(minimal, exact, plan, diagnostic));
    --exact.external_memory.maximum_bytes;
    core::ResourcePlan none{0};
    const auto denied = read(minimal, exact, none, diagnostic);
    REQUIRE_FALSE(denied);
    CHECK(denied.error().code == core::FailureCode::ResourceExceeded);
    if (!diagnostic || !diagnostic->parser_plan) {
        FAIL_CHECK("Missing parser admission diagnostic");
        return;
    }
    CHECK(diagnostic->reason == InputDiagnosticReason::ExternalMemoryLimit);
    CHECK(diagnostic->parser_plan->estimate_bytes() == memory.estimate_bytes());
    CHECK_FALSE(diagnostic->digest);
    CHECK(none.allocated_bytes() == 0);
    for (const auto cap :
         {UINT64_C(0), memory.estimate_bytes(), std::numeric_limits<std::uint64_t>::max()}) {
        auto strict = limits();
        strict.external_memory = {ParserAllocationMode::RequireProvenBound, cap};
        const auto result = read(minimal, strict, none, diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::Unsupported);
        CHECK(none.allocated_bytes() == 0);
    }
    auto invalid = limits();
    invalid.external_memory.mode = static_cast<ParserAllocationMode>(255);
    const auto unknown = read(minimal, invalid, none, diagnostic);
    REQUIRE_FALSE(unknown);
    CHECK(unknown.error().code == core::FailureCode::InvalidInput);
    invalid = limits();
    invalid.maximum_nodes = std::numeric_limits<std::uint64_t>::max();
    const auto overflow = read(minimal, invalid, none, diagnostic);
    REQUIRE_FALSE(overflow);
    CHECK(overflow.error().code == core::FailureCode::ResourceExceeded);
    core::ResourcePlan measured{budget};
    const auto sampled = read(minimal, limits(), measured, diagnostic);
    REQUIRE(sampled);
    core::ResourcePlan equality{measured.allocated_bytes()};
    CHECK(read(minimal, limits(), equality, diagnostic));
    core::ResourcePlan shortage{measured.allocated_bytes() - 1};
    CHECK_FALSE(read(minimal, limits(), shortage, diagnostic));
    CHECK(shortage.allocated_bytes() == 0);
}

TEST_CASE("Source hashes diagnostic copies and immutable documents survive all input owners",
          "[verification][property][negative][IO4-01][N4-03][N4-08]") {
    auto kept = [] {
        std::string input{minimal};
        core::ResourcePlan plan{budget};
        std::optional<InputDiagnostic> diagnostic;
        auto result = take(read(input, limits(), plan, diagnostic,
                                SourceId{std::numeric_limits<std::uint64_t>::max()}));
        input.assign(input.size(), 'x');
        return result;
    }();
    const auto before = take(kept.view());
    CHECK(before.identity.source.value() == std::numeric_limits<std::uint64_t>::max());
    CHECK(hex(before.identity.digest) ==
          "3db8e7560eaace99fac2b0bdcf8ab118a5275dfc4d15a4efeebb16de5e6f2934");
    CHECK(before.source_bytes.size() == minimal.size());
    CHECK(std::equal(before.source_bytes.begin(), before.source_bytes.end(),
                     std::as_bytes(std::span{minimal.data(), minimal.size()}).begin()));
    const auto spaced = document(" " + std::string{minimal});
    CHECK(take(spaced.view()).identity.source.value() == 0);
    CHECK_FALSE(take(spaced.view()).identity.digest == before.identity.digest);
    const core::FailureDetails context{core::ModelGeneration{1},
                                       core::AcceptedStateId{2},
                                       core::TrialId{3},
                                       core::NumericalGeneration{4},
                                       core::ComponentId{5},
                                       core::ProblemId{core::PanelId{6}},
                                       core::ResourceUsage{core::ResourceKind::Iterations, 7, 8},
                                       core::AcceptedStateId{9}};
    core::ResourcePlan plan{budget};
    std::optional<InputDiagnostic> diagnostic;
    const auto bad = read("[", limits(), plan, diagnostic, SourceId{7}, context);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().details == context);
    if (!diagnostic) {
        FAIL_CHECK("Missing syntax diagnostic");
        return;
    }
    CHECK(diagnostic->byte_offset == 1);
    CHECK(diagnostic->digest);
    const auto saved = *diagnostic;
    const auto duplicate = read(R"({"a":0,"a":1})", limits(), plan, diagnostic);
    REQUIRE_FALSE(duplicate);
    if (!diagnostic) {
        FAIL_CHECK("Missing duplicate diagnostic");
        return;
    }
    CHECK_FALSE(diagnostic->byte_offset);
    CHECK(diagnostic->event_ordinal);
    CHECK(saved.source.value() == 7);
    CHECK(saved.byte_offset == 1);
    const auto good = read(minimal, limits(), plan, diagnostic);
    REQUIRE(good);
    CHECK_FALSE(diagnostic);
    CHECK(take(kept.view()).identity == before.identity);
    auto moved = std::move(kept);
    // NOLINTNEXTLINE(bugprone-use-after-move): moved-from queries are explicitly supported.
    CHECK_FALSE(kept.view());
    CHECK(take(moved.view()).identity == before.identity);
    std::array<bool, 2> same{};
    std::array<std::jthread, 2> workers;
    for (std::size_t i = 0; i < 2; ++i)
        workers[i] = std::jthread([&, i] {
            core::ResourcePlan local{budget};
            std::optional<InputDiagnostic> own;
            const auto result = read(minimal, limits(), local, own, SourceId{9});
            same[i] =
                result && !own && take(result->view()).identity.digest == before.identity.digest;
        });
    for (auto& worker : workers)
        worker.join();
    CHECK(same[0]);
    CHECK(same[1]);
}

template <class T>
concept RvalueView = requires(const T&& owner) { std::move(owner).view(); };
template <class T>
concept RvalueString = requires(const T&& owner) { std::move(owner).string(); };
template <class T>
concept RvalueSources = requires(const T&& owner) { std::move(owner).sources(); };
template <class T>
concept RvalueArray = requires(const T&& owner) { std::move(owner).array(); };
template <class T>
concept RvalueObject = requires(const T&& owner) { std::move(owner).object(); };
template <class T>
concept RvalueDecimal = requires(const T&& owner) { std::move(owner).decimal_token(); };
template <class T>
concept RvalueAt = requires(const T&& owner) { std::move(owner).at(0); };
template <class T>
concept RvalueFind = requires(const T&& owner) { std::move(owner).find("x"); };
static_assert(!RvalueView<InputDocument>);
static_assert(!RvalueString<InputNodeView>);
static_assert(!RvalueSources<ParserMemoryPlan>);
static_assert(!RvalueArray<InputNodeView>);
static_assert(!RvalueObject<InputNodeView>);
static_assert(!RvalueDecimal<InputNodeView>);
static_assert(!RvalueAt<InputArrayView>);
static_assert(!RvalueAt<InputObjectView>);
static_assert(!RvalueFind<InputObjectView>);
static_assert(!std::is_default_constructible_v<InputDocument>);
static_assert(!std::is_copy_constructible_v<InputDocument>);
static_assert(!std::is_constructible_v<InputNodeOrdinal, core::Index>);
static_assert(!std::is_constructible_v<InputNodeView, const void*, core::Index>);
static_assert(!std::is_default_constructible_v<ParserMemoryPlan>);
static_assert(!std::is_convertible_v<SourceId, core::DofId>);
static_assert(!std::is_convertible_v<SourceId, core::Index>);
static_assert(std::is_nothrow_copy_constructible_v<InputDiagnostic>);
static_assert(noexcept(read_input_document({}, SourceId{0}, std::declval<const InputLimits&>(),
                                           std::declval<core::ResourcePlan&>(),
                                           std::declval<std::optional<InputDiagnostic>&>())));
