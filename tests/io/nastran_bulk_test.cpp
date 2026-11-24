#include "nastran_bulk.hpp"
#include "source_digest.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
using namespace uvlas;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using enum core::FailureCode;

constexpr NastranBulkLimits limits{10000, 1024, 128, 32, 256, 50000};
constexpr std::uint64_t memory = 50000;

template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::logic_error("Expected a value, code=" +
                               std::to_string(static_cast<unsigned>(value.error().code)));
    return std::move(*value);
}

template <class T>
T take(std::optional<T> value) {
    if (!value)
        throw std::logic_error("Expected a view");
    return *value;
}

class Reader final : public ReferenceFileReader {
public:
    explicit Reader(std::span<const std::byte> bytes) noexcept : bytes_(bytes) {}

    core::Result<ReferenceFileInfo> info() noexcept override {
        return ReferenceFileInfo{{0, {}}, bytes_.size(), {}, true};
    }

    core::Result<std::size_t> read(std::span<std::byte> out) noexcept override {
        const auto count = std::min(bytes_.size() - offset_, out.size());
        std::ranges::copy(bytes_.subspan(offset_, count), out.begin());
        offset_ += count;
        return count;
    }

private:
    std::span<const std::byte> bytes_;
    std::size_t offset_ = 0;
};

AcquiredRootSnapshot snapshot(std::string_view text,
                              ReferenceRole role = ReferenceRole::NastranDeck) {
    const auto bytes = std::as_bytes(std::span{text});
    Reader reader{bytes};
    core::ResourcePlan plan{text.size()};
    auto file =
        read_reference_file(reader, source_digest(bytes),
                            {text.size(), text.size() + 1, 3, text.size()}, plan, {}, nullptr);
    return {{{{SourceId{90}, 100, {}}, 1, 1}, SourceId{99}, role, 31}, take(std::move(file)), {}};
}

std::string fixed(std::string_view head, std::initializer_list<std::string_view> data,
                  std::size_t width, std::string_view field_ten = {}) {
    std::string result{head};
    result.resize(8, ' ');
    for (const auto field : data) {
        REQUIRE(field.size() <= width);
        result.append(field);
        result.append(width - field.size(), ' ');
    }
    result.resize(72, ' ');
    result.append(field_ten);
    return result;
}

TEST_CASE("Nastran lexer retains formats, spelling and source positions in one snapshot",
          "[IO4-03][IO4-01][verification]") {
    const auto small_head =
        fixed("CORD2R", {"42", "", "0.0", "-0.0", "0.0", "1.0", "0.0", "0.0"}, 8, "xA1");
    const auto small_named = fixed("+A1", {"1.0", "1.0", "0.0"}, 8);
    const auto large_head =
        fixed("GRID*", {"12345678901", "", "1234567890.12D-2", "-1.25+3"}, 16, "+L7");
    const auto large_named = fixed("*L7", {"0.0", "", "136", ""}, 16, "note");
    const auto automatic = fixed("", {"91", "92"}, 8);
    const auto text = "$ project note\r\nBEGIN BULK\r\n" + small_named + "\r\n" + large_head +
                      "\r\n$ another comment\r\n" + small_head + "\r\n" + large_named + "\r\n" +
                      fixed("SPC1", {"3", "12", "11", "12", "13", "14", "15", "16"}, 8) + "\r\n" +
                      automatic +
                      "\r\nRBE2,1,2,123,3,4,5,6,7,8,9\r\nSPC1,2,123,3,4,5,6,7,8\r\n,10,,+11 $ "
                      "unchanged\r\nENDDATA";
    const char* original_address = nullptr;
    SourceDigest digest{};
    const auto deck = [&] {
        auto original = snapshot(text);
        const auto input = take(original.file.view());
        original_address = reinterpret_cast<const char*>(input.bytes.data());
        digest = input.digest;
        core::ResourcePlan plan{memory};
        auto prepared = take(NastranBulkDeck::prepare(std::move(original), limits, plan, {}));
        // NOLINTNEXTLINE(bugprone-use-after-move): success moves this snapshot into the deck.
        CHECK_FALSE(original.file.view());
        return prepared;
    }();
    const auto result = take(deck.view());
    CHECK(result.bytes.data() == original_address);
    CHECK(result.bytes == text);
    CHECK(result.digest == digest);
    CHECK(result.location.parent.source.source == SourceId{90});
    CHECK(result.location.source == SourceId{99});
    CHECK(result.location.declaration_node == 31);
    CHECK(result.work > text.size());
    REQUIRE(result.cards.size() == 5);

    const auto& large = result.cards[0];
    CHECK(large.format == NastranBulkFormat::Large);
    CHECK(result.text(large.name) == "GRID");
    const auto large_fields = result.data_fields(large);
    REQUIRE(large_fields.size() == 8);
    CHECK(result.text(large_fields[0]) == "12345678901");
    CHECK(result.text(large_fields[1]).empty());
    CHECK(result.text(large_fields[2]) == "1234567890.12D-2");
    CHECK(result.text(large_fields[3]) == "-1.25+3");
    CHECK(take(result.integer(large_fields[0])) == 12345678901);
    const auto explicit_exponent = take(result.real(large_fields[2]));
    const auto implicit_exponent = take(result.real(large_fields[3]));
    CHECK(explicit_exponent.decimal() == "1234567890.12e-2");
    CHECK(implicit_exponent.decimal() == "-1.25e+3");
    CHECK_FALSE(result.real(large_fields[0]));
    CHECK_FALSE(result.integer(large_fields[2]));
    CHECK(large_fields[2].position == NastranBulkPosition{text.find("1234567890.12D-2"), 4, 41});
    CHECK(large_fields[4].position == NastranBulkPosition{text.find("*L7") + 8, 7, 9});
    CHECK(result.text(large_fields[6]) == "136");

    const auto& small = result.cards[1];
    CHECK(small.format == NastranBulkFormat::Small);
    CHECK(result.text(small.name) == "CORD2R");
    const auto small_fields = result.data_fields(small);
    REQUIRE(small_fields.size() == 16);
    CHECK(result.text(small_fields[0]) == "42");
    CHECK(result.text(small_fields[3]) == "-0.0");
    CHECK(result.text(small_fields[8]) == "1.0");
    CHECK(small_fields[8].position == NastranBulkPosition{text.find("+A1") + 8, 3, 9});
    CHECK(small_fields[15].position == NastranBulkPosition{text.find("+A1") + 64, 3, 65});

    const auto auto_fields = result.data_fields(result.cards[2]);
    REQUIRE(auto_fields.size() == 16);
    CHECK(result.text(auto_fields[7]) == "16");
    CHECK(result.text(auto_fields[8]) == "91");
    CHECK(auto_fields[8].position.line == 9);

    const auto& free = result.cards[3];
    CHECK(free.format == NastranBulkFormat::Free);
    CHECK(result.text(free.name) == "RBE2");
    const auto free_fields = result.data_fields(free);
    REQUIRE(free_fields.size() == 10);
    CHECK(result.text(free_fields[9]) == "9");

    const auto& free_continued = result.cards[4];
    CHECK(free_continued.format == NastranBulkFormat::Free);
    CHECK(result.text(free_continued.name) == "SPC1");
    const auto continued = result.data_fields(free_continued);
    REQUIRE(continued.size() == 11);
    CHECK(result.text(continued[7]) == "8");
    CHECK(result.text(continued[8]) == "10");
    CHECK(result.text(continued[9]).empty());
    CHECK(continued[9].position == NastranBulkPosition{text.find(",,+11") + 1, 12, 5});
    CHECK(result.text(continued[10]) == "+11");
    CHECK(take(result.integer(continued[10])) == 11);
    CHECK_FALSE(result.integer(continued[9]));
}

TEST_CASE("The checked lexical limits count padding and preserve the snapshot on failure",
          "[IO4-03][N4-02][negative]") {
    const auto text = std::string{"GRID    42\nSPC1,1, 12,,45,"};
    auto initial = snapshot(text);
    core::ResourcePlan plan{memory};
    const auto deck = take(NastranBulkDeck::prepare(std::move(initial), limits, plan, {}));
    const auto view = take(deck.view());
    REQUIRE(view.cards.size() == 2);
    const auto first = view.data_fields(view.cards[0]);
    const auto second = view.data_fields(view.cards[1]);
    REQUIRE(first.size() == 8);
    REQUIRE(second.size() == 5);
    CHECK(first[1].position == NastranBulkPosition{10, 1, 17});
    CHECK(first[7].position == NastranBulkPosition{10, 1, 65});
    CHECK(view.text(first[1]).empty());
    CHECK(view.text(second[1]) == "12");
    CHECK(view.text(second[2]).empty());
    CHECK(view.text(second[4]).empty());

    auto caps = limits;
    caps.file_bytes = text.size();
    caps.physical_lines = 2;
    caps.cards = view.cards.size();
    caps.fields = view.fields.size();
    caps.work = view.work;
    core::ResourcePlan exact{plan.allocated_bytes()};
    auto repeated = snapshot(text);
    const auto replay = NastranBulkDeck::prepare(std::move(repeated), caps, exact, {});
    CHECK(replay);
    CHECK(exact.allocated_bytes() == plan.allocated_bytes());

    auto amount = exact.allocated_bytes();
    SECTION("file size") {
        --caps.file_bytes;
    }
    SECTION("physical lines") {
        --caps.physical_lines;
    }
    SECTION("line bytes") {
        caps.physical_line_bytes = 9;
    }
    SECTION("card count") {
        --caps.cards;
    }
    SECTION("field count") {
        --caps.fields;
    }
    SECTION("work") {
        --caps.work;
    }
    SECTION("storage") {
        --amount;
    }
    auto input = snapshot(text);
    core::ResourcePlan restricted{amount};
    std::optional<NastranBulkDiagnostic> diagnostic;
    const auto rejected =
        NastranBulkDeck::prepare(std::move(input), caps, restricted, {}, &diagnostic);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == ResourceExceeded);
    // NOLINTNEXTLINE(bugprone-use-after-move): rejected calls leave the input snapshot intact.
    CHECK(take(input.file.view()).bytes.size() == text.size());
    REQUIRE(diagnostic);
    if (diagnostic)
        CHECK(diagnostic->reason == NastranBulkFailure::Limit);
}

TEST_CASE("Ambiguous forms and inconsistent continuation graphs fail visibly",
          "[IO4-03][IO4-01][negative]") {
    std::string text;
    auto code = InvalidInput;
    auto line = std::size_t{1};
    SECTION("orphan automatic small field") {
        text = "        12\n";
    }
    SECTION("orphan named small field") {
        text = "+LOST   12\n";
    }
    SECTION("large format has two named halves at minimum") {
        text = "GRID*   12\n";
    }
    SECTION("explicit tail marker cannot disappear") {
        text = fixed("GRID", {"12"}, 8, "+LOST");
        code = Unsupported;
    }
    SECTION("one outgoing label cannot be owned twice") {
        text = fixed("GRID", {"10"}, 8, "xYES") + "\n" + fixed("SPC1", {"11"}, 8, "yYES") + "\n" +
               fixed("+YES", {"2"}, 8);
        line = 2;
    }
    SECTION("cycle without a card head") {
        text = fixed("+TWO", {"10"}, 8, "+ONE") + "\n" + fixed("+ONE", {"11"}, 8, "+TWO");
    }
    SECTION("fixed mixed-width handoff") {
        text = fixed("GRID*", {"11"}, 16, "*NEXT") + "\n" + fixed("+NEXT", {"2"}, 8);
        line = 2;
        code = Unsupported;
    }
    SECTION("free continuation after an incomplete group") {
        text = "CORD2R,1,0,0.0\n,1,2,3\n";
        line = 2;
        code = Unsupported;
    }
    SECTION("named free continuation") {
        text = "CORD2R,1,0,0.0\n+N,1,2,3\n";
        line = 2;
        code = Unsupported;
    }
    SECTION("patched fields cannot contain embedded spaces") {
        text = fixed("GRID", {"11", "1 .0"}, 8);
    }
    SECTION("an INCLUDE needs its own adopted extraction") {
        text = "INCLUDE 'part.dat'";
        code = Unsupported;
    }
    SECTION("outside Bulk Data") {
        text = "SOL 101\nCEND\nBEGIN BULK\nGRID,1\nENDDATA";
        code = Unsupported;
    }
    SECTION("data after ENDDATA") {
        text = "GRID,1\nENDDATA\nGRID,2\n";
        line = 3;
    }
    SECTION("NUL in source") {
        text = "GRID,";
        text.push_back('\0');
    }
    auto input = snapshot(text);
    core::ResourcePlan plan{memory};
    std::optional<NastranBulkDiagnostic> diagnostic;
    const auto rejected = NastranBulkDeck::prepare(std::move(input), limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == code);
    CHECK(rejected.error().operation == core::Operation::PrepareModel);
    REQUIRE(diagnostic);
    if (diagnostic) {
        CHECK(diagnostic->position.line == line);
        CHECK(diagnostic->position.byte < text.size());
    }
    // NOLINTNEXTLINE(bugprone-use-after-move): rejected calls leave the input snapshot intact.
    CHECK(input.file.view());
}

TEST_CASE("Nastran lexer cannot consume other role snapshots and observes stop",
          "[IO4-03][IO4-01][negative]") {
    auto other = snapshot("GRID,1\n", ReferenceRole::NastranSidecar);
    core::ResourcePlan plan{0};
    std::optional<NastranBulkDiagnostic> diagnostic;
    const auto wrong_role =
        NastranBulkDeck::prepare(std::move(other), limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(wrong_role);
    CHECK(wrong_role.error().code == InvalidInput);
    REQUIRE(diagnostic);
    std::stop_source stop; // NOLINT(misc-const-correctness): MSVC request_stop is non-const.
    REQUIRE(stop.request_stop());
    auto input = snapshot("GRID,1\n");
    const auto halted =
        NastranBulkDeck::prepare(std::move(input), limits, plan, stop.get_token(), &diagnostic);
    REQUIRE_FALSE(halted);
    CHECK(halted.error().code == Cancelled);
    CHECK_FALSE(diagnostic);
    CHECK(plan.allocated_bytes() == 0);
}
} // namespace
