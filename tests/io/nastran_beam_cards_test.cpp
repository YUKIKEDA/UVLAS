#include "nastran_beam_cards.hpp"
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
using Reason = NastranBeamFailure;

constexpr std::uint64_t memory = 200000;
constexpr NastranBulkLimits lexical{10000, 1024, 128, 32, 256, 50000};
constexpr NastranCoordinateLimits limits{8, 8, 8, 10000};

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

NastranBulkDeck lexical_deck(std::string_view text) {
    const auto bytes = std::as_bytes(std::span{text});
    Reader reader{bytes};
    core::ResourcePlan plan{memory};
    auto file = take(read_reference_file(reader, source_digest(bytes),
                                         {text.size(), text.size() + 1, 3, text.size()}, plan, {},
                                         nullptr));
    AcquiredRootSnapshot root{
        {{{SourceId{90}, 100, {}}, 1, 1}, SourceId{99}, ReferenceRole::NastranDeck, 31},
        std::move(file),
        {}};
    return take(NastranBulkDeck::prepare(std::move(root), lexical, plan, {}));
}

constexpr NastranBeamLimits beam_limits{16, 16, 20000};
constexpr std::string_view points = "GRID,11\nGRID,12\nGRID,13\nGRID,42\n";
constexpr std::string_view property_card = "PBEAM,100\n";

NastranShellCards shell_cards(std::string_view source) {
    auto deck = lexical_deck(source);
    core::ResourcePlan plan{memory};
    auto coordinates = take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
    return take(NastranShellCards::prepare(std::move(coordinates), {16, 16, 20000}, plan, {}));
}

std::string fixed(std::string_view name, std::initializer_list<std::string_view> values,
                  std::size_t width = 8, std::string_view label = {}) {
    std::string line{name};
    line.resize(8, ' ');
    for (const auto field : values) {
        REQUIRE(field.size() <= width);
        line.append(field);
        line.append(width - field.size(), ' ');
    }
    line.resize(72, ' ');
    return line + std::string{label} + "\n";
}

TEST_CASE("CBEAM preserves ordered ends, exact offsets, releases and PBEAM references",
          "[IO4-03][IO4-01][verification]") {
    const auto source = std::string{points} + "CBEAM,100,,11,12,,2.0D-3,-.25,BOG\n"
                                              ",513,46,-.10,0.0,,1.25D+1,,3.0\n"
                                              ",42,0\n"
                                              "CBEAM,101,90,12,11,13,,,GGG\n"
                                              "CBEAM,102,90,12,11,13\n"
                                              "PSHELL,90\nPBEAM,90\nPBEAM,100\n";
    const auto parent = shell_cards(source);
    const auto input = take(parent.view());
    core::ResourcePlan plan{memory};
    auto prepared = take(NastranBeamCards::prepare(parent, beam_limits, plan, {}));
    const auto view = take(prepared.view());
    REQUIRE(view.beams.size() == 3);
    REQUIRE(view.properties.size() == 2);
    CHECK(view.source.digest == source_digest(std::as_bytes(std::span{source})));
    CHECK(view.source.location.source == input.input.bulk.location.source);
    const auto& bulk = input.input.bulk;
    const auto& first = view.beams[0];
    CHECK(first.id == NastranElementId{100});
    CHECK(bulk.text(bulk.cards[first.card_ordinal].name) == "CBEAM");
    CHECK(first.property.id == NastranPropertyId{100});
    CHECK(take(first.property.field).size == 0);
    CHECK(view.properties[first.property.declaration_ordinal].id == first.property.id);
    const auto& property = view.properties[first.property.declaration_ordinal];
    CHECK(bulk.text(bulk.cards[property.card_ordinal].name) == "PBEAM");
    CHECK(first.a.point.id == NastranGridId{11});
    CHECK(first.b.point.id == NastranGridId{12});
    CHECK(input.input.grids[first.a.point.grid_ordinal].id == first.a.point.id);
    CHECK(input.input.grids[first.b.point.grid_ordinal].id == first.b.point.id);
    const auto& vector = std::get<NastranBeamVector>(first.orientation);
    CHECK_FALSE(vector[0]);
    const auto vy = take(vector[1]);
    CHECK(vy.value.decimal() == "2.0e-3");
    CHECK(bulk.text(vy.field) == "2.0D-3");
    const auto vz = take(vector[2]);
    CHECK(vz.value.decimal() == "-0.25");
    CHECK(first.interpretation.orientation == NastranBeamOrientationBasis::Basic);
    CHECK(first.interpretation.a == NastranBeamOffsetBasis::OffsetSystem);
    CHECK(first.interpretation.b == NastranBeamOffsetBasis::GridDisplacement);
    CHECK(bulk.text(take(first.interpretation.field)) == "BOG");
    CHECK(first.a.releases.element_components == 0b010101);
    CHECK(first.b.releases.element_components == 0b101000);
    CHECK(bulk.text(take(first.a.releases.field)) == "513");
    const auto a_offset = take(first.a.offset[0]);
    CHECK(a_offset.value.decimal() == "-0.10");
    CHECK(a_offset.field.position == NastranBulkPosition{source.find("-.10"), 6, 9});
    const auto a_zero = take(first.a.offset[1]);
    CHECK(a_zero.value.decimal() == "0.0");
    CHECK_FALSE(first.a.offset[2]);
    const auto b_offset = take(first.b.offset[0]);
    CHECK(b_offset.value.decimal() == "1.25e+1");
    CHECK_FALSE(first.b.offset[1]);
    const auto sa = take(first.a.warping);
    CHECK(sa.id == NastranGridId{42});
    CHECK(input.input.grids[take(sa.grid_ordinal)].id == sa.id);
    const auto sb = take(first.b.warping);
    CHECK(sb.id == NastranGridId{0});
    CHECK_FALSE(sb.grid_ordinal);
    CHECK(bulk.text(sb.field) == "0");
    const auto& second = view.beams[1];
    CHECK(second.property.id == NastranPropertyId{90});
    const auto& g0 = std::get<NastranBeamGridReference>(second.orientation);
    CHECK(g0.id == NastranGridId{13});
    CHECK(input.input.grids[g0.grid_ordinal].id == g0.id);
    CHECK(bulk.text(take(second.interpretation.field)) == "GGG");
    CHECK(second.interpretation.orientation == NastranBeamOrientationBasis::GridDisplacement);
    const auto& third = view.beams[2];
    CHECK_FALSE(third.interpretation.field);
    CHECK(third.interpretation.orientation == NastranBeamOrientationBasis::GridDisplacement);
    CHECK_FALSE(third.a.warping);
    CHECK(take(parent.view()).input.bulk.bytes.data() == bulk.bytes.data());
    CHECK(take(parent.view()).input.bulk.bytes == source);

    const auto moved = std::move(prepared);
    // NOLINTNEXTLINE(bugprone-use-after-move): moving clears the previous view.
    CHECK_FALSE(prepared.view());
    CHECK(take(moved.view()).beams.size() == 3);
}

TEST_CASE("CBEAM large continuations preserve an empty offsets row and both warping columns",
          "[IO4-03][verification]") {
    const auto source = std::string{points} + std::string{property_card} +
                        fixed("CBEAM*", {"1", "100", "11", "12"}, 16, "*B1") +
                        fixed("*B1", {"13", "", "", "BGO"}, 16, "*B2") +
                        fixed("*B2", {}, 16, "*B3") + fixed("*B3", {}, 16, "*B4") +
                        fixed("*B4", {"0", "42"}, 16, "*B5") + fixed("*B5", {}, 16);
    core::ResourcePlan plan{memory};
    const auto prepared = [&] {
        const auto parent = shell_cards(source);
        return take(NastranBeamCards::prepare(parent, beam_limits, plan, {}));
    }();
    const auto view = take(prepared.view());
    REQUIRE(view.beams.size() == 1);
    const auto& row = view.beams[0];
    const auto& g0 = std::get<NastranBeamGridReference>(row.orientation);
    CHECK(g0.id == NastranGridId{13});
    CHECK(row.interpretation.b == NastranBeamOffsetBasis::OffsetSystem);
    CHECK(row.a.releases.element_components == 0);
    CHECK_FALSE(row.a.offset[0]);
    CHECK_FALSE(row.b.offset[2]);
    const auto sa = take(row.a.warping);
    const auto sb = take(row.b.warping);
    CHECK(sa.id == NastranGridId{0});
    CHECK(sb.id == NastranGridId{42});
    CHECK(sa.field.position.line == 10);
    CHECK(sb.field.position.column == 25);
}

TEST_CASE("CBEAM small fields preserve all eight MSC offset interpretation combinations",
          "[IO4-03][verification]") {
    using V = NastranBeamOrientationBasis;
    using O = NastranBeamOffsetBasis;
    struct Offt {
        std::string_view text;
        V vector;
        O a, b;
    };
    constexpr auto vg = V::GridDisplacement;
    constexpr auto vb = V::Basic;
    constexpr auto g = O::GridDisplacement;
    constexpr auto o = O::OffsetSystem;
    constexpr std::array variants{
        Offt{"GGG", vg, g, g}, Offt{"BGG", vb, g, g}, Offt{"GGO", vg, g, o}, Offt{"BGO", vb, g, o},
        Offt{"GOG", vg, o, g}, Offt{"BOG", vb, o, g}, Offt{"GOO", vg, o, o}, Offt{"BOO", vb, o, o}};
    for (const auto& expected : variants) {
        CAPTURE(expected.text);
        const auto source = std::string{points} + std::string{property_card} +
                            fixed("CBEAM", {"1", "100", "11", "12", "1.0", "", "", expected.text});
        const auto parent = shell_cards(source);
        core::ResourcePlan plan{memory};
        const auto result = take(NastranBeamCards::prepare(parent, beam_limits, plan, {}));
        const auto view = take(result.view());
        REQUIRE(view.beams.size() == 1);
        const auto& mode = view.beams[0].interpretation;
        CHECK(mode.orientation == expected.vector);
        CHECK(mode.a == expected.a);
        CHECK(mode.b == expected.b);
    }
}

TEST_CASE("CBEAM rejects invalid field combinations at their source location",
          "[IO4-03][negative][verification]") {
    struct FailureCase {
        std::string_view body;
        Reason reason;
        core::FailureCode code;
        std::size_t line;
        std::size_t column;
    };
    constexpr std::array cases{
        FailureCase{"CBEAM,1,100,11,99,13\n", Reason::UndefinedGrid, InvalidInput, 6, 16},
        FailureCase{"CBEAM,1,100,11,11,13\n", Reason::DuplicateEndpoint, InvalidInput, 6, 16},
        FailureCase{"CBEAM,1,100,11,12,99\n", Reason::UndefinedGrid, InvalidInput, 6, 19},
        FailureCase{"CBEAM,1,100,11,12,11\n", Reason::OrientationEndpoint, InvalidInput, 6, 19},
        FailureCase{"CBEAM,1,100,11,12,12\n", Reason::OrientationEndpoint, InvalidInput, 6, 19},
        FailureCase{"CBEAM,1,100,11,12,13,1.0\n", Reason::ExtraField, InvalidInput, 6, 22},
        FailureCase{"CBEAM,1,100,11,12,13,,1.0\n", Reason::ExtraField, InvalidInput, 6, 23},
        FailureCase{"CBEAM,1,100,11,12\n", Reason::MissingField, InvalidInput, 6, 1},
        FailureCase{"CBEAM,1,100,11,12,13,,,OOO\n", Reason::UnsupportedField, Unsupported, 6, 24},
        FailureCase{"CBEAM,1,100,11,12,13,,,GEG\n", Reason::UnsupportedField, Unsupported, 6, 24},
        FailureCase{"CBEAM,1,100,11,12,13,,,0.5\n", Reason::UnsupportedField, Unsupported, 6, 24},
        FailureCase{"CBEAM,1,90,11,12,13\n", Reason::UndefinedProperty, InvalidInput, 6, 9},
        FailureCase{"CBEAM,1,90,11,12,13\nPSHELL,90\n", Reason::WrongPropertyKind, InvalidInput, 6,
                    9},
        FailureCase{"PBEAM,100\n", Reason::DuplicateProperty, InvalidInput, 6, 7},
        FailureCase{"CBEAM,1,000000100,11,12,13\n", Reason::InvalidField, InvalidInput, 6, 9},
        FailureCase{"CBEAM,1,100,11,12,13,,,GGG\n,111\n", Reason::InvalidField, InvalidInput, 7, 2},
        FailureCase{"CBEAM,1,100,11,12,13,,,GGG\n,123456\n", Reason::InvalidField, InvalidInput, 7,
                    2},
        FailureCase{"CBEAM,1,100,11,12,13,,,GGG\n,,,,,,,,\n,99\n", Reason::UndefinedWarpingPoint,
                    InvalidInput, 8, 2},
        FailureCase{"CBEAM,1,100,11,12,13,,,GGG\n,,,,,,,,\n,0,-1\n", Reason::InvalidField,
                    InvalidInput, 8, 4},
        FailureCase{"CBEAM,1,100,11,12,13,,,GGG\n,,,,,,,,\n,0,0,2\n", Reason::ExtraField,
                    InvalidInput, 8, 6}};
    for (const auto& entry : cases) {
        CAPTURE(entry.body);
        const auto source =
            std::string{points} + std::string{property_card} + std::string{entry.body};
        const auto parent = shell_cards(source);
        core::ResourcePlan plan{memory};
        std::optional<NastranBeamDiagnostic> diagnostic;
        const auto result = NastranBeamCards::prepare(parent, beam_limits, plan, {}, &diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == entry.code);
        const auto actual = take(diagnostic);
        CHECK(actual.reason == entry.reason);
        CHECK(actual.position.line == entry.line);
        CHECK(actual.position.column == entry.column);
        CHECK(take(parent.view()).input.bulk.bytes == source);
    }
}

TEST_CASE("CBEAM budgets bound declarations and indexes while leaving the deck usable",
          "[IO4-03][N4-02][negative][verification]") {
    const auto source = std::string{points} + std::string{property_card} + "CBEAM,100,,11,12,13\n";
    const auto parent = shell_cards(source);
    core::ResourcePlan enough{memory};
    const auto no_beams = NastranBeamCards::prepare(parent, {0, 1, beam_limits.work}, enough, {});
    REQUIRE_FALSE(no_beams);
    CHECK(no_beams.error().code == ResourceExceeded);
    const auto no_properties =
        NastranBeamCards::prepare(parent, {1, 0, beam_limits.work}, enough, {});
    REQUIRE_FALSE(no_properties);
    CHECK(no_properties.error().code == ResourceExceeded);
    CHECK(enough.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness) MSVC request_stop() is non-const.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    const auto cancelled =
        NastranBeamCards::prepare(parent, beam_limits, enough, cancel.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK(enough.allocated_bytes() == 0);
    const auto measured = NastranBeamCards::prepare(parent, beam_limits, enough, {});
    REQUIRE(measured);
    const auto work = take(measured->view()).work;
    const auto bytes = enough.allocated_bytes();
    REQUIRE(work > 0);
    REQUIRE(bytes > 0);
    core::ResourcePlan exact{bytes};
    CHECK(NastranBeamCards::prepare(parent, {1, 1, work}, exact, {}));
    CHECK(exact.remaining_bytes() == 0);
    core::ResourcePlan short_memory{bytes - 1};
    const auto byte_rejected = NastranBeamCards::prepare(parent, {1, 1, work}, short_memory, {});
    REQUIRE_FALSE(byte_rejected);
    CHECK(byte_rejected.error().code == ResourceExceeded);
    CHECK(short_memory.allocated_bytes() == 0);
    core::ResourcePlan full_memory{bytes};
    const auto work_rejected = NastranBeamCards::prepare(parent, {1, 1, work - 1}, full_memory, {});
    REQUIRE_FALSE(work_rejected);
    CHECK(work_rejected.error().code == ResourceExceeded);
    CHECK(take(parent.view()).input.bulk.bytes == source);
}

} // namespace
