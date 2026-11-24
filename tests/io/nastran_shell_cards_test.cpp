#include "nastran_shell_cards.hpp"
#include "nastran_support_cards.hpp"
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
using Reason = NastranShellFailure;

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

constexpr NastranShellLimits shell_limits{16, 16, 5000};
constexpr std::string_view grids = "GRID,17\nGRID,8\nGRID,41\nGRID,22\n";
constexpr std::string_view triangle = "CTRIA3,10,20,8,41,17\n";

NastranCoordinateCards coordinates(std::string_view source) {
    auto deck = lexical_deck(source);
    core::ResourcePlan plan{memory};
    return take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
}

std::string fixed(std::string_view name, std::initializer_list<std::string_view> values,
                  std::size_t width = 8, std::string_view tail = {}) {
    std::string s{name};
    s.resize(8, ' ');
    for (const auto value : values) {
        REQUIRE(value.size() <= width);
        s.append(value);
        s.append(width - value.size(), ' ');
    }
    s.resize(72, ' ');
    return s + std::string{tail} + "\n";
}

TEST_CASE("Nastran shells preserve ordered corner identities, PSHELL and PCOMP across field forms",
          "[IO4-03][IO4-01][verification]") {
    const auto source =
        fixed("ctria3", {"10", "", "8", "41", "17", "-1.5D+1", "-.125"}, 8, "+TRI") +
        fixed("+TRI", {"", "0", "0.0", "1.25-3", ""}) +
        fixed("CQUAD4*", {"11", "40", "41", "22"}, 16, "*Q1") +
        fixed("*Q1", {"17", "8", "9", ""}, 16, "*Q2") + fixed("*Q2", {"", "", "", ""}, 16, "*Q3") +
        fixed("*Q3", {}, 16) +
        "GRID,17\nGRID,8\nGRID,41\n"
        "GRID*,22,,1.25D+1,-2.0,\n*,3.0\n"
        "PSHELL,10,31,0.2\nPCOMP,40\n"
        "CORD2R,9,,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n";
    core::ResourcePlan plan{memory};
    const char* original_bytes = nullptr;
    const auto shells = [&] {
        auto input = coordinates(source);
        original_bytes = take(input.view()).bulk.bytes.data();
        auto parsed = take(NastranShellCards::prepare(std::move(input), shell_limits, plan, {}));
        // NOLINTNEXTLINE(bugprone-use-after-move): success leaves the input view absent.
        CHECK_FALSE(input.view());
        return parsed;
    }();
    const auto view = take(shells.view());
    CHECK(view.input.bulk.bytes == source);
    CHECK(view.input.bulk.bytes.data() == original_bytes);
    CHECK(view.input.bulk.digest == source_digest(std::as_bytes(std::span{source})));
    REQUIRE(view.elements.size() == 2);
    REQUIRE(view.properties.size() == 2);
    const auto& tri = view.elements[0];
    CHECK(tri.id == NastranElementId{10});
    CHECK(tri.property.id == NastranPropertyId{10});
    CHECK(take(tri.property.field).size == 0);
    CHECK(view.properties[tri.property.property_ordinal].kind == NastranShellPropertyKind::Pshell);
    CHECK(view.input.bulk.text(view.input.bulk.cards[tri.card_ordinal].name) == "ctria3");
    const auto& tri_corners = std::get<std::array<NastranShellCorner, 3>>(tri.corners);
    CHECK(tri_corners[0].id == NastranGridId{8});
    CHECK(tri_corners[1].id == NastranGridId{41});
    CHECK(tri_corners[2].id == NastranGridId{17});
    for (const auto& corner : tri_corners)
        CHECK(view.input.grids[corner.grid_ordinal].id == corner.id);
    const auto angle = std::get<NastranShellReal>(tri.orientation);
    CHECK(angle.value.decimal() == "-1.5e+1");
    CHECK(view.input.bulk.text(angle.field) == "-1.5D+1");
    CHECK(angle.field.position == NastranBulkPosition{source.find("-1.5D+1"), 1, 49});
    const auto offset = take(tri.offset);
    CHECK(offset.value.decimal() == "-0.125");
    CHECK(tri.thickness_mode == NastranShellThicknessMode::Absolute);
    const auto flag = take(tri.thickness_mode_field);
    CHECK(flag.position.line == 2);
    CHECK(flag.position.column == 17);
    const auto zero = take(tri_corners[0].thickness);
    CHECK(zero.value.decimal() == "0.0");
    const auto nonzero = take(tri_corners[1].thickness);
    CHECK(nonzero.value.decimal() == "1.25e-3");
    CHECK_FALSE(tri_corners[2].thickness);
    const auto& quad = view.elements[1];
    CHECK(quad.id == NastranElementId{11});
    CHECK(quad.property.id == NastranPropertyId{40});
    CHECK(view.properties[quad.property.property_ordinal].kind == NastranShellPropertyKind::Pcomp);
    const auto& quad_corners = std::get<std::array<NastranShellCorner, 4>>(quad.corners);
    CHECK(quad_corners[0].id == NastranGridId{41});
    CHECK(quad_corners[1].id == NastranGridId{22});
    REQUIRE(quad_corners[1].grid_ordinal < view.input.grids.size());
    const auto& large_point = view.input.grids[quad_corners[1].grid_ordinal];
    CHECK(large_point.id == NastranGridId{22});
    CHECK(large_point.position[0].decimal() == "1.25e+1");
    CHECK(large_point.position[2].decimal() == "3.0");
    CHECK(quad_corners[2].id == NastranGridId{17});
    CHECK(quad_corners[3].id == NastranGridId{8});
    const auto axis = std::get<NastranShellMaterialAxis>(quad.orientation);
    CHECK(axis.id == NastranCoordinateId{9});
    CHECK(view.input.coordinates[take(axis.coordinate_ordinal)].id == axis.id);
    const auto free_values = view.input.bulk.data_fields(
        view.input.bulk.cards[view.properties[quad.property.property_ordinal].card_ordinal]);
    REQUIRE_FALSE(free_values.empty());
    CHECK(view.input.bulk.text(free_values[0]) == "40");
}

TEST_CASE("Nastran free shells distinguish implicit and explicit axes and fractional thickness",
          "[IO4-03][verification]") {
    const auto source = std::string{grids} + "PSHELL,20,30,0.1\n"
                                             "CTRIA3,10,20,8,41,17\n"
                                             "CTRIA3,11,20,8,41,17,0\n"
                                             "CTRIA3,12,20,8,41,17,0.0\n"
                                             "CQUAD4,13,20,8,41,17,22,45.0,0.0\n"
                                             ",,1,1.25,,0.5,1.0\n";
    auto input = coordinates(source);
    core::ResourcePlan plan{memory};
    const auto output = take(NastranShellCards::prepare(std::move(input), shell_limits, plan, {}));
    const auto view = take(output.view());
    REQUIRE(view.elements.size() == 4);
    CHECK(std::holds_alternative<std::monostate>(view.elements[0].orientation));
    const auto axis = std::get<NastranShellMaterialAxis>(view.elements[1].orientation);
    CHECK(axis.id == NastranCoordinateId{0});
    CHECK_FALSE(axis.coordinate_ordinal);
    const auto angle = std::get<NastranShellReal>(view.elements[2].orientation);
    CHECK(angle.value.decimal() == "0.0");
    const auto& row = view.elements[3];
    CHECK(row.thickness_mode == NastranShellThicknessMode::Fraction);
    const auto offset = take(row.offset);
    CHECK(offset.value.decimal() == "0.0");
    const auto& q = std::get<std::array<NastranShellCorner, 4>>(row.corners);
    const auto first = take(q[0].thickness);
    const auto last = take(q[3].thickness);
    CHECK(first.value.decimal() == "1.25");
    CHECK_FALSE(q[1].thickness);
    CHECK(last.value.decimal() == "1.0");
    CHECK(last.field.position.line == 10);
    CHECK(view.input.bulk.text(last.field) == "1.0");
}

TEST_CASE("Nastran shell failures preserve the coordinate owner and report the source field",
          "[IO4-03][N4-02][negative][verification]") {
    struct FailureCase {
        std::string_view body;
        Reason reason;
        std::size_t relative_line;
        std::size_t column;
    };
    constexpr std::array cases{
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,8\n", Reason::DuplicateGrid, 2, 19},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,80\n", Reason::UndefinedGrid, 2, 19},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,\n", Reason::MissingField, 2, 19},
        FailureCase{"PSHELL,21\nCTRIA3,10,20,8,41,17\n", Reason::UndefinedProperty, 2, 11},
        FailureCase{"PBEAM,20\nCTRIA3,10,20,8,41,17\n", Reason::WrongPropertyKind, 2, 11},
        FailureCase{"PSHELL,20\nPCOMP,20\n", Reason::DuplicateProperty, 2, 7},
        FailureCase{"CTRIA3,10,20,8,41,17\nCBEAM,10\n", Reason::DuplicateElement, 2, 7},
        FailureCase{"CTRIA3,10,20,8,41,17\nCTRIA3,10,20,8,41,17\n", Reason::DuplicateElement, 2, 8},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,99\n", Reason::UndefinedCoordinate, 2, 22},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,-1\n", Reason::InvalidField, 2, 22},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,TOP\n", Reason::InvalidField, 2, 23},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,,1\n", Reason::ExtraField, 2, 24},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,,\n,1\n", Reason::ExtraField, 3, 2},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,,\n,,2\n", Reason::InvalidField, 3, 3},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,,\n,,,0.0,0.0,-0.0\n", Reason::InvalidField,
                    3, 12},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,,\n,,,,-.1\n", Reason::InvalidField, 3, 5},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,,\n,,1,,,0.0\n", Reason::InvalidField, 3, 7},
        FailureCase{"PSHELL,20\nCTRIA3,10,20,8,41,17,,,\n,,,,,,9\n", Reason::ExtraField, 3, 7},
        FailureCase{"PCOMP,10000000\n", Reason::InvalidField, 1, 7},
        FailureCase{"CTRIA3,100000000,20,8,41,17\n", Reason::InvalidField, 1, 8}};
    for (const auto& test : cases) {
        CAPTURE(test.body);
        const auto source = std::string{grids} + std::string{test.body};
        auto owner = coordinates(source);
        core::ResourcePlan plan{memory};
        std::optional<NastranShellDiagnostic> diagnostic;
        const auto result =
            NastranShellCards::prepare(std::move(owner), shell_limits, plan, {}, &diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == InvalidInput);
        const auto failure = take(diagnostic);
        CHECK(failure.reason == test.reason);
        CHECK(failure.position.line == 4 + test.relative_line);
        CHECK(failure.position.column == test.column);
        // NOLINTNEXTLINE(bugprone-use-after-move): rejection leaves the original owner intact.
        CHECK(take(owner.view()).bulk.bytes == source);
        CHECK(plan.remaining_bytes() < memory);
    }
}

TEST_CASE("Nastran shell identities, work and memory obey the supplied plan",
          "[IO4-03][N4-02][negative][verification]") {
    const auto source = std::string{grids} + std::string{triangle} +
                        "PSHELL,20\nCBEAM,21\nCONM2,22\nRBE2,23\nPBEAM,24\n";
    const auto attempt = [&](NastranShellLimits bounds, core::ResourcePlan& plan,
                             const std::stop_token& stop = {}) {
        auto input = coordinates(source);
        auto result = NastranShellCards::prepare(std::move(input), bounds, plan, stop);
        if (!result) {
            // NOLINTNEXTLINE(bugprone-use-after-move): rejected calls preserve input.
            CHECK(take(input.view()).bulk.bytes == source);
        }
        return result;
    };
    core::ResourcePlan unused{0};
    const auto no_memory = attempt(shell_limits, unused);
    REQUIRE_FALSE(no_memory);
    CHECK(no_memory.error().code == ResourceExceeded);
    core::ResourcePlan full{memory};
    const auto too_many_elements = attempt({3, 2, shell_limits.work}, full);
    REQUIRE_FALSE(too_many_elements);
    CHECK(too_many_elements.error().code == ResourceExceeded);
    const auto too_many_properties = attempt({4, 1, shell_limits.work}, full);
    REQUIRE_FALSE(too_many_properties);
    CHECK(too_many_properties.error().code == ResourceExceeded);

    // NOLINTNEXTLINE(misc-const-correctness) MSVC request_stop() is non-const.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    const auto cancelled = attempt(shell_limits, full, cancel.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK(full.remaining_bytes() == memory);
    const auto accepted = attempt({4, 2, shell_limits.work}, full);
    REQUIRE(accepted);
    const auto work = take(accepted->view()).work;
    REQUIRE(work > 0);
    const auto bytes = memory - full.remaining_bytes();
    REQUIRE(bytes > 0);
    core::ResourcePlan exact{bytes};
    const auto verified = attempt({4, 2, work}, exact);
    CHECK(verified.has_value());
    CHECK(exact.remaining_bytes() == 0);
    core::ResourcePlan short_memory{bytes - 1};
    const auto byte_rejected = attempt({4, 2, work}, short_memory);
    REQUIRE_FALSE(byte_rejected);
    CHECK(byte_rejected.error().code == ResourceExceeded);
    CHECK(short_memory.remaining_bytes() == bytes - 1);
    core::ResourcePlan short_work{memory};
    const auto rejected = attempt({4, 2, work - 1}, short_work);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == ResourceExceeded);
    CHECK(short_work.remaining_bytes() < memory);
}

TEST_CASE("Nastran structural supports preserve prescribed motion, THRU and set membership",
          "[IO4-03][verification]") {
    const auto source = std::string{grids} + "SPC,10,17,135,-2.5D+1,8,46,\n"
                                             "SPC1,30,26,7,THRU,24\n"
                                             "SPC1,30,126,41,22\n"
                                             "SPC1,100,1,17\nSPCADD,100,30,10\n";
    const auto parent = coordinates(source);
    const auto parent_view = take(parent.view());
    core::ResourcePlan plan{memory};
    auto resolved = take(NastranSupportCards::prepare(parent, {6, 2, 10000}, plan, {}));
    const auto initial = take(resolved.view());
    CHECK(initial.source.digest == parent_view.bulk.digest);
    CHECK(initial.source.location.source == parent_view.bulk.location.source);
    REQUIRE(initial.terms.size() == 6);
    REQUIRE(initial.union_members.size() == 2);
    CHECK(initial.terms[0].set == NastranSupportSetId{10});
    CHECK(initial.terms[0].components == 0b010101);
    const auto& point = std::get<NastranSupportPoint>(initial.terms[0].target);
    CHECK(point.grid == NastranGridId{17});
    CHECK(parent_view.grids[point.grid_ordinal].id == point.grid);
    const auto prescribed = take(initial.terms[0].prescribed_motion);
    CHECK(prescribed.decimal() == "-2.5e+1");
    CHECK(parent_view.bulk.text(take(initial.terms[0].motion_field)) == "-2.5D+1");
    CHECK(initial.terms[1].components == 0b101000);
    CHECK_FALSE(initial.terms[1].prescribed_motion);
    CHECK(parent_view.bulk.text(take(initial.terms[1].motion_field)).empty());
    const auto& range = std::get<NastranSupportRange>(initial.terms[2].target);
    CHECK(range.first == NastranGridId{7});
    CHECK(range.last == NastranGridId{24});
    CHECK(parent_view.bulk.text(range.thru_field) == "THRU");
    CHECK(std::get<NastranSupportPoint>(initial.terms[3].target).grid == NastranGridId{41});
    CHECK(std::get<NastranSupportPoint>(initial.terms[4].target).grid == NastranGridId{22});
    CHECK(initial.terms[5].set == NastranSupportSetId{100});
    CHECK(initial.union_members[0].combination == NastranSupportSetId{100});
    CHECK(initial.union_members[0].member == NastranSupportSetId{30});
    CHECK(initial.union_members[1].member == NastranSupportSetId{10});
    CHECK(take(parent.view()).bulk.bytes == source);

    const auto moved = std::move(resolved);
    // NOLINTNEXTLINE(bugprone-use-after-move): moving invalidates the old support view.
    CHECK_FALSE(resolved.view());
    CHECK(take(moved.view()).terms.size() == 6);

    core::ResourcePlan exact{plan.allocated_bytes()};
    CHECK(NastranSupportCards::prepare(parent, {6, 2, initial.work}, exact, {}));
    core::ResourcePlan constrained{plan.allocated_bytes() - 1};
    const auto short_memory =
        NastranSupportCards::prepare(parent, {6, 2, initial.work}, constrained, {});
    REQUIRE_FALSE(short_memory);
    CHECK(short_memory.error().code == ResourceExceeded);
    core::ResourcePlan enough{memory};
    const auto short_work =
        NastranSupportCards::prepare(parent, {6, 2, initial.work - 1}, enough, {});
    REQUIRE_FALSE(short_work);
    CHECK(short_work.error().code == ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop is non-const.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    core::ResourcePlan no_memory{0};
    const auto cancelled =
        NastranSupportCards::prepare(parent, {6, 2, initial.work}, no_memory, stop.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK(no_memory.allocated_bytes() == 0);
    CHECK(take(parent.view()).bulk.bytes == source);
}

TEST_CASE("Nastran structural supports diagnose invalid references before returning a set",
          "[IO4-03][negative]") {
    using SupportIssue = NastranSupportFailure;
    std::string_view input;
    SupportIssue reason{};
    std::size_t expected_line = 5;
    std::optional<std::size_t> expected_column;
    SECTION("unknown explicit point") {
        input = "SPC1,10,2,12";
        reason = SupportIssue::UnknownGrid;
    }
    SECTION("repeated component digit") {
        input = "SPC,10,17,122,0.0";
        reason = SupportIssue::InvalidField;
    }
    SECTION("invalid range direction") {
        input = "SPC1,10,123,24,THRU,7";
        reason = SupportIssue::InvalidField;
    }
    SECTION("invalid prescribed real") {
        input = "SPC,10,17,1,TRU";
        reason = SupportIssue::InvalidField;
    }
    SECTION("small-width free integer fields have at most eight characters") {
        input = "SPC1,123456789,1,17";
        reason = SupportIssue::InvalidField;
        expected_column = 6;
    }
    SECTION("unsupported nesting") {
        input = "SPCADD,30,20\nSPCADD,20,10\nSPC1,10,2,17";
        reason = SupportIssue::NestedCombination;
    }
    SECTION("unresolved union member") {
        input = "SPCADD,30,10";
        reason = SupportIssue::UnknownSet;
    }
    SECTION("duplicate union member") {
        input = "SPCADD,30,10,10\nSPC1,10,2,17";
        reason = SupportIssue::DuplicateMember;
        expected_column = 14;
    }
    SECTION("duplicate set combination") {
        input = "SPCADD,30,20\nSPCADD,30,10\nSPC1,10,2,17\nSPC1,20,2,8";
        reason = SupportIssue::DuplicateCombination;
        expected_line = 6;
    }
    const auto source = std::string{grids} + std::string{input};
    const auto parent = coordinates(source);
    core::ResourcePlan plan{memory};
    std::optional<NastranSupportDiagnostic> diagnostic;
    const auto result = NastranSupportCards::prepare(parent, {8, 8, 10000}, plan, {}, &diagnostic);
    REQUIRE_FALSE(result);
    const auto recorded = take(diagnostic);
    CHECK(recorded.reason == reason);
    CHECK(recorded.position.line == expected_line);
    if (expected_column)
        CHECK(recorded.position.column == *expected_column);
    CHECK(take(parent.view()).bulk.bytes == source);
}

} // namespace
