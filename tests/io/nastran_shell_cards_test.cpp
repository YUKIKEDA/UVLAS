#include "../../src/structure/exact_polynomial_detail.hpp"
#include "nastran_lumped_mass_geometry.hpp"
#include "nastran_rigid_attachment.hpp"
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

TEST_CASE("CONM2 retains local offsets, basic centres and named inertia products",
          "[IO4-03][verification]") {
    const auto source =
        "CONM2*,101,17,9,12.5\n*,1.25D+1,-2.0,3.0\n" +
        fixed("CONM2", {"102", "41", "-1", "2.75", "4.0", "-5.0", "6.0"}, 8, "+M2") +
        fixed("+M2", {"11.0", "-21.0", "22.0", "31.0", "-32.0", "33.0"}, 8) +
        "CONM2,103,8,,0.5\nCORD2C,9,0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n" + std::string{grids};
    core::ResourcePlan shell_plan{memory};
    const auto parent =
        take(NastranShellCards::prepare(coordinates(source), shell_limits, shell_plan, {}));
    const auto shell = take(parent.view());
    core::ResourcePlan plan{memory};
    auto parsed = take(NastranLumpedMassCards::prepare(parent, {3, 10000}, plan, {}));
    const auto held = take(parsed.view());
    CHECK(held.source.digest == shell.input.bulk.digest);
    CHECK(held.source.location.source == shell.input.bulk.location.source);
    REQUIRE(held.masses.size() == 3);
    const auto& framed = held.masses[0];
    CHECK(framed.id == NastranElementId{101});
    CHECK(framed.grid == NastranGridId{17});
    CHECK(shell.input.grids[framed.grid_ordinal].id == framed.grid);
    const auto& local = std::get<NastranFramedMassOffset>(framed.frame);
    CHECK(local.frame == NastranCoordinateId{9});
    CHECK(shell.input.coordinates[local.frame_ordinal].kind == NastranCoordinateKind::Cylindrical);
    const auto x = take(framed.position[0].supplied);
    CHECK(x.decimal() == "1.25e+1");
    CHECK_FALSE(framed.inertia.i11.supplied);
    CHECK_FALSE(framed.inertia.i11.field);
    const auto& absolute = held.masses[1];
    CHECK(std::holds_alternative<NastranBasicMassCentre>(absolute.frame));
    CHECK(shell.input.bulk.text(take(absolute.frame_field)) == "-1");
    const auto supplied_mass = take(absolute.mass.supplied);
    CHECK(supplied_mass.decimal() == "2.75");
    const auto& i = absolute.inertia;
    const std::array inertia{&i.i11, &i.i21, &i.i22, &i.i31, &i.i32, &i.i33};
    const std::array<std::string_view, 6> expected{"11.0", "-21.0", "22.0",
                                                   "31.0", "-32.0", "33.0"};
    for (std::size_t n = 0; n < inertia.size(); ++n) {
        const auto value = take(inertia[n]->supplied);
        const auto field = take(inertia[n]->field);
        CHECK(value.decimal() == expected[n]);
        CHECK(field.position == NastranBulkPosition{source.find(expected[n]), 4, 9 + 8 * n});
    }
    const auto& basic = held.masses[2];
    CHECK(std::holds_alternative<NastranBasicMassOffset>(basic.frame));
    CHECK(shell.input.bulk.text(take(basic.frame_field)).empty());
    CHECK_FALSE(basic.position[0].supplied);

    const auto moved = std::move(parsed);
    // NOLINTNEXTLINE(bugprone-use-after-move): moving invalidates the old mass view.
    CHECK_FALSE(parsed.view());
    CHECK(take(moved.view()).masses.size() == 3);
    core::ResourcePlan exact{plan.allocated_bytes()};
    CHECK(NastranLumpedMassCards::prepare(parent, {3, held.work}, exact, {}));
    core::ResourcePlan short_memory{plan.allocated_bytes() - 1};
    const auto memory_limited =
        NastranLumpedMassCards::prepare(parent, {3, held.work}, short_memory, {});
    REQUIRE_FALSE(memory_limited);
    CHECK(memory_limited.error().code == ResourceExceeded);
    const auto work_limited = NastranLumpedMassCards::prepare(parent, {3, held.work - 1}, plan, {});
    REQUIRE_FALSE(work_limited);
    CHECK(work_limited.error().code == ResourceExceeded);
    CHECK(take(parent.view()).input.bulk.bytes == source);
}

TEST_CASE("CONM2 identifies the field that blocks typed mass preparation", "[IO4-03][negative]") {
    using MassIssue = NastranLumpedMassFailure;
    std::string_view line;
    MassIssue expected{};
    SECTION("unknown GRID") {
        line = "CONM2,1,99,0,1.0";
        expected = MassIssue::UnknownGrid;
    }
    SECTION("missing mass") {
        line = "CONM2,1,17,0,";
        expected = MassIssue::MissingField;
    }
    SECTION("unknown coordinate") {
        line = "CONM2,1,17,99,1.0";
        expected = MassIssue::UnknownFrame;
    }
    SECTION("invalid negative coordinate") {
        line = "CONM2,1,17,-2,1.0";
        expected = MassIssue::InvalidField;
    }
    SECTION("an integer spelling is not a mass value") {
        line = "CONM2,1,17,0,10";
        expected = MassIssue::InvalidField;
    }
    SECTION("a word in the reserved first row") {
        line = "CONM2,1,17,0,1.0,,,,3.0";
        expected = MassIssue::ExtraField;
    }
    SECTION("a word beyond the last inertia slot") {
        line = "CONM2,1,17,0,1.0,,,,,1.0,,2.0,,,3.0,4.0";
        expected = MassIssue::ExtraField;
    }
    const auto source = std::string{grids} + std::string{line};
    core::ResourcePlan shell_plan{memory}, mass_plan{memory};
    const auto parent =
        take(NastranShellCards::prepare(coordinates(source), shell_limits, shell_plan, {}));
    std::optional<NastranLumpedMassDiagnostic> diagnostic;
    const auto result =
        NastranLumpedMassCards::prepare(parent, {3, 10000}, mass_plan, {}, &diagnostic);
    REQUIRE_FALSE(result);
    const auto found = take(diagnostic);
    CHECK(found.reason == expected);
    CHECK(found.position.line == 5);
    CHECK(take(parent.view()).input.bulk.bytes == source);
}

TEST_CASE("RBE2 separates the ordered GRID list from thermal terms", "[IO4-03][verification]") {
    const auto source = std::string{"RBE2,101,17,136,22,8,6.5-6,-273.15\n"} +
                        fixed("RBE2", {"102", "17", "45", "41", "8", "22"}, 8) +
                        fixed("RBE2", {"103", "17", "2", "22", "8", "41", "5", "6"}, 8, "+R1") +
                        fixed("+R1", {"0.0", "15.0"}, 8) + std::string{grids} + "GRID,5\nGRID,6\n";
    core::ResourcePlan shell_plan{memory}, plan{memory};
    const auto parent =
        take(NastranShellCards::prepare(coordinates(source), shell_limits, shell_plan, {}));
    const auto shell = take(parent.view());
    auto prepared = take(NastranRigidLinkCards::prepare(parent, {3, 10, 10000}, plan, {}));
    const auto view = take(prepared.view());
    CHECK(view.source.digest == shell.input.bulk.digest);
    CHECK(view.source.location.source == shell.input.bulk.location.source);
    REQUIRE(view.links.size() == 3);
    REQUIRE(view.dependents.size() == 10);
    const auto& first = view.links[0];
    CHECK(first.id == NastranElementId{101});
    CHECK(first.dependent_components == 0b100101);
    CHECK(first.independent.id == NastranGridId{17});
    CHECK(shell.input.grids[first.independent.ordinal].id == first.independent.id);
    CHECK(first.first_dependent == 0);
    CHECK(first.dependent_count == 2);
    CHECK(view.dependents[0].id == NastranGridId{22});
    CHECK(view.dependents[1].id == NastranGridId{8});
    const auto alpha = take(first.alpha);
    const auto tref = take(first.reference_temperature);
    CHECK(alpha.value.decimal() == "6.5e-6");
    CHECK(tref.value.decimal() == "-273.15");
    CHECK(shell.input.bulk.text(alpha.field) == "6.5-6");
    CHECK(view.links[1].first_dependent == 2);
    CHECK(view.links[1].dependent_count == 3);
    CHECK(view.links[1].independent.id == NastranGridId{17});
    CHECK_FALSE(view.links[1].alpha);
    CHECK_FALSE(view.links[1].reference_temperature);
    CHECK(view.links[2].first_dependent == 5);
    CHECK(view.links[2].dependent_count == 5);
    CHECK(view.links[2].dependent_components == 0b000010);
    CHECK(view.dependents[8].id == NastranGridId{5});
    CHECK(view.dependents[9].id == NastranGridId{6});
    const auto fixed_tref = take(view.links[2].reference_temperature);
    CHECK(fixed_tref.field.position == NastranBulkPosition{source.find("15.0"), 4, 17});
    CHECK(fixed_tref.value.decimal() == "15.0");
    for (const auto& point : view.dependents)
        CHECK(shell.input.grids[point.ordinal].id == point.id);

    const auto moved = std::move(prepared);
    // NOLINTNEXTLINE(bugprone-use-after-move): moving invalidates the old rigid-link view.
    CHECK_FALSE(prepared.view());
    CHECK(take(moved.view()).links.size() == 3);
    core::ResourcePlan exact{plan.allocated_bytes()};
    CHECK(NastranRigidLinkCards::prepare(parent, {3, 10, view.work}, exact, {}));
    core::ResourcePlan short_memory{plan.allocated_bytes() - 1};
    const auto constrained =
        NastranRigidLinkCards::prepare(parent, {3, 10, view.work}, short_memory, {});
    REQUIRE_FALSE(constrained);
    CHECK(constrained.error().code == ResourceExceeded);
    const auto short_members = NastranRigidLinkCards::prepare(parent, {3, 9, view.work}, plan, {});
    REQUIRE_FALSE(short_members);
    CHECK(short_members.error().code == ResourceExceeded);
    CHECK(take(parent.view()).input.bulk.bytes == source);
}

TEST_CASE("RBE2 rejects ambiguous thermal fields and inconsistent GRID references",
          "[IO4-03][negative]") {
    using RigidIssue = NastranRigidLinkFailure;
    std::string_view line;
    RigidIssue expected{};
    auto code = InvalidInput;
    SECTION("unknown independent GRID") {
        line = "RBE2,11,99,1,22";
        expected = RigidIssue::UnknownGrid;
    }
    SECTION("unknown dependent GRID") {
        line = "RBE2,11,17,1,99";
        expected = RigidIssue::UnknownGrid;
    }
    SECTION("independent GRID reused as a dependent") {
        line = "RBE2,11,17,1,17";
        expected = RigidIssue::RepeatedGrid;
    }
    SECTION("dependent GRID repeated") {
        line = "RBE2,11,17,1,22,8,22";
        expected = RigidIssue::RepeatedGrid;
    }
    SECTION("component repeated") {
        line = "RBE2,11,17,66,22";
        expected = RigidIssue::InvalidField;
    }
    SECTION("no dependent GRID") {
        line = "RBE2,11,17,1,0.0";
        expected = RigidIssue::MissingField;
    }
    SECTION("temperature without a spelled ALPHA") {
        line = "RBE2,11,17,1,22,,15.0";
        expected = RigidIssue::UnsupportedThermalForm;
        code = Unsupported;
    }
    SECTION("more than two thermal fields") {
        line = "RBE2,11,17,1,22,0.0,1.0,2.0";
        expected = RigidIssue::ExtraField;
    }
    const auto source = std::string{grids} + std::string{line};
    core::ResourcePlan shell_plan{memory}, plan{memory};
    const auto parent =
        take(NastranShellCards::prepare(coordinates(source), shell_limits, shell_plan, {}));
    std::optional<NastranRigidLinkDiagnostic> diagnostic;
    const auto result =
        NastranRigidLinkCards::prepare(parent, {3, 10, 10000}, plan, {}, &diagnostic);
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
    const auto failure = take(diagnostic);
    CHECK(failure.reason == expected);
    CHECK(failure.position.line == 5);
    CHECK(take(parent.view()).input.bulk.bytes == source);
}

using MassQ = structure::detail::Rational;
using MassVector = std::array<MassQ, 3>;
using MassGeometryReason = NastranMassGeometryFailure;
constexpr auto mass_geometry_budget = std::uint64_t{1} << 28U;

[[nodiscard]] NastranMassGeometryLimits mass_geometry_limits() noexcept {
    const structure::ExactLimits exact{512,
                                       8,
                                       32,
                                       8,
                                       8,
                                       4096,
                                       8192,
                                       structure::ExactEstimateAdmission{mass_geometry_budget},
                                       structure::ExactTimingMode::Disabled};
    return {{exact, 30000, 30000, {}}, 64000};
}

[[nodiscard]] std::string mass_geometry_fixture() {
    return "CORD2R,9,0,100.0,200.0,300.0,100.0,200.0,301.0,101.0,200.0,300.0\n"
           "GRID,17,9,0.0,0.0,0.0,9\n" +
           fixed("CONM2", {"101", "17", "0", "2000.0", "300.0", "-400.0", "100.0"}, 8, "+OFFSET") +
           fixed("+OFFSET", {"1.4E8", "2.0E7", "1.0E8", "", "-3.0E7", "6.0E7"}) +
           fixed("CONM2", {"102", "17", "-1", "2000.0", "400.0", "-200.0", "400.0"}, 8, "+CENTRE") +
           fixed("+CENTRE", {"1.4E8", "2.0E7", "1.0E8", "", "-3.0E7", "6.0E7"}) +
           "CONM2,103,17,,500.0\nCONM2,104,17,9,1.0\nCONM2,105,17,0,-1.0\n" +
           fixed("CONM2", {"106", "17", "0", "1.0"}, 8, "+SLACK") +
           fixed("+SLACK", {"10.0", "", "1.0", "", "", "1.0"}) +
           fixed("CONM2", {"107", "17", "0", "1.0"}, 8, "+SPD") +
           fixed("+SPD", {"1.0", "2.0", "1.0", "", "", "1.0"}) +
           "CONM2,108,17,-1,1.0,400.0,,400.0\n";
}

constexpr SiUnits mass_geometry_units{SiLengthUnit::Millimetre, SiMassUnit::Gram,
                                      SiTimeUnit::Second, SiAngleUnit::Degree};

[[nodiscard]] const MassQ& mass_q(const structure::ExactRational& r) {
    return structure::detail::ExactFactory::value(r);
}

void si_contains(const NastranCartesianInterval& enclosure, const MassQ& exact) {
    const auto& lower = mass_q(enclosure.lower);
    const auto& upper = mass_q(enclosure.upper);
    CHECK(lower <= exact);
    CHECK(upper >= exact);
    CHECK(upper - lower < MassQ("1/10000000000000000000000000000000000000000"));
}

[[nodiscard]] MassVector mv_add(const MassVector& a, const MassVector& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

[[nodiscard]] MassVector mv_cross(const MassVector& a, const MassVector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

[[nodiscard]] MassQ mv_dot(const MassVector& a, const MassVector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

TEST_CASE("CONM2 basic offset and absolute centre preserve independent rigid kinetic energy",
          "[IO4-03][ST-07][R-05][integration][verification]") {
    const auto source = mass_geometry_fixture();
    const auto projections = [&] {
        core::ResourcePlan plan{mass_geometry_budget};
        const auto deck =
            take(NastranShellCards::prepare(coordinates(source), shell_limits, plan, {}));
        const auto parsed = take(NastranLumpedMassCards::prepare(deck, {8, 10000}, plan, {}));
        const auto prepare = [&](std::uint64_t id) {
            return take(prepare_nastran_basic_mass(deck, parsed, NastranElementId{id},
                                                   mass_geometry_units, mass_geometry_limits(),
                                                   plan));
        };
        return std::array{prepare(101), prepare(102), prepare(103)};
    }();

    const MassVector grid{MassQ("1/10"), MassQ("1/5"), MassQ("3/10")};
    const MassVector centre{MassQ("2/5"), MassQ("-1/5"), MassQ("2/5")};
    const MassVector lever{MassQ("3/10"), MassQ("-2/5"), MassQ("1/10")};
    // Four equal point masses about their centre. The second pair's yz product
    // is negative, so the two nonzero card products have opposite signs.
    const std::array particles{MassVector{MassQ("1/10"), MassQ("1/5"), MassQ{0}},
                               MassVector{MassQ("-1/10"), MassQ("-1/5"), MassQ{0}},
                               MassVector{MassQ{0}, MassQ("1/10"), MassQ("-3/10")},
                               MassVector{MassQ{0}, MassQ("-1/10"), MassQ("3/10")}};
    const MassVector omega{MassQ{2}, MassQ{-1}, MassQ{3}};
    const MassVector velocity{MassQ("-1/5"), MassQ("2/5"), MassQ("1/10")};
    const MassVector acceleration{MassQ("1/10"), MassQ("1/5"), MassQ("3/10")};
    const MassVector alpha{MassQ{0}, MassQ{2}, MassQ{1}};
    MassQ particle_energy{}, particle_work{};
    for (const auto& relative : particles) {
        const auto position = mv_add(lever, relative);
        const auto v = mv_add(velocity, mv_cross(omega, position));
        const auto a = mv_add(acceleration, mv_cross(alpha, position));
        particle_energy += mv_dot(v, v) / 4; // Each particle has mass 1/2 kg.
        particle_work += mv_dot(v, a) / 2;
    }
    const std::array expected_inertia{MassQ("7/50"), MassQ("-1/50"), MassQ("1/10"),
                                      MassQ{0},      MassQ("3/100"), MassQ("3/50")};
    CHECK(std::holds_alternative<NastranBasicMassOffset>(projections[0].card.frame));
    CHECK(std::holds_alternative<NastranBasicMassCentre>(projections[1].card.frame));
    for (const auto& row : std::span{projections}.first(2)) {
        CHECK(row.kind == NastranBasicMassKind::Rigid);
        CHECK(row.source.digest == source_digest(std::as_bytes(std::span{source})));
        CHECK(row.card.grid == NastranGridId{17});
        CHECK(row.input_units.mass == SiMassUnit::Gram);
        CHECK(mass_q(row.mass_kg) == 2);
        const auto original_i21 = take(row.card.inertia.i21.supplied);
        const auto original_i32 = take(row.card.inertia.i32.supplied);
        CHECK(original_i21.decimal() == "2.0e7");
        CHECK(original_i32.decimal() == "-3.0e7");
        for (std::size_t i = 0; i < grid.size(); ++i) {
            si_contains(row.grid_metres[i], grid[i]);
            si_contains(row.centre_metres[i], centre[i]);
            si_contains(row.grid_to_centre_metres[i], lever[i]);
        }
        for (std::size_t i = 0; i < expected_inertia.size(); ++i) {
            CHECK(mass_q(row.centre_inertia_kg_m2[i]) == expected_inertia[i]);
        }

        // Compare the mass bilinear form with particle velocities/accelerations,
        // both on the diagonal (twice the energy) and with an independent pair.
        const auto compare_grid_work = [&](const MassVector& linear, const MassVector& angular,
                                           const MassQ& expected) {
            const auto translation =
                2 * (mv_dot(velocity, linear) + mv_dot(velocity, mv_cross(angular, lever)) +
                     mv_dot(linear, mv_cross(omega, lever)));
            MassQ lower = translation, upper = translation;
            std::size_t entry = 0;
            for (std::size_t i = 0; i < omega.size(); ++i) {
                for (std::size_t j = 0; j <= i; ++j) {
                    const MassQ weight = i == j ? omega[i] * angular[i]
                                                : omega[i] * angular[j] + omega[j] * angular[i];
                    const auto& value = row.grid_inertia_kg_m2[entry++];
                    lower += weight * mass_q(weight < 0 ? value.upper : value.lower);
                    upper += weight * mass_q(weight < 0 ? value.lower : value.upper);
                }
            }
            CHECK(lower <= expected);
            CHECK(upper >= expected);
            CHECK(upper - lower < MassQ("1/10000000000000000000000000000000000000000"));
        };
        compare_grid_work(velocity, omega, 2 * particle_energy);
        compare_grid_work(acceleration, alpha, particle_work);
    }
    const auto& point_mass = projections[2];
    CHECK(point_mass.kind == NastranBasicMassKind::Particle);
    CHECK(mass_q(point_mass.mass_kg) == MassQ("1/2"));
    for (const auto& value : point_mass.grid_inertia_kg_m2) {
        si_contains(value, MassQ{0});
    }
}

TEST_CASE("CONM2 SI preparation rejects unsupported frames and unphysical inertia",
          "[IO4-03][ST-07][integration][negative][verification]") {
    core::ResourcePlan plan{mass_geometry_budget};
    const auto source = mass_geometry_fixture();
    const auto deck = take(NastranShellCards::prepare(coordinates(source), shell_limits, plan, {}));
    const auto masses = take(NastranLumpedMassCards::prepare(deck, {8, 10000}, plan, {}));
    const auto cards = take(masses.view());
    std::optional<NastranMassGeometryDiagnostic> diagnostic;
    const auto prepare = [&](std::uint64_t id, const NastranMassGeometryLimits& bounded) {
        return prepare_nastran_basic_mass(deck, masses, NastranElementId{id}, mass_geometry_units,
                                          bounded, plan, &diagnostic);
    };
    const auto framed = prepare(104, mass_geometry_limits());
    REQUIRE_FALSE(framed);
    CHECK(framed.error().code == Unsupported);
    CHECK(take(diagnostic).reason == MassGeometryReason::Frame);
    CHECK(take(diagnostic).position == take(cards.masses[3].frame_field).position);
    CHECK(take(diagnostic).coordinate == NastranCoordinateId{9});
    const auto negative = prepare(105, mass_geometry_limits());
    REQUIRE_FALSE(negative);
    CHECK(negative.error().code == InvalidInput);
    CHECK(take(diagnostic).reason == MassGeometryReason::Mass);
    CHECK(take(diagnostic).position == take(cards.masses[4].mass.field).position);
    for (const auto id : {106U, 107U}) {
        const auto impossible = prepare(id, mass_geometry_limits());
        REQUIRE_FALSE(impossible);
        CHECK(impossible.error().code == InvalidInput);
        CHECK(take(diagnostic).reason == MassGeometryReason::Inertia);
    }

    const auto incomplete_centre = prepare(108, mass_geometry_limits());
    REQUIRE_FALSE(incomplete_centre);
    CHECK(incomplete_centre.error().code == Unsupported);
    CHECK(take(diagnostic).reason == MassGeometryReason::Position);
    CHECK(take(diagnostic).position == take(cards.masses[7].position[1].field).position);

    const auto measured = take(prepare(101, mass_geometry_limits()));
    auto short_work = mass_geometry_limits();
    short_work.maximum_directed_work = measured.directed_work - 1;
    const auto denied = prepare(101, short_work);
    REQUIRE_FALSE(denied);
    CHECK(denied.error().code == ResourceExceeded);
    ++short_work.maximum_directed_work;
    REQUIRE(prepare(101, short_work));
    CHECK_FALSE(diagnostic);
    CHECK(take(deck.view()).input.bulk.bytes == source);

    const auto other = take(NastranShellCards::prepare(
        coordinates(source + "$ different snapshot\n"), shell_limits, plan, {}));
    const auto stale =
        prepare_nastran_basic_mass(other, masses, NastranElementId{101}, mass_geometry_units,
                                   mass_geometry_limits(), plan, &diagnostic);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == StaleGeneration);
    CHECK(take(diagnostic).reason == MassGeometryReason::Source);
}

constexpr std::string_view rigid_geometry_fixture =
    "CORD2C,20,0,1000.0,2000.0,3000.0,1000.0,2000.0,3100.0,1100.0,2000.0,3000.0\n"
    "GRID,8,0,1000.0,2000.0,3000.0\n"
    "GRID,10,20,2000.0,90.0,0.0,20\n"
    "GRID,12,0,1000.0,2000.0,4000.0\n"
    "GRID,13,20,0.0,40.0,2000.0,20\n"
    "RBE2,200,8,123456,10,12,0.0,293.15\n"
    "RBE2,201,8,26,10,12\n"
    "RBE2,202,8,1,10,0.00001,293.15\n"
    "RBE2,203,8,36,13\n"
    "RBE2,204,8,13,13\n";

[[nodiscard]] NastranRigidAttachmentLimits rigid_limits() noexcept {
    const auto base = mass_geometry_limits();
    return {base.arithmetic, 150000, 4};
}

[[nodiscard]] NastranRigidAttachment rigid_reference(std::uint64_t id, core::ResourcePlan& plan) {
    const auto shells = take(
        NastranShellCards::prepare(coordinates(rigid_geometry_fixture), shell_limits, plan, {}));
    const auto links = take(NastranRigidLinkCards::prepare(shells, {8, 16, 20000}, plan, {}));
    return take(NastranRigidAttachment::prepare(shells, links, NastranElementId{id},
                                                SiLengthUnit::Millimetre, rigid_limits(), plan));
}

struct RigidFixture {
    core::ResourcePlan& plan;

    [[nodiscard]] structure::ExactRational integer(std::string_view text) const {
        return take(structure::ExactRational::from_integer_pair(
            text, "1", rigid_limits().arithmetic.algebra, plan));
    }

    [[nodiscard]] NastranRigidExactVector vector(std::string_view x, std::string_view y,
                                                 std::string_view z) const {
        return {integer(x), integer(y), integer(z)};
    }

    [[nodiscard]] std::array<structure::ExactRational, 9> rotation() const {
        return {integer("0"), integer("-1"), integer("0"), integer("1"), integer("0"),
                integer("0"), integer("0"),  integer("0"), integer("1")};
    }
};

void rigid_contains(const NastranRigidIntervalVector& values, std::array<int, 3> expected) {
    for (std::size_t i = 0; i < expected.size(); ++i) {
        si_contains(values[i], MassQ{expected[i]});
    }
}

void rigid_lanes(const std::array<std::optional<NastranCartesianInterval>, 3>& values,
                 const std::array<std::optional<int>, 3>& expected) {
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(values[i].has_value() == expected[i].has_value());
        if (const auto& actual = values[i]; actual && expected[i]) {
            si_contains(*actual, MassQ{*expected[i]});
        }
    }
}

TEST_CASE("RBE2 projects all selected point motions and their world wrench transpose",
          "[IO4-03][CN-D02][CN-D03][R-05][integration][verification]") {
    core::ResourcePlan plan{mass_geometry_budget};
    const auto owned = rigid_reference(200, plan);
    const auto reference = take(owned.view());
    CHECK(reference.source.digest ==
          source_digest(std::as_bytes(std::span{rigid_geometry_fixture})));
    CHECK(reference.card.dependent_components == 0b111111);
    const auto temperature = take(reference.card.reference_temperature);
    CHECK(temperature.value.decimal() == "293.15");
    rigid_contains(reference.master_metres, {1, 2, 3});
    REQUIRE(reference.dependents.size() == 2);
    CHECK(reference.dependents[0].source.id == NastranGridId{10});
    CHECK(reference.dependents[1].source.id == NastranGridId{12});
    rigid_contains(reference.dependents[0].offset_from_master_metres, {0, 2, 0});
    rigid_contains(reference.dependents[1].offset_from_master_metres, {0, 0, 1});

    const RigidFixture r{plan};
    const auto origin = r.vector("10", "20", "30");
    const auto v = r.vector("1", "2", "3");
    const auto w = r.vector("2", "-1", "4");
    const auto rotation = r.rotation();
    const NastranRigidMotion motion{origin, rotation, v, w};
    const std::array loads{
        NastranRigidPointLoad{NastranGridId{10}, r.vector("2", "3", "4"), r.vector("5", "6", "7")},
        NastranRigidPointLoad{NastranGridId{12}, r.vector("1", "-2", "3"),
                              r.vector("-4", "5", "2")}};
    const auto projected =
        take(evaluate_nastran_rigid_attachment(owned, motion, loads, rigid_limits(), plan));
    CHECK(projected.link == NastranElementId{200});
    const auto rows = projected.dependents.values();
    REQUIRE(rows.size() == 2);
    rigid_contains(rows[0].rigid_target_metres, {8, 20, 30});
    rigid_contains(rows[1].rigid_target_metres, {10, 20, 31});
    // CD at GRID 10 is (Y,-X,Z). The current lever is (-2,0,0), so
    // omega x lever is (0,-8,-2). GRID 12 uses basic axes.
    rigid_lanes(rows[0].translation_target, {20, -8, 30});
    rigid_lanes(rows[0].translational_velocity, {-6, -1, 1});
    rigid_lanes(rows[0].spin, {-1, -2, 4});
    rigid_lanes(rows[1].translation_target, {10, 20, 31});
    rigid_lanes(rows[1].translational_velocity, {0, 0, 3});
    rigid_lanes(rows[1].spin, {2, -1, 4});
    if (const auto& orientation = rows[0].rigid_target_orientation) {
        rigid_contains((*orientation)[0], {-1, 0, 0});
        rigid_contains((*orientation)[1], {0, -1, 0});
        rigid_contains((*orientation)[2], {0, 0, 1});
    } else {
        FAIL("A complete rotational mask lost its finite target frame");
    }
    // Independent sum: F=(-3,2,4)+(1,-2,3); M=(-6,13,3)+(-2,6,2).
    rigid_contains(projected.world_force_n, {-2, 0, 7});
    rigid_contains(projected.world_moment_about_master_nm, {-8, 19, 5});

    // Compare the master-side power enclosure to the independent sum of the
    // six dependent per-GRID component products, first 0 and then 4 watts.
    MassQ lower{}, upper{};
    const auto accumulate = [&](const NastranRigidIntervalVector& wrench,
                                const NastranRigidExactVector& motion_vector) {
        for (std::size_t i = 0; i < wrench.size(); ++i) {
            const auto& q = mass_q(motion_vector[i]);
            lower += q * mass_q(q < 0 ? wrench[i].upper : wrench[i].lower);
            upper += q * mass_q(q < 0 ? wrench[i].lower : wrench[i].upper);
        }
    };
    accumulate(projected.world_force_n, v);
    accumulate(projected.world_moment_about_master_nm, w);
    CHECK(lower <= 4);
    CHECK(upper >= 4);
    CHECK(upper - lower < MassQ("1/10000000000000000000000000000000000000000"));
}

TEST_CASE("RBE2 partial selection preserves independent directions and diagnoses unsupported input",
          "[IO4-03][CN-D03][integration][negative][verification]") {
    core::ResourcePlan plan{mass_geometry_budget};
    const auto partial = rigid_reference(201, plan);
    const RigidFixture r{plan};
    const auto origin = r.vector("10", "20", "30");
    const auto v = r.vector("1", "2", "3");
    const auto w = r.vector("2", "-1", "4");
    const auto rotation = r.rotation();
    const NastranRigidMotion motion{origin, rotation, v, w};
    std::array loads{
        NastranRigidPointLoad{NastranGridId{10}, r.vector("0", "3", "0"), r.vector("0", "0", "7")},
        NastranRigidPointLoad{NastranGridId{12}, r.vector("0", "-2", "0"),
                              r.vector("0", "0", "2")}};
    const auto good =
        take(evaluate_nastran_rigid_attachment(partial, motion, loads, rigid_limits(), plan));
    const auto& first = good.dependents.values().front();
    rigid_lanes(first.translation_target, {std::nullopt, -8, std::nullopt});
    rigid_lanes(first.translational_velocity, {std::nullopt, -1, std::nullopt});
    rigid_lanes(first.spin, {std::nullopt, std::nullopt, 4});
    CHECK_FALSE(first.rigid_target_orientation);
    rigid_contains(good.world_force_n, {-3, -2, 0});
    rigid_contains(good.world_moment_about_master_nm, {2, 0, 9});

    std::optional<NastranAttachmentDiagnostic> problem;
    loads[0].force_n[0] = r.integer("1");
    const auto off_component =
        evaluate_nastran_rigid_attachment(partial, motion, loads, rigid_limits(), plan, &problem);
    REQUIRE_FALSE(off_component);
    CHECK(off_component.error().code == InvalidInput);
    CHECK(take(problem).reason == NastranAttachmentFailure::LoadSelection);
    CHECK(take(problem).grid == NastranGridId{10});
    CHECK(take(problem).position == take(partial.view()).card.components_field.position);
    loads[0].force_n[0] = r.integer("0");
    const std::array reflection{r.integer("-1"), r.integer("0"), r.integer("0"),
                                r.integer("0"),  r.integer("1"), r.integer("0"),
                                r.integer("0"),  r.integer("0"), r.integer("1")};
    const auto improper = evaluate_nastran_rigid_attachment(partial, {origin, reflection, v, w},
                                                            loads, rigid_limits(), plan, &problem);
    REQUIRE_FALSE(improper);
    CHECK(improper.error().code == GeometryInvalid);
    CHECK(take(problem).reason == NastranAttachmentFailure::Pose);

    const auto shell = take(
        NastranShellCards::prepare(coordinates(rigid_geometry_fixture), shell_limits, plan, {}));
    const auto links = take(NastranRigidLinkCards::prepare(shell, {8, 16, 20000}, plan, {}));
    const auto prepare = [&](std::uint64_t link, const NastranRigidAttachmentLimits& input_limits) {
        return NastranRigidAttachment::prepare(shell, links, NastranElementId{link},
                                               SiLengthUnit::Millimetre, input_limits, plan,
                                               &problem);
    };
    const auto hot = prepare(202, rigid_limits());
    REQUIRE_FALSE(hot);
    CHECK(hot.error().code == Unsupported);
    CHECK(take(problem).reason == NastranAttachmentFailure::Thermal);
    CHECK(take(take(problem).position).line == 8);
    const auto axial = take(prepare(203, rigid_limits()));
    if (const auto& axis = take(axial.view()).dependents[0].component_axes[2]) {
        rigid_contains(*axis, {0, 0, 1});
    } else {
        FAIL("The cylindrical axial component is defined at the pole");
    }
    const auto pole = prepare(204, rigid_limits());
    REQUIRE_FALSE(pole);
    CHECK(pole.error().code == Unsupported);
    CHECK(take(problem).reason == NastranAttachmentFailure::DisplacementAxes);
    CHECK(take(problem).coordinate == NastranCoordinateId{20});
    CHECK(take(take(problem).position).line == 5);
    auto bounded = rigid_limits();
    bounded.maximum_dependents = 1;
    const auto limited = prepare(200, bounded);
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == ResourceExceeded);
    CHECK(take(shell.view()).input.bulk.bytes == rigid_geometry_fixture);
    const auto other = take(NastranShellCards::prepare(
        coordinates(std::string{rigid_geometry_fixture} + "$ other\n"), shell_limits, plan, {}));
    const auto stale =
        NastranRigidAttachment::prepare(other, links, NastranElementId{200},
                                        SiLengthUnit::Millimetre, rigid_limits(), plan, &problem);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == StaleGeneration);
}

} // namespace
