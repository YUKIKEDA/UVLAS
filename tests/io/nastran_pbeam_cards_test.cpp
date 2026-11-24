#include "nastran_pbeam_cards.hpp"
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
using Reason = NastranPbeamFailure;

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

NastranBulkDeck lexical_deck(std::string_view text, const RootSnapshotLocation& location) {
    const auto bytes = std::as_bytes(std::span{text});
    Reader reader{bytes};
    core::ResourcePlan plan{memory};
    auto file = take(read_reference_file(reader, source_digest(bytes),
                                         {text.size(), text.size() + 1, 3, text.size()}, plan, {},
                                         nullptr));
    AcquiredRootSnapshot root{location, std::move(file), {}};
    return take(NastranBulkDeck::prepare(std::move(root), lexical, plan, {}));
}

constexpr NastranPbeamLimits property_limits{16, 24, 50000};
constexpr std::string_view beam_source =
    "GRID,11\nGRID,12\nGRID,13\nCBEAM,90,,11,12,13\nMAT1,6,2.0,,.3\n";
constexpr RootSnapshotLocation deck_location{
    {{SourceId{90}, 100, {}}, 1, 1}, SourceId{99}, ReferenceRole::NastranDeck, 31};

struct Parents {
    NastranIsotropicShellCards materials;
    NastranBeamCards beams;
};

Parents parents(std::string_view text, const RootSnapshotLocation& location = deck_location) {
    auto deck = lexical_deck(text, location);
    core::ResourcePlan plan{memory};
    auto coordinates = take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
    auto shells =
        take(NastranShellCards::prepare(std::move(coordinates), {16, 16, 20000}, plan, {}));
    auto beams = take(NastranBeamCards::prepare(shells, {16, 16, 20000}, plan, {}));
    auto materials =
        take(NastranIsotropicShellCards::prepare(std::move(shells), {16, 20000}, plan, {}));
    return {std::move(materials), std::move(beams)};
}

std::string fixed(std::string_view name, std::initializer_list<std::string_view> values,
                  std::size_t width = 8, std::string_view label = {}) {
    std::string line{name};
    line.resize(8, ' ');
    for (const auto value : values) {
        REQUIRE(value.size() <= width);
        line.append(value);
        line.append(width - value.size(), ' ');
    }
    line.resize(72, ' ');
    return line + std::string{label} + "\n";
}

std::string group(std::initializer_list<std::string_view> values, std::string_view label = {}) {
    return fixed({}, values, 8, label);
}

TEST_CASE("PBEAM keeps declared stations and distinct shear-center attributes under CBEAM identity",
          "[IO4-03][IO4-01][verification]") {
    const auto source = std::string{beam_source} + "PBEAM,41,6,1.1,2.2,3.3\n"
                                                   "PBEAM,90,6,2.9,3.5,5.97,-.10,0.0,.15\n"
                                                   ",.11,-.12,.21,-.22,.31,-.32,.41,-.42\n"
                                                   ",YESA,.1D+1,5.3,,,,,\n"
                                                   ",NO,.25,,,,,,\n"
                                                   ",YES,5.0-1,4.0,4.2,5.0,.00,,.21\n"
                                                   ",.51,.52,.61,.62,.71,.72,.81,.82\n"
                                                   ",1.1,2.1,.21,,.32,.33,.0,.45\n"
                                                   ",.51,.52,.53,.54,-.55,-.56,-.57,-.58\n";
    const auto input = parents(source);
    const auto incoming = take(input.materials.view());
    const auto links = take(input.beams.view());
    core::ResourcePlan plan{memory};
    auto qualified =
        take(NastranPbeamCards::prepare(input.materials, input.beams, property_limits, plan, {}));
    const auto view = take(qualified.view());
    const auto& bulk = incoming.input.input.bulk;
    REQUIRE(view.properties.size() == 2);
    REQUIRE(links.beams.size() == 1);
    REQUIRE(view.stations.size() == 3);
    CHECK(view.source.digest == source_digest(std::as_bytes(std::span{source})));
    CHECK(view.source.location.source == links.source.location.source);
    CHECK(view.properties[0].declaration.id == NastranPropertyId{41});
    CHECK_FALSE(view.properties[0].end_b_station);
    CHECK(view.sections(view.properties[0]).empty());
    const auto& cbeam = links.beams[0];
    const auto& section = view.properties[cbeam.property.declaration_ordinal];
    CHECK(section.declaration.id == cbeam.property.id);
    CHECK(section.declaration.card_ordinal == links.properties[1].card_ordinal);
    const auto& mat = incoming.materials[section.material.declaration_ordinal];
    CHECK(mat.id == NastranMaterialId{6});
    CHECK(mat.kind == NastranMaterialKind::Mat1);
    CHECK(take(section.end_b_station) == 0);
    const auto check = [&](const NastranPbeamReal& spec, std::string_view original,
                           std::string_view expected) {
        const auto real = take(spec);
        CHECK(real.value.decimal() == expected);
        CHECK(bulk.text(real.field) == original);
    };
    check(section.end_a.area, "2.9", "2.9");
    check(section.end_a.i1, "3.5", "3.5");
    check(section.end_a.i2, "5.97", "5.97");
    check(section.end_a.i12, "-.10", "-0.10");
    check(section.end_a.torsion, "0.0", "0.0");
    check(section.end_a.nonstructural_mass, ".15", "0.15");
    const auto recovery = take(section.recovery_a);
    check(recovery.locations[1], "-.12", "-0.12");
    check(recovery.locations[7], "-.42", "-0.42");
    const auto stations = view.sections(section);
    REQUIRE(stations.size() == 3);
    CHECK(stations[0].output == NastranPbeamOutput::EndAPoints);
    CHECK_FALSE(stations[0].recovery);
    check(stations[0].fraction, ".1D+1", "0.1e+1");
    check(stations[0].section.area, "5.3", "5.3");
    CHECK_FALSE(stations[0].section.i1);
    CHECK(stations[1].output == NastranPbeamOutput::No);
    CHECK_FALSE(stations[1].section.area);
    check(stations[1].fraction, ".25", "0.25");
    CHECK(stations[2].output == NastranPbeamOutput::SectionPoints);
    check(stations[2].fraction, "5.0-1", "5.0e-1");
    const auto fraction = take(stations[2].fraction);
    CHECK(fraction.field.position == NastranBulkPosition{source.find("5.0-1"), 11, 6});
    const auto at_station = take(stations[2].recovery);
    check(at_station.locations[2], ".61", "0.61");
    const auto fields = take(section.modifiers);
    check(fields.shear_factor_1, "1.1", "1.1");
    check(fields.shear_factor_2, "2.1", "2.1");
    check(fields.taper_relief_1, ".21", "0.21");
    CHECK_FALSE(fields.taper_relief_2);
    check(fields.nonstructural_inertia[0], ".32", "0.32");
    check(fields.nonstructural_inertia[1], ".33", "0.33");
    check(fields.warping[0], ".0", "0.0");
    check(fields.warping[1], ".45", "0.45");
    const auto positions = take(section.offsets);
    check(positions.nonstructural_mass_center[0], ".51", "0.51");
    check(positions.nonstructural_mass_center[3], ".54", "0.54");
    check(positions.neutral_axis[0], "-.55", "-0.55");
    check(positions.neutral_axis[3], "-.58", "-0.58");
    CHECK(take(input.materials.view()).input.input.bulk.bytes == source);
    CHECK(take(input.beams.view()).beams[0].property.declaration_ordinal == 1);
    const auto moved = std::move(qualified);
    // NOLINTNEXTLINE(bugprone-use-after-move): moving clears the previous view.
    CHECK_FALSE(qualified.view());
    CHECK(take(moved.view()).stations.size() == 3);
}

TEST_CASE("PBEAM distinguishes omitted end-A recovery and explicit empty YES recovery",
          "[IO4-03][verification]") {
    const auto source = std::string{beam_source} +
                        fixed("PBEAM*", {"90", "6", "1.23456789012D+7", "9.0"}, 16, "*P1") +
                        fixed("*P1", {"12.0"}, 16, "*P2") +
                        fixed("*P2", {"YES", "1.0"}, 16, "*P3") + fixed("*P3", {}, 16, "*P4") +
                        fixed("*P4", {}, 16, "*P5") + fixed("*P5", {}, 16);
    core::ResourcePlan plan{memory};
    const auto typed = [&] {
        const auto input = parents(source);
        return take(
            NastranPbeamCards::prepare(input.materials, input.beams, property_limits, plan, {}));
    }();
    const auto view = take(typed.view());
    REQUIRE(view.properties.size() == 1);
    REQUIRE(view.stations.size() == 1);
    const auto& property = view.properties[0];
    CHECK_FALSE(property.recovery_a);
    CHECK_FALSE(property.modifiers);
    const auto value = take(property.end_a.area);
    CHECK(value.value.decimal() == "1.23456789012e+7");
    CHECK(value.field.position.column == 41);
    CHECK(take(property.end_b_station) == 0);
    const auto recovery = take(view.stations[0].recovery);
    CHECK_FALSE(recovery.locations[0]);
    CHECK(recovery.position.line == 10);
}

TEST_CASE("PBEAM numeric rows without stations start with the end-A recovery group",
          "[IO4-03][verification]") {
    const auto header =
        std::string{beam_source} + fixed("PBEAM", {"90", "6", "2.0", "3.0", "4.0"}, 8, "+A");
    const auto recovery_only = parents(header + fixed("+A", {"1.1", "2.1"}));
    core::ResourcePlan first_plan{memory};
    const auto parsed = take(NastranPbeamCards::prepare(
        recovery_only.materials, recovery_only.beams, property_limits, first_plan, {}));
    const auto first_view = take(parsed.view());
    REQUIRE(first_view.properties.size() == 1);
    const auto& uniform_property = first_view.properties[0];
    CHECK(first_view.sections(uniform_property).empty());
    CHECK_FALSE(uniform_property.modifiers);
    const auto recovery = take(uniform_property.recovery_a);
    const auto first_location = take(recovery.locations[0]);
    CHECK(first_location.value.decimal() == "1.1");

    const auto with_attributes =
        parents(header + fixed("+A", {}, 8, "+K") + fixed("+K", {"1.1", "2.1"}, 8, "+M") +
                fixed("+M", {".25", "-.3"}));
    core::ResourcePlan second_plan{memory};
    const auto typed = take(NastranPbeamCards::prepare(
        with_attributes.materials, with_attributes.beams, property_limits, second_plan, {}));
    const auto second_view = take(typed.view());
    REQUIRE(second_view.properties.size() == 1);
    const auto& property = second_view.properties[0];
    CHECK_FALSE(property.end_b_station);
    CHECK(second_view.sections(property).empty());
    const auto empty_recovery = take(property.recovery_a);
    CHECK_FALSE(empty_recovery.locations[0]);
    const auto modifiers = take(property.modifiers);
    const auto shear = take(modifiers.shear_factor_1);
    CHECK(shear.value.decimal() == "1.1");
    CHECK(shear.field.position.line == 8);
    const auto offsets = take(property.offsets);
    const auto axis = take(offsets.nonstructural_mass_center[1]);
    CHECK(axis.value.decimal() == "-0.3");
    CHECK(axis.field.position.line == 9);
}

TEST_CASE("PBEAM preserves the full allowed station count without assuming declaration order",
          "[IO4-03][verification]") {
    const auto source_prefix =
        std::string{beam_source} + fixed("PBEAM", {"90", "6", "2.0", "3.0", "4.0"});
    auto source = source_prefix + group({"NO", "1.0"});
    for (const auto fraction : {".01", ".09", ".08", ".03", ".07", ".05", ".04", ".02", ".06"})
        source += group({"NO", fraction});
    const auto input = parents(source);
    core::ResourcePlan plan{memory};
    const auto accepted =
        take(NastranPbeamCards::prepare(input.materials, input.beams, {1, 10, 50000}, plan, {}));
    const auto view = take(accepted.view());
    REQUIRE(view.properties.size() == 1);
    const auto rows = view.sections(view.properties[0]);
    REQUIRE(rows.size() == 10);
    CHECK(take(view.properties[0].end_b_station) == 0);
    const auto first_mid = take(rows[1].fraction);
    CHECK(first_mid.value.decimal() == "0.01");
    const auto next_mid = take(rows[2].fraction);
    CHECK(next_mid.value.decimal() == "0.09");
    core::ResourcePlan excess_plan{memory};
    const auto excess_input = parents(source + group({"NO", ".11"}));
    std::optional<NastranPbeamDiagnostic> diagnostic;
    const auto excess = NastranPbeamCards::prepare(excess_input.materials, excess_input.beams,
                                                   property_limits, excess_plan, {}, &diagnostic);
    REQUIRE_FALSE(excess);
    CHECK(take(diagnostic).reason == Reason::ExtraField);
    CHECK(excess_plan.allocated_bytes() == 0);
}

TEST_CASE("PBEAM end-A torsion rejects exact negatives but preserves blank and zero fields",
          "[IO4-03][negative][verification]") {
    constexpr std::string_view prefix = "PBEAM,90,6,2.9,3.5,5.97,,";
    const auto deck = [&](std::string_view j) {
        return std::string{beam_source} + std::string{prefix} + std::string{j} + "\n";
    };
    for (const std::string_view j : {"-1.0", "-1.0D-999", "-1.0,BOGUS"}) {
        CAPTURE(j);
        const auto source = deck(j);
        const auto input = parents(source);
        core::ResourcePlan plan{memory};
        std::optional<NastranPbeamDiagnostic> diagnostic;
        const auto rejected = NastranPbeamCards::prepare(input.materials, input.beams,
                                                         property_limits, plan, {}, &diagnostic);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == InvalidInput);
        const auto error = take(diagnostic);
        CHECK(error.reason == Reason::InvalidField);
        CHECK(error.position == NastranBulkPosition{beam_source.size() + prefix.size(), 6, 26});
    }

    for (const std::string_view j : {"", "0.0", "-0.0D-999"}) {
        CAPTURE(j);
        const auto input = parents(deck(j));
        core::ResourcePlan plan{memory};
        const auto parsed = take(
            NastranPbeamCards::prepare(input.materials, input.beams, property_limits, plan, {}));
        const auto view = take(parsed.view());
        REQUIRE(view.properties.size() == 1);
        const auto& torsion = view.properties[0].end_a.torsion;
        if (j.empty()) {
            CHECK_FALSE(torsion);
        } else {
            const auto value = take(torsion);
            const auto source = take(input.materials.view());
            CHECK(source.input.input.bulk.text(value.field) == j);
            CHECK(value.value.decimal() == (j == "0.0" ? "0.0" : "-0.0e-999"));
        }
    }
}

TEST_CASE("PBEAM diagnoses invalid section fields, unsupported station modes and references",
          "[IO4-03][negative][verification]") {
    struct FailureCase {
        std::string_view part;
        Reason reason;
        core::FailureCode code;
        std::size_t line;
    };
    constexpr std::array cases{
        FailureCase{"PBEAM,90,,2.9,3.5,5.97\n", Reason::MissingField, InvalidInput, 6},
        FailureCase{"PBEAM,90,61,2.9,3.5,5.97\n", Reason::UndefinedMaterial, InvalidInput, 6},
        FailureCase{"PBEAM,90,7,2.9,3.5,5.97\nMAT8,7\n", Reason::WrongMaterialKind, InvalidInput,
                    6},
        FailureCase{"PBEAM,90,6,0.0,3.5,5.97\n", Reason::InvalidField, InvalidInput, 6},
        FailureCase{"PBEAM,90,6,2.9,,5.97\n", Reason::MissingField, InvalidInput, 6},
        FailureCase{"PBEAM,90,6,2.9,3.5,-5.97\n", Reason::InvalidField, InvalidInput, 6},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,MAYBE,1.0\n", Reason::UnsupportedField,
                    Unsupported, 7},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,NO,0.0\n", Reason::InvalidField, InvalidInput, 7},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,NO,1.01\n", Reason::InvalidField, InvalidInput,
                    7},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,NO,.5\n", Reason::MissingEndB, InvalidInput, 6},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,YES,1.0\n", Reason::MissingField, InvalidInput,
                    6},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,NO,1.0,,,,,,\n,YESA,.10D+1\n",
                    Reason::DuplicateStation, InvalidInput, 8},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,NO,.2D-18446744073709551615,,,,,,\n"
                    ",NO,2.0D-18446744073709551616\n",
                    Reason::DuplicateStation, InvalidInput, 8},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,NO,.9D+18446744073709551616\n",
                    Reason::InvalidField, InvalidInput, 7},
        FailureCase{"PBEAM,90,6,2.9,3.5,5.97,,,\n,NO,1.0,,,,,,\n"
                    ",0.0,,,,,,,\n,0.0,,,,,,,\n,0.0\n",
                    Reason::ExtraField, InvalidInput, 10}};
    for (const auto& entry : cases) {
        CAPTURE(entry.part);
        const auto source = std::string{beam_source} + std::string{entry.part};
        const auto input = parents(source);
        core::ResourcePlan plan{memory};
        std::optional<NastranPbeamDiagnostic> diagnostic;
        const auto result = NastranPbeamCards::prepare(input.materials, input.beams,
                                                       property_limits, plan, {}, &diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == entry.code);
        const auto failure = take(diagnostic);
        CHECK(failure.reason == entry.reason);
        CHECK(failure.position.line == entry.line);
        CHECK(take(input.materials.view()).input.input.bulk.bytes == source);
    }
}

TEST_CASE("PBEAM checks source identity and exact resource limits before returning properties",
          "[IO4-03][N4-02][negative][verification]") {
    const auto source = std::string{beam_source} +
                        fixed("PBEAM", {"90", "6", "2.9", "3.5", "5.97"}) + group({"NO", "1.0"});
    const auto input = parents(source);
    core::ResourcePlan plan{memory};
    std::optional<NastranPbeamDiagnostic> diagnostic;
    auto location = deck_location;
    location.source = SourceId{101};
    const auto other = parents(source, location);
    const auto mismatch = NastranPbeamCards::prepare(input.materials, other.beams, property_limits,
                                                     plan, {}, &diagnostic);
    REQUIRE_FALSE(mismatch);
    CHECK(take(diagnostic).reason == Reason::SourceMismatch);
    location = deck_location;
    ++location.declaration_node;
    const auto distinct_declaration = parents(source, location);
    const auto misplaced = NastranPbeamCards::prepare(input.materials, distinct_declaration.beams,
                                                      property_limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(misplaced);
    CHECK(take(diagnostic).reason == Reason::SourceMismatch);
    location = deck_location;
    ++location.parent.contract_revision;
    const auto distinct_revision = parents(source, location);
    const auto old_revision = NastranPbeamCards::prepare(input.materials, distinct_revision.beams,
                                                         property_limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(old_revision);
    CHECK(take(diagnostic).reason == Reason::SourceMismatch);
    const auto different = parents(source + "$ additional bytes\n");
    const auto differs = NastranPbeamCards::prepare(input.materials, different.beams,
                                                    property_limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(differs);
    CHECK(take(diagnostic).reason == Reason::SourceMismatch);
    CHECK(plan.allocated_bytes() == 0);
    const auto limited =
        NastranPbeamCards::prepare(input.materials, input.beams, {1, 0, 50000}, plan, {});
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == ResourceExceeded);
    CHECK(plan.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness) MSVC request_stop() is non-const.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    const auto cancelled = NastranPbeamCards::prepare(input.materials, input.beams, property_limits,
                                                      plan, stop.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK(plan.allocated_bytes() == 0);
    const auto baseline =
        take(NastranPbeamCards::prepare(input.materials, input.beams, property_limits, plan, {}));
    const auto work = take(baseline.view()).work;
    const auto bytes = plan.allocated_bytes();
    REQUIRE(work > 0);
    REQUIRE(bytes > 0);
    core::ResourcePlan exact{bytes};
    CHECK(NastranPbeamCards::prepare(input.materials, input.beams, {1, 1, work}, exact, {}));
    CHECK(exact.remaining_bytes() == 0);
    core::ResourcePlan short_memory{bytes - 1};
    const auto byte_rejected =
        NastranPbeamCards::prepare(input.materials, input.beams, {1, 1, work}, short_memory, {});
    REQUIRE_FALSE(byte_rejected);
    CHECK(byte_rejected.error().code == ResourceExceeded);
    CHECK(short_memory.allocated_bytes() == 0);
    core::ResourcePlan enough{bytes};
    const auto work_rejected =
        NastranPbeamCards::prepare(input.materials, input.beams, {1, 1, work - 1}, enough, {});
    REQUIRE_FALSE(work_rejected);
    CHECK(work_rejected.error().code == ResourceExceeded);
    CHECK(take(input.materials.view()).input.input.bulk.bytes == source);
}

} // namespace
