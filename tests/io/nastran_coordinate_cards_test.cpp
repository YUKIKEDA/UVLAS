#include "../../src/structure/exact_polynomial_detail.hpp"
#include "nastran_position_enclosures.hpp"
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
using Reason = NastranCoordinateFailure;

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

std::string cord(std::string_view name, std::string_view id, std::string_view rid) {
    return std::string{name} + "," + std::string{id} + "," + std::string{rid} +
           ",0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n";
}

std::string fixed(std::string_view name, std::initializer_list<std::string_view> values,
                  std::string_view tail = {}) {
    std::string s{name};
    s.resize(8, ' ');
    for (const auto value : values) {
        REQUIRE(value.size() <= 16);
        s.append(value);
        s.append(16 - value.size(), ' ');
    }
    s.resize(72, ' ');
    return s + std::string{tail} + "\n";
}

TEST_CASE("Nastran coordinate cards own exact tuples and distinct position/displacement references",
          "[IO4-03][IO4-01][verification]") {
    const auto source = cord("CORD2s", "13", "12") + "GRID,98,12,1.0D+3,,-.0,13,316,0\n" +
                        cord("CORD2C", "12", "11") + cord("CORD2R", "11", "") + "GRID,11\n" +
                        "GRID,99999999,0\n" + "PSHELL,20,10,1.0\n";
    const auto typed = [&] {
        auto deck = lexical_deck(source);
        core::ResourcePlan plan{memory};
        return take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
    }();

    const auto result = take(typed.view());
    CHECK(result.bulk.bytes == source);
    CHECK(result.bulk.location.source == SourceId{99});
    CHECK(result.bulk.digest == source_digest(std::as_bytes(std::span{source})));
    CHECK(result.deferred_cards == 1);
    REQUIRE(result.grids.size() == 3);
    REQUIRE(result.coordinates.size() == 3);
    CHECK(std::ranges::equal(result.dependency_first, std::array<std::size_t, 3>{2, 1, 0}));
    const auto& point = result.grids[0];
    CHECK(point.id == NastranGridId{98});
    CHECK(point.position_system.id == NastranCoordinateId{12});
    CHECK(point.displacement_system.id == NastranCoordinateId{13});
    CHECK(point.position[0].decimal() == "1.0e+3");
    CHECK(point.position[1].decimal() == "0.0");
    CHECK_FALSE(point.position[1].supplied);
    CHECK(point.position[2].decimal() == "-0.0");
    const auto location = take(point.position[0].field);
    CHECK(result.bulk.text(location) == "1.0D+3");
    CHECK(location.position == NastranBulkPosition{source.find("1.0D+3"), 2, 12});
    CHECK(point.permanent_components == 37);
    CHECK(result.bulk.text(take(point.permanent_field)) == "316");
    CHECK(result.bulk.text(take(point.superelement_field)) == "0");
    CHECK(result.grids[1].id == NastranGridId{11});
    CHECK(result.grids[1].position_system.id == NastranCoordinateId{0});
    CHECK(result.grids[1].displacement_system.id == NastranCoordinateId{0});
    CHECK_FALSE(result.grids[1].position_system.field);
    CHECK(result.grids[1].position[0].decimal() == "0.0");
    CHECK(result.grids[2].id == NastranGridId{99999999});
    CHECK(result.bulk.text(take(result.grids[2].position_system.field)) == "0");
    CHECK(result.coordinates[0].kind == NastranCoordinateKind::Spherical);
    CHECK(result.coordinates[0].reference.id == NastranCoordinateId{12});
    CHECK(result.coordinates[1].kind == NastranCoordinateKind::Cylindrical);
    CHECK(result.coordinates[1].reference.id == NastranCoordinateId{11});
    CHECK(result.coordinates[2].kind == NastranCoordinateKind::Rectangular);
    CHECK(result.coordinates[2].reference.id == NastranCoordinateId{0});
    CHECK(take(result.coordinates[2].reference.field).size == 0);
    CHECK(result.coordinates[0].a[0].decimal() == "0.0");
    CHECK(result.coordinates[0].b[2].decimal() == "1.0");
    CHECK(result.coordinates[0].c[0].decimal() == "1.0");
    CHECK(result.bulk.text(result.coordinates[0].name) == "CORD2s");
}

TEST_CASE("A large CORD supports exact local IDs and its required physical continuations",
          "[IO4-03][verification]") {
    const auto source = fixed("CORD2R*", {"1000000000000", "", "1.23456789012D-9", "0.0"}, "*ONE") +
                        fixed("*ONE", {"0.0", "0.0", "0.0", "1.0"}, "*TWO") +
                        fixed("*TWO", {"1.0", "0.0", "0.0"}, "*THREE") + fixed("*THREE", {}) +
                        fixed("GRID*", {"2", "1000000000000", "0.125-9", "0.0"}, "*GRID") +
                        fixed("*GRID", {"0.0", "0", "0", "0"});
    auto deck = lexical_deck(source);
    core::ResourcePlan plan{memory};
    const auto typed = take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
    const auto result = take(typed.view());
    REQUIRE(result.coordinates.size() == 1);
    REQUIRE(result.grids.size() == 1);
    CHECK(result.coordinates[0].id == NastranCoordinateId{1000000000000});
    CHECK(result.grids[0].position_system.id == result.coordinates[0].id);
    CHECK(result.coordinates[0].a[0].decimal() == "1.23456789012e-9");
    CHECK(result.grids[0].position[0].decimal() == "0.125e-9");
    CHECK(result.grids[0].permanent_components == 0);
}

TEST_CASE("Nastran coordinate qualification diagnoses unsupported and inconsistent declarations",
          "[IO4-03][negative][verification]") {
    std::string source;
    auto reason = Reason::InvalidField;
    auto code = InvalidInput;
    auto line = std::size_t{1};
    SECTION("duplicate grid IDs") {
        source = "GRID,1\nGRID,1\n";
        reason = Reason::DuplicateGrid;
        line = 2;
    }
    SECTION("duplicate CORD IDs across chart types") {
        source = cord("CORD2R", "3", "") + cord("CORD2C", "3", "");
        reason = Reason::DuplicateCoordinate;
        line = 2;
    }
    SECTION("missing GRID position reference") {
        source = "GRID,1,9\n";
        reason = Reason::UndefinedCoordinate;
    }
    SECTION("missing GRID displacement reference") {
        source = "GRID,1,0,,,,9\n";
        reason = Reason::UndefinedCoordinate;
    }
    SECTION("missing CORD reference") {
        source = cord("CORD2S", "3", "17");
        reason = Reason::UndefinedCoordinate;
    }
    SECTION("cycles follow typed RID edges") {
        source = cord("CORD2R", "3", "4") + cord("CORD2C", "4", "3");
        reason = Reason::CoordinateCycle;
        line = 2;
    }
    SECTION("self reference is a cycle") {
        source = cord("CORD2R", "3", "3");
        reason = Reason::CoordinateCycle;
    }
    SECTION("GRID is below the normal QRG ID limit") {
        source = fixed("GRID*", {"100000000"}, "*NEXT") + fixed("*NEXT", {});
    }
    SECTION("negative CP is invalid") {
        source = "GRID,1,-1\n";
    }
    SECTION("free integer over eight bytes") {
        source = cord("CORD2C", "123456789", "");
    }
    SECTION("the basic system is implicit") {
        source = cord("CORD2C", "0", "");
    }
    SECTION("missing CORD continuation") {
        source = "CORD2R,1,0,0.0,0.0,0.0,0.0,0.0,1.0\n";
        reason = Reason::MissingField;
    }
    SECTION("required real cannot become an integer") {
        source = "CORD2R,1,0,0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n";
    }
    SECTION("data after the last field is not used as data in the record") {
        source = "GRID,1,0,0.0,0.0,0.0,0,0,0,EXTRA\n";
        reason = Reason::ExtraField;
    }
    SECTION("permanent fields list components one through six") {
        source = "GRID,1,0,0.0,0.0,0.0,0,17\n";
    }
    SECTION("fluid displacement sentinel") {
        source = "GRID,1,0,,,,-1\n";
        code = Unsupported;
        reason = Reason::UnsupportedField;
    }
    SECTION("named superelement") {
        source = "GRID,1,0,0.0,0.0,0.0,0,0,4\n";
        code = Unsupported;
        reason = Reason::UnsupportedField;
    }
    SECTION("GRID defaults are not applied when GRDSET is present") {
        source = "GRID,1\nGRDSET,,17\n";
        code = Unsupported;
        reason = Reason::UnsupportedCard;
        line = 2;
    }
    SECTION("other coordinate definitions cannot be discarded") {
        source = "CORD1R,1,2,3,4\n";
        code = Unsupported;
        reason = Reason::UnsupportedCard;
    }
    auto deck = lexical_deck(source);
    core::ResourcePlan plan{memory};
    std::optional<NastranCoordinateDiagnostic> diagnostic;
    const auto error =
        NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}, &diagnostic);
    REQUIRE_FALSE(error);
    CHECK(error.error().code == code);
    CHECK(error.error().operation == core::Operation::PrepareModel);
    REQUIRE(diagnostic);
    if (diagnostic) {
        CHECK(diagnostic->reason == reason);
        CHECK(diagnostic->position.line == line);
    }
    // NOLINTNEXTLINE(bugprone-use-after-move): failed qualification preserves the lexical owner.
    CHECK(take(deck.view()).bytes == source);
}

TEST_CASE("Nastran dependency depth counts an already visited parent path",
          "[IO4-03][N4-02][negative]") {
    const auto first = cord("CORD2R", "1", "");
    const auto second = cord("CORD2S", "2", "1");
    const auto third = cord("CORD2C", "3", "2");
    const auto forward = first + second + third;
    const auto reverse = third + second + first;
    for (const auto* source : {&forward, &reverse}) {
        auto input = lexical_deck(*source);
        auto caps = limits;
        caps.depth = 3;
        core::ResourcePlan plan{memory};
        const auto object = take(NastranCoordinateCards::prepare(std::move(input), caps, plan, {}));
        const auto result = take(object.view());
        REQUIRE(result.dependency_first.size() == 3);
        for (std::size_t i = 0; i < result.dependency_first.size(); ++i)
            CHECK(result.coordinates[result.dependency_first[i]].id == NastranCoordinateId{i + 1});

        --caps.depth;
        auto short_input = lexical_deck(*source);
        core::ResourcePlan constrained{memory};
        std::optional<NastranCoordinateDiagnostic> diagnostic;
        const auto denied = NastranCoordinateCards::prepare(std::move(short_input), caps,
                                                            constrained, {}, &diagnostic);
        REQUIRE_FALSE(denied);
        CHECK(denied.error().code == ResourceExceeded);
        REQUIRE(diagnostic);
        if (diagnostic)
            CHECK(diagnostic->reason == Reason::Limit);
    }
}

TEST_CASE("Nastran typed allocation and work caps have exact negative boundaries",
          "[IO4-03][N4-02][negative]") {
    const auto source = cord("CORD2C", "9", "") + "GRID,1,9\n";
    core::ResourcePlan baseline{memory};
    auto seed = lexical_deck(source);
    const auto checked =
        take(NastranCoordinateCards::prepare(std::move(seed), limits, baseline, {}));
    auto caps = limits;
    caps.grids = 1;
    caps.coordinate_systems = 1;
    caps.depth = 1;
    auto measured = lexical_deck(source);
    core::ResourcePlan measurement{memory};
    const auto exact_baseline =
        take(NastranCoordinateCards::prepare(std::move(measured), caps, measurement, {}));
    caps.work = take(exact_baseline.view()).work;
    auto amount = measurement.allocated_bytes();
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() requires a mutable source.
    std::stop_source stop;
    bool accepts = false;
    SECTION("exact work and capacity") {
        accepts = true;
    }
    SECTION("too little work") {
        --caps.work;
    }
    SECTION("no room for grids") {
        --caps.grids;
    }
    SECTION("no room for coordinate systems") {
        --caps.coordinate_systems;
    }
    SECTION("one byte too little") {
        --amount;
    }
    SECTION("cancelled") {
        REQUIRE(stop.request_stop());
    }
    auto sample = lexical_deck(source);
    core::ResourcePlan plan{amount};
    const auto value =
        NastranCoordinateCards::prepare(std::move(sample), caps, plan, stop.get_token());
    if (accepts) {
        REQUIRE(value);
        CHECK(plan.allocated_bytes() == measurement.allocated_bytes());
    } else {
        REQUIRE_FALSE(value);
        CHECK(value.error().code == (stop.stop_requested() ? Cancelled : ResourceExceeded));
        if (stop.stop_requested())
            CHECK(plan.allocated_bytes() == 0);
    }
    CHECK(take(checked.view()).grids.size() == 1);
}

using Q = uvlas::structure::detail::Rational;
using PositionReason = NastranPositionFailure;
constexpr std::uint64_t position_budget = std::uint64_t{1} << 28U;

[[nodiscard]] NastranPositionLimits position_limits() noexcept {
    return {{512, 8, 32, 8, 8, 4096, 8192, structure::ExactEstimateAdmission{position_budget},
             structure::ExactTimingMode::Disabled},
            8,
            64000,
            {}};
}

[[nodiscard]] const Q& q(const structure::ExactRational& value) {
    return uvlas::structure::detail::ExactFactory::value(value);
}

void contains(const NastranCartesianInterval& value, const Q& expected) {
    const auto& lo = q(value.lower);
    const auto& hi = q(value.upper);
    CHECK(lo <= expected);
    CHECK(hi >= expected);
    CHECK(hi - lo < Q("1/10000000000000000000000000000000000000000"));
}

void contains_root(const NastranCartesianInterval& value, const Q& offset, const Q& scale,
                   const Q& radicand) {
    const Q lower = (q(value.lower) - offset) / scale;
    const Q upper = (q(value.upper) - offset) / scale;
    CHECK(lower > 0);
    CHECK(lower * lower < radicand);
    CHECK(upper * upper > radicand);
    CHECK(q(value.upper) - q(value.lower) < Q("1/10000000000000000000000000000000000000000"));
}

constexpr std::string_view coordinate_fixture =
    "CORD2R,10,0,0.0,0.0,0.0,1.0,0.0,1.0,1.0,0.0,-1.0\n"
    "CORD2S,30,20,0.0,723.0,200.0,0.0,-120.0,201.0,1.0,0.0,200.0\n"
    "CORD2R,40,30,0.0,45.0,-999.0,10.0,0.0,120.0,10.0,90.0,0.0\n"
    "CORD2C,20,0,100.0,-200.0,50.0,100.0,-200.0,51.0,101.0,-200.0,50.0\n"
    "GRID,1,0,-4.0,5.0,6.0,30\n"
    "GRID,2,10,2.0,0.3,4.0\n"
    "GRID,3,20,3.0,-360000000000000000000000000000030.0,5.0\n"
    "GRID,4,30,10.0,60.0,30.0\n"
    "GRID,5,40,1.0,2.0,3.0\n"
    "GRID,6,20,2.0,90.0000000000000000000000000001,5.0\n";

TEST_CASE("The Nastran CP position envelope owns rectangular cylindrical and spherical results",
          "[IO4-03][IO4-01][N4-02][integration][verification]") {
    constexpr std::array selection{NastranGridId{5}, NastranGridId{3}, NastranGridId{4},
                                   NastranGridId{2}, NastranGridId{1}, NastranGridId{6}};
    const auto output = [&] {
        auto deck = lexical_deck(coordinate_fixture);
        core::ResourcePlan plan{position_budget};
        auto typed = take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
        return take(NastranPositionEnclosures::prepare(typed, SiLengthUnit::Centimetre, selection,
                                                       position_limits(), plan));
    }();

    const auto read = take(output.view());
    REQUIRE(read.grids.size() == selection.size());
    CHECK(read.work > 0);
    CHECK(read.source.source == SourceId{99});
    CHECK(read.digest == source_digest(std::as_bytes(std::span{coordinate_fixture})));
    for (std::size_t i = 0; i < selection.size(); ++i) {
        CHECK(read.grids[i].grid == selection[i]);
        CHECK(read.grids[i].cp_field.has_value());
        CHECK(read.grids[i].tuple_fields[0].has_value());
    }

    // The mixed R->S->C path retains Cartesian coordinates when moving to the parent.
    contains(read.grids[0].metres[0], Q("101/100"));
    contains(read.grids[0].metres[1], Q("-99/50"));
    contains(read.grids[0].metres[2], Q("253/100"));
    CHECK(read.grids[0].cp == NastranCoordinateId{40});
    // -30 degrees represented after an exact large negative number of turns.
    contains_root(read.grids[1].metres[0], Q{1}, Q("3/200"), Q{3});
    contains(read.grids[1].metres[1], Q("-403/200"));
    contains(read.grids[1].metres[2], Q("11/20"));
    // Spherical polar angle 60, azimuth 30 and an origin defined in cylindrical RID.
    contains(read.grids[2].metres[0], Q("43/40"));
    contains_root(read.grids[2].metres[1], Q{-2}, Q("1/40"), Q{3});
    contains(read.grids[2].metres[2], Q("51/20"));
    // Rectangular normalization stretches the bounds to cover two irrational components.
    contains_root(read.grids[3].metres[0], Q{0}, Q("3/100"), Q{2});
    contains(read.grids[3].metres[1], Q("3/1000"));
    contains_root(read.grids[3].metres[2], Q{0}, Q("1/100"), Q{2});
    // CD on the basic GRID gives its vector basis and cannot move the position.
    contains(read.grids[4].metres[0], Q("-1/25"));
    contains(read.grids[4].metres[1], Q("1/20"));
    contains(read.grids[4].metres[2], Q("3/50"));
    // Just after 90 degrees the x coordinate must be strictly below the frame origin.
    CHECK(q(read.grids[5].metres[0].upper) < 1);
    CHECK(q(read.grids[5].metres[0].lower) > Q("9999999999999999999999/10000000000000000000000"));
}

TEST_CASE("Nastran position conversion distinguishes degenerate frames and admission failures",
          "[IO4-03][IO4-01][N4-02][integration][negative][verification]") {
    constexpr std::string_view fixture = "CORD2R,10,0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,2.0\n"
                                         "GRID,1,10,1.0,2.0,3.0\nGRID,2,0,4.0,5.0,6.0\n";
    auto deck = lexical_deck(fixture);
    core::ResourcePlan plan{position_budget};
    auto coordinates = take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
    const auto input = take(coordinates.view());
    std::optional<NastranPositionDiagnostic> info;
    const auto read = [&](std::span<const NastranGridId> ids, const NastranPositionLimits& bounds,
                          core::ResourcePlan& payload) {
        return NastranPositionEnclosures::prepare(coordinates, SiLengthUnit::Metre, ids, bounds,
                                                  payload, &info);
    };
    const std::array degenerate{NastranGridId{1}};
    const std::array good{NastranGridId{2}};
    const auto invalid = read(degenerate, position_limits(), plan);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == GeometryInvalid);
    const auto diagnostic = take(info);
    CHECK(diagnostic.reason == PositionReason::FramePlane);
    CHECK(diagnostic.grid == NastranGridId{1});
    CHECK(diagnostic.coordinate == NastranCoordinateId{10});
    CHECK(diagnostic.position == take(input.coordinates[0].c[0].field).position);

    auto capped = position_limits();
    capped.maximum_work = 0;
    core::ResourcePlan no_work{position_budget};
    const auto resource = read(good, capped, no_work);
    REQUIRE_FALSE(resource);
    CHECK(resource.error().code == ResourceExceeded);
    CHECK(no_work.allocated_bytes() == 0);
    const auto empty = take(read({}, capped, no_work));
    CHECK(take(empty.view()).grids.empty());
    CHECK_FALSE(info);
    const std::array missing{NastranGridId{91}};
    const auto absent = read(missing, position_limits(), plan);
    REQUIRE_FALSE(absent);
    CHECK(absent.error().code == InvalidInput);
    CHECK(take(info).reason == PositionReason::Selection);
    CHECK(take(info).grid == NastranGridId{91});
    CHECK_FALSE(take(info).position);
    const auto later = take(read(good, position_limits(), plan));
    CHECK_FALSE(info);
    CHECK(q(take(later.view()).grids.front().metres[0].lower) == 4);
    CHECK(q(take(later.view()).grids.front().metres[0].upper) == 4);
    CHECK(take(coordinates.view()).bulk.bytes == fixture);

    const std::string ambiguous_input = cord("CORD2C", "11", "0") +
                                        "CORD2R,12,11,1.0,30.0,0.0,1.0,390.0,0.0,2.0,0.0,0.0\n"
                                        "GRID,1,12,1.0,2.0,3.0\n";
    auto raw = lexical_deck(ambiguous_input);
    const auto prepared = take(NastranCoordinateCards::prepare(std::move(raw), limits, plan, {}));
    // The equal angles give identical mathematical points, but subtracting their
    // non-singleton enclosures cannot prove equality of the defining points.
    const auto unresolved = NastranPositionEnclosures::prepare(
        prepared, SiLengthUnit::Metre, degenerate, position_limits(), plan, &info);
    REQUIRE_FALSE(unresolved);
    CHECK(unresolved.error().code == PrecisionUnresolved);
    CHECK(take(info).reason == PositionReason::FrameAxis);
    CHECK(take(info).coordinate == NastranCoordinateId{12});

    const auto shells =
        take(NastranShellCards::prepare(std::move(coordinates), {8, 8, 10000}, plan, {}));
    const auto bridged = take(NastranPositionEnclosures::prepare(shells, SiLengthUnit::Metre, good,
                                                                 position_limits(), plan, &info));
    CHECK_FALSE(info);
    const auto shell_view = take(bridged.view());
    CHECK(shell_view.digest == take(later.view()).digest);
    CHECK(q(shell_view.grids.front().metres[0].lower) == 4);
}

constexpr std::string_view spc_fixture =
    "CORD2R,10,0,100.0,200.0,300.0,101.0,200.0,301.0,101.0,200.0,299.0\n"
    "CORD2C,20,10,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n"
    "CORD2S,30,0,11.0,21.0,31.0,11.0,21.0,32.0,12.0,21.0,31.0\n"
    "GRID,1,,0.0,0.0,0.0,10\n"
    "GRID,2,10,0.0,2.0,3.0,20\n"
    "GRID,3,,12.0,23.0,33.0,30\n"
    "GRID,4,,11.0,21.0,35.0,30\n"
    "GRID,6,20,1.0,30.0,1.0,20\n"
    "GRID,7,10,0.0,0.0,2.0,20\n"
    "SPC,10,1,14,0.125\n"
    "SPC1,20,123,2,THRU,7\n"
    "SPC,30,3,123456,0.0\n"
    "SPC,40,4,2,0.0\n"
    "GRID,8,20,0.0,30.0,2.0,20\n"
    "SPC1,50,36,7,THRU,8\n";

[[nodiscard]] NastranSpcDirectionLimits spc_limits() noexcept {
    const auto base = position_limits();
    return {base.exact, base.maximum_work, {}};
}

void contains_signed_root(const NastranCartesianInterval& entry, int sign, const Q& squared) {
    const Q lower = sign > 0 ? q(entry.lower) : -q(entry.upper);
    const Q upper = sign > 0 ? q(entry.upper) : -q(entry.lower);
    CHECK(lower > 0);
    CHECK(lower * lower < squared);
    CHECK(upper * upper > squared);
    CHECK(upper - lower < Q("1/10000000000000000000000000000000000000000"));
}

TEST_CASE("SPC components use the matched GRID CD basis and retain their deck fields",
          "[IO4-03][IO4-01][N4-02][integration][verification]") {
    const auto owners = [&] {
        auto raw = lexical_deck(spc_fixture);
        core::ResourcePlan plan{position_budget};
        const auto input = take(NastranCoordinateCards::prepare(std::move(raw), limits, plan, {}));
        const auto supports = take(NastranSupportCards::prepare(input, {8, 8, 10000}, plan, {}));
        const auto read = [&](std::size_t term, std::uint64_t grid) {
            return take(NastranSpcDirections::prepare(input, supports, term, NastranGridId{grid},
                                                      spc_limits(), plan));
        };
        return std::array{read(0, 1), read(1, 2), read(2, 3), read(1, 6)};
    }();

    using Kind = NastranSpcComponentKind;
    const auto& rect = owners[0];
    const auto* source = rect.origin();
    REQUIRE(source);
    CHECK(source->set == NastranSupportSetId{10});
    CHECK(source->source.location.source == SourceId{99});
    CHECK(source->source.digest == source_digest(std::as_bytes(std::span{spc_fixture})));
    CHECK(source->cd == NastranCoordinateId{10});
    CHECK(take(source->cd_field).line == 4);
    CHECK(source->component_field.line == 10);
    CHECK(take(source->motion_field).line == 10);
    if (const auto& value = source->prescribed_motion) {
        CHECK(value->decimal() == "0.125");
    } else {
        FAIL("Missing source SPC value");
    }
    CHECK_FALSE(rect.direction(0));
    CHECK_FALSE(rect.direction(2));
    CHECK_FALSE(rect.direction(7));
    CHECK(take(rect.direction(1)).kind == Kind::Translation);
    CHECK(take(rect.direction(4)).kind == Kind::Rotation);
    const auto& local_x = take(rect.direction(1)).basic;
    CHECK(&local_x == &take(rect.direction(4)).basic);
    contains_signed_root(local_x[0], 1, Q("1/2"));
    contains(local_x[1], Q{0});
    contains_signed_root(local_x[2], -1, Q("1/2"));

    const auto& cylinder = owners[1];
    // Polar azimuth at GRID 2 is +Y in the translated, rotated CD20 frame.
    const auto& radial = take(cylinder.direction(1)).basic;
    const auto& tangential = take(cylinder.direction(2)).basic;
    const auto& axial = take(cylinder.direction(3)).basic;
    contains(radial[0], Q{0});
    contains(radial[1], Q{1});
    contains(radial[2], Q{0});
    contains_signed_root(tangential[0], -1, Q("1/2"));
    contains(tangential[1], Q{0});
    contains_signed_root(tangential[2], 1, Q("1/2"));
    contains_signed_root(axial[0], 1, Q("1/2"));
    contains(axial[1], Q{0});
    contains_signed_root(axial[2], 1, Q("1/2"));

    // In the spherical CD30 frame the relative position is (1,2,2).
    // The three directions are (1,2,2)/3, (2,4,-5)/(3sqrt(5)), (-2,1,0)/sqrt(5).
    const auto& sphere = owners[2];
    const auto& sphere_r = take(sphere.direction(1)).basic;
    const auto& sphere_theta = take(sphere.direction(2)).basic;
    const auto& sphere_phi = take(sphere.direction(6)).basic;
    contains(sphere_r[0], Q("1/3"));
    contains(sphere_r[1], Q("2/3"));
    contains(sphere_r[2], Q("2/3"));
    contains_signed_root(sphere_theta[0], 1, Q("4/45"));
    contains_signed_root(sphere_theta[1], 1, Q("16/45"));
    contains_signed_root(sphere_theta[2], -1, Q("5/9"));
    contains_signed_root(sphere_phi[0], -1, Q("4/5"));
    contains_signed_root(sphere_phi[1], 1, Q("1/5"));
    contains(sphere_phi[2], Q{0});

    const auto& shared_cp = take(owners[3].direction(1)).basic;
    contains_signed_root(shared_cp[0], 1, Q("3/8"));
    contains(shared_cp[1], Q("1/2"));
    contains_signed_root(shared_cp[2], -1, Q("3/8"));
}

TEST_CASE("SPC direction preparation diagnoses polar axes and stale or absent source selections",
          "[IO4-03][IO4-01][N4-02][integration][negative][verification]") {
    core::ResourcePlan plan{position_budget};
    auto raw = lexical_deck(spc_fixture);
    const auto input = take(NastranCoordinateCards::prepare(std::move(raw), limits, plan, {}));
    const auto support = take(NastranSupportCards::prepare(input, {8, 8, 10000}, plan, {}));
    const auto deck = take(input.view());
    std::optional<NastranPositionDiagnostic> diagnostic;
    const auto read = [&](std::size_t term, std::uint64_t grid) {
        return NastranSpcDirections::prepare(input, support, term, NastranGridId{grid},
                                             spc_limits(), plan, &diagnostic);
    };
    const auto pole = read(3, 4);
    REQUIRE_FALSE(pole);
    CHECK(pole.error().code == Unsupported);
    CHECK(take(diagnostic).reason == PositionReason::DisplacementAxis);
    CHECK(take(diagnostic).position == take(deck.grids[3].displacement_system.field).position);
    CHECK(take(diagnostic).coordinate == NastranCoordinateId{30});

    const auto unresolved = read(1, 7);
    REQUIRE_FALSE(unresolved);
    CHECK(unresolved.error().code == PrecisionUnresolved);
    CHECK(take(diagnostic).reason == PositionReason::DisplacementAxis);
    CHECK(take(diagnostic).grid == NastranGridId{7});
    CHECK(take(diagnostic).position == take(deck.grids[5].displacement_system.field).position);

    // The cylindrical axial direction is unique even at the polar axis.
    for (const auto grid : {7U, 8U}) {
        const auto axial = take(read(4, grid));
        CHECK_FALSE(axial.direction(1));
        const auto& direction = take(axial.direction(3)).basic;
        CHECK(&direction == &take(axial.direction(6)).basic);
        contains_signed_root(direction[0], 1, Q("1/2"));
        contains(direction[1], Q{0});
        contains_signed_root(direction[2], 1, Q("1/2"));
    }

    const auto missing = read(1, 5);
    REQUIRE_FALSE(missing);
    CHECK(missing.error().code == InvalidInput);
    CHECK(take(diagnostic).reason == PositionReason::Selection);
    const auto outside = read(0, 2);
    REQUIRE_FALSE(outside);
    CHECK(outside.error().code == InvalidInput);
    CHECK(take(diagnostic).reason == PositionReason::SpcTerm);
    CHECK(take(diagnostic).position == take(support.view()).terms[0].component_field.position);

    auto changed_raw = lexical_deck(std::string{spc_fixture} + "$ differing source bytes\n");
    const auto changed =
        take(NastranCoordinateCards::prepare(std::move(changed_raw), limits, plan, {}));
    const auto stale = NastranSpcDirections::prepare(changed, support, 0, NastranGridId{1},
                                                     spc_limits(), plan, &diagnostic);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == StaleGeneration);
    CHECK(take(diagnostic).reason == PositionReason::Source);

    auto first = take(read(0, 1));
    auto exact_work = spc_limits();
    exact_work.maximum_work = first.work() - 1;
    auto limited =
        NastranSpcDirections::prepare(input, support, 0, NastranGridId{1}, exact_work, plan);
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == ResourceExceeded);
    ++exact_work.maximum_work;
    REQUIRE(NastranSpcDirections::prepare(input, support, 0, NastranGridId{1}, exact_work, plan));
    const auto moved = std::move(first);
    CHECK(moved.origin());
    // The owning type explicitly permits observing an empty moved-from origin.
    // NOLINTNEXTLINE(bugprone-use-after-move)
    CHECK_FALSE(first.origin());
    CHECK(take(input.view()).bulk.bytes == spc_fixture);
}

} // namespace
