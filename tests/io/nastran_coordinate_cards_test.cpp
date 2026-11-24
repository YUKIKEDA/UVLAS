#include "nastran_coordinate_cards.hpp"
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

} // namespace
