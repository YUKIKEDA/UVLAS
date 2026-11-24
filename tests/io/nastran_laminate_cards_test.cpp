#include "nastran_laminate_cards.hpp"
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
using Reason = NastranLaminateFailure;

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

constexpr NastranLaminateLimits laminate_limits{24, 25000};
constexpr std::string_view basic_mat = "MAT8,81,130.0,9.0,.28\n";
constexpr std::string_view simple_pcomp = "PCOMP,30,,,,,,,\n,81,.01\n";

NastranIsotropicShellCards material_cards(std::string_view source) {
    auto deck = lexical_deck(source);
    core::ResourcePlan plan{memory};
    auto coords = take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
    auto shell = take(NastranShellCards::prepare(std::move(coords), {16, 16, 20000}, plan, {}));
    return take(NastranIsotropicShellCards::prepare(std::move(shell), {16, 20000}, plan, {}));
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

TEST_CASE("PCOMP keeps ordered declared plies and exact MAT8 fields from the shell snapshot",
          "[IO4-03][IO4-01][verification]") {
    const std::string source = "GRID,1\nGRID,2\nGRID,3\nCTRIA3,10,30,1,2,3,45.0\n"
                               "PCOMP,30,-.0025,.14,2.5D+6,STRN,77.0,.04,SYM\n"
                               ",81,.001,45.0,YES,,,0.0,NO\n"
                               ",52,.002,-30.0,,81,,15.0,YES\n"
                               "PCOMP,31,,,,,,,\n,52,.05\n"
                               "MAT1,52,1.0\n"
                               "MAT8,81,13.0D+10,9.0D+9,.28,6.0+9,-0.0,,1.56+3\n"
                               ",4.0-7,3.5-5,110.0,2.1+9,1.1+9,3.3+7,8.85+7,6.6+7\n"
                               ",.02,-.001,.1D+1\n";
    core::ResourcePlan plan{memory};
    const char* original = nullptr;
    const auto cards = [&] {
        auto input = material_cards(source);
        original = take(input.view()).input.input.bulk.bytes.data();
        auto typed =
            take(NastranLaminateCards::prepare(std::move(input), laminate_limits, plan, {}));
        // NOLINTNEXTLINE(bugprone-use-after-move): success leaves the input view absent.
        CHECK_FALSE(input.view());
        return typed;
    }();
    const auto view = take(cards.view());
    const auto& bulk = view.input.input.input.bulk;
    CHECK(bulk.bytes.data() == original);
    CHECK(bulk.bytes == source);
    CHECK(bulk.digest == source_digest(std::as_bytes(std::span{source})));
    REQUIRE(view.mat8.size() == 1);
    REQUIRE(view.pcomp.size() == 2);
    REQUIRE(view.plies.size() == 5);
    REQUIRE(view.input.input.elements.size() == 1);
    const auto& lam = view.pcomp[0];
    CHECK(lam.property_ordinal == view.input.input.elements[0].property.property_ordinal);
    CHECK(view.input.input.properties[lam.property_ordinal].id == NastranPropertyId{30});
    CHECK(lam.mode == NastranLaminateMode::SymmetricHalf);
    CHECK(bulk.text(take(lam.mode_field)) == "SYM");
    const auto theory = take(lam.failure_theory);
    CHECK(theory.theory == NastranPcompFailureTheory::MaximumStrain);
    CHECK(bulk.text(theory.field) == "STRN");
    const auto real = [&](const std::optional<NastranShellReal>& source_field, std::string_view raw,
                          std::string_view decimal) {
        const auto value = take(source_field);
        CHECK(value.value.decimal() == decimal);
        CHECK(bulk.text(value.field) == raw);
    };
    real(lam.bottom_offset, "-.0025", "-0.0025");
    real(lam.nonstructural_mass, ".14", "0.14");
    real(lam.bond_shear_limit, "2.5D+6", "2.5e+6");
    real(lam.reference_temperature, "77.0", "77.0");
    real(lam.damping, ".04", "0.04");
    const auto& mat = view.mat8[0];
    const auto material = view.input.materials[mat.declaration_ordinal];
    CHECK(material.id == NastranMaterialId{81});
    CHECK(material.kind == NastranMaterialKind::Mat8);
    real(mat.longitudinal_young, "13.0D+10", "13.0e+10");
    real(mat.transverse_young, "9.0D+9", "9.0e+9");
    real(mat.poisson_12, ".28", "0.28");
    real(mat.shear_12, "6.0+9", "6.0e+9");
    real(mat.shear_1z, "-0.0", "-0.0");
    CHECK_FALSE(mat.shear_2z);
    real(mat.density, "1.56+3", "1.56e+3");
    real(mat.expansion_1, "4.0-7", "4.0e-7");
    real(mat.expansion_2, "3.5-5", "3.5e-5");
    real(mat.reference_temperature, "110.0", "110.0");
    real(mat.longitudinal_tension, "2.1+9", "2.1e+9");
    real(mat.longitudinal_compression, "1.1+9", "1.1e+9");
    real(mat.transverse_tension, "3.3+7", "3.3e+7");
    real(mat.transverse_compression, "8.85+7", "8.85e+7");
    real(mat.shear_limit, "6.6+7", "6.6e+7");
    real(mat.damping, ".02", "0.02");
    real(mat.tsai_wu_interaction, "-.001", "-0.001");
    real(mat.strain_allowables, ".1D+1", "0.1e+1");

    const auto plies = view.layers(lam);
    REQUIRE(plies.size() == 4);
    CHECK(take(plies[0].material).declaration_ordinal == mat.declaration_ordinal);
    CHECK(plies[0].material_source_ply == 0);
    CHECK(plies[0].thickness_source_ply == 0);
    CHECK(take(plies[0].stress_output).requested);
    real(plies[0].thickness, ".001", "0.001");
    real(plies[0].angle, "45.0", "45.0");
    CHECK_FALSE(plies[1].material);
    CHECK_FALSE(plies[1].thickness);
    CHECK(plies[1].material_source_ply == 0);
    CHECK(plies[1].thickness_source_ply == 0);
    CHECK(plies[1].source.position == NastranBulkPosition{source.find(",,,0.0") + 3, 6, 21});
    CHECK_FALSE(take(plies[1].stress_output).requested);
    CHECK(take(plies[2].material).id == NastranMaterialId{52});
    CHECK(plies[2].thickness_source_ply == 2);
    CHECK_FALSE(plies[2].stress_output);
    real(plies[2].angle, "-30.0", "-30.0");
    CHECK(take(plies[3].material).id == NastranMaterialId{81});
    CHECK(plies[3].material_source_ply == 3);
    CHECK(plies[3].thickness_source_ply == 2);
    CHECK_FALSE(plies[3].thickness);
    const auto& full = view.pcomp[1];
    CHECK(full.mode == NastranLaminateMode::FullStack);
    CHECK_FALSE(full.mode_field);
    CHECK_FALSE(full.bottom_offset);
    const auto single = view.layers(full);
    REQUIRE(single.size() == 1);
    CHECK(single[0].material_source_ply == 0);
    CHECK(single[0].thickness_source_ply == 0);
    CHECK_FALSE(single[0].angle);
}

TEST_CASE("PCOMP preserves a one-ply fixed row, skips its blank group and reads LAM distinctions",
          "[IO4-03][verification]") {
    using Lam = NastranLaminateMode;
    using Ft = NastranPcompFailureTheory;
    struct Option {
        std::string_view lam, ft;
        Lam mode;
        Ft theory;
    };
    constexpr std::array options{Option{"SYM", "HILL", Lam::SymmetricHalf, Ft::Hill},
                                 Option{"MEM", "HOFF", Lam::MembraneOnly, Ft::Hoffman},
                                 Option{"BEND", "TSAI", Lam::BendingOnly, Ft::TsaiWu},
                                 Option{"SMEAR", "HFAIL", Lam::Smeared, Ft::Hashin},
                                 Option{"SMCORE", "HTAPE", Lam::SmearedWithCore, Ft::HashinTape},
                                 Option{"", "HFABR", Lam::FullStack, Ft::HashinFabric}};
    for (const auto& option : options) {
        CAPTURE(option.lam);
        const auto source =
            std::string{basic_mat} +
            fixed("PCOMP", {"30", "", "", "2.0", option.ft, "", "", option.lam}, 8, "+PLY1") +
            fixed("+PLY1", {"81", ".04", "10.0", "YES"}, 8, "+PLY2") +
            fixed("+PLY2", {"", "", "-10.0", "NO"});
        auto input = material_cards(source);
        core::ResourcePlan plan{memory};
        const auto typed =
            take(NastranLaminateCards::prepare(std::move(input), laminate_limits, plan, {}));
        const auto view = take(typed.view());
        REQUIRE(view.pcomp.size() == 1);
        const auto& row = view.pcomp[0];
        CHECK(row.mode == option.mode);
        CHECK(take(row.failure_theory).theory == option.theory);
        const auto layers = view.layers(row);
        REQUIRE(layers.size() == 2);
        CHECK(layers[1].material_source_ply == 0);
        CHECK(layers[1].thickness_source_ply == 0);
        const auto angle = take(layers[1].angle);
        CHECK(angle.value.decimal() == "-10.0");
        CHECK(angle.field.position == NastranBulkPosition{source.find("-10.0"), 4, 25});
        CHECK_FALSE(view.mat8[0].shear_12);
        CHECK_FALSE(view.mat8[0].density);
    }
}

TEST_CASE("Large-field MAT8 and PCOMP keep exact decimals and four fields per ply",
          "[IO4-03][verification]") {
    const auto source = fixed("MAT8*", {"81", "1.23456789012D+7", "9.0", ".64"}, 16, "*M1") +
                        fixed("*M1", {"3.0", "4.0", "5.0", "1.0"}, 16) +
                        fixed("PCOMP*", {"30", "-.02", "", ""}, 16, "*P1") +
                        fixed("*P1", {"", "", "", "SYM"}, 16, "*P2") +
                        fixed("*P2", {"81", ".04", "30.0", "YES"}, 16, "*P3") +
                        fixed("*P3", {"", "", "-30.0", "NO"}, 16);
    auto input = material_cards(source);
    core::ResourcePlan plan{memory};
    const auto typed =
        take(NastranLaminateCards::prepare(std::move(input), laminate_limits, plan, {}));
    const auto view = take(typed.view());
    REQUIRE(view.mat8.size() == 1);
    REQUIRE(view.pcomp.size() == 1);
    const auto young = take(view.mat8[0].longitudinal_young);
    CHECK(young.value.decimal() == "1.23456789012e+7");
    CHECK(young.field.position == NastranBulkPosition{source.find("1.23456789012D+7"), 1, 25});
    const auto poisson = take(view.mat8[0].poisson_12);
    CHECK(poisson.value.decimal() == "0.64");
    const auto layers = view.layers(view.pcomp[0]);
    REQUIRE(layers.size() == 2);
    CHECK(take(layers[0].material).id == NastranMaterialId{81});
    CHECK(layers[1].material_source_ply == 0);
    CHECK(layers[1].thickness_source_ply == 0);
    CHECK_FALSE(layers[1].thickness);
    const auto angle = take(layers[1].angle);
    CHECK(angle.value.decimal() == "-30.0");
    CHECK(angle.field.position == NastranBulkPosition{source.find("-30.0"), 6, 41});
}

TEST_CASE(
    "MAT8 and PCOMP reject unsupported options and invalid ply references at the source field",
    "[IO4-03][negative][verification]") {
    struct FailureCase {
        std::string_view entry;
        Reason reason;
        core::FailureCode code;
        std::size_t line;
        std::size_t column;
    };
    constexpr std::array cases{
        FailureCase{"MAT8,82,,2.0,.3\n", Reason::MissingField, InvalidInput, 2, 9},
        FailureCase{"MAT8,82,1.0,2.0\n", Reason::MissingField, InvalidInput, 2, 1},
        FailureCase{"MAT8,82,-1.0,2.0,.3\n", Reason::InvalidField, InvalidInput, 2, 9},
        FailureCase{"MAT8,82,1.0,0.0,.3\n", Reason::InvalidField, InvalidInput, 2, 13},
        FailureCase{"MAT8,82,1.0,2.0,.3,0.0\n", Reason::InvalidField, InvalidInput, 2, 20},
        FailureCase{"MAT8,82,1.0,2.0,.3,1.0,-.01\n", Reason::InvalidField, InvalidInput, 2, 24},
        FailureCase{"MAT8,82,1.0,2.0,.3,,,,\n,0.0,,,,,,,\n,,,0.99999999999\n", Reason::InvalidField,
                    InvalidInput, 4, 4},
        FailureCase{"MAT8,82,1.0,2.0,.3,,,,\n,0.0,,,,,,,\n,,,1.0,,,,,\n,HFAIL,2.0\n",
                    Reason::UnsupportedField, Unsupported, 5, 2},
        FailureCase{"PCOMP,30\n", Reason::MissingField, InvalidInput, 2, 1},
        FailureCase{"PCOMP,30,,,,,,,\n,,.01,10.0\n", Reason::MissingField, InvalidInput, 3, 2},
        FailureCase{"PCOMP,30,,,,,,,\n,81,,10.0\n", Reason::MissingField, InvalidInput, 3, 5},
        FailureCase{"PCOMP,30,,,,,,,\n,82,.01\n", Reason::UndefinedMaterial, InvalidInput, 3, 2},
        FailureCase{"PCOMP,30,,,,,,,\n,81,.01,10.0,MAYBE\n", Reason::UnsupportedField, Unsupported,
                    3, 14},
        FailureCase{"PCOMP,30,,,,,,,SYMMETRIC\n,81,.01\n", Reason::UnsupportedField, Unsupported, 2,
                    16},
        FailureCase{"PCOMP,30,,,2.0,FAIL,,,\n,81,.01\n", Reason::UnsupportedField, Unsupported, 2,
                    16},
        FailureCase{"PCOMP,30,,,,HILL,,,\n,81,.01\n", Reason::MissingField, InvalidInput, 2, 12},
        FailureCase{"PCOMP,30,,,0.0,HILL,,,\n,81,.01\n", Reason::InvalidField, InvalidInput, 2, 12},
        FailureCase{"PCOMP,30,,,,,,,\n,81,.01,10.0,,000000081,.02\n", Reason::InvalidField,
                    InvalidInput, 3, 15}};
    for (const auto& entry : cases) {
        CAPTURE(entry.entry);
        const auto source = std::string{basic_mat} + std::string{entry.entry};
        auto input = material_cards(source);
        core::ResourcePlan plan{memory};
        std::optional<NastranLaminateDiagnostic> diagnostic;
        const auto result =
            NastranLaminateCards::prepare(std::move(input), laminate_limits, plan, {}, &diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == entry.code);
        const auto actual = take(diagnostic);
        CHECK(actual.reason == entry.reason);
        CHECK(actual.position.line == entry.line);
        CHECK(actual.position.column == entry.column);
        // NOLINTNEXTLINE(bugprone-use-after-move): rejection preserves the input owner.
        CHECK(take(input.view()).input.input.bulk.bytes == source);
    }
}

TEST_CASE("Nastran laminate counts declared plies and bounds memory, work and cancellation",
          "[IO4-03][N4-02][negative][verification]") {
    const auto source = std::string{basic_mat} + std::string{simple_pcomp};
    const auto run = [&](NastranLaminateLimits capacity, core::ResourcePlan& plan,
                         const std::stop_token& stop = {}) {
        auto input = material_cards(source);
        auto result = NastranLaminateCards::prepare(std::move(input), capacity, plan, stop);
        if (!result) {
            // NOLINTNEXTLINE(bugprone-use-after-move): rejected calls preserve input.
            CHECK(take(input.view()).input.input.bulk.bytes == source);
        }
        return result;
    };
    core::ResourcePlan enough{memory};
    const auto too_many_plies = run({0, laminate_limits.work}, enough);
    REQUIRE_FALSE(too_many_plies);
    CHECK(too_many_plies.error().code == ResourceExceeded);
    CHECK(enough.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness) MSVC request_stop() is non-const.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    const auto cancelled = run(laminate_limits, enough, cancel.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK(enough.allocated_bytes() == 0);
    const auto prepared = run({1, laminate_limits.work}, enough);
    REQUIRE(prepared);
    const auto work = take(prepared->view()).work;
    const auto bytes = enough.allocated_bytes();
    REQUIRE(work > 0);
    REQUIRE(bytes > 0);
    core::ResourcePlan exact{bytes};
    const auto exact_fit = run({1, work}, exact);
    CHECK(exact_fit.has_value());
    CHECK(exact.remaining_bytes() == 0);
    core::ResourcePlan short_memory{bytes - 1};
    const auto byte_rejected = run({1, work}, short_memory);
    REQUIRE_FALSE(byte_rejected);
    CHECK(byte_rejected.error().code == ResourceExceeded);
    CHECK(short_memory.allocated_bytes() == 0);
    core::ResourcePlan full_memory{bytes};
    const auto work_rejected = run({1, work - 1}, full_memory);
    REQUIRE_FALSE(work_rejected);
    CHECK(work_rejected.error().code == ResourceExceeded);
}

} // namespace
