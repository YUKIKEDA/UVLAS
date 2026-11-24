#include "../structure/pg_partition_test_support.hpp"
#include "case_schema_test_adapter.hpp"
#include "external_mesh_quantity.hpp"
#include "nastran_laminate_cards.hpp"
#include "nastran_pcomp_physical.hpp"
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

namespace a = uvlas::test::schema_adapter;
constexpr std::uint64_t reconciliation_budget = 256ULL * 1024 * 1024;
constexpr SiUnits si{SiLengthUnit::Metre, SiMassUnit::Kilogram, SiTimeUnit::Second,
                     SiAngleUnit::Radian};
constexpr SiUnits deck_units{SiLengthUnit::Millimetre, SiMassUnit::Tonne, SiTimeUnit::Second,
                             SiAngleUnit::Degree};

structure::ExactLimits comparison_limits() noexcept {
    return {256,
            0,
            0,
            0,
            0,
            4096,
            8192,
            structure::ExactEstimateAdmission{reconciliation_budget},
            structure::ExactTimingMode::Disabled};
}

SchemaCase checked_case(const a::JsonToken& input) {
    const auto wire = input.wire();
    const InputLimits syntax{
        1000000, 256,
        100000,  100000,
        1000000, 1000000,
        1000000, {ParserAllocationMode::EstimateAdmission, reconciliation_budget}};
    core::ResourcePlan plan{reconciliation_budget};
    std::optional<InputDiagnostic> parsing;
    const auto doc = take(
        read_input_document(std::as_bytes(std::span{wire}), SourceId{83}, syntax, plan, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(SchemaQuota::Count)> caps{};
    caps.fill(reconciliation_budget);
    SchemaDiagnostic diagnostic;
    auto result = validate_case_schema(doc, SchemaLimits{caps}, plan, &diagnostic);
    INFO("Schema fixture field " << diagnostic.path().value_or(""));
    return take(std::move(result));
}

TEST_CASE("External explicit material fields must agree in exact SI with the case",
          "[IO4-03][IO4-04][IO4-01][N4-02][verification][negative]") {
    auto input = a::root();
    auto& materials = input.field("model").field("materials").elements;
    auto& isotropic = materials.at(0);
    isotropic.set("young", a::parse("200000000000"));
    isotropic.set("poisson", a::parse("0.25"));
    isotropic.set("density", a::parse("7850"));
    auto orthotropic = a::fragment("Material.orthotropic3d");
    orthotropic.set("id", a::string("67"));
    orthotropic.set("young", a::parse("[130000000000,9000000000,11000000000]"));
    orthotropic.set("shear", a::parse("[4000000000,5000000000,6000000000]"));
    orthotropic.set("poisson", a::parse("[0.28,0.15,0.11]"));
    orthotropic.set("density", a::parse("1560"));
    materials.push_back(std::move(orthotropic));
    const auto schema = checked_case(input);
    const auto case_model = take(schema.root()).model.record();
    const auto case_mat1 = take(case_model.materials.at(0)).record();
    const auto case_mat8 = take(case_model.materials.at(1)).record();

    core::ResourcePlan plan{reconciliation_budget};
    const auto check = [&](std::string_view mat1_e, std::string_view mat8_g13) {
        const auto deck = std::string{"MAT1,9,"} + std::string{mat1_e} +
                          ",80000.0,.25,7.85D-9\nMAT8,10,130000.0,9000.0,.28,6000.0," +
                          std::string{mat8_g13} + ",4000.0,1.56D-9\n";
        auto card_input = material_cards(deck);
        return take(
            NastranLaminateCards::prepare(std::move(card_input), laminate_limits, plan, {}));
    };

    const auto same_deck = check("200000.0", "5000.0");
    const auto same = take(same_deck.view());
    std::optional<ImportedQuantityMismatch> diagnostic;
    REQUIRE(check_imported_material(same.input.mat1.front(), case_mat1, deck_units, si,
                                    comparison_limits(), plan, {}, &diagnostic));
    REQUIRE(check_imported_material(same.mat8.front(), case_mat8, deck_units, si,
                                    comparison_limits(), plan, {}, &diagnostic));
    CHECK_FALSE(diagnostic);

    const auto changed_deck = check("200000.00000000001", "5100.0");
    const auto changed = take(changed_deck.view());
    const auto differing_young =
        check_imported_material(changed.input.mat1.front(), case_mat1, deck_units, si,
                                comparison_limits(), plan, {}, &diagnostic);
    REQUIRE_FALSE(differing_young);
    CHECK(differing_young.error().code == InvalidInput);
    const auto young_detail = take(diagnostic);
    CHECK(young_detail.quantity == ImportedQuantityField::Young1);
    CHECK(young_detail.card_position == take(changed.input.mat1.front().young).field.position);

    const auto differing_shear =
        check_imported_material(changed.mat8.front(), case_mat8, deck_units, si,
                                comparison_limits(), plan, {}, &diagnostic);
    REQUIRE_FALSE(differing_shear);
    CHECK(differing_shear.error().code == InvalidInput);
    const auto shear_detail = take(diagnostic);
    CHECK(shear_detail.quantity == ImportedQuantityField::Shear13);
    CHECK(shear_detail.card_position == take(changed.mat8.front().shear_1z).field.position);
}

TEST_CASE("A declared constant whole-region full thickness compares with every case PG vertex",
          "[IO4-03][IO4-04][N4-02][verification][negative]") {
    auto input = a::root();
    auto& values =
        input.field("model").field("geometry").elements.at(0).field("thickness").field("values");
    for (auto& value : values.elements) {
        value = a::parse("0.00254");
    }
    const auto constant = checked_case(input);
    const auto thickness = [](const SchemaCase& schema) {
        const auto model = take(schema.root()).model.record();
        const auto geometry = take(model.geometry.at(0)).record();
        return std::get<schema::GeometryPg>(geometry).thickness.record();
    };
    const auto expected = thickness(constant);
    constexpr ImportedUniformThickness declared{"0.1", SiLengthUnit::Inch};
    std::optional<ImportedQuantityMismatch> diagnostic;
    core::ResourcePlan plan{reconciliation_budget};
    REQUIRE(check_imported_uniform_thickness(declared, expected, si.length, expected.values.size(),
                                             comparison_limits(), plan, {}, &diagnostic));
    CHECK_FALSE(diagnostic);
    const auto insufficient =
        check_imported_uniform_thickness(declared, expected, si.length, expected.values.size() - 1,
                                         comparison_limits(), plan, {}, &diagnostic);
    REQUIRE_FALSE(insufficient);
    CHECK(insufficient.error().code == ResourceExceeded);
    CHECK_FALSE(diagnostic);

    values.elements.back() = a::parse("0.00254000000000000000000000000001");
    const auto varying = checked_case(input);
    const auto vertices = thickness(varying);
    const auto disagreement =
        check_imported_uniform_thickness(declared, vertices, si.length, vertices.values.size(),
                                         comparison_limits(), plan, {}, &diagnostic);
    REQUIRE_FALSE(disagreement);
    CHECK(disagreement.error().code == InvalidInput);
    const auto sidecar_detail = take(diagnostic);
    CHECK(sidecar_detail.quantity == ImportedQuantityField::FullThickness);
    CHECK(sidecar_detail.vertex == static_cast<core::Index>(vertices.values.size() - 1));

    const auto pshell = [](std::string_view thick, bool corner = false) {
        const auto element = corner ? "CTRIA3,50,91,1,2,3,,,\n,,,3.0\n" : "CTRIA3,50,91,1,2,3\n";
        const auto deck = std::string{"GRID,1\nGRID,2\nGRID,3\n"} + element + "PSHELL,91,9," +
                          std::string{thick} + "\nMAT1,9,1.0\n";
        return material_cards(deck);
    };
    const auto candidate = pshell("2.54");
    const auto source = take(candidate.view());
    REQUIRE(check_imported_pshell_thickness(
        source.input.elements.front(), source.pshell.front(), SiLengthUnit::Millimetre, expected,
        si.length, expected.values.size(), comparison_limits(), plan, {}, &diagnostic));
    const auto wrong = pshell("2.5400000000000000000000000000001");
    const auto wrong_source = take(wrong.view());
    const auto wrong_result = check_imported_pshell_thickness(
        wrong_source.input.elements.front(), wrong_source.pshell.front(), SiLengthUnit::Millimetre,
        expected, si.length, expected.values.size(), comparison_limits(), plan, {}, &diagnostic);
    REQUIRE_FALSE(wrong_result);
    CHECK(wrong_result.error().code == InvalidInput);
    const auto deck_detail = take(diagnostic);
    CHECK(deck_detail.card_position == take(wrong_source.pshell.front().thickness).field.position);
    CHECK(deck_detail.vertex == 0);

    const auto corner = pshell("2.54", true);
    const auto corner_source = take(corner.view());
    const auto unqualified = check_imported_pshell_thickness(
        corner_source.input.elements.front(), corner_source.pshell.front(),
        SiLengthUnit::Millimetre, expected, si.length, expected.values.size(), comparison_limits(),
        plan, {}, &diagnostic);
    REQUIRE_FALSE(unqualified);
    CHECK(unqualified.error().code == Unsupported);
    CHECK_FALSE(diagnostic);
}

[[nodiscard]] structure::PgReferenceGeometry pcomp_pg(partition_test::Fixture& input,
                                                      std::span<const std::string_view> endpoints,
                                                      std::string_view full_numerator,
                                                      std::string_view full_denominator,
                                                      bool varying, core::ResourcePlan& plan) {
    auto number = [&](std::string_view value) {
        const auto slash = value.find('/');
        if (slash != std::string_view::npos) {
            return input.r(value.substr(0, slash), value.substr(slash + 1));
        }
        return input.r(value);
    };
    std::vector<structure::ExactRational> zeta;
    zeta.reserve(endpoints.size());
    for (const auto value : endpoints) {
        zeta.push_back(number(value));
    }
    input.layers.clear();
    for (std::size_t i = 0; i + 1 < zeta.size(); ++i) {
        input.layers.emplace_back(core::LayerId{11U + i}, zeta[i], zeta[i + 1]);
    }
    for (auto& height : input.heights) {
        height = input.r(full_numerator, full_denominator);
    }
    if (varying) {
        input.heights[3] = input.r("41", "100");
    }
    const auto exact = partition_test::exact_limits();
    const structure::PgReferenceLimits reference{exact, 20'000'000, 10'000'000, {}};
    const auto one = input.r("1");
    const std::array origin{input.r("0"), input.r("0"), input.r("0")};
    const std::array rotation{input.r("1"), input.r("0"), input.r("0"), input.r("0"), input.r("1"),
                              input.r("0"), input.r("0"), input.r("0"), input.r("1")};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto dictionary = take(structure::PgCoefficientDictionary::prepare(
        input.chart_graph(), {exact, 2000, 20'000'000, 100'000, {}}, workspace, plan));
    auto partition = take(structure::PgIntegrationPartition::prepare(
        input.request(), partition_test::limits(), workspace, plan));
    input.layers.clear();
    structure::PgReferenceInputs assembled{std::move(input.chart_graph()), std::move(dictionary),
                                           std::move(partition)};
    return take(structure::PgReferenceGeometry::prepare(std::move(assembled),
                                                        {rotation, origin, one}, reference, plan));
}

TEST_CASE("PCOMP layers retain exact thickness and reverse the physical ply order",
          "[IO4-03][FD-03][R-05][N4-02][integration][verification][negative]") {
    constexpr auto forward = ExternalPgWinding::Forward;
    constexpr auto reverse = ExternalPgWinding::Reversed;
    constexpr std::array<std::string_view, 3> forward_intervals{"-1", "-1/2", "1"};
    constexpr std::array<std::string_view, 3> reverse_intervals{"-1", "1/2", "1"};
    constexpr std::array<std::string_view, 5> split_intervals{"-1", "-3/4", "0", "3/4", "1"};
    constexpr std::array<std::string_view, 4> joined_intervals{"-1", "-3/4", "3/4", "1"};
    auto endpoints = std::span<const std::string_view>{forward_intervals};
    std::string_view z0 = "-200.", mode, first = "100.", second = "300.";
    ExternalPgWinding direction = forward;
    bool variable_case = false;
    using Failure = NastranPcompGeometryFailure;
    std::optional<Failure> failure;
    SECTION("Explicit centered offset with asymmetric plies") {}
    SECTION("Reversed face uses the original top ply first") {
        direction = reverse;
        endpoints = reverse_intervals;
        z0 = {};
    }
    SECTION("An unchanged Case ply order cannot match a reversed asymmetric stack") {
        direction = reverse;
        failure = Failure::Interface;
    }
    SECTION("The exact layer interface catches a changed split with the same total") {
        first = "101.";
        second = "299.";
        failure = Failure::Interface;
    }
    SECTION("A sub-ulp deck thickness difference fails") {
        first = "100.00000000000000000001";
        z0 = {};
        failure = Failure::FullThickness;
    }
    SECTION("The actual positive-area PG thickness is compared") {
        variable_case = true;
        failure = Failure::FullThickness;
    }
    SECTION("An explicit bottom shift cannot match the prepared reference plane") {
        z0 = "-199.";
        failure = Failure::ReferencePlane;
    }
    SECTION("Both halves of SYM can occupy independently declared Case layers") {
        mode = "SYM";
        endpoints = split_intervals;
        z0 = "-400.";
    }
    SECTION("SYM matches one declared center layer spanning the center plane") {
        mode = "SYM";
        endpoints = joined_intervals;
        z0 = {};
    }
    SECTION("A smeared stiffness option does not certify a physical stack") {
        mode = "SMEAR";
        z0 = {};
        failure = Failure::Mode;
    }
    SECTION("A PCOMP stack does not merge unmatched Case layers") {
        endpoints = joined_intervals;
        failure = Failure::LayerCount;
    }
    const auto head =
        std::string{"PCOMP,30,"} + std::string{z0} + ",,,,,," + std::string{mode} + "\n";
    const auto bulk = head + ",81," + std::string{first} + ",,,52," + std::string{second} +
                      "\nMAT8,81,130.0,9.0,.28\nMAT1,52,1.0\n";
    core::ResourcePlan plan{partition_test::cap};
    const auto cards =
        take(NastranLaminateCards::prepare(material_cards(bulk), laminate_limits, plan, {}));
    const auto deck = take(cards.view());
    const auto expected_hash = source_digest(std::as_bytes(std::span{bulk}));
    const auto property = deck.pcomp.front().property_ordinal;
    partition_test::Fixture input;
    const bool symmetric = mode == "SYM";
    const auto pg = pcomp_pg(input, endpoints, symmetric ? "4" : "2", "5", variable_case, plan);
    const NastranPcompGeometryLimits bounds{comparison_limits(), 4, 128};
    std::optional<NastranPcompGeometryDiagnostic> diagnostic;
    const auto compare = [&](const NastranPcompGeometryLimits& quota,
                             const std::stop_token& stop = {}) {
        return compare_nastran_pcomp_geometry(cards, property, SiLengthUnit::Millimetre, direction,
                                              pg, quota, plan, stop, &diagnostic);
    };
    const auto result = compare(bounds);
    if (failure) {
        REQUIRE_FALSE(result);
        CHECK(result.error().operation == core::Operation::PrepareModel);
        CHECK(result.error().code == (failure == Failure::Mode ? Unsupported : InvalidInput));
        const auto invalid = take(diagnostic);
        CHECK(invalid.reason == *failure);
        CHECK(invalid.field);
        if (failure == Failure::Interface) {
            CHECK(invalid.case_layer == core::LayerId{11});
        }
        if (variable_case) {
            CHECK(invalid.thickness_vertex == 3);
        }
        return;
    }
    REQUIRE(result);
    CHECK_FALSE(diagnostic);
    CHECK(result->source.source == SourceId{99});
    CHECK(result->digest == expected_hash);
    CHECK(result->pg == take(pg.identity()));
    CHECK(result->winding == direction);
    CHECK(result->full_thickness_metres == input.r(symmetric ? "4" : "2", "5"));
    const auto match = result->layers.values();
    REQUIRE(match.size() == endpoints.size() - 1);
    for (std::size_t i = 0; i < match.size(); ++i) {
        CHECK(match[i].layer == core::LayerId{11U + i});
        CHECK(match[i].case_ordinal == i);
        const auto declared = match[i].declared_ply;
        REQUIRE(declared < 2);
        const auto ply = deck.layers(deck.pcomp.front());
        CHECK(match[i].material.id == take(ply[declared].material).id);
        CHECK(match[i].thickness == take(ply[declared].thickness).field.position);
    }
    if (symmetric) {
        CHECK(match.front().material.id == NastranMaterialId{81});
        CHECK(match.back().material.id == NastranMaterialId{81});
        CHECK(match[1].ply_thickness_metres == input.r(match.size() == 3 ? "6" : "3", "10"));
    } else {
        CHECK(match.front().material.id == NastranMaterialId{direction == forward ? 81U : 52U});
        CHECK(match.front().ply_thickness_metres ==
              input.r(direction == forward ? "1" : "3", "10"));
    }

    auto exact = bounds;
    exact.work = result->work;
    REQUIRE(compare(exact));
    --exact.work;
    const auto exhausted = compare(exact);
    REQUIRE_FALSE(exhausted);
    CHECK(exhausted.error().code == ResourceExceeded);
    exact = bounds;
    exact.thickness_vertices = 3;
    const auto vertices = compare(exact);
    REQUIRE_FALSE(vertices);
    CHECK(vertices.error().code == ResourceExceeded);
    const auto cancelled = compare(bounds, partition_test::requested_stop_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK_FALSE(diagnostic);
    CHECK(compare(bounds));
}
} // namespace
