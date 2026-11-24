#include "nastran_isotropic_shell_cards.hpp"
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
using Reason = NastranIsotropicShellFailure;

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

constexpr NastranIsotropicShellLimits material_limits{16, 20000};
constexpr std::string_view isotropic = "MAT1,52,2.10D+11,,.31\n";
constexpr std::string_view grid_shell = "GRID,1\nGRID,2\nGRID,3\nCTRIA3,72,20,1,2,3,,-.01\n";

NastranShellCards shells(std::string_view source) {
    auto deck = lexical_deck(source);
    core::ResourcePlan plan{memory};
    auto coords = take(NastranCoordinateCards::prepare(std::move(deck), limits, plan, {}));
    return take(NastranShellCards::prepare(std::move(coords), {16, 16, 20000}, plan, {}));
}

TEST_CASE("MAT1 and PSHELL keep exact fields and shell-to-material identities from one snapshot",
          "[IO4-03][IO4-01][verification]") {
    const auto source = std::string{grid_shell} +
                        "PSHELL,20,52,2.5D-3,52,1.0,52,.833333,.1\n"
                        ",-.002,.003\n"
                        "PCOMP,21\n"
                        "PSHELL,22,52,.002,55,1.25,54,.75,0.0\n"
                        ",-.002,.003,56\n"
                        "mat1,52,2.10D+11,8.0+10,.3125,7.850+3,1.2-5,293.15,.02\n"
                        ",2.0+8,1.8+8,1.0+8,41\n"
                        "MAT8,54\nMAT1,55,3.0\nMAT1,56,,4.0\n"
                        "CORD2R,41,,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0\n";
    core::ResourcePlan plan{memory};
    const char* original = nullptr;
    const auto cards = [&] {
        auto input = shells(source);
        original = take(input.view()).input.bulk.bytes.data();
        auto typed =
            take(NastranIsotropicShellCards::prepare(std::move(input), material_limits, plan, {}));
        // NOLINTNEXTLINE(bugprone-use-after-move): success leaves the input view absent.
        CHECK_FALSE(input.view());
        return typed;
    }();
    const auto view = take(cards.view());
    const auto& bulk = view.input.input.bulk;
    CHECK(bulk.bytes.data() == original);
    CHECK(bulk.bytes == source);
    CHECK(bulk.digest == source_digest(std::as_bytes(std::span{source})));
    REQUIRE(view.materials.size() == 4);
    REQUIRE(view.mat1.size() == 3);
    REQUIRE(view.pshell.size() == 2);
    REQUIRE(view.input.elements.size() == 1);
    const auto& used = view.input.elements.front();
    const auto& pshell = view.pshell[0];
    CHECK(used.property.property_ordinal == pshell.property_ordinal);
    CHECK(view.input.properties[pshell.property_ordinal].id == NastranPropertyId{20});
    CHECK(view.input.properties[view.pshell[1].property_ordinal].id == NastranPropertyId{22});
    const auto membrane = take(pshell.membrane);
    CHECK(membrane.id == NastranMaterialId{52});
    const auto& material = view.materials[membrane.declaration_ordinal];
    CHECK(material.id == membrane.id);
    CHECK(material.kind == NastranMaterialKind::Mat1);
    CHECK(bulk.text(bulk.cards[material.card_ordinal].name) == "mat1");
    const auto& mat = view.mat1[take(material.mat1_ordinal)];
    CHECK(mat.declaration_ordinal == membrane.declaration_ordinal);
    const auto check_real = [&](const std::optional<NastranShellReal>& field,
                                std::string_view spelling, std::string_view decimal) {
        const auto value = take(field);
        CHECK(bulk.text(value.field) == spelling);
        CHECK(value.value.decimal() == decimal);
    };
    check_real(mat.young, "2.10D+11", "2.10e+11");
    check_real(mat.shear, "8.0+10", "8.0e+10");
    check_real(mat.poisson, ".3125", "0.3125");
    check_real(mat.density, "7.850+3", "7.850e+3");
    check_real(mat.expansion, "1.2-5", "1.2e-5");
    check_real(mat.reference_temperature, "293.15", "293.15");
    check_real(mat.damping, ".02", "0.02");
    check_real(mat.tension_limit, "2.0+8", "2.0e+8");
    check_real(mat.compression_limit, "1.8+8", "1.8e+8");
    check_real(mat.shear_limit, "1.0+8", "1.0e+8");
    const auto axis = take(mat.material_coordinate);
    CHECK(axis.id == NastranCoordinateId{41});
    CHECK(view.input.input.coordinates[take(axis.coordinate_ordinal)].id == axis.id);
    CHECK(bulk.text(axis.field) == "41");
    check_real(pshell.thickness, "2.5D-3", "2.5e-3");
    check_real(pshell.bending_ratio, "1.0", "1.0");
    check_real(pshell.shear_ratio, ".833333", "0.833333");
    check_real(pshell.nonstructural_mass, ".1", "0.1");
    check_real(pshell.first_fiber, "-.002", "-0.002");
    check_real(pshell.second_fiber, ".003", "0.003");
    const auto fiber = take(pshell.first_fiber);
    CHECK(fiber.field.position == NastranBulkPosition{source.find("-.002"), 6, 2});
    CHECK(std::get<NastranMaterialReference>(pshell.bending).id == membrane.id);
    CHECK(take(pshell.transverse_shear).id == membrane.id);
    CHECK_FALSE(pshell.distinct_materials);
    CHECK_FALSE(pshell.orthotropic_material);
    CHECK_FALSE(pshell.coupling);
    const auto& mixed = view.pshell[1];
    CHECK(mixed.distinct_materials);
    CHECK(mixed.orthotropic_material);
    CHECK(take(mixed.coupling).id == NastranMaterialId{56});
    const auto shear = take(mixed.transverse_shear);
    CHECK(view.materials[shear.declaration_ordinal].kind == NastranMaterialKind::Mat8);
    CHECK_FALSE(view.materials[shear.declaration_ordinal].mat1_ordinal);
}

TEST_CASE("MAT1 large halves and PSHELL small continuation keep the documented field positions",
          "[IO4-03][verification]") {
    const auto fixed = [](std::string_view name, std::size_t width,
                          std::initializer_list<std::string_view> values,
                          std::string_view label = {}) {
        std::string line{name};
        line.resize(8, ' ');
        for (const auto field : values) {
            REQUIRE(field.size() <= width);
            line.append(field);
            line.append(width - field.size(), ' ');
        }
        line.resize(72, ' ');
        return line + std::string{label} + "\n";
    };
    const auto source =
        fixed("MAT1*", 16, {"53", "1.23456789012D+1", "", ".30"}, "*M1") +
        fixed("*M1", 16, {"7.800+3", "", "293.15", ""}, "*M2") +
        fixed("*M2", 16, {"8.0+2", "", "3.0+2", ""}, "*M3") + fixed("*M3", 16, {}) +
        fixed("PSHELL", 8, {"20", "53", ".125", "53", "1.0", "53", ".833333", ".0"}, "+P1") +
        fixed("+P1", 8, {"-.0625", ".0625"});
    auto parent = shells(source);
    core::ResourcePlan plan{memory};
    const auto typed =
        take(NastranIsotropicShellCards::prepare(std::move(parent), material_limits, plan, {}));
    const auto view = take(typed.view());
    REQUIRE(view.mat1.size() == 1);
    REQUIRE(view.pshell.size() == 1);
    const auto young = take(view.mat1[0].young);
    CHECK(young.value.decimal() == "1.23456789012e+1");
    CHECK(young.field.position == NastranBulkPosition{source.find("1.23456789012D+1"), 1, 25});
    const auto shear_limit = take(view.mat1[0].shear_limit);
    CHECK(shear_limit.value.decimal() == "3.0e+2");
    CHECK(shear_limit.field.position.line == 3);
    const auto fiber = take(view.pshell[0].first_fiber);
    CHECK(fiber.value.decimal() == "-0.0625");
    CHECK(fiber.field.position == NastranBulkPosition{source.find("-.0625"), 6, 9});
}

TEST_CASE("Optional values and a PSHELL bending designation do not invent a 3D material",
          "[IO4-03][verification]") {
    const auto source = std::string{isotropic} + "PSHELL,21,52\n"
                                                 "PSHELL,22,52,,-1\n"
                                                 "PSHELL,23,52,,0,,52\n"
                                                 "MAT1,53,,2.0,-.99999999999999\n"
                                                 "MAT1,54,5.0,,4.9-1\n"
                                                 "MAT1,55,5.0,,7.0D-18446744073709551616\n"
                                                 "MAT1,56,5.0,,-0.0+18446744073709551616\n";
    auto input = shells(source);
    core::ResourcePlan plan{memory};
    const auto typed =
        take(NastranIsotropicShellCards::prepare(std::move(input), material_limits, plan, {}));
    const auto view = take(typed.view());
    REQUIRE(view.mat1.size() == 5);
    REQUIRE(view.pshell.size() == 3);
    CHECK_FALSE(view.mat1[0].shear);
    CHECK_FALSE(view.mat1[0].density);
    CHECK_FALSE(view.mat1[0].expansion);
    CHECK_FALSE(view.mat1[0].reference_temperature);
    CHECK_FALSE(view.mat1[0].material_coordinate);
    CHECK_FALSE(view.mat1[1].young);
    const auto tiny = take(view.mat1[3].poisson);
    CHECK(tiny.value.decimal() == "7.0e-18446744073709551616");
    CHECK(std::holds_alternative<std::monostate>(view.pshell[0].bending));
    CHECK_FALSE(view.pshell[0].thickness);
    CHECK_FALSE(view.pshell[0].bending_ratio);
    CHECK_FALSE(view.pshell[0].shear_ratio);
    CHECK_FALSE(view.pshell[0].first_fiber);
    CHECK_FALSE(view.pshell[0].second_fiber);
    const auto plane = std::get<NastranBendingDesignation>(view.pshell[1].bending);
    CHECK(plane.kind == NastranSpecialBending::PlaneStrain);
    CHECK(view.input.input.bulk.text(plane.field) == "-1");
    const auto zero = std::get<NastranBendingDesignation>(view.pshell[2].bending);
    CHECK(zero.kind == NastranSpecialBending::Zero);
    CHECK(view.input.input.bulk.text(zero.field) == "0");
    CHECK(take(view.pshell[2].transverse_shear).id == NastranMaterialId{52});
}

TEST_CASE("MAT1 and PSHELL diagnostics preserve the shell owner and locate invalid fields",
          "[IO4-03][negative][verification]") {
    struct FailureCase {
        std::string_view text;
        Reason result;
        std::size_t line;
        std::size_t column;
    };
    constexpr std::array cases{
        FailureCase{"MAT1,53,,,,7.850+3\n", Reason::MissingField, 2, 9},
        FailureCase{"MAT1,53,0.0\n", Reason::InvalidField, 2, 9},
        FailureCase{"MAT1,53,,0.0\n", Reason::InvalidField, 2, 10},
        FailureCase{"MAT1,53,1.0,,.500000000001\n", Reason::InvalidField, 2, 14},
        FailureCase{"MAT1,53,1.0,,0.050+1\n", Reason::InvalidField, 2, 14},
        FailureCase{"MAT1,53,1.0,,0.05D+18446744073709551616\n", Reason::InvalidField, 2, 14},
        FailureCase{"MAT1,53,1.0,,-1.00-0\n", Reason::InvalidField, 2, 14},
        FailureCase{"MAT1,53,1.0,,-.100+1\n", Reason::InvalidField, 2, 14},
        FailureCase{"MAT1,53,1.0,,,,,,\n,,,,12\n", Reason::UndefinedCoordinate, 3, 5},
        FailureCase{"MAT1,53,1.0,,,,,,\n,,,,,12\n", Reason::ExtraField, 3, 6},
        FailureCase{"MAT8,52\n", Reason::DuplicateMaterial, 2, 6},
        FailureCase{"MAT8,100000000\n", Reason::InvalidField, 2, 6},
        FailureCase{"PSHELL,20,77\n", Reason::UndefinedMaterial, 2, 11},
        FailureCase{"PSHELL,20,52,,77\n", Reason::UndefinedMaterial, 2, 15},
        FailureCase{"PSHELL,20,52,,52,,77\n", Reason::UndefinedMaterial, 2, 19},
        FailureCase{"PSHELL,20,52,,52,,,,\n,,,77\n", Reason::UndefinedMaterial, 3, 4},
        FailureCase{"PSHELL,20,52,,-2\n", Reason::InvalidField, 2, 15},
        FailureCase{"PSHELL,20,,-1,-1\n", Reason::InvalidField, 2, 12},
        FailureCase{"PSHELL,20,,,-1\n", Reason::IncompatibleMaterial, 2, 11},
        FailureCase{"MAT8,53\nPSHELL,20,53,,-1\n", Reason::IncompatibleMaterial, 3, 11},
        FailureCase{"PSHELL,20,52,,,52,52\n", Reason::InvalidField, 2, 16},
        FailureCase{"PSHELL,20,52,,,,52\n", Reason::IncompatibleMaterial, 2, 17},
        FailureCase{"PSHELL,20,52,,-1,,52\n", Reason::IncompatibleMaterial, 2, 19},
        FailureCase{"PSHELL,20,52,,52,0.0\n", Reason::InvalidField, 2, 18},
        FailureCase{"PSHELL,20,52,,52,,52,-0.8\n", Reason::InvalidField, 2, 22},
        FailureCase{"PSHELL,20,52,,52,,,,\n,,,52\n", Reason::IncompatibleMaterial, 3, 4},
        FailureCase{"PSHELL,20,52,,52,,,,\n,,,,5\n", Reason::ExtraField, 3, 5}};
    for (const auto& entry : cases) {
        CAPTURE(entry.text);
        const auto source = std::string{isotropic} + std::string{entry.text};
        auto parent = shells(source);
        core::ResourcePlan plan{memory};
        std::optional<NastranIsotropicShellDiagnostic> diagnostic;
        const auto result = NastranIsotropicShellCards::prepare(std::move(parent), material_limits,
                                                                plan, {}, &diagnostic);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == InvalidInput);
        const auto actual = take(diagnostic);
        CHECK(actual.reason == entry.result);
        CHECK(actual.position.line == entry.line);
        CHECK(actual.position.column == entry.column);
        // NOLINTNEXTLINE(bugprone-use-after-move): rejected calls preserve the input.
        CHECK(take(parent.view()).input.bulk.bytes == source);
    }
}

TEST_CASE("Nastran material qualification accounts for retained rows, indexes and cancellation",
          "[IO4-03][N4-02][negative][verification]") {
    const auto source = std::string{isotropic} + "MAT8,53\nPSHELL,20,52,0.1,52\n";
    const auto run = [&](NastranIsotropicShellLimits capacity, core::ResourcePlan& plan,
                         const std::stop_token& stop = {}) {
        auto parent = shells(source);
        auto result = NastranIsotropicShellCards::prepare(std::move(parent), capacity, plan, stop);
        if (!result) {
            // NOLINTNEXTLINE(bugprone-use-after-move): rejected calls preserve input.
            CHECK(take(parent.view()).input.bulk.bytes == source);
        }
        return result;
    };
    core::ResourcePlan full{memory};
    const auto limited = run({1, material_limits.work}, full);
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == ResourceExceeded);
    CHECK(full.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness) MSVC request_stop() is non-const.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    const auto cancelled = run(material_limits, full, cancel.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Cancelled);
    CHECK(full.allocated_bytes() == 0);
    const auto accepted = run({2, material_limits.work}, full);
    REQUIRE(accepted);
    const auto work = take(accepted->view()).work;
    REQUIRE(work > 0);
    const auto bytes = full.allocated_bytes();
    REQUIRE(bytes > 0);
    core::ResourcePlan exact{bytes};
    const auto verified = run({2, work}, exact);
    CHECK(verified.has_value());
    CHECK(exact.remaining_bytes() == 0);
    core::ResourcePlan short_memory{bytes - 1};
    const auto byte_rejected = run({2, work}, short_memory);
    REQUIRE_FALSE(byte_rejected);
    CHECK(byte_rejected.error().code == ResourceExceeded);
    CHECK(short_memory.allocated_bytes() == 0);
    core::ResourcePlan short_work{bytes};
    const auto work_rejected = run({2, work - 1}, short_work);
    REQUIRE_FALSE(work_rejected);
    CHECK(work_rejected.error().code == ResourceExceeded);
}

} // namespace
