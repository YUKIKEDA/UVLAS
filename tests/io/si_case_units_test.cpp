#include "case_schema_snapshot.hpp"
#include "case_schema_test_adapter.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <uvlas/io/si_case.hpp>
#include <variant>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
constexpr std::uint64_t budget = std::uint64_t{1} << 50U;
template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::runtime_error("SI unit control rejected, code=" +
                                 std::to_string(static_cast<unsigned>(value.error().code)));
    return std::move(*value);
}
template <class T>
T si_phase(core::Result<T> value, std::string_view phase, const io::SiDiagnostic& diag) {
    if (!value) {
        std::string message{phase};
        message += " code=" + std::to_string(static_cast<unsigned>(value.error().code));
        const auto detail = diag.detail();
        if (detail)
            message += " reason=" + std::to_string(static_cast<unsigned>(detail->reason));
        const auto path = diag.path();
        if (path)
            message += " path=" + std::string(*path);
        throw std::runtime_error(message);
    }
    return std::move(*value);
}
io::SiLimits unit_limits() {
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> q{};
    q.fill(1000000);
    for (const auto quota : {io::SiQuota::OwnedBytes, io::SiQuota::ScratchBytes,
                             io::SiQuota::CumulativeBytes, io::SiQuota::Work})
        q[static_cast<std::size_t>(quota)] = budget;
    q[static_cast<std::size_t>(io::SiQuota::UnitPower)] = 64;
    return {q,
            {100000, 100000, 100000, 64, 64, 1000000, 2000000,
             structure::ExactEstimateAdmission{budget}, structure::ExactTimingMode::Disabled}};
}
io::SchemaCase checked(const JsonToken& source) {
    const auto wire = source.wire();
    const io::InputLimits syntax{
        10000000, 256,      1000000,  1000000,
        10000000, 10000000, 10000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    core::ResourcePlan plan{budget};
    std::optional<io::InputDiagnostic> diag;
    const auto doc = take(io::read_input_document(
        std::as_bytes(std::span{wire.data(), wire.size()}), io::SourceId{716}, syntax, plan, diag));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> caps{};
    caps.fill(1000000000);
    io::SchemaDiagnostic schema;
    auto result = io::validate_case_schema(doc, io::SchemaLimits{caps}, plan, &schema);
    if (!result) {
        std::string message =
            "schema code=" + std::to_string(static_cast<unsigned>(result.error().code));
        const auto detail = schema.detail();
        if (detail)
            message += " reason=" + std::to_string(static_cast<unsigned>(detail->reason));
        const auto path = schema.path();
        if (path)
            message += " path=" + std::string(*path);
        throw std::runtime_error(message);
    }
    return std::move(*result);
}
std::vector<io::SiRequirementRow> requirements(const io::SiCatalog& catalog,
                                               std::optional<io::SiQuantityKey> exact = {}) {
    std::vector<io::SiRequirementRow> rows;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto entry = take(catalog.entry(i));
        const bool zero = exact && *exact == entry.key;
        rows.push_back({entry.key,
                        {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                        io::SiDecimalInput{zero ? "0" : "1e-25"},
                        io::SiDecimalInput{zero ? "0" : "1e-25"}});
    }
    return rows;
}
io::SiCase converted(const JsonToken& source, bool exact_length = false) {
    const auto input = checked(source);
    core::ResourcePlan plan{budget};
    io::SiDiagnostic diag;
    const auto limits = unit_limits();
    const auto catalog = si_phase(io::make_si_catalog(input, limits, plan, &diag), "catalog", diag);
    std::optional<io::SiQuantityKey> exact;
    if (exact_length) {
        for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
            const auto entry = take(catalog.entry(i));
            if (entry.field == "reference_length")
                exact = entry.key;
        }
        if (!exact)
            throw std::runtime_error("Independent length fixture did not reach catalog");
    }
    const auto rows = requirements(catalog, exact);
    const auto request =
        si_phase(io::make_si_request(catalog, rows, limits, plan, &diag), "request", diag);
    return si_phase(io::convert_case_to_si(input, catalog, request, limits, plan, &diag), "convert",
                    diag);
}
// Decode each IEEE word by its bits; never add high/low or rescale in double.
Rat binary(double word) {
    const auto bits = std::bit_cast<std::uint64_t>(word);
    const auto exponent = (bits >> 52U) & 0x7ffU;
    if (exponent == 0x7ffU)
        throw std::runtime_error("Nonfinite SI result");
    Big significand = bits & 0xfffffffffffffULL;
    if (exponent != 0)
        significand += Big{1} << 52U;
    if (significand == 0)
        return Rat{0};
    Rat value{significand};
    const auto power = exponent == 0 ? -1074 : static_cast<int>(exponent) - 1075;
    if (power >= 0)
        value *= Big{1} << static_cast<unsigned>(power);
    else
        value /= Big{1} << static_cast<unsigned>(-power);
    if ((bits >> 63U) != 0)
        return -value;
    return value;
}
std::pair<Rat, Rat> interval(const io::SiQuantityView& quantity) {
    return std::visit(
        [](const auto& x) {
            using X = std::decay_t<decltype(x)>;
            Rat center;
            if constexpr (std::is_same_v<X, io::ScaledBinary64>)
                center = binary(x.center);
            else
                center = binary(x.high) + binary(x.low);
            const Rat radius = binary(x.radius);
            CHECK(radius >= 0);
            REQUIRE(x.exponent2 > -4096);
            REQUIRE(x.exponent2 < 4096);
            Rat scale = 1;
            if (x.exponent2 >= 0)
                scale *= Big{1} << static_cast<unsigned>(x.exponent2);
            else
                scale /= Big{1} << static_cast<unsigned>(-x.exponent2);
            return std::pair<Rat, Rat>{(center - radius) * scale, (center + radius) * scale};
        },
        quantity.enclosure);
}
void contains(const io::SiQuantityView& quantity, const Rat& lower, const Rat& upper) {
    const auto bound = interval(quantity);
    CHECK(bound.first <= lower);
    CHECK(bound.second >= upper);
    // Independently enforce the request's SI scale=1 and a=r=10^-25.
    const Rat magnitude_lower =
        lower <= 0 && upper >= 0 ? Rat{0} : Rat{std::min(abs(lower), abs(upper))};
    const Rat allowed = (1 + magnitude_lower) / boost::multiprecision::pow(Big{10}, 25);
    CHECK((bound.second - bound.first) / 2 <= allowed);
    if (quantity.recipe.pi_power == 0) {
        REQUIRE(lower == upper);
        core::ResourcePlan plan{budget};
        const auto config = unit_limits();
        const auto rational = take(io::materialize_rational(quantity.recipe, config.exact(), plan));
        const auto expected = take(structure::ExactRational::from_integer_pair(
            numerator(lower).str(), denominator(lower).str(), config.exact(), plan));
        CHECK(rational == expected);
    }
}
void exact_value(const io::SiQuantityView& quantity, const Rat& value) {
    INFO("dims " << quantity.recipe.dimensions.exponents[0] << ','
                 << quantity.recipe.dimensions.exponents[1] << ','
                 << quantity.recipe.dimensions.exponents[2] << ','
                 << quantity.recipe.dimensions.exponents[3]);
    const auto description = std::visit(
        [](const auto& x) {
            using X = std::decay_t<decltype(x)>;
            std::string text;
            if constexpr (std::is_same_v<X, io::ScaledBinary64>)
                text = "center=" + std::to_string(std::bit_cast<std::uint64_t>(x.center));
            else
                text = "hi=" + std::to_string(std::bit_cast<std::uint64_t>(x.high)) +
                       " low=" + std::to_string(std::bit_cast<std::uint64_t>(x.low));
            return text + " radius=" + std::to_string(std::bit_cast<std::uint64_t>(x.radius)) +
                   " q=" + std::to_string(x.exponent2);
        },
        quantity.enclosure);
    INFO(description);
    const auto bound = interval(quantity);
    CHECK(numerator(bound.first) * denominator(value) ==
          numerator(value) * denominator(bound.first));
    CHECK(numerator(bound.second) * denominator(value) ==
          numerator(value) * denominator(bound.second));
    CHECK(bound.first == value);
    CHECK(bound.second == value);
}
io::schema::Case case_root(const io::SiCase& value) {
    return take(take(value.source_case()).get().root());
}
io::schema::GeometryPg geometry(const io::SiCase& value) {
    const auto payload = take(case_root(value).model.record().geometry.at(0)).record();
    return std::get<io::schema::GeometryPg>(payload);
}
io::schema::MaterialIsotropic material(const io::SiCase& value) {
    return std::get<io::schema::MaterialIsotropic>(
        take(case_root(value).model.record().materials.at(0)).record());
}
std::pair<Rat, Rat> pi_bounds() {
    // Independently exact rational alternating sums, with the next-term theorem.
    const auto atan_bounds = [](unsigned d) {
        Rat sum = 0;
        for (unsigned j = 0; j < 160; ++j) {
            const Rat term = Rat{1} / ((2 * j + 1) * boost::multiprecision::pow(Big{d}, 2 * j + 1));
            sum += (j & 1U) ? -term : term;
        }
        const Rat remainder = Rat{1} / (321 * boost::multiprecision::pow(Big{d}, 321));
        return std::pair<Rat, Rat>{sum, sum + remainder};
    };
    const auto a = atan_bounds(5), b = atan_bounds(239);
    return {16 * a.first - 4 * b.second, 16 * a.second - 4 * b.first};
}
JsonToken beam_case() {
    auto value = root();
    value.field("model").field("beams").elements.push_back(fragment("Beam"));
    auto support = fragment("Support.elastic");
    support.field("target").field("target").set("kind", string("beam"));
    support.field("target").field("target").set("id", string("1"));
    support.set("id", string("2"));
    value.field("model").field("supports").elements.push_back(std::move(support));
    // Adding a beam requires its translation/rotation and moment balance scales.
    // Elastic support adds no orientation multiplier block; retain the original
    // shell trace constraint's position-only constraint scale policy.
    for (const auto field : {"beam_translation", "beam_rotation", "moment_balance"})
        value.field("numerics").field("scales").elements.push_back(scale(field));
    return value;
}
} // namespace

TEST_CASE("Every declared SI base unit uses independent exact physical factors",
          "[IO4-01][verification][property]") {
    const std::array<std::string_view, 6> lengths{"m", "mm", "cm", "km", "in", "ft"};
    const std::array<Rat, 6> lfactor{Rat{1},    Rat{1} / 1000,   Rat{1} / 100,
                                     Rat{1000}, Rat{127} / 5000, Rat{381} / 1250};
    for (std::size_t i = 0; i < lengths.size(); ++i) {
        INFO(lengths[i]);
        auto source = root();
        source.field("units").set("length", string(lengths[i]));
        source.field("model").field("geometry").elements[0].set("reference_length", parse("37"));
        source.field("model").field("materials").elements[0].set("density", parse("311"));
        const auto value = converted(source);
        contains(take(value.quantity(geometry(value).reference_length)), 37 * lfactor[i],
                 37 * lfactor[i]);
        const Rat density = Rat{311} / (lfactor[i] * lfactor[i] * lfactor[i]);
        contains(take(value.quantity(material(value).density)), density, density);
        // Chart geometry remains dimensionless even when physical length changes.
        const auto width = take(value.quantity(geometry(value).width));
        contains(width, Rat{1}, Rat{1});
        CHECK(width.recipe.pi_power == 0);
    }
    const std::array<std::string_view, 4> masses{"kg", "g", "t", "lbm"};
    const std::array<Rat, 4> mfactor{Rat{1}, Rat{1} / 1000, Rat{1000}, Rat{45359237} / 100000000};
    for (std::size_t i = 0; i < masses.size(); ++i) {
        INFO(masses[i]);
        auto source = root();
        source.field("units").set("mass", string(masses[i]));
        source.field("model").field("materials").elements[0].set("density", parse("311"));
        const auto value = converted(source);
        const Rat expected = 311 * mfactor[i];
        contains(take(value.quantity(material(value).density)), expected, expected);
    }
    const std::array<std::string_view, 3> times{"s", "ms", "min"};
    const std::array<Rat, 3> tfactor{Rat{1}, Rat{1} / 1000, Rat{60}};
    for (std::size_t i = 0; i < times.size(); ++i) {
        INFO(times[i]);
        auto source = root("coupled_transient");
        source.field("units").set("time", string(times[i]));
        source.field("analysis").field("intervals").elements[0].set("end", parse("37"));
        source.field("model").field("materials").elements[0].set("young", parse("17"));
        const auto value = converted(source);
        const Rat pressure = Rat{17} / (tfactor[i] * tfactor[i]);
        contains(take(value.quantity(material(value).young)), pressure, pressure);
        const auto analysis =
            std::get<io::schema::AnalysisTransient>(case_root(value).analysis.record());
        const auto interval_payload = take(analysis.intervals.at(0)).record();
        contains(take(value.quantity(interval_payload.end)), 37 * tfactor[i], 37 * tfactor[i]);
    }
}
TEST_CASE("Exact zero budgets accept unit-cancelled dyadics without an interval floor",
          "[IO4-01][verification][property]") {
    const std::array<std::pair<std::uint64_t, Rat>, 7> words{
        {{0, Rat{0}},
         {0x8000000000000000ULL, Rat{0}},
         {0x3ff9000000000000ULL, Rat{25} / 16},
         {0x3fe0000000000000ULL, Rat{1} / 2},
         {0xbff4000000000000ULL, Rat{-5} / 4},
         {0x3fb999999999999aULL, Rat{3602879701896397LL} / Big{36028797018963968LL}},
         {1, Rat{1} / (Big{1} << 1074U)}}};
    for (const auto& [word, expected] : words) {
        const auto decoded = binary(std::bit_cast<double>(word));
        CHECK(numerator(decoded) * denominator(expected) ==
              numerator(expected) * denominator(decoded));
        CHECK(decoded == expected);
    }
    for (const auto choice : {0, 1}) {
        auto source = root();
        source.field("units").set("length", string(choice == 0 ? "km" : "mm"));
        source.field("model")
            .field("geometry")
            .elements[0]
            .set("reference_length", parse(choice == 0 ? "0.1" : "500"));
        const auto value = converted(source, true);
        exact_value(take(value.quantity(geometry(value).reference_length)),
                    choice == 0 ? Rat{100} : Rat{1} / 2);
    }
}
TEST_CASE("Beam curvature and support matrices preserve mixed dimensions and inverse pi",
          "[IO4-01][verification][property]") {
    auto source = beam_case();
    auto& units = source.field("units");
    units.set("length", string("mm"));
    units.set("mass", string("g"));
    units.set("time", string("ms"));
    units.set("angle", string("deg"));
    auto& beam = source.field("model").field("beams").elements[0];
    beam.set("initial_strain", parse("[2,2,2,2,2,2]"));
    auto& support = source.field("model").field("supports").elements[1];
    support.field("stiffness").elements[0].elements[3] = parse("7");
    support.field("stiffness").elements[3].elements[0] = parse("7");
    support.field("reference").field("rotation_vector").set("value", parse("[37,0,0]"));
    const auto value = converted(source);
    const auto model = case_root(value).model.record();
    const auto b = take(model.beams.at(0)).record();
    for (core::Index i = 0; i < 6; ++i) {
        const auto q = take(value.quantity(take(b.initial_strain.at(i))));
        const Rat expected = i < 3 ? Rat{2} : Rat{2000};
        contains(q, expected, expected);
        CHECK(q.recipe.pi_power == 0);
        const io::SchemaDimension expected_dim{
            {static_cast<std::int16_t>(i < 3 ? 0 : -1), 0, 0, 0}};
        CHECK(q.recipe.dimensions == expected_dim);
    }
    const auto k = std::get<io::schema::SectionMatrix>(b.section.record()).stiffness;

    for (core::Index i = 0; i < 6; ++i) {
        const auto q = take(value.quantity(take(take(k.at(i)).at(i))));
        const Rat expected = i < 3 ? Rat{1} : Rat{1} / 1000000;
        contains(q, expected, expected);
        CHECK(q.recipe.pi_power == 0);
    }
    const auto e = std::get<io::schema::SupportElastic>(take(model.supports.at(1)).record());
    const auto pi = pi_bounds();
    contains(take(value.quantity(take(take(e.stiffness.at(0)).at(0)))), Rat{1000}, Rat{1000});
    const auto cross = take(value.quantity(take(take(e.stiffness.at(0)).at(3))));
    CHECK(cross.recipe.pi_power == -1);
    contains(cross, Rat{1260} / pi.second, Rat{1260} / pi.first);
    const auto rotational = take(value.quantity(take(take(e.stiffness.at(3)).at(3))));
    CHECK(rotational.recipe.pi_power == -2);
    contains(rotational, Rat{162} / (5 * pi.second * pi.second),
             Rat{162} / (5 * pi.first * pi.first));
    const auto rotation =
        std::get<io::schema::TimeSeriesConstant>(e.reference.record().rotation_vector.record());
    const auto angle = take(value.quantity(take(rotation.value.at(0))));
    CHECK(angle.recipe.pi_power == 1);
    contains(angle, 37 * pi.first / 180, 37 * pi.second / 180);
    core::ResourcePlan plan{budget};
    const auto config = unit_limits();
    const auto rational = io::materialize_rational(angle.recipe, config.exact(), plan);
    REQUIRE_FALSE(rational);
    CHECK(rational.error().code == core::FailureCode::PrecisionUnsupported);
    const auto zero_angle = take(value.quantity(take(rotation.value.at(1))));
    exact_value(zero_angle, Rat{0});
    const auto zero_rational =
        take(io::materialize_rational(zero_angle.recipe, config.exact(), plan));
    CHECK(zero_rational.sign() == 0);
}
TEST_CASE("SI host requests reject missing duplicate dimensions and invalid signs",
          "[IO4-01][negative][property]") {
    const auto input = checked(root());
    const auto limits = unit_limits();
    core::ResourcePlan plan{budget};
    io::SiDiagnostic diag;
    const auto catalog = si_phase(io::make_si_catalog(input, limits, plan, &diag), "catalog", diag);
    const auto correct = requirements(catalog);
    REQUIRE(correct.size() > 1);
    for (const auto defect : {0, 1, 2, 3, 4}) {
        INFO(defect);
        auto rows = correct;
        if (defect == 0)
            rows.pop_back();
        if (defect == 1)
            rows.back() = rows.front();
        if (defect == 2)
            rows.front().scale.dimensions.exponents[0] = 17;
        if (defect == 3)
            rows.front().scale.value = io::SiDecimalInput{"0"};
        if (defect == 4)
            rows.front().absolute = io::SiDecimalInput{"-1"};
        const auto request = io::make_si_request(catalog, rows, limits, plan, &diag);
        REQUIRE_FALSE(request);
        CHECK(request.error().code == core::FailureCode::InvalidInput);
    }
    const auto control = io::make_si_request(catalog, correct, limits, plan, &diag);
    CHECK(control.has_value());
}

TEST_CASE(
    "All five frozen roots preserve their complete typed payload through public SI conversion",
    "[IO4-01][verification][property]") {
    std::uint64_t count = 0;
    for (const auto& example : test::case_schema_corpus::roots) {
        INFO(example.name);
        const auto source = parse(example.json);
        const auto original = checked(source);
        const auto expected = test::schema_snapshot::serialize(take(original.root()));
        const auto converted_case = converted(source);
        const auto owned = take(converted_case.source_case());
        const auto actual = test::schema_snapshot::serialize(take(owned.get().root()));
        CHECK(actual == expected);
        ++count;
    }
    CHECK(count == 5);
}
