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
io::schema::Case case_root(const io::SiCase& value) {
    return take(take(value.source_case()).get().root());
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

void expect_schema_reject(const JsonToken& source) {
    try {
        (void)checked(source);
    } catch (const std::runtime_error& error) {
        CHECK(std::string_view{error.what()}.starts_with("schema code=0"));
        return;
    }
    FAIL("The adopted schema-negative input was accepted");
}
void dim(const io::SiQuantityView& q, io::SchemaDimension expected) {
    CHECK(q.recipe.dimensions == expected);
}
template <io::detail::FieldType T>
std::uint64_t ordinal(io::schema::Field<T> f) {
    if constexpr (io::detail::field_is_signed[static_cast<std::size_t>(T)])
        return static_cast<std::uint64_t>(f.signed_value());
    else
        return f.unsigned_value();
}
JsonToken expression_source(unsigned consumer) {
    auto source = root(consumer == 1 ? "coupled_transient" : "structural_static");
    source.field("units").set("length", string("mm"));
    source.field("units").set("mass", string("g"));
    source.field("units").set("time", string("ms"));
    auto expression = parse(
        R"({"domain":{},"nodes":[{"op":"variable","axis":0},{"op":"constant","value":37,"dimensions":[0,0,0,0]},{"op":"multiply","first":0,"second":1}],"outputs":[2]})");
    if (consumer == 0) {
        expression.set(
            "domain",
            parse(
                R"({"kind":"reference","coordinates":"chart","lower":[0,0,-1],"upper":[1,1,1]})"));
        expression.set("outputs", parse("[2,2,2,2,2,2]"));
        auto field = parse(R"({"kind":"expression","expression":{},"error":[0,0,0,0,0,0]})");
        field.set("expression", expression);
        source.field("model").field("shells").elements[0].field("layers").elements[0].set(
            "initial_green", field);
    } else if (consumer == 1) {
        expression.set("domain", parse(R"({"kind":"time","lower":0,"upper":1,"events":[0,1]})"));
        expression.field("nodes").elements[1].set("dimensions", parse("[1,0,-2,0]"));
        expression.set("outputs", parse("[2,2,2]"));
        auto field = parse(R"({"kind":"expression","expression":{}})");
        field.set("expression", expression);
        source.field("flow").set("uniform", field);
    } else {
        expression.set(
            "domain",
            parse(
                R"({"kind":"current_position_time","position_lower":[-1,-2,-3],"position_upper":[4,5,6],"time_lower":0,"time_upper":1,"events":[]})"));
        expression.field("nodes").elements[1].set("dimensions", parse("[-2,1,-2,0]"));
        auto load = fragment("Load.current_pressure");
        load.set("expression", expression);
        source.field("analysis").field("loads").elements.push_back(load);
    }
    return source;
}
JsonToken& expression_json(JsonToken& source, unsigned consumer) {
    if (consumer == 0)
        return source.field("model")
            .field("shells")
            .elements[0]
            .field("layers")
            .elements[0]
            .field("initial_green")
            .field("expression");
    if (consumer == 1)
        return source.field("flow").field("uniform").field("expression");
    return source.field("analysis").field("loads").elements[0].field("expression");
}
io::schema::Expression expression_record(const io::SiCase& value, unsigned consumer) {
    const auto root = case_root(value);
    if (consumer == 0) {
        const auto shell = take(root.model.record().shells.at(0)).record();
        const auto layer = take(shell.layers.at(0)).record();
        return std::get<io::schema::ReferenceFieldExpression>(layer.initial_green.record())
            .expression.record();
    }
    if (consumer == 1) {
        const auto flow = std::get<io::schema::FlowKnown>(root.flow.record());
        return std::get<io::schema::TimeSeriesExpression>(flow.uniform.record())
            .expression.record();
    }
    const auto analysis = std::get<io::schema::AnalysisStatic>(root.analysis.record());
    return std::get<io::schema::LoadCurrentExpression>(take(analysis.loads.at(0)).record())
        .expression.record();
}
void polynomial(const io::SiCase& value, const io::schema::VectorPolynomial& p,
                const std::vector<Rat>& z, const Rat& multiplier) {
    REQUIRE(p.terms.size() == z.size());
    for (std::size_t j = 0; j < z.size(); ++j) {
        const auto term = take(p.terms.at(static_cast<core::Index>(j))).record();
        for (core::Index component = 0; component < 3; ++component) {
            const Rat expected = component == 2 ? Rat{z[j] * multiplier} : Rat{0};
            const auto q = take(value.quantity(take(term.value.at(component))));
            contains(q, expected, expected);
            dim(q, io::SchemaDimension{
                       {1, 0, static_cast<std::int16_t>(multiplier == Rat{1} ? -1 : 0), 0}});
        }
    }
}
void point_checks(const io::SiCase& value, const io::schema::PerturbationPhysical& p,
                  const std::vector<Rat>& dz, const std::vector<Rat>& vz, bool rigid) {
    for (unsigned velocity = 0; velocity < 2; ++velocity) {
        const auto list = velocity ? p.velocities : p.displacements;
        const auto& expected_z = velocity ? vz : dz;
        REQUIRE(list.size() == expected_z.size());
        for (std::size_t j = 0; j < expected_z.size(); ++j) {
            const auto point = take(list.at(static_cast<core::Index>(j))).record();
            for (core::Index component = 0; component < 3; ++component) {
                const Rat expected = component == 2
                                         ? Rat{expected_z[j] * (velocity ? Rat{1} : Rat{1} / 1000)}
                                         : Rat{0};
                const auto value_q = take(value.quantity(take(point.value.at(component))));
                contains(value_q, expected, expected);
                dim(value_q,
                    io::SchemaDimension{{1, 0, static_cast<std::int16_t>(velocity ? -1 : 0), 0}});
                const Rat tolerance =
                    Rat{1} / boost::multiprecision::pow(Big{10}, velocity ? 12 : 15);
                const auto tolerance_q = take(value.quantity(take(point.tolerance.at(component))));
                contains(tolerance_q, tolerance, tolerance);
                dim(tolerance_q,
                    io::SchemaDimension{{1, 0, static_cast<std::int16_t>(velocity ? -1 : 0), 0}});
            }
            const auto selected = point.point.record();
            if (rigid) {
                const auto coords = std::get<io::schema::InitialPointRigid>(selected).coordinates;
                for (core::Index k = 0; k < 3; ++k) {
                    const Rat expected = k == 0 ? Rat{1} / 1000 : Rat{0};
                    const auto coordinate = take(value.quantity(take(coords.at(k))));
                    contains(coordinate, expected, expected);
                    dim(coordinate, io::SchemaDimension{{1, 0, 0, 0}});
                }
            } else {
                const auto shell = std::get<io::schema::InitialPointShell>(selected);
                CHECK(shell.side_s.text() == "lower");
                CHECK(shell.side_t.text() == "lower");
                contains(take(value.quantity(take(shell.coordinates.at(0)))), Rat{1}, Rat{1});
                contains(take(value.quantity(take(shell.coordinates.at(1)))), Rat{1} / 2,
                         Rat{1} / 2);
                const Rat zeta = expected_z.size() == 2 ? (j == 0 ? Rat{1} : Rat{-1}) : Rat{0};
                const auto coordinate = take(value.quantity(take(shell.coordinates.at(2))));
                contains(coordinate, zeta, zeta);
                dim(coordinate, io::SchemaDimension{{0, 0, 0, 0}});
            }
        }
    }
}
} // namespace
TEST_CASE("SI viscosity and damping convert every nonzero mixed matrix component",
          "[IO4-01][verification][property]") {
    auto source = beam_case();
    auto& units = source.field("units");
    units.set("length", string("mm"));
    units.set("mass", string("g"));
    units.set("time", string("ms"));
    units.set("angle", string("deg"));
    auto& viscosity = source.field("model").field("beams").elements[0].field("viscosity");
    auto& damping = source.field("model").field("supports").elements[1].field("damping");
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) {
            viscosity.elements[i].elements[j] = parse("2");
            damping.elements[i].elements[j] = parse("3");
        }
    const auto value = converted(source);
    const auto model = case_root(value).model.record();
    const auto beam = take(model.beams.at(0)).record();
    const auto support = std::get<io::schema::SupportElastic>(take(model.supports.at(1)).record());
    const auto pi = pi_bounds();
    for (core::Index i = 0; i < 6; ++i)
        for (core::Index j = 0; j < 6; ++j) {
            const auto rotational = static_cast<unsigned>(i >= 3) + static_cast<unsigned>(j >= 3);
            const Rat expected = Rat{2} / boost::multiprecision::pow(Big{1000}, rotational + 1);
            const auto v = take(value.quantity(take(take(beam.viscosity.at(i)).at(j))));
            contains(v, expected, expected);
            CHECK(v.recipe.pi_power == 0);
            const auto d = take(value.quantity(take(take(support.damping.at(i)).at(j))));
            CHECK(d.recipe.pi_power == -static_cast<int>(rotational));
            if (rotational == 0)
                contains(d, Rat{3}, Rat{3});
            else if (rotational == 1)
                contains(d, Rat{27} / (50 * pi.second), Rat{27} / (50 * pi.first));
            else
                contains(d, Rat{243} / (2500 * pi.second * pi.second),
                         Rat{243} / (2500 * pi.first * pi.first));
        }
}
TEST_CASE("SI expression constants and domains convert once while the DAG remains unevaluated",
          "[IO4-01][verification][negative][property]") {
    for (unsigned consumer = 0; consumer < 3; ++consumer) {
        INFO(consumer);
        const auto source = expression_source(consumer);
        const auto value = converted(source);
        const auto expression = expression_record(value, consumer);
        REQUIRE(expression.nodes.size() == 3);
        const auto variable =
            std::get<io::schema::ExprNodeVariable>(take(expression.nodes.at(0)).record());
        CHECK(variable.axis.unsigned_value() == 0);
        const auto constant =
            std::get<io::schema::ExprNodeConstant>(take(expression.nodes.at(1)).record());
        CHECK(constant.value.number().digits == "37");
        const Rat expected =
            consumer == 0 ? Rat{37} : (consumer == 1 ? Rat{37000} : Rat{37000000000LL});
        contains(take(value.quantity(constant.value)), expected, expected);
        const auto multiply =
            std::get<io::schema::ExprNodeBinary>(take(expression.nodes.at(2)).record());
        CHECK(multiply.op.text() == "multiply");
        CHECK(multiply.first.unsigned_value() == 0);
        CHECK(multiply.second.unsigned_value() == 1);
        for (core::Index j = 0; j < static_cast<core::Index>(expression.outputs.size()); ++j)
            CHECK(take(expression.outputs.at(j)).unsigned_value() == 2);
        const auto domain = expression.domain.record();
        if (consumer == 0) {
            const auto d = std::get<io::schema::ExpressionDomainReference>(domain);
            CHECK(d.coordinates.text() == "chart");
            contains(take(value.quantity(take(d.upper.at(0)))), Rat{1}, Rat{1});
        } else if (consumer == 1) {
            const auto d = std::get<io::schema::ExpressionDomainTime>(domain);
            contains(take(value.quantity(d.upper)), Rat{1} / 1000, Rat{1} / 1000);
            contains(take(value.quantity(take(d.events.at(1)))), Rat{1} / 1000, Rat{1} / 1000);
        } else {
            const auto d = std::get<io::schema::ExpressionDomainCurrentPositionTime>(domain);
            contains(take(value.quantity(take(d.position_lower.at(1)))), Rat{-2} / 1000,
                     Rat{-2} / 1000);
            contains(take(value.quantity(take(d.position_upper.at(2)))), Rat{6} / 1000,
                     Rat{6} / 1000);
            contains(take(value.quantity(d.time_upper)), Rat{1} / 1000, Rat{1} / 1000);
        }
        auto bad = source;
        expression_json(bad, consumer)
            .field("nodes")
            .elements[1]
            .set("dimensions", parse(consumer == 0 ? "[1,0,0,0]" : "[0,0,0,0]"));
        expect_schema_reject(bad);
    }
}
TEST_CASE(
    "All four full physical initial fields preserve coefficients states and point checks in SI",
    "[IO4-01][verification][property]") {
    std::uint64_t cases = 0;
    for (const auto& example : test::case_schema_corpus::physical) {
        INFO(example.name);
        auto source = physical_case(example.context_json, example.perturbation_json);
        source.field("units").set("length", string("mm"));
        source.field("units").set("mass", string("g"));
        source.field("units").set("time", string("ms"));
        const auto value = converted(source);
        const auto root = case_root(value);
        const auto initial_field = root.initial;
        if (!initial_field.has_value()) {
            FAIL_CHECK("Complete physical input lost its initial state");
            continue;
        }
        const auto initial =
            std::get<io::schema::InitialEquilibrium>(initial_field.value().record());
        const auto physical =
            std::get<io::schema::PerturbationPhysical>(initial.perturbation.record());
        if (example.name == "beam_rigid_nonzero") {
            REQUIRE(physical.beams.size() == 1);
            REQUIRE(physical.rigids.size() == 1);
            const auto b = take(physical.beams.at(0)).record();
            const auto rigid = take(physical.rigids.at(0)).record().state.record();
            for (core::Index node = 0; node < 2; ++node) {
                const auto state = take(b.nodes.at(node)).record();
                for (core::Index k = 0; k < 3; ++k) {
                    const Rat position = node == 1 && k == 0 ? Rat{1} / 1000 : Rat{0};
                    contains(take(value.quantity(take(state.position.at(k)))), position, position);
                    const Rat omega = node == 1 && k == 2 ? Rat{1000} : Rat{0};
                    contains(take(value.quantity(take(state.angular_velocity.at(k)))), omega,
                             omega);
                    contains(take(value.quantity(take(state.linear_velocity.at(k)))), Rat{0},
                             Rat{0});
                    const Rat rp = k == 0 ? Rat{1} / 1000 : Rat{0};
                    contains(take(value.quantity(take(rigid.position.at(k)))), rp, rp);
                }
                const std::array<std::array<int, 3>, 3> rotation =
                    node == 0
                        ? std::array<std::array<int, 3>, 3>{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}
                        : std::array<std::array<int, 3>, 3>{{{1, 0, 0}, {0, 0, -1}, {0, 1, 0}}};
                for (core::Index i = 0; i < 3; ++i)
                    for (core::Index j = 0; j < 3; ++j) {
                        const Rat expected{
                            rotation[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]};
                        contains(take(value.quantity(take(take(state.rotation.at(i)).at(j)))),
                                 expected, expected);
                        if (node == 1)
                            contains(take(value.quantity(take(take(rigid.rotation.at(i)).at(j)))),
                                     expected, expected);
                    }
            }
            for (core::Index k = 0; k < 3; ++k) {
                const Rat expected = k == 2 ? Rat{1000} : Rat{0};
                contains(take(value.quantity(take(rigid.angular_velocity.at(k)))), expected,
                         expected);
                contains(take(value.quantity(take(rigid.linear_velocity.at(k)))), Rat{0}, Rat{0});
            }
            point_checks(value, physical, {Rat{0}}, {Rat{1}}, true);
        } else {
            REQUIRE(physical.shells.size() == 1);
            const auto shell = take(physical.shells.at(0)).record();
            for (unsigned velocity = 0; velocity < 2; ++velocity) {
                const auto field = (velocity ? shell.velocity : shell.displacement).record();
                const Rat multiplier = velocity ? Rat{1} : Rat{1} / 1000;
                if (example.name == "flat_half_cells") {
                    const auto cells = std::get<io::schema::FullShellFieldHalfCells>(field).cells;
                    REQUIRE(cells.size() == 4);
                    for (core::Index cell = 0; cell < 4; ++cell) {
                        const auto c = take(cells.at(cell)).record();
                        CHECK(ordinal(take(c.cell.at(0))) == static_cast<std::uint64_t>(cell / 2));
                        CHECK(ordinal(take(c.cell.at(1))) == static_cast<std::uint64_t>(cell % 2));
                        const Rat vscale = velocity ? Rat{2} : Rat{1};
                        const std::vector<Rat> z =
                            cell < 2 ? std::vector<Rat>{vscale / 1600, vscale / 800, vscale / 1600}
                                     : std::vector<Rat>{9 * vscale / 1600, 3 * vscale / 800,
                                                        vscale / 1600};
                        polynomial(value, c.r.record(), z, multiplier);
                        polynomial(value, c.p.record(), {}, multiplier);
                        polynomial(value, c.c.record(), {}, multiplier);
                    }
                } else {
                    const auto global = std::get<io::schema::FullShellFieldGlobalPolynomial>(field);
                    if (example.name == "flat_clamped_global") {
                        polynomial(value, global.r.record(),
                                   {velocity ? Rat{1} / 50 : Rat{1} / 100}, multiplier);
                        polynomial(value, global.p.record(), {}, multiplier);
                        polynomial(value, global.c.record(), {}, multiplier);
                    } else {
                        polynomial(value, global.r.record(), {}, multiplier);
                        polynomial(value, global.p.record(),
                                   {velocity ? Rat{1} / 50 : Rat{1} / 100}, multiplier);
                        polynomial(value, global.c.record(),
                                   {velocity ? Rat{1} / 100 : Rat{1} / 200}, multiplier);
                    }
                }
            }
            if (example.name == "flat_full_thickness")
                point_checks(value, physical, {Rat{9} / 1600, Rat{-7} / 1600},
                             {Rat{9} / 800, Rat{-7} / 800}, false);
            else
                point_checks(value, physical, {Rat{1} / 100}, {Rat{1} / 50}, false);
        }
        ++cases;
    }
    CHECK(cases == 4);
}
TEST_CASE("Frozen physical negative stages reject thirteen schemas and retain nine SI-only states",
          "[IO4-01][verification][negative]") {
    std::uint64_t schema_count = 0, prepare_count = 0;
    for (const auto& example : test::case_schema_corpus::negatives) {
        INFO(example.name);
        const auto source = physical_case(example.context_json, example.perturbation_json);
        if (example.stage == test::case_schema_corpus::ExpectedStage::Schema) {
            expect_schema_reject(source);
            ++schema_count;
        } else {
            const auto expected = checked(source);
            const auto bytes = test::schema_snapshot::serialize(take(expected.root()));
            const auto result = converted(source);
            const auto retained = take(take(result.source_case()).get().root());
            CHECK(test::schema_snapshot::serialize(retained) == bytes);
            ++prepare_count;
        }
    }
    CHECK(schema_count == 13);
    CHECK(prepare_count == 9);
}

TEST_CASE("SI requirement views own all rational and decimal request tokens after caller mutation",
          "[IO4-01][verification][property]") {
    const auto value = [] {
        auto original = root();
        original.field("units").set("length", string("mm"));
        const auto source = checked(original);
        const auto limits = unit_limits();
        core::ResourcePlan plan{budget};
        io::SiDiagnostic diagnostic;
        const auto catalog =
            si_phase(io::make_si_catalog(source, limits, plan, &diagnostic), "catalog", diagnostic);
        std::string numerator = "3", denominator = "2", absolute = "1e-25", relative_n = "1",
                    relative_d = "10000000000000000000000000";
        auto rows = requirements(catalog);
        for (auto& row : rows) {
            row.scale.value = io::SiRationalInput{numerator, denominator};
            row.absolute = io::SiDecimalInput{absolute};
            row.relative = io::SiRationalInput{relative_n, relative_d};
        }
        const auto request = si_phase(io::make_si_request(catalog, rows, limits, plan, &diagnostic),
                                      "request", diagnostic);
        numerator = "7";
        denominator = "5";
        absolute = "9e-20";
        relative_n = "9";
        relative_d = "1";
        return si_phase(io::convert_case_to_si(source, catalog, request, limits, plan, &diagnostic),
                        "convert", diagnostic);
    }();
    const auto geometry = std::get<io::schema::GeometryPg>(
        take(case_root(value).model.record().geometry.at(0)).record());
    const auto q = take(value.quantity(geometry.reference_length));
    const auto scale = std::get<io::SiRationalInput>(q.requirement.scale.value);
    CHECK(scale.numerator == "3");
    CHECK(scale.denominator == "2");
    CHECK(std::get<io::SiDecimalInput>(q.requirement.absolute).token == "1e-25");
    const auto relative = std::get<io::SiRationalInput>(q.requirement.relative);
    CHECK(relative.numerator == "1");
    CHECK(relative.denominator == "10000000000000000000000000");
    const io::SchemaDimension length{{1, 0, 0, 0}};
    CHECK(q.requirement.scale.dimensions == length);
}
