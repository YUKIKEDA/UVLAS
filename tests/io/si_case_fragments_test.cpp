#include "case_schema_expected_fields.hpp"
#include "case_schema_payload.hpp"
#include "si_case_detail.hpp"
#include "si_case_fragment_context.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <set>
#include <stdexcept>
#include <string>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
using namespace uvlas::test::si_fragment_adapter;
namespace si = io::detail::si;
using Big = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
constexpr std::uint64_t budget = std::uint64_t{1} << 40;
template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw std::runtime_error("Scoped SI control failed, code=" +
                                 std::to_string(static_cast<unsigned>(result.error().code)));
    return std::move(*result);
}
Rat decimal(std::string_view token) {
    bool negative = false;
    if (token.front() == '-' || token.front() == '+') {
        negative = token.front() == '-';
        token.remove_prefix(1);
    }
    Big coefficient = 0;
    std::int64_t fractional = 0;
    bool dot = false;
    const auto e = token.find_first_of("eE");
    const auto mantissa = token.substr(0, e);
    for (const char c : mantissa) {
        if (c == '.') {
            dot = true;
            continue;
        }
        REQUIRE(c >= '0');
        REQUIRE(c <= '9');
        coefficient = 10 * coefficient + (c - '0');
        if (dot)
            ++fractional;
    }
    const auto exponent =
        e == std::string_view::npos ? 0 : std::stoll(std::string{token.substr(e + 1)});
    const auto power = exponent - fractional;
    REQUIRE(power > -10000);
    REQUIRE(power < 10000);
    if (coefficient == 0)
        return 0;
    Rat value{negative ? -coefficient : coefficient};
    if (power >= 0)
        value *= boost::multiprecision::pow(Big{10}, static_cast<unsigned>(power));
    else
        value /= boost::multiprecision::pow(Big{10}, static_cast<unsigned>(-power));
    return value;
}
Rat word(double value) {
    const auto raw = std::bit_cast<std::uint64_t>(value);
    const auto exponent = (raw >> 52) & 0x7ffU;
    REQUIRE(exponent != 0x7ffU);
    Big significand = raw & 0xfffffffffffffULL;
    if (exponent)
        significand += Big{1} << 52;
    if (significand == 0)
        return 0;
    Rat result{significand};
    const auto power = exponent ? static_cast<int>(exponent) - 1075 : -1074;
    if (power >= 0)
        result *= Big{1} << static_cast<unsigned>(power);
    else
        result /= Big{1} << static_cast<unsigned>(-power);
    return (raw >> 63) ? -result : result;
}
Rat power(Rat base, std::int16_t exponent) {
    auto remaining = static_cast<unsigned>(exponent < 0 ? -static_cast<int>(exponent) : exponent);
    if (exponent < 0)
        base = 1 / base;
    Rat value = 1;
    while (remaining) {
        if (remaining & 1U)
            value *= base;
        remaining >>= 1;
        if (remaining)
            base *= base;
    }
    return value;
}
std::pair<Rat, Rat> pi_bounds() {
    const auto atan = [](unsigned q) {
        Rat sum = 0;
        for (unsigned j = 0; j < 160; ++j) {
            const Rat term = Rat{1} / ((2 * j + 1) * boost::multiprecision::pow(Big{q}, 2 * j + 1));
            sum += (j & 1) ? -term : term;
        }
        const Rat next = Rat{1} / (321 * boost::multiprecision::pow(Big{q}, 321));
        return std::pair<Rat, Rat>{sum, sum + next};
    };
    const auto a = atan(5), b = atan(239);
    return {16 * a.first - 4 * b.second, 16 * a.second - 4 * b.first};
}
std::optional<std::uint64_t> member(const io::detail::InputDocumentData& source, std::uint64_t node,
                                    std::string_view name) {
    const auto& object = std::get<io::detail::ObjectValue>(
        source.nodes.values()[static_cast<std::size_t>(node)].value);
    auto index = object.members.head;
    while (index >= 0) {
        const auto& item = source.members.values()[static_cast<std::size_t>(index)];
        if (std::string_view{source.text.values().data() + item.key.offset, item.key.size} == name)
            return static_cast<std::uint64_t>(item.child);
        index = item.next;
    }
    return std::nullopt;
}
Rat original(const io::detail::InputDocumentData& source, std::uint64_t node) {
    const auto& value = source.nodes.values()[static_cast<std::size_t>(node)].value;
    if (const auto* integer = std::get_if<std::int64_t>(&value))
        return Rat{*integer};
    if (const auto* integer = std::get_if<std::uint64_t>(&value))
        return Rat{*integer};
    if (std::holds_alternative<io::detail::DecimalValue>(value))
        return decimal(io::detail::node_text(source, node));
    const auto numerator = member(source, node, "numerator"),
               denominator = member(source, node, "denominator");
    if (!numerator || !denominator)
        throw std::runtime_error("A qualified rational must retain both original integer tokens");
    return decimal(io::detail::node_text(source, *numerator)) /
           decimal(io::detail::node_text(source, *denominator));
}
bool descendant(std::span<const io::detail::SchemaNode> nodes, std::uint64_t root,
                std::uint64_t node) {
    while (node != root) {
        const auto parent = nodes[static_cast<std::size_t>(node)].parent;
        if (!parent)
            return false;
        node = *parent;
    }
    return true;
}
struct Coverage {
    std::set<std::string> records, fields;
    std::uint64_t numbers = 0, native = 0;
};
void exercise(const FragmentInput& input, bool accepted, Coverage& coverage) {
    auto wrapped = input.context;
    if (input.type == "Case")
        wrapped = input.payload;
    else
        wrapped.set("units", input.payload);
    const auto wire = wrapped.wire();
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    core::ResourcePlan source_plan{budget};
    std::optional<io::InputDiagnostic> parser_diagnostic;
    const auto document =
        take(io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                     io::SourceId{371}, syntax, source_plan, parser_diagnostic));
    const auto view = take(document.view());
    const auto object = view.root.object();
    if (!object)
        throw std::runtime_error("Scoped context root is an object");
    const auto payload = input.type == "Case" ? std::optional{view.root} : object->find("units");
    if (!payload)
        throw std::runtime_error("Scoped payload exists");
    const auto types = io::detail::schema_types();
    const auto selected =
        std::ranges::find_if(types, [&](const auto& type) { return type.spelling == input.type; });
    REQUIRE(selected != types.end());
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(1000000);
    schema_caps[static_cast<std::size_t>(io::SchemaQuota::TypedBytes)] = budget;
    schema_caps[static_cast<std::size_t>(io::SchemaQuota::ScratchBytes)] = budget;
    schema_caps[static_cast<std::size_t>(io::SchemaQuota::Work)] = budget;
    schema_caps[static_cast<std::size_t>(io::SchemaQuota::ExactDigits)] = 4096;
    io::SchemaDiagnostic schema_diagnostic;
    const auto checked = io::detail::validate_schema_payload(
        document,
        {*payload, static_cast<io::detail::FieldType>(selected - types.begin()), input.components,
         input.dimensions, input.type != "Case"},
        io::SchemaLimits{schema_caps}, source_plan, &schema_diagnostic);
    CHECK(checked.has_value() == accepted);
    if (!checked) {
        if (accepted)
            FAIL("Required scoped schema control was rejected");
        CHECK(checked.error().code == core::FailureCode::InvalidInput);
        return;
    }
    const auto& data = checked->data.values()[0];
    const auto& source = data.snapshot.values()[0];
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> caps{};
    caps.fill(budget);
    const io::SiLimits limits{caps,
                              {100000, 100000, 100000, 64, 64, 1000000, 2000000,
                               structure::ExactEstimateAdmission{budget},
                               structure::ExactTimingMode::Disabled}};
    core::ResourcePlan plan{budget};
    io::SiDiagnostic diagnostic;
    io::detail::SiContext context{limits, plan, &diagnostic};
    context.program = context.allocate<si::ProgramScratch>(1, true);
    // A scoped numeric test is not a public Case/SiCase qualification. Its caller
    // supplies an explicit non-SI unit recipe to the same production scalar kernel.
    const io::SiUnits units{io::SiLengthUnit::Inch, io::SiMassUnit::PoundMass,
                            io::SiTimeUnit::Minute, io::SiAngleUnit::Degree};
    const auto pi = si::programmed_pi(context.program.values()[0], context.work);
    static const auto pi_reference = pi_bounds();
    for (std::uint64_t node = 0; node < data.nodes.size(); ++node) {
        const auto& typed = data.nodes.values()[static_cast<std::size_t>(node)];
        if (!typed.assigned || !descendant(data.nodes.values(), checked->root, node))
            continue;
        const auto kind = io::detail::field_kinds[static_cast<std::size_t>(typed.type)];
        const auto spelling =
            io::detail::schema_types()[static_cast<std::size_t>(typed.type)].spelling;
        const bool numeric_spelling = spelling == "Number" || spelling == "Rational" ||
                                      spelling == "ExactQuantity" || spelling.starts_with("Q<");
        CHECK(io::detail::si_numeric_node(typed) == numeric_spelling);
        if (io::detail::node_kind(source, node) == io::InputNodeKind::Object) {
            const auto& record = io::detail::schema_records()[typed.record_kind];
            coverage.records.emplace(record.name);
            for (const auto& field : record.fields)
                if (member(source, node, field.name))
                    coverage.fields.insert(std::string{record.name} + "/" +
                                           std::string{field.name});
        }
        if (kind == io::detail::FieldKind::Integer || kind == io::detail::FieldKind::Identifier ||
            kind == io::detail::FieldKind::Reference) {
            CHECK_FALSE(io::detail::si_numeric_node(typed));
            ++coverage.native;
            continue;
        }
        if (!io::detail::si_numeric_node(typed))
            continue;
        ++coverage.numbers;
        const auto recipe =
            io::detail::si_recipe(source, data.nodes.values(), data.digits.values(), node, units);
        const io::SiRequirementView requirement{{io::SiDecimalInput{"1"}, recipe.dimensions},
                                                io::SiDecimalInput{"1e-25"},
                                                io::SiDecimalInput{"1e-25"}};
        const auto actual = io::detail::si_evaluate_quantity(context, recipe, requirement, pi);
        const auto& d = recipe.dimensions.exponents;
        const Rat coefficient = original(source, node) * power(Rat{127} / 5000, d[0]) *
                                power(Rat{45359237} / 100000000, d[1]) * power(Rat{60}, d[2]) *
                                power(Rat{1} / 180, d[3]);
        Rat lower = coefficient * power(d[3] < 0 ? pi_reference.second : pi_reference.first, d[3]);
        Rat upper = coefficient * power(d[3] < 0 ? pi_reference.first : pi_reference.second, d[3]);
        if (lower > upper)
            std::swap(lower, upper);
        std::visit(
            [&](const auto& x) {
                using T = std::decay_t<decltype(x)>;
                Rat center;
                if constexpr (std::is_same_v<T, io::ScaledBinary64>)
                    center = word(x.center);
                else
                    center = word(x.high) + word(x.low);
                REQUIRE(x.exponent2 > -10000);
                REQUIRE(x.exponent2 < 10000);
                Rat factor = 1;
                if (x.exponent2 >= 0)
                    factor *= Big{1} << static_cast<unsigned>(x.exponent2);
                else
                    factor /= Big{1} << static_cast<unsigned>(-x.exponent2);
                const Rat radius = word(x.radius) * factor;
                center *= factor;
                CHECK(center - radius <= lower);
                CHECK(center + radius >= upper);
                const Rat magnitude =
                    lower <= 0 && upper >= 0 ? Rat{0} : Rat{std::min(abs(lower), abs(upper))};
                CHECK(radius <= (1 + magnitude) / boost::multiprecision::pow(Big{10}, 25));
            },
            actual);
    }
}
} // namespace
TEST_CASE("Frozen dictionary numeric leaves share the production SI scalar kernel",
          "[IO4-01][verification][property]") {
    Coverage coverage;
    std::uint64_t visited = 0, accepted = 0, rejected = 0;
    for (const auto& example : test::case_schema_corpus::fragments) {
        INFO(example.name);
        auto input = input_for(example.name, parse(example.json));
        const bool normal = !example.superseded && example.name != "Geometry.gmsh";
        exercise(input, normal, coverage);
        ++visited;
        if (normal)
            ++accepted;
        else
            ++rejected;
        if (example.name == "Geometry.gmsh") {
            input.payload.field("reference").field("file").set("source_id", string("3"));
            exercise(input, true, coverage);
        }
    }
    CHECK(visited == 144);
    CHECK(accepted == 142);
    CHECK(rejected == 2);
    for (const auto& physical : test::case_schema_corpus::physical) {
        INFO(physical.name);
        auto value = physical_case(physical.context_json, physical.perturbation_json);
        const auto context = parse(physical.context_json);
        if (!context.field("beam_ids").elements.empty()) {
            auto check = parse(
                R"({"point":{"target":{"kind":"beam","id":"1"},"coordinates":[0,0,0]},"value":[0,0,0],"tolerance":[0,0,0]})");
            check.field("point").field("target").set("id", context.field("beam_ids").elements[0]);
            value.field("initial")
                .field("perturbation")
                .field("displacements")
                .elements.push_back(std::move(check));
        }
        exercise({root(), std::move(value), "Case", 1, {{0, 0, 0, 0}}}, true, coverage);
    }
    exercise(
        input_for(
            "SurfaceWake",
            parse(
                R"({"surface":"1","rows":1,"points":[[0,0,0],[1,0,0]],"coefficients":[0],"row_ids":["1"],"point_ids":["1","2"],"band_ids":["1"],"birth_times":[0]})")),
        true, coverage);
    exercise(input_for("EigenCluster", parse(R"({"members":[0],"expected_dimension":1})")), true,
             coverage);
    auto value = root("coupled_flutter");
    value.set("restart", fragment("Restart"));
    value.set("provenance",
              parse(R"({"description":"explicit SI control","references":[],"assumptions":[]})"));
    auto& graph = value.field("model").field("geometry").elements[0];
    graph.set("degree", parse("[1,0]"));
    graph.set("coefficients",
              parse(R"([{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}])"));
    auto other = graph;
    other.set("id", string("2"));
    value.field("model").field("geometry").elements.push_back(std::move(other));
    value.field("model")
        .field("geometry")
        .elements[0]
        .set(
            "continuity",
            parse(
                R"([{"other":"2","kind":"smooth","own_boundary":0,"other_boundary":0,"orientation":"same"}])"));
    auto& surfaces = value.field("model").field("aero_surfaces").elements;
    auto second = surfaces[0];
    second.set("id", string("2"));
    second.set("geometry", string("2"));
    surfaces.push_back(std::move(second));
    surfaces[0].set(
        "sharing",
        parse(R"([{"other":"2","own_boundary":0,"other_boundary":0,"orientation":"same"}])"));
    exercise({root(), std::move(value), "Case", 1, {{0, 0, 0, 0}}}, true, coverage);
    for (const auto& expected : test::schema_field_inventory::fields) {
        INFO(expected.record << "/" << expected.name);
        CHECK(coverage.records.contains(std::string{expected.record}));
        CHECK(coverage.fields.contains(std::string{expected.record} + "/" +
                                       std::string{expected.name}));
    }
    CHECK(coverage.records.size() == 191);
    CHECK(coverage.fields.size() == 806);
    CHECK(coverage.numbers > 0);
    CHECK(coverage.native > 0);
}

TEST_CASE("SI exact fallback retains sparse low words and refuses unrepresentable tails",
          "[IO4-01][property][negative]") {
    for (const unsigned tail : {1000U, 1074U, 1100U}) {
        INFO(tail);
        std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> caps{};
        caps.fill(budget);
        const io::SiLimits limits{caps,
                                  {100000, 100000, 100000, 64, 64, 1000000, 2000000,
                                   structure::ExactEstimateAdmission{budget},
                                   structure::ExactTimingMode::Disabled}};
        core::ResourcePlan plan{budget};
        io::detail::SiContext context{limits, plan, nullptr};
        context.program = context.allocate<si::ProgramScratch>(1, true);
        // Independent finite decimal: (1 + 2^-k) = (2^k+1)*5^k*10^-k.
        const Big coefficient = ((Big{1} << tail) + 1) * boost::multiprecision::pow(Big{5}, tail);
        const auto digits = coefficient.convert_to<std::string>();
        const io::SchemaDimension dimension{{0, 0, 0, 0}};
        const io::SiRecipeView recipe{
            io::SiDecimalView{1, digits, -static_cast<std::int64_t>(tail)},
            {},
            {io::SiLengthUnit::Metre, io::SiMassUnit::Kilogram, io::SiTimeUnit::Second,
             io::SiAngleUnit::Radian},
            dimension,
            0};
        const io::SiRequirementView wanted{
            {io::SiDecimalInput{"1"}, dimension}, io::SiDecimalInput{"0"}, io::SiDecimalInput{"0"}};
        if (tail <= 1074U) {
            const auto output =
                io::detail::si_evaluate_quantity(context, recipe, wanted, si::integer(0));
            const auto* paired = std::get_if<io::ScaledTwoWord>(&output);
            REQUIRE(paired != nullptr);
            CHECK(paired->exponent2 == 0);
            CHECK(paired->radius == 0);
            CHECK(word(paired->high) + word(paired->low) == Rat{1} + Rat{1} / (Big{1} << tail));
        } else {
            try {
                (void)io::detail::si_evaluate_quantity(context, recipe, wanted, si::integer(0));
                FAIL("An exact zero budget must not discard an unrepresentable sparse tail");
            } catch (const core::Failure& failure) {
                CHECK(failure.code == core::FailureCode::PrecisionUnresolved);
            }
        }
    }
}
