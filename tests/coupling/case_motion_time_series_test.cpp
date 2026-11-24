#include "coupling/case_motion_time_series.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>

namespace {
using namespace uvlas;
namespace c = coupling::detail;
namespace p = aero::detail::pchip;
namespace a = uvlas::test::schema_adapter;
using enum p::JetOrder;
using enum p::RequestedSide;

constexpr std::uint64_t budget = std::uint64_t{1} << 48U;
constexpr p::PchipInputStamp stamp{p::PchipIssuerId{71}, p::PchipInputGeneration{8}};

structure::ExactLimits exact_limits() noexcept {
    return {100000,
            100000,
            100000,
            64,
            64,
            100000,
            200000,
            structure::ExactEstimateAdmission{budget},
            structure::ExactTimingMode::Disabled};
}

const c::CaseMotionTimeLimits limits{
    exact_limits(),
    {{128, 100000, budget, 10000}, numerics::linear::EstimateAdmission{budget}},
    100,
    1000};

template <class T>
T take(core::Result<T> result) {
    if (!result)
        FAIL("Unexpected Motion failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}

template <class T>
void rejected(const core::Result<T>& result, core::FailureCode code) {
    if (result)
        FAIL("Unexpected Motion success");
    CHECK(result.error().code == code);
}

io::SiCase checked_case(const a::JsonToken& input) {
    core::ResourcePlan plan{budget};
    const auto wire = input.wire();
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    std::optional<io::InputDiagnostic> parsing;
    const auto doc =
        take(io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                     io::SourceId{19}, syntax, plan, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(budget);
    io::SchemaDiagnostic schema_diagnostic;
    auto checked =
        io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, plan, &schema_diagnostic);
    if (!checked)
        FAIL("Motion schema fixture at " << schema_diagnostic.path().value_or(""));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(budget);
    const io::SiLimits bounds{si_caps, exact_limits()};
    const auto catalog = take(io::make_si_catalog(*checked, bounds, plan, nullptr));
    std::vector<io::SiRequirementRow> requirements;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto entry = take(catalog.entry(i));
        requirements.push_back({entry.key,
                                {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                                io::SiDecimalInput{"1e-25"},
                                io::SiDecimalInput{"1e-25"}});
    }
    const auto request = take(io::make_si_request(catalog, requirements, bounds, plan, nullptr));
    return take(io::convert_case_to_si(*checked, catalog, request, bounds, plan, nullptr));
}

c::CaseMotionField motion_field(const io::SiCase& source) {
    const auto root = take(take(source.source_case()).get().root());
    const auto model = root.model.record();
    const auto object = take(model.aero_surfaces.at(0)).record().attachment.record();
    if (const auto* motion = std::get_if<io::schema::AeroAttachmentPrescribed>(&object))
        return motion->motion;
    FAIL("Missing prescribed attachment");
}

a::JsonToken& motion_json(a::JsonToken& root) {
    return root.field("model")
        .field("aero_surfaces")
        .elements[0]
        .field("attachment")
        .field("motion");
}

structure::ExactRational exact_time(std::string_view decimal, core::ResourcePlan& plan) {
    return take(structure::ExactRational::from_decimal(decimal, exact_limits(), plan));
}

void exact(core::Result<p::RationalView> result, std::uint32_t n, std::uint32_t d = 1,
           int sign = 1) {
    const auto view = take(std::move(result));
    const std::array numerator{n};
    const std::array denominator{d};
    CHECK(view.sign == (n == 0 ? 0 : sign));
    const auto expected = n == 0 ? std::span<const std::uint32_t>{} : std::span{numerator}.first(1);
    CHECK(std::ranges::equal(view.numerator, expected));
    CHECK(std::ranges::equal(view.denominator, denominator));
}
} // namespace

TEST_CASE("Typed Motion owns polynomial jets and a common exact coefficient of pi",
          "[AE-07][IO4-01][N4-02][N4-03][property][negative]") {
    auto input = a::root("aerodynamic_prescribed_transient");
    input.field("units").set("time", a::string("ms"));
    input.field("units").set("angle", a::string("deg"));
    auto& motion = motion_json(input);
    motion.set(
        "translation",
        a::parse(
            R"({"kind":"polynomial","origin":1000,"scale":2000,"domain":[0,3000],"coefficients":[[3,0,0],[4,0,0],[5,0,0]]})"));
    motion.set("rotation_vector", a::parse(R"({"kind":"constant","value":[90,-45,0]})"));
    motion.set("events", a::parse("[]"));

    p::CurveWorkspace workspace;
    const auto field = [&] {
        const auto si = checked_case(input);
        core::ResourcePlan plan{budget};
        const auto begin = exact_time("0", plan);
        const auto end = exact_time("3", plan);
        return take(c::prepare_case_motion_time_series(si, motion_field(si), {begin, end}, stamp,
                                                       workspace, limits, plan));
    }();
    // A separately converted case with the same source identity can serve as the current owner.
    const auto si = checked_case(input);
    core::ResourcePlan plan{budget};
    const auto zero = exact_time("0", plan);
    const auto end = exact_time("3", plan);
    const auto current = take(c::prepare_case_motion_time_series(si, motion_field(si), {zero, end},
                                                                 stamp, workspace, limits, plan));
    const auto t = exact_time("1.5", plan);
    const auto jet = take(c::evaluate_case_motion_time_series(
        field, current, {stamp, t, Second, Automatic}, workspace, limits, plan));
    CHECK(field.events().empty());
    CHECK(jet.rotation_pi_power() == 1);
    exact(jet.translation(0, Value), 69, 16);
    exact(jet.translation(0, First), 13, 4);
    exact(jet.translation(0, Second), 5, 2);
    exact(jet.rotation_coefficient(0, Value), 1, 2);
    exact(jet.rotation_coefficient(1, Value), 1, 4, -1);
    exact(jet.rotation_coefficient(0, Second), 0);
    CHECK_FALSE(jet.translation(3, Value));
    CHECK(take(jet.translation_location()).selected == p::JetCandidate::Only);
    const auto value = take(c::evaluate_case_motion_time_series(
        field, current, {stamp, t, Value, Automatic}, workspace, limits, plan));
    CHECK_FALSE(value.translation(0, First));

    const auto outside = exact_time("3.01", plan);
    rejected(c::evaluate_case_motion_time_series(field, current, {stamp, outside}, workspace,
                                                 limits, plan),
             core::FailureCode::InvalidInput);
    auto other = stamp;
    other.generation = p::PchipInputGeneration{9};
    rejected(
        c::evaluate_case_motion_time_series(field, current, {other, t}, workspace, limits, plan),
        core::FailureCode::StaleGeneration);
    auto capped = limits;
    capped.maximum_exact_numbers = 15;
    rejected(c::prepare_case_motion_time_series(si, motion_field(si), {zero, end}, stamp, workspace,
                                                capped, plan),
             core::FailureCode::ResourceExceeded);
}

TEST_CASE("Typed Motion PCHIP preserves both accelerations and a shared physical query",
          "[AE-07][IO4-01][N4-02][N4-03][property][negative]") {
    auto input = a::root("aerodynamic_prescribed_transient");
    input.field("units").set("angle", a::string("deg"));
    input.field("analysis")
        .set(
            "intervals",
            a::parse(
                R"([{"start":0,"end":1,"step":1,"restart":"initial"},{"start":1,"end":2,"step":1,"restart":"event"}])"));
    auto& motion = motion_json(input);
    motion.set(
        "translation",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[1,0,0],[3,0,0]],"outside":{"kind":"reject"}})"));
    motion.set(
        "rotation_vector",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[90,0,0],[270,0,0]],"outside":{"kind":"reject"}})"));
    motion.set("events", a::parse("[0,1,2]"));
    const auto si = checked_case(input);
    core::ResourcePlan plan{budget};
    p::CurveWorkspace workspace;
    const auto begin = exact_time("0", plan);
    const auto end = exact_time("2", plan);
    const auto at = exact_time("1", plan);
    const auto field = take(c::prepare_case_motion_time_series(si, motion_field(si), {begin, end},
                                                               stamp, workspace, limits, plan));
    REQUIRE(field.events().size() == 3);
    const auto evaluate = [&](p::JetOrder order, p::RequestedSide side) {
        return c::evaluate_case_motion_time_series(field, field, {stamp, at, order, side},
                                                   workspace, limits, plan);
    };
    const auto common = take(evaluate(First, Automatic));
    exact(common.translation(0, Value), 1);
    exact(common.translation(0, First), 4, 3);
    exact(common.rotation_coefficient(0, First), 2, 3);
    rejected(evaluate(Second, Automatic), core::FailureCode::InvalidInput);
    const auto left = take(evaluate(Second, Left));
    const auto right = take(evaluate(Second, Right));
    exact(left.translation(0, Second), 1, 3);
    exact(right.translation(0, Second), 5, 3);
    exact(left.rotation_coefficient(0, Second), 1, 6);
    exact(right.rotation_coefficient(0, Second), 5, 6);
    CHECK(take(left.translation_location()).selected == p::JetCandidate::Left);
    CHECK(take(right.rotation_location()).selected == p::JetCandidate::Right);
    CHECK(left.rotation_pi_power() == 1);
    CHECK(right.rotation_pi_power() == 1);
}
