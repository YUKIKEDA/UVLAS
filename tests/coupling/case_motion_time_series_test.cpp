#include "coupling/case_flow_clock.hpp"
#include "coupling/case_motion.hpp"
#include "coupling/case_motion_time_series.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>

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

io::schema::AeroAttachmentPrescribed attachment_field(const io::SiCase& source) {
    const auto root = take(take(source.source_case()).get().root());
    const auto model = root.model.record();
    const auto object = take(model.aero_surfaces.at(0)).record().attachment.record();
    if (const auto* motion = std::get_if<io::schema::AeroAttachmentPrescribed>(&object))
        return *motion;
    FAIL("Missing prescribed attachment");
}

c::CaseMotionField motion_field(const io::SiCase& source) {
    return attachment_field(source).motion;
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

void same_enclosure(core::Result<io::SiEnclosure> a, core::Result<io::SiEnclosure> b) {
    const auto left = std::get<io::ScaledTwoWord>(take(std::move(a)));
    const auto right = std::get<io::ScaledTwoWord>(take(std::move(b)));
    CHECK(std::bit_cast<std::uint64_t>(left.high) == std::bit_cast<std::uint64_t>(right.high));
    CHECK(std::bit_cast<std::uint64_t>(left.low) == std::bit_cast<std::uint64_t>(right.low));
    CHECK(std::bit_cast<std::uint64_t>(left.radius) == std::bit_cast<std::uint64_t>(right.radius));
    CHECK(left.exponent2 == right.exponent2);
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

    capped.maximum_exact_numbers = 16;
    CHECK(c::prepare_case_motion_time_series(si, motion_field(si), {zero, end}, stamp, workspace,
                                             capped, plan));

    auto radians_input = input;
    radians_input.field("units").set("angle", a::string("rad"));
    const auto radians_case = checked_case(radians_input);
    const auto radians = take(c::prepare_case_motion_time_series(
        radians_case, motion_field(radians_case), {zero, end}, stamp, workspace, limits, plan));
    const auto radians_jet = take(c::evaluate_case_motion_time_series(
        radians, radians, {stamp, t, Value}, workspace, limits, plan));
    CHECK(radians_jet.rotation_pi_power() == 0);
    exact(radians_jet.rotation_coefficient(0, Value), 90);
    exact(radians_jet.rotation_coefficient(1, Value), 45, 1, -1);
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

TEST_CASE("Typed Motion retains Hold and periodic knot accelerations in SI time",
          "[AE-07][IO4-01][N4-02][N4-03][property][negative]") {
    auto input = a::root("aerodynamic_prescribed_transient");
    input.field("units").set("angle", a::string("deg"));
    input.field("analysis")
        .set("intervals", a::parse(R"([{"start":0,"end":1,"step":1,"restart":"initial"}])"));
    auto& motion = motion_json(input);
    motion.set(
        "translation",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[1,0,0],[5,0,0]],"outside":{"kind":"hold"}})"));
    motion.set(
        "rotation_vector",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2,3],"values":[[0,0,0],[90,0,0],[-90,0,0],[0,0,0]],"outside":{"kind":"periodic","period":3}})"));
    motion.set("events", a::parse("[0,1,2,3]"));
    const auto si = checked_case(input);
    core::ResourcePlan plan{budget};
    p::CurveWorkspace workspace;
    const auto begin = exact_time("-3", plan);
    const auto end = exact_time("1", plan);
    const auto zero = exact_time("0", plan);
    const auto field = take(c::prepare_case_motion_time_series(si, motion_field(si), {begin, end},
                                                               stamp, workspace, limits, plan));
    REQUIRE(field.events().size() == 4);
    const auto evaluate = [&](const structure::ExactRational& time, p::JetOrder order,
                              p::RequestedSide side) {
        return c::evaluate_case_motion_time_series(field, field, {stamp, time, order, side},
                                                   workspace, limits, plan);
    };
    const auto before = take(evaluate(zero, Second, Left));
    const auto after = take(evaluate(zero, Second, Right));
    exact(before.translation(0, Second), 0);
    exact(after.translation(0, Second), 14, 5);
    rejected(evaluate(zero, Second, Automatic), core::FailureCode::InvalidInput);
    const auto repeated_left = take(evaluate(begin, Second, Left));
    const auto repeated_right = take(evaluate(begin, Second, Right));
    exact(repeated_left.rotation_coefficient(0, Value), 0);
    exact(repeated_left.rotation_coefficient(0, First), 5, 4);
    exact(repeated_left.rotation_coefficient(0, Second), 2);
    exact(repeated_right.rotation_coefficient(0, Second), 2, 1, -1);
    exact(repeated_right.translation(0, Second), 0);
}

TEST_CASE("Clocked prescribed inputs share one exact node and event side",
          "[AE-07][FD-02][FD-06][FD-07][IO4-06][N4-03][N4-06][verification][negative]") {
    auto input = a::root("aerodynamic_prescribed_transient");
    input.field("analysis")
        .set(
            "intervals",
            a::parse(
                R"([{"start":0,"end":1,"step":0.5,"restart":"initial"},{"start":1,"end":2,"step":0.5,"restart":"event"}])"));
    input.field("analysis")
        .field("wake")
        .set(
            "convection",
            a::parse(
                R"({"kind":"prescribed","velocity":{"kind":"polynomial","origin":0,"scale":1,"domain":[0,2],"coefficients":[[4,0,0],[2,0,0]]}})"));
    input.field("flow").set(
        "uniform",
        a::parse(
            R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,2],"coefficients":[[1,0,0],[3,0,0]]})"));
    auto& motion = motion_json(input);
    motion.set(
        "translation",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[1,0,0],[3,0,0]],"outside":{"kind":"reject"}})"));
    motion.set(
        "rotation_vector",
        a::parse(
            R"({"kind":"pchip","times":[0,2],"values":[[0,0,0],[0,0,1]],"outside":{"kind":"reject"}})"));
    motion.set("events", a::parse("[0,1,2]"));

    core::ResourcePlan memory{budget};
    aero::detail::UniformFlowWorkspace flow_workspace;
    p::CurveWorkspace motion_workspace;
    const c::CaseUniformFlowLimits uniform{exact_limits(), limits.curve, 100};
    const c::CaseClockedPrescribedLimits clock_limits{
        {{uniform, 8, 8}, {uniform, 2048, 50'000}},
        {{exact_limits(), 10'000, c::case_motion_frame_interval_scratch_bytes},
         limits,
         2048,
         50'000},
        budget};
    const auto prepare = [&](const a::JsonToken& json, const c::CaseClockedPrescribedLimits& caps) {
        const auto source = checked_case(json);
        return c::prepare_case_clocked_prescribed_inputs(source, attachment_field(source), stamp,
                                                         flow_workspace, motion_workspace, caps,
                                                         memory);
    };
    const std::array<std::uint32_t, 1> one_word{1};
    const std::array<std::uint32_t, 2> denominator{0, 1U << 12U};
    const p::RationalView one{1, one_word, one_word}, zero{0, {}, one_word};
    const c::CaseFlowComponentBudget component{one, {1, one_word, denominator}, zero};
    std::array<c::CaseFlowComponentBudget, 36> pose{};
    std::array<c::CaseFlowComponentBudget, 9> gradient{};
    pose.fill(component);
    gradient.fill(component);
    const std::array vector{component, component, component};
    const c::CaseKnownFlowProjection ambient{core::Precision::TwoWord, vector, vector, gradient};
    const std::optional<c::CaseFlowProjectionRequest> wake{
        std::in_place, core::Precision::TwoWord, vector, vector,
        c::CaseFlowProjectionLimits{2048, 50'000}};
    const c::CaseMotionProjection projection{core::Precision::TwoWord, pose, {}};
    const auto evaluate = [&](const c::CaseClockedPrescribedInputs& ready,
                              const c::CaseClockedPrescribedInputs& current,
                              c::CaseFlowClockQuery query,
                              const std::optional<c::CaseFlowProjectionRequest>& prescribed) {
        return c::evaluate_case_clocked_prescribed_inputs(
            ready, current, query, {zero, zero, zero}, ambient, prescribed, Second, projection,
            flow_workspace, motion_workspace, clock_limits, memory);
    };

    SECTION("the arrival and departure agree in time and retain different accelerations") {
        const auto ready = take(prepare(input, clock_limits));
        const auto current = take(prepare(input, clock_limits));
        const auto middle = take(evaluate(ready, current, {stamp, 0, 1, Automatic}, wake));
        exact(middle.ambient.time(), 1, 2);
        const auto left = take(evaluate(ready, current, {stamp, 0, 2, Left}, wake));
        const auto right = take(evaluate(ready, current, {stamp, 1, 0, Right}, wake));
        const auto verify = [](const c::CaseClockedPrescribedSample& sample, p::RequestedSide side,
                               std::uint32_t acceleration) {
            const auto& time = take(sample.motion.time()).get();
            exact(sample.ambient.time(), 1);
            exact(time.time(), 1);
            exact(time.translation(0, First), 4, 3);
            exact(time.translation(0, Second), acceleration, 3);
            exact(time.rotation_coefficient(2, First), 1, 2);
            exact(time.rotation_coefficient(2, Second), 0);
            REQUIRE(sample.prescribed);
            exact(sample.prescribed->jet().time(), 1);
            exact(sample.prescribed->jet().velocity(0), 6);
            exact(sample.prescribed->jet().time_partial(0), 2);
            const auto acceleration_word = take(sample.ambient.time_partial(0));
            const auto& dt = std::get<io::ScaledTwoWord>(acceleration_word);
            CHECK(std::ldexp(dt.high + dt.low, static_cast<int>(dt.exponent2)) == 3.0);
            const auto velocity = take(sample.ambient.velocity(0));
            const auto& enclosed = std::get<io::ScaledTwoWord>(velocity);
            CHECK(std::ldexp(enclosed.high + enclosed.low, static_cast<int>(enclosed.exponent2)) ==
                  4.0);
            CHECK(time.requested_side() == side);
            CHECK(sample.ambient.requested_side() == side);
            CHECK(sample.prescribed->requested_side() == side);
        };
        verify(left, Left, 1);
        verify(right, Right, 5);
        rejected(evaluate(ready, current, {stamp, 0, 2, Right}, wake),
                 core::FailureCode::InvalidInput);
        const p::PchipInputStamp stale{p::PchipIssuerId{71}, p::PchipInputGeneration{9}};
        rejected(evaluate(ready, current, {stale, 0, 2, Left}, wake),
                 core::FailureCode::StaleGeneration);
    }

    SECTION("the owned snapshot restores every checked input at both sides of the event") {
        const auto current = take(prepare(input, clock_limits));
        const auto restored = [&] {
            const auto archive = [&] {
                core::ResourcePlan local{budget};
                aero::detail::UniformFlowWorkspace local_fields;
                p::CurveWorkspace local_motion;
                const auto original = checked_case(input);
                const auto ready = take(c::prepare_case_clocked_prescribed_inputs(
                    original, attachment_field(original), stamp, local_fields, local_motion,
                    clock_limits, local));
                return take(c::snapshot_case_clocked_prescribed_inputs(
                    ready, local_fields, local_motion, clock_limits, local));
            }();
            const auto requalify = [&](const c::CaseClockedPrescribedInputs& live) {
                return c::restore_case_clocked_prescribed_inputs(
                    archive, live, flow_workspace, motion_workspace, clock_limits, memory);
            };
            const auto refuse_case = [&](const a::JsonToken& json, core::Operation operation) {
                const auto changed = take(prepare(json, clock_limits));
                const auto result = requalify(changed);
                if (result)
                    FAIL("Unexpected clocked restore success");
                CHECK(result.error().code == core::FailureCode::StaleGeneration);
                CHECK(result.error().operation == operation);
            };
            auto alternate = input;
            alternate.field("analysis")
                .field("intervals")
                .elements[1]
                .set("step", a::parse("0.25"));
            refuse_case(alternate, core::Operation::EvaluateAero);
            alternate = input;
            alternate.field("analysis")
                .field("intervals")
                .elements[1]
                .set("restart", a::string("step_change"));
            refuse_case(alternate, core::Operation::EvaluateAero);
            alternate = input;
            alternate.field("flow").field("uniform").set("coefficients",
                                                         a::parse("[[1,0,0],[4,0,0]]"));
            refuse_case(alternate, core::Operation::EvaluateAero);
            alternate = input;
            alternate.field("analysis")
                .field("wake")
                .field("convection")
                .field("velocity")
                .set("coefficients", a::parse("[[4,0,0],[3,0,0]]"));
            refuse_case(alternate, core::Operation::EvaluateAero);
            alternate = input;
            motion_json(alternate)
                .field("translation")
                .set("values", a::parse("[[0,0,0],[1,0,0],[4,0,0]]"));
            refuse_case(alternate, core::Operation::PrepareModel);
            alternate = input;
            alternate.field("analysis")
                .field("wake")
                .set("convection", a::parse(R"({"kind":"free"})"));
            const auto free = take(prepare(alternate, clock_limits));
            core::ResourcePlan mismatch{0};
            rejected(c::restore_case_clocked_prescribed_inputs(
                         archive, free, flow_workspace, motion_workspace, clock_limits, mismatch),
                     core::FailureCode::StaleGeneration);
            CHECK(mismatch.allocated_bytes() == 0);
            const auto free_archive = take(c::snapshot_case_clocked_prescribed_inputs(
                free, flow_workspace, motion_workspace, clock_limits, memory));
            const auto current_free = take(prepare(alternate, clock_limits));
            const auto restored_free = take(c::restore_case_clocked_prescribed_inputs(
                free_archive, current_free, flow_workspace, motion_workspace, clock_limits,
                memory));
            const auto free_sample =
                take(evaluate(restored_free, current_free, {stamp, 0, 2, Left}, {}));
            CHECK_FALSE(free_sample.prescribed);
            exact(free_sample.ambient.time(), 1);
            rejected(c::restore_case_clocked_prescribed_inputs(free_archive, current,
                                                               flow_workspace, motion_workspace,
                                                               clock_limits, mismatch),
                     core::FailureCode::StaleGeneration);
            return take(requalify(current));
        }();
        for (const auto address : {c::CaseFlowClockQuery{stamp, 0, 2, Left},
                                   c::CaseFlowClockQuery{stamp, 1, 0, Right}}) {
            const auto uninterrupted = take(evaluate(current, current, address, wake));
            const auto resumed = take(evaluate(restored, current, address, wake));
            CHECK(resumed.address.side == address.side);
            CHECK(resumed.restart == uninterrupted.restart);
            CHECK(resumed.restart == (address.side == Left ? c::CaseClockRestart::ItemInitial
                                                           : c::CaseClockRestart::ItemEvent));
            const auto step = resumed.step.view();
            exact(p::RationalView{step.sign, step.numerator, step.denominator}, 1, 2);
            const auto& motion_jet = take(resumed.motion.time()).get();
            exact(motion_jet.time(), 1);
            exact(motion_jet.translation(0, First), 4, 3);
            exact(motion_jet.translation(0, Second), address.side == Left ? 1U : 5U, 3);
            exact(motion_jet.rotation_coefficient(2, First), 1, 2);
            exact(resumed.ambient.time(), 1);
            CHECK(motion_jet.requested_side() == address.side);
            CHECK(resumed.ambient.requested_side() == address.side);
            const auto& saved_wake = resumed.prescribed;
            const auto& original_wake = uninterrupted.prescribed;
            if (!saved_wake || !original_wake)
                FAIL("The restored prescribed mode lost its velocity input");
            CHECK(saved_wake->requested_side() == address.side);
            exact(saved_wake->jet().time(), 1);
            exact(saved_wake->jet().velocity(0), 6);
            exact(saved_wake->jet().time_partial(0), 2);
            for (std::size_t k = 0; k < 3; ++k) {
                same_enclosure(resumed.ambient.velocity(k), uninterrupted.ambient.velocity(k));
                same_enclosure(resumed.ambient.time_partial(k),
                               uninterrupted.ambient.time_partial(k));
                same_enclosure(saved_wake->velocity(k), original_wake->velocity(k));
                same_enclosure(saved_wake->time_partial(k), original_wake->time_partial(k));
                for (std::size_t j = 0; j < 3; ++j)
                    same_enclosure(resumed.ambient.spatial_partial(k, j),
                                   uninterrupted.ambient.spatial_partial(k, j));
            }
            for (const auto order : {Value, First, Second}) {
                for (std::size_t i = 0; i < 3; ++i) {
                    same_enclosure(resumed.motion.centre(i, order),
                                   uninterrupted.motion.centre(i, order));
                    for (std::size_t j = 0; j < 3; ++j)
                        same_enclosure(resumed.motion.rotation(i, j, order),
                                       uninterrupted.motion.rotation(i, j, order));
                }
            }
        }
    }

    SECTION("direct evaluation refuses independently changed ready input") {
        const auto ready = take(prepare(input, clock_limits));
        input.field("analysis").field("intervals").elements[1].set("step", a::parse("0.25"));
        const auto changed_grid = take(prepare(input, clock_limits));
        rejected(evaluate(ready, changed_grid, {stamp, 0, 2, Left}, wake),
                 core::FailureCode::StaleGeneration);
        input.field("analysis").field("intervals").elements[1].set("step", a::parse("0.5"));
        motion.set("rotation_vector", a::parse(R"({"kind":"constant","value":[0,0,0]})"));
        const auto changed_motion = take(prepare(input, clock_limits));
        rejected(evaluate(ready, changed_motion, {stamp, 0, 2, Left}, wake),
                 core::FailureCode::StaleGeneration);
        input.field("analysis").field("wake").set("convection", a::parse(R"({"kind":"free"})"));
        const auto free = take(prepare(input, clock_limits));
        const auto free_sample = take(evaluate(free, free, {stamp, 0, 2, Left}, std::nullopt));
        exact(free_sample.ambient.time(), 1);
        CHECK_FALSE(free_sample.prescribed);
        rejected(evaluate(free, free, {stamp, 0, 2, Left}, wake), core::FailureCode::InvalidInput);
        rejected(evaluate(ready, free, {stamp, 0, 2, Left}, wake),
                 core::FailureCode::StaleGeneration);
    }

    SECTION("a negative periodic extension cannot cross an undeclared event") {
        motion.set(
            "translation",
            a::parse(
                R"({"kind":"pchip","times":[0,1,2,3],"values":[[0,0,0],[1,0,0],[-1,0,0],[0,0,0]],"outside":{"kind":"periodic","period":3}})"));
        motion.set("events", a::parse("[0,1,2,3]"));
        motion.set("rotation_vector", a::parse(R"({"kind":"constant","value":[0,0,0]})"));
        input.field("flow").set("uniform", a::parse(R"({"kind":"constant","value":[1,0,0]})"));
        input.field("analysis")
            .field("wake")
            .set("convection",
                 a::parse(
                     R"({"kind":"prescribed","velocity":{"kind":"constant","value":[4,0,0]}})"));
        auto& analysis = input.field("analysis");
        analysis.set("intervals",
                     a::parse(R"([{"start":-3.5,"end":-2.5,"step":0.5,"restart":"initial"}])"));
        rejected(prepare(input, clock_limits), core::FailureCode::InvalidInput);
        analysis.set(
            "intervals",
            a::parse(
                R"([{"start":-3.5,"end":-3,"step":0.5,"restart":"initial"},{"start":-3,"end":-2.5,"step":0.5,"restart":"event"}])"));
        const auto ready = take(prepare(input, clock_limits));
        const auto left = take(evaluate(ready, ready, {stamp, 0, 1, Left}, wake));
        const auto right = take(evaluate(ready, ready, {stamp, 1, 0, Right}, wake));
        exact(take(left.motion.time()).get().translation(0, Second), 4);
        exact(take(right.motion.time()).get().translation(0, Second), 4, 1, -1);
        auto capped = clock_limits;
        capped.maximum_event_work = 1;
        rejected(prepare(input, capped), core::FailureCode::ResourceExceeded);
    }
}

TEST_CASE("Prepared Motion projects a rotated Frame and PCHIP side accelerations together",
          "[AE-07][IO4-01][N4-03][N4-06][verification][negative]") {
    auto input = a::root("aerodynamic_prescribed_transient");
    input.field("units").set("time", a::string("s"));
    input.field("units").set("length", a::string("m"));
    input.field("units").set("angle", a::string("deg"));
    input.field("analysis")
        .set(
            "intervals",
            a::parse(
                R"([{"start":0,"end":1,"step":1,"restart":"initial"},{"start":1,"end":2,"step":1,"restart":"event"}])"));
    auto& motion = motion_json(input);
    motion.set("origin", a::parse(R"({"origin":[3,4,5],"rotation":[[1,0,0],[0,0,-1],[0,1,0]]})"));
    motion.set(
        "translation",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[0,1,0],[0,3,0]],"outside":{"kind":"reject"}})"));
    motion.set(
        "rotation_vector",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[0,0,90],[0,0,270]],"outside":{"kind":"reject"}})"));
    motion.set("events", a::parse("[0,1,2]"));

    core::ResourcePlan plan{budget};
    p::CurveWorkspace workspace;
    const auto begin = exact_time("0", plan);
    const auto end = exact_time("2", plan);
    const auto at = exact_time("1", plan);
    const c::CaseMotionLimits motion_limits{
        {exact_limits(), 10'000, c::case_motion_frame_interval_scratch_bytes},
        limits,
        2048,
        50'000};
    const auto prepare = [&](const a::JsonToken& json) {
        const auto si = checked_case(json);
        return take(c::prepare_case_motion(si, motion_field(si), {begin, end}, stamp, workspace,
                                           motion_limits, plan));
    };
    const auto ready = prepare(input);
    const auto current = prepare(input);
    REQUIRE(ready.events().size() == 3);
    const std::array<std::uint32_t, 1> one{1};
    const std::array<std::uint32_t, 3> denominator{0, 0, 1U << 17U};
    const p::RationalView unit{1, one, one}, absolute{1, one, denominator}, zero{0, {}, one};
    std::array<c::CaseFlowComponentBudget, 36> budgets{};
    budgets.fill({unit, absolute, zero});
    const std::array<std::uint32_t, 1> three{3}, four{4}, five{5};
    const p::RationalView x0{1, three, one}, y0{1, four, one}, z0{1, five, one};
    const auto point_budgets = std::span<const c::CaseFlowComponentBudget>{budgets}.first(9);
    const std::array points{c::CaseMotionReferencePoint{{x0, y0, z0}, point_budgets},
                            c::CaseMotionReferencePoint{{y0, y0, z0}, point_budgets}};
    const c::CaseMotionProjection projection{core::Precision::TwoWord, budgets, points};
    const auto evaluate = [&](p::RequestedSide side) {
        return c::evaluate_case_motion(ready, current, {stamp, at, Second, side}, projection,
                                       workspace, motion_limits, plan);
    };
    rejected(evaluate(Automatic), core::FailureCode::InvalidInput);
    const auto left = take(evaluate(Left));
    const auto right = take(evaluate(Right));

    using Dec = boost::multiprecision::cpp_dec_float_100;
    // Independent decimal enclosure of pi, with an error far smaller than the
    // requested projection. The Hermite accelerations below come from the
    // slopes 1/2, 4/3 and 5/2 of the respective unscaled translation ordinates.
    const Dec pi{"3.141592653589793238462643383279502884197169399375105820974944592307816406"};
    const auto expect = [&](core::Result<io::SiEnclosure> supplied, const Dec& oracle) {
        const auto enclosure = take(std::move(supplied));
        const auto* words = std::get_if<io::ScaledTwoWord>(&enclosure);
        REQUIRE(words != nullptr);
        const auto exp = static_cast<int>(words->exponent2);
        const Dec center = ldexp(Dec{words->high} + Dec{words->low}, exp);
        const Dec radius = ldexp(Dec{words->radius}, exp);
        CHECK(abs(center - oracle) <= radius + Dec{"1e-65"});
        CHECK(radius <= ldexp(Dec{1}, -81));
    };
    for (const auto* evaluated : {&left, &right}) {
        CHECK(evaluated->point_count() == 2);
        expect(evaluated->centre(0, Value), 3);
        expect(evaluated->centre(1, Value), 4);
        expect(evaluated->centre(2, Value), 6);
        expect(evaluated->centre(2, First), Dec{4} / 3);
        expect(evaluated->rotation(0, 0, Value), 0);
        expect(evaluated->rotation(1, 0, Value), 0);
        expect(evaluated->rotation(2, 0, Value), 1);
        expect(evaluated->rotation(0, 0, First), -2 * pi / 3);
        expect(evaluated->rotation(2, 0, Second), -4 * pi * pi / 9);
        expect(evaluated->point(0, 2, Value), 6);
        expect(evaluated->point(1, 0, Value), 3);
        expect(evaluated->point(1, 1, Value), 4);
        expect(evaluated->point(1, 2, Value), 7);
        expect(evaluated->point(1, 0, First), -2 * pi / 3);
        expect(evaluated->point(1, 2, First), Dec{4} / 3);
        CHECK_FALSE(evaluated->point(2, 0, Value));
    }
    expect(left.centre(2, Second), Dec{1} / 3);
    expect(right.centre(2, Second), Dec{5} / 3);
    expect(left.rotation(0, 0, Second), -pi / 6);
    expect(right.rotation(0, 0, Second), -5 * pi / 6);
    expect(left.point(1, 2, Second), Dec{1} / 3 - 4 * pi * pi / 9);
    expect(right.point(1, 2, Second), Dec{5} / 3 - 4 * pi * pi / 9);
    CHECK(take(left.time()).get().requested_side() == Left);
    CHECK(take(right.time()).get().requested_side() == Right);
    CHECK(left.interval_work() != 0);

    auto changed = input;
    motion_json(changed).field("origin").set("origin", a::parse("[3,4,6]"));
    const auto other = prepare(changed);
    rejected(c::evaluate_case_motion(ready, other, {stamp, at, Second, Left}, projection, workspace,
                                     motion_limits, plan),
             core::FailureCode::StaleGeneration);
    auto strict = budgets;
    strict[12 + 3].absolute = zero;
    rejected(c::evaluate_case_motion(ready, current, {stamp, at, Second, Left},
                                     {core::Precision::TwoWord, strict}, workspace, motion_limits,
                                     plan),
             core::FailureCode::PrecisionUnresolved);
}

TEST_CASE("Prescribed Motion restores owned exact inputs and requalifies saved coefficients",
          "[AE-07][IO4-06][FD-05][FD-06][N4-06][verification][negative]") {
    auto input = a::root("aerodynamic_prescribed_transient");
    input.field("units").set("time", a::string("s"));
    input.field("units").set("length", a::string("m"));
    input.field("units").set("angle", a::string("deg"));
    input.field("analysis")
        .set(
            "intervals",
            a::parse(
                R"([{"start":0,"end":1,"step":1,"restart":"initial"},{"start":1,"end":2,"step":1,"restart":"event"}])"));
    auto& motion = motion_json(input);
    motion.set("origin", a::parse(R"({"origin":[3,4,5],"rotation":[[1,0,0],[0,1,0],[0,0,1]]})"));
    motion.set(
        "translation",
        a::parse(
            R"({"kind":"pchip","times":[0,1,2],"values":[[0,0,0],[1,0,0],[3,0,0]],"outside":{"kind":"reject"}})"));
    motion.set(
        "rotation_vector",
        a::parse(
            R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,2],"coefficients":[[0,0,0],[0,0,90]]})"));
    motion.set("events", a::parse("[0,1,2]"));

    const c::CaseMotionLimits motion_limits{
        {exact_limits(), 10'000, c::case_motion_frame_interval_scratch_bytes},
        limits,
        2048,
        50'000};
    const auto prepare = [&](p::CurveWorkspace& workspace, core::ResourcePlan& plan) {
        const auto source = checked_case(input);
        const auto begin = exact_time("0", plan);
        const auto end = exact_time("2", plan);
        return take(c::prepare_case_motion(source, motion_field(source), {begin, end}, stamp,
                                           workspace, motion_limits, plan));
    };
    core::ResourcePlan memory{budget};
    p::CurveWorkspace workspace;
    const auto current = prepare(workspace, memory);
    const auto restored = [&] {
        const auto archive = [&] {
            core::ResourcePlan local{budget};
            p::CurveWorkspace builder;
            const auto original = prepare(builder, local);
            return take(c::snapshot_case_motion(original, builder, motion_limits, local));
        }();
        const auto saved = take(archive.view());
        const auto requalify = [&](const c::CaseMotionSnapshotView& record) {
            return c::restore_case_motion(record, current, workspace, motion_limits, memory);
        };
        const std::array<std::uint32_t, 1> one{1}, four{4}, nine{9};
        const p::RationalView changed_frame{1, four, one}, changed_coefficient{1, nine, one};
        const auto* pchip = std::get_if<p::PchipSnapshotView>(&saved.time.translation);
        const auto* formula = std::get_if<c::CaseMotionFormulaSnapshotView>(&saved.time.rotation);
        if (pchip == nullptr || formula == nullptr)
            FAIL("Expected a PCHIP translation and polynomial rotation");
        REQUIRE_FALSE(pchip->coefficients.empty());
        REQUIRE_FALSE(formula->coefficients.empty());
        REQUIRE(saved.time.events.size() == 3);
        REQUIRE(saved.time.pi_power == 1);
        exact(saved.frame.components[0], 3);
        exact(saved.time.events[1], 1);

        auto invalid = saved;
        invalid.formula_version += 1;
        rejected(requalify(invalid), core::FailureCode::CheckpointIncompatible);
        invalid = saved;
        invalid.frame.components[0] = changed_frame;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);
        const p::RationalView positive{1, one, one}, negative{-1, one, one}, zero{0, {}, one};
        const std::array rotated{zero, negative, zero, positive, zero, zero, zero, zero, positive};
        invalid = saved;
        std::ranges::copy(rotated, invalid.frame.components.begin() + 3);
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);
        invalid = saved;
        invalid.time.pi_power = 0;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);

        const auto nearby = exact_time("1.000000000000000001", memory);
        const auto nearby_words = take(io::detail::copy_rational_words(nearby, memory, {}));
        const auto near = nearby_words.view();
        const p::RationalView distinct_event{near.sign, near.numerator, near.denominator};
        auto events = std::array{saved.time.events[0], distinct_event, saved.time.events[2]};
        invalid = saved;
        invalid.time.events = events;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);
        invalid = saved;
        invalid.time.end = distinct_event;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);

        std::vector<p::RationalView> raw{pchip->input.samples.values.begin(),
                                         pchip->input.samples.values.end()};
        REQUIRE_FALSE(raw.empty());
        raw[0] = changed_coefficient;
        auto altered = *pchip;
        altered.input.samples.values = raw;
        invalid = saved;
        invalid.time.translation = altered;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);
        altered = *pchip;
        altered.input.boundary.policy = p::HoldBoundary{};
        invalid.time.translation = altered;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);

        std::vector<p::RationalView> generated{pchip->coefficients.begin(),
                                               pchip->coefficients.end()};
        generated[0] = changed_coefficient;
        altered = *pchip;
        altered.coefficients = generated;
        invalid.time.translation = altered;
        rejected(requalify(invalid), core::FailureCode::InvalidInput);
        invalid = saved;
        std::vector<p::RationalView> polynomial{formula->coefficients.begin(),
                                                formula->coefficients.end()};
        polynomial[0] = changed_coefficient;
        auto altered_formula = *formula;
        altered_formula.coefficients = polynomial;
        invalid.time.rotation = altered_formula;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);
        altered_formula = *formula;
        if (const auto& domain = altered_formula.domain) {
            auto changed_domain = *domain;
            changed_domain[3] = changed_frame;
            altered_formula.domain = changed_domain;
        } else {
            FAIL("The saved rotation lost its polynomial domain");
        }
        invalid = saved;
        invalid.time.rotation = altered_formula;
        rejected(requalify(invalid), core::FailureCode::StaleGeneration);

        auto zero_work = motion_limits;
        zero_work.maximum_interval_work = 0;
        core::ResourcePlan capped{budget};
        rejected(c::restore_case_motion(saved, current, workspace, zero_work, capped),
                 core::FailureCode::ResourceExceeded);
        CHECK(capped.allocated_bytes() == 0);

        core::ResourcePlan refused{0};
        rejected(c::restore_case_motion(saved, current, workspace, motion_limits, refused),
                 core::FailureCode::ResourceExceeded);
        CHECK(refused.allocated_bytes() == 0);
        // Kept until all independent corruptions have been refused; these calls may not consume it.
        return take(requalify(saved));
    }();

    const auto at = exact_time("1", memory);
    const std::array<std::uint32_t, 1> one{1};
    const std::array<std::uint32_t, 3> denominator{0, 0, 1U << 17U};
    const p::RationalView unit{1, one, one}, absolute{1, one, denominator}, zero{0, {}, one};
    std::array<c::CaseFlowComponentBudget, 36> budgets{};
    budgets.fill({unit, absolute, zero});
    const c::CaseMotionProjection projection{core::Precision::TwoWord, budgets};
    for (const auto side : {Left, Right}) {
        const auto evaluate = [&](const c::CaseMotionReady& prepared) {
            return take(c::evaluate_case_motion(prepared, current, {stamp, at, Second, side},
                                                projection, workspace, motion_limits, memory));
        };
        const auto uninterrupted = evaluate(current);
        const auto resumed = evaluate(restored);
        const auto& time = take(resumed.time()).get();
        exact(time.time(), 1);
        exact(time.translation(0, First), 4, 3);
        exact(time.translation(0, Second), side == Left ? 1U : 5U, 3);
        exact(time.rotation_coefficient(2, First), 1, 2);
        CHECK(time.rotation_pi_power() == 1);
        CHECK(time.requested_side() == side);
        for (const auto order : {Value, First, Second}) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                same_enclosure(resumed.centre(axis, order), uninterrupted.centre(axis, order));
                for (std::size_t col = 0; col < 3; ++col)
                    same_enclosure(resumed.rotation(axis, col, order),
                                   uninterrupted.rotation(axis, col, order));
            }
        }
    }
}

TEST_CASE("The checked case clock includes knots of every declared Motion consumer",
          "[AE-07][FD-02][FD-06][FD-07][N4-06][verification][negative]") {
    auto base = a::root("coupled_transient");
    const auto fixed_motion = a::parse(
        R"({"origin":{"origin":[0,0,0],"rotation":[[1,0,0],[0,1,0],[0,0,1]]},"translation":{"kind":"constant","value":[0,0,0]},"rotation_vector":{"kind":"constant","value":[0,0,0]},"events":[]})");
    auto& baseline_surfaces = base.field("model").field("aero_surfaces").elements;
    auto selected = baseline_surfaces.front();
    selected.set("id", a::string("2"));
    auto attachment = a::parse(R"({"kind":"prescribed_motion"})");
    attachment.set("motion", fixed_motion);
    selected.set("attachment", attachment);
    baseline_surfaces.insert(baseline_surfaces.begin(), selected);
    const auto periodic = a::parse(
        R"({"kind":"pchip","times":[0,1,2,3],"values":[[0,0,0],[1,0,0],[-1,0,0],[0,0,0]],"outside":{"kind":"periodic","period":3}})");
    auto recurring = fixed_motion;
    recurring.set("translation", periodic);
    // 1.25 is an additional absolute event, not a periodic knot.
    recurring.set("events", a::parse("[0,1,1.25,2,3]"));
    const auto split = a::parse(
        R"([{"start":-3.5,"end":-3,"step":0.5,"restart":"initial"},{"start":-3,"end":-2.5,"step":0.5,"restart":"event"}])");
    const auto unsplit = a::parse(R"([{"start":-3.5,"end":-2.5,"step":0.5,"restart":"initial"}])");

    const c::CaseUniformFlowLimits uniform{exact_limits(), limits.curve, 100};
    const c::CaseClockedPrescribedLimits clock_limits{
        {{uniform, 8, 8}, {uniform, 2048, 50'000}},
        {{exact_limits(), 10'000, c::case_motion_frame_interval_scratch_bytes},
         limits,
         2048,
         50'000},
        budget};
    core::ResourcePlan memory{budget};
    aero::detail::UniformFlowWorkspace fields;
    p::CurveWorkspace motion_workspace;
    const auto prepare = [&](const a::JsonToken& json) {
        const auto source = checked_case(json);
        return c::prepare_case_clocked_prescribed_inputs(source, attachment_field(source), stamp,
                                                         fields, motion_workspace, clock_limits,
                                                         memory);
    };
    const std::array<std::uint32_t, 1> one_word{1};
    const std::array<std::uint32_t, 2> denominator{0, 1U << 12U};
    const p::RationalView one{1, one_word, one_word}, zero{0, {}, one_word};
    const c::CaseFlowComponentBudget component{one, {1, one_word, denominator}, zero};
    std::array<c::CaseFlowComponentBudget, 36> pose{};
    std::array<c::CaseFlowComponentBudget, 9> gradient{};
    pose.fill(component);
    gradient.fill(component);
    const std::array vector{component, component, component};
    const c::CaseKnownFlowProjection ambient{core::Precision::TwoWord, vector, vector, gradient};
    const c::CaseMotionProjection projection{core::Precision::TwoWord, pose};
    const auto evaluate = [&](const c::CaseClockedPrescribedInputs& ready,
                              const c::CaseClockedPrescribedInputs& current,
                              c::CaseFlowClockQuery node) {
        return c::evaluate_case_clocked_prescribed_inputs(ready, current, node, {zero, zero, zero},
                                                          ambient, {}, Second, projection, fields,
                                                          motion_workspace, clock_limits, memory);
    };

    constexpr std::array<std::string_view, 5> consumers{"another aero attachment", "point support",
                                                        "trace support", "elastic reference",
                                                        "gravity"};
    for (const auto consumer : consumers) {
        CAPTURE(consumer);
        auto input = base;
        auto& model = input.field("model");
        auto& analysis = input.field("analysis");
        auto& supports = model.field("supports").elements;
        if (consumer == "another aero attachment") {
            auto other = selected;
            other.set("id", a::string("3"));
            other.field("attachment").set("motion", recurring);
            model.field("aero_surfaces").elements.push_back(std::move(other));
        } else if (consumer == "point support") {
            auto support = a::parse(
                R"({"kind":"constraint","id":"2","target":{"target":{"kind":"shell","id":"1"},"coordinates":[0,0,0],"side":"interior"},"components":[0]})");
            support.set("motion", recurring);
            supports.push_back(std::move(support));
        } else if (consumer == "trace support") {
            auto with_rotation = fixed_motion;
            with_rotation.set("rotation_vector", periodic);
            with_rotation.set("events", recurring.field("events"));
            supports.front().set("motion", std::move(with_rotation));
        } else if (consumer == "elastic reference") {
            auto elastic = a::parse(
                R"({"kind":"elastic","id":"2","target":{"target":{"kind":"shell","id":"1"},"coordinates":[0,0,0],"side":"interior"},"stiffness":[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]],"damping":[[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]],"frame":{"origin":[0,0,0],"rotation":[[1,0,0],[0,1,0],[0,0,1]]}})");
            elastic.set("reference", recurring);
            supports.push_back(std::move(elastic));
        } else {
            auto load = a::parse(R"({"kind":"gravity","id":"1"})");
            load.set("acceleration", periodic);
            analysis.field("loads").elements.push_back(std::move(load));
        }

        analysis.set("intervals", unsplit);
        rejected(prepare(input), core::FailureCode::InvalidInput);
        analysis.set("intervals", split);
        const auto ready = take(prepare(input));
        const auto current = take(prepare(input));
        const auto archive = take(c::snapshot_case_clocked_prescribed_inputs(
            ready, fields, motion_workspace, clock_limits, memory));
        const auto restored = take(c::restore_case_clocked_prescribed_inputs(
            archive, current, fields, motion_workspace, clock_limits, memory));
        for (const auto address : {c::CaseFlowClockQuery{stamp, 0, 1, Left},
                                   c::CaseFlowClockQuery{stamp, 1, 0, Right}}) {
            const auto sample = take(evaluate(restored, current, address));
            const auto& jet = take(sample.motion.time()).get();
            exact(sample.ambient.time(), 3, 1, -1);
            exact(jet.time(), 3, 1, -1);
            CHECK(jet.requested_side() == address.side);
            CHECK(sample.ambient.requested_side() == address.side);
            exact(jet.translation(0, Second), 0);
        }
        rejected(evaluate(restored, current, {stamp, 0, 1, Right}),
                 core::FailureCode::InvalidInput);

        if (consumer == "trace support") {
            auto different = input;
            different.field("model")
                .field("supports")
                .elements.front()
                .field("motion")
                .field("rotation_vector")
                .set("values", a::parse("[[0,0,0],[2,0,0],[-2,0,0],[0,0,0]]"));
            const auto altered = take(prepare(different));
            rejected(c::restore_case_clocked_prescribed_inputs(
                         archive, altered, fields, motion_workspace, clock_limits, memory),
                     core::FailureCode::StaleGeneration);

            supports.front().field("motion").set("events", a::parse("[-3.25,0,1,1.25,2,3]"));
            analysis.set(
                "intervals",
                a::parse(
                    R"([{"start":-3.5,"end":-3.25,"step":0.25,"restart":"initial"},{"start":-3.25,"end":-3,"step":0.25,"restart":"event"},{"start":-3,"end":-2.5,"step":0.5,"restart":"event"}])"));
            const auto with_absolute = take(prepare(input));
            for (const auto node : {c::CaseFlowClockQuery{stamp, 0, 1, Left},
                                    c::CaseFlowClockQuery{stamp, 1, 0, Right}}) {
                const auto actual = take(evaluate(with_absolute, with_absolute, node));
                exact(actual.ambient.time(), 13, 4, -1);
                CHECK(take(actual.motion.time()).get().requested_side() == node.side);
            }
        }

        // Here -2 is a periodic knot at the first endpoint. The extra absolute
        // event 1.25 would be inside the interval if incorrectly repeated as -1.75.
        analysis.set("intervals",
                     a::parse(R"([{"start":-2,"end":-1.5,"step":0.5,"restart":"initial"}])"));
        REQUIRE(prepare(input));
    }
}
