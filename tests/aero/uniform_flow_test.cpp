#include "aero/uniform_flow.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <limits>

namespace {
namespace a = uvlas::aero::detail;
namespace p = a::pchip;
namespace c = uvlas::core;
namespace backend = uvlas::numerics::linear;

constexpr auto max = std::numeric_limits<std::uint64_t>::max();
const p::CurveLimits cap{{96, 100000, max, 10000}, backend::EstimateAdmission{max}};
const p::PchipInputStamp stamp{p::PchipIssuerId{7}, p::PchipInputGeneration{13}};

class Fraction {
public:
    Fraction(std::int32_t num, std::uint32_t den = 1)
        : magnitude_{static_cast<std::uint32_t>(num < 0 ? -static_cast<std::int64_t>(num) : num)},
          denominator_{den},
          sign_((num > 0) - (num < 0)) {}

    [[nodiscard]] p::RationalView view() const& {
        return {sign_, sign_ == 0 ? std::span<const std::uint32_t>{} : magnitude_, denominator_};
    }
    [[nodiscard]] p::RationalView view() const&& = delete;

private:
    std::array<std::uint32_t, 1> magnitude_, denominator_;
    int sign_;
};

template <std::size_t N>
std::array<p::RationalView, N> views(const std::array<Fraction, N>& source) {
    std::array<p::RationalView, N> result{};
    std::ranges::transform(source, result.begin(),
                           [](const Fraction& value) { return value.view(); });
    return result;
}

template <class T>
T take(c::Result<T> result) {
    if (!result)
        FAIL("Unexpected uniform-flow failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}

template <class T>
void reject(const c::Result<T>& result, c::FailureCode expected) {
    if (result)
        FAIL("Unexpected uniform-flow success");
    CHECK(result.error().code == expected);
    CHECK(result.error().operation == c::Operation::EvaluateAero);
}

void exact(p::RationalView actual, const Fraction& literal) {
    const auto expected = literal.view();
    CHECK(actual.sign == expected.sign);
    CHECK(std::ranges::equal(actual.numerator, expected.numerator));
    CHECK(std::ranges::equal(actual.denominator, expected.denominator));
}

void check_jet(const a::ExactUniformFlowJet& jet, const Fraction& time,
               const std::array<Fraction, 3>& value, const std::array<Fraction, 3>& derivative) {
    CHECK(take(jet.stamp()) == stamp);
    exact(take(jet.time()), time);
    for (std::size_t i = 0; i < 3; ++i) {
        CAPTURE(i);
        exact(take(jet.velocity(i)), value[i]);
        exact(take(jet.time_partial(i)), derivative[i]);
        for (std::size_t j = 0; j < 3; ++j)
            exact(take(jet.spatial_partial(i, j)), Fraction{0});
    }
    reject(jet.spatial_partial(0, 3), c::FailureCode::InvalidInput);
}

a::UniformFlowReady prepare(const a::UniformFlowInput& input, a::UniformFlowWorkspace& workspace) {
    const auto forecast = take(a::forecast_uniform_flow_prepare(input, cap));
    auto limits = cap;
    limits.arithmetic.steps = forecast.reserved_steps;
    limits.metadata = backend::EstimateAdmission{forecast.metadata_estimate_bytes};
    c::ResourcePlan exact_memory{forecast.payload_bytes};
    auto result = take(a::prepare_uniform_flow(input, workspace, limits, exact_memory));
    CHECK(result.used_steps() <= forecast.reserved_steps);
    CHECK(exact_memory.allocated_bytes() == forecast.payload_bytes);
    REQUIRE(forecast.payload_bytes > 0);
    c::ResourcePlan too_small{forecast.payload_bytes - 1};
    reject(a::prepare_uniform_flow(input, workspace, limits, too_small),
           c::FailureCode::ResourceExceeded);
    CHECK(too_small.allocated_bytes() == 0);
    --limits.arithmetic.steps;
    c::ResourcePlan refused{forecast.payload_bytes};
    reject(a::prepare_uniform_flow(input, workspace, limits, refused),
           c::FailureCode::ResourceExceeded);
    CHECK(refused.allocated_bytes() == 0);
    return result;
}

p::SharedPreparationInput current_pchip(const a::UniformFlowInput& input) {
    const auto* series = std::get_if<a::PchipUniformFlow>(&input.series);
    if (!series)
        FAIL("Expected PCHIP current input");
    constexpr uvlas::numerics::SiDimension velocity{{1, 0, -1, 0, 0, 0, 0}};
    constexpr std::array dimensions{velocity, velocity, velocity};
    p::PreparationWorkspace workspace;
    p::SharedPreparationWorkspace shared;
    c::ResourcePlan memory{max};
    auto original = take(p::make_preparation_input(
        {series->samples, dimensions, series->boundary, input.events, input.stamp}, workspace,
        cap.arithmetic, memory));
    return take(p::share_preparation_input(std::move(original), shared,
                                           {cap.metadata, cap.arithmetic.steps}, memory));
}

a::ExactUniformFlowJet evaluate(const a::UniformFlowReady& field, const Fraction& time,
                                a::UniformFlowWorkspace& workspace,
                                p::RequestedSide side = p::RequestedSide::Automatic,
                                const p::SharedPreparationInput* expected = nullptr) {
    a::UniformFlowQuery request{stamp, time.view(), side};
    if (expected)
        request.expected_pchip.emplace(*expected);
    const auto forecast = take(a::forecast_uniform_flow_evaluation(field, request, cap));
    auto limits = cap;
    limits.arithmetic.steps = forecast.reserved_steps;
    limits.metadata = backend::EstimateAdmission{forecast.metadata_estimate_bytes};
    c::ResourcePlan exact_memory{forecast.payload_bytes};
    auto result = take(a::evaluate_uniform_flow(field, request, workspace, limits, exact_memory));
    CHECK(result.used_steps() <= forecast.reserved_steps);
    CHECK(exact_memory.allocated_bytes() == forecast.payload_bytes);
    REQUIRE(forecast.payload_bytes > 0);
    c::ResourcePlan too_small{forecast.payload_bytes - 1};
    reject(a::evaluate_uniform_flow(field, request, workspace, limits, too_small),
           c::FailureCode::ResourceExceeded);
    CHECK(too_small.allocated_bytes() == 0);
    --limits.arithmetic.steps;
    c::ResourcePlan refused{forecast.payload_bytes};
    reject(a::evaluate_uniform_flow(field, request, workspace, limits, refused),
           c::FailureCode::ResourceExceeded);
    CHECK(refused.allocated_bytes() == 0);
    return result;
}

a::UniformFlowSnapshot snapshot(const a::UniformFlowReady& field, a::UniformFlowWorkspace& work) {
    const auto plan = take(a::forecast_uniform_flow_snapshot(field, cap));
    c::ResourcePlan memory{plan.payload_bytes};
    auto result = take(a::snapshot_uniform_flow(field, work, cap, memory));
    CHECK(memory.allocated_bytes() == plan.payload_bytes);
    return result;
}

a::UniformFlowReady restore(const a::UniformFlowSnapshotView& saved,
                            const a::UniformFlowReady& current,
                            a::UniformFlowWorkspace& workspace) {
    const auto forecast = take(a::forecast_uniform_flow_restore(saved, current, cap));
    auto limits = cap;
    limits.arithmetic.steps = forecast.reserved_steps;
    limits.metadata = backend::EstimateAdmission{forecast.metadata_estimate_bytes};
    c::ResourcePlan memory{forecast.payload_bytes};
    auto result = take(a::restore_uniform_flow(saved, current, workspace, limits, memory));
    CHECK(result.used_steps() <= forecast.reserved_steps);
    CHECK(memory.allocated_bytes() == forecast.payload_bytes);
    --limits.arithmetic.steps;
    c::ResourcePlan refused{forecast.payload_bytes};
    reject(a::restore_uniform_flow(saved, current, workspace, limits, refused),
           c::FailureCode::ResourceExceeded);
    CHECK(refused.allocated_bytes() == 0);
    return result;
}
} // namespace

TEST_CASE("Exact uniform flow has spatially constant velocity and survives both source lifetimes",
          "[AE-07][FD-06][N4-02][N4-03][N4-08][N4-09][property]") {
    a::UniformFlowWorkspace workspace;
    const std::array velocity{Fraction{9, 7}, Fraction{-2, 3}, Fraction{0}};
    const std::array zeros{Fraction{0}, Fraction{0}, Fraction{0}};
    const Fraction time{-11, 5};
    const auto saving = [&] {
        auto source = velocity;
        const a::UniformFlowInput request{stamp, {}, a::ConstantUniformFlow{views(source)}};
        const auto field = prepare(request, workspace);
        source[0] = Fraction{3};
        return snapshot(field, workspace);
    }();
    const auto actual = take(saving.view());
    const auto result = [&] {
        const a::UniformFlowInput request{stamp, {}, a::ConstantUniformFlow{views(velocity)}};
        const auto independently_prepared = prepare(request, workspace);
        const auto restored = restore(actual, independently_prepared, workspace);
        return evaluate(restored, time, workspace);
    }();
    check_jet(result, time, velocity, zeros);
}

TEST_CASE("Uniform polynomial gives the exact SI chain rule on its closed domain",
          "[AE-07][FD-06][N4-02][N4-03][N4-08][N4-09][property][negative]") {
    const Fraction origin{1, 3}, scale{2}, begin{-1}, end{3}, time{4, 3}, outside{13, 4};
    const std::array coefficients{Fraction{1},  Fraction{-2}, Fraction{0},  Fraction{2},
                                  Fraction{0},  Fraction{3},  Fraction{-1}, Fraction{4},
                                  Fraction{-2}, Fraction{3},  Fraction{0},  Fraction{1}};
    const auto records = views(coefficients);
    const std::array events{Fraction{-1, 2}, Fraction{3, 2}};
    const auto schedule = views(events);
    const a::UniformFlowInput request{
        stamp, schedule,
        a::PolynomialUniformFlow{origin.view(), scale.view(), begin.view(), end.view(), records}};

    a::UniformFlowWorkspace work;
    const auto field = prepare(request, work);
    const auto current = prepare(request, work);
    const auto original = take(field.input());
    REQUIRE(original.events.size() == 2);
    exact(original.events[0], events[0]);
    exact(original.events[1], events[1]);

    // u=1/2. Expanding the cubic and differentiating each term gives these fractions.
    const std::array value{Fraction{17, 8}, Fraction{-1}, Fraction{9, 8}};
    const std::array derivative{Fraction{13, 8}, Fraction{2}, Fraction{7, 8}};
    const auto jet = evaluate(field, time, work);
    check_jet(jet, time, value, derivative);
    const auto saved = snapshot(field, work);
    const auto saved_view = take(saved.view());
    const auto restored = restore(saved_view, current, work);
    const auto replay = evaluate(restored, time, work);
    check_jet(replay, time, value, derivative);
    CHECK(evaluate(field, begin, work).active());
    CHECK(evaluate(field, end, work, p::RequestedSide::Left).active());

    c::ResourcePlan memory{max};
    reject(a::evaluate_uniform_flow(field, {stamp, outside.view()}, work, cap, memory),
           c::FailureCode::InvalidInput);
    const p::PchipInputStamp old{stamp.issuer, p::PchipInputGeneration{12}};
    reject(a::evaluate_uniform_flow(field, {old, time.view()}, work, cap, memory),
           c::FailureCode::StaleGeneration);
    const auto* source = std::get_if<a::UniformFormulaSnapshotView>(&saved_view);
    if (!source)
        FAIL("Expected owned polynomial replay");
    auto changed = *source;
    changed.input.events = {};
    const auto before = memory.allocated_bytes();
    reject(a::restore_uniform_flow(changed, current, work, cap, memory),
           c::FailureCode::StaleGeneration);
    CHECK(memory.allocated_bytes() == before);

    check_jet(jet, time, value, derivative);
}

TEST_CASE("Uniform PCHIP replays the actual qualified coefficients and physical time derivative",
          "[AE-07][FD-06][N4-02][N4-03][N4-08][N4-09][property][negative]") {
    const std::array times{Fraction{0}, Fraction{1}, Fraction{2}};
    const std::array samples{Fraction{0},  Fraction{-2}, Fraction{1},  Fraction{1}, Fraction{-2},
                             Fraction{-1}, Fraction{4},  Fraction{-2}, Fraction{-7}};
    const auto x = views(times);
    const auto y = views(samples);
    const a::UniformFlowInput request{
        stamp, {}, a::PchipUniformFlow{{x, y, 3}, {x.front(), x.back(), p::RejectBoundary{}}}};
    const Fraction time{1, 2};
    a::UniformFlowWorkspace work;
    const auto field = prepare(request, work);
    const auto current = prepare(request, work);
    const auto owned_save = snapshot(field, work);
    const auto saved_view = take(owned_save.view());
    const auto restored = restore(saved_view, current, work);
    const auto expected = current_pchip(request);
    // On [0,1] the first component is 3t²/2-t³/2. The third is 1-2 times the first.
    const std::array value{Fraction{5, 16}, Fraction{-2}, Fraction{3, 8}};
    const std::array derivative{Fraction{9, 8}, Fraction{0}, Fraction{-9, 4}};
    const auto jet = evaluate(restored, time, work, p::RequestedSide::Automatic, &expected);
    check_jet(jet, time, value, derivative);

    for (auto side : {p::RequestedSide::Left, p::RequestedSide::Right}) {
        const auto at_knot = evaluate(field, times[1], work, side, &expected);
        check_jet(at_knot, times[1], {Fraction{1}, Fraction{-2}, Fraction{-1}},
                  {Fraction{3, 2}, Fraction{0}, Fraction{-3}});
    }

    const auto* saved_pchip = std::get_if<p::PchipSnapshotView>(&saved_view);
    if (!saved_pchip)
        FAIL("Expected actual PCHIP coefficient snapshot");
    REQUIRE(saved_pchip->coefficients.size() == 24);
    std::array<p::RationalView, 24> edited{};
    std::ranges::copy(saved_pchip->coefficients, edited.begin());
    edited.front() = times[1].view();
    auto forged = *saved_pchip;
    forged.coefficients = edited;
    c::ResourcePlan memory{max};
    reject(a::restore_uniform_flow(forged, current, work, cap, memory),
           c::FailureCode::InvalidInput);
    check_jet(jet, time, value, derivative);
}

TEST_CASE("Uniform PCHIP honours declared endpoint hold and periodic normalization",
          "[AE-07][FD-06][N4-02][N4-03][property][negative]") {
    const std::array times{Fraction{0}, Fraction{1}, Fraction{2}, Fraction{3}, Fraction{4}};
    const std::array samples{
        Fraction{0},  Fraction{1},  Fraction{-1}, Fraction{0},  Fraction{1},
        Fraction{-1}, Fraction{2},  Fraction{1},  Fraction{-3}, Fraction{0},
        Fraction{1},  Fraction{-1}, Fraction{0},  Fraction{1},  Fraction{-1},
    };
    const auto x = views(times);
    const auto y = views(samples);
    const Fraction outside{-2};
    const std::array zeros{Fraction{0}, Fraction{0}, Fraction{0}};
    a::UniformFlowWorkspace work;
    for (const auto& [policy, value] :
         std::array<std::pair<decltype(p::BoundaryRequest::policy), std::array<Fraction, 3>>, 2>{{
             {p::HoldBoundary{}, {Fraction{0}, Fraction{1}, Fraction{-1}}},
             {p::PeriodicBoundary{x.back()}, {Fraction{2}, Fraction{1}, Fraction{-3}}},
         }}) {
        const a::UniformFlowInput request{
            stamp, {}, a::PchipUniformFlow{{x, y, 3}, {x.front(), x.back(), policy}}};
        const auto field = prepare(request, work);
        const auto expected = current_pchip(request);
        const auto jet = evaluate(field, outside, work, p::RequestedSide::Automatic, &expected);
        check_jet(jet, outside, value, zeros);
    }
}

TEST_CASE("Uniform PCHIP forecast and evaluation check the independent current source",
          "[AE-07][N4-02][N4-08][N4-09][negative]") {
    const std::array times{Fraction{0}, Fraction{1}, Fraction{2}};
    const std::array samples{Fraction{0},  Fraction{-2}, Fraction{1},  Fraction{1}, Fraction{-2},
                             Fraction{-1}, Fraction{4},  Fraction{-2}, Fraction{-7}};
    const auto x = views(times);
    const auto y = views(samples);
    const a::PchipUniformFlow source{{x, y, 3}, {x.front(), x.back(), p::RejectBoundary{}}};
    const a::UniformFlowInput input{stamp, {}, source};
    const Fraction time{1, 2};
    a::UniformFlowWorkspace work;
    const auto field = prepare(input, work);
    const auto expected = current_pchip(input);
    const auto good = evaluate(field, time, work, p::RequestedSide::Automatic, &expected);
    const std::array value{Fraction{5, 16}, Fraction{-2}, Fraction{3, 8}};
    const std::array derivative{Fraction{9, 8}, Fraction{0}, Fraction{-9, 4}};
    check_jet(good, time, value, derivative);

    auto edited_values = y;
    edited_values[4] = times[1].view();
    auto edited = source;
    edited.samples.values = edited_values;
    const auto changed_samples = current_pchip({stamp, {}, edited});
    const std::array events{time.view()};
    const auto changed_events = current_pchip({stamp, events, source});
    const auto check_refusal = [&](a::UniformFlowQuery query, c::FailureCode code) {
        c::ResourcePlan memory{max};
        reject(a::forecast_uniform_flow_evaluation(field, query, cap), code);
        reject(a::evaluate_uniform_flow(field, query, work, cap, memory), code);
        CHECK(memory.allocated_bytes() == 0);
    };
    check_refusal({stamp, time.view()}, c::FailureCode::InvalidInput);
    for (const auto& current : {std::cref(changed_samples), std::cref(changed_events)})
        check_refusal({stamp, time.view(), p::RequestedSide::Automatic, current},
                      c::FailureCode::StaleGeneration);
    check_jet(good, time, value, derivative);
}

TEST_CASE("Uniform flow rejects malformed input and cancelled admission",
          "[AE-07][FD-06][N4-02][N4-08][N4-09][negative]") {
    const std::array coefficient{Fraction{0}, Fraction{1}, Fraction{2}};
    const auto row = views(coefficient);
    const Fraction zero{0}, upper{2}, one{1}, noncanonical{2, 4};
    const auto base =
        a::PolynomialUniformFlow{zero.view(), one.view(), zero.view(), upper.view(), row};
    a::UniformFlowWorkspace work;
    const auto check_invalid = [&](const a::UniformFlowInput& request) {
        c::ResourcePlan memory{max};
        reject(a::prepare_uniform_flow(request, work, cap, memory), c::FailureCode::InvalidInput);
    };
    auto bad = base;
    bad.scale = zero.view();
    check_invalid({stamp, {}, bad});
    bad = base;
    bad.begin = upper.view();
    check_invalid({stamp, {}, bad});
    const std::array repeated{one.view(), one.view()};
    check_invalid({stamp, repeated, base});
    auto nonreduced = row;
    nonreduced[0] = noncanonical.view();
    check_invalid({stamp, {}, a::ConstantUniformFlow{nonreduced}});

    const a::UniformFlowInput good{stamp, {}, base};
    const auto field = prepare(good, work);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 declares request_stop non-const.
    std::stop_source cancel;
    cancel.request_stop();
    c::ResourcePlan refused{max};
    reject(a::prepare_uniform_flow(good, work, cap, refused, cancel.get_token()),
           c::FailureCode::Cancelled);
    reject(a::evaluate_uniform_flow(field, {stamp, one.view()}, work, cap, refused,
                                    cancel.get_token()),
           c::FailureCode::Cancelled);
    CHECK(refused.allocated_bytes() == 0);
}
