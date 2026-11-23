#include "aero/pchip_curve.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <fstream>
#include <future>
#include <map>
#include <vector>

namespace {
namespace p = uvlas::aero::detail::pchip;
namespace core = uvlas::core;
namespace linear = uvlas::numerics::linear;
using Rat = boost::multiprecision::cpp_rational;
using Big = boost::multiprecision::cpp_int;

template <class T>
T take(core::Result<T> result) {
    if (!result)
        FAIL("Unexpected failure " << static_cast<int>(result.error().code));
    return std::move(*result);
}

template <class T>
void failed(const core::Result<T>& result, core::FailureCode code) {
    if (result)
        FAIL("Unexpected success");
    CHECK(result.error().code == code);
    CHECK(result.error().operation == core::Operation::EvaluateAero);
}

struct WordValue {
    std::vector<std::uint32_t> numerator, denominator;
    int sign = 0;

    explicit WordValue(const Rat& value) {
        Big n = boost::multiprecision::numerator(value);
        sign = n == 0 ? 0 : n < 0 ? -1 : 1;
        if (n < 0)
            n = -n;
        if (n != 0)
            boost::multiprecision::export_bits(n, std::back_inserter(numerator), 32, false);
        const auto d = boost::multiprecision::denominator(value);
        boost::multiprecision::export_bits(d, std::back_inserter(denominator), 32, false);
    }

    [[nodiscard]] p::RationalView view() const noexcept {
        return {sign, numerator, denominator};
    }
};

Rat decode(p::RationalView value) {
    Big n, d;
    for (std::size_t i = value.numerator.size(); i != 0; --i)
        n = (n << 32) + value.numerator[i - 1];
    for (std::size_t i = value.denominator.size(); i != 0; --i)
        d = (d << 32) + value.denominator[i - 1];
    REQUIRE(d > 0);
    REQUIRE(value.denominator.back() != 0);
    if (n == 0) {
        CHECK(value.numerator.empty());
        CHECK(value.sign == 0);
        CHECK(d == 1);
    } else {
        CHECK(value.numerator.back() != 0);
        CHECK((value.sign == -1 || value.sign == 1));
    }
    const Rat exact = Rat{value.sign < 0 ? -n : n} / d;
    CHECK(boost::multiprecision::denominator(exact) == d);
    return exact;
}

struct Input {
    using Sequence = std::vector<p::RationalView>;
    std::vector<WordValue> words;
    Sequence times, values, events;
    std::vector<uvlas::numerics::SiDimension> dimensions;

    Input(const std::vector<Rat>& t, const std::vector<Rat>& y, std::size_t components,
          const std::vector<Rat>& e = {})
        : dimensions(components) {
        words.reserve(t.size() + y.size() + e.size());
        const auto load = [&](const auto& values_to_load, auto& views) {
            for (const auto& value : values_to_load)
                views.push_back(words.emplace_back(value).view());
        };
        load(t, times);
        load(y, values);
        load(e, events);
        for (auto& dimension : dimensions) {
            dimension.powers[0] = 1;
            dimension.powers[2] = -1;
        }
    }

    [[nodiscard]] p::PreparationRequest request() const noexcept {
        return {{times, values, dimensions.size()},
                dimensions,
                {times.front(), times.back(), p::RejectBoundary{}},
                events,
                {p::PchipIssuerId{4}, p::PchipInputGeneration{6}}};
    }
};

p::CurveLimits limits() {
    return {{4096, 1000000, UINT64_MAX, 100000}, linear::EstimateAdmission{4096}};
}
p::SharedPreparationInput host_input(const p::PreparationRequest& input) {
    const auto cap = limits();
    core::ResourcePlan memory{UINT64_MAX};
    p::PreparationWorkspace preparation;
    p::SharedPreparationWorkspace sharing;
    auto owned = take(p::make_preparation_input(input, preparation, cap.arithmetic, memory));
    return take(p::share_preparation_input(std::move(owned), sharing,
                                           {cap.metadata, cap.arithmetic.steps}, memory));
}

p::PchipCurveReady prepare(const p::PreparationRequest& input) {
    p::CurveWorkspace workspace;
    auto cap = limits();
    const auto plan = take(p::forecast_pchip_prepare(input, cap));
    cap.arithmetic.steps = plan.reserved_steps;
    cap.metadata = linear::EstimateAdmission{plan.metadata_estimate_bytes};
    core::ResourcePlan memory{plan.payload_bytes};
    auto curve = take(p::prepare_pchip_curve(input, workspace, cap, memory));
    CHECK(memory.allocated_bytes() <= plan.payload_bytes);
    CHECK(curve.used_steps() <= plan.reserved_steps);
    return curve;
}

p::ExactPchipJet evaluate(const p::PchipCurveReady& curve, const p::SharedPreparationInput& input,
                          const Rat& time, p::JetOrder order,
                          p::RequestedSide side = p::RequestedSide::Automatic) {
    const WordValue query{time};
    const auto request = p::PchipQuery{query.view(), order, side};
    p::CurveWorkspace workspace;
    auto cap = limits();
    const auto plan = take(p::forecast_pchip_evaluation(curve, input, request, cap));
    cap.arithmetic.steps = plan.reserved_steps;
    cap.metadata = linear::EstimateAdmission{plan.metadata_estimate_bytes};
    core::ResourcePlan memory{plan.payload_bytes};
    auto jet = take(p::evaluate_pchip_curve(curve, input, request, workspace, cap, memory));
    CHECK(memory.allocated_bytes() <= plan.payload_bytes);
    CHECK(jet.used_steps() <= plan.reserved_steps);
    CHECK(take(jet.order()) == order);
    CHECK(take(jet.location()).requested == side);
    CHECK(decode(take(jet.original_time())) == time);
    return jet;
}

struct ReferenceCurve {
    p::PchipCurveReady curve;
    p::SharedPreparationInput input;
};
} // namespace

TEST_CASE("Qualified PCHIP curves reproduce exact independent coefficient and Hermite fixtures",
          "[AE-07][FD-02][FD-06][N4-02][N4-03][N4-06][verification]") {
    std::ifstream coefficient_stream{UVLAS_PCHIP_COEFFICIENT_REFERENCE};
    std::ifstream jet_stream{UVLAS_PCHIP_INTERVAL_REFERENCE};
    REQUIRE(coefficient_stream.good());
    REQUIRE(jet_stream.good());
    const auto coefficient_reference = nlohmann::json::parse(coefficient_stream);
    const auto jet_reference = nlohmann::json::parse(jet_stream);
    std::map<std::string, ReferenceCurve> sources;
    for (const auto& record : coefficient_reference.at("cases")) {
        CAPTURE(record.at("name"));
        std::vector<Rat> times, values;
        for (const auto& v : record.at("times"))
            times.emplace_back(v.get<std::string>());
        for (const auto& v : record.at("values"))
            values.emplace_back(v.get<std::string>());
        const auto components = record.at("components").get<std::uint64_t>();
        const Input input{times, values, components};
        const auto request = input.request();
        auto curve = prepare(request);
        const auto& coefficients = record.at("coefficients");
        for (std::size_t index = 0; index < coefficients.size(); ++index)
            for (std::size_t power = 0; power < 4; ++power)
                CHECK(decode(
                          take(curve.coefficient(index / components, index % components, power))) ==
                      Rat{coefficients.at(index).at(power).get<std::string>()});
        sources.try_emplace(record.at("name").get<std::string>(),
                            ReferenceCurve{std::move(curve), host_input(request)});
    }
    const std::array<const char*, 3> names{"value", "dt", "dtt"};
    for (const auto& record : jet_reference.at("cases")) {
        CAPTURE(record.at("source_name"), record.at("interval"), record.at("u"));
        const auto& prepared = sources.at(record.at("source_name").get<std::string>());
        const auto time = Rat{record.at("time").get<std::string>()};
        const auto side = record.at("u") == "1" ? p::RequestedSide::Left : p::RequestedSide::Right;
        for (const auto order : {p::JetOrder::Value, p::JetOrder::First, p::JetOrder::Second}) {
            const auto jet = evaluate(prepared.curve, prepared.input, time, order, side);
            CHECK(take(jet.location()).interval == record.at("interval").get<std::uint64_t>());
            CHECK(decode(take(jet.normalized_time())) == time);
            for (std::uint64_t component = 0; component < jet.components(); ++component) {
                for (std::uint64_t d = 0; d < names.size(); ++d) {
                    const auto value = jet.derivative(component, static_cast<p::JetOrder>(d));
                    if (d <= static_cast<std::uint64_t>(order))
                        CHECK(decode(take(value)) == Rat{record.at("values")
                                                             .at(component)
                                                             .at(names[static_cast<std::size_t>(d)])
                                                             .get<std::string>()});
                    else
                        failed(value, core::FailureCode::InvalidInput);
                }
            }
        }
    }
}

TEST_CASE("Qualified PCHIP policies select the physical side and enforce C1 on preparation",
          "[AE-07][N4-03][N4-06][N4-08][property][negative]") {
    const Input identity{{0, 1, 3}, {0, 1, 3}, 1};
    const auto reject = identity.request();
    const auto reject_curve = prepare(reject);
    const auto reject_input = host_input(reject);
    const WordValue before{-1}, after{4}, first{0}, end{3}, negative{Rat{-1} / 2};
    p::CurveWorkspace workspace;
    core::ResourcePlan memory{UINT64_MAX};
    failed(p::evaluate_pchip_curve(reject_curve, reject_input, {before.view(), p::JetOrder::First},
                                   workspace, limits(), memory),
           core::FailureCode::InvalidInput);
    const auto linear_join = evaluate(reject_curve, reject_input, 1, p::JetOrder::Second);
    CHECK(decode(take(linear_join.derivative(0, p::JetOrder::Second))) == 0);
    const Input corner{{0, 1, 3}, {0, 1, 5}, 1};
    const auto corner_request = corner.request();
    const auto corner_curve = prepare(corner_request);
    const auto corner_input = host_input(corner_request);
    const WordValue join{1};
    failed(p::evaluate_pchip_curve(corner_curve, corner_input, {join.view(), p::JetOrder::Second},
                                   workspace, limits(), memory),
           core::FailureCode::InvalidInput);
    const auto corner_left =
        evaluate(corner_curve, corner_input, 1, p::JetOrder::Second, p::RequestedSide::Left);
    const auto corner_right =
        evaluate(corner_curve, corner_input, 1, p::JetOrder::Second, p::RequestedSide::Right);
    CHECK(decode(take(corner_left.derivative(0, p::JetOrder::First))) ==
          decode(take(corner_right.derivative(0, p::JetOrder::First))));
    CHECK(decode(take(corner_left.derivative(0, p::JetOrder::Second))) !=
          decode(take(corner_right.derivative(0, p::JetOrder::Second))));

    const WordValue half{Rat{1} / 2}, two_half{Rat{5} / 2};
    auto hold = identity.request();
    hold.boundary.policy = p::HoldBoundary{};
    failed(p::prepare_pchip_curve(hold, workspace, limits(), memory),
           core::FailureCode::InvalidInput);
    hold.boundary.begin = half.view();
    hold.boundary.end = two_half.view();
    const auto hold_curve = prepare(hold);
    const auto hold_input = host_input(hold);
    const auto outside = evaluate(hold_curve, hold_input, 4, p::JetOrder::Second);
    CHECK_FALSE(take(outside.location()).interval);
    CHECK(decode(take(outside.normalized_time())) == 4);
    CHECK(decode(take(outside.derivative(0, p::JetOrder::Value))) == 3);
    CHECK(decode(take(outside.derivative(0, p::JetOrder::First))) == 0);
    failed(p::evaluate_pchip_curve(hold_curve, hold_input, {first.view(), p::JetOrder::First},
                                   workspace, limits(), memory),
           core::FailureCode::InvalidInput);
    const auto held =
        evaluate(hold_curve, hold_input, 0, p::JetOrder::Second, p::RequestedSide::Left);
    const auto inside =
        evaluate(hold_curve, hold_input, 0, p::JetOrder::Second, p::RequestedSide::Right);
    CHECK(decode(take(held.derivative(0, p::JetOrder::First))) == 0);
    CHECK(decode(take(inside.derivative(0, p::JetOrder::First))) == 1);

    const Input flat_start{{0, 1, 2}, {0, 1, 5}, 1};
    auto hold_flat = flat_start.request();
    hold_flat.boundary.policy = p::HoldBoundary{};
    hold_flat.boundary.end = join.view();
    const auto hold_flat_curve = prepare(hold_flat);
    const auto hold_flat_input = host_input(hold_flat);
    const auto admitted_first = evaluate(hold_flat_curve, hold_flat_input, 0, p::JetOrder::First);
    CHECK(decode(take(admitted_first.derivative(0, p::JetOrder::First))) == 0);
    failed(p::evaluate_pchip_curve(hold_flat_curve, hold_flat_input,
                                   {first.view(), p::JetOrder::Second}, workspace, limits(),
                                   memory),
           core::FailureCode::InvalidInput);

    const Input cyclic{{0, 1, 2, 4, 5}, {0, 0, 1, 0, 0}, 1};
    const WordValue cycle{5};
    auto periodic = cyclic.request();
    periodic.boundary.policy = p::PeriodicBoundary{cycle.view()};
    const auto periodic_curve = prepare(periodic);
    const auto periodic_input = host_input(periodic);
    const auto wrapped = evaluate(periodic_curve, periodic_input, Rat{-1} / 2, p::JetOrder::Second);
    CHECK(decode(take(wrapped.normalized_time())) == Rat{9} / 2);
    CHECK(decode(take(wrapped.derivative(0, p::JetOrder::Value))) == 0);
    const auto seam = evaluate(periodic_curve, periodic_input, -5, p::JetOrder::Second);
    CHECK(decode(take(seam.normalized_time())) == 0);
    CHECK(take(seam.location()).interval == 0);
    const auto seam_left =
        evaluate(periodic_curve, periodic_input, -5, p::JetOrder::Second, p::RequestedSide::Left);
    CHECK(take(seam_left.location()).interval == 3);

    const Input asymmetric_seam{{0, 1, 2, 3, 5}, {0, 1, 0, Rat{-4} / 3, 0}, 1};
    auto asymmetric = asymmetric_seam.request();
    asymmetric.boundary.policy = p::PeriodicBoundary{cycle.view()};
    const auto asymmetric_curve = prepare(asymmetric);
    const auto asymmetric_input = host_input(asymmetric);
    const auto asymmetric_left = evaluate(asymmetric_curve, asymmetric_input, -5,
                                          p::JetOrder::Second, p::RequestedSide::Left);
    const auto asymmetric_right = evaluate(asymmetric_curve, asymmetric_input, -5,
                                           p::JetOrder::Second, p::RequestedSide::Right);
    CHECK(decode(take(asymmetric_left.derivative(0, p::JetOrder::First))) == 2);
    CHECK(decode(take(asymmetric_right.derivative(0, p::JetOrder::First))) == 2);
    CHECK(decode(take(asymmetric_left.derivative(0, p::JetOrder::Second))) == 2);
    CHECK(decode(take(asymmetric_right.derivative(0, p::JetOrder::Second))) == -2);
    const WordValue negative_seam{-5};
    failed(p::evaluate_pchip_curve(asymmetric_curve, asymmetric_input,
                                   {negative_seam.view(), p::JetOrder::Second}, workspace, limits(),
                                   memory),
           core::FailureCode::InvalidInput);
    const Input slope_mismatch{{0, 1, 2, 3, 5}, {0, 1, 0, -1, 0}, 1};
    auto bad_slope = slope_mismatch.request();
    bad_slope.boundary.policy = p::PeriodicBoundary{cycle.view()};
    failed(p::prepare_pchip_curve(bad_slope, workspace, limits(), memory),
           core::FailureCode::InvalidInput);
    auto bad_period = periodic;
    bad_period.boundary.policy = p::PeriodicBoundary{after.view()};
    failed(p::prepare_pchip_curve(bad_period, workspace, limits(), memory),
           core::FailureCode::InvalidInput);
    auto bad_seam = identity.request();
    bad_seam.boundary.policy = p::PeriodicBoundary{end.view()};
    failed(p::prepare_pchip_curve(bad_seam, workspace, limits(), memory),
           core::FailureCode::InvalidInput);
    failed(p::evaluate_pchip_curve(periodic_curve, periodic_input,
                                   {negative.view(), static_cast<p::JetOrder>(3)}, workspace,
                                   limits(), memory),
           core::FailureCode::InvalidInput);
}

TEST_CASE("PCHIP snapshot retains original descriptor and restores only requalified coefficients",
          "[AE-07][FD-02][FD-06][N4-06][N4-08][N4-09][verification][negative]") {
    const Input definition{{-1, 0, 2}, {-2, 5, 8, 0, 4, 9, 4, 2, 11}, 3, {-3, 0, 4}};
    const auto current = host_input(definition.request());
    auto snapshot = [&] {
        const auto source = prepare(definition.request());
        p::CurveWorkspace workspace;
        const auto saving = take(p::forecast_pchip_snapshot(source, limits()));
        core::ResourcePlan memory{saving.payload_bytes};
        return take(p::snapshot_pchip_curve(source, workspace, limits(), memory));
    }();
    const auto saved = take(snapshot.view());
    REQUIRE(saved.input.events.size() == 3);
    CHECK(decode(saved.input.events.front()) == -3);
    CHECK(decode(saved.input.events.back()) == 4);
    CHECK(saved.input.stamp == definition.request().stamp);
    const auto planned = take(p::forecast_pchip_restore(saved, current, limits()));
    p::CurveWorkspace workspace;
    core::ResourcePlan memory{planned.payload_bytes};
    auto cap = limits();
    cap.arithmetic.steps = planned.reserved_steps;
    cap.metadata = linear::EstimateAdmission{planned.metadata_estimate_bytes};
    auto restored = take(p::restore_pchip_curve(saved, current, workspace, cap, memory));
    CHECK(restored.used_steps() <= planned.reserved_steps);
    CHECK(memory.allocated_bytes() <= planned.payload_bytes);
    const auto expected = host_input(definition.request());
    // The caller uses its own original owner, and subsequently an independently identical owner.
    const auto same = evaluate(restored, current, 1, p::JetOrder::Second);
    const auto different = evaluate(restored, expected, 1, p::JetOrder::Second);
    for (std::uint64_t component = 0; component < 3; ++component)
        for (const auto d : {p::JetOrder::Value, p::JetOrder::First, p::JetOrder::Second})
            CHECK(decode(take(same.derivative(component, d))) ==
                  decode(take(different.derivative(component, d))));
    CHECK(decode(take(different.derivative(0, p::JetOrder::Value))) == 2);
    CHECK(decode(take(different.derivative(1, p::JetOrder::First))) == -1);
    CHECK(decode(take(different.derivative(2, p::JetOrder::Second))) == 0);

    const WordValue query{1};
    auto older = definition.request();
    older.stamp.generation = p::PchipInputGeneration{5};
    const auto outdated = host_input(older);
    core::ResourcePlan fail_memory{UINT64_MAX};
    failed(p::restore_pchip_curve(saved, outdated, workspace, limits(), fail_memory),
           core::FailureCode::StaleGeneration);
    failed(p::evaluate_pchip_curve(restored, outdated, {query.view()}, workspace, limits(),
                                   fail_memory),
           core::FailureCode::StaleGeneration);
    CHECK(fail_memory.allocated_bytes() == 0);
    auto other_issuer = definition.request();
    other_issuer.stamp.issuer = p::PchipIssuerId{5};
    const auto issued_elsewhere = host_input(other_issuer);
    failed(p::restore_pchip_curve(saved, issued_elsewhere, workspace, limits(), fail_memory),
           core::FailureCode::StaleGeneration);

    auto altered = definition.request();
    const WordValue small{Rat{1} / (Big{1} << 200)};
    auto values = definition.values;
    values[3] = small.view();
    altered.samples.values = values;
    const auto content_changed = host_input(altered);
    failed(p::restore_pchip_curve(saved, content_changed, workspace, limits(), fail_memory),
           core::FailureCode::StaleGeneration);
    failed(p::evaluate_pchip_curve(restored, content_changed, {query.view()}, workspace, limits(),
                                   fail_memory),
           core::FailureCode::StaleGeneration);
    altered = definition.request();
    altered.events = altered.events.subspan(1);
    failed(p::restore_pchip_curve(saved, host_input(altered), workspace, limits(), fail_memory),
           core::FailureCode::StaleGeneration);

    auto corrupt = saved;
    ++corrupt.formula_version;
    failed(p::restore_pchip_curve(corrupt, current, workspace, limits(), fail_memory),
           core::FailureCode::InvalidInput);
    corrupt = saved;
    corrupt.input.samples.values = corrupt.input.samples.values.subspan(1);
    failed(p::restore_pchip_curve(corrupt, current, workspace, limits(), fail_memory),
           core::FailureCode::InvalidInput);
    corrupt = saved;
    corrupt.coefficients = corrupt.coefficients.subspan(1);
    failed(p::restore_pchip_curve(corrupt, current, workspace, limits(), fail_memory),
           core::FailureCode::InvalidInput);
    CHECK(fail_memory.allocated_bytes() == 0);

    const WordValue wrong{17};
    std::vector<p::RationalView> coefficients(saved.coefficients.begin(), saved.coefficients.end());
    coefficients[0] = wrong.view(); // A canonical number in the correct shape still cannot restore.
    corrupt.coefficients = coefficients;
    failed(p::restore_pchip_curve(corrupt, current, workspace, limits(), fail_memory),
           core::FailureCode::InvalidInput);
    const std::array<std::uint32_t, 1> numerator{4}, denominator{2};
    coefficients[0] = {-1, numerator, denominator}; // Numerically original, noncanonical.
    failed(p::restore_pchip_curve(corrupt, current, workspace, limits(), fail_memory),
           core::FailureCode::InvalidInput);
    CHECK(decode(take(restored.coefficient(0, 0, 0))) == -2);

    auto jet = [&] {
        auto moved_curve = std::move(restored);
        // NOLINTNEXTLINE(bugprone-use-after-move): inactivity is the tested contract.
        CHECK_FALSE(restored.active());
        failed(restored.input(), core::FailureCode::InvalidInput);
        return evaluate(moved_curve, expected, 1, p::JetOrder::Second);
    }();
    auto moved = std::move(jet);
    // NOLINTNEXTLINE(bugprone-use-after-move): inactivity is the tested contract.
    CHECK_FALSE(jet.active());
    failed(jet.original_time(), core::FailureCode::InvalidInput);
    CHECK(take(moved.stamp()) == saved.input.stamp);
    CHECK(take(moved.formula_version()) == saved.formula_version);
    CHECK(decode(take(moved.derivative(2, p::JetOrder::Value))) == 10);

    const auto move_saved = std::move(snapshot);
    // NOLINTNEXTLINE(bugprone-use-after-move): inactivity is the tested contract.
    failed(snapshot.view(), core::FailureCode::InvalidInput);
    CHECK(take(move_saved.view()).coefficients.size() == 24);
}

TEST_CASE(
    "Exact PCHIP results and work agree with independent workspaces on one two and four threads",
    "[AE-07][N4-06][N4-09][property]") {
    const Input fixture{{0, 1, 3}, {0, 1, 3}, 1};
    const auto input = host_input(fixture.request());
    const auto curve = prepare(fixture.request());
    const WordValue time{1};
    const p::PchipQuery query{time.view(), p::JetOrder::Second};
    const auto run = [&]() {
        const auto plan = p::forecast_pchip_evaluation(curve, input, query, limits());
        if (!plan)
            return core::Result<p::ExactPchipJet>{std::unexpected(plan.error())};
        core::ResourcePlan memory{plan->payload_bytes};
        p::CurveWorkspace workspace;
        return p::evaluate_pchip_curve(curve, input, query, workspace, limits(), memory);
    };
    const auto baseline = take(run());
    for (const std::size_t peers : {std::size_t{1}, std::size_t{2}, std::size_t{4}}) {
        std::vector<std::future<core::Result<p::ExactPchipJet>>> results;
        results.reserve(peers);
        for (std::size_t i = 0; i < peers; ++i)
            results.push_back(std::async(std::launch::async, run));
        for (auto& each : results) {
            const auto jet = take(each.get());
            CHECK(jet.used_steps() == baseline.used_steps());
            CHECK(take(jet.location()).selected == take(baseline.location()).selected);
            for (const auto derivative :
                 {p::JetOrder::Value, p::JetOrder::First, p::JetOrder::Second})
                CHECK(decode(take(jet.derivative(0, derivative))) ==
                      decode(take(baseline.derivative(0, derivative))));
        }
    }
}
