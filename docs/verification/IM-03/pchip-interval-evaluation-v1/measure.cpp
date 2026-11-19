#include "pchip_interval_jet.hpp"

#include <array>
#include <chrono>
#include <iostream>
#include <uvlas/core/build_info.hpp>
#include <uvlas/core/instrumentation.hpp>
namespace p = uvlas::aero::detail::pchip;
namespace c = uvlas::core;
template <class T>
T take(c::Result<T> value) {
    if (!value)
        std::terminate();
    return std::move(*value);
}
int main() {
    constexpr std::uint64_t repeats = 1000;
    const std::array<std::uint32_t, 1> one{1}, two{2}, four{4};
    const std::array<p::RationalView, 3> times{{{0, {}, one}, {1, one, one}, {1, two, one}}},
        values{{{0, {}, one}, {1, one, one}, {1, four, one}}};
    const p::RationalLimits limits{64, 100000, UINT64_MAX, 10000};
    const p::SampleGrid grid{times, values, 1};
    const auto source_plan = take(p::forecast_coefficients(grid, limits));
    c::ResourcePlan source_memory{source_plan.payload_bytes()};
    p::CoefficientWorkspace source_workspace;
    const auto source =
        take(p::prepare_coefficients(grid, source_workspace, limits, source_memory));
    const p::IntervalRequest query{0, {1, one, two}, p::JetOrder::Second};
    const auto plan = take(p::forecast_interval(source, query, limits));
    p::IntervalWorkspace workspace;
    c::ResourcePlan held_memory{plan.payload_bytes()};
    const auto held = take(p::evaluate_interval(source, query, workspace, limits, held_memory));
    const auto before = take(c::sample_process_memory(c::Operation::EvaluateAero));
    std::uint64_t sink = 0;
    const auto begin = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < repeats; ++i)
        sink += take(p::forecast_interval(source, query, limits)).reserved_steps();
    const auto counted = std::chrono::steady_clock::now();
    const auto after_count = take(c::sample_process_memory(c::Operation::EvaluateAero));
    const auto start = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < repeats; ++i) {
        c::ResourcePlan memory{plan.payload_bytes()};
        const auto result = take(p::evaluate_interval(source, query, workspace, limits, memory));
        const std::array<std::uint32_t, 3> numerator{5, 9, 3}, denominator{16, 8, 2};
        for (std::uint64_t d = 0; d < 3; ++d) {
            const auto value = take(result.derivative(0, static_cast<p::JetOrder>(d)));
            const auto j = static_cast<std::size_t>(d);
            if (value.sign != 1 || value.numerator.size() != 1 || value.denominator.size() != 1 ||
                value.numerator[0] != numerator[j] || value.denominator[0] != denominator[j])
                return 2;
        }
        sink += result.used_steps();
    }
    const auto end = std::chrono::steady_clock::now();
    const auto after = take(c::sample_process_memory(c::Operation::EvaluateAero));
    const auto ns = [](auto a, auto b) {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
    };
    std::cout << "{\n\"repeats\":" << repeats << ",\n\"count_total_ns\":" << ns(begin, counted)
              << ",\n\"evaluate_destroy_total_ns\":" << ns(start, end)
              << ",\n\"payload_bytes\":" << plan.payload_bytes()
              << ",\n\"reserved_steps\":" << plan.reserved_steps()
              << ",\n\"actual_steps\":" << held.used_steps()
              << ",\n\"reserved_slots\":" << plan.arithmetic().slots()
              << ",\n\"actual_slots\":" << held.used_slots()
              << ",\n\"bits\":" << plan.arithmetic().integer_bits()
              << ",\n\"words\":" << plan.arithmetic().words()
              << ",\n\"resident_before\":" << before.resident_bytes
              << ",\n\"resident_after_count\":" << after_count.resident_bytes
              << ",\n\"resident_after\":" << after.resident_bytes
              << ",\n\"lifetime_peak_after\":" << after.peak_resident_bytes
              << ",\n\"sink\":" << sink << ",\n\"build_info\":" << uvlas::build_info_json()
              << "\n}\n";
}
