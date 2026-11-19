#include "pchip_coefficients.hpp"

#include <array>
#include <chrono>
#include <cstdlib>
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
    const p::SampleGrid grid{times, values, 1};
    const p::RationalLimits limits{64, 100000, UINT64_MAX, 10000};
    const auto plan = take(p::forecast_coefficients(grid, limits));
    p::CoefficientWorkspace workspace;
    c::ResourcePlan held_memory{plan.payload_bytes()};
    const auto held = take(p::prepare_coefficients(grid, workspace, limits, held_memory));
    const auto before = take(c::sample_process_memory(c::Operation::EvaluateAero));
    std::uint64_t sink = 0;
    const auto start = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < repeats; ++i)
        sink += take(p::forecast_coefficients(grid, limits)).reserved_steps();
    const auto counted = std::chrono::steady_clock::now();
    const auto after_count = take(c::sample_process_memory(c::Operation::EvaluateAero));
    const auto prepare_start = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < repeats; ++i) {
        c::ResourcePlan memory{plan.payload_bytes()};
        const auto block = take(p::prepare_coefficients(grid, workspace, limits, memory));
        const std::array<int, 8> numerators{0, 0, 3, -1, 1, 3, 2, -1},
            denominators{1, 1, 2, 2, 1, 2, 1, 2};
        for (std::uint64_t interval = 0; interval < 2; ++interval)
            for (std::uint64_t power = 0; power < 4; ++power) {
                const auto view = take(block.coefficient(interval, 0, power));
                const auto index = static_cast<std::size_t>(4 * interval + power);
                const auto expected = numerators[index];
                const auto sign = expected < 0 ? -1 : expected > 0 ? 1 : 0;
                if (view.sign != sign || view.denominator.size() != 1 ||
                    view.denominator[0] != static_cast<std::uint32_t>(denominators[index]))
                    return 2;
                if (expected == 0) {
                    if (!view.numerator.empty())
                        return 3;
                } else if (view.numerator.size() != 1 ||
                           view.numerator[0] !=
                               static_cast<std::uint32_t>(expected < 0 ? -expected : expected))
                    return 4;
            }
        sink += block.used_steps();
    }
    const auto finished = std::chrono::steady_clock::now();
    const auto after = take(c::sample_process_memory(c::Operation::EvaluateAero));
    const auto ns = [](auto a, auto b) {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
    };
    std::cout << "{\n\"repeats\":" << repeats << ",\n\"count_total_ns\":" << ns(start, counted)
              << ",\n\"prepare_destroy_total_ns\":" << ns(prepare_start, finished)
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
