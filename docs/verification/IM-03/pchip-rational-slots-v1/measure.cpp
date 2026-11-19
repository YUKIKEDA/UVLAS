#include "pchip_rational_slots.hpp"

#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <utility>
#include <uvlas/core/instrumentation.hpp>
namespace p = uvlas::aero::detail::pchip;
namespace c = uvlas::core;
template <class T>
T take(c::Result<T> value) {
    if (!value)
        std::terminate();
    return std::move(*value);
}
p::RationalPlan prepare() {
    p::Counter counter{{64, 256, UINT64_MAX, 16}};
    auto a = take(counter.input({3, 2}));
    auto b = take(counter.input({3, 4}));
    auto sum = take(counter.binary(p::RationalOp::Add, a, b));
    take(counter.floor(sum));
    return take(counter.finish());
}
int main() {
    constexpr std::uint64_t repeats = 10000;
    std::uint64_t sink = 0;
    auto before = take(c::sample_process_memory(c::Operation::EvaluateAero));
    const auto start = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < repeats; ++i)
        sink += prepare().reserved_steps();
    const auto counted = std::chrono::steady_clock::now();
    const auto plan = prepare();
    for (std::uint64_t i = 0; i < repeats; ++i) {
        c::ResourcePlan memory{plan.payload_bytes()};
        auto workspace = take(p::Workspace::create(plan, memory));
        sink += memory.allocated_bytes();
    }
    const auto allocated = std::chrono::steady_clock::now();
    c::ResourcePlan memory{plan.payload_bytes()};
    auto workspace = take(p::Workspace::create(plan, memory));
    auto held = take(c::sample_process_memory(c::Operation::EvaluateAero));
    std::uint64_t actual = 0;
    const std::array<std::uint32_t, 1> n{7}, d{3}, bn{5}, bd{11};
    const auto evaluation_start = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < repeats; ++i) {
        auto session = take(workspace.begin());
        auto a = take(session.input({-1, n, d}));
        auto b = take(session.input({1, bn, bd}));
        auto sum = take(session.binary(p::RationalOp::Add, a, b));
        auto result = take(session.floor(sum));
        auto view = take(session.view(result));
        if (view.sign != -1 || view.numerator.size() != 1 || view.numerator[0] != 2 ||
            view.denominator[0] != 1)
            return 2;
        actual = take(session.usage()).steps;
        sink += actual;
    }
    const auto finish = std::chrono::steady_clock::now();
    auto after = take(c::sample_process_memory(c::Operation::EvaluateAero));
    const auto ns = [](auto a, auto b) {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
    };
    std::cout << "{\n  \"repeats\":" << repeats << ",\n  \"count_total_ns\":" << ns(start, counted)
              << ",\n  \"allocate_destroy_total_ns\":" << ns(counted, allocated)
              << ",\n  \"evaluate_total_ns\":" << ns(evaluation_start, finish)
              << ",\n  \"bits\":" << plan.integer_bits() << ",\n  \"words\":" << plan.words()
              << ",\n  \"slots\":" << plan.slots()
              << ",\n  \"payload_bytes\":" << plan.payload_bytes()
              << ",\n  \"reserved_steps\":" << plan.reserved_steps()
              << ",\n  \"actual_steps\":" << actual
              << ",\n  \"resident_before\":" << before.resident_bytes
              << ",\n  \"resident_held\":" << held.resident_bytes
              << ",\n  \"resident_after\":" << after.resident_bytes
              << ",\n  \"lifetime_peak_after\":" << after.peak_resident_bytes
              << ",\n  \"sink\":" << sink << "\n}\n";
}
