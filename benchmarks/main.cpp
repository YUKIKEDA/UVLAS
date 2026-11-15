#include <benchmark/benchmark.h>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <thread>
#include <uvlas/core/build_info.hpp>

namespace {
[[nodiscard]] std::string environment_variable(const char* name) {
#if defined(_MSC_VER)
    std::size_t required = 0;
    if (getenv_s(&required, nullptr, 0, name) != 0 || required == 0) {
        return {};
    }
    std::string value(required, '\0');
    if (getenv_s(&required, value.data(), value.size(), name) != 0) {
        return {};
    }
    value.resize(required - 1);
    return value;
#else
    const char* value = std::getenv(name);
    return value == nullptr ? std::string{} : std::string{value};
#endif
}

int run(int argc, char* argv[]) {
    const std::string machine = environment_variable("UVLAS_BENCH_MACHINE");
    const std::string cpu = environment_variable("UVLAS_BENCH_CPU");
    const std::string memory_text = environment_variable("UVLAS_BENCH_MEMORY_BYTES");
    if (machine.empty() || cpu.empty() || memory_text.empty()) {
        std::cerr << "Set UVLAS_BENCH_MACHINE, UVLAS_BENCH_CPU and UVLAS_BENCH_MEMORY_BYTES.\n";
        return 2;
    }
    const std::string_view memory{memory_text};
    std::uint64_t memory_bytes = 0;
    const auto parsed = std::from_chars(memory.data(), memory.data() + memory.size(), memory_bytes);
    if (parsed.ec != std::errc{} || parsed.ptr != memory.data() + memory.size() ||
        memory_bytes == 0) {
        std::cerr << "Physical memory must be a positive decimal byte count.\n";
        return 2;
    }
    const auto logical_processors = std::thread::hardware_concurrency();
    if (logical_processors == 0) {
        std::cerr << "Logical processor count could not be determined.\n";
        return 2;
    }
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 2;
    }
    const auto info = uvlas::build_info();
    benchmark::AddCustomContext("uvlas.machine", machine);
    benchmark::AddCustomContext("uvlas.cpu", cpu);
    benchmark::AddCustomContext("uvlas.physical_memory_bytes", std::to_string(memory_bytes));
    benchmark::AddCustomContext("uvlas.logical_processors", std::to_string(logical_processors));
    benchmark::AddCustomContext("uvlas.version", std::string{info.version});
    benchmark::AddCustomContext("uvlas.commit", std::string{info.commit});
    benchmark::AddCustomContext("uvlas.dirty", info.uncommitted_changes ? "true" : "false");
    benchmark::AddCustomContext("uvlas.compiler", std::string{info.compiler});
    benchmark::AddCustomContext("uvlas.platform", std::string{info.platform});
    benchmark::AddCustomContext("uvlas.platform_key", std::string{info.platform_key});
    benchmark::AddCustomContext("uvlas.floating_point", std::string{info.floating_point});
    benchmark::AddCustomContext("uvlas.configuration", std::string{info.configuration});
    benchmark::AddCustomContext("uvlas.instruction_set", std::string{info.instruction_set});
    benchmark::AddCustomContext("uvlas.threads", "1");
    benchmark::AddCustomContext(
        "uvlas.settings",
        R"({"checked_extents":[5000,64,3],"dense_orders":[16,64],"dense_rhs_columns":1,"dense_update_rank":0,"dense_refinements_max":5,"dense_owning_output_included":true})");
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
} // namespace

int main(int argc, char* argv[]) {
    try {
        return run(argc, argv);
    } catch (...) {
        return 1;
    }
}
