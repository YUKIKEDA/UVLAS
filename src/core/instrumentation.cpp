#include <nlohmann/json.hpp>

#include <limits>
#include <uvlas/core/build_info.hpp>
#include <uvlas/core/instrumentation.hpp>

#if defined(_WIN32)
#include <windows.h>

// PSAPI declarations require the Windows types above.
#include <psapi.h>
#elif defined(__linux__)
#include <cstdio>
#include <sys/resource.h>
#include <unistd.h>
#endif

namespace uvlas::core {
namespace {
constexpr std::array<std::string_view, static_cast<std::size_t>(MeasurementStage::Count)>
    stage_names{"preparation",     "structure", "aerodynamics",   "linear_solve",
                "nonlinear_solve", "time_step", "flutter_search", "checkpoint"};
}

Result<void> StageRecorder::record_sample(MeasurementStage stage, std::uint64_t elapsed) noexcept {
    if (!enabled_) {
        return {};
    }
    if (failure_) {
        return std::unexpected(*failure_);
    }
    const auto index = static_cast<std::size_t>(stage);
    if (index >= timings_.size()) {
        failure_ = Failure{FailureCode::InvalidInput, operation_};
        return std::unexpected(*failure_);
    }
    auto& value = timings_[index];
    if (value.calls == std::numeric_limits<std::uint64_t>::max() ||
        elapsed > std::numeric_limits<std::uint64_t>::max() - value.elapsed_nanoseconds) {
        failure_ = Failure{FailureCode::ResourceExceeded, operation_};
        return std::unexpected(*failure_);
    }
    ++value.calls;
    value.elapsed_nanoseconds += elapsed;
    return {};
}

Result<std::array<StageTiming, static_cast<std::size_t>(MeasurementStage::Count)>>
StageRecorder::snapshot() const noexcept {
    if (failure_) {
        return std::unexpected(*failure_);
    }
    return timings_;
}

StageScope::StageScope(StageRecorder& recorder, MeasurementStage stage) noexcept
    : recorder_(recorder), stage_(stage) {
    if (recorder_.enabled()) {
        start_ = std::chrono::steady_clock::now();
    }
}

StageScope::~StageScope() noexcept {
    if (!start_) {
        return;
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                             std::chrono::steady_clock::now() - *start_)
                             .count();
    // A negative elapsed value indicates an invalid duration; route it through the
    // recorder's retained failure rather than manufacturing a successful measurement.
    const auto recorded =
        elapsed < 0 ? recorder_.record_sample(MeasurementStage::Count, 0)
                    : recorder_.record_sample(stage_, static_cast<std::uint64_t>(elapsed));
    if (!recorded) {
        return;
    }
}

Result<ProcessMemory> sample_process_memory(Operation operation) noexcept {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS counters{};
    counters.cb = sizeof(counters);
    if (!GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))) {
        return std::unexpected(Failure{FailureCode::Unsupported, operation});
    }
    return ProcessMemory{static_cast<std::uint64_t>(counters.WorkingSetSize),
                         static_cast<std::uint64_t>(counters.PeakWorkingSetSize)};
#elif defined(__linux__)
    unsigned long long total_pages = 0;
    unsigned long long resident_pages = 0;
    std::FILE* file = std::fopen("/proc/self/statm", "r");
    if (file == nullptr) {
        return std::unexpected(Failure{FailureCode::Unsupported, operation});
    }
    const int fields = std::fscanf(file, "%llu %llu", &total_pages, &resident_pages);
    const int closed = std::fclose(file);
    const long page_size = sysconf(_SC_PAGESIZE);
    rusage usage{};
    if (fields != 2 || closed != 0 || page_size <= 0 || getrusage(RUSAGE_SELF, &usage) != 0 ||
        usage.ru_maxrss < 0) {
        return std::unexpected(Failure{FailureCode::Unsupported, operation});
    }
    const auto page_bytes = static_cast<std::uint64_t>(page_size);
    const auto peak_kib = static_cast<std::uint64_t>(usage.ru_maxrss);
    if (resident_pages > std::numeric_limits<std::uint64_t>::max() / page_bytes ||
        peak_kib > std::numeric_limits<std::uint64_t>::max() / 1024) {
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
    }
    return ProcessMemory{static_cast<std::uint64_t>(resident_pages) * page_bytes, peak_kib * 1024};
#else
    return std::unexpected(Failure{FailureCode::Unsupported, operation});
#endif
}

Result<std::string> measurement_report_json(const MeasurementConditions& conditions,
                                            const StageRecorder& recorder,
                                            const ResourcePlan& resources, ProcessMemory memory,
                                            Operation operation) noexcept {
    if (conditions.machine.label.empty() || conditions.machine.cpu.empty() ||
        conditions.machine.logical_processors == 0 ||
        conditions.machine.physical_memory_bytes == 0 || conditions.case_id.empty() ||
        conditions.threads == 0) {
        return std::unexpected(Failure{FailureCode::InvalidInput, operation});
    }
    const auto times = recorder.snapshot();
    if (!times) {
        return std::unexpected(times.error());
    }
    return guard_entry(operation, FailureCode::InvalidInput, {}, [&]() -> Result<std::string> {
        auto settings = nlohmann::json::parse(conditions.settings_json, nullptr, false);
        if (!settings.is_object()) {
            return std::unexpected(Failure{FailureCode::InvalidInput, operation});
        }
        nlohmann::json stages = nlohmann::json::array();
        for (std::size_t i = 0; i < times->size(); ++i) {
            stages.push_back({{"stage", stage_names[i]},
                              {"calls", (*times)[i].calls},
                              {"elapsed_nanoseconds", (*times)[i].elapsed_nanoseconds}});
        }
        // N4-09: payload estimates and measured resident high-water are different quantities.
        const nlohmann::json report{
            {"schema", "uvlas.measurement/1"},
            {"machine",
             {{"label", conditions.machine.label},
              {"cpu", conditions.machine.cpu},
              {"logical_processors", conditions.machine.logical_processors},
              {"physical_memory_bytes", conditions.machine.physical_memory_bytes}}},
            {"case_id", conditions.case_id},
            {"threads", conditions.threads},
            {"settings", std::move(settings)},
            {"build", nlohmann::json::parse(uvlas::build_info_json())},
            {"timing_enabled", recorder.enabled()},
            {"stages", std::move(stages)},
            {"memory",
             {{"resident_bytes", memory.resident_bytes},
              {"peak_resident_bytes", memory.peak_resident_bytes},
              {"planned_payload_bytes", resources.allocated_bytes()},
              {"payload_limit_bytes", resources.byte_limit()}}}};
        return report.dump(2);
    });
}
} // namespace uvlas::core
