#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <uvlas/core/failure.hpp>
#include <uvlas/core/resource_plan.hpp>

namespace uvlas::core {

/// Coarse computation phases; Count is a sentinel, not a measurable phase.
enum class MeasurementStage : std::uint8_t {
    Preparation,    ///< Geometry, topology and workspace preparation.
    Structure,      ///< Structural evaluation.
    Aerodynamics,   ///< Aerodynamic evaluation.
    LinearSolve,    ///< Linear solution.
    NonlinearSolve, ///< Nonlinear solution.
    TimeStep,       ///< One time step.
    FlutterSearch,  ///< Stability search.
    Checkpoint,     ///< State serialization and recovery.
    Count           ///< Number of supported stages.
};

/// Accumulated wall time for one phase; overlapping scopes count independently.
struct StageTiming {
    /// Number of recorded scope completions or manually supplied samples.
    std::uint64_t calls = 0;
    /// Sum of elapsed wall-clock durations in nanoseconds.
    std::uint64_t elapsed_nanoseconds = 0;
};

/// Fixed-size phase counters. No operation allocates, and disabled recording is a no-op.
/// Use a separate recorder per concurrent worker; do not share mutable counters.
class StageRecorder {
public:
    /// Enable or disable recording for this recorder's lifetime.
    explicit StageRecorder(bool enabled, Operation operation) noexcept
        : enabled_(enabled), operation_(operation) {}
    /// Whether scopes should read the clock and record samples.
    [[nodiscard]] bool enabled() const noexcept {
        return enabled_;
    }
    /// Add an externally measured duration; disabled recorders ignore the sample.
    /// Invalid phases return InvalidInput and counter overflow returns ResourceExceeded.
    /// The first failure is retained, so a failed sample cannot produce a valid report.
    [[nodiscard]] Result<void> record_sample(MeasurementStage stage,
                                             std::uint64_t elapsed_nanoseconds) noexcept;
    /// Return owned counters, or the first recording failure; no references are borrowed.
    [[nodiscard]] Result<std::array<StageTiming, static_cast<std::size_t>(MeasurementStage::Count)>>
    snapshot() const noexcept;

private:
    bool enabled_;
    Operation operation_;
    std::array<StageTiming, static_cast<std::size_t>(MeasurementStage::Count)> timings_{};
    std::optional<Failure> failure_;
};

/// Measure a phase's wall time until scope exit, including exception unwinding.
/// A disabled recorder causes no clock read or counter update. The recorder must outlive
/// this scope. Copying and moving are forbidden to prevent duplicate completion samples.
class StageScope {
public:
    /// Start a scope; stage must name a phase rather than the Count sentinel.
    StageScope(StageRecorder& recorder, MeasurementStage stage) noexcept;
    /// Record completion without allocating or throwing; inspect recorder.snapshot() for errors.
    ~StageScope() noexcept;
    /// A scope has exactly one owner.
    StageScope(const StageScope&) = delete;
    /// A scope cannot be reassigned.
    StageScope& operator=(const StageScope&) = delete;
    /// A scope cannot transfer its completion responsibility.
    StageScope(StageScope&&) = delete;
    /// A scope cannot be move-assigned.
    StageScope& operator=(StageScope&&) = delete;

private:
    StageRecorder& recorder_;
    MeasurementStage stage_;
    std::optional<std::chrono::steady_clock::time_point> start_;
};

/// Process resident memory, distinct from planned array payload and allocation capacity.
struct ProcessMemory {
    /// Resident bytes observed near the sample time.
    std::uint64_t resident_bytes;
    /// Process lifetime high-water resident bytes, not a phase-local peak.
    std::uint64_t peak_resident_bytes;
};

/// Sample resident memory on Windows or Linux; return Unsupported if OS sampling fails.
/// This function may perform OS I/O and must be called outside allocation-free kernels.
[[nodiscard]] Result<ProcessMemory> sample_process_memory(Operation operation) noexcept;

/// Explicit machine characteristics used to interpret measurements.
/// Strings borrow storage only while a report is serialized.
struct MachineDescription {
    /// Nonempty machine label without user names or local paths.
    std::string_view label;
    /// Nonempty processor model or architecture description.
    std::string_view cpu;
    /// Number of logical processors; must be positive.
    std::uint32_t logical_processors;
    /// Installed physical memory in bytes; must be positive.
    std::uint64_t physical_memory_bytes;
};

/// Caller-supplied reproducibility conditions; strings are borrowed only while reporting.
struct MeasurementConditions {
    /// Explicit machine characteristics, rather than a hostname used as a proxy.
    MachineDescription machine;
    /// Nonempty input/case identity.
    std::string_view case_id;
    /// Number of application workers used by the measured computation; must be positive.
    std::uint32_t threads;
    /// JSON object describing input sizes and solver settings, not an unstructured log.
    std::string_view settings_json;
};

/// Serialize owned JSON with conditions, build identity, phase times, process RSS and payload.
/// Invalid conditions return InvalidInput; recorder failures propagate and bad_alloc becomes
/// ResourceExceeded. Serialization allocates and belongs outside measured numerical scopes.
[[nodiscard]] Result<std::string> measurement_report_json(const MeasurementConditions& conditions,
                                                          const StageRecorder& recorder,
                                                          const ResourcePlan& resources,
                                                          ProcessMemory memory,
                                                          Operation operation) noexcept;

} // namespace uvlas::core
