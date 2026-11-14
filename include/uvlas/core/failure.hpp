#pragma once

#include <cstdint>
#include <expected>
#include <functional>
#include <new>
#include <optional>
#include <type_traits>
#include <utility>
#include <uvlas/core/identifier.hpp>
#include <variant>

namespace uvlas::core {

// N4-08: names are the contract; ordinal values are not a file format.
enum class FailureCode : std::uint8_t {
    InvalidInput,
    Unsupported,
    GeometryInvalid,
    RankUnresolved,
    ChartSingular,
    Nonfinite,
    FactorizationFailed,
    LinearUnconverged,
    NonlinearUnconverged,
    EigenCoverageUnresolved,
    PrecisionUnsupported,
    PrecisionUnresolved,
    ResourceExceeded,
    Cancelled,
    StaleGeneration,
    CheckpointCorrupt,
    CheckpointIncompatible
};

enum class Operation : std::uint8_t {
    PrepareModel,
    EvaluateStructure,
    EvaluateAero,
    Initialize,
    EvaluateStep,
    SolveLinear,
    AdvanceTrial,
    Commit,
    LinearizeEquilibrium,
    SearchFlutter,
    Recover,
    CertifyError
};

enum class ResourceKind : std::uint8_t { Bytes, Iterations, OperatorApplications, Nanoseconds };

struct ResourceUsage {
    ResourceKind kind;
    std::uint64_t required;
    std::uint64_t observed;
    constexpr ResourceUsage(ResourceKind unit, std::uint64_t required_count,
                            std::uint64_t observed_count) noexcept
        : kind(unit), required(required_count), observed(observed_count) {}
    bool operator==(const ResourceUsage&) const = default;
};

using ProblemId =
    std::variant<DofId, PanelId, WakeVertexId, WakeFaceId, LayerId, ObservationId, ScopeId>;

struct FailureDetails {
    std::optional<ModelGeneration> model;
    std::optional<AcceptedStateId> parent;
    std::optional<TrialId> trial;
    std::optional<NumericalGeneration> numerical_generation;
    std::optional<ComponentId> component;
    std::optional<ProblemId> problem;
    std::optional<ResourceUsage> resources;
    std::optional<AcceptedStateId> last_accepted;
    bool operator==(const FailureDetails&) const = default;
};

struct Failure {
    FailureCode code;
    Operation operation;
    FailureDetails details;

    Failure(FailureCode failure_code, Operation process, FailureDetails context = {}) noexcept
        : code(failure_code), operation(process), details(context) {}
    bool operator==(const Failure&) const = default;
};

template <class T>
using Result = std::expected<T, Failure>;

namespace detail {
template <class T>
inline constexpr bool is_result = false;
template <class T>
inline constexpr bool is_result<Result<T>> = true;
} // namespace detail

// No allocating diagnostic is constructed while translating an allocation failure.
static_assert(std::is_nothrow_copy_constructible_v<Failure>);

template <class F>
    requires detail::is_result<std::invoke_result_t<F>>
[[nodiscard]] auto guard_entry(Operation operation, FailureCode exception_code,
                               const FailureDetails& details, F&& action) noexcept
    -> std::invoke_result_t<F> {
    try {
        return std::invoke(std::forward<F>(action));
    } catch (const std::bad_alloc&) {
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation, details});
    } catch (...) {
        // The entry chooses the operation-specific code; existing expected failures
        // returned by action are preserved, not reclassified by this guard.
        return std::unexpected(Failure{exception_code, operation, details});
    }
}

} // namespace uvlas::core
