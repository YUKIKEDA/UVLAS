#include "result_documents.hpp"

#include "source_digest.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <string>
#include <utility>
#include <uvlas/core/build_info.hpp>

namespace uvlas::io::detail {
namespace {

using Json = nlohmann::ordered_json;
using Code = core::FailureCode;
using Operation = core::Operation;
constexpr auto write = Operation::WriteResults;

[[nodiscard]] auto fail(Code code) noexcept {
    return std::unexpected(core::Failure{code, write});
}

#define UVLAS_RESULT_NAME(type, value)                                                             \
    std::pair {                                                                                    \
        type::value, std::string_view {                                                            \
            #value                                                                                 \
        }                                                                                          \
    }
constexpr std::array code_names{UVLAS_RESULT_NAME(Code, InvalidInput),
                                UVLAS_RESULT_NAME(Code, Unsupported),
                                UVLAS_RESULT_NAME(Code, GeometryInvalid),
                                UVLAS_RESULT_NAME(Code, RankUnresolved),
                                UVLAS_RESULT_NAME(Code, ChartSingular),
                                UVLAS_RESULT_NAME(Code, Nonfinite),
                                UVLAS_RESULT_NAME(Code, FactorizationFailed),
                                UVLAS_RESULT_NAME(Code, LinearUnconverged),
                                UVLAS_RESULT_NAME(Code, NonlinearUnconverged),
                                UVLAS_RESULT_NAME(Code, EigenCoverageUnresolved),
                                UVLAS_RESULT_NAME(Code, PrecisionUnsupported),
                                UVLAS_RESULT_NAME(Code, PrecisionUnresolved),
                                UVLAS_RESULT_NAME(Code, ResourceExceeded),
                                UVLAS_RESULT_NAME(Code, Cancelled),
                                UVLAS_RESULT_NAME(Code, StaleGeneration),
                                UVLAS_RESULT_NAME(Code, CheckpointCorrupt),
                                UVLAS_RESULT_NAME(Code, CheckpointIncompatible)};
constexpr std::array operation_names{UVLAS_RESULT_NAME(Operation, PrepareModel),
                                     UVLAS_RESULT_NAME(Operation, EvaluateStructure),
                                     UVLAS_RESULT_NAME(Operation, EvaluateAero),
                                     UVLAS_RESULT_NAME(Operation, Initialize),
                                     UVLAS_RESULT_NAME(Operation, EvaluateStep),
                                     UVLAS_RESULT_NAME(Operation, SolveLinear),
                                     UVLAS_RESULT_NAME(Operation, AdvanceTrial),
                                     UVLAS_RESULT_NAME(Operation, Commit),
                                     UVLAS_RESULT_NAME(Operation, LinearizeEquilibrium),
                                     UVLAS_RESULT_NAME(Operation, SearchFlutter),
                                     UVLAS_RESULT_NAME(Operation, Recover),
                                     UVLAS_RESULT_NAME(Operation, CertifyError),
                                     UVLAS_RESULT_NAME(Operation, ExtractEigen),
                                     UVLAS_RESULT_NAME(Operation, SaveCheckpoint),
                                     UVLAS_RESULT_NAME(Operation, WriteResults)};
using Resource = core::ResourceKind;
constexpr std::array resource_names{
    UVLAS_RESULT_NAME(Resource, Bytes), UVLAS_RESULT_NAME(Resource, Iterations),
    UVLAS_RESULT_NAME(Resource, OperatorApplications), UVLAS_RESULT_NAME(Resource, Nanoseconds)};
#undef UVLAS_RESULT_NAME

template <class Enum, std::size_t N>
[[nodiscard]] std::string_view
name_of(Enum value, const std::array<std::pair<Enum, std::string_view>, N>& names) noexcept {
    const auto entry = std::ranges::find(names, value, &std::pair<Enum, std::string_view>::first);
    return entry == names.end() ? std::string_view{} : entry->second;
}

template <class Id>
[[nodiscard]] std::string_view problem_domain(Id /*id*/) noexcept {
    if constexpr (std::is_same_v<Id, core::DofId>) {
        return "dof";
    } else if constexpr (std::is_same_v<Id, core::PanelId>) {
        return "panel";
    } else if constexpr (std::is_same_v<Id, core::WakeVertexId>) {
        return "wake_vertex";
    } else if constexpr (std::is_same_v<Id, core::WakeFaceId>) {
        return "wake_face";
    } else if constexpr (std::is_same_v<Id, core::LayerId>) {
        return "layer";
    } else if constexpr (std::is_same_v<Id, core::ObservationId>) {
        return "observation";
    } else {
        static_assert(std::is_same_v<Id, core::ScopeId>);
        return "scope";
    }
}

[[nodiscard]] bool known_exit(FailedResultExit exit) noexcept {
    switch (exit) {
    case FailedResultExit::Input:
    case FailedResultExit::Geometry:
    case FailedResultExit::Numerical:
    case FailedResultExit::Resources:
    case FailedResultExit::Persistence:
    case FailedResultExit::Interrupted:
        return true;
    }
    return false;
}

[[nodiscard]] std::optional<AcceptedResultIdentity>
snapshot(const core::AcceptedState* accepted) noexcept {
    if (accepted == nullptr) {
        return std::nullopt;
    }
    return AcceptedResultIdentity{accepted->id(), accepted->model().generation(),
                                  accepted->numerical_generation(),
                                  accepted->model().fingerprint()};
}

template <class Id>
[[nodiscard]] Json identifier(std::optional<Id> id) {
    if (id) {
        return std::to_string(id->value());
    }
    return nullptr;
}

[[nodiscard]] std::string hex_digest(const std::array<std::byte, 32>& bytes) {
    constexpr std::string_view hex = "0123456789abcdef";
    std::string text(2 * bytes.size(), '0');
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        const auto byte = std::to_integer<unsigned>(bytes[i]);
        text[2 * i] = hex[byte >> 4];
        text[2 * i + 1] = hex[byte & 15U];
    }
    return text;
}

[[nodiscard]] Json failure_json(const core::Failure& failure) {
    const auto& details = failure.details;
    Json result = {{"code", name_of(failure.code, code_names)},
                   {"operation", name_of(failure.operation, operation_names)},
                   {"details",
                    {{"model", identifier(details.model)},
                     {"parent", identifier(details.parent)},
                     {"trial", identifier(details.trial)},
                     {"numerical_generation", identifier(details.numerical_generation)},
                     {"component", identifier(details.component)},
                     {"problem", nullptr},
                     {"problem_domain", nullptr},
                     {"resources", nullptr},
                     {"last_accepted", identifier(details.last_accepted)}}}};
    auto& out = result["details"];
    if (details.problem) {
        std::visit(
            [&out](auto id) {
                out["problem"] = std::to_string(id.value());
                out["problem_domain"] = problem_domain(id);
            },
            *details.problem);
    }
    if (details.resources) {
        const auto& resources = *details.resources;
        out["resources"] = {{"kind", name_of(resources.kind, resource_names)},
                            {"required", std::to_string(resources.required)},
                            {"observed", std::to_string(resources.observed)}};
    }
    return result;
}

[[nodiscard]] Json last_accepted(const ResultStatus& status,
                                 const std::optional<core::Failure>& failure) {
    if (const auto accepted = status.accepted()) {
        return {{"id", std::to_string(accepted->id.value())},
                {"model_generation", std::to_string(accepted->model_generation.value())},
                {"numerical_generation", std::to_string(accepted->numerical_generation.value())},
                {"model_sha256", hex_digest(accepted->model.bytes)}};
    }
    if (failure && failure->details.last_accepted) {
        return {{"id", identifier(failure->details.last_accepted)}};
    }
    return nullptr;
}

[[nodiscard]] std::string provenance_json(const ResultProvenance& provenance) {
    Json result = {{"build", Json::parse(uvlas::build_info_json())},
                   {"original_case", nullptr},
                   {"resolved_case_sha256", nullptr}};
    if (provenance.original_case) {
        const auto& source = *provenance.original_case;
        result["original_case"] = {{"source_id", std::to_string(source.source.value())},
                                   {"byte_count", std::to_string(source.byte_count)},
                                   {"sha256", hex_digest(source.digest.bytes)}};
    }
    if (provenance.resolved_case) {
        result["resolved_case_sha256"] = hex_digest(provenance.resolved_case->bytes);
    }
    return result.dump() + '\n';
}

[[nodiscard]] std::string status_json(const ResultStatus& status, SourceDigest provenance) {
    const auto failure = status.failure();
    const Json result = {{"completed", !failure},
                         {"exit_code", status.exit_code()},
                         {"failure", failure ? failure_json(*failure) : Json(nullptr)},
                         {"last_accepted", last_accepted(status, failure)},
                         {"provenance_sha256", hex_digest(provenance.bytes)}};
    return result.dump() + '\n';
}

} // namespace

ResultStatus::ResultStatus(std::uint8_t exit_code, std::optional<core::Failure> failure,
                           std::optional<AcceptedResultIdentity> accepted) noexcept
    : exit_code_(exit_code), failure_(failure), accepted_(accepted) {}

core::Result<ResultStatus> ResultStatus::completed(const core::AcceptedState* last) noexcept {
    if (last != nullptr && !last->active()) {
        return fail(Code::StaleGeneration);
    }
    return ResultStatus{0, std::nullopt, snapshot(last)};
}

core::Result<ResultStatus> ResultStatus::failed(core::Failure failure, FailedResultExit exit,
                                                const core::AcceptedState* last) noexcept {
    const auto& details = failure.details;
    if (!known_exit(exit) || name_of(failure.code, code_names).empty() ||
        name_of(failure.operation, operation_names).empty() ||
        (details.resources && name_of(details.resources->kind, resource_names).empty()) ||
        ((failure.code == Code::Cancelled) != (exit == FailedResultExit::Interrupted)) ||
        ((failure.code == Code::ResourceExceeded) != (exit == FailedResultExit::Resources))) {
        return fail(Code::InvalidInput);
    }
    if (last != nullptr &&
        (!last->active() || (details.last_accepted && *details.last_accepted != last->id()) ||
         (details.model && *details.model != last->model().generation()))) {
        return fail(Code::StaleGeneration);
    }
    return ResultStatus{static_cast<std::uint8_t>(exit), failure, snapshot(last)};
}

ResultDocuments::ResultDocuments(core::PlannedBuffer<char> bytes, std::size_t status_bytes) noexcept
    : bytes_(std::move(bytes)), status_bytes_(status_bytes) {}

std::string_view ResultDocuments::status() const& noexcept {
    const auto data = bytes_.values();
    return data.empty() ? std::string_view{} : std::string_view{data.data(), status_bytes_};
}

std::string_view ResultDocuments::provenance() const& noexcept {
    const auto data = bytes_.values();
    return data.empty()
               ? std::string_view{}
               : std::string_view{data.data() + status_bytes_, data.size() - status_bytes_};
}

core::Result<ResultDocuments> write_result_documents(const ResultStatus& status,
                                                     const ResultProvenance& provenance,
                                                     std::uint64_t maximum_output_bytes,
                                                     core::ResourcePlan& plan,
                                                     const std::stop_token& stop) noexcept {
    return core::guard_entry(write, Code::InvalidInput, {}, [&]() -> core::Result<ResultDocuments> {
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        if (provenance.resolved_case && !provenance.original_case) {
            return fail(Code::InvalidInput);
        }
        // Both temporaries contain fixed-size metadata, never physical arrays or user text.
        const auto source = provenance_json(provenance);
        const auto encoded = status_json(status, source_digest(std::as_bytes(std::span{source})));
        const auto budget = std::min(maximum_output_bytes, plan.remaining_bytes());
        if (source.size() > budget || encoded.size() > budget - source.size()) {
            return fail(Code::ResourceExceeded);
        }
        const auto count = source.size() + encoded.size();
        if (!std::in_range<core::Index>(count)) {
            return fail(Code::ResourceExceeded);
        }
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        auto buffer = plan.allocate<char>(static_cast<core::Index>(count), write);
        if (!buffer) {
            return std::unexpected(buffer.error());
        }
        auto cursor = std::ranges::copy(encoded, buffer->values().begin()).out;
        std::ranges::copy(source, cursor);
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        return ResultDocuments{std::move(*buffer), encoded.size()};
    });
}

} // namespace uvlas::io::detail
