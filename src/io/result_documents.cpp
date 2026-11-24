#include "result_documents.hpp"

#include "json_object_reader.hpp"
#include "source_digest.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <charconv>
#include <string>
#include <type_traits>
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

[[nodiscard]] bool matching_exit(Code code, FailedResultExit exit) noexcept {
    return known_exit(exit) &&
           ((code == Code::Cancelled) == (exit == FailedResultExit::Interrupted)) &&
           ((code == Code::ResourceExceeded) == (exit == FailedResultExit::Resources));
}

[[nodiscard]] bool matching_accepted(const core::FailureDetails& details, core::AcceptedStateId id,
                                     core::ModelGeneration model) noexcept {
    return (!details.last_accepted || *details.last_accepted == id) &&
           (!details.model || *details.model == model);
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

constexpr auto read = Operation::Recover;

[[noreturn]] void bad_document(Code code = Code::InvalidInput, core::FailureDetails details = {}) {
    throw core::Failure{code, read, details};
}

void require_field(bool matches) {
    if (!matches)
        bad_document();
}

template <class T>
[[nodiscard]] T field_value(std::optional<T> value) {
    if (!value)
        bad_document();
    return std::move(*value);
}

[[nodiscard]] InputNodeView member(const InputObjectView& object, std::string_view key) {
    return field_value(object.find(key));
}

[[nodiscard]] InputObjectView fields(const InputNodeView& node, core::Index count) {
    const auto object = field_value(node.object());
    require_field(object.size() == count);
    return object;
}

[[nodiscard]] std::string_view string_value(const InputNodeView& node) {
    return field_value(node.string());
}

[[nodiscard]] std::uint64_t identifier_number(const InputNodeView& node) {
    const auto token = string_value(node);
    require_field(!token.empty() && token.size() <= 20 &&
                  (token.size() == 1 || token.front() != '0'));
    std::uint64_t result = 0;
    const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), result);
    require_field(error == std::errc{} && end == token.data() + token.size());
    return result;
}

template <class Id>
[[nodiscard]] std::optional<Id> optional_id(const InputObjectView& object, std::string_view key) {
    const auto node = member(object, key);
    if (node.kind() == InputNodeKind::Null)
        return std::nullopt;
    return Id{identifier_number(node)};
}

[[nodiscard]] SourceDigest declared_digest(const InputNodeView& node) {
    const auto token = string_value(node);
    require_field(token.size() == 64);
    constexpr std::string_view hex = "0123456789abcdef";
    SourceDigest value{};
    for (std::size_t i = 0; i < value.bytes.size(); ++i) {
        const auto high = hex.find(token[2 * i]);
        const auto low = hex.find(token[2 * i + 1]);
        require_field(high != hex.npos && low != hex.npos);
        value.bytes[i] = static_cast<std::byte>((high << 4) | low);
    }
    return value;
}

template <class Enum, std::size_t N>
[[nodiscard]] Enum enum_value(const InputNodeView& node,
                              const std::array<std::pair<Enum, std::string_view>, N>& names) {
    const auto name = string_value(node);
    const auto entry = std::ranges::find(names, name, &std::pair<Enum, std::string_view>::second);
    if (entry == names.end())
        bad_document();
    return entry->first;
}

template <class... Id>
[[nodiscard]] core::ProblemId problem_value(std::string_view name, std::uint64_t value,
                                            std::type_identity<std::variant<Id...>>) {
    std::optional<core::ProblemId> result;
    const auto match = [&]<class Identifier>(Identifier id) {
        if (name == problem_domain(id))
            result = id;
    };
    (match(Id{value}), ...);
    return field_value(result);
}

[[nodiscard]] core::Failure read_failure(const InputNodeView& node) {
    const auto object = fields(node, 3);
    const auto code = enum_value(member(object, "code"), code_names);
    const auto operation = enum_value(member(object, "operation"), operation_names);
    const auto details = fields(member(object, "details"), 9);
    core::FailureDetails value{};
    value.model = optional_id<core::ModelGeneration>(details, "model");
    value.parent = optional_id<core::AcceptedStateId>(details, "parent");
    value.trial = optional_id<core::TrialId>(details, "trial");
    value.numerical_generation =
        optional_id<core::NumericalGeneration>(details, "numerical_generation");
    value.component = optional_id<core::ComponentId>(details, "component");
    value.last_accepted = optional_id<core::AcceptedStateId>(details, "last_accepted");
    const auto problem = member(details, "problem");
    const auto domain = member(details, "problem_domain");
    if (problem.kind() != InputNodeKind::Null) {
        value.problem = problem_value(string_value(domain), identifier_number(problem),
                                      std::type_identity<core::ProblemId>{});
    } else {
        require_field(domain.kind() == InputNodeKind::Null);
    }
    const auto resource = member(details, "resources");
    if (resource.kind() != InputNodeKind::Null) {
        const auto usage = fields(resource, 3);
        value.resources.emplace(enum_value(member(usage, "kind"), resource_names),
                                identifier_number(member(usage, "required")),
                                identifier_number(member(usage, "observed")));
    }
    return {code, operation, value};
}

struct DecodedStatus {
    InspectedResultStatus status;
    SourceDigest digest;
};

[[nodiscard]] DecodedStatus decode_status(const InputNodeView& node) {
    const auto object = fields(node, 5);
    const auto exit_node = member(object, "exit_code");
    const auto exit = field_value(exit_node.unsigned_integer());
    require_field(std::in_range<std::uint8_t>(exit));
    InspectedResultStatus result{static_cast<std::uint8_t>(exit), {}, {}, {}};
    const auto complete_node = member(object, "completed");
    const auto completed = field_value(complete_node.boolean());
    const auto failure = member(object, "failure");
    if (failure.kind() != InputNodeKind::Null)
        result.failure = read_failure(failure);
    require_field(completed == (exit == 0) && completed == !result.failure);
    if (const auto failed = result.failure) {
        require_field(matching_exit(failed->code, static_cast<FailedResultExit>(exit)));
    }

    const auto last = member(object, "last_accepted");
    if (last.kind() != InputNodeKind::Null) {
        const auto accepted = field_value(last.object());
        require_field(accepted.size() == 1 || accepted.size() == 4);
        const core::AcceptedStateId id{identifier_number(member(accepted, "id"))};
        result.last_accepted_id = id;
        if (accepted.size() == 1) {
            require_field(result.failure && result.failure->details.last_accepted == id);
        } else {
            const core::ModelGeneration model{
                identifier_number(member(accepted, "model_generation"))};
            const core::NumericalGeneration numerical{
                identifier_number(member(accepted, "numerical_generation"))};
            const auto model_digest = declared_digest(member(accepted, "model_sha256"));
            result.accepted = AcceptedResultIdentity{id, model, numerical,
                                                     core::ModelFingerprint{model_digest.bytes}};
            if (const auto failed = result.failure)
                require_field(matching_accepted(failed->details, id, model));
        }
    } else if (const auto failed = result.failure) {
        require_field(!failed->details.last_accepted);
    }
    return {result, declared_digest(member(object, "provenance_sha256"))};
}

void validate_build(const InputNodeView& node) {
    const auto object = field_value(node.object());
    if (string_value(member(object, "schema")) != "uvlas.build-info/1")
        bad_document(Code::Unsupported);
    require_field(object.size() == 12);
    for (const auto name : {"version", "commit", "compiler", "platform", "configuration",
                            "floating_point", "instruction_set", "platform_key"})
        static_cast<void>(string_value(member(object, name)));
    const auto dirty = member(object, "uncommitted_changes");
    require_field(dirty.boolean().has_value());
    const auto dependencies = fields(member(object, "dependencies"), 4);
    for (const auto name : {"eigen", "spectra", "boost", "nlohmann_json"})
        static_cast<void>(string_value(member(dependencies, name)));
    const auto repairs = fields(member(object, "dependency_repairs"), 1);
    const auto allocation = member(repairs, "nlohmann_json_sax_allocation");
    require_field(allocation.unsigned_integer().has_value());
}

[[nodiscard]] ResultProvenance decode_provenance(const InputNodeView& node) {
    const auto object = fields(node, 3);
    validate_build(member(object, "build"));
    ResultProvenance result;
    const auto original = member(object, "original_case");
    if (original.kind() != InputNodeKind::Null) {
        const auto source = fields(original, 3);
        result.original_case =
            SourceIdentity{SourceId{identifier_number(member(source, "source_id"))},
                           identifier_number(member(source, "byte_count")),
                           declared_digest(member(source, "sha256"))};
    }
    const auto resolved = member(object, "resolved_case_sha256");
    if (resolved.kind() != InputNodeKind::Null) {
        require_field(result.original_case.has_value());
        result.resolved_case = declared_digest(resolved);
    }
    return result;
}

void poll_read(const std::stop_token& stop) {
    if (stop.stop_requested())
        bad_document(Code::Cancelled);
}

[[nodiscard]] InputDocument read_object(std::span<const std::byte> bytes, SourceId id,
                                        const InputLimits& limits, core::ResourcePlan& plan,
                                        const std::stop_token& stop) {
    poll_read(stop);
    std::optional<InputDiagnostic> diagnostic;
    auto result = read_json_object(bytes, id, limits, plan, diagnostic);
    poll_read(stop);
    if (!result) {
        const auto& error = result.error();
        bad_document(error.code == Code::Nonfinite ? Code::InvalidInput : error.code,
                     error.details);
    }
    return std::move(*result);
}

[[nodiscard]] InputDocumentView parsed_view(const InputDocument& document) {
    const auto view = document.view();
    if (!view)
        bad_document(view.error().code, view.error().details);
    return *view;
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
    if (!matching_exit(failure.code, exit) || name_of(failure.code, code_names).empty() ||
        name_of(failure.operation, operation_names).empty() ||
        (details.resources && name_of(details.resources->kind, resource_names).empty())) {
        return fail(Code::InvalidInput);
    }
    if (last != nullptr &&
        (!last->active() || !matching_accepted(details, last->id(), last->model().generation()))) {
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

InspectedResultDocuments::InspectedResultDocuments(InspectedResultStatus status,
                                                   ResultProvenance provenance,
                                                   InputDocument build_owner) noexcept
    : status_(status), provenance_(provenance), build_owner_(std::move(build_owner)) {}

std::optional<InputObjectView> InspectedResultDocuments::build_info() const& noexcept {
    const auto doc = build_owner_.view();
    if (!doc)
        return std::nullopt;
    const auto envelope = doc->root.object();
    if (!envelope)
        return std::nullopt;
    const auto build = envelope->find("build");
    return build ? build->object() : std::nullopt;
}

core::Result<InspectedResultDocuments> read_result_documents(std::span<const std::byte> status,
                                                             std::span<const std::byte> provenance,
                                                             const InputLimits& limits,
                                                             core::ResourcePlan& plan,
                                                             const std::stop_token& stop) noexcept {
    return core::guard_entry(
        read, Code::InvalidInput, {}, [&]() -> core::Result<InspectedResultDocuments> {
            try {
                const auto checked = [&] {
                    const auto status_document =
                        read_object(status, SourceId{0}, limits, plan, stop);
                    return decode_status(parsed_view(status_document).root);
                }();
                auto provenance_document = read_object(provenance, SourceId{1}, limits, plan, stop);
                const auto view = parsed_view(provenance_document);
                require_field(view.identity.digest == checked.digest);
                const auto checked_provenance = decode_provenance(view.root);
                poll_read(stop);
                return InspectedResultDocuments{checked.status, checked_provenance,
                                                std::move(provenance_document)};
            } catch (const core::Failure& error) {
                return std::unexpected(error);
            }
        });
}
} // namespace uvlas::io::detail
