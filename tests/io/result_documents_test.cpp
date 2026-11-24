#include "result_documents.hpp"
#include "source_digest.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <uvlas/core/build_info.hpp>

namespace {

using namespace uvlas::core;
using namespace uvlas::io;
using namespace uvlas::io::detail;
using Json = nlohmann::json;

constexpr auto maximum = std::numeric_limits<std::uint64_t>::max();
constexpr auto largest = "18446744073709551615";
constexpr std::uint64_t capacity = 8192;
constexpr std::uint64_t parser_budget = 1U << 20U;
const InputLimits read_limits{
    capacity, 12, 250,      250,
    512,      64, capacity, {ParserAllocationMode::EstimateAdmission, parser_budget}};

[[nodiscard]] std::span<const std::byte> bytes(std::string_view text) noexcept {
    return std::as_bytes(std::span{text});
}

template <class T>
T take(Result<T> value) {
    REQUIRE(value);
    if (!value) {
        throw std::logic_error("Expected a result value");
    }
    return std::move(*value);
}

AcceptedState accept(ResourcePlan& plan, ModelGeneration generation = ModelGeneration{maximum}) {
    const AcceptancePlan requirements{{{ObservationId{1}, ScopeId{20}, 1, Precision::Binary64},
                                       {ObservationId{2}, ScopeId{20}, 1, Precision::Binary64},
                                       {ObservationId{3}, ScopeId{20}, 1, Precision::Binary64},
                                       {ObservationId{4}, ScopeId{20}, 1, Precision::Binary64}}};
    const auto model = take(Model::create({}, generation, requirements, plan));
    auto trial = take(TrialState::initial(model, TrialId{1}, NumericalGeneration{0}, {}, plan));
    const auto assessment = [&](std::size_t i) {
        const auto& requirement = requirements[i];
        return take(ErrorCertificate::assess({requirement.observation, requirement.scope,
                                              NumericalGeneration{0}, requirement.representation},
                                             requirement.budget, 0.25));
    };
    const auto report = take(EvaluationReport::create(
        model, trial.target(), {assessment(0), assessment(1), assessment(2), assessment(3)}));
    auto packet = take(qualify(trial, report));
    return take(commit_initial(trial, packet, AcceptedStateId{maximum}));
}

std::string hex(const auto& bytes) {
    constexpr std::string_view digits = "0123456789abcdef";
    std::string s;
    for (const auto value : bytes) {
        const auto byte = std::to_integer<unsigned>(value);
        s.push_back(digits[byte / 16]);
        s.push_back(digits[byte % 16]);
    }
    return s;
}

TEST_CASE("Owned result documents preserve a real accepted state and link provenance",
          "[IO4-05][QA-06]") {
    ResourcePlan state_plan{capacity};
    const auto status = [&] {
        const auto accepted = accept(state_plan);
        return take(ResultStatus::completed(&accepted));
    }();
    const auto accepted = status.accepted();
    REQUIRE(accepted);
    SourceDigest original{};
    original.bytes.front() = std::byte{0xab};
    SourceDigest resolved{};
    resolved.bytes.back() = std::byte{0xcd};
    ResourcePlan output_plan{capacity};
    auto documents = take(write_result_documents(
        status, {SourceIdentity{SourceId{maximum}, maximum, original}, resolved}, capacity,
        output_plan, {}));
    auto moved = std::move(documents);
    // A moved-from owner expressly provides empty views.
    // NOLINTNEXTLINE(bugprone-use-after-move)
    CHECK(documents.status().empty());
    CHECK(documents.provenance().empty());
    const auto run = Json::parse(moved.status());
    const auto provenance = Json::parse(moved.provenance());
    CHECK(run.at("completed") == true);
    CHECK(run.at("exit_code") == 0);
    CHECK(run.at("failure").is_null());
    CHECK(run.at("last_accepted").at("id") == largest);
    CHECK(run.at("last_accepted").at("model_generation") == largest);
    CHECK(run.at("last_accepted").at("numerical_generation") == "0");
    if (accepted) {
        CHECK(run.at("last_accepted").at("model_sha256") == hex(accepted->model.bytes));
    }
    CHECK(provenance.at("build") == Json::parse(uvlas::build_info_json()));
    CHECK(provenance.at("original_case").at("source_id") == largest);
    CHECK(provenance.at("original_case").at("byte_count") == largest);
    CHECK(provenance.at("original_case").at("sha256") == ("ab" + std::string(62, '0')));
    CHECK(provenance.at("resolved_case_sha256") == std::string(62, '0') + "cd");
    CHECK(run.at("provenance_sha256") ==
          hex(source_digest(std::as_bytes(std::span{moved.provenance()})).bytes));
    CHECK(output_plan.allocated_bytes() == moved.status().size() + moved.provenance().size());
}

TEST_CASE("Unconverged results retain failure details without reporting completion",
          "[IO4-05][N4-08]") {
    const Failure failure{FailureCode::NonlinearUnconverged,
                          Operation::AdvanceTrial,
                          {.model = ModelGeneration{0},
                           .parent = AcceptedStateId{maximum - 1},
                           .trial = TrialId{maximum},
                           .numerical_generation = NumericalGeneration{2},
                           .component = ComponentId{3},
                           .problem = WakeFaceId{maximum},
                           .resources = ResourceUsage{ResourceKind::Iterations, maximum, 0},
                           .last_accepted = AcceptedStateId{0}}};
    const auto status = take(ResultStatus::failed(failure, FailedResultExit::Numerical));
    CHECK(status.failure() == failure);
    ResourcePlan plan{capacity};
    const auto documents = take(write_result_documents(status, {}, capacity, plan, {}));
    const auto run = Json::parse(documents.status());
    const auto provenance = Json::parse(documents.provenance());
    CHECK(run.at("completed") == false);
    CHECK(run.at("exit_code") == 4);
    CHECK(run.at("last_accepted") == Json::parse(R"({"id":"0"})"));
    CHECK(run.at("failure") == Json::parse(R"({
        "code":"NonlinearUnconverged","operation":"AdvanceTrial",
        "details":{"model":"0","parent":"18446744073709551614","trial":"18446744073709551615",
        "numerical_generation":"2","component":"3","problem":"18446744073709551615",
        "problem_domain":"wake_face","resources":{"kind":"Iterations",
        "required":"18446744073709551615","observed":"0"},"last_accepted":"0"}
    })"));
    CHECK(provenance.at("original_case").is_null());
    CHECK(provenance.at("resolved_case_sha256").is_null());
}

TEST_CASE("Result statuses reject invalid termination and moved or inconsistent parents",
          "[IO4-05][N4-08]") {
    const Failure cancelled{FailureCode::Cancelled, Operation::EvaluateStep};
    const Failure resource{FailureCode::ResourceExceeded, Operation::WriteResults};
    CHECK_FALSE(ResultStatus::failed(cancelled, FailedResultExit::Numerical));
    CHECK_FALSE(ResultStatus::failed(resource, FailedResultExit::Input));
    CHECK_FALSE(ResultStatus::failed(cancelled, static_cast<FailedResultExit>(0)));

    ResourcePlan plan{capacity};
    const auto cancelled_status =
        take(ResultStatus::failed(cancelled, FailedResultExit::Interrupted));
    const auto interrupted = take(write_result_documents(cancelled_status, {}, capacity, plan, {}));
    const auto run = Json::parse(interrupted.status());
    CHECK(run.at("exit_code") == 130);
    CHECK(run.at("failure").at("code") == "Cancelled");
    CHECK(run.at("last_accepted").is_null());
    CHECK(run.at("failure").at("details").at("last_accepted").is_null());
    const auto limited_status = take(ResultStatus::failed(resource, FailedResultExit::Resources));
    const auto limited = take(write_result_documents(limited_status, {}, capacity, plan, {}));
    CHECK(Json::parse(limited.status()).at("exit_code") == 5);
    CHECK(Json::parse(limited.status()).at("failure").at("operation") == "WriteResults");

    auto old = accept(plan);
    const auto accepted = std::move(old);
    // AcceptedState expressly exposes an inactive moved-from value.
    // NOLINTNEXTLINE(bugprone-use-after-move)
    CHECK_FALSE(ResultStatus::completed(&old));
    FailureDetails mismatched_context{};
    mismatched_context.last_accepted = AcceptedStateId{0};
    const Failure mismatched{FailureCode::Nonfinite, Operation::EvaluateStep, mismatched_context};
    const auto stale = ResultStatus::failed(mismatched, FailedResultExit::Numerical, &accepted);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == FailureCode::StaleGeneration);
    CHECK(stale.error().operation == Operation::WriteResults);

    const auto other = accept(plan, ModelGeneration{maximum - 1});
    REQUIRE(other.id() == accepted.id());
    FailureDetails same_id_context{};
    same_id_context.model = ModelGeneration{maximum};
    same_id_context.last_accepted = accepted.id();
    const Failure same_id{FailureCode::Nonfinite, Operation::EvaluateStep, same_id_context};
    CHECK(ResultStatus::failed(same_id, FailedResultExit::Numerical, &accepted));
    const auto stale_model = ResultStatus::failed(same_id, FailedResultExit::Numerical, &other);
    REQUIRE_FALSE(stale_model);
    CHECK(stale_model.error().code == FailureCode::StaleGeneration);
    CHECK(stale_model.error().operation == Operation::WriteResults);

    const Failure nonfinite{FailureCode::Nonfinite, Operation::EvaluateStep};
    CHECK(take(ResultStatus::failed(nonfinite, FailedResultExit::Numerical)).exit_code() == 4);
    CHECK(take(ResultStatus::failed(nonfinite, FailedResultExit::Input)).exit_code() == 2);
}

TEST_CASE("Result document admission is exact and failure delivers no partial bytes", "[IO4-05]") {
    const auto status = take(ResultStatus::completed());
    ResourcePlan first{capacity};
    const auto reference = take(write_result_documents(status, {}, capacity, first, {}));
    const auto size = first.allocated_bytes();
    ResourcePlan exact{size};
    const auto same = take(write_result_documents(status, {}, size, exact, {}));
    CHECK(reference.status() == same.status());
    CHECK(reference.provenance() == same.provenance());

    ResourcePlan short_plan{size - 1};
    const auto over_plan = write_result_documents(status, {}, size, short_plan, {});
    REQUIRE_FALSE(over_plan);
    CHECK(over_plan.error().code == FailureCode::ResourceExceeded);
    CHECK(over_plan.error().operation == Operation::WriteResults);
    CHECK(short_plan.allocated_bytes() == 0);
    ResourcePlan short_output{size};
    const auto over_output = write_result_documents(status, {}, size - 1, short_output, {});
    REQUIRE_FALSE(over_output);
    CHECK(over_output.error().code == FailureCode::ResourceExceeded);
    CHECK(short_output.allocated_bytes() == 0);
    const auto inconsistent =
        write_result_documents(status, {std::nullopt, SourceDigest{}}, size, short_output, {});
    REQUIRE_FALSE(inconsistent);
    CHECK(inconsistent.error().code == FailureCode::InvalidInput);
    // MSVC request_stop is non-const.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source cancel;
    static_cast<void>(cancel.request_stop());
    const auto cancelled =
        write_result_documents(status, {}, size, short_output, cancel.get_token());
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == FailureCode::Cancelled);
    CHECK(short_output.allocated_bytes() == 0);
}

TEST_CASE("The adopted failed exit categories are emitted as numbers", "[IO4-05]") {
    struct Entry {
        FailureCode code;
        FailedResultExit exit;
        int number;
    };
    constexpr std::array entries{
        Entry{FailureCode::InvalidInput, FailedResultExit::Input, 2},
        Entry{FailureCode::GeometryInvalid, FailedResultExit::Geometry, 3},
        Entry{FailureCode::LinearUnconverged, FailedResultExit::Numerical, 4},
        Entry{FailureCode::ResourceExceeded, FailedResultExit::Resources, 5},
        Entry{FailureCode::CheckpointIncompatible, FailedResultExit::Persistence, 6},
        Entry{FailureCode::Cancelled, FailedResultExit::Interrupted, 130}};
    for (const auto& [code, exit, number] : entries) {
        ResourcePlan plan{capacity};
        const auto status = take(ResultStatus::failed({code, Operation::PrepareModel}, exit));
        const auto documents = take(write_result_documents(status, {}, capacity, plan, {}));
        const auto run = Json::parse(documents.status());
        CHECK(run.at("exit_code") == number);
        CHECK(run.at("completed") == false);
        ResourcePlan reading{parser_budget};
        const auto recovered = take(read_result_documents(
            bytes(documents.status()), bytes(documents.provenance()), read_limits, reading));
        CHECK(recovered.status().exit_code == number);
        CHECK(recovered.status().failure == Failure{code, Operation::PrepareModel});
    }
}

TEST_CASE("Result inspection owns historical outcome, build and case identities",
          "[IO4-05][QA-06][verification]") {
    SourceDigest original{};
    original.bytes.front() = std::byte{0xab};
    const ResultProvenance expected{SourceIdentity{SourceId{maximum}, maximum, original},
                                    SourceDigest{}};
    FailureDetails context{};
    context.last_accepted = AcceptedStateId{maximum};
    context.model = ModelGeneration{maximum};
    context.trial = TrialId{0};
    context.problem = ObservationId{maximum};
    context.resources = ResourceUsage{ResourceKind::OperatorApplications, maximum, 0};
    const Failure report{FailureCode::NonlinearUnconverged, Operation::AdvanceTrial, context};
    const auto inspected = [&] {
        ResourcePlan state_plan{capacity};
        const auto accepted = accept(state_plan);
        const auto status =
            take(ResultStatus::failed(report, FailedResultExit::Numerical, &accepted));
        ResourcePlan plan{parser_budget};
        const auto output = take(write_result_documents(status, expected, capacity, plan, {}));
        return take(read_result_documents(bytes(output.status()), bytes(output.provenance()),
                                          read_limits, plan));
    }();
    const auto status = inspected.status();
    CHECK(status.exit_code == 4);
    CHECK(status.failure == report);
    CHECK(status.last_accepted_id == AcceptedStateId{maximum});
    REQUIRE(status.accepted);
    if (const auto accepted = status.accepted) {
        CHECK(accepted->id == AcceptedStateId{maximum});
        CHECK(accepted->model_generation == ModelGeneration{maximum});
        CHECK(accepted->numerical_generation == NumericalGeneration{0});
        ResourcePlan state_plan{capacity};
        const auto state = accept(state_plan);
        CHECK(accepted->model == state.model().fingerprint());
    }
    const auto provenance = inspected.provenance();
    CHECK(provenance.original_case == expected.original_case);
    CHECK(provenance.resolved_case == expected.resolved_case);
    const auto build = inspected.build_info();
    REQUIRE(build);
    if (build) {
        const auto schema = build->find("schema");
        if (schema) {
            CHECK(schema->string() == "uvlas.build-info/1");
        } else {
            FAIL("Missing build schema");
        }
        const auto version = build->find("version");
        if (version) {
            CHECK(version->string() == uvlas::build_info().version);
        } else {
            FAIL("Missing historical build version");
        }
    }

    ResourcePlan plan{parser_budget};
    context.last_accepted = AcceptedStateId{0};
    const auto id_only = take(ResultStatus::failed({report.code, report.operation, context},
                                                   FailedResultExit::Numerical));
    const auto output = take(write_result_documents(id_only, {}, capacity, plan, {}));
    const auto id_read = take(read_result_documents(bytes(output.status()),
                                                    bytes(output.provenance()), read_limits, plan));
    CHECK(id_read.status().last_accepted_id == AcceptedStateId{0});
    CHECK_FALSE(id_read.status().accepted);

    auto earlier_provenance = Json::parse(output.provenance());
    earlier_provenance["build"]["commit"] = std::string(40, 'a');
    const auto earlier_text = earlier_provenance.dump();
    auto linked = Json::parse(output.status());
    linked["provenance_sha256"] = hex(source_digest(bytes(earlier_text)).bytes);
    const auto earlier =
        take(read_result_documents(bytes(linked.dump()), bytes(earlier_text), read_limits, plan));
    const auto earlier_build = earlier.build_info();
    REQUIRE(earlier_build);
    if (earlier_build) {
        const auto commit = earlier_build->find("commit");
        if (commit) {
            CHECK(commit->string() == std::string(40, 'a'));
        } else {
            FAIL("Missing older build identity");
        }
    }
}

TEST_CASE("Result inspection rejects inconsistent content even under a matching digest",
          "[IO4-05][QA-06][negative][verification]") {
    ResourcePlan writer{capacity};
    const auto accepted = accept(writer);
    FailureDetails context{};
    context.model = accepted.model().generation();
    context.last_accepted = accepted.id();
    const auto outcome =
        take(ResultStatus::failed({FailureCode::LinearUnconverged, Operation::SolveLinear, context},
                                  FailedResultExit::Numerical, &accepted));
    const auto documents = take(write_result_documents(outcome, {}, capacity, writer, {}));
    auto status = Json::parse(documents.status());
    auto provenance = Json::parse(documents.provenance());
    auto expected = FailureCode::InvalidInput;
    auto consistent_digest = true;
    auto duplicate_status = false;
    SECTION("edited provenance bytes") {
        provenance["build"]["version"] = "older";
        consistent_digest = false;
    }
    SECTION("duplicate decoded status key") {
        duplicate_status = true;
    }
    SECTION("completed contradicts the recorded failure") {
        status["completed"] = true;
    }
    SECTION("exit must be an integer rather than a string") {
        status["exit_code"] = "4";
    }
    SECTION("exit must agree with a failure category") {
        status["exit_code"] = 130;
    }
    SECTION("model generation must match the reported accepted state") {
        status["last_accepted"]["model_generation"] = "0";
    }
    SECTION("distinct accepted IDs cannot agree through binary rounding") {
        status["last_accepted"]["id"] = "18446744073709551614";
    }
    SECTION("an ID must be a canonical uint64 string") {
        status["failure"]["details"]["parent"] = "00";
    }
    SECTION("an overflowing ID is unrepresentable") {
        status["failure"]["details"]["parent"] = "18446744073709551616";
    }
    SECTION("a failure's present last ID cannot become absent") {
        status["last_accepted"] = nullptr;
    }
    SECTION("a problem domain requires an id") {
        status["failure"]["details"]["problem_domain"] = "panel";
    }
    SECTION("a provenance digest has exactly 64 lowercase characters") {
        provenance["original_case"] = {
            {"source_id", "0"}, {"byte_count", "0"}, {"sha256", std::string(64, 'A')}};
    }
    SECTION("resolved digest requires an original source") {
        provenance["resolved_case_sha256"] = std::string(64, 'a');
    }
    SECTION("build fields are typed even if the document hash matches") {
        provenance["build"]["uncommitted_changes"] = 0;
    }
    SECTION("new build schemas are unsupported") {
        provenance["build"]["schema"] = "uvlas.build-info/2";
        expected = FailureCode::Unsupported;
    }
    const auto provenance_text = provenance.dump();
    if (consistent_digest)
        status["provenance_sha256"] = hex(source_digest(bytes(provenance_text)).bytes);
    auto status_text = status.dump();
    if (duplicate_status)
        status_text.insert(1, R"("compl\u0065ted":false,)");
    ResourcePlan plan{parser_budget};
    const auto result =
        read_result_documents(bytes(status_text), bytes(provenance_text), read_limits, plan);
    REQUIRE_FALSE(result);
    CHECK(result.error().code == expected);
    CHECK(result.error().operation == Operation::Recover);
}

TEST_CASE("Result inspection checks exact source admission, plan and cancellation",
          "[IO4-05][N4-02][negative]") {
    ResourcePlan writer{capacity};
    const auto completed = take(ResultStatus::completed());
    const auto docs = take(write_result_documents(completed, {}, capacity, writer, {}));
    const auto status = bytes(docs.status());
    const auto provenance = bytes(docs.provenance());
    ResourcePlan first{parser_budget};
    const auto baseline = take(read_result_documents(status, provenance, read_limits, first));
    CHECK(baseline.status().exit_code == 0);
    CHECK_FALSE(baseline.status().failure);
    CHECK_FALSE(baseline.status().last_accepted_id);
    CHECK_FALSE(baseline.provenance().original_case);
    auto caps = read_limits;
    caps.maximum_source_bytes = std::max(status.size(), provenance.size());
    ResourcePlan exact{first.allocated_bytes()};
    CHECK(read_result_documents(status, provenance, caps, exact));
    CHECK(exact.allocated_bytes() == first.allocated_bytes());
    ResourcePlan short_plan{exact.allocated_bytes() - 1};
    const auto over_plan = read_result_documents(status, provenance, caps, short_plan);
    REQUIRE_FALSE(over_plan);
    CHECK(over_plan.error().code == FailureCode::ResourceExceeded);
    --caps.maximum_source_bytes;
    ResourcePlan short_input{parser_budget};
    const auto over_input = read_result_documents(status, provenance, caps, short_input);
    REQUIRE_FALSE(over_input);
    CHECK(over_input.error().code == FailureCode::ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC request_stop() requires a mutable source.
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    ResourcePlan cancelled{0};
    const auto rejected =
        read_result_documents(status, provenance, caps, cancelled, stop.get_token());
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == FailureCode::Cancelled);
    CHECK(cancelled.allocated_bytes() == 0);
}
} // namespace
