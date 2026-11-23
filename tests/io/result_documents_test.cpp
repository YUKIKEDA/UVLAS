#include "result_documents.hpp"
#include "source_digest.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

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

template <class T>
T take(Result<T> value) {
    REQUIRE(value);
    if (!value) {
        throw std::logic_error("Expected a result value");
    }
    return std::move(*value);
}

AcceptedState accept(ResourcePlan& plan) {
    const AcceptancePlan requirements{{{ObservationId{1}, ScopeId{20}, 1, Precision::Binary64},
                                       {ObservationId{2}, ScopeId{20}, 1, Precision::Binary64},
                                       {ObservationId{3}, ScopeId{20}, 1, Precision::Binary64},
                                       {ObservationId{4}, ScopeId{20}, 1, Precision::Binary64}}};
    const auto model = take(Model::create({}, ModelGeneration{maximum}, requirements, plan));
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
    }
}

} // namespace
