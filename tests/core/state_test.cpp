#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <string_view>
#include <thread>
#include <type_traits>
#include <uvlas/core/state.hpp>

namespace {
using namespace uvlas::core;
template <class T>
concept RvaluePayload = requires(T&& value) { std::move(value).payload(); };
template <class T>
concept RvalueModel = requires(T&& value) { std::move(value).model(); };
template <class T>
concept RvalueDefinition = requires(T&& value) { std::move(value).definition(); };
template <class T>
concept RvalueRequirements = requires(T&& value) { std::move(value).requirements(); };
template <class T>
concept RvalueTarget = requires(T&& value) { std::move(value).target(); };
template <class T>
concept RvalueAssessments = requires(T&& value) { std::move(value).assessments(); };
constexpr std::array model_bytes{std::byte{'m'}, std::byte{'o'}, std::byte{'d'}, std::byte{'e'},
                                 std::byte{'l'}};
constexpr std::array initial_bytes{std::byte{1}, std::byte{2}, std::byte{3}};
constexpr std::array changed_bytes{std::byte{4}, std::byte{5}};

AcceptancePlan requirements() {
    return {{{ObservationId{10}, ScopeId{100}, 1.0, Precision::Binary64},
             {ObservationId{11}, ScopeId{101}, 1.0, Precision::Binary64},
             {ObservationId{12}, ScopeId{102}, 1.0, Precision::Binary64},
             {ObservationId{13}, ScopeId{103}, 1.0, Precision::Binary64}}};
}

Model model(ResourcePlan& plan, ModelGeneration generation = ModelGeneration{1}) {
    const auto result = Model::create(model_bytes, generation, requirements(), plan);
    REQUIRE(result);
    return *result;
}

ErrorCertificate assessment(const TrialState& trial, std::size_t condition,
                            std::optional<double> bound = 0.5) {
    const auto& requirement = trial.model().requirements()[condition];
    const CertificateContext context{requirement.observation, requirement.scope,
                                     trial.target().numerical_generation,
                                     requirement.representation};
    auto result = bound ? ErrorCertificate::assess(context, requirement.budget, *bound)
                        : ErrorCertificate::unevaluated(context, requirement.budget);
    REQUIRE(result);
    return *result;
}

std::array<ErrorCertificate, 4> assessments(const TrialState& trial) {
    return {assessment(trial, 0), assessment(trial, 1), assessment(trial, 2), assessment(trial, 3)};
}

QualificationPacket qualify_trial(const TrialState& trial) {
    const auto evaluated_target = trial.target();
    const auto report =
        EvaluationReport::create(trial.model(), evaluated_target, assessments(trial));
    REQUIRE(report);
    auto packet = qualify(trial, *report);
    REQUIRE(packet);
    return std::move(*packet);
}

AcceptedState initial(const Model& prepared, ResourcePlan& plan,
                      AcceptedStateId id = AcceptedStateId{42},
                      NumericalGeneration generation = NumericalGeneration{1},
                      std::span<const std::byte> bytes = initial_bytes) {
    auto trial = TrialState::initial(prepared, TrialId{1}, generation, bytes, plan);
    REQUIRE(trial);
    auto packet = qualify_trial(*trial);
    auto result = commit_initial(*trial, packet, id);
    REQUIRE(result);
    return std::move(*result);
}

template <class T>
void check_failure(const Result<T>& result, FailureCode code) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
}

std::array<std::byte, 32> digest(std::string_view hex) {
    std::array<std::byte, 32> result{};
    auto nibble = [](char character) {
        return character <= '9' ? character - '0' : character - 'a' + 10;
    };
    REQUIRE(hex.size() == 64);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = static_cast<std::byte>(16 * nibble(hex[2 * i]) + nibble(hex[2 * i + 1]));
    }
    return result;
}
} // namespace

TEST_CASE("Models own immutable definitions and content fingerprints", "[unit][N4-03]") {
    ResourcePlan plan{1024};
    auto definition = model_bytes;
    const auto prepared = Model::create(definition, ModelGeneration{1}, requirements(), plan);
    REQUIRE(prepared);
    definition[0] = std::byte{0};
    CHECK(prepared->definition()[0] == std::byte{'m'});
    // Independent Python hashlib of the documented explicit model encoding.
    CHECK(prepared->fingerprint().bytes ==
          digest("5a87a80088356be4b752c35502ff30be36417655231d85a791b12d6f612c6631"));
    const auto same = model(plan, ModelGeneration{2});
    CHECK(same.fingerprint() == prepared->fingerprint());
    auto altered = requirements();
    altered[0].budget = 0.5;
    const auto different = Model::create(model_bytes, ModelGeneration{1}, altered, plan);
    REQUIRE(different);
    CHECK(different->fingerprint() != prepared->fingerprint());
    altered[1].observation = altered[0].observation;
    check_failure(Model::create(model_bytes, ModelGeneration{1}, altered, plan),
                  FailureCode::InvalidInput);
    altered = requirements();
    altered[0].budget = std::numeric_limits<double>::quiet_NaN();
    check_failure(Model::create(model_bytes, ModelGeneration{1}, altered, plan),
                  FailureCode::Nonfinite);
}

TEST_CASE("Qualified commits transfer owned state once and preserve the parent", "[unit][N4-03]") {
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    const auto parent = initial(prepared, plan);
    const auto parent_hash =
        digest("039058c6f2c0cb492c533b0a4d14ef77cc0f78abccced5287d84a1a2011cfb81");
    auto candidate =
        TrialState::begin(parent, TrialId{2}, NumericalGeneration{2}, changed_bytes, plan);
    REQUIRE(candidate);
    if (const auto parent_stamp = candidate->target().parent; parent_stamp) {
        CHECK(parent_stamp->payload.bytes == parent_hash);
    } else {
        FAIL("Accepted parent identity was omitted");
    }
    auto packet = qualify_trial(*candidate);
    const auto original_budget = plan.allocated_bytes();
    const auto accepted = commit(parent, *candidate, packet);
    REQUIRE(accepted);
    CHECK(accepted->id() == AcceptedStateId{43});
    CHECK(accepted->payload()[0] == changed_bytes[0]);
    CHECK(parent.id() == AcceptedStateId{42});
    CHECK(std::ranges::equal(parent.payload(), initial_bytes));
    CHECK(plan.allocated_bytes() == original_budget);
    CHECK_FALSE(candidate->active());
    CHECK(packet.consumed());
    check_failure(commit(parent, *candidate, packet), FailureCode::StaleGeneration);
    check_failure(candidate->replace_payload(initial_bytes, plan), FailureCode::StaleGeneration);
}

TEST_CASE("Every unachieved category prevents qualification", "[negative][N4-03]") {
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    const auto candidate =
        TrialState::initial(prepared, TrialId{4}, NumericalGeneration{5}, initial_bytes, plan);
    REQUIRE(candidate);
    const std::array exceeded{FailureCode::NonlinearUnconverged, FailureCode::GeometryInvalid,
                              FailureCode::PrecisionUnresolved, FailureCode::PrecisionUnresolved};
    for (std::size_t condition = 0; condition < 4; ++condition) {
        for (const bool evaluated : {false, true}) {
            auto certificates = assessments(*candidate);
            certificates[condition] = assessment(
                *candidate, condition, evaluated ? std::optional<double>{2.0} : std::nullopt);
            const auto report =
                EvaluationReport::create(candidate->model(), candidate->target(), certificates);
            REQUIRE(report);
            check_failure(qualify(*candidate, *report),
                          evaluated ? exceeded[condition] : FailureCode::PrecisionUnresolved);
            CHECK(candidate->active());
        }
    }
}

TEST_CASE("Assessment metadata must match the model and numerical generation",
          "[negative][N4-03]") {
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    const auto candidate =
        TrialState::initial(prepared, TrialId{4}, NumericalGeneration{5}, initial_bytes, plan);
    REQUIRE(candidate);
    for (const int altered_field : {0, 1, 2, 3, 4}) {
        auto certificates = assessments(*candidate);
        CertificateContext context{ObservationId{10}, ScopeId{100}, NumericalGeneration{5},
                                   Precision::Binary64};
        double budget = 1.0;
        if (altered_field == 0)
            context.observation = ObservationId{99};
        if (altered_field == 1)
            context.scope = ScopeId{99};
        if (altered_field == 2)
            context.numerical_generation = NumericalGeneration{99};
        if (altered_field == 3)
            context.representation = Precision::TwoWord;
        if (altered_field == 4)
            budget = 2.0;
        const auto wrong = ErrorCertificate::assess(context, budget, 0.5);
        REQUIRE(wrong);
        certificates[0] = *wrong;
        check_failure(
            EvaluationReport::create(candidate->model(), candidate->target(), certificates),
            FailureCode::StaleGeneration);
    }
}

TEST_CASE("Reports and packets cannot authorize another trial or model generation",
          "[negative][N4-03]") {
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    const auto same_content_new_model = model(plan, ModelGeneration{2});
    auto first =
        TrialState::initial(prepared, TrialId{4}, NumericalGeneration{5}, initial_bytes, plan);
    auto different_trial =
        TrialState::initial(prepared, TrialId{6}, NumericalGeneration{5}, initial_bytes, plan);
    auto different_model = TrialState::initial(same_content_new_model, TrialId{4},
                                               NumericalGeneration{5}, initial_bytes, plan);
    REQUIRE(first);
    REQUIRE(different_trial);
    REQUIRE(different_model);
    const auto report =
        EvaluationReport::create(first->model(), first->target(), assessments(*first));
    REQUIRE(report);
    check_failure(qualify(*different_trial, *report), FailureCode::StaleGeneration);
    check_failure(qualify(*different_model, *report), FailureCode::StaleGeneration);
    auto packet = qualify_trial(*first);
    check_failure(commit_initial(*different_trial, packet, AcceptedStateId{0}),
                  FailureCode::StaleGeneration);
    CHECK_FALSE(packet.consumed());
}

TEST_CASE("Replacement invalidates packets even when old candidate bytes are restored",
          "[negative][N4-03]") {
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    const auto parent = initial(prepared, plan);
    auto candidate =
        TrialState::begin(parent, TrialId{2}, NumericalGeneration{2}, changed_bytes, plan);
    REQUIRE(candidate);
    const auto report =
        EvaluationReport::create(candidate->model(), candidate->target(), assessments(*candidate));
    REQUIRE(report);
    auto packet = qualify_trial(*candidate);
    REQUIRE(candidate->replace_payload(initial_bytes, plan));
    REQUIRE(candidate->replace_payload(changed_bytes, plan));
    check_failure(qualify(*candidate, *report), FailureCode::StaleGeneration);
    check_failure(commit(parent, *candidate, packet), FailureCode::StaleGeneration);
    CHECK_FALSE(packet.consumed());
    CHECK(std::ranges::equal(parent.payload(), initial_bytes));
    auto fresh = qualify_trial(*candidate);
    REQUIRE(commit(parent, *candidate, fresh));
}

TEST_CASE("Evaluation results retain the target captured before evaluation", "[negative][N4-03]") {
    ResourcePlan plan{4096};
    const auto prepared = model(plan);
    auto candidate =
        TrialState::initial(prepared, TrialId{2}, NumericalGeneration{5}, initial_bytes, plan);
    auto other =
        TrialState::initial(prepared, TrialId{3}, NumericalGeneration{5}, initial_bytes, plan);
    REQUIRE(candidate);
    REQUIRE(other);
    const auto started_target = candidate->target();
    const auto completed_assessments = assessments(*candidate);
    REQUIRE(candidate->replace_payload(changed_bytes, plan));
    // Constructing the report after mutation must preserve the original target,
    // not attach completed evidence to the candidate that happens to be current.
    const auto report = EvaluationReport::create(prepared, started_target, completed_assessments);
    REQUIRE(report);
    CHECK(report->target() == started_target);
    check_failure(qualify(*candidate, *report), FailureCode::StaleGeneration);
    check_failure(qualify(*other, *report), FailureCode::StaleGeneration);
    auto corrupted_target = other->target();
    corrupted_target.payload.bytes[0] ^= std::byte{1};
    const auto corrupted_report =
        EvaluationReport::create(prepared, corrupted_target, assessments(*other));
    REQUIRE(corrupted_report);
    check_failure(qualify(*other, *corrupted_report), FailureCode::StaleGeneration);
}

TEST_CASE("Wrong parent content generation and model are rejected before ownership changes",
          "[negative][N4-03]") {
    ResourcePlan plan{4096};
    const auto prepared = model(plan);
    const auto parent = initial(prepared, plan);
    auto candidate =
        TrialState::begin(parent, TrialId{2}, NumericalGeneration{2}, changed_bytes, plan);
    REQUIRE(candidate);
    auto packet = qualify_trial(*candidate);
    const auto wrong_id = initial(prepared, plan, AcceptedStateId{99});
    const auto wrong_content =
        initial(prepared, plan, AcceptedStateId{42}, NumericalGeneration{1}, changed_bytes);
    const auto wrong_generation =
        initial(prepared, plan, AcceptedStateId{42}, NumericalGeneration{3});
    const auto different_model = model(plan, ModelGeneration{2});
    const auto wrong_model = initial(different_model, plan);
    for (const AcceptedState* other :
         {&wrong_id, &wrong_content, &wrong_generation, &wrong_model}) {
        const auto rejected = commit(*other, *candidate, packet);
        check_failure(rejected, FailureCode::StaleGeneration);
        CHECK(rejected.error().details.last_accepted == other->id());
        CHECK(candidate->active());
        CHECK_FALSE(packet.consumed());
    }
    REQUIRE(commit(parent, *candidate, packet));
}

TEST_CASE("Cancellation resource failure and identity exhaustion preserve accepted state",
          "[negative][N4-03]") {
    ResourcePlan plan{4096};
    const auto prepared = model(plan);
    const auto parent = initial(prepared, plan);
    for (const auto old_generation : {NumericalGeneration{0}, NumericalGeneration{1}}) {
        check_failure(TrialState::begin(parent, TrialId{99}, old_generation, changed_bytes, plan),
                      FailureCode::StaleGeneration);
    }
    auto candidate =
        TrialState::begin(parent, TrialId{2}, NumericalGeneration{2}, changed_bytes, plan);
    REQUIRE(candidate);
    const auto target = candidate->target();
    ResourcePlan insufficient{0};
    check_failure(candidate->replace_payload(initial_bytes, insufficient),
                  FailureCode::ResourceExceeded);
    CHECK(candidate->target() == target);
    auto packet = qualify_trial(*candidate);
    CHECK(candidate->cancel().code == FailureCode::Cancelled);
    check_failure(commit(parent, *candidate, packet), FailureCode::StaleGeneration);
    CHECK(std::ranges::equal(parent.payload(), initial_bytes));
    const auto exhausted =
        initial(prepared, plan, AcceptedStateId{std::numeric_limits<std::uint64_t>::max()});
    auto next =
        TrialState::begin(exhausted, TrialId{3}, NumericalGeneration{2}, changed_bytes, plan);
    REQUIRE(next);
    auto next_packet = qualify_trial(*next);
    check_failure(commit(exhausted, *next, next_packet), FailureCode::ResourceExceeded);
    CHECK(next->active());
    CHECK_FALSE(next_packet.consumed());
    auto generation_exhausted = TrialState::initial(
        prepared, TrialId{7}, NumericalGeneration{std::numeric_limits<std::uint64_t>::max()},
        initial_bytes, plan);
    REQUIRE(generation_exhausted);
    const auto before = generation_exhausted->target();
    check_failure(generation_exhausted->replace_payload(changed_bytes, plan),
                  FailureCode::ResourceExceeded);
    CHECK(generation_exhausted->target() == before);
}

TEST_CASE("All cache generations are checked even when only the right hand side changes",
          "[unit][N4-03]") {
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    const CacheKey key{prepared.fingerprint(),  ModelGeneration{1},        GeometryGeneration{2},
                       ChartGeneration{3},      SparsityGeneration{4},     NumericalGeneration{5},
                       ScaleGeneration{6},      PrecisionGeneration{7},    Precision::Binary64,
                       WakePolicyGeneration{8}, TimeIntervalGeneration{9}, BackendGeneration{10}};
    CHECK(can_reuse(key, {key, RhsGeneration{1}}));
    CHECK(can_reuse(key, {key, RhsGeneration{2}}));
    const std::array<void (*)(CacheKey&), 12> change{
        [](CacheKey& value) { value.fingerprint.bytes[0] ^= std::byte{1}; },
        [](CacheKey& value) { value.model = ModelGeneration{99}; },
        [](CacheKey& value) { value.geometry = GeometryGeneration{99}; },
        [](CacheKey& value) { value.chart = ChartGeneration{99}; },
        [](CacheKey& value) { value.sparsity = SparsityGeneration{99}; },
        [](CacheKey& value) { value.matrix = NumericalGeneration{99}; },
        [](CacheKey& value) { value.scales = ScaleGeneration{99}; },
        [](CacheKey& value) { value.precision = PrecisionGeneration{99}; },
        [](CacheKey& value) { value.representation = Precision::TwoWord; },
        [](CacheKey& value) { value.wake = WakePolicyGeneration{99}; },
        [](CacheKey& value) { value.interval = TimeIntervalGeneration{99}; },
        [](CacheKey& value) { value.backend = BackendGeneration{99}; }};
    for (const auto modify : change) {
        auto other = key;
        modify(other);
        CHECK_FALSE(can_reuse(key, {other, RhsGeneration{1}}));
        CHECK_FALSE(can_reuse(key, {other, RhsGeneration{2}}));
    }
    auto unsupported = key;
    unsupported.representation = static_cast<Precision>(255);
    CHECK_FALSE(can_reuse(unsupported, {unsupported, RhsGeneration{1}}));
}

TEST_CASE("Workspace leases reject reentry and concurrent mutable sharing", "[negative][N4-03]") {
    ResourcePlan plan{1024};
    auto workspace = Workspace::create(4, plan);
    REQUIRE(workspace);
    {
        auto lease = (*workspace)->acquire();
        REQUIRE(lease);
        check_failure((*workspace)->acquire(), FailureCode::InvalidInput);
        std::optional<FailureCode> other_failure;
        std::thread other{[&] {
            const auto concurrent = (*workspace)->acquire();
            if (!concurrent)
                other_failure = concurrent.error().code;
        }};
        other.join();
        CHECK(other_failure == FailureCode::InvalidInput);
        lease->bytes()[0] = std::byte{7};
    }
    auto again = (*workspace)->acquire();
    REQUIRE(again);
    CHECK(again->bytes()[0] == std::byte{7});
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<Workspace>);
    STATIC_REQUIRE_FALSE(std::is_move_constructible_v<Workspace>);
}

TEST_CASE("Returned states never borrow reusable workspace storage", "[unit][N4-03]") {
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    auto workspace = Workspace::create(3, plan);
    REQUIRE(workspace);
    auto candidate = [&]() {
        auto lease = (*workspace)->acquire();
        REQUIRE(lease);
        std::ranges::copy(initial_bytes, lease->bytes().begin());
        return TrialState::initial(prepared, TrialId{1}, NumericalGeneration{1}, lease->bytes(),
                                   plan);
    }();
    REQUIRE(candidate);
    {
        auto lease = (*workspace)->acquire();
        REQUIRE(lease);
        std::ranges::fill(lease->bytes(), std::byte{99});
    }
    workspace->reset();
    CHECK(std::ranges::equal(candidate->payload(), initial_bytes));
    auto packet = qualify_trial(*candidate);
    const auto accepted = commit_initial(*candidate, packet, AcceptedStateId{0});
    REQUIRE(accepted);
    CHECK(std::ranges::equal(accepted->payload(), initial_bytes));
}

TEST_CASE("Exception unwinding releases exclusive workspace access", "[negative][N4-03]") {
    ResourcePlan plan{16};
    auto workspace = Workspace::create(4, plan);
    REQUIRE(workspace);
    const auto failed =
        guard_entry(Operation::EvaluateStep, FailureCode::Cancelled, {}, [&]() -> Result<void> {
            auto lease = (*workspace)->acquire();
            if (!lease)
                return std::unexpected(lease.error());
            throw 7;
        });
    check_failure(failed, FailureCode::Cancelled);
    REQUIRE((*workspace)->acquire());
}

TEST_CASE("Empty model and state payloads are valid owned values", "[unit][N4-03]") {
    ResourcePlan plan{0};
    const auto prepared = Model::create({}, ModelGeneration{0}, requirements(), plan);
    REQUIRE(prepared);
    const auto parent = initial(*prepared, plan, AcceptedStateId{0}, NumericalGeneration{0}, {});
    CHECK(parent.payload().empty());
    auto candidate = TrialState::begin(parent, TrialId{0}, NumericalGeneration{1}, {}, plan);
    REQUIRE(candidate);
    auto packet = qualify_trial(*candidate);
    const auto accepted = commit(parent, *candidate, packet);
    REQUIRE(accepted);
    CHECK(accepted->id() == AcceptedStateId{1});
    CHECK(accepted->payload().empty());
    CHECK(plan.allocated_bytes() == 0);
}

TEST_CASE("Borrowing APIs reject direct temporaries and ownership moves consume the source",
          "[unit][N4-03]") {
    STATIC_REQUIRE_FALSE(RvaluePayload<TrialState>);
    STATIC_REQUIRE_FALSE(RvaluePayload<AcceptedState>);
    STATIC_REQUIRE_FALSE(RvalueModel<TrialState>);
    STATIC_REQUIRE_FALSE(RvalueModel<AcceptedState>);
    STATIC_REQUIRE_FALSE(RvalueDefinition<Model>);
    STATIC_REQUIRE_FALSE(RvalueRequirements<Model>);
    STATIC_REQUIRE_FALSE(RvalueTarget<EvaluationReport>);
    STATIC_REQUIRE_FALSE(RvalueTarget<QualificationPacket>);
    STATIC_REQUIRE_FALSE(RvalueAssessments<EvaluationReport>);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<TrialState>);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<AcceptedState>);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<QualificationPacket>);
    ResourcePlan plan{1024};
    const auto prepared = model(plan);
    auto candidate =
        TrialState::initial(prepared, TrialId{1}, NumericalGeneration{1}, initial_bytes, plan);
    REQUIRE(candidate);
    TrialState moved{std::move(*candidate)};
    CHECK_FALSE(candidate->active());
    CHECK(candidate->payload().empty());
    auto packet = qualify_trial(moved);
    QualificationPacket transferred{std::move(packet)};
    // NOLINTNEXTLINE(bugprone-use-after-move): the packet move contract consumes the source.
    CHECK(packet.consumed());
    check_failure(commit_initial(moved, packet, AcceptedStateId{0}), FailureCode::StaleGeneration);
    auto accepted = commit_initial(moved, transferred, AcceptedStateId{0});
    REQUIRE(accepted);
    const AcceptedState final_state{std::move(*accepted)};
    CHECK_FALSE(accepted->active());
    CHECK(accepted->payload().empty());
    CHECK(std::ranges::equal(final_state.payload(), initial_bytes));
}
