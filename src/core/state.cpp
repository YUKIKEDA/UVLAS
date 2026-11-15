#include <boost/hash2/sha2.hpp>

#include <algorithm>
#include <bit>
#include <limits>
#include <utility>
#include <uvlas/core/state.hpp>

namespace uvlas::core {
namespace detail {
struct ModelData {
    ModelData(PlannedBuffer<std::byte> definition_bytes, ModelGeneration model_generation,
              AcceptancePlan plan, ModelFingerprint content) noexcept
        : definition(std::move(definition_bytes)),
          generation(model_generation),
          requirements(plan),
          fingerprint(content) {}
    PlannedBuffer<std::byte> definition;
    ModelGeneration generation;
    AcceptancePlan requirements;
    ModelFingerprint fingerprint;
};
} // namespace detail
namespace {
using Digest = std::array<std::byte, 32>;

Digest finish(boost::hash2::sha2_256& hash) noexcept {
    const auto result = hash.result();
    Digest digest{};
    for (std::size_t i = 0; i < digest.size(); ++i) {
        digest[i] = static_cast<std::byte>(result[i]);
    }
    return digest;
}

StateFingerprint payload_hash(std::span<const std::byte> payload) noexcept {
    boost::hash2::sha2_256 hash;
    // Empty spans may carry null data; do not pass them through an implementation
    // that might form pointer arithmetic on a null input (N4-03).
    if (!payload.empty()) {
        hash.update(payload.data(), payload.size());
    }
    return {finish(hash)};
}

void append_u64(boost::hash2::sha2_256& hash, std::uint64_t value) noexcept {
    std::array<unsigned char, 8> encoded{};
    for (std::size_t i = 0; i < encoded.size(); ++i) {
        encoded[i] = static_cast<unsigned char>(value & 0xffU);
        value >>= 8U;
    }
    hash.update(encoded.data(), encoded.size());
}

ModelFingerprint model_hash(std::span<const std::byte> definition,
                            const AcceptancePlan& requirements) noexcept {
    boost::hash2::sha2_256 hash;
    constexpr std::array<unsigned char, 8> domain{'U', 'V', 'L', 'A', 'S', 'M', '0', '1'};
    hash.update(domain.data(), domain.size());
    append_u64(hash, definition.size());
    if (!definition.empty()) {
        hash.update(definition.data(), definition.size());
    }
    for (const auto& requirement : requirements) {
        append_u64(hash, requirement.observation.value());
        append_u64(hash, requirement.scope.value());
        append_u64(hash, std::bit_cast<std::uint64_t>(requirement.budget));
        const auto representation = static_cast<unsigned char>(requirement.representation);
        hash.update(&representation, 1);
    }
    return {finish(hash)};
}

Result<PlannedBuffer<std::byte>> copy_payload(std::span<const std::byte> bytes, ResourcePlan& plan,
                                              Operation operation) noexcept {
    if (!std::in_range<Index>(bytes.size())) {
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation});
    }
    auto storage = plan.allocate<std::byte>(static_cast<Index>(bytes.size()), operation);
    if (!storage) {
        return std::unexpected(storage.error());
    }
    std::ranges::copy(bytes, storage->values().begin());
    return storage;
}

Failure target_failure(FailureCode code, Operation operation,
                       const EvaluationTarget& target) noexcept {
    FailureDetails details;
    details.model = target.model_generation;
    if (target.parent) {
        details.parent = target.parent->id;
        details.last_accepted = target.parent->id;
    }
    details.trial = target.trial;
    details.numerical_generation = target.numerical_generation;
    return {code, operation, details};
}

Failure trial_failure(FailureCode code, Operation operation, const TrialState& trial) noexcept {
    return target_failure(code, operation, trial.target());
}

bool packet_matches(const TrialState& trial, const QualificationPacket& packet) noexcept {
    return trial.active() && !packet.consumed() && trial.target() == packet.target();
}
} // namespace

Result<Model> Model::create(std::span<const std::byte> definition, ModelGeneration generation,
                            AcceptancePlan requirements, ResourcePlan& plan) noexcept {
    for (std::size_t i = 0; i < requirements.size(); ++i) {
        const auto& requirement = requirements[i];
        const auto valid =
            ErrorCertificate::unevaluated({requirement.observation, requirement.scope,
                                           NumericalGeneration{0}, requirement.representation},
                                          requirement.budget);
        if (!valid) {
            auto failure = valid.error();
            failure.operation = Operation::PrepareModel;
            failure.details.model = generation;
            failure.details.numerical_generation.reset();
            return std::unexpected(failure);
        }
        for (std::size_t j = 0; j < i; ++j) {
            if (requirements[j].observation == requirement.observation) {
                FailureDetails details;
                details.model = generation;
                details.problem = requirement.observation;
                return std::unexpected(
                    Failure{FailureCode::InvalidInput, Operation::PrepareModel, details});
            }
        }
    }
    auto bytes = copy_payload(definition, plan, Operation::PrepareModel);
    if (!bytes) {
        auto failure = bytes.error();
        failure.details.model = generation;
        return std::unexpected(failure);
    }
    const auto fingerprint = model_hash(bytes->values(), requirements);
    FailureDetails details;
    details.model = generation;
    return guard_entry(Operation::PrepareModel, FailureCode::InvalidInput, details,
                       [&]() -> Result<Model> {
                           return Model{std::make_shared<const detail::ModelData>(
                               std::move(*bytes), generation, requirements, fingerprint)};
                       });
}

std::span<const std::byte> Model::definition() const& noexcept {
    return data_->definition.values();
}
ModelGeneration Model::generation() const noexcept {
    return data_->generation;
}
ModelFingerprint Model::fingerprint() const noexcept {
    return data_->fingerprint;
}
const AcceptancePlan& Model::requirements() const& noexcept {
    return data_->requirements;
}

AcceptedState::AcceptedState(const Model& model, AcceptedStateId id, NumericalGeneration generation,
                             PlannedBuffer<std::byte> payload) noexcept
    : model_(model), id_(id), generation_(generation), payload_(std::move(payload)) {}

AcceptedState::AcceptedState(AcceptedState&& other) noexcept
    : model_(other.model_),
      id_(other.id_),
      generation_(other.generation_),
      payload_(std::move(other.payload_)),
      active_(std::exchange(other.active_, false)) {}

AcceptedState& AcceptedState::operator=(AcceptedState&& other) noexcept {
    if (this != &other) {
        model_ = other.model_;
        id_ = other.id_;
        generation_ = other.generation_;
        payload_ = std::move(other.payload_);
        active_ = std::exchange(other.active_, false);
    }
    return *this;
}

TrialState::TrialState(const Model& model, TrialId id, NumericalGeneration generation,
                       std::optional<ParentIdentity> parent,
                       PlannedBuffer<std::byte> payload) noexcept
    : model_(model),
      id_(id),
      generation_(generation),
      parent_(parent),
      payload_(std::move(payload)) {}

TrialState::TrialState(TrialState&& other) noexcept
    : model_(other.model_),
      id_(other.id_),
      generation_(other.generation_),
      parent_(other.parent_),
      payload_(std::move(other.payload_)),
      revision_(other.revision_),
      active_(std::exchange(other.active_, false)) {}

TrialState& TrialState::operator=(TrialState&& other) noexcept {
    if (this != &other) {
        model_ = other.model_;
        id_ = other.id_;
        generation_ = other.generation_;
        parent_ = other.parent_;
        payload_ = std::move(other.payload_);
        revision_ = other.revision_;
        active_ = std::exchange(other.active_, false);
    }
    return *this;
}

Result<TrialState> TrialState::initial(const Model& model, TrialId id,
                                       NumericalGeneration generation,
                                       std::span<const std::byte> payload,
                                       ResourcePlan& plan) noexcept {
    auto bytes = copy_payload(payload, plan, Operation::Initialize);
    if (!bytes) {
        auto failure = bytes.error();
        failure.details.model = model.generation();
        failure.details.trial = id;
        failure.details.numerical_generation = generation;
        return std::unexpected(failure);
    }
    return TrialState{model, id, generation, std::nullopt, std::move(*bytes)};
}

Result<TrialState> TrialState::begin(const AcceptedState& parent, TrialId id,
                                     NumericalGeneration generation,
                                     std::span<const std::byte> payload,
                                     ResourcePlan& plan) noexcept {
    if (!parent.active_ || generation.value() <= parent.generation_.value()) {
        FailureDetails details;
        details.model = parent.model_.generation();
        details.parent = parent.id_;
        details.last_accepted = parent.id_;
        details.trial = id;
        details.numerical_generation = generation;
        return std::unexpected(
            Failure{FailureCode::StaleGeneration, Operation::AdvanceTrial, details});
    }
    auto bytes = copy_payload(payload, plan, Operation::AdvanceTrial);
    if (!bytes) {
        auto failure = bytes.error();
        failure.details.model = parent.model_.generation();
        failure.details.parent = parent.id_;
        failure.details.last_accepted = parent.id_;
        failure.details.trial = id;
        failure.details.numerical_generation = generation;
        return std::unexpected(failure);
    }
    return TrialState{
        parent.model_, id, generation,
        ParentIdentity{parent.id_, parent.generation_, payload_hash(parent.payload())},
        std::move(*bytes)};
}

Result<void> TrialState::replace_payload(std::span<const std::byte> payload,
                                         ResourcePlan& plan) noexcept {
    if (!active_) {
        return std::unexpected(
            trial_failure(FailureCode::StaleGeneration, Operation::AdvanceTrial, *this));
    }
    // revision starts at zero and grows with generation, so revision <= generation.
    // Checking generation also prevents revision overflow without an unreachable branch.
    if (generation_.value() == std::numeric_limits<std::uint64_t>::max()) {
        return std::unexpected(
            trial_failure(FailureCode::ResourceExceeded, Operation::AdvanceTrial, *this));
    }
    auto replacement = copy_payload(payload, plan, Operation::AdvanceTrial);
    if (!replacement) {
        auto failure = trial_failure(replacement.error().code, Operation::AdvanceTrial, *this);
        failure.details.resources = replacement.error().details.resources;
        return std::unexpected(failure);
    }
    payload_ = std::move(*replacement);
    generation_ = NumericalGeneration{generation_.value() + 1};
    ++revision_;
    return {};
}

Failure TrialState::cancel() noexcept {
    const auto failure =
        trial_failure(active_ ? FailureCode::Cancelled : FailureCode::StaleGeneration,
                      Operation::AdvanceTrial, *this);
    active_ = false;
    return failure;
}

EvaluationTarget TrialState::target() const noexcept {
    return {model_.fingerprint(),   model_.generation(), parent_, id_, generation_, revision_,
            payload_hash(payload())};
}

Result<EvaluationReport>
EvaluationReport::create(const Model& model, EvaluationTarget evaluated_target,
                         std::array<ErrorCertificate, 4> assessments) noexcept {
    if (evaluated_target.model != model.fingerprint() ||
        evaluated_target.model_generation != model.generation()) {
        return std::unexpected(target_failure(FailureCode::StaleGeneration, Operation::AdvanceTrial,
                                              evaluated_target));
    }
    const auto& requirements = model.requirements();
    for (std::size_t i = 0; i < assessments.size(); ++i) {
        const auto& expected = requirements[i];
        const auto& actual = assessments[i];
        if (actual.context().observation != expected.observation ||
            actual.context().scope != expected.scope ||
            actual.context().numerical_generation != evaluated_target.numerical_generation ||
            actual.context().representation != expected.representation ||
            std::bit_cast<std::uint64_t>(actual.budget()) !=
                std::bit_cast<std::uint64_t>(expected.budget)) {
            auto failure = target_failure(FailureCode::StaleGeneration, Operation::AdvanceTrial,
                                          evaluated_target);
            failure.details.problem = actual.context().observation;
            return std::unexpected(failure);
        }
    }
    return EvaluationReport{evaluated_target, assessments};
}

QualificationPacket::QualificationPacket(QualificationPacket&& other) noexcept
    : target_(other.target_), consumed_(std::exchange(other.consumed_, true)) {}

QualificationPacket& QualificationPacket::operator=(QualificationPacket&& other) noexcept {
    if (this != &other) {
        target_ = other.target_;
        consumed_ = std::exchange(other.consumed_, true);
    }
    return *this;
}

Result<QualificationPacket> qualify(const TrialState& trial,
                                    const EvaluationReport& report) noexcept {
    if (!trial.active() || trial.target() != report.target()) {
        return std::unexpected(
            trial_failure(FailureCode::StaleGeneration, Operation::AdvanceTrial, trial));
    }
    constexpr std::array failures{FailureCode::NonlinearUnconverged, FailureCode::GeometryInvalid,
                                  FailureCode::PrecisionUnresolved,
                                  FailureCode::PrecisionUnresolved};
    for (std::size_t i = 0; i < report.assessments().size(); ++i) {
        const auto& assessment = report.assessments()[i].assessment();
        if (!std::holds_alternative<Achieved>(assessment)) {
            auto failure = trial_failure(std::holds_alternative<Exceeded>(assessment)
                                             ? failures[i]
                                             : FailureCode::PrecisionUnresolved,
                                         Operation::AdvanceTrial, trial);
            failure.details.problem = report.assessments()[i].context().observation;
            return std::unexpected(failure);
        }
    }
    return QualificationPacket{trial.target()};
}

Result<AcceptedState> commit_initial(TrialState& trial, QualificationPacket& packet,
                                     AcceptedStateId initial_id) noexcept {
    if (!packet_matches(trial, packet) || trial.parent_) {
        return std::unexpected(
            trial_failure(FailureCode::StaleGeneration, Operation::Commit, trial));
    }
    AcceptedState result{trial.model_, initial_id, trial.generation_, std::move(trial.payload_)};
    trial.active_ = false;
    packet.consumed_ = true;
    return result;
}

Result<AcceptedState> commit(const AcceptedState& parent, TrialState& trial,
                             QualificationPacket& packet) noexcept {
    if (!parent.active_ || !packet_matches(trial, packet) ||
        trial.parent_ !=
            ParentIdentity{parent.id_, parent.generation_, payload_hash(parent.payload())} ||
        trial.model_.fingerprint() != parent.model_.fingerprint() ||
        trial.model_.generation() != parent.model_.generation()) {
        auto failure = trial_failure(FailureCode::StaleGeneration, Operation::Commit, trial);
        failure.details.last_accepted = parent.id_;
        return std::unexpected(failure);
    }
    if (parent.id_.value() == std::numeric_limits<std::uint64_t>::max()) {
        return std::unexpected(
            trial_failure(FailureCode::ResourceExceeded, Operation::Commit, trial));
    }
    // All validation precedes ownership transfer; the parent is never modified (N4-03).
    AcceptedState result{trial.model_, AcceptedStateId{parent.id_.value() + 1}, trial.generation_,
                         std::move(trial.payload_)};
    trial.active_ = false;
    packet.consumed_ = true;
    return result;
}

bool can_reuse(const CacheKey& cached, const CacheRequest& request) noexcept {
    return (cached.representation == Precision::Binary64 ||
            cached.representation == Precision::TwoWord) &&
           cached == request.key;
}

Result<std::unique_ptr<Workspace>> Workspace::create(Index bytes, ResourcePlan& plan) noexcept {
    auto scratch = plan.allocate<std::byte>(bytes, Operation::PrepareModel);
    if (!scratch) {
        return std::unexpected(scratch.error());
    }
    return guard_entry(Operation::PrepareModel, FailureCode::InvalidInput, {},
                       [&]() -> Result<std::unique_ptr<Workspace>> {
                           return std::make_unique<Workspace>(std::move(*scratch));
                       });
}

Result<WorkspaceLease> Workspace::acquire() noexcept {
    if (busy_.test_and_set(std::memory_order_acquire)) {
        return std::unexpected(Failure{FailureCode::InvalidInput, Operation::EvaluateStep});
    }
    return WorkspaceLease{*this};
}

WorkspaceLease::WorkspaceLease(WorkspaceLease&& other) noexcept
    : owner_(std::exchange(other.owner_, nullptr)) {}

WorkspaceLease& WorkspaceLease::operator=(WorkspaceLease&& other) noexcept {
    if (this != &other) {
        if (owner_) {
            owner_->busy_.clear(std::memory_order_release);
        }
        owner_ = std::exchange(other.owner_, nullptr);
    }
    return *this;
}

WorkspaceLease::~WorkspaceLease() {
    if (owner_) {
        owner_->busy_.clear(std::memory_order_release);
    }
}

std::span<std::byte> WorkspaceLease::bytes() & noexcept {
    return owner_ ? owner_->scratch_.values() : std::span<std::byte>{};
}

} // namespace uvlas::core
