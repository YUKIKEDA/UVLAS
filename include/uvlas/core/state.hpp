#pragma once

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <uvlas/core/error_certificate.hpp>
#include <uvlas/core/resource_plan.hpp>

namespace uvlas::core {

/// SHA-256 identity of a complete immutable model definition and its acceptance requirements.
struct ModelFingerprint {
    /// Digest bytes in SHA-256 output order; no host object representation is serialized.
    std::array<std::byte, 32> bytes;
    /// Compare all digest bytes exactly.
    bool operator==(const ModelFingerprint&) const = default;
};

/// SHA-256 identity of the canonical bytes owned by a state candidate.
struct StateFingerprint {
    /// Digest bytes in SHA-256 output order.
    std::array<std::byte, 32> bytes;
    /// Compare all digest bytes exactly.
    bool operator==(const StateFingerprint&) const = default;
};

/// Required categories evaluated before a candidate can be accepted.
enum class AcceptanceCondition : std::uint8_t {
    Residual, ///< Governing residual checks.
    Geometry, ///< Geometric admissibility checks.
    Time,     ///< Time-discretization checks.
    Precision ///< Retained and consumed precision checks.
};

/// One model-defined error requirement; quantities use the same units or normalization as their
/// bound.
struct ObservationRequirement {
    /// Unique observable for this acceptance category.
    ObservationId observation;
    /// Domain or subset over which the bound must hold.
    ScopeId scope;
    /// Nonnegative finite permitted error.
    double budget;
    /// Arithmetic representation actually required for evaluation.
    Precision representation;
};

/// Requirements in Residual, Geometry, Time, Precision order; all four are mandatory.
using AcceptancePlan = std::array<ObservationRequirement, 4>;

namespace detail {
struct ModelData;
}

/// Cheap owning reference to a model that cannot change after construction.
/// Definition bytes must be a complete, deterministic encoding supplied by model preparation;
/// this protocol does not define a physical model or state serialization format.
class Model {
public:
    /// Copy an owning reference to immutable storage; no mutable alias is exposed.
    Model(const Model&) noexcept = default;
    /// Replace the owning reference; borrowed views into the old model may expire.
    Model& operator=(const Model&) noexcept = default;

    /// Copy canonical definition bytes, validate four distinct observations and compute their
    /// digest.
    /// @return An immutable owning model, or a finite-value/input/representation/resource failure.
    /// Array payload uses plan; shared-reference metadata is ordinary allocation overhead.
    [[nodiscard]] static Result<Model> create(std::span<const std::byte> definition,
                                              ModelGeneration generation,
                                              AcceptancePlan requirements,
                                              ResourcePlan& plan) noexcept;
    /// Read the immutable definition while this owning reference remains alive and unchanged.
    [[nodiscard]] std::span<const std::byte> definition() const& noexcept;
    /// A temporary model cannot provide a borrowed definition.
    [[nodiscard]] std::span<const std::byte> definition() const&& = delete;
    /// Model generation supplied during preparation.
    [[nodiscard]] ModelGeneration generation() const noexcept;
    /// Content identity, including all acceptance requirements.
    [[nodiscard]] ModelFingerprint fingerprint() const noexcept;
    /// Borrow requirements while this reference remains alive and unchanged.
    [[nodiscard]] const AcceptancePlan& requirements() const& noexcept;
    /// A temporary model cannot provide borrowed requirements.
    [[nodiscard]] const AcceptancePlan& requirements() const&& = delete;

private:
    explicit Model(std::shared_ptr<const detail::ModelData> data) noexcept
        : data_(std::move(data)) {}
    std::shared_ptr<const detail::ModelData> data_;
};

class TrialState;
class QualificationPacket;

/// Owns a previously accepted canonical payload; processing a trial never mutates it.
class AcceptedState {
public:
    /// Transfer payload ownership; the source can no longer act as an accepted parent.
    AcceptedState(AcceptedState&& other) noexcept;
    /// Transfer payload ownership and release the old payload.
    AcceptedState& operator=(AcceptedState&& other) noexcept;
    /// Accepted payloads are not implicitly copied.
    AcceptedState(const AcceptedState&) = delete;
    /// Accepted payloads are not implicitly copied.
    AcceptedState& operator=(const AcceptedState&) = delete;
    /// Release this owned payload, independently of any workspace.
    ~AcceptedState() = default;
    /// Identity of this accepted state.
    [[nodiscard]] AcceptedStateId id() const noexcept {
        return id_;
    }
    /// Generation of the accepted numerical values.
    [[nodiscard]] NumericalGeneration numerical_generation() const noexcept {
        return generation_;
    }
    /// Borrow the immutable model reference while this state lives.
    [[nodiscard]] const Model& model() const& noexcept {
        return model_;
    }
    /// A temporary state cannot provide a borrowed model reference.
    [[nodiscard]] const Model& model() const&& = delete;
    /// Borrow immutable payload bytes until this state is moved from, assigned to or destroyed.
    [[nodiscard]] std::span<const std::byte> payload() const& noexcept {
        return payload_.values();
    }
    /// A temporary state cannot provide a borrowed payload.
    [[nodiscard]] std::span<const std::byte> payload() const&& = delete;
    /// Whether this state still owns an accepted parent, including a valid empty payload.
    [[nodiscard]] bool active() const noexcept {
        return active_;
    }

private:
    friend class TrialState;
    friend Result<AcceptedState> commit_initial(TrialState&, QualificationPacket&,
                                                AcceptedStateId) noexcept;
    friend Result<AcceptedState> commit(const AcceptedState&, TrialState&,
                                        QualificationPacket&) noexcept;
    AcceptedState(const Model& model, AcceptedStateId id, NumericalGeneration generation,
                  PlannedBuffer<std::byte> payload) noexcept;
    Model model_;
    AcceptedStateId id_;
    NumericalGeneration generation_;
    PlannedBuffer<std::byte> payload_;
    bool active_ = true;
};

/// Accepted-parent identity and content, kept together to prevent partial parent stamps.
struct ParentIdentity {
    /// Accepted state number.
    AcceptedStateId id;
    /// Numerical generation of that accepted state.
    NumericalGeneration generation;
    /// Digest of the actual accepted payload.
    StateFingerprint payload;
    /// Compare all parent metadata and content exactly.
    bool operator==(const ParentIdentity&) const = default;
};

/// Complete immutable target to which evaluation reports and qualifications are bound.
struct EvaluationTarget {
    /// Model content identity.
    ModelFingerprint model;
    /// Model definition generation, distinct even when model content coincides.
    ModelGeneration model_generation;
    /// Complete accepted parent stamp, absent only for initialization.
    std::optional<ParentIdentity> parent;
    /// Candidate identity.
    TrialId trial;
    /// Candidate numerical generation.
    NumericalGeneration numerical_generation;
    /// Revision changes on every successful payload replacement, including restoration of old
    /// bytes.
    std::uint64_t revision;
    /// Digest computed from the candidate's actual owned payload.
    StateFingerprint payload;
    /// Compare every identity, revision and digest exactly.
    bool operator==(const EvaluationTarget&) const = default;
};

/// Owns an unaccepted candidate, separate from its read-only accepted parent.
class TrialState {
public:
    /// Transfer candidate ownership; the source cannot be edited, qualified or committed afterward.
    TrialState(TrialState&& other) noexcept;
    /// Transfer candidate ownership and release the previous candidate.
    TrialState& operator=(TrialState&& other) noexcept;
    /// Candidate ownership cannot be copied.
    TrialState(const TrialState&) = delete;
    /// Candidate ownership cannot be copied.
    TrialState& operator=(const TrialState&) = delete;
    /// Release owned candidate storage.
    ~TrialState() = default;

    /// Copy an initial candidate with no parent; numerical initialization checks are still
    /// required.
    [[nodiscard]] static Result<TrialState> initial(const Model& model, TrialId id,
                                                    NumericalGeneration generation,
                                                    std::span<const std::byte> payload,
                                                    ResourcePlan& plan) noexcept;
    /// Copy a candidate without modifying parent; generation must be newer than parent's
    /// generation.
    /// @return Owned candidate, StaleGeneration for a consumed parent or old generation,
    /// or a resource failure. Input bytes need only live for this synchronous call.
    [[nodiscard]] static Result<TrialState> begin(const AcceptedState& parent, TrialId id,
                                                  NumericalGeneration generation,
                                                  std::span<const std::byte> payload,
                                                  ResourcePlan& plan) noexcept;
    /// Replace payload transactionally and increment numerical generation and revision.
    /// Any prior report/packet becomes stale, even if the bytes are restored later.
    /// On failure this candidate is unchanged; consumed candidates return StaleGeneration.
    [[nodiscard]] Result<void> replace_payload(std::span<const std::byte> payload,
                                               ResourcePlan& plan) noexcept;
    /// Invalidate this candidate and return Cancelled without changing its parent.
    /// An already consumed candidate instead returns StaleGeneration.
    [[nodiscard]] Failure cancel() noexcept;
    /// Borrow immutable candidate bytes until replacement, move, commit, assignment or destruction.
    [[nodiscard]] std::span<const std::byte> payload() const& noexcept {
        return payload_.values();
    }
    /// A temporary candidate cannot provide a borrowed payload.
    [[nodiscard]] std::span<const std::byte> payload() const&& = delete;
    /// Borrow the immutable owning model reference while this candidate lives.
    [[nodiscard]] const Model& model() const& noexcept {
        return model_;
    }
    /// A temporary candidate cannot provide a borrowed model reference.
    [[nodiscard]] const Model& model() const&& = delete;
    /// Construct the current target by hashing the actual payload without allocating.
    [[nodiscard]] EvaluationTarget target() const noexcept;
    /// Whether this candidate can still be edited or evaluated.
    [[nodiscard]] bool active() const noexcept {
        return active_;
    }

private:
    friend Result<AcceptedState> commit_initial(TrialState&, QualificationPacket&,
                                                AcceptedStateId) noexcept;
    friend Result<AcceptedState> commit(const AcceptedState&, TrialState&,
                                        QualificationPacket&) noexcept;
    TrialState(const Model& model, TrialId id, NumericalGeneration generation,
               std::optional<ParentIdentity> parent, PlannedBuffer<std::byte> payload) noexcept;
    Model model_;
    TrialId id_;
    NumericalGeneration generation_;
    std::optional<ParentIdentity> parent_;
    PlannedBuffer<std::byte> payload_;
    std::uint64_t revision_ = 0;
    bool active_ = true;
};

/// Validated metadata and four explicit assessments, bound to one actual candidate revision.
/// The physical evaluator must compute conservative bounds for these exact bytes and scopes;
/// this type does not evaluate governing equations or certify physics by itself.
class EvaluationReport {
public:
    /// Preserve the target captured before evaluation with assessments in
    /// Residual, Geometry, Time, Precision order. Carry that target with the results;
    /// never replace it with a newer candidate's target after evaluation completes.
    /// @return A report if identities, generations, representations and exact budget bits match
    /// the model. Unevaluated/exceeded assessments remain explicit and cannot qualify a trial.
    [[nodiscard]] static Result<EvaluationReport>
    create(const Model& model, EvaluationTarget evaluated_target,
           std::array<ErrorCertificate, 4> assessments) noexcept;
    /// The exact candidate revision and content evaluated by this report.
    [[nodiscard]] const EvaluationTarget& target() const& noexcept {
        return target_;
    }
    /// A temporary report cannot provide a borrowed target.
    [[nodiscard]] const EvaluationTarget& target() const&& = delete;
    /// Borrow all four assessments while the report lives and is unchanged.
    [[nodiscard]] const std::array<ErrorCertificate, 4>& assessments() const& noexcept {
        return assessments_;
    }
    /// A temporary report cannot provide borrowed assessments.
    [[nodiscard]] const std::array<ErrorCertificate, 4>& assessments() const&& = delete;

private:
    EvaluationReport(EvaluationTarget target, std::array<ErrorCertificate, 4> assessments) noexcept
        : target_(target), assessments_(assessments) {}
    EvaluationTarget target_;
    std::array<ErrorCertificate, 4> assessments_;
};

/// A noncopyable single-use qualification for one candidate and accepted parent.
class QualificationPacket {
public:
    /// Transfer qualification; the source becomes consumed.
    QualificationPacket(QualificationPacket&& other) noexcept;
    /// Transfer qualification and invalidate the previous target.
    QualificationPacket& operator=(QualificationPacket&& other) noexcept;
    /// Qualifications cannot be duplicated.
    QualificationPacket(const QualificationPacket&) = delete;
    /// Qualifications cannot be duplicated.
    QualificationPacket& operator=(const QualificationPacket&) = delete;
    /// Whether a successful commit has consumed this packet, or ownership has moved away.
    [[nodiscard]] bool consumed() const noexcept {
        return consumed_;
    }
    /// Borrow the sealed target while this packet lives and is unchanged.
    [[nodiscard]] const EvaluationTarget& target() const& noexcept {
        return target_;
    }
    /// A temporary packet cannot provide a borrowed target.
    [[nodiscard]] const EvaluationTarget& target() const&& = delete;

private:
    friend Result<QualificationPacket> qualify(const TrialState&, const EvaluationReport&) noexcept;
    friend Result<AcceptedState> commit_initial(TrialState&, QualificationPacket&,
                                                AcceptedStateId) noexcept;
    friend Result<AcceptedState> commit(const AcceptedState&, TrialState&,
                                        QualificationPacket&) noexcept;
    explicit QualificationPacket(EvaluationTarget target) noexcept : target_(target) {}
    EvaluationTarget target_;
    bool consumed_ = false;
};

/// Issue a packet only when the report still matches the actual candidate and all four bounds
/// are achieved. Old/different targets return StaleGeneration; unevaluated bounds return
/// PrecisionUnresolved; exceeded residual/geometry/time/precision checks return respectively
/// NonlinearUnconverged, GeometryInvalid, PrecisionUnresolved, PrecisionUnresolved.
[[nodiscard]] Result<QualificationPacket> qualify(const TrialState& trial,
                                                  const EvaluationReport& report) noexcept;

/// Accept a qualified initial candidate exactly once, without allocating payload storage.
/// Parent-bearing candidates or consumed/stale targets return StaleGeneration.
[[nodiscard]] Result<AcceptedState> commit_initial(TrialState& trial, QualificationPacket& packet,
                                                   AcceptedStateId initial_id) noexcept;

/// Recheck parent identity/content/model and actual candidate digest, then transfer owned
/// candidate storage to a new accepted state. Parent never changes. The new identity is
/// parent.id()+1; exhaustion returns ResourceExceeded. Success consumes trial and packet;
/// every rejection leaves parent, trial and packet unchanged.
[[nodiscard]] Result<AcceptedState> commit(const AcceptedState& parent, TrialState& trial,
                                           QualificationPacket& packet) noexcept;

namespace detail {
struct GeometryGenerationTag;
struct ChartGenerationTag;
struct SparsityGenerationTag;
struct ScaleGenerationTag;
struct PrecisionGenerationTag;
struct WakePolicyGenerationTag;
struct TimeIntervalGenerationTag;
struct BackendGenerationTag;
struct RhsGenerationTag;
} // namespace detail
/// Generation of evaluated geometry.
using GeometryGeneration = Identifier<detail::GeometryGenerationTag>;
/// Generation of the constraint coordinate chart.
using ChartGeneration = Identifier<detail::ChartGenerationTag>;
/// Generation of the compressed sparse pattern.
using SparsityGeneration = Identifier<detail::SparsityGenerationTag>;
/// Generation of component scaling.
using ScaleGeneration = Identifier<detail::ScaleGenerationTag>;
/// Generation of the precision policy.
using PrecisionGeneration = Identifier<detail::PrecisionGenerationTag>;
/// Generation of the wake treatment policy.
using WakePolicyGeneration = Identifier<detail::WakePolicyGenerationTag>;
/// Identity of the time interval to which cached operations apply.
using TimeIntervalGeneration = Identifier<detail::TimeIntervalGenerationTag>;
/// Generation of the numerical backend policy.
using BackendGeneration = Identifier<detail::BackendGenerationTag>;
/// Generation of a right-hand side, separate from matrix factorization identity.
using RhsGeneration = Identifier<detail::RhsGenerationTag>;

/// All generations that determine whether a prepared numerical operation can be reused.
struct CacheKey {
    /// Model content identity.
    ModelFingerprint fingerprint;
    /// Model definition generation.
    ModelGeneration model;
    /// Geometry generation.
    GeometryGeneration geometry;
    /// Constraint chart generation.
    ChartGeneration chart;
    /// Compressed sparse pattern generation.
    SparsityGeneration sparsity;
    /// Numerical matrix generation, checked even for a changed right-hand side.
    NumericalGeneration matrix;
    /// Component-scale generation.
    ScaleGeneration scales;
    /// Precision-policy generation.
    PrecisionGeneration precision;
    /// Actual arithmetic representation.
    Precision representation;
    /// Wake-policy generation.
    WakePolicyGeneration wake;
    /// Time-interval generation.
    TimeIntervalGeneration interval;
    /// Backend generation.
    BackendGeneration backend;
    /// Compare every field; neither addresses nor dimensions establish identity.
    bool operator==(const CacheKey&) const = default;
};

/// Request for a prepared factorization with an independently identified right-hand side.
struct CacheRequest {
    /// Required operator/cache identity.
    CacheKey key;
    /// Right-hand-side generation; changing it never bypasses matrix identity checks.
    RhsGeneration rhs;
};

/// Reuse only when every operator generation and representation matches. A different
/// right-hand side may reuse factors, but never makes a different matrix reusable.
[[nodiscard]] bool can_reuse(const CacheKey& cached, const CacheRequest& request) noexcept;

class WorkspaceLease;

/// Caller-owned scratch storage with one active synchronous user at a time.
/// Concurrent or reentrant acquisition returns InvalidInput instead of sharing mutable data.
class Workspace {
public:
    /// Transfer previously planned scratch storage into a new exclusive workspace.
    explicit Workspace(PlannedBuffer<std::byte> scratch) noexcept : scratch_(std::move(scratch)) {}
    /// Allocate scratch payload through plan; owner metadata is ordinary allocation overhead.
    [[nodiscard]] static Result<std::unique_ptr<Workspace>> create(Index bytes,
                                                                   ResourcePlan& plan) noexcept;
    /// Workspace identity cannot be copied.
    Workspace(const Workspace&) = delete;
    /// Workspace identity cannot be reassigned.
    Workspace& operator=(const Workspace&) = delete;
    /// Atomic lease state cannot be moved while borrowed.
    Workspace(Workspace&&) = delete;
    /// Atomic lease state cannot be moved while borrowed.
    Workspace& operator=(Workspace&&) = delete;
    /// Release scratch storage; the caller must let every lease expire first.
    ~Workspace() = default;
    /// Acquire exclusive scratch access without allocation; Workspace must outlive the lease.
    [[nodiscard]] Result<WorkspaceLease> acquire() noexcept;

private:
    friend class WorkspaceLease;
    PlannedBuffer<std::byte> scratch_;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};

/// Exclusive borrowed workspace access; releasing or moving from it invalidates its views.
/// Results that must survive release must copy into their own storage before returning.
class WorkspaceLease {
public:
    /// Transfer exclusive access; the source becomes empty.
    WorkspaceLease(WorkspaceLease&& other) noexcept;
    /// Release old access and transfer the new lease.
    WorkspaceLease& operator=(WorkspaceLease&& other) noexcept;
    /// A lease cannot duplicate mutable access.
    WorkspaceLease(const WorkspaceLease&) = delete;
    /// A lease cannot duplicate mutable access.
    WorkspaceLease& operator=(const WorkspaceLease&) = delete;
    /// Release the workspace for its next user.
    ~WorkspaceLease();
    /// Borrow writable scratch until this lease is moved, assigned to or destroyed.
    /// A moved-from lease returns an empty view.
    [[nodiscard]] std::span<std::byte> bytes() & noexcept;
    /// Temporary leases cannot provide borrowed views.
    [[nodiscard]] std::span<std::byte> bytes() && = delete;

private:
    friend class Workspace;
    explicit WorkspaceLease(Workspace& owner) noexcept : owner_(&owner) {}
    Workspace* owner_;
};

} // namespace uvlas::core
