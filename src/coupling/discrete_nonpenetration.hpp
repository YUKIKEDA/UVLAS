#pragma once

#include "aero/p1_shared_action.hpp"
#include "case_known_flow.hpp"

#include <uvlas/numerics/linear_solve.hpp>

namespace uvlas::coupling::detail {
using NonpenetrationVector = std::array<aero::detail::ScalarWords, 3>;
using ExplicitCollocationTangents = std::array<NonpenetrationVector, 2>;

struct CaseWingNormalProjection {
    std::array<io::SiEnclosure, 3> unit;
    std::uint64_t used_steps{};
};

/// Enclose the oriented Case chart normal and qualify it against both retained child facets.
[[nodiscard]] core::Result<CaseWingNormalProjection> project_case_wing_normal(
    const aero::detail::OwnedGeometry&, std::span<const aero::detail::FacetId, 2>,
    const std::array<io::SiEnclosure, 3>& chord, const std::array<io::SiEnclosure, 3>& span,
    bool reversed, core::Precision, const std::array<CaseFlowComponentBudget, 3>&,
    CaseFlowProjectionLimits, std::stop_token) noexcept;

/// Retained current coordinates of an oriented facet. The tangents are the
/// exact differences 1-0 and 2-0, evaluated without intermediate projection.
struct RoundedFacetTangents {
    std::array<NonpenetrationVector, 3> points;
};

/// Explicit discrete collocation in the current geometry. The tangents orient its
/// normal and the velocity is in m/s. No physical PG mapping is inferred.
struct DiscreteCollocation {
    DiscreteCollocation(aero::detail::FacetId source, const aero::detail::ResolvedQuery& point,
                        const ExplicitCollocationTangents& frame,
                        const NonpenetrationVector& motion, const aero::detail::ScalarWords& gamma,
                        const aero::detail::ScalarWords& direction) noexcept
        : facet(source),
          query(point),
          tangents(frame),
          velocity(motion),
          circulation(gamma),
          circulation_direction(direction) {}
    DiscreteCollocation(const RoundedFacetTangents& frame, aero::detail::FacetId source,
                        const aero::detail::ResolvedQuery& point,
                        const NonpenetrationVector& motion, const aero::detail::ScalarWords& gamma,
                        const aero::detail::ScalarWords& direction) noexcept
        : facet(source),
          query(point),
          tangents(frame),
          velocity(motion),
          circulation(gamma),
          circulation_direction(direction) {}
    aero::detail::FacetId facet; ///< Either child of the bound physical quadrilateral.
    const aero::detail::ResolvedQuery& query;
    std::variant<ExplicitCollocationTangents, RoundedFacetTangents> tangents;
    NonpenetrationVector velocity;
    aero::detail::ScalarWords circulation;
    aero::detail::ScalarWords circulation_direction;
};

struct AttachedWakeFront {};
/// Fixed steady strip. The downstream trace uses the signed attached front column.
struct SteadyWakeClosure {};
struct DiscreteWakePanel {
    DiscreteWakePanel(aero::detail::FacetId source,
                      std::variant<AttachedWakeFront, aero::detail::P1StrengthId> left,
                      std::variant<aero::detail::P1StrengthId, SteadyWakeClosure> right) noexcept
        : facet(source), leading(left), trailing(right) {}
    aero::detail::FacetId facet; ///< Either child of the retained wake quadrilateral.
    std::variant<AttachedWakeFront, aero::detail::P1StrengthId> leading;
    std::variant<aero::detail::P1StrengthId, SteadyWakeClosure> trailing;
};

struct FrozenWakeTrace {
    aero::detail::P1StrengthId id;
    aero::detail::ScalarWords strength{aero::detail::SingleWord{0}};
};
struct FrozenWakeTraces {
    FrozenWakeTraces(core::AcceptedStateId accepted,
                     std::span<const FrozenWakeTrace> trace_columns) noexcept
        : parent(accepted), columns(trace_columns) {}
    core::AcceptedStateId parent;
    std::span<const FrozenWakeTrace> columns; ///< Strictly increasing independent history IDs.
};

struct DiscreteNonpenetrationInput {
    DiscreteNonpenetrationInput(const aero::detail::GlobalContactState& qualified,
                                std::span<const DiscreteCollocation> rows,
                                std::span<const DiscreteWakePanel> wake_panels,
                                FrozenWakeTraces frozen, aero::detail::P1ActionStamp trial,
                                aero::detail::pchip::PchipInputStamp known_flow,
                                aero::detail::pchip::RequestedSide side) noexcept
        : contact(qualified),
          collocations(rows),
          wake(wake_panels),
          history(frozen),
          trial_stamp(trial),
          flow_stamp(known_flow),
          flow_side(side) {}
    const aero::detail::GlobalContactState& contact;
    std::span<const DiscreteCollocation> collocations;
    std::span<const DiscreteWakePanel> wake;
    FrozenWakeTraces history;
    aero::detail::P1ActionStamp trial_stamp;
    aero::detail::pchip::PchipInputStamp flow_stamp;
    aero::detail::pchip::RequestedSide flow_side;
};

struct DiscreteNonpenetrationLimits {
    DiscreteNonpenetrationLimits(aero::detail::LocalAccuracy packet,
                                 aero::detail::P1Accuracy action,
                                 const aero::detail::P1AssemblyLimits& assembly,
                                 std::uint64_t per_action, std::uint64_t total) noexcept
        : geometry_accuracy(packet),
          induction_accuracy(action),
          induction(assembly),
          action_steps(per_action),
          steps(total) {}
    aero::detail::LocalAccuracy geometry_accuracy;
    aero::detail::P1Accuracy induction_accuracy;
    aero::detail::P1AssemblyLimits induction;
    std::uint64_t action_steps;
    std::uint64_t steps;
};

struct DiscreteNonpenetrationRow {
    aero::detail::PhysicalFaceId face{0};
    aero::detail::ScaledEnclosure residual;           ///< m/s, including the original affine terms.
    aero::detail::ScaledEnclosure circulation_action; ///< A_eff δΓ in m/s.
    bool
        includes_wake_quadrature{}; ///< The finite wake formula carries the existing P1 diagnostic.
};
struct NonpenetrationStorage;
struct NonpenetrationBuilder;

/// Own the finite discrete row values, original flow source and checked geometry.
class DiscreteNonpenetration {
public:
    DiscreteNonpenetration(DiscreteNonpenetration&&) noexcept;
    DiscreteNonpenetration& operator=(DiscreteNonpenetration&&) noexcept;
    DiscreteNonpenetration(const DiscreteNonpenetration&) = delete;
    DiscreteNonpenetration& operator=(const DiscreteNonpenetration&) = delete;
    ~DiscreteNonpenetration();
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] core::Result<DiscreteNonpenetrationRow> row(std::size_t) const noexcept;
    [[nodiscard]] const ProjectedCaseUniformFlow& flow() const& noexcept;
    [[nodiscard]] const ProjectedCaseUniformFlow& flow() const&& = delete;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;
    [[nodiscard]] aero::detail::GeometryContext geometry_context() const noexcept;
    [[nodiscard]] aero::detail::P1ActionStamp trial_stamp() const noexcept;

private:
    friend struct NonpenetrationBuilder;
    DiscreteNonpenetration(core::PlannedBuffer<NonpenetrationStorage>, ProjectedCaseUniformFlow&&,
                           std::optional<aero::detail::GlobalContactState>,
                           aero::detail::GeometryContext, aero::detail::P1ActionStamp,
                           std::uint64_t) noexcept;
    core::PlannedBuffer<NonpenetrationStorage> rows_;
    ProjectedCaseUniformFlow flow_;
    std::optional<aero::detail::GlobalContactState> geometry_owner_;
    aero::detail::GeometryContext geometry_context_;
    aero::detail::P1ActionStamp trial_stamp_;
    std::uint64_t used_steps_;
};

/// Evaluate the original residual and its affine circulation action with history
/// frozen. Resolve all attached L columns from the actual directed TE rows.
/// Consume the supplied flow only after every row succeeds.
[[nodiscard]] core::Result<DiscreteNonpenetration>
evaluate_discrete_nonpenetration(const DiscreteNonpenetrationInput&, ProjectedCaseUniformFlow&&,
                                 const DiscreteNonpenetrationLimits&,
                                 aero::detail::SharedP1Workspace&,
                                 aero::detail::LocalPacketWorkspace&, aero::detail::P1Workspace&,
                                 core::ResourcePlan&, const std::stop_token& = {}) noexcept;
/// Explicit capacities and budgets for the real dense reference route. Matrix and
/// RHS budgets bound their binary64 input projections in 1/m and m/s respectively.
struct DiscreteCirculationRequest {
    DiscreteCirculationRequest(numerics::linear::Scales physical_scales,
                               const numerics::linear::Limits& dense_limits,
                               numerics::linear::Accuracy algebraic_accuracy,
                               core::RhsGeneration rhs, core::NumericalGeneration retained,
                               double matrix_error, double rhs_error, double absolute,
                               double relative) noexcept
        : scales(physical_scales),
          solver_limits(dense_limits),
          algebraic(algebraic_accuracy),
          rhs_generation(rhs),
          retained_generation(retained),
          matrix_absolute(matrix_error),
          rhs_absolute(rhs_error),
          original_absolute(absolute),
          original_relative(relative) {}
    numerics::linear::Scales scales; ///< Ordered by ascending bound physical face ID.
    numerics::linear::Limits solver_limits;
    numerics::linear::Accuracy algebraic;
    core::RhsGeneration rhs_generation;
    core::NumericalGeneration retained_generation;
    double matrix_absolute;
    double rhs_absolute;
    double original_absolute; ///< Multiplier of the supplied velocity scales.
    double original_relative; ///< Multiplier of the unprojected affine RHS magnitude.
    /// Use outward Binary64 contour assembly for a primal fixed steady closure.
    bool certified_closed_contour{};
    /// Admit one to four independent workspaces for the steady contour solve.
    /// Each parallel range receives its own share of the remaining step budget.
    std::size_t contour_workers{1};
};

struct CirculationResidualRow {
    aero::detail::PhysicalFaceId face{0};
    aero::detail::ScaledEnclosure original_residual; ///< m/s from the native finite expression.
    bool includes_wake_quadrature{};
};

struct NativeCirculationRows;

/// Own both the LU-qualified retained words and the independent native residual
/// check. No physical mapping or continuous wake-integral accuracy is inferred.
class DiscreteCirculation {
public:
    DiscreteCirculation(DiscreteCirculation&&) noexcept;
    DiscreteCirculation& operator=(DiscreteCirculation&&) noexcept;
    ~DiscreteCirculation();
    DiscreteCirculation(const DiscreteCirculation&) = delete;
    DiscreteCirculation& operator=(const DiscreteCirculation&) = delete;
    [[nodiscard]] const numerics::linear::Solution<double>& solution() const& noexcept {
        return solution_;
    }
    [[nodiscard]] const numerics::linear::Solution<double>& solution() const&& = delete;
    [[nodiscard]] core::Result<CirculationResidualRow> row(std::size_t) const noexcept;
    [[nodiscard]] const ProjectedCaseUniformFlow& flow() const& noexcept {
        return flow_;
    }
    [[nodiscard]] const ProjectedCaseUniformFlow& flow() const&& = delete;
    [[nodiscard]] const aero::detail::PreparedTopologyOwner& source() const& noexcept {
        return source_;
    }
    [[nodiscard]] const aero::detail::PreparedTopologyOwner& source() const&& = delete;
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return steps_;
    }

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculation(numerics::linear::Solution<double>, ProjectedCaseUniformFlow&&,
                        core::PlannedBuffer<std::optional<CirculationResidualRow>>,
                        core::PlannedBuffer<NativeCirculationRows>,
                        aero::detail::PreparedTopologyOwner, std::uint64_t) noexcept;
    numerics::linear::Solution<double> solution_;
    ProjectedCaseUniformFlow flow_;
    core::PlannedBuffer<std::optional<CirculationResidualRow>> rows_;
    core::PlannedBuffer<NativeCirculationRows> operator_rows_;
    aero::detail::PreparedTopologyOwner source_;
    std::uint64_t steps_;
};

/// Synchronous replay witness for a qualified real operator. The P1 route compares
/// dense reference coefficients. The closed-contour route replays the same projection
/// and compares retained circulation and residual words without a dense interval copy.
struct CirculationRowReplay {
    std::size_t rows, coefficients;
    std::uint64_t replay_row_bytes, reference_dense_bytes, allocated_bytes, used_steps;
};

[[nodiscard]] core::Result<CirculationRowReplay> verify_discrete_circulation_row_replay(
    const DiscreteNonpenetrationInput&, const DiscreteCirculation&,
    const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest&,
    aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
    aero::detail::P1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

/// Assemble and solve all bound rows, then check the actual returned circulation
/// with the native original residual. Trial input circulation/direction is not
/// used as a starting guess. Consume the flow only on complete success.
[[nodiscard]] core::Result<DiscreteCirculation>
solve_discrete_circulation(const DiscreteNonpenetrationInput&, ProjectedCaseUniformFlow&&,
                           const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest&,
                           aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
                           aero::detail::P1Workspace&, core::ResourcePlan&,
                           const std::stop_token& = {}) noexcept;
/// The coefficient transpose uses the real dual pairing for a dimensionless
/// scalar: its seed has units s/m² and its solution has units s/m.
enum class CirculationSense : std::uint8_t { Tangent, Transpose };

struct CirculationLinearRhs {
    CirculationLinearRhs(aero::detail::PhysicalFaceId entity,
                         const aero::detail::ScalarWords& words) noexcept
        : face(entity), value(words) {}
    aero::detail::PhysicalFaceId face;
    aero::detail::ScalarWords value; ///< Fixed-circulation δR (m/s), or the dual seed (s/m²).
};

struct CirculationLinearRow {
    aero::detail::PhysicalFaceId face;
    aero::detail::ScaledEnclosure original_residual; ///< m/s in Tangent, s/m² in Transpose.
};

/// A certified solve on the fixed finite-discrete operator. The retained solution
/// and row enclosures own their data independently of the primal and RHS.
class DiscreteCirculationLinearization {
public:
    DiscreteCirculationLinearization(DiscreteCirculationLinearization&&) noexcept = default;
    DiscreteCirculationLinearization&
    operator=(DiscreteCirculationLinearization&&) noexcept = default;
    DiscreteCirculationLinearization(const DiscreteCirculationLinearization&) = delete;
    DiscreteCirculationLinearization& operator=(const DiscreteCirculationLinearization&) = delete;
    [[nodiscard]] const numerics::linear::Solution<double>& solution() const& noexcept {
        return solution_;
    }
    [[nodiscard]] const numerics::linear::Solution<double>& solution() const&& = delete;
    [[nodiscard]] core::Result<CirculationLinearRow> row(std::size_t) const noexcept;
    [[nodiscard]] CirculationSense sense() const noexcept {
        return sense_;
    }
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return steps_;
    }

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationLinearization(numerics::linear::Solution<double>,
                                     core::PlannedBuffer<std::optional<CirculationLinearRow>>,
                                     CirculationSense, std::uint64_t) noexcept;
    numerics::linear::Solution<double> solution_;
    core::PlannedBuffer<std::optional<CirculationLinearRow>> rows_;
    CirculationSense sense_;
    std::uint64_t steps_;
};

/// Solve A_eff δΓ = −δR or A_effᵀ λ = seed. The latter pairs with −δR.
/// RHS rows and component scales follow the primal's ascending face IDs. The
/// transpose has distinct matrix and scale generations and its own dual layout IDs.
/// Reusing the primal scale generation requires its certified scale values and units.
/// rhs_absolute uses the RHS's units and matrix_absolute uses 1/m in both senses.
/// max_steps covers native row preparation, projection, solve and original residual.
[[nodiscard]] core::Result<DiscreteCirculationLinearization>
solve_discrete_circulation_linearization(const DiscreteCirculation&, CirculationSense,
                                         std::span<const CirculationLinearRhs>,
                                         const DiscreteCirculationRequest&, std::uint64_t max_steps,
                                         core::ResourcePlan&, const std::stop_token& = {}) noexcept;
struct CirculationCollocationDirection {
    CirculationCollocationDirection(aero::detail::PhysicalFaceId id,
                                    const std::optional<NonpenetrationVector>& query,
                                    const std::array<NonpenetrationVector, 2>& frame,
                                    const NonpenetrationVector& motion) noexcept
        : face(id), free_query(query), tangents(frame), velocity(motion) {}
    aero::detail::PhysicalFaceId face;
    std::optional<NonpenetrationVector> free_query; ///< Present only for a free query.
    std::array<NonpenetrationVector, 2> tangents;   ///< Variations of the original two tangents.
    NonpenetrationVector velocity;                  ///< Explicit variation of wing velocity, m/s.
};

struct DiscreteCirculationDirection {
    DiscreteCirculationDirection(std::span<const aero::detail::PointId> vertices,
                                 std::span<const NonpenetrationVector> positions,
                                 std::span<const CirculationCollocationDirection> collocations,
                                 const NonpenetrationVector& flow) noexcept
        : points(vertices), displacements(positions), rows(collocations), uniform_flow(flow) {}
    std::span<const aero::detail::PointId> points;       ///< All canonical IDs in ascending order.
    std::span<const NonpenetrationVector> displacements; ///< SI position directions in point order.
    std::span<const CirculationCollocationDirection> rows; ///< Ascending bound physical face IDs.
    NonpenetrationVector uniform_flow; ///< Uniform velocity variation at fixed time, m/s.
};

struct CirculationForcingRow {
    aero::detail::PhysicalFaceId face{0};
    aero::detail::ScaledEnclosure residual_direction; ///< Fixed-circulation δR, in m/s.
};

/// Own the original solution, directional forcing, and qualified implicit tangent.
class DiscreteCirculationDirectionalSolve {
public:
    DiscreteCirculationDirectionalSolve(DiscreteCirculationDirectionalSolve&&) noexcept = default;
    DiscreteCirculationDirectionalSolve&
    operator=(DiscreteCirculationDirectionalSolve&&) noexcept = default;
    DiscreteCirculationDirectionalSolve(const DiscreteCirculationDirectionalSolve&) = delete;
    DiscreteCirculationDirectionalSolve&
    operator=(const DiscreteCirculationDirectionalSolve&) = delete;
    [[nodiscard]] const DiscreteCirculation& primal() const& noexcept {
        return primal_;
    }
    [[nodiscard]] const DiscreteCirculation& primal() const&& = delete;
    [[nodiscard]] const DiscreteCirculationLinearization& tangent() const& noexcept {
        return tangent_;
    }
    [[nodiscard]] const DiscreteCirculationLinearization& tangent() const&& = delete;
    [[nodiscard]] core::Result<CirculationForcingRow> forcing(std::size_t) const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return steps_;
    }

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationDirectionalSolve(DiscreteCirculation&&, DiscreteCirculationLinearization&&,
                                        core::PlannedBuffer<std::optional<CirculationForcingRow>>,
                                        std::uint64_t) noexcept;
    DiscreteCirculation primal_;
    DiscreteCirculationLinearization tangent_;
    core::PlannedBuffer<std::optional<CirculationForcingRow>> forcing_;
    std::uint64_t steps_;
};

/// Solve the primal and its fixed-contact-chart direction in one transaction.
/// The strength used for the geometry JVP is the attained circulation, including
/// signed attached wake fronts and fixed history. The unprojected δR enclosure
/// is checked with the attained tangent. Both use limits.steps cumulatively.
/// Time, radius and old wake strengths are fixed. No PG mapping is inferred.
/// The primal and tangent must have distinct RHS and retained-value generations.
/// Consume the supplied flow only once both original residual checks succeed.
[[nodiscard]] core::Result<DiscreteCirculationDirectionalSolve>
solve_discrete_circulation_direction(
    const DiscreteNonpenetrationInput&, const DiscreteCirculationDirection&,
    ProjectedCaseUniformFlow&&, const DiscreteNonpenetrationLimits&,
    const DiscreteCirculationRequest& primal, const DiscreteCirculationRequest& tangent,
    aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
    aero::detail::P1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

using CirculationVectorCovector = std::array<aero::detail::ScaledEnclosure, 3>;

struct CirculationPointCovector {
    aero::detail::PointId point;
    CirculationVectorCovector position; ///< Dual to that canonical SI position.
};

struct CirculationCollocationCovector {
    aero::detail::PhysicalFaceId face;
    std::optional<CirculationVectorCovector> free_query; ///< Present only for a free query.
    std::array<CirculationVectorCovector, 2> tangents;   ///< Duals to the two input tangents.
    CirculationVectorCovector velocity;                  ///< Dual to the face velocity, s/m.
};

/// Pullback for a dimensionless scalar objective with a covector on circulation.
/// All covectors enclose the action of the actual retained transpose solution.
class DiscreteCirculationPullback {
public:
    DiscreteCirculationPullback(DiscreteCirculationPullback&&) noexcept = default;
    DiscreteCirculationPullback& operator=(DiscreteCirculationPullback&&) noexcept = default;
    DiscreteCirculationPullback(const DiscreteCirculationPullback&) = delete;
    DiscreteCirculationPullback& operator=(const DiscreteCirculationPullback&) = delete;
    [[nodiscard]] const DiscreteCirculation& primal() const& noexcept {
        return primal_;
    }
    [[nodiscard]] const DiscreteCirculation& primal() const&& = delete;
    [[nodiscard]] const DiscreteCirculationLinearization& transpose() const& noexcept {
        return transpose_;
    }
    [[nodiscard]] const DiscreteCirculationLinearization& transpose() const&& = delete;
    [[nodiscard]] CirculationVectorCovector uniform_flow() const noexcept {
        return uniform_flow_;
    }
    [[nodiscard]] std::size_t point_count() const noexcept {
        return points_.size();
    }
    [[nodiscard]] std::size_t row_count() const noexcept {
        return rows_.size();
    }
    [[nodiscard]] core::Result<CirculationPointCovector> point(std::size_t) const noexcept;
    [[nodiscard]] core::Result<CirculationCollocationCovector> row(std::size_t) const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return steps_;
    }

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationPullback(DiscreteCirculation&&, DiscreteCirculationLinearization&&,
                                core::PlannedBuffer<std::optional<CirculationPointCovector>>,
                                core::PlannedBuffer<std::optional<CirculationCollocationCovector>>,
                                const CirculationVectorCovector&, std::uint64_t) noexcept;
    DiscreteCirculation primal_;
    DiscreteCirculationLinearization transpose_;
    core::PlannedBuffer<std::optional<CirculationPointCovector>> points_;
    core::PlannedBuffer<std::optional<CirculationCollocationCovector>> rows_;
    CirculationVectorCovector uniform_flow_;
    std::uint64_t steps_;
};

/// Solve A_effᵀ λ=seed in the same discrete chart and return −(dR/dinputs)ᵀ λ
/// using the attained primal and transpose words. Material query contributions
/// accumulate in their canonical point covectors. Time, core radius and inherited
/// strength are fixed. The supplied flow is consumed only on complete success.
[[nodiscard]] core::Result<DiscreteCirculationPullback> solve_discrete_circulation_pullback(
    const DiscreteNonpenetrationInput&, std::span<const CirculationLinearRhs> seed,
    ProjectedCaseUniformFlow&&, const DiscreteNonpenetrationLimits&,
    const DiscreteCirculationRequest& primal, const DiscreteCirculationRequest& transpose,
    aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
    aero::detail::P1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

using DiscreteLoadVector = std::array<aero::detail::ScaledEnclosure, 3>;

struct EdgeMidpointVelocity {
    EdgeMidpointVelocity(aero::detail::EdgeId id, const NonpenetrationVector& v) noexcept
        : edge(id), velocity(v) {}
    aero::detail::EdgeId edge;     ///< Canonical non-TE bound edge.
    NonpenetrationVector velocity; ///< Physical velocity at its exact midpoint, m/s.
};

struct BoundCirculationRate {
    BoundCirculationRate(aero::detail::PhysicalFaceId id,
                         const aero::detail::ScalarWords& value) noexcept
        : face(id), rate(value) {}
    aero::detail::PhysicalFaceId face;
    aero::detail::ScalarWords rate; ///< Externally supplied Γdot, m²/s², for this discrete face.
};

class FixedBdf2Circulation;

enum class DiscreteLoadRetention : std::uint8_t { Aggregates, OriginalRows };

struct DiscreteVectorLoadInput {
    DiscreteVectorLoadInput(aero::detail::P1ActionStamp trial,
                            core::NumericalGeneration circulation_generation,
                            std::span<const EdgeMidpointVelocity> edge_velocities,
                            std::span<const BoundCirculationRate> circulation_rates,
                            const aero::detail::ScalarWords& rho,
                            const NonpenetrationVector& origin) noexcept
        : stamp(trial),
          circulation(circulation_generation),
          edges(edge_velocities),
          rates(circulation_rates),
          density(rho),
          moment_origin(origin) {}
    aero::detail::P1ActionStamp stamp;
    core::NumericalGeneration circulation; ///< Generation of the accompanying circulation trial.
    std::span<const EdgeMidpointVelocity> edges; ///< Ascending canonical force edge IDs.
    std::span<const BoundCirculationRate> rates; ///< Ascending physical bound face IDs.
    aero::detail::ScalarWords density;           ///< kg/m³, finite and strictly positive.
    NonpenetrationVector moment_origin;          ///< Reference point in the geometry's SI frame.
    const FixedBdf2Circulation* fixed_bdf2 = nullptr; ///< Alternative to independent rates.
    DiscreteLoadRetention retention = DiscreteLoadRetention::Aggregates;
};

struct DiscreteCirculatoryLoad {
    aero::detail::EdgeId edge;
    std::array<aero::detail::PointId, 2> endpoints;
    DiscreteLoadVector position, force;
    std::uint64_t excluded_bound_occurrences;
    bool includes_wake_quadrature;
};

struct DiscreteNoncirculatoryLoad {
    aero::detail::PhysicalFaceId face;
    aero::detail::FacetId facet;
    std::array<aero::detail::PointId, 3> vertices;
    DiscreteLoadVector position, force;
};

struct DiscreteNodalLoad {
    aero::detail::PointId point;
    DiscreteLoadVector force;
};

/// One already weighted contribution to an original facet vertex. The canonical
/// point is for global work checks; face, facet and local vertex select the original row.
struct DiscreteOriginalLoadContribution {
    aero::detail::PhysicalFaceId face;
    aero::detail::FacetId facet;
    std::uint8_t local_vertex{};
    aero::detail::PointId canonical;
    std::optional<aero::detail::EdgeId> edge; ///< Empty for triangle area loads.
    DiscreteLoadVector force;
};

struct RetainedBdf2Rate {
    aero::detail::PhysicalFaceId face{0};
    aero::detail::ScaledEnclosure rate; ///< Same D_h Gamma used in the load, m²/s².
};

/// Exact adopted trial identity retained with assembled nodal forces.
struct DiscreteLoadSource {
    aero::detail::GeometryContext context;
    std::array<std::byte, 32> geometry;
    aero::detail::P1ActionStamp stamp;
};

/// Own one certified circulation and its vector loads, in SI units.
class DiscreteCirculationLoads {
public:
    DiscreteCirculationLoads(DiscreteCirculationLoads&&) noexcept = default;
    DiscreteCirculationLoads& operator=(DiscreteCirculationLoads&&) noexcept = default;
    DiscreteCirculationLoads(const DiscreteCirculationLoads&) = delete;
    DiscreteCirculationLoads& operator=(const DiscreteCirculationLoads&) = delete;
    [[nodiscard]] const DiscreteCirculation& primal() const& noexcept {
        return primal_;
    }
    [[nodiscard]] const DiscreteCirculation& primal() const&& = delete;
    [[nodiscard]] DiscreteLoadVector total_force() const noexcept {
        return force_;
    }
    [[nodiscard]] DiscreteLoadVector total_moment() const noexcept {
        return moment_;
    }
    [[nodiscard]] DiscreteLoadVector moment_origin() const noexcept {
        return origin_;
    }
    [[nodiscard]] std::size_t edge_count() const noexcept {
        return edges_.size();
    }
    [[nodiscard]] std::size_t facet_count() const noexcept {
        return facets_.size();
    }
    [[nodiscard]] std::size_t point_count() const noexcept {
        return points_.size();
    }
    [[nodiscard]] core::Result<DiscreteCirculatoryLoad> edge(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteNoncirculatoryLoad> facet(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteNodalLoad> point(std::size_t) const noexcept;
    [[nodiscard]] std::size_t original_count() const noexcept;
    [[nodiscard]] core::Result<DiscreteOriginalLoadContribution>
    original(std::size_t) const noexcept;
    [[nodiscard]] std::size_t bdf2_rate_count() const noexcept {
        return bdf2_rates_.size();
    }
    [[nodiscard]] core::Result<RetainedBdf2Rate> bdf2_rate(std::size_t) const noexcept;
    [[nodiscard]] DiscreteLoadSource source() const noexcept {
        return source_;
    }
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return steps_;
    }

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationLoads(DiscreteCirculation&&,
                             core::PlannedBuffer<std::optional<DiscreteCirculatoryLoad>>,
                             core::PlannedBuffer<std::optional<DiscreteNoncirculatoryLoad>>,
                             core::PlannedBuffer<std::optional<DiscreteNodalLoad>>,
                             const DiscreteLoadVector& force, const DiscreteLoadVector& moment,
                             const DiscreteLoadVector& origin,
                             core::PlannedBuffer<std::optional<RetainedBdf2Rate>>,
                             core::PlannedBuffer<std::optional<DiscreteOriginalLoadContribution>>,
                             std::size_t original_count, DiscreteLoadSource,
                             std::uint64_t) noexcept;
    DiscreteCirculation primal_;
    core::PlannedBuffer<std::optional<DiscreteCirculatoryLoad>> edges_;
    core::PlannedBuffer<std::optional<DiscreteNoncirculatoryLoad>> facets_;
    core::PlannedBuffer<std::optional<DiscreteNodalLoad>> points_;
    DiscreteLoadVector force_, moment_, origin_;
    core::PlannedBuffer<std::optional<RetainedBdf2Rate>> bdf2_rates_;
    core::PlannedBuffer<std::optional<DiscreteOriginalLoadContribution>> original_;
    std::size_t original_count_{};
    DiscreteLoadSource source_;
    std::uint64_t steps_;
};

/// Evaluate adopted midpoint and triangle loads on the same solved discrete trial.
/// The current geometry must have a positive regularization radius for incident
/// midpoint queries. The caller supplies the rate and midpoint velocities;
/// the same state, exact edge IDs and face order are required. No history is
/// reconstructed. The flow is consumed only on complete success.
[[nodiscard]] core::Result<DiscreteCirculationLoads>
solve_discrete_circulation_loads(const DiscreteNonpenetrationInput&, const DiscreteVectorLoadInput&,
                                 ProjectedCaseUniformFlow&&, const DiscreteNonpenetrationLimits&,
                                 const DiscreteCirculationRequest&,
                                 aero::detail::GeometryWorkspace&, aero::detail::SharedP1Workspace&,
                                 aero::detail::LocalPacketWorkspace&, aero::detail::P1Workspace&,
                                 core::ResourcePlan&, const std::stop_token& = {}) noexcept;

struct BoundLoadDirection {
    BoundLoadDirection(aero::detail::PhysicalFaceId id, const aero::detail::ScalarWords& gamma,
                       const aero::detail::ScalarWords& rate) noexcept
        : face(id), circulation(gamma), circulation_rate(rate) {}
    aero::detail::PhysicalFaceId face;
    aero::detail::ScalarWords circulation, circulation_rate; ///< m²/s and m²/s² directions.
};

/// Fixed-chart load direction at the same nominal circulation and trial.
struct DiscreteVectorLoadDirection {
    DiscreteVectorLoadDirection(aero::detail::P1ActionStamp trial,
                                core::NumericalGeneration circulation_generation,
                                std::span<const aero::detail::PointId> point_ids,
                                std::span<const NonpenetrationVector> dx,
                                std::span<const EdgeMidpointVelocity> dv,
                                std::span<const BoundLoadDirection> gamma,
                                const NonpenetrationVector& ambient) noexcept
        : stamp(trial),
          circulation(circulation_generation),
          points(point_ids),
          displacements(dx),
          edges(dv),
          bound(gamma),
          uniform_flow(ambient) {}
    aero::detail::P1ActionStamp stamp;
    core::NumericalGeneration circulation;
    std::span<const aero::detail::PointId> points;       ///< Ascending canonical geometry IDs.
    std::span<const NonpenetrationVector> displacements; ///< SI point position directions.
    std::span<const EdgeMidpointVelocity>
        edges; ///< Midpoint-velocity directions in value edge order.
    std::span<const BoundLoadDirection> bound; ///< Ascending bound face IDs.
    NonpenetrationVector uniform_flow;         ///< Uniform velocity direction at fixed time, m/s.
};

struct DiscreteLoadSampleDirection {
    DiscreteLoadVector position, force; ///< Direction of the same indexed value load sample.
};

/// Value loads and the matching forward direction. Indices follow the value result.
class DiscreteCirculationLoadJvp {
public:
    DiscreteCirculationLoadJvp(DiscreteCirculationLoadJvp&&) noexcept = default;
    DiscreteCirculationLoadJvp& operator=(DiscreteCirculationLoadJvp&&) noexcept = default;
    DiscreteCirculationLoadJvp(const DiscreteCirculationLoadJvp&) = delete;
    DiscreteCirculationLoadJvp& operator=(const DiscreteCirculationLoadJvp&) = delete;
    [[nodiscard]] const DiscreteCirculationLoads& loads() const& noexcept {
        return loads_;
    }
    [[nodiscard]] const DiscreteCirculationLoads& loads() const&& = delete;
    [[nodiscard]] DiscreteLoadVector total_force_direction() const noexcept {
        return force_;
    }
    [[nodiscard]] DiscreteLoadVector total_moment_direction() const noexcept {
        return moment_;
    }
    [[nodiscard]] core::Result<DiscreteLoadSampleDirection> edge(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadSampleDirection> facet(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteNodalLoad> point(std::size_t) const noexcept;

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationLoadJvp(DiscreteCirculationLoads&&,
                               core::PlannedBuffer<std::optional<DiscreteLoadSampleDirection>>,
                               core::PlannedBuffer<std::optional<DiscreteLoadSampleDirection>>,
                               core::PlannedBuffer<std::optional<DiscreteNodalLoad>>,
                               const DiscreteLoadVector& force,
                               const DiscreteLoadVector& moment) noexcept;
    DiscreteCirculationLoads loads_;
    core::PlannedBuffer<std::optional<DiscreteLoadSampleDirection>> edges_, facets_;
    core::PlannedBuffer<std::optional<DiscreteNodalLoad>> points_;
    DiscreteLoadVector force_, moment_;
};

/// Use supplied circulation/rate directions at the attained value circulation.
/// Shared induction requires the complete canonical point order. Time, density,
/// core radius, inherited wake strength, and moment origin remain fixed.
/// Both value and direction must succeed before the flow is consumed.
[[nodiscard]] core::Result<DiscreteCirculationLoadJvp> solve_discrete_circulation_load_direction(
    const DiscreteNonpenetrationInput&, const DiscreteVectorLoadInput&,
    const DiscreteVectorLoadDirection&, ProjectedCaseUniformFlow&&,
    const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest&,
    aero::detail::GeometryWorkspace&, aero::detail::SharedP1Workspace&,
    aero::detail::LocalPacketWorkspace&, aero::detail::P1Workspace&, core::ResourcePlan&,
    const std::stop_token& = {}) noexcept;

/// Shared fixed-chart circulation directions, with independent load-only inputs.
struct DiscreteImplicitVectorLoadDirection {
    DiscreteImplicitVectorLoadDirection(aero::detail::P1ActionStamp trial,
                                        core::NumericalGeneration generation,
                                        const DiscreteCirculationDirection& variation,
                                        std::span<const EdgeMidpointVelocity> midpoint,
                                        std::span<const BoundCirculationRate> rate) noexcept
        : stamp(trial), circulation(generation), state(variation), edges(midpoint), rates(rate) {}
    aero::detail::P1ActionStamp stamp;
    core::NumericalGeneration circulation; ///< Retained primal generation.
    DiscreteCirculationDirection
        state; ///< The points and ambient direction are shared by both stages.
    std::span<const EdgeMidpointVelocity> edges; ///< Midpoint velocity directions, m/s.
    std::span<const BoundCirculationRate> rates; ///< Independent Γdot directions, m²/s².
};

class DiscreteCirculationImplicitLoadJvp {
public:
    DiscreteCirculationImplicitLoadJvp(DiscreteCirculationImplicitLoadJvp&&) noexcept = default;
    DiscreteCirculationImplicitLoadJvp&
    operator=(DiscreteCirculationImplicitLoadJvp&&) noexcept = default;
    DiscreteCirculationImplicitLoadJvp(const DiscreteCirculationImplicitLoadJvp&) = delete;
    DiscreteCirculationImplicitLoadJvp&
    operator=(const DiscreteCirculationImplicitLoadJvp&) = delete;
    [[nodiscard]] const DiscreteCirculationLoadJvp& direction() const& noexcept {
        return direction_;
    }
    [[nodiscard]] const DiscreteCirculationLoadJvp& direction() const&& = delete;
    [[nodiscard]] const DiscreteCirculationLinearization& tangent() const& noexcept {
        return tangent_;
    }
    [[nodiscard]] const DiscreteCirculationLinearization& tangent() const&& = delete;

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationImplicitLoadJvp(DiscreteCirculationLoadJvp&& direction,
                                       DiscreteCirculationLinearization&& tangent) noexcept
        : direction_(std::move(direction)), tangent_(std::move(tangent)) {}
    DiscreteCirculationLoadJvp direction_;
    DiscreteCirculationLinearization tangent_;
};

/// Solve one primal, then A_eff δΓ = −δR using its fixed-chart native forcing.
/// Pass both retained tangent words to the load JVP. The tangent must use the
/// primal operator identity and distinct RHS and retained generations. A shared
/// scale generation requires the certified primal scales. Flow moves on success.
[[nodiscard]] core::Result<DiscreteCirculationImplicitLoadJvp>
solve_discrete_circulation_implicit_load_direction(
    const DiscreteNonpenetrationInput&, const DiscreteVectorLoadInput&,
    const DiscreteImplicitVectorLoadDirection&, ProjectedCaseUniformFlow&&,
    const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest& primal,
    const DiscreteCirculationRequest& tangent, aero::detail::GeometryWorkspace&,
    aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
    aero::detail::P1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

struct DiscreteLoadSampleWeights {
    DiscreteLoadSampleWeights(const NonpenetrationVector& x, const NonpenetrationVector& f) noexcept
        : position(x), force(f) {}
    NonpenetrationVector position, force; ///< 1/m and 1/N for a dimensionless objective.
};

/// Explicit weights on all value load outputs, in the corresponding result order.
struct DiscreteVectorLoadWeights {
    DiscreteVectorLoadWeights(aero::detail::P1ActionStamp trial,
                              core::NumericalGeneration circulation_generation,
                              std::span<const DiscreteLoadSampleWeights> edge_rows,
                              std::span<const DiscreteLoadSampleWeights> facet_rows,
                              std::span<const NonpenetrationVector> point_rows,
                              const NonpenetrationVector& total_force,
                              const NonpenetrationVector& total_moment) noexcept
        : stamp(trial),
          circulation(circulation_generation),
          edges(edge_rows),
          facets(facet_rows),
          points(point_rows),
          force(total_force),
          moment(total_moment) {}
    aero::detail::P1ActionStamp stamp;
    core::NumericalGeneration circulation;
    std::span<const DiscreteLoadSampleWeights> edges, facets;
    std::span<const NonpenetrationVector> points; ///< Weights on the assembled nodal forces, 1/N.
    NonpenetrationVector force, moment;           ///< Weights on the totals, 1/N and 1/(N m).
};

struct LoadEdgeVelocityCovector {
    aero::detail::EdgeId edge;
    DiscreteLoadVector velocity; ///< Dual to midpoint velocity, s/m.
};

struct LoadCirculationCovector {
    aero::detail::PhysicalFaceId face;
    aero::detail::ScaledEnclosure circulation, circulation_rate; ///< s/m² and s²/m².
};

/// Transpose of the explicit load variables at the attained value circulation.
class DiscreteCirculationLoadVjp {
public:
    DiscreteCirculationLoadVjp(DiscreteCirculationLoadVjp&&) noexcept = default;
    DiscreteCirculationLoadVjp& operator=(DiscreteCirculationLoadVjp&&) noexcept = default;
    DiscreteCirculationLoadVjp(const DiscreteCirculationLoadVjp&) = delete;
    DiscreteCirculationLoadVjp& operator=(const DiscreteCirculationLoadVjp&) = delete;
    [[nodiscard]] const DiscreteCirculationLoads& loads() const& noexcept {
        return loads_;
    }
    [[nodiscard]] const DiscreteCirculationLoads& loads() const&& = delete;
    [[nodiscard]] std::size_t point_count() const noexcept {
        return points_.size();
    }
    [[nodiscard]] core::Result<CirculationPointCovector> point(std::size_t) const noexcept;
    [[nodiscard]] core::Result<LoadEdgeVelocityCovector>
    midpoint_velocity(std::size_t) const noexcept;
    [[nodiscard]] core::Result<LoadCirculationCovector> bound(std::size_t) const noexcept;
    [[nodiscard]] DiscreteLoadVector uniform_flow() const noexcept {
        return flow_;
    }

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationLoadVjp(DiscreteCirculationLoads&&,
                               core::PlannedBuffer<std::optional<CirculationPointCovector>>,
                               core::PlannedBuffer<std::optional<LoadEdgeVelocityCovector>>,
                               core::PlannedBuffer<std::optional<LoadCirculationCovector>>,
                               const DiscreteLoadVector& flow) noexcept;
    DiscreteCirculationLoads loads_;
    core::PlannedBuffer<std::optional<CirculationPointCovector>> points_;
    core::PlannedBuffer<std::optional<LoadEdgeVelocityCovector>> edges_;
    core::PlannedBuffer<std::optional<LoadCirculationCovector>> bound_;
    DiscreteLoadVector flow_;
};

/// Pull all supplied load duals back to the fixed-chart inputs used by the JVP.
/// Point covectors are ascending canonical IDs. With no midpoint evaluations,
/// only points of the loaded triangles can receive covectors.
/// Time, radius, density, moment origin and inherited wake strength are fixed.
/// Consume the flow only after the value and pullback both succeed.
[[nodiscard]] core::Result<DiscreteCirculationLoadVjp> solve_discrete_circulation_load_pullback(
    const DiscreteNonpenetrationInput&, const DiscreteVectorLoadInput&,
    const DiscreteVectorLoadWeights&, ProjectedCaseUniformFlow&&,
    const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest&,
    aero::detail::GeometryWorkspace&, aero::detail::SharedP1Workspace&,
    aero::detail::LocalPacketWorkspace&, aero::detail::P1Workspace&, core::ResourcePlan&,
    const std::stop_token& = {}) noexcept;

struct LoadCirculationRateCovector {
    aero::detail::PhysicalFaceId face{0};
    aero::detail::ScaledEnclosure rate; ///< Dual to the independent Γdot, s²/m².
};

/// Fixed-chart transpose of the solved load outputs. Point and flow duals
/// include the response of retained circulation. With fixed BDF2 they include
/// the current-rate chain while the older solves and step size stay fixed.
class DiscreteCirculationImplicitLoadVjp {
public:
    DiscreteCirculationImplicitLoadVjp(DiscreteCirculationImplicitLoadVjp&&) noexcept = default;
    DiscreteCirculationImplicitLoadVjp&
    operator=(DiscreteCirculationImplicitLoadVjp&&) noexcept = default;
    DiscreteCirculationImplicitLoadVjp(const DiscreteCirculationImplicitLoadVjp&) = delete;
    DiscreteCirculationImplicitLoadVjp&
    operator=(const DiscreteCirculationImplicitLoadVjp&) = delete;
    [[nodiscard]] const DiscreteCirculationLoads& loads() const& noexcept {
        return values_.loads();
    }
    [[nodiscard]] const DiscreteCirculationLoads& loads() const&& = delete;
    [[nodiscard]] const DiscreteCirculationLinearization& transpose() const& noexcept {
        return transpose_;
    }
    [[nodiscard]] const DiscreteCirculationLinearization& transpose() const&& = delete;
    [[nodiscard]] DiscreteLoadVector uniform_flow() const noexcept {
        return values_.uniform_flow();
    }
    [[nodiscard]] std::size_t point_count() const noexcept {
        return values_.point_count();
    }
    [[nodiscard]] core::Result<CirculationPointCovector> point(std::size_t i) const noexcept {
        return values_.point(i);
    }
    [[nodiscard]] core::Result<LoadEdgeVelocityCovector>
    midpoint_velocity(std::size_t i) const noexcept {
        return values_.midpoint_velocity(i);
    }
    [[nodiscard]] core::Result<LoadCirculationRateCovector>
    circulation_rate(std::size_t) const noexcept;
    [[nodiscard]] core::Result<CirculationCollocationCovector>
    collocation(std::size_t) const noexcept;

private:
    friend struct NonpenetrationBuilder;
    DiscreteCirculationImplicitLoadVjp(
        DiscreteCirculationLoadVjp&& values, DiscreteCirculationLinearization&& transpose,
        core::PlannedBuffer<std::optional<CirculationCollocationCovector>>&& rows) noexcept
        : values_(std::move(values)), transpose_(std::move(transpose)), rows_(std::move(rows)) {}
    DiscreteCirculationLoadVjp values_;
    DiscreteCirculationLinearization transpose_;
    core::PlannedBuffer<std::optional<CirculationCollocationCovector>> rows_;
};

/// Solve the primal once and solve A_effᵀ λ = ∂L/∂Γ against its native influence
/// and the native load seed. Add −(dR/dinputs)ᵀ λ to the direct load duals.
/// Each transpose layout and its matrix, scale, RHS and retained generations
/// must be distinct from the corresponding primal data. Flow is consumed on success.
[[nodiscard]] core::Result<DiscreteCirculationImplicitLoadVjp>
solve_discrete_circulation_implicit_load_pullback(
    const DiscreteNonpenetrationInput&, const DiscreteVectorLoadInput&,
    const DiscreteVectorLoadWeights&, ProjectedCaseUniformFlow&&,
    const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest& primal,
    const DiscreteCirculationRequest& transpose, aero::detail::GeometryWorkspace&,
    aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
    aero::detail::P1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

struct DiscreteKnownFlowInput {
    const CaseKnownFlowReady& source;
    const CaseKnownFlowReady& current;
    CaseUniformFlowQuery time;
    CaseKnownFlowProjection projection;
    CaseKnownFlowLimits flow_limits; ///< Provider caps per collocation.
    aero::detail::GeometryLimits coordinate_limits;
    aero::detail::GeometryWorkspace& coordinates;
    aero::detail::UniformFlowWorkspace& workspace;
};

struct KnownCirculationTangent {
    const DiscreteCirculationDirection& state;
    aero::detail::ScalarWords flow_time; ///< Variation of ambient time only, in seconds.
    const DiscreteCirculationRequest& request;
};

struct KnownCirculationAdjoint {
    std::span<const CirculationLinearRhs> seed;
    const DiscreteCirculationRequest& request;
};

using KnownCirculationDerivative =
    std::variant<std::monostate, KnownCirculationTangent, KnownCirculationAdjoint>;
struct KnownCirculationRun;

/// Real discrete solve with the known flow sampled once at each qualified collocation.
/// Holds the native original residual, input jets, and the selected fixed-chart derivative.
class DiscreteKnownCirculation {
public:
    DiscreteKnownCirculation(DiscreteKnownCirculation&&) noexcept;
    DiscreteKnownCirculation& operator=(DiscreteKnownCirculation&&) noexcept;
    ~DiscreteKnownCirculation();
    DiscreteKnownCirculation(const DiscreteKnownCirculation&) = delete;
    DiscreteKnownCirculation& operator=(const DiscreteKnownCirculation&) = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const numerics::linear::Solution<double>>>
    solution() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const numerics::linear::Solution<double>>>
    solution() const&& = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteCirculationLinearization>>
    linearization() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteCirculationLinearization>>
    linearization() const&& = delete;
    [[nodiscard]] core::Result<CirculationResidualRow> row(std::size_t) const noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseKnownFlowJet>>
    flow(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseKnownFlowJet>>
    flow(std::size_t) const&& = delete;
    [[nodiscard]] core::Result<CirculationForcingRow> forcing(std::size_t) const noexcept;
    [[nodiscard]] core::Result<CirculationPointCovector> point(std::size_t) const noexcept;
    [[nodiscard]] core::Result<CirculationCollocationCovector>
    collocation(std::size_t) const noexcept;
    [[nodiscard]] core::Result<CirculationVectorCovector> uniform_component_dual() const noexcept;
    /// Covector to a common ambient time shift, in inverse seconds.
    [[nodiscard]] core::Result<aero::detail::ScaledEnclosure> flow_time_dual() const noexcept;

private:
    friend struct NonpenetrationBuilder;
    explicit DiscreteKnownCirculation(core::PlannedBuffer<KnownCirculationRun>&&) noexcept;
    core::PlannedBuffer<KnownCirculationRun> data_;
};

/// Requires at least one bound row. Geometry, attached-wake chart, core radius and wake history
/// remain those of this trial.
/// The time direction differentiates only the analytic known-flow input at its requested side.
[[nodiscard]] core::Result<DiscreteKnownCirculation> solve_discrete_known_circulation(
    const DiscreteNonpenetrationInput&, const DiscreteKnownFlowInput&,
    const KnownCirculationDerivative&, const DiscreteNonpenetrationLimits&,
    const DiscreteCirculationRequest&, aero::detail::SharedP1Workspace&,
    aero::detail::LocalPacketWorkspace&, aero::detail::P1Workspace&, core::ResourcePlan&,
    const std::stop_token& = {}) noexcept;

struct KnownLoadTangent {
    const DiscreteImplicitVectorLoadDirection& state;
    aero::detail::ScalarWords flow_time; ///< Ambient-time direction only, in seconds.
    const DiscreteCirculationRequest& request;
};
struct KnownLoadAdjoint {
    const DiscreteVectorLoadWeights& weights;
    const DiscreteCirculationRequest& request;
};
using KnownLoadDerivative = std::variant<std::monostate, KnownLoadTangent, KnownLoadAdjoint>;

struct KnownLoadRun;

/// Own the solved circulation, the flow at all used midpoints, and the matching loads.
class DiscreteKnownLoads {
public:
    DiscreteKnownLoads(DiscreteKnownLoads&&) noexcept;
    DiscreteKnownLoads& operator=(DiscreteKnownLoads&&) noexcept;
    DiscreteKnownLoads(const DiscreteKnownLoads&) = delete;
    DiscreteKnownLoads& operator=(const DiscreteKnownLoads&) = delete;
    ~DiscreteKnownLoads();
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteKnownCirculation>>
    circulation() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteKnownCirculation>>
    circulation() const&& = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseKnownFlowJet>>
    midpoint_flow(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseKnownFlowJet>>
    midpoint_flow(std::size_t) const&& = delete;
    [[nodiscard]] core::Result<DiscreteCirculatoryLoad> edge(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteNoncirculatoryLoad> facet(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteNodalLoad> point(std::size_t) const noexcept;
    [[nodiscard]] std::size_t original_count() const noexcept;
    [[nodiscard]] core::Result<DiscreteOriginalLoadContribution>
    original(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadVector> total_force() const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadVector> total_moment() const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadVector> moment_origin() const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadSampleDirection>
    edge_direction(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadSampleDirection>
    facet_direction(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteNodalLoad> point_direction(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadVector> force_direction() const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadVector> moment_direction() const noexcept;
    [[nodiscard]] core::Result<CirculationPointCovector> point_dual(std::size_t) const noexcept;
    [[nodiscard]] core::Result<LoadEdgeVelocityCovector>
    midpoint_velocity_dual(std::size_t) const noexcept;
    /// Dual to independent Γdot. The solved Γ response is included in the other covectors.
    [[nodiscard]] core::Result<LoadCirculationRateCovector>
    circulation_dual(std::size_t) const noexcept;
    [[nodiscard]] core::Result<CirculationCollocationCovector>
    collocation_dual(std::size_t) const noexcept;
    [[nodiscard]] core::Result<DiscreteLoadVector> uniform_component_dual() const noexcept;
    [[nodiscard]] core::Result<aero::detail::ScaledEnclosure> flow_time_dual() const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend struct NonpenetrationBuilder;
    explicit DiscreteKnownLoads(core::PlannedBuffer<KnownLoadRun>&&) noexcept;
    core::PlannedBuffer<KnownLoadRun> data_;
};

/// One nominal circulation solve, then fixed-chart loads and the selected implicit action.
/// Circulation rate and midpoint wing velocity remain explicit independent inputs.
[[nodiscard]] core::Result<DiscreteKnownLoads>
solve_discrete_known_loads(const DiscreteNonpenetrationInput&, const DiscreteVectorLoadInput&,
                           const DiscreteKnownFlowInput&, const KnownLoadDerivative&,
                           const DiscreteNonpenetrationLimits&, const DiscreteCirculationRequest&,
                           aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
                           aero::detail::P1Workspace&, core::ResourcePlan&,
                           const std::stop_token& = {}) noexcept;

struct DiscreteFreeWakeDirection {
    std::span<const aero::detail::PointId> points; ///< Ascending canonical point IDs.
    std::span<const NonpenetrationVector> displacements;
    std::span<const FrozenWakeTrace> history; ///< Old trace directions in source ID order.
    NonpenetrationVector uniform_flow;
    aero::detail::ScalarWords flow_time; ///< Only the input field's time coordinate.
    aero::detail::ScalarWords core_radius;
};

struct DiscreteFreeWakeLimits {
    DiscreteNonpenetrationLimits action;
    std::uint64_t max_particles{};
};

struct FreeWakeVelocitySource {
    aero::detail::GeometryContext geometry;
    aero::detail::P1ActionStamp stamp;
    core::AcceptedStateId parent{0};
};

struct FreeWakeVelocityRow {
    aero::detail::PointId particle{0};
    aero::detail::PointId owner{0};
    DiscreteLoadVector velocity;
    std::optional<DiscreteLoadVector> direction;
    std::uint64_t source_facets{}, far_facets{}, cancelled_analytic_pairs{};
    bool includes_wake_quadrature{};
};

struct FreeWakeStrengthCovector {
    aero::detail::P1StrengthId id; ///< Bound circulation or independent old wake trace.
    aero::detail::ScaledEnclosure value;
};
struct FreeWakeFieldCovector {
    DiscreteLoadVector uniform_flow;
    aero::detail::ScaledEnclosure flow_time, core_radius;
};
struct DiscreteFreeWakeRun;

class DiscreteFreeWakeVelocity {
public:
    DiscreteFreeWakeVelocity(DiscreteFreeWakeVelocity&&) noexcept;
    DiscreteFreeWakeVelocity& operator=(DiscreteFreeWakeVelocity&&) noexcept;
    DiscreteFreeWakeVelocity(const DiscreteFreeWakeVelocity&) = delete;
    DiscreteFreeWakeVelocity& operator=(const DiscreteFreeWakeVelocity&) = delete;
    ~DiscreteFreeWakeVelocity();
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;
    [[nodiscard]] core::Result<FreeWakeVelocitySource> source() const noexcept;
    [[nodiscard]] core::Result<FreeWakeVelocityRow> row(std::size_t) const noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseKnownFlowJet>>
    ambient(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseKnownFlowJet>>
    ambient(std::size_t) const&& = delete;
    [[nodiscard]] std::size_t point_count() const noexcept;
    [[nodiscard]] std::size_t strength_count() const noexcept;
    [[nodiscard]] core::Result<CirculationPointCovector> point_dual(std::size_t) const noexcept;
    [[nodiscard]] core::Result<FreeWakeStrengthCovector> strength_dual(std::size_t) const noexcept;
    [[nodiscard]] core::Result<FreeWakeFieldCovector> field_dual() const noexcept;

private:
    friend struct NonpenetrationBuilder;
    explicit DiscreteFreeWakeVelocity(core::PlannedBuffer<DiscreteFreeWakeRun>&&) noexcept;
    core::PlannedBuffer<DiscreteFreeWakeRun> data_;
};

/// Compose the velocity pullback directly into a signed residual.
enum class FreeWakeSeedScale : std::uint8_t { Positive, Negative, NegativeHalf };

/// Use a fixed contact-qualified source state and its positive common SI radius.
/// Each selected particle belongs to a retained wake face. The listed IDs can be
/// explicit aliases of wake points, with one evaluation per canonical owner.
/// The circulation words (and directions, when requested) are the supplied
/// bound rows. The attached trace follows T; history entries remain independent.
/// Seeds are empty for no adjoint, otherwise one real velocity seed per particle.
/// No material identity is released or advanced here.
[[nodiscard]] core::Result<DiscreteFreeWakeVelocity> evaluate_discrete_free_wake_velocity(
    const DiscreteNonpenetrationInput&, const DiscreteKnownFlowInput&,
    std::span<const aero::detail::PointId> particles, const DiscreteFreeWakeDirection*,
    std::span<const NonpenetrationVector> seeds, const DiscreteFreeWakeLimits&,
    aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
    aero::detail::P1Workspace&, core::ResourcePlan&, const std::stop_token& = {},
    FreeWakeSeedScale = FreeWakeSeedScale::Positive) noexcept;

} // namespace uvlas::coupling::detail
