#pragma once

#include "aero/p1_shared_action.hpp"
#include "case_uniform_flow_projection.hpp"

namespace uvlas::coupling::detail {
using NonpenetrationVector = std::array<aero::detail::ScalarWords, 3>;

/// Explicit discrete collocation in the current geometry. The tangents orient its
/// normal and the velocity is in m/s. No physical PG mapping is inferred.
struct DiscreteCollocation {
    DiscreteCollocation(aero::detail::FacetId source, const aero::detail::ResolvedQuery& point,
                        const std::array<NonpenetrationVector, 2>& frame,
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
    std::array<NonpenetrationVector, 2> tangents;
    NonpenetrationVector velocity;
    aero::detail::ScalarWords circulation;
    aero::detail::ScalarWords circulation_direction;
};

struct AttachedWakeFront {};
struct DiscreteWakePanel {
    DiscreteWakePanel(aero::detail::FacetId source,
                      std::variant<AttachedWakeFront, aero::detail::P1StrengthId> left,
                      aero::detail::P1StrengthId right) noexcept
        : facet(source), leading(left), trailing(right) {}
    aero::detail::FacetId facet; ///< Either child of the retained wake quadrilateral.
    std::variant<AttachedWakeFront, aero::detail::P1StrengthId> leading;
    aero::detail::P1StrengthId trailing;
};

struct FrozenWakeTrace {
    aero::detail::P1StrengthId id;
    aero::detail::ScalarWords strength;
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
    aero::detail::PhysicalFaceId face;
    aero::detail::ScaledEnclosure residual;           ///< m/s, including the original affine terms.
    aero::detail::ScaledEnclosure circulation_action; ///< A_eff δΓ in m/s.
    bool includes_wake_quadrature; ///< The finite wake formula carries the existing P1 diagnostic.
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
                           std::optional<aero::detail::P1ActionValue>,
                           aero::detail::GeometryContext, aero::detail::P1ActionStamp,
                           std::uint64_t) noexcept;
    core::PlannedBuffer<NonpenetrationStorage> rows_;
    ProjectedCaseUniformFlow flow_;
    std::optional<aero::detail::P1ActionValue> geometry_owner_;
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
} // namespace uvlas::coupling::detail
