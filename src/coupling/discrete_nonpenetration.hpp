#pragma once

#include "aero/p1_shared_action.hpp"
#include "case_uniform_flow_projection.hpp"

namespace uvlas::coupling::detail {
using NonpenetrationVector = std::array<aero::detail::ScalarWords, 3>;

/// Explicit discrete collocation in the current geometry. The tangents orient its
/// normal and the velocity is in m/s. No physical PG mapping is inferred.
struct DiscreteCollocation {
    aero::detail::FacetId facet; ///< Either child of the bound physical quadrilateral.
    const aero::detail::ResolvedQuery& query;
    std::array<NonpenetrationVector, 2> tangents;
    NonpenetrationVector velocity;
    aero::detail::ScalarWords circulation;
    aero::detail::ScalarWords circulation_direction;
};

struct AttachedWakeFront {};
struct DiscreteWakePanel {
    aero::detail::FacetId facet; ///< Either child of the retained wake quadrilateral.
    std::variant<AttachedWakeFront, aero::detail::P1StrengthId> leading;
    aero::detail::P1StrengthId trailing;
};

struct FrozenWakeTrace {
    aero::detail::P1StrengthId id;
    aero::detail::ScalarWords strength;
};
struct FrozenWakeTraces {
    core::AcceptedStateId parent;
    std::span<const FrozenWakeTrace> columns; ///< Strictly increasing independent history IDs.
};

struct DiscreteNonpenetrationInput {
    const aero::detail::GlobalContactState& contact;
    std::span<const DiscreteCollocation> collocations;
    std::span<const DiscreteWakePanel> wake;
    FrozenWakeTraces history;
    aero::detail::P1ActionStamp trial_stamp;
    aero::detail::pchip::PchipInputStamp flow_stamp;
    aero::detail::pchip::RequestedSide flow_side;
};

struct DiscreteNonpenetrationLimits {
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
