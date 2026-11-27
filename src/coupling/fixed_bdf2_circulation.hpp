#pragma once

#include "case_flow_clock.hpp"

namespace uvlas::aero::detail::p1 {
struct Work;
}

namespace uvlas::coupling::detail {

/// A checked same-interval stencil of numerical circulation solves. The two
/// prior values own the actual retained solution words and their physical IDs.
/// Only the time-step owner can establish that those historical solves were accepted.
class FixedBdf2Circulation {
public:
    FixedBdf2Circulation(FixedBdf2Circulation&&) noexcept;
    FixedBdf2Circulation& operator=(FixedBdf2Circulation&&) noexcept;
    FixedBdf2Circulation(const FixedBdf2Circulation&) = delete;
    FixedBdf2Circulation& operator=(const FixedBdf2Circulation&) = delete;
    ~FixedBdf2Circulation();

    /// Read the newest to oldest checked time nodes. Copy only fully qualified
    /// circulation values, and take the destination node on success.
    [[nodiscard]] static core::Result<FixedBdf2Circulation>
    prepare(CaseClockedFlowSample&& destination, const CaseClockedFlowSample& previous,
            const CaseClockedFlowSample& older, const DiscreteCirculation& previous_solution,
            const DiscreteCirculation& older_solution, std::uint64_t maximum_integer_bits,
            std::uint64_t maximum_steps, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

    /// Associate numerically qualified known-flow histories with the trial owner's parent ID.
    /// The owner supplies accepted-state provenance. This operation qualifies Case, clock,
    /// topology and retained values without publishing an accepted state.
    [[nodiscard]] static core::Result<FixedBdf2Circulation>
    prepare(CaseClockedKnownFlowSample&& destination, const CaseClockedKnownFlowSample& previous,
            const CaseClockedKnownFlowSample& older,
            const DiscreteKnownCirculation& previous_solution,
            const DiscreteKnownCirculation& older_solution, core::AcceptedStateId parent,
            std::uint64_t maximum_integer_bits, std::uint64_t maximum_steps, core::ResourcePlan&,
            const std::stop_token& = {}) noexcept;

private:
    friend struct NonpenetrationBuilder;
    struct Data;
    explicit FixedBdf2Circulation(core::PlannedBuffer<Data>&&) noexcept;
    core::PlannedBuffer<Data> data_;

    template <class Sample, class Solve>
    [[nodiscard]] static core::Result<FixedBdf2Circulation>
    prepare_history(Sample&& destination, const Sample& previous, const Sample& older,
                    const Solve& previous_solution, const Solve& older_solution,
                    std::optional<core::AcceptedStateId> parent, std::uint64_t maximum_integer_bits,
                    std::uint64_t maximum_steps, core::ResourcePlan&,
                    const std::stop_token&) noexcept;

    [[nodiscard]] const Data& validate_source(const aero::detail::PreparedTopologyOwner&,
                                              core::NumericalGeneration, std::size_t,
                                              aero::detail::p1::Work&) const;
    /// These internal actions use the load routine's same bounded arithmetic.
    void bind(const aero::detail::PreparedTopologyOwner&, const ProjectedCaseUniformFlow&,
              core::NumericalGeneration, std::size_t, aero::detail::p1::Work&) const;
    void bind(const aero::detail::PreparedTopologyOwner&, const DiscreteKnownFlowInput&,
              core::AcceptedStateId, core::NumericalGeneration, std::size_t,
              aero::detail::p1::Work&) const;
    [[nodiscard]] aero::detail::p1::Interval value(std::size_t, aero::detail::PhysicalFaceId,
                                                   const aero::detail::p1::Interval&,
                                                   aero::detail::p1::Work&) const;
    [[nodiscard]] aero::detail::p1::Interval current_action(const aero::detail::p1::Interval&,
                                                            aero::detail::p1::Work&) const;
};

} // namespace uvlas::coupling::detail
