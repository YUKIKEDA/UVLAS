#pragma once

#include <uvlas/structure/pg_sparse_chart.hpp>

namespace uvlas::structure {

struct PgSupportTimeKnot {
    const ExactRational& seconds;
    std::span<const PgCartesianExactPosition> value, rate; ///< All of a support's original rows.
};

struct PgSupportTimeTrack {
    PgConstraintId id;
    std::span<const PgSupportTimeKnot> knots;
};

struct PgSupportTimeWindow {
    const ExactRational& begin_seconds;
    const ExactRational& end_seconds;
};

enum class PgSupportTimeSide : std::uint8_t { Interior, Left, Right };
enum class PgSupportTimeOrder : std::uint8_t { Value, Rate, Acceleration };

struct PgSupportHistoryLimits {
    PgReferenceLimits arithmetic; ///< Budget for a preparation/evaluation outside chart solves.
    /// Per coefficient solve, at most four per interval. Cancellation uses arithmetic.cancellation.
    PgReferenceLimits chart_verification;
    std::uint64_t maximum_input_knots{}, maximum_coefficient_triples{};
};

/// Exact C1 original-row right-hand sides from complete spatially mapped support samples.
class PgSupportHistory {
public:
    PgSupportHistory(PgSupportHistory&&) noexcept;
    PgSupportHistory& operator=(PgSupportHistory&&) noexcept;
    PgSupportHistory(const PgSupportHistory&) = delete;
    PgSupportHistory& operator=(const PgSupportHistory&) = delete;
    ~PgSupportHistory();

    /// A support absent from the track list keeps the chart's original constant RHS.
    [[nodiscard]] static core::Result<PgSupportHistory>
    prepare(PgSparseChart&&, PgSupportTimeWindow, std::span<const PgSupportTimeTrack>,
            const PgSupportHistoryLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    [[nodiscard]] const PgSparseChart& chart() const& noexcept;
    [[nodiscard]] const PgSparseChart& chart() const&& = delete;
    [[nodiscard]] std::span<const ExactRational> events_seconds() const& noexcept;
    [[nodiscard]] std::span<const ExactRational> events_seconds() const&& = delete;

    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    evaluate(const ExactRational& seconds, PgSupportTimeOrder, PgSupportTimeSide,
             const PgReferenceLimits&, core::Workspace&, core::ResourcePlan&) const noexcept;

private:
    struct Data;
    PgSupportHistory(PgSparseChart&&, core::PlannedBuffer<Data>) noexcept;
    PgSparseChart chart_;
    core::PlannedBuffer<Data> data_;
};

} // namespace uvlas::structure
