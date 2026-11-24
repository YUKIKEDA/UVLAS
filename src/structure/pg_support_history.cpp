#include "pg_reference_arithmetic.hpp"

#include <algorithm>
#include <uvlas/structure/pg_support_history.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;
using Position = PgCartesianExactPosition;
using Cubic = std::array<ExactRational, 4>;

template <class T>
auto allocate(ReferenceArithmetic& math, std::uint64_t n) {
    return partition_take(math.plan.allocate<T>(partition_index(n), partition_operation));
}

struct Constants {
    ExactRational two, three, six;
};

Cubic coefficients(ReferenceArithmetic& math, const Constants& k, const ExactRational& h,
                   const ExactRational& x0, const ExactRational& x1, const ExactRational& v0,
                   const ExactRational& v1) {
    const auto delta = math.sub(x1, x0);
    auto a0 = math.copy(x0);
    auto a1 = math.mul(h, v0);
    auto a2 = math.sub(math.mul(k.three, delta), math.mul(h, math.add(math.mul(k.two, v0), v1)));
    auto a3 = math.sub(math.mul(h, math.add(v0, v1)), math.mul(k.two, delta));
    return {std::move(a0), std::move(a1), std::move(a2), std::move(a3)};
}

std::array<ExactRational, 2> value_and_rate(ReferenceArithmetic& math, const Constants& k,
                                            const ExactRational& u, const ExactRational& h,
                                            const Cubic& a) {
    auto value =
        math.add(a[0], math.mul(u, math.add(a[1], math.mul(u, math.add(a[2], math.mul(u, a[3]))))));
    auto rate =
        math.div(math.add(a[1], math.mul(u, math.add(math.mul(k.two, a[2]),
                                                     math.mul(u, math.mul(k.three, a[3]))))),
                 h);
    return {std::move(value), std::move(rate)};
}

template <typename At>
std::size_t lower_index(ReferenceArithmetic& math, std::size_t size, At at,
                        const ExactRational& x) {
    std::size_t first = 0;
    while (size != 0) {
        math.visit(1);
        const auto half = size / 2;
        const auto mid = first + half;
        if (math.compare(at(mid), x) < 0) {
            first = mid + 1;
            size -= half + 1;
        } else {
            size = half;
        }
    }
    return first;
}

std::size_t interval_begin(ReferenceArithmetic& math, std::span<const PgSupportTimeKnot> knots,
                           const ExactRational& x) {
    const auto i = lower_index(
        math, knots.size(), [&](std::size_t j) -> const auto& { return knots[j].seconds; }, x);
    if (i >= knots.size())
        partition_fail(core::FailureCode::InvalidInput);
    if (i + 1 < knots.size() && math.compare(x, knots[i].seconds) == 0)
        return i;
    if (i == 0)
        partition_fail(core::FailureCode::InvalidInput);
    return i - 1;
}

ExactRational at_time(ReferenceArithmetic& math, const Constants& k, const ExactRational& u,
                      const ExactRational& h, PgSupportTimeOrder order, const ExactRational& a0,
                      const ExactRational& a1, const ExactRational& a2, const ExactRational& a3) {
    switch (order) {
    case PgSupportTimeOrder::Value:
        return math.add(a0, math.mul(u, math.add(a1, math.mul(u, math.add(a2, math.mul(u, a3))))));
    case PgSupportTimeOrder::Rate:
        return math.div(math.add(a1, math.mul(u, math.add(math.mul(k.two, a2),
                                                          math.mul(u, math.mul(k.three, a3))))),
                        h);
    case PgSupportTimeOrder::Acceleration:
        return math.div(math.add(math.mul(k.two, a2), math.mul(u, math.mul(k.six, a3))),
                        math.square(h));
    }
    partition_fail(core::FailureCode::InvalidInput);
}

std::size_t selected_interval(ReferenceArithmetic& math, std::span<const ExactRational> events,
                              const ExactRational& t, PgSupportTimeSide side) {
    const auto i = lower_index(
        math, events.size(), [&](std::size_t j) -> const auto& { return events[j]; }, t);
    if (i == events.size())
        partition_fail(core::FailureCode::InvalidInput);
    if (math.compare(t, events[i]) == 0) {
        if (i != 0 && side == PgSupportTimeSide::Left)
            return i - 1;
        if (i + 1 < events.size() && side == PgSupportTimeSide::Right)
            return i;
    } else if (i != 0 && side == PgSupportTimeSide::Interior) {
        return i - 1;
    }
    partition_fail(core::FailureCode::InvalidInput);
}
} // namespace

struct PgSupportHistory::Data {
    core::PlannedBuffer<ExactRational> events;
    core::PlannedBuffer<Position> rhs;
    std::size_t event_count{}, rows{};

    [[nodiscard]] static std::uint64_t live(std::size_t capacity, std::size_t count) {
        return partition_add(capacity, partition_mul(3, count));
    }

    [[nodiscard]] std::span<const ExactRational> time() const noexcept {
        return events.values().first(event_count);
    }

    [[nodiscard]] std::size_t offset(std::size_t interval, std::size_t power) const noexcept {
        return (4 * interval + power) * rows;
    }
};

PgSupportHistory::PgSupportHistory(PgSparseChart&& chart, core::PlannedBuffer<Data> data) noexcept
    : chart_(std::move(chart)), data_(std::move(data)) {}
PgSupportHistory::PgSupportHistory(PgSupportHistory&&) noexcept = default;
PgSupportHistory& PgSupportHistory::operator=(PgSupportHistory&&) noexcept = default;
PgSupportHistory::~PgSupportHistory() = default;

core::Result<PgSupportHistory> PgSupportHistory::prepare(PgSparseChart&& chart,
                                                         PgSupportTimeWindow window,
                                                         std::span<const PgSupportTimeTrack> tracks,
                                                         const PgSupportHistoryLimits& limits,
                                                         core::Workspace& workspace,
                                                         core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgSupportHistory> {
        ReferenceArithmetic math{limits.arithmetic, plan};
        math.cancel();
        const auto identity = partition_take(chart.identity());
        const auto sources = chart.original().sources();
        const auto original_rows = chart.original().original_rows();
        if (tracks.size() > sources.size())
            partition_fail(core::FailureCode::InvalidInput);

        core::PlannedBuffer<Data> owner;
        {
            const auto lease = partition_take(workspace.acquire());
            static_cast<void>(lease);
            std::uint64_t knot_capacity = 2;
            for (const auto& track : tracks) {
                math.visit(1);
                if (track.knots.size() < 2)
                    partition_fail(core::FailureCode::InvalidInput);
                knot_capacity = partition_add(knot_capacity, track.knots.size());
            }
            if (knot_capacity - 2 > limits.maximum_input_knots)
                partition_fail(core::FailureCode::ResourceExceeded);
            math.admit(partition_add(reference_live_rationals, knot_capacity));
            if (math.compare(window.begin_seconds, window.end_seconds) >= 0)
                partition_fail(core::FailureCode::InvalidInput);

            owner = allocate<Data>(math, 1);
            auto& data = owner.values().front();
            data.rows = original_rows.size();
            data.events = allocate<ExactRational>(math, knot_capacity);
            auto events = data.events.values();
            events[0] = math.copy(window.begin_seconds);
            events[1] = math.copy(window.end_seconds);
            std::size_t events_used = 2;
            auto by_source = allocate<const PgSupportTimeTrack*>(math, sources.size());
            for (const auto& track : tracks) {
                const auto s = std::ranges::find_if(sources, [&](const auto& source) {
                    math.visit(1);
                    return source.id == track.id;
                });
                if (s == sources.end() || s->kind != PgTraceKind::Support)
                    partition_fail(core::FailureCode::InvalidInput);
                auto& mapped = by_source.values()[static_cast<std::size_t>(s - sources.begin())];
                if (mapped != nullptr)
                    partition_fail(core::FailureCode::InvalidInput);
                mapped = &track;
                const auto knots = track.knots;
                if (math.compare(knots.front().seconds, window.begin_seconds) > 0 ||
                    math.compare(knots.back().seconds, window.end_seconds) < 0)
                    partition_fail(core::FailureCode::InvalidInput);
                for (std::size_t i = 0; i < knots.size(); ++i) {
                    math.visit(1);
                    const auto& x = knots[i];
                    if (x.value.size() != static_cast<std::size_t>(s->rows) ||
                        x.rate.size() != x.value.size() ||
                        (i != 0 && math.compare(knots[i - 1].seconds, x.seconds) >= 0))
                        partition_fail(core::FailureCode::InvalidInput);
                    if (math.compare(x.seconds, window.begin_seconds) > 0 &&
                        math.compare(x.seconds, window.end_seconds) < 0)
                        events[events_used++] = math.copy(x.seconds);
                }
            }
            const auto active = events.first(events_used);
            std::ranges::sort(active, [&](const auto& a, const auto& b) {
                math.visit(1);
                return math.compare(a, b) < 0;
            });
            const auto duplicate = std::ranges::unique(active, [&](const auto& a, const auto& b) {
                math.visit(1);
                return math.compare(a, b) == 0;
            });
            for (auto& extra : duplicate) {
                math.visit(1);
                extra = {};
            }
            data.event_count = static_cast<std::size_t>(duplicate.begin() - active.begin());
            const auto count = partition_mul(partition_mul(data.event_count - 1, 4), data.rows);
            if (count > limits.maximum_coefficient_triples)
                partition_fail(core::FailureCode::ResourceExceeded);
            // The combined allowance already includes events. Replace the leaf envelope
            // from the original limits while preserving the whole call's work counters.
            ReferenceArithmetic combined{limits.arithmetic, plan};
            combined.admit(
                partition_add(reference_live_rationals, Data::live(events.size(), count)));
            math.leaf_limits = combined.leaf_limits;
            data.rhs = allocate<Position>(math, count);
            const Constants k{math.constant("2"), math.constant("3"), math.constant("6")};

            for (std::size_t interval = 0; interval + 1 < data.event_count; ++interval) {
                math.visit(1);
                const auto& start = events[interval];
                const auto& end = events[interval + 1];
                const auto h = math.sub(end, start);
                auto rhs = data.rhs.values();
                for (std::size_t source = 0; source < sources.size(); ++source) {
                    math.visit(1);
                    const auto& s = sources[source];
                    const auto first = static_cast<std::size_t>(s.first_row);
                    const auto n = static_cast<std::size_t>(s.rows);
                    const auto* track = by_source.values()[source];
                    if (track == nullptr) {
                        for (std::size_t row = first; row < first + n; ++row) {
                            math.visit(1);
                            for (std::size_t c = 0; c < 3; ++c)
                                rhs[data.offset(interval, 0) + row][c] =
                                    math.copy(original_rows[row].prescribed[c]);
                        }
                        continue;
                    }
                    const auto j = interval_begin(math, track->knots, start);
                    const auto& a = track->knots[j];
                    const auto& b = track->knots[j + 1];
                    const auto width = math.sub(b.seconds, a.seconds);
                    const auto u0 = math.div(math.sub(start, a.seconds), width);
                    const auto u1 = math.div(math.sub(end, a.seconds), width);
                    for (std::size_t row = 0; row < n; ++row) {
                        math.visit(1);
                        for (std::size_t c = 0; c < 3; ++c) {
                            const auto base =
                                coefficients(math, k, width, a.value[row][c], b.value[row][c],
                                             a.rate[row][c], b.rate[row][c]);
                            const auto left = value_and_rate(math, k, u0, width, base);
                            const auto right = value_and_rate(math, k, u1, width, base);
                            auto piece =
                                coefficients(math, k, h, left[0], right[0], left[1], right[1]);
                            for (std::size_t power = 0; power < piece.size(); ++power)
                                rhs[data.offset(interval, power) + first + row][c] =
                                    std::move(piece[power]);
                        }
                    }
                }
            }
        }
        const auto& data = owner.values().front();
        auto checking = limits.chart_verification;
        checking.cancellation = limits.arithmetic.cancellation;
        math.reserve_nested_rationals(checking.algebra,
                                      Data::live(data.events.size(), data.rhs.size()));
        const auto independent = allocate<Position>(math, chart.free_columns().size());
        for (std::size_t i = 0; i + 1 < data.event_count; ++i) {
            for (std::size_t k = 0; k < 4; ++k) {
                math.visit(1);
                const auto rhs = data.rhs.values().subspan(data.offset(i, k), data.rows);
                static_cast<void>(partition_take(chart.lift_rhs(
                    {identity, independent.values()}, {identity, rhs}, checking, workspace, plan)));
            }
        }
        math.cancel();
        return PgSupportHistory{std::move(chart), std::move(owner)};
    });
}

const PgSparseChart& PgSupportHistory::chart() const& noexcept {
    return chart_;
}

std::span<const ExactRational> PgSupportHistory::events_seconds() const& noexcept {
    return data_.size() == 0 ? std::span<const ExactRational>{} : data_.values().front().time();
}

core::Result<core::PlannedBuffer<Position>>
PgSupportHistory::evaluate(const ExactRational& seconds, PgSupportTimeOrder order,
                           PgSupportTimeSide side, const PgReferenceLimits& limits,
                           core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<Position>> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        if (data_.size() == 0)
            partition_fail(core::FailureCode::StaleGeneration);
        if (order != PgSupportTimeOrder::Value && order != PgSupportTimeOrder::Rate &&
            order != PgSupportTimeOrder::Acceleration)
            partition_fail(core::FailureCode::InvalidInput);
        const auto& data = data_.values().front();
        math.admit(partition_add(reference_live_rationals, partition_mul(3, data.rows)));
        const auto time = data.time();
        const auto i = selected_interval(math, time, seconds, side);
        const auto h = math.sub(time[i + 1], time[i]);
        const auto u = math.div(math.sub(seconds, time[i]), h);
        const Constants k{math.constant("2"), math.constant("3"), math.constant("6")};
        const auto a = data.rhs.values();
        auto result = allocate<Position>(math, data.rows);
        for (std::size_t row = 0; row < data.rows; ++row) {
            math.visit(1);
            for (std::size_t c = 0; c < 3; ++c)
                result.values()[row][c] =
                    at_time(math, k, u, h, order, a[data.offset(i, 0) + row][c],
                            a[data.offset(i, 1) + row][c], a[data.offset(i, 2) + row][c],
                            a[data.offset(i, 3) + row][c]);
        }
        math.cancel();
        return result;
    });
}

} // namespace uvlas::structure
