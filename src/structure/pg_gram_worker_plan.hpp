#pragma once

#include "pg_rank_storage.hpp"
#include "pg_reference_arithmetic.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <uvlas/structure/pg_discrete_gram.hpp>
#include <variant>

namespace uvlas::structure::detail::parallel_gram {
using Count = std::uint64_t;
struct Budget {
    Count bytes{}, scalar{}, visits{}, integer{};
    void add(const Budget& b) {
        bytes = partition_add(bytes, b.bytes);
        scalar = partition_add(scalar, b.scalar);
        visits = partition_add(visits, b.visits);
        integer = partition_add(integer, b.integer);
    }
};
struct Job {
    std::size_t first{}, points{};
    Budget cap;
};
struct Forecast {
    std::array<Job, 4> jobs;
    Budget sum;
    Count external_per_worker{}, worker_live{};
};

template <std::size_t workers, class Interval, class Numerators, class Work>
std::optional<Forecast>
forecast(const PgReferenceBasis& basis, const PgIntegrationPartition& partition,
         std::span<const PgOwnedDiscreteVolumePoint> rule, std::size_t depths, Count entries,
         Count rank, Count held_rationals, ReferenceArithmetic& math,
         const PgDiscreteGramLimits& limits, const Work& used) {
    static_assert(workers == 4);
    constexpr Count native = std::numeric_limits<RankWord>::digits;
    constexpr Count M = 8192 / native, W = (1536 + native - 1) / native;
    const auto accum_words = (Count{3} * 1536 + partition_log(rule.size()) + native - 1) / native;
    // Fixed integer work is charged per operation/limb product. Square-root
    // The initial relative excess is at most one. The real Newton step
    // squares it and divides it by at least two; flooring cannot increase it.
    // Once the absolute excess is below one, the next decreasing step
    // reaches the integer root. Two extra iterations cover the floor and exit.
    constexpr Count root_iterations = std::bit_width(std::uint64_t{8192}) + 2;
    constexpr Count root_work = root_iterations * (1 + M * M) + 6 * M * M + 4 * M + 5;
    constexpr Count radical_work = 20 + 10 * M + 6 * M * M;
    const auto point_integer =
        partition_add(partition_add(root_work, partition_mul(rank + 1, radical_work)),
                      partition_add(partition_mul(entries, 1 + 8 * W * W + 2 * accum_words),
                                    partition_mul(rank, 2 * W * W)));
    const auto shared_external = math.external_integer_payload(partition_mul(2, held_rationals));
    const auto* external =
        std::get_if<ExactEstimateAdmission>(&limits.arithmetic.algebra.external_memory);
    if (!external || shared_external > external->bytes)
        return std::nullopt;
    Forecast result{};
    result.external_per_worker = (external->bytes - shared_external) / workers;
    constexpr auto power_cache =
        PgRankMonomialPowers::retained_rationals + PgRankMonomialPowers::retained_mixed;
    result.worker_live = partition_add(
        basis_live_rationals + 16,
        partition_add(partition_mul(rank, 3), power_cache + PgRankMonomialPowers::retained_common));
    auto preflight_limits = limits.arithmetic;
    preflight_limits.algebra.external_memory = ExactEstimateAdmission{result.external_per_worker};
    ReferenceArithmetic child_preflight{preflight_limits, math.plan};
    child_preflight.admit(result.worker_live);
    // Standard rational leaf operations hold at most two borrowed rationals.
    // Provide their external allowance and a spare slot above the guarded
    // maximum already checked by admit(). No worker can enter solely because
    // the shared rule was counted as though it lived in just one child.
    constexpr Count integer_limb = IntegerBackend::limb_bits;
    const auto allowed_limbs = limits.arithmetic.algebra.maximum_work_bits / integer_limb;
    partition_take(exact_rational_admission(child_preflight.leaf_limits,
                                            (allowed_limbs - 2) * integer_limb,
                                            child_preflight.external_integer_payload(8)));
    const auto* selected = partition_take(basis.rank());
    const auto rows = partition_take(selected->rows());
    const auto groups = rule.size() / depths;
    const auto universal_row = rows.size();
    // Confirm using only extents that scanning can fit alongside every
    // possible serial visit. If the other estimates fail, the same visit
    // budget can still admit the unchanged serial kernel.
    const auto worst_search = 3 * (static_cast<Count>(std::bit_width(universal_row)) + 1);
    const auto serial_visits =
        partition_add(partition_mul(rule.size(), rank + 2),
                      partition_mul(groups, partition_mul(rank, 1 + worst_search + universal_row)));
    const auto planner_visits =
        partition_add(8, partition_add(partition_mul(rank, partition_add(1, universal_row)),
                                       partition_mul(groups, 3)));
    if (partition_add(serial_visits, planner_visits) >
        limits.arithmetic.maximum_partition_visits - math.partition_visits)
        return std::nullopt;
    math.visit(8);
    const auto graph = partition_take(partition_take(basis.source()).dictionary.graph());
    const auto axis = static_cast<std::size_t>(2 * graph.n);
    std::array<Count, 64> term_scalar{}, term_visits{};
    const auto cell_slot = [axis](PgCellIndex c) {
        return static_cast<std::size_t>(c.first) * axis + static_cast<std::size_t>(c.second);
    };
    Count all_column_searches = 0;
    for (std::size_t i = 0; i < rank; ++i) {
        const auto indices = partition_take(selected->basis_indices({partition_index(i)}));
        const auto search_width = static_cast<Count>(std::bit_width(indices.size()));
        all_column_searches = partition_add(all_column_searches, 2 + 3 * (search_width + 1));
        math.visit(partition_add(1, indices.size()));
        for (auto index : indices) {
            math.cancel();
            const auto& row = rows[index];
            const auto slot = cell_slot(row.cell);
            term_scalar[slot] = partition_add(
                term_scalar[slot], 3 + static_cast<Count>(row.powers.first + row.powers.second));
            term_visits[slot] = partition_add(term_visits[slot], 1);
        }
    }
    const auto volumes = partition_take(partition.volumes());
    const auto triangles = partition_take(partition.triangles());
    const auto pieces = partition_take(partition.pieces());
    // <= 10 charged operations per height term per jet order. The straight
    // line reference geometry/determinant, domain blend and Gram measure
    // use below 512 excluding the height polynomial and per-column calls.
    const Count point_scalar = 512 + 150 * graph.height.terms().size() + 3 * rank;
    const auto scratch = partition_add(partition_mul(entries, sizeof(Interval)),
                                       partition_mul(rank, sizeof(Interval) + sizeof(Numerators)));
    for (std::size_t task = 0; task < workers; ++task) {
        auto& job = result.jobs[task];
        const auto gbegin = (groups / workers) * task + std::min(groups % workers, task);
        const auto gend = (groups / workers) * (task + 1) + std::min(groups % workers, task + 1);
        job.first = gbegin * depths;
        job.points = (gend - gbegin) * depths;
        auto& cap = job.cap;
        cap.scalar = partition_mul(job.points, point_scalar);
        cap.visits = partition_mul(job.points, 2);
        cap.integer = partition_mul(job.points, point_integer);
        for (std::size_t i = job.first; i < job.first + job.points; i += depths) {
            math.cancel();
            const auto& volume = volumes[static_cast<std::size_t>(rule[i].volume)];
            const auto& triangle = triangles[static_cast<std::size_t>(volume.triangle)];
            const auto cell = cell_slot(pieces[static_cast<std::size_t>(triangle.piece)].cell);
            cap.scalar = partition_add(cap.scalar, partition_add(power_cache, term_scalar[cell]));
            cap.visits = partition_add(
                cap.visits, partition_add(partition_add(all_column_searches, term_visits[cell]),
                                          partition_mul(depths - 1, rank)));
            math.visit(3);
        }
        // Each ReferenceArithmetic scalar charge creates at most one owned
        // RationalStorage; common_sum charges >= its single materialization.
        cap.bytes = partition_add(scratch, partition_mul(cap.scalar, sizeof(RationalStorage)));
        result.sum.add(cap);
    }
    const Count finalize_integer =
        partition_mul(entries, partition_add(2 * workers * accum_words, 3 + 2 * M));
    auto all = result.sum;
    all.add(
        {partition_mul(2 * entries, sizeof(RationalStorage)), 2 * entries, 0, finalize_integer});
    const auto avail_visits = limits.arithmetic.maximum_partition_visits - math.partition_visits;
    const auto avail_scalar = limits.arithmetic.maximum_scalar_operations - math.scalar_operations;
    const auto avail_int = limits.maximum_integer_steps - used.steps;
    if (all.bytes > math.plan.remaining_bytes() || all.scalar > avail_scalar ||
        all.visits > avail_visits || all.integer > avail_int)
        return std::nullopt;
    return result;
}
} // namespace uvlas::structure::detail::parallel_gram
