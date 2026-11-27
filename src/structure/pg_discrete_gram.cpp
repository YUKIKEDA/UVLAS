#include "../core/resource_plan_worker.hpp"
#include "pg_conditioning_factor.hpp"
#include "pg_gram_worker_plan.hpp"
#include "pg_reference_arithmetic.hpp"

#include <new>
#include <stop_token>
#include <thread>
#include <uvlas/structure/pg_discrete_gram.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;

struct IntegerWork {
    ReferenceArithmetic& math;
    std::uint64_t limit = 0;
    std::uint64_t steps = 0;

    void charge(std::uint64_t count = 1) {
        math.cancel();
        if (count > limit - steps)
            partition_fail(core::FailureCode::ResourceExceeded);
        steps += count;
    }
};

struct Interval {
    RankWide lower;
    RankWide upper;
};

[[nodiscard]] RankWide shift_magnitude(IntegerWork& work, const RankWide& value,
                                       std::uint32_t bits) {
    const auto size = partition_add(rank_bits(value), bits);
    if (size > 8192)
        partition_fail(core::FailureCode::ResourceExceeded);
    work.charge(1 + size / std::numeric_limits<RankWord>::digits);
    return rank_abs(value) << bits;
}

[[nodiscard]] Interval outward_quotient(IntegerWork& work, const RankWide& numerator,
                                        const RankWide& denominator, std::uint32_t bits) {
    const auto scaled = shift_magnitude(work, numerator, bits);
    const auto quotient = conditioning_divide(work, scaled, denominator);
    const auto inexact = rank_product(work, quotient, denominator) != scaled;
    if (numerator < 0)
        return {-quotient - (inexact ? 1 : 0), -quotient};
    return {quotient, quotient + (inexact ? 1 : 0)};
}

[[nodiscard]] RankFraction<RankWide> fraction(IntegerWork& work, const ExactRational& value) {
    const auto n = value.numerator_bits(), d = value.denominator_bits();
    if (n > 8192 || d > 8192)
        partition_fail(core::FailureCode::ResourceExceeded);
    constexpr auto limb_bits = std::numeric_limits<RankWord>::digits;
    work.charge(2 + n / limb_bits + d / limb_bits);
    return rank_fraction<RankWide>(value);
}

[[nodiscard]] Interval rational_bounds(IntegerWork& work, const ExactRational& value,
                                       std::uint32_t bits) {
    work.charge();
    if (value.sign() == 0)
        return {};
    const auto exact = fraction(work, value);
    return outward_quotient(work, exact.numerator, exact.denominator, bits);
}

[[nodiscard]] Interval root_bounds(IntegerWork& work, const ExactRational& g, std::uint32_t bits) {
    const auto exact = fraction(work, g);
    const auto scaled = shift_magnitude(work, exact.numerator, 2 * bits);
    const auto quotient = conditioning_divide(work, scaled, exact.denominator);
    auto root = conditioning_root(work, quotient);
    auto square = rank_product(work, root, root);
    if (square > quotient) {
        --root;
        square = rank_product(work, root, root);
    }
    const bool same = rank_product(work, square, exact.denominator) == scaled;
    return {root, root + (same ? 0 : 1)};
}

[[nodiscard]] RankWide multiply_endpoint(IntegerWork& work, const RankWide& left,
                                         const RankWide& right) {
    using DoubleWord = boost::multiprecision::double_limb_type;
    constexpr auto limb_bits = std::numeric_limits<RankWord>::digits;
    static_assert(std::numeric_limits<DoubleWord>::digits == 2 * limb_bits);
    const auto& a = left.backend();
    const auto& b = right.backend();
    work.charge(partition_mul(a.size(), b.size()));
    const auto size = a.size() + b.size();
    if (size > 8192 / limb_bits)
        partition_fail(core::FailureCode::ResourceExceeded);

    RankWide product;
    auto& target = product.backend();
    target.resize(size, size);
    const std::span<RankWord> output{target.limbs(), size};
    const std::span<const RankWord> lhs{a.limbs(), a.size()}, rhs{b.limbs(), b.size()};
    std::ranges::fill(output, RankWord{0});
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        DoubleWord carry = 0;
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            // A limb product plus an existing limb and carry fits in two limbs.
            const DoubleWord word =
                static_cast<DoubleWord>(lhs[i]) * rhs[j] + output[i + j] + carry;
            output[i + j] = static_cast<RankWord>(word);
            carry = word >> limb_bits;
        }
        output[i + rhs.size()] = static_cast<RankWord>(carry);
    }
    target.normalize();
    target.sign(a.sign() != b.sign());
    return product;
}

/// Multiply an arbitrary interval by a nonnegative interval on the same grids.
[[nodiscard]] Interval with_positive(IntegerWork& work, const Interval& value,
                                     const Interval& positive) {
    return {
        multiply_endpoint(work, value.lower, value.lower < 0 ? positive.upper : positive.lower),
        multiply_endpoint(work, value.upper, value.upper < 0 ? positive.lower : positive.upper)};
}

[[nodiscard]] Interval interval_product(IntegerWork& work, const Interval& left,
                                        const Interval& right) {
    if (right.lower >= 0)
        return with_positive(work, left, right);
    if (left.lower >= 0)
        return with_positive(work, right, left);
    if (left.upper <= 0) {
        const auto& minimum = right.upper <= 0 ? left.upper : left.lower;
        return {multiply_endpoint(work, minimum, right.upper),
                multiply_endpoint(work, left.lower, right.lower)};
    }
    if (right.upper <= 0)
        return {multiply_endpoint(work, left.upper, right.lower),
                multiply_endpoint(work, left.lower, right.lower)};
    const std::array endpoints{multiply_endpoint(work, left.lower, right.lower),
                               multiply_endpoint(work, left.lower, right.upper),
                               multiply_endpoint(work, left.upper, right.lower),
                               multiply_endpoint(work, left.upper, right.upper)};
    const auto bounds = std::ranges::minmax_element(endpoints);
    return {*bounds.min, *bounds.max};
}

[[nodiscard]] Interval rescale(IntegerWork& work, const Interval& value, std::uint32_t drop) {
    if (drop == 0)
        return value;
    const auto lower = rank_abs(value.lower) >> drop;
    const auto upper = rank_abs(value.upper) >> drop;
    const auto lower_remainder = shift_magnitude(work, lower, drop) != rank_abs(value.lower);
    const auto upper_remainder = shift_magnitude(work, upper, drop) != rank_abs(value.upper);
    return {value.lower < 0 ? RankWide{-lower - (lower_remainder ? 1 : 0)} : lower,
            value.upper < 0 ? RankWide{-upper} : RankWide{upper + (upper_remainder ? 1 : 0)}};
}

void check_slot(const Interval& value) {
    if (rank_bits(value.lower) > 1536 || rank_bits(value.upper) > 1536)
        partition_fail(core::FailureCode::ResourceExceeded);
}

[[nodiscard]] Interval radical_bounds(IntegerWork& work, const PgReferenceScalar& value,
                                      const Interval& root, std::uint32_t bits) {
    const auto a = rational_bounds(work, value.rational, bits);
    const auto b = rational_bounds(work, value.radical, bits);
    // These limits also bound the intermediate product before reducing it.
    check_slot(a);
    check_slot(b);
    auto radical = rescale(work, with_positive(work, b, root), bits);
    radical.lower += a.lower;
    radical.upper += a.upper;
    check_slot(radical);
    return radical;
}

void validate_rule(ReferenceArithmetic& math, std::span<const PgDiscreteVolumePoint> points,
                   std::size_t volumes) {
    const auto count = partition_index(volumes);
    core::Index previous = -1;
    for (const auto& point : points) {
        math.visit(1);
        const auto current = point.point.volume;
        if (current < 0 || current >= count || current < previous || current > previous + 1 ||
            point.weight.sign() <= 0)
            partition_fail(core::FailureCode::InvalidInput);
        previous = current;
    }
    if (count == 0 || previous != count - 1)
        partition_fail(core::FailureCode::InvalidInput);
}

[[nodiscard]] std::uint64_t product_count(PgVolumeProductRule input, std::size_t volumes) {
    if (volumes == 0 || input.radial.empty() || input.transverse.empty() || input.depth.empty())
        partition_fail(core::FailureCode::InvalidInput);
    return partition_mul(partition_mul(volumes, input.radial.size()),
                         partition_mul(input.transverse.size(), input.depth.size()));
}

void validate_factor(ReferenceArithmetic& math, std::span<const PgLineQuadratureNode> nodes,
                     const ExactRational& one, bool triangle_jacobian) {
    ExactRational mass, first_moment;
    const ExactRational* previous = nullptr;
    for (const auto& node : nodes) {
        math.visit(1);
        if (node.coordinate.sign() <= 0 || math.compare(node.coordinate, one) >= 0 ||
            node.weight.sign() <= 0 ||
            (previous != nullptr && math.compare(*previous, node.coordinate) >= 0))
            partition_fail(core::FailureCode::InvalidInput);
        previous = &node.coordinate;
        mass = math.add(mass, node.weight);
        if (triangle_jacobian)
            first_moment = math.add(first_moment, math.mul(node.coordinate, node.weight));
    }
    if (mass != one || (triangle_jacobian && math.add(first_moment, first_moment) != one))
        partition_fail(core::FailureCode::InvalidInput);
}

void fill_product(ReferenceArithmetic& math, PgVolumeProductRule input, core::Index volumes,
                  std::span<PgOwnedDiscreteVolumePoint> output) {
    const auto one = math.constant("1");
    std::size_t index = 0;
    for (core::Index volume = 0; volume < volumes; ++volume) {
        for (const auto& radial : input.radial) {
            const auto side = math.sub(one, radial.coordinate);
            const auto mapped_weight = math.mul(radial.weight, side);
            for (const auto& transverse : input.transverse) {
                const auto v = math.mul(side, transverse.coordinate);
                const auto surface_weight = math.mul(mapped_weight, transverse.weight);
                for (const auto& depth : input.depth) {
                    math.visit(1);
                    output[index++] = {volume, math.copy(radial.coordinate), math.copy(v),
                                       math.copy(depth.coordinate),
                                       math.mul(surface_weight, depth.weight)};
                }
            }
        }
    }
}

void accumulate_sample(IntegerWork& work, std::span<Interval> entries,
                       std::span<const Interval> values, const Interval& weight) {
    for (std::size_t j = 0; j < values.size(); ++j) {
        const auto& right = values[j];
        // No rounding occurs on these integer grids. Interval multiplication
        // is associative, so the weighted column can be shared by the row.
        const auto weighted =
            right.lower == 0 && right.upper == 0 ? Interval{} : with_positive(work, right, weight);
        for (std::size_t i = 0; i <= j; ++i) {
            work.charge();
            const auto& left = values[i];
            if ((left.lower == 0 && left.upper == 0) ||
                (weighted.lower == 0 && weighted.upper == 0))
                continue;
            const auto contribution = interval_product(work, left, weighted);
            auto& entry = entries[packed(i, j)];
            work.charge(
                std::max(entry.lower.backend().size(), contribution.lower.backend().size()) +
                std::max(entry.upper.backend().size(), contribution.upper.backend().size()));
            entry.lower += contribution.lower;
            entry.upper += contribution.upper;
        }
    }
}

} // namespace

PgFixedDiscreteGramCertificate::PgFixedDiscreteGramCertificate(
    PgConditioningCertificate&& conditioning,
    core::PlannedBuffer<PgOwnedDiscreteVolumePoint> rule) noexcept
    : conditioning_(std::move(conditioning)), rule_(std::move(rule)) {}

PgDiscreteGram::PgDiscreteGram(PgReferenceBasis&& basis, core::PlannedBuffer<PgGramEntry> entries,
                               std::uint32_t bits,
                               core::PlannedBuffer<PgOwnedDiscreteVolumePoint> rule,
                               PgDiscreteGramInfo info) noexcept
    : basis_(std::move(basis)),
      entries_(std::move(entries)),
      fractional_bits_(bits),
      rule_(std::move(rule)),
      info_(info) {}

core::Result<PgDiscreteGram>
PgDiscreteGram::prepare(PgReferenceBasis&& basis, std::span<const PgDiscreteVolumePoint> points,
                        PgDiscreteGramGrid grid, const PgDiscreteGramLimits& limits,
                        core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return prepare_rule(std::move(basis), points, grid, limits, workspace, plan);
}

core::Result<PgDiscreteGram>
PgDiscreteGram::prepare_product(PgReferenceBasis&& basis, PgVolumeProductRule product,
                                PgDiscreteGramGrid grid, const PgDiscreteGramLimits& limits,
                                core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return prepare_rule(std::move(basis), product, grid, limits, workspace, plan);
}

core::Result<PgDiscreteGram> PgDiscreteGram::prepare_rule(PgReferenceBasis&& basis, RuleInput input,
                                                          PgDiscreteGramGrid grid,
                                                          const PgDiscreteGramLimits& limits,
                                                          core::Workspace& workspace,
                                                          core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgDiscreteGram> {
        try {
            auto math = ReferenceArithmetic{limits.arithmetic, plan};
            math.cancel();
            const auto lease = partition_take(workspace.acquire());
            static_cast<void>(lease);
            const auto bits = grid.working_fractional_bits;
            if (bits > 1535 || grid.output_fractional_bits > bits)
                partition_fail(core::FailureCode::InvalidInput);
            const auto rank = partition_take(basis.size());
            if (rank == 0 || rank > rank_max_columns)
                partition_fail(core::FailureCode::InvalidInput);
            const auto source = partition_take(basis.source());
            const auto volumes = partition_take(source.partition.volumes());
            const auto* points = std::get_if<std::span<const PgDiscreteVolumePoint>>(&input);
            const auto* product = std::get_if<PgVolumeProductRule>(&input);
            const auto size =
                points != nullptr ? points->size() : product_count(*product, volumes.size());
            const auto point_count = partition_index(size);
            if (points != nullptr)
                validate_rule(math, *points, volumes.size());
            const auto entries = packed(0, static_cast<std::size_t>(rank));
            const auto storage =
                partition_add(partition_mul(entries, sizeof(Interval) + sizeof(PgGramEntry)),
                              partition_mul(rank, sizeof(Interval)));
            const auto rule_bytes = partition_mul(size, sizeof(PgOwnedDiscreteVolumePoint));
            if (partition_add(storage, rule_bytes) > plan.remaining_bytes())
                partition_fail(core::FailureCode::ResourceExceeded);
            // The exact rule and both endpoints of every result entry can coexist.
            const auto held_rationals =
                partition_add(partition_mul(4, size), partition_mul(2, entries));
            const auto live_rationals = partition_add(basis_live_rationals + 16, held_rationals);
            math.admit(live_rationals);
            if (product != nullptr) {
                const auto one = math.constant("1");
                validate_factor(math, product->radial, one, true);
                validate_factor(math, product->transverse, one, false);
                validate_factor(math, product->depth, one, false);
            }
            auto rule = partition_take(
                plan.allocate<PgOwnedDiscreteVolumePoint>(point_count, partition_operation));
            if (points != nullptr) {
                for (std::size_t i = 0; i < points->size(); ++i) {
                    math.visit(1);
                    const auto& p = (*points)[i];
                    rule.values()[i] = {p.point.volume, math.copy(p.point.u), math.copy(p.point.v),
                                        math.copy(p.point.depth), math.copy(p.weight)};
                }
            } else {
                fill_product(math, *product, partition_index(volumes.size()), rule.values());
            }
            auto samples =
                partition_take(plan.allocate<Interval>(partition_index(rank), partition_operation));
            auto sums = partition_take(
                plan.allocate<Interval>(partition_index(entries), partition_operation));
            auto result = partition_take(
                plan.allocate<PgGramEntry>(partition_index(entries), partition_operation));
            using Numerators = std::array<ExactRational, 3>;
            core::PlannedBuffer<Numerators> cached;
            std::size_t depth_group = 1;
            using PowerReuse = PgRankMonomialPowers::Reuse;
            std::optional<PowerReuse> power_reuse;
            if (product != nullptr && product->depth.size() > 1) {
                auto cache_arithmetic = ReferenceArithmetic{limits.arithmetic, plan};
                const auto cache_rationals = partition_add(live_rationals, partition_mul(rank, 3));
                try {
                    cache_arithmetic.admit(cache_rationals);
                    cached = partition_take(
                        plan.allocate<Numerators>(partition_index(rank), partition_operation));
                    math.leaf_limits = cache_arithmetic.leaf_limits;
                    depth_group = product->depth.size();
                } catch (const PartitionAbort& unavailable) {
                    if (unavailable.failure.code != core::FailureCode::ResourceExceeded)
                        throw;
                }
                if (depth_group > 1) {
                    const auto axis_slots =
                        partition_add(cache_rationals, PgRankMonomialPowers::retained_rationals);
                    const auto mixed_slots =
                        partition_add(axis_slots, PgRankMonomialPowers::retained_mixed);
                    const std::array preferences{
                        std::pair{
                            PowerReuse::CommonDenominator,
                            partition_add(mixed_slots, PgRankMonomialPowers::retained_common)},
                        std::pair{PowerReuse::Mixed, mixed_slots},
                        std::pair{PowerReuse::Axis, axis_slots}};
                    for (const auto& [reuse, slots] : preferences) {
                        cache_arithmetic.leaf_limits = limits.arithmetic.algebra;
                        try {
                            cache_arithmetic.admit(slots);
                            math.leaf_limits = cache_arithmetic.leaf_limits;
                            power_reuse = reuse;
                            break;
                        } catch (const PartitionAbort& unavailable) {
                            if (unavailable.failure.code != core::FailureCode::ResourceExceeded)
                                throw;
                        }
                    }
                }
            }
            IntegerWork work{math, limits.maximum_integer_steps};
            const auto& length = basis.length_metres();
            const auto length_cube = math.mul(math.square(length), length);
            const auto two = math.constant("2");
            const auto two_length = math.mul(two, length);
            const auto evaluate_range = [&](std::span<const PgOwnedDiscreteVolumePoint> range,
                                            ReferenceArithmetic& point_math,
                                            IntegerWork& point_work,
                                            std::span<Interval> point_samples,
                                            std::span<Numerators> point_cache,
                                            std::span<Interval> output) {
                std::size_t depth_slot = 0;
                for (const auto& p : range) {
                    point_math.visit(1);
                    const auto located = basis.locate(
                        point_math, PgReferenceVolumePoint{p.volume, p.u, p.v, p.depth});
                    const auto& geometry = located.reference;
                    const auto& jacobian =
                        std::get<PgReferenceVolumeMeasure>(geometry.measure).jacobian;
                    const auto root = root_bounds(point_work, geometry.g, bits);
                    check_slot(root);
                    ReferenceRadicals field{point_math, geometry.g};
                    const auto scale = point_math.div(p.weight, length_cube);
                    auto measure =
                        radical_bounds(point_work, field.scale(jacobian, scale), root, bits);
                    // The certified Jacobian and the supplied scalar weight are positive.
                    measure.lower = std::max(RankWide{}, measure.lower);
                    const auto z = point_math.div(
                        point_math.mul(geometry.zeta, geometry.full_thickness_metres), two_length);
                    const auto half_z_squared = point_math.div(point_math.square(z), two);
                    const auto director_radical = point_math.div(z, point_math.square(geometry.g));
                    const auto evaluate = [&](const Numerators& terms) {
                        auto a = point_math.add(terms[0], point_math.mul(half_z_squared, terms[2]));
                        auto b = point_math.mul(director_radical, terms[1]);
                        return radical_bounds(point_work, {std::move(a), std::move(b)}, root, bits);
                    };
                    std::optional<PgRankMonomialPowers> monomials;
                    if (depth_slot == 0 && power_reuse &&
                        PgRankMonomialPowers::fits(located.u, located.v, point_math.leaf_limits))
                        monomials.emplace(located.u, located.v, *power_reuse);
                    auto values = point_samples;
                    for (std::size_t i = 0; i < values.size(); ++i) {
                        point_math.visit(1);
                        const auto numerator_values = [&]() -> Numerators {
                            auto jet = basis.numerator_jet(
                                point_math, {partition_index(i)}, located,
                                /*derivatives=*/false, monomials ? &*monomials : nullptr);
                            return {std::move(jet[0][0]), std::move(jet[1][0]),
                                    std::move(jet[2][0])};
                        };
                        if (point_cache.size() == 0) {
                            values[i] = evaluate(numerator_values());
                            continue;
                        }
                        auto& terms = point_cache[i];
                        if (depth_slot == 0)
                            terms = numerator_values();
                        values[i] = evaluate(terms);
                    }
                    accumulate_sample(point_work, output, values, measure);
                    // fill_product emits all depths consecutively for one volume and (u, v).
                    if (++depth_slot == depth_group)
                        depth_slot = 0;
                }
            };
            using namespace detail::parallel_gram;
            constexpr std::size_t count = 4;
            const auto concurrency =
                std::min<std::size_t>(count, std::thread::hardware_concurrency());
            bool merged = false;
            // The recorded crossover control contains 192 points per fixed task.
            constexpr std::size_t minimum_points_per_task = 192;
            if (concurrency > 1 && depth_group > 1 &&
                rule.values().size() >= count * minimum_points_per_task &&
                rule.values().size() / depth_group >= count &&
                rule.values().size() % depth_group == 0) {
                const auto suggest = [&]() -> std::optional<Forecast> {
                    try {
                        return forecast<count, Interval, Numerators>(
                            basis, source.partition, rule.values(), depth_group, entries, rank,
                            held_rationals, math, limits, work);
                    } catch (const PartitionAbort& failure) {
                        if (failure.failure.code != core::FailureCode::ResourceExceeded)
                            throw;
                        return std::nullopt;
                    }
                };
                const auto admitted_jobs = suggest();
                if (admitted_jobs) {
                    core::detail::WorkerPayloadReservation reservation{plan};
                    if (!reservation.begin(admitted_jobs->sum.bytes))
                        partition_fail(core::FailureCode::ResourceExceeded);
                    struct WorkerResult {
                        core::PlannedBuffer<Interval> sums;
                        Budget used;
                        std::optional<core::Failure> failure;
                    };
                    std::array<WorkerResult, count> children;
                    std::stop_source stop_source;
                    const std::stop_callback caller_stop{limits.arithmetic.cancellation,
                                                         [&] { stop_source.request_stop(); }};
                    const auto worker = [&](std::size_t i) {
                        const auto& job = admitted_jobs->jobs[i];
                        auto& child = children[i];
                        core::ResourcePlan local_plan{job.cap.bytes};
                        core::Workspace local_workspace{core::PlannedBuffer<std::byte>{}};
                        auto local_limits = limits.arithmetic;
                        local_limits.maximum_scalar_operations = job.cap.scalar;
                        local_limits.maximum_partition_visits = job.cap.visits;
                        local_limits.algebra.external_memory =
                            ExactEstimateAdmission{admitted_jobs->external_per_worker};
                        local_limits.cancellation = stop_source.get_token();
                        ReferenceArithmetic local_math{local_limits, local_plan};
                        IntegerWork local_work{local_math, job.cap.integer};
                        try {
                            local_math.cancel();
                            const auto local_lease = partition_take(local_workspace.acquire());
                            local_math.admit(admitted_jobs->worker_live);
                            auto child_sums = partition_take(local_plan.allocate<Interval>(
                                partition_index(entries), partition_operation));
                            auto child_samples = partition_take(local_plan.allocate<Interval>(
                                partition_index(rank), partition_operation));
                            auto child_cache = partition_take(local_plan.allocate<Numerators>(
                                partition_index(rank), partition_operation));
                            evaluate_range(rule.values().subspan(job.first, job.points), local_math,
                                           local_work, child_samples.values(), child_cache.values(),
                                           child_sums.values());
                            local_math.cancel();
                            child.sums = std::move(child_sums);
                        } catch (const PartitionAbort& failure) {
                            child.failure = failure.failure;
                        } catch (const std::overflow_error&) {
                            child.failure = core::Failure{core::FailureCode::ResourceExceeded,
                                                          partition_operation};
                        } catch (const std::bad_alloc&) {
                            child.failure = core::Failure{core::FailureCode::ResourceExceeded,
                                                          partition_operation};
                        } catch (...) {
                            child.failure =
                                core::Failure{core::FailureCode::InvalidInput, partition_operation};
                        }
                        child.used = {local_plan.allocated_bytes(), local_math.scalar_operations,
                                      local_math.partition_visits, local_work.steps};
                    };
                    std::optional<core::Failure> launch_failure;
                    {
                        std::array<std::jthread, count> threads;
                        try {
                            for (std::size_t i = 0; i < concurrency; ++i) {
                                threads[i] = std::jthread([&, i] {
                                    for (std::size_t id = i; id < count; id += concurrency)
                                        worker(id);
                                });
                            }
                        } catch (...) {
                            stop_source.request_stop();
                            launch_failure = core::Failure{core::FailureCode::ResourceExceeded,
                                                           partition_operation};
                        }
                    }
                    Budget receipt{};
                    for (const auto& child : children)
                        receipt.add(child.used);
                    if (!reservation.settle(receipt.bytes))
                        partition_fail(core::FailureCode::ResourceExceeded);
                    math.scalar_operations = partition_add(math.scalar_operations, receipt.scalar);
                    math.partition_visits = partition_add(math.partition_visits, receipt.visits);
                    work.steps = partition_add(work.steps, receipt.integer);
                    if (launch_failure)
                        throw PartitionAbort{*launch_failure};
                    math.cancel();
                    for (const auto& child : children)
                        if (child.failure)
                            throw PartitionAbort{*child.failure};
                    for (const auto& child : children)
                        for (std::size_t i = 0; i < entries; ++i) {
                            auto& dst = sums.values()[i];
                            const auto& src = child.sums.values()[i];
                            work.charge(
                                std::max(dst.lower.backend().size(), src.lower.backend().size()) +
                                std::max(dst.upper.backend().size(), src.upper.backend().size()));
                            dst.lower += src.lower;
                            dst.upper += src.upper;
                        }
                    merged = true;
                }
            }
            if (!merged) {
                evaluate_range(rule.values(), math, work, samples.values(), cached.values(),
                               sums.values());
            }
            const auto drop = 3 * bits - grid.output_fractional_bits;
            for (std::size_t i = 0; i < entries; ++i) {
                work.charge();
                const auto output = rescale(work, sums.values()[i], drop);
                check_slot(output);
                math.charge(2);
                result.values()[i] = {
                    rank_materialize_dyadic(output.lower, grid.output_fractional_bits,
                                            math.leaf_limits, plan),
                    rank_materialize_dyadic(output.upper, grid.output_fractional_bits,
                                            math.leaf_limits, plan)};
            }
            math.cancel();
            return PgDiscreteGram{std::move(basis),
                                  std::move(result),
                                  grid.output_fractional_bits,
                                  std::move(rule),
                                  {math.info(), work.steps, static_cast<std::size_t>(point_count)}};
        } catch (const std::overflow_error&) {
            partition_fail(core::FailureCode::ResourceExceeded);
        }
    });
}

core::Result<PgDiscreteGramEnclosure> PgDiscreteGram::enclosure() const& noexcept {
    const auto count = basis_.size();
    if (!count)
        return std::unexpected(count.error());
    return PgDiscreteGramEnclosure{entries_.values(), fractional_bits_};
}

core::Result<PgFixedDiscreteGramCertificate>
PgDiscreteGram::generate_and_certify(PgConditioningGeneration generation,
                                     const PgConditioningLimits& limits, core::Workspace& workspace,
                                     core::ResourcePlan& plan) && noexcept {
    const auto input = enclosure();
    if (!input)
        return std::unexpected(input.error());
    auto accepted = PgConditioningCertificate::generate_and_certify(
        std::move(basis_), *input, generation, limits, workspace, plan);
    if (!accepted)
        return std::unexpected(accepted.error());
    entries_ = {};
    return PgFixedDiscreteGramCertificate{std::move(*accepted), std::move(rule_)};
}

} // namespace uvlas::structure
