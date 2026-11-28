#include "pg_composed_basis.hpp"
#include "pg_conditioned_volume_row.hpp"
#include "pg_conditioning_factor.hpp"
#include "pg_reference_arithmetic.hpp"
#include "pg_reference_hash.hpp"

#include <uvlas/structure/pg_conditioning.hpp>

namespace uvlas::structure::detail {
struct PgConditioningState {
    PgConditioningInfo info{};
    std::uint32_t transform_scale = 0;
    std::uint32_t congruence_scale = 0;
    RankWide bound{};
    core::PlannedBuffer<RankCell> transform;
};
} // namespace uvlas::structure::detail

namespace uvlas::structure {
namespace {
using namespace detail;
constexpr auto operation = core::Operation::CertifyError;
constexpr std::uint32_t maximum_fractional_bits = 1535;

[[noreturn]] void fail(core::FailureCode code) {
    throw PartitionAbort{{code, operation}};
}

struct ConditioningBudget {
    const PgConditioningLimits& limits;
    std::uint64_t steps = 0;

    void charge(std::uint64_t units = 1) {
        if (limits.cancellation.stop_requested())
            fail(core::FailureCode::Cancelled);
        if (units > limits.maximum_steps - steps)
            fail(core::FailureCode::ResourceExceeded);
        steps += units;
    }
};

template <class I>
struct EnclosedInteger {
    I lower{}, upper{};
};

[[nodiscard]] std::size_t symmetric(std::size_t i, std::size_t j) noexcept {
    return i <= j ? packed(i, j) : packed(j, i);
}

[[nodiscard]] RankCell dyadic_units(const ExactRational& r, std::uint32_t scale,
                                    ConditioningBudget& work) {
    work.charge();
    if (r.numerator_bits() > 1536 || r.denominator_bits() > 1536)
        fail(core::FailureCode::ResourceExceeded);
    const auto& backend = ExactFactory::value(r).backend();
    work.charge(std::uint64_t{backend.num().size()} + backend.denom().size());
    const auto value = rank_fraction<RankCell>(r);
    const auto exponent = static_cast<unsigned>(rank_bits(value.denominator) - 1);
    if (exponent > scale || value.denominator != (RankCell{1} << exponent))
        fail(core::FailureCode::InvalidInput);
    auto units = rank_abs(value.numerator);
    units <<= (scale - exponent);
    return value.numerator < 0 ? -units : units;
}

template <class I>
void accumulate(EnclosedInteger<RankWide>& out, const EnclosedInteger<I>& entry,
                const RankCell& factor, ConditioningBudget& work) {
    const bool positive = factor > 0;
    const RankWide widened = factor;
    out.lower += rank_product(work, widened, RankWide{positive ? entry.lower : entry.upper});
    out.upper += rank_product(work, widened, RankWide{positive ? entry.upper : entry.lower});
}

[[nodiscard]] RankWide certify_congruence(std::span<const EnclosedInteger<RankCell>> gram,
                                          std::span<const RankCell> transform,
                                          std::span<EnclosedInteger<RankWide>> intermediate,
                                          std::span<RankWide> rows, const RankWide& one,
                                          ConditioningBudget& work) {
    for (std::size_t j = 0; j < rows.size(); ++j) {
        work.charge(intermediate.size());
        std::ranges::fill(intermediate, EnclosedInteger<RankWide>{});
        for (std::size_t k = 0; k <= j; ++k) {
            work.charge();
            const auto& factor = transform[packed(k, j)];
            if (factor == 0)
                continue;
            for (std::size_t i = 0; i < rows.size(); ++i)
                accumulate(intermediate[i], gram[symmetric(i, k)], factor, work);
        }
        for (std::size_t i = 0; i <= j; ++i) {
            EnclosedInteger<RankWide> value{};
            for (std::size_t k = 0; k <= i; ++k) {
                work.charge();
                const auto& factor = transform[packed(k, i)];
                if (factor != 0)
                    accumulate(value, intermediate[k], factor, work);
            }
            if (i == j) {
                value.lower -= one;
                value.upper -= one;
            }
            work.charge();
            const auto bound = std::max(rank_abs(value.lower), rank_abs(value.upper));
            rows[i] += bound;
            if (i != j)
                rows[j] += bound;
        }
    }
    work.charge(rows.size());
    return *std::ranges::max_element(rows);
}

[[nodiscard]] const PgConditioningState&
state(const core::PlannedBuffer<PgConditioningState>& buffer) {
    if (buffer.size() != 1)
        fail(core::FailureCode::InvalidInput);
    return buffer.values().front();
}

[[nodiscard]] ExactRational entry(ReferenceArithmetic& math, const PgConditioningState& source,
                                  std::size_t i, std::size_t j) {
    math.charge();
    return rank_materialize_dyadic(source.transform.values()[packed(i, j)], source.transform_scale,
                                   math.leaf_limits, math.plan);
}

} // namespace

PgConditioningCertificate::PgConditioningCertificate(
    PgReferenceBasis&& basis, core::PlannedBuffer<PgConditioningState> storage) noexcept
    : basis_(std::move(basis)), state_(std::move(storage)) {}
PgConditioningCertificate::PgConditioningCertificate(PgConditioningCertificate&&) noexcept =
    default;
PgConditioningCertificate&
PgConditioningCertificate::operator=(PgConditioningCertificate&&) noexcept = default;
PgConditioningCertificate::~PgConditioningCertificate() = default;

core::Result<PgConditioningCertificate>
PgConditioningCertificate::certify(PgReferenceBasis&& basis, PgDiscreteGramEnclosure enclosure,
                                   PgConditioningCandidate candidate,
                                   const PgConditioningLimits& limits, core::Workspace& workspace,
                                   core::ResourcePlan& plan) noexcept {
    return certify_impl(std::move(basis), enclosure, candidate, limits, workspace, plan);
}

core::Result<PgConditioningCertificate> PgConditioningCertificate::generate_and_certify(
    PgReferenceBasis&& basis, PgDiscreteGramEnclosure enclosure,
    PgConditioningGeneration generation, const PgConditioningLimits& limits,
    core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return certify_impl(std::move(basis), enclosure, generation, limits, workspace, plan);
}

core::Result<PgConditioningCertificate> PgConditioningCertificate::certify_impl(
    PgReferenceBasis&& basis, PgDiscreteGramEnclosure enclosure,
    std::variant<PgConditioningCandidate, PgConditioningGeneration> input,
    const PgConditioningLimits& limits, core::Workspace& workspace,
    core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::InvalidInput, {},
        [&]() -> core::Result<PgConditioningCertificate> {
            try {
                const auto lease = partition_take(workspace.acquire());
                static_cast<void>(lease);
                auto work = ConditioningBudget{limits};
                work.charge(0);
                const auto* rank = partition_take(basis.rank());
                const auto info = partition_take(rank->preparation());
                if (info.rank == 0 || info.rank > rank_max_columns ||
                    info.independent_required > info.rank)
                    fail(core::FailureCode::RankUnresolved);
                const auto count = static_cast<std::size_t>(info.rank);
                const auto terms = packed(0, count);
                const auto* generation = std::get_if<PgConditioningGeneration>(&input);
                const auto* candidate = std::get_if<PgConditioningCandidate>(&input);
                const auto transform_scale = generation != nullptr
                                                 ? generation->candidate_fractional_bits
                                                 : candidate->fractional_bits;
                if (enclosure.upper_triangle.size() != terms ||
                    (candidate != nullptr && candidate->upper_triangle.size() != terms))
                    fail(core::FailureCode::InvalidInput);
                if (enclosure.fractional_bits > maximum_fractional_bits ||
                    transform_scale > maximum_fractional_bits ||
                    (generation != nullptr &&
                     generation->factor_fractional_bits > maximum_fractional_bits))
                    fail(core::FailureCode::ResourceExceeded);

                RankSizer forecast;
                forecast.add<PgConditioningState>(1);
                forecast.add<RankCell>(terms);
                forecast.add<EnclosedInteger<RankCell>>(terms);
                forecast.add<EnclosedInteger<RankWide>>(count);
                forecast.add<RankWide>(count);
                if (generation != nullptr) {
                    forecast.add<PgConditioningFactor>(terms);
                    forecast.add<PgConditioningFactor>(count);
                }
                if (forecast.bytes > limits.maximum_payload_bytes ||
                    forecast.bytes > plan.remaining_bytes())
                    fail(core::FailureCode::ResourceExceeded);
                work.charge(0);
                auto storage = partition_take(plan.allocate<PgConditioningState>(1, operation));
                auto& out = storage.values().front();
                work.charge(0);
                out.transform =
                    partition_take(plan.allocate<RankCell>(partition_index(terms), operation));
                work.charge(0);
                auto gram = partition_take(
                    plan.allocate<EnclosedInteger<RankCell>>(partition_index(terms), operation));
                work.charge(0);
                auto column = partition_take(
                    plan.allocate<EnclosedInteger<RankWide>>(partition_index(count), operation));
                work.charge(0);
                auto rows =
                    partition_take(plan.allocate<RankWide>(partition_index(count), operation));

                const auto transform = out.transform.values();
                const auto bounds = gram.values();
                for (std::size_t j = 0; j < count; ++j) {
                    for (std::size_t i = 0; i <= j; ++i) {
                        const auto k = packed(i, j);
                        const auto& interval = enclosure.upper_triangle[k];
                        bounds[k].lower =
                            dyadic_units(interval.lower, enclosure.fractional_bits, work);
                        bounds[k].upper =
                            dyadic_units(interval.upper, enclosure.fractional_bits, work);
                        if (candidate != nullptr)
                            transform[k] =
                                dyadic_units(candidate->upper_triangle[k], transform_scale, work);
                        if (bounds[k].lower > bounds[k].upper)
                            fail(core::FailureCode::InvalidInput);
                    }
                    if (candidate != nullptr && transform[packed(j, j)] <= 0)
                        fail(core::FailureCode::RankUnresolved);
                }

                if (generation != nullptr) {
                    work.charge(0);
                    auto factor = partition_take(
                        plan.allocate<PgConditioningFactor>(partition_index(terms), operation));
                    work.charge(0);
                    auto inverse = partition_take(
                        plan.allocate<PgConditioningFactor>(partition_index(count), operation));
                    conditioning_factor(work, std::span<const EnclosedInteger<RankCell>>{bounds},
                                        enclosure.fractional_bits,
                                        generation->factor_fractional_bits, count, factor.values());
                    conditioning_inverse(work, factor.values(), generation->factor_fractional_bits,
                                         transform_scale, inverse.values(), transform);
                    for (std::size_t j = 0; j < count; ++j) {
                        work.charge();
                        if (transform[packed(j, j)] <= 0)
                            fail(core::FailureCode::RankUnresolved);
                    }
                }

                out.transform_scale = transform_scale;
                out.congruence_scale = enclosure.fractional_bits + 2 * transform_scale;
                const RankWide one = RankWide{1} << out.congruence_scale;
                out.bound = certify_congruence(bounds, transform, column.values(), rows.values(),
                                               one, work);
                const RankWide ten_billion{10'000'000'000ULL};
                const auto tolerance_denominator = rank_product(work, ten_billion, ten_billion);
                if (rank_product(work, out.bound, tolerance_denominator) > one)
                    fail(core::FailureCode::PrecisionUnresolved);
                work.charge(0);
                out.info = {info.rank, info.independent_required, forecast.bytes, work.steps};
                return PgConditioningCertificate{std::move(basis), std::move(storage)};
            } catch (const PartitionAbort& failure) {
                return std::unexpected(failure.failure);
            } catch (const std::overflow_error&) {
                return std::unexpected(
                    core::Failure{core::FailureCode::ResourceExceeded, operation});
            } catch (const std::range_error&) {
                return std::unexpected(
                    core::Failure{core::FailureCode::ResourceExceeded, operation});
            }
        });
}

core::Result<PgConditioningInfo> PgConditioningCertificate::info() const noexcept {
    return partition_entry(
        [&]() -> core::Result<PgConditioningInfo> { return state(state_).info; });
}

core::Result<ExactRational>
PgConditioningCertificate::congruence_bound(const ExactLimits& limits,
                                            core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        const auto& accepted = state(state_);
        return rank_materialize_dyadic(accepted.bound, accepted.congruence_scale, limits, plan);
    });
}

core::Result<const PgReferenceBasis*> PgConditioningCertificate::reference_basis() const& noexcept {
    return partition_entry([&]() -> core::Result<const PgReferenceBasis*> {
        static_cast<void>(state(state_));
        return &basis_;
    });
}

core::Result<ExactRational>
PgConditioningCertificate::transform_coefficient(PgRankBasisColumn row, PgRankBasisColumn column,
                                                 const ExactLimits& limits,
                                                 core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        const auto& accepted = state(state_);
        const auto i = static_cast<std::uint64_t>(row.value);
        const auto j = static_cast<std::uint64_t>(column.value);
        if (i >= accepted.info.rank || j >= accepted.info.rank)
            fail(core::FailureCode::InvalidInput);
        if (i > j)
            return ExactRational{};
        return rank_materialize_dyadic(
            accepted.transform
                .values()[packed(static_cast<std::size_t>(i), static_cast<std::size_t>(j))],
            accepted.transform_scale, limits, plan);
    });
}

core::Result<PgReferenceBasisValue>
PgConditioningCertificate::evaluate(const PgReferencePoint& query, PgRankBasisColumn column,
                                    const PgReferenceLimits& limits,
                                    core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgReferenceBasisValue> {
        auto math = ReferenceArithmetic{limits, plan};
        // Besides the ordinary point path: nine accumulated values, an exact
        // weight, and the two arithmetic operands/result during accumulation.
        math.admit(basis_live_rationals + 16);
        const auto& accepted = state(state_);
        if (column.value < 0 || static_cast<std::uint64_t>(column.value) >= accepted.info.rank)
            fail(core::FailureCode::InvalidInput);
        auto point = basis_.locate(math, query);
        std::array<std::array<ExactRational, 3>, 3> combined{};
        const auto j = static_cast<std::size_t>(column.value);
        for (std::size_t i = 0; i <= j; ++i) {
            math.visit(1);
            if (accepted.transform.values()[packed(i, j)] == 0)
                continue;
            const auto weight = entry(math, accepted, i, j);
            const auto source = basis_.numerator_jet(math, {partition_index(i)}, point);
            for (std::size_t field = 0; field < combined.size(); ++field)
                for (std::size_t d = 0; d < combined[field].size(); ++d)
                    if (source[field][d].sign() != 0)
                        combined[field][d] =
                            math.add(combined[field][d], math.mul(weight, source[field][d]));
        }
        return basis_.assemble(math, std::move(point), column, combined);
    });
}

core::Result<detail::PgConditionedVolumeRow> detail::PgConditionedVolumeRow::evaluate(
    const PgConditioningCertificate& certificate, const PgReferenceVolumePoint& query,
    const PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgConditionedVolumeRow> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto& accepted = state(certificate.state_);
        const auto rank = accepted.info.rank;
        math.admit(partition_add(basis_live_rationals + 16 + 74, partition_mul(2, rank)));
        auto values = partition_take(
            plan.allocate<PgReferenceScalar>(partition_index(rank), partition_operation));
        const auto& basis = certificate.basis_;
        auto point = basis.locate(math, query);
        const auto two = math.constant("2");
        const auto z = math.mul(
            math.div(point.reference.full_thickness_metres, math.mul(two, basis.length_metres())),
            point.reference.zeta);
        const auto half_z_squared = math.div(math.square(z), two);
        const auto g_squared = math.square(point.reference.g);
        auto row = values.values();
        for (std::size_t i = 0; i < row.size(); ++i) {
            math.visit(1);
            const auto [r, p, c] =
                basis.numerator_jet(math, {partition_index(i)}, point, /*derivatives=*/false);
            row[i] = {math.add(r[0], math.mul(half_z_squared, c[0])),
                      math.mul(z, math.div(p[0], g_squared))};
        }
        ReferenceRadicals field{math, point.reference.g};
        for (std::size_t end = row.size(); end != 0; --end) {
            const auto j = end - 1;
            PgReferenceScalar combined{};
            for (std::size_t i = 0; i <= j; ++i) {
                math.visit(1);
                if (accepted.transform.values()[packed(i, j)] == 0)
                    continue;
                const auto weight = entry(math, accepted, i, j);
                combined = field.add(combined, field.scale(row[i], weight));
            }
            row[j] = std::move(combined);
        }
        math.cancel();
        point.reference.call_info = math.info();
        return PgConditionedVolumeRow{std::move(point.reference), std::move(values), math.info()};
    });
}

namespace {
[[nodiscard]] std::size_t composed_group(const PgRankRow& row) noexcept {
    return 3 * (2 * static_cast<std::size_t>(row.cell.first) +
                static_cast<std::size_t>(row.cell.second)) +
           static_cast<std::size_t>(row.field);
}

template <class I>
[[nodiscard]] I exact_divide(ReferenceArithmetic& math, const I& a, const I& d) {
    math.charge(partition_mul(a.backend().size(), d.backend().size()));
    return a / d;
}
} // namespace

std::optional<detail::PgComposedBasis>
detail::PgComposedBasis::prepare(const PgConditioningCertificate& certificate,
                                 ReferenceArithmetic& math) {
    math.cancel();
    const auto& accepted = state(certificate.state_);
    const auto source = partition_take(certificate.basis_.source());
    const auto graph = partition_take(source.dictionary.graph());
    // Fixed composed scratch stays inside the caller's exact work-bit budget.
    if (graph.n != 1 || accepted.transform_scale >= coefficient_bits ||
        math.limits.algebra.maximum_work_bits < std::numeric_limits<RankCell>::digits)
        return std::nullopt;
    const auto* rank = partition_take(certificate.basis_.rank());
    const auto rows = partition_take(rank->rows());
    const auto n = static_cast<std::size_t>(accepted.info.rank);
    const auto transform_slots = accepted.transform.values().size();
    if (partition_add(n, transform_slots) >
        math.limits.maximum_partition_visits - math.partition_visits)
        return std::nullopt;
    std::uint64_t used = 0, updates = 0;
    for (std::size_t i = 0; i < n; ++i) {
        math.visit(1);
        const auto indices = partition_take(rank->basis_indices({partition_index(i)}));
        used = partition_add(used, indices.size());
        updates = partition_add(updates, partition_mul(indices.size(), n - i));
    }
    constexpr auto limbs = coefficient_bits / std::numeric_limits<RankWord>::digits;
    const auto storage = partition_mul(partition_mul(rows.size(), n), sizeof(Units));
    const auto estimated_contract =
        partition_add(partition_mul(updates, 1 + limbs + limbs * limbs),
                      partition_mul(used, 1 + 2 * limbs + 2 * limbs * limbs));
    if (storage > math.plan.remaining_bytes() ||
        estimated_contract > math.limits.maximum_scalar_operations - math.scalar_operations ||
        partition_add(transform_slots, partition_mul(2, used)) >
            math.limits.maximum_partition_visits - math.partition_visits)
        return std::nullopt;
    std::uint64_t t_bits = 0;
    for (const auto& value : accepted.transform.values()) {
        math.visit(1);
        t_bits = std::max(t_bits, rank_bits(value));
    }
    if (t_bits > coefficient_bits)
        return std::nullopt;

    PgComposedBasis result;
    result.owner_ = &certificate;
    result.basis_ = &certificate.basis_;
    result.rows_ = rows;
    result.rank_ = n;
    result.denominators_.fill(Units{1});
    std::uint64_t b_bits = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto column = source_column(*rank, {partition_index(i)});
        for (std::size_t term = 0; term < column.rows.size(); ++term) {
            math.visit(1);
            math.charge(1 + 2 * limbs);
            const auto& value = column.values[term];
            const auto numerator_bits = rank_bits(value.numerator);
            const auto denominator_bits = rank_bits(value.denominator);
            auto& common = result.denominators_[composed_group(rows[column.rows[term]])];
            if (numerator_bits > coefficient_bits ||
                rank_bits(common) + denominator_bits > coefficient_bits)
                return std::nullopt;
            b_bits = std::max(b_bits, numerator_bits);
            common = rank_lcm(math, common, rank_copy_integer<Units>(value.denominator.backend()));
        }
    }
    for (const auto& common : result.denominators_) {
        if (rank_bits(common) + accepted.transform_scale > coefficient_bits ||
            b_bits + rank_bits(common) + t_bits + static_cast<std::uint64_t>(std::bit_width(n)) >
                coefficient_bits)
            return std::nullopt;
    }

    if (estimated_contract > math.limits.maximum_scalar_operations - math.scalar_operations)
        return std::nullopt;

    result.entries_ = partition_take(math.plan.allocate<Units>(
        partition_index(partition_mul(rows.size(), n)), partition_operation));
    for (std::size_t i = 0; i < n; ++i) {
        const auto column = source_column(*rank, {partition_index(i)});
        for (std::size_t term = 0; term < column.rows.size(); ++term) {
            math.visit(1);
            math.charge(1 + 2 * limbs);
            const auto index = static_cast<std::size_t>(column.rows[term]);
            const auto& fraction = column.values[term];
            const auto& common = result.denominators_[composed_group(rows[index])];
            const auto denominator = rank_copy_integer<Units>(fraction.denominator.backend());
            const auto scaled =
                rank_product(math, rank_copy_integer<Units>(fraction.numerator.backend()),
                             exact_divide(math, common, denominator));
            auto out = result.entries_.values().subspan(index * n, n);
            for (std::size_t j = i; j < n; ++j) {
                math.charge(1 + limbs);
                const auto& transform = accepted.transform.values()[packed(i, j)];
                if (transform != 0)
                    out[j] +=
                        rank_product(math, scaled, rank_copy_integer<Units>(transform.backend()));
            }
        }
    }
    for (auto& common : result.denominators_)
        common <<= accepted.transform_scale;
    return result;
}

std::optional<PgConditionedReferenceBatch>
detail::PgComposedBasis::evaluate(const PgConditioningCertificate& certificate,
                                  const PgReferenceVolumePoint& query,
                                  ReferenceArithmetic& math) const {
    math.cancel();
    if (&certificate != owner_)
        partition_fail(core::FailureCode::StaleGeneration);
    math.admit(partition_add(basis_live_rationals + 16 + 74, partition_mul(8, rank_)));
    auto point = basis_->locate(math, query);
    if (std::max({point.u.numerator_bits(), point.u.denominator_bits(), point.v.numerator_bits(),
                  point.v.denominator_bits()}) > coordinate_bits)
        return std::nullopt;
    const auto u = rank_fraction<RankCell>(point.u), v = rank_fraction<RankCell>(point.v);
    const auto common = rank_lcm(math, u.denominator, v.denominator);
    if (rank_bits(common) > coordinate_bits)
        return std::nullopt;
    const auto x = rank_product(math, u.numerator, exact_divide(math, common, u.denominator));
    const auto y = rank_product(math, v.numerator, exact_divide(math, common, v.denominator));
    if (rank_abs(x) > common || rank_abs(y) > common)
        return std::nullopt;

    const auto powers = [&](const RankCell& value) {
        std::array<RankCell, degree + 1> out;
        out.front() = 1;
        for (std::size_t j = 1; j < out.size(); ++j)
            out[j] = rank_product(math, out[j - 1], value);
        return out;
    };
    const auto xp = powers(x), yp = powers(y), dp = powers(common);
    using Jet = std::array<std::array<RankCell, 3>, 3>;
    // The temporary is allocated through the caller's bounded plan and its
    // signed accumulators are at most 512 + 21*32 + ceil(log2 rows) + 7 bits.
    auto sums =
        partition_take(math.plan.allocate<Jet>(partition_index(rank_), partition_operation));
    const auto cell = point.reference.cell;
    const auto base =
        3 * (2 * static_cast<std::size_t>(cell.first) + static_cast<std::size_t>(cell.second));
    constexpr auto width = std::numeric_limits<RankCell>::digits;
    constexpr auto limbs = width / std::numeric_limits<RankWord>::digits;
    for (std::size_t i = 0; i < rows_.size(); ++i) {
        math.visit(1);
        const auto& row = rows_[i];
        if (row.cell != cell)
            continue;
        const auto field = static_cast<std::size_t>(row.field);
        const auto a = static_cast<std::size_t>(row.powers.first);
        const auto b = static_cast<std::size_t>(row.powers.second);
        const auto monomial = [&](std::size_t s, std::size_t t) {
            return rank_product(math, rank_product(math, xp[s], yp[t]), dp[degree - s - t]);
        };
        const std::array scales{
            monomial(a, b),
            a == 0 ? RankCell{} : rank_product(math, RankCell{4 * a}, monomial(a - 1, b)),
            b == 0 ? RankCell{} : rank_product(math, RankCell{4 * b}, monomial(a, b - 1))};
        const auto entries = entries_.values().subspan(i * rank_, rank_);
        for (std::size_t j = 0; j < entries.size(); ++j) {
            math.charge();
            if (entries[j] == 0)
                continue;
            const auto coefficient = rank_copy_integer<RankCell>(entries[j].backend());
            for (std::size_t d = 0; d < scales.size(); ++d) {
                if (scales[d] == 0)
                    continue;
                math.charge(limbs);
                sums.values()[j][field][d] += rank_product(math, coefficient, scales[d]);
            }
        }
    }
    std::array<RankCell, 3> output_denominators;
    for (std::size_t f = 0; f < output_denominators.size(); ++f)
        output_denominators[f] = rank_product(
            math, rank_copy_integer<RankCell>(denominators_[base + f].backend()), dp.back());

    auto points = partition_take(math.plan.allocate<PgReferenceValue>(1, partition_operation));
    auto columns = partition_take(math.plan.allocate<PgConditionedReferenceColumn>(
        partition_index(rank_), partition_operation));
    for (std::size_t j = 0; j < rank_; ++j) {
        std::array<std::array<ExactRational, 3>, 3> jet;
        for (std::size_t f = 0; f < jet.size(); ++f)
            for (std::size_t d = 0; d < jet[f].size(); ++d) {
                const auto& numerator = sums.values()[j][f][d];
                if (numerator == 0)
                    continue;
                const auto& denominator = output_denominators[f];
                const auto divisor = rank_gcd(math, numerator, denominator);
                const RankFraction<RankCell> reduced{exact_divide(math, numerator, divisor),
                                                     exact_divide(math, denominator, divisor)};
                math.charge();
                jet[f][d] = rank_exact(reduced, math.leaf_limits, math.plan);
            }
        auto value = basis_->sample(math, point, jet);
        columns.values()[j] = {std::move(value.value),
                               std::move(value.physical_gradient_per_metre)};
    }
    math.cancel();
    point.reference.call_info = math.info();
    points.values().front() = std::move(point.reference);
    return PgConditionedReferenceBatch{rank_, std::move(points), std::move(columns), math.info()};
}

core::Result<PgConditionedReferenceBatch>
PgConditioningCertificate::evaluate_reference(std::span<const PgReferencePoint> queries,
                                              const PgReferenceLimits& limits,
                                              core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgConditionedReferenceBatch> {
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto& accepted = state(state_);
        if (queries.empty())
            return PgConditionedReferenceBatch{accepted.info.rank, {}, {}, math.info()};

        const auto point_count = partition_index(queries.size());
        const auto column_count =
            partition_index(partition_mul(queries.size(), accepted.info.rank));
        const auto point_bytes = partition_take(
            core::checked_byte_count(point_count, sizeof(PgReferenceValue), partition_operation));
        const auto column_bytes = partition_take(core::checked_byte_count(
            column_count, sizeof(PgConditionedReferenceColumn), partition_operation));
        if (partition_add(point_bytes, column_bytes) > plan.remaining_bytes())
            partition_fail(core::FailureCode::ResourceExceeded);
        // Per point: 22 chart rationals, 44 Cartesian/determinant rationals,
        // and up to 8 rationals for an area measure. Per column: 2 for N, 6 for G.
        const auto retained =
            partition_add(partition_mul(74, queries.size()),
                          partition_mul(8, static_cast<std::uint64_t>(column_count)));
        math.admit(partition_add(basis_live_rationals + 16, retained));
        math.visit(partition_add(queries.size(), static_cast<std::uint64_t>(column_count)));
        auto points =
            partition_take(plan.allocate<PgReferenceValue>(point_count, partition_operation));
        auto columns = partition_take(
            plan.allocate<PgConditionedReferenceColumn>(column_count, partition_operation));
        const auto rank = static_cast<std::size_t>(accepted.info.rank);
        auto unfilled = columns.values();
        for (std::size_t p = 0; p < queries.size(); ++p) {
            const auto start = math.info();
            math.visit(1);
            auto point = basis_.locate(math, queries[p]);
            auto row = unfilled.first(rank);
            unfilled = unfilled.subspan(rank);
            for (std::size_t i = 0; i < rank; ++i) {
                math.visit(1);
                const auto jet = basis_.numerator_jet(math, {partition_index(i)}, point);
                auto value = basis_.sample(math, point, jet);
                row[i] = {std::move(value.value), std::move(value.physical_gradient_per_metre)};
            }
            // N and G are linear in the selected basis. Traverse downwards so
            // each column can replace its input without losing any untransformed
            // row i <= j needed by the other columns of the upper triangular T.
            auto field = ReferenceRadicals{math, point.reference.g};
            for (std::size_t end = rank; end != 0; --end) {
                const auto j = end - 1;
                PgConditionedReferenceColumn combined{};
                for (std::size_t i = 0; i <= j; ++i) {
                    math.visit(1);
                    if (accepted.transform.values()[packed(i, j)] == 0)
                        continue;
                    const auto weight = entry(math, accepted, i, j);
                    combined.value = field.add(combined.value, field.scale(row[i].value, weight));
                    auto& gradient = combined.physical_gradient_per_metre;
                    for (std::size_t d = 0; d < gradient.size(); ++d)
                        gradient[d] =
                            field.add(gradient[d],
                                      field.scale(row[i].physical_gradient_per_metre[d], weight));
                }
                row[j] = std::move(combined);
            }
            const auto after = math.info();
            point.reference.call_info = {after.scalar_operations - start.scalar_operations,
                                         after.partition_visits - start.partition_visits,
                                         after.allocated_owned_payload_bytes -
                                             start.allocated_owned_payload_bytes};
            points.values()[p] = std::move(point.reference);
        }
        math.cancel();
        return PgConditionedReferenceBatch{accepted.info.rank, std::move(points),
                                           std::move(columns), math.info()};
    });
}

core::Result<ExactRational>
PgConditioningCertificate::coefficient(PgRankBasisColumn column, core::Index row,
                                       const PgReferenceLimits& limits,
                                       core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        auto math = ReferenceArithmetic{limits, plan};
        math.admit();
        const auto& accepted = state(state_);
        const auto* rank = partition_take(basis_.rank());
        const auto rows = partition_take(rank->rows());
        if (column.value < 0 || static_cast<std::uint64_t>(column.value) >= accepted.info.rank ||
            row < 0 || static_cast<std::uint64_t>(row) >= rows.size())
            fail(core::FailureCode::InvalidInput);
        auto result = ExactRational{};
        const auto j = static_cast<std::size_t>(column.value);
        for (std::size_t i = 0; i <= j; ++i) {
            math.visit(1);
            if (accepted.transform.values()[packed(i, j)] == 0)
                continue;
            math.visit(static_cast<std::uint64_t>(std::bit_width(rows.size())));
            math.charge();
            const auto value =
                partition_take(rank->basis({partition_index(i)}, row, math.leaf_limits, plan));
            result = math.add(result, math.mul(value, entry(math, accepted, i, j)));
        }
        math.cancel();
        return result;
    });
}

core::PlannedBuffer<ExactRational>
PgConditioningCertificate::materialize_required(ReferenceArithmetic& math,
                                                PgRankRequiredColumn target) const {
    const auto& accepted = state(state_);
    if (target.value < 0 || target.value >= rank_required)
        fail(core::FailureCode::InvalidInput);
    const auto* rank = partition_take(basis_.rank());
    auto output = partition_take(
        math.plan.allocate<ExactRational>(partition_index(accepted.info.rank), operation));
    auto values = output.values();
    const auto prefix = static_cast<std::size_t>(accepted.info.protected_rank);
    for (std::size_t offset = prefix; offset > 0; --offset) {
        const auto i = offset - 1;
        math.charge();
        auto residual = partition_take(
            rank->required({partition_index(i)}, target, math.leaf_limits, math.plan));
        for (std::size_t j = i + 1; j < prefix; ++j) {
            math.visit(1);
            if (accepted.transform.values()[packed(i, j)] != 0)
                residual = math.sub(residual, math.mul(entry(math, accepted, i, j), values[j]));
        }
        values[i] = math.div(residual, entry(math, accepted, i, i));
    }
    math.cancel();
    return output;
}

void PgConditioningCertificate::append_basis_identity(ReferenceArithmetic& math,
                                                      ReferenceHash& stream) const {
    const auto& accepted = state(state_);
    const auto geometry = partition_take(basis_.reference_identity());
    stream.bytes(geometry.bytes.data(), geometry.bytes.size());
    const auto* rank = partition_take(basis_.rank());
    const auto recipe = partition_take(rank->witness());
    stream.word(recipe.recipe_version);
    stream.word(recipe.reduction_version);
    stream.word(accepted.info.rank);
    stream.word(accepted.info.protected_rank);
    stream.word(accepted.transform_scale);
    for (const auto ordinal : partition_take(rank->selected())) {
        math.visit(1);
        stream.word(ordinal);
    }
    for (const auto& value : accepted.transform.values()) {
        math.visit(1);
        stream.signed_magnitude(value.backend(), rank_bits(value));
    }
}

core::Result<core::PlannedBuffer<ExactRational>> PgConditioningCertificate::required_coordinates(
    PgRankRequiredColumn target, const PgReferenceLimits& limits, core::Workspace& workspace,
    core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<ExactRational>> {
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        auto math = ReferenceArithmetic{limits, plan};
        math.cancel();
        const auto& accepted = state(state_);
        if (target.value < 0 || target.value >= rank_required)
            fail(core::FailureCode::InvalidInput);
        math.admit(reference_live_rationals + accepted.info.rank);
        return materialize_required(math, target);
    });
}

} // namespace uvlas::structure
