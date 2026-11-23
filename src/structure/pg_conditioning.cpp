#include "pg_rank_storage.hpp"
#include "pg_reference_arithmetic.hpp"

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

[[nodiscard]] std::size_t packed(std::size_t i, std::size_t j) noexcept {
    return j * (j + 1) / 2 + i;
}

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
[[nodiscard]] ExactRational materialize(const I& units, std::uint32_t scale,
                                        const ExactLimits& limits, core::ResourcePlan& plan) {
    if (units == 0)
        return {};
    const auto magnitude = rank_abs(units);
    const auto trailing = std::min<std::uint64_t>(scale, boost::multiprecision::lsb(magnitude));
    auto numerator = I{magnitude >> trailing};
    if (units < 0)
        numerator = -numerator;
    return rank_exact(RankFraction<I>{std::move(numerator), I{I{1} << (scale - trailing)}}, limits,
                      plan);
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
    return materialize(source.transform.values()[packed(i, j)], source.transform_scale,
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
                if (enclosure.upper_triangle.size() != terms ||
                    candidate.upper_triangle.size() != terms)
                    fail(core::FailureCode::InvalidInput);
                if (enclosure.fractional_bits > maximum_fractional_bits ||
                    candidate.fractional_bits > maximum_fractional_bits)
                    fail(core::FailureCode::ResourceExceeded);

                RankSizer forecast;
                forecast.add<PgConditioningState>(1);
                forecast.add<RankCell>(terms);
                forecast.add<EnclosedInteger<RankCell>>(terms);
                forecast.add<EnclosedInteger<RankWide>>(count);
                forecast.add<RankWide>(count);
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
                        transform[k] = dyadic_units(candidate.upper_triangle[k],
                                                    candidate.fractional_bits, work);
                        if (bounds[k].lower > bounds[k].upper)
                            fail(core::FailureCode::InvalidInput);
                    }
                    if (transform[packed(j, j)] <= 0)
                        fail(core::FailureCode::RankUnresolved);
                }

                out.transform_scale = candidate.fractional_bits;
                out.congruence_scale = enclosure.fractional_bits + 2 * candidate.fractional_bits;
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
        return materialize(accepted.bound, accepted.congruence_scale, limits, plan);
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
        return materialize(accepted.transform.values()[packed(static_cast<std::size_t>(i),
                                                              static_cast<std::size_t>(j))],
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
        const auto* rank = partition_take(basis_.rank());
        math.admit(reference_live_rationals + accepted.info.rank);
        auto output = partition_take(
            plan.allocate<ExactRational>(partition_index(accepted.info.rank), operation));
        auto values = output.values();
        const auto prefix = static_cast<std::size_t>(accepted.info.protected_rank);
        for (std::size_t offset = prefix; offset > 0; --offset) {
            const auto i = offset - 1;
            math.charge();
            auto residual = partition_take(
                rank->required({partition_index(i)}, target, math.leaf_limits, plan));
            for (std::size_t j = i + 1; j < prefix; ++j) {
                math.visit(1);
                if (accepted.transform.values()[packed(i, j)] != 0)
                    residual = math.sub(residual, math.mul(entry(math, accepted, i, j), values[j]));
            }
            values[i] = math.div(residual, entry(math, accepted, i, i));
        }
        math.cancel();
        return output;
    });
}

} // namespace uvlas::structure
