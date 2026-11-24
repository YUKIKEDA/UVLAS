#include "pg_rank_storage.hpp"
#include "pg_reference_polynomials.hpp"

#include <tuple>
#include <uvlas/structure/pg_flat_gram.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;

struct GramWork {
    ReferenceArithmetic& math;
    std::uint64_t remaining;

    void charge(std::uint64_t units = 1) {
        math.cancel();
        if (units > remaining)
            partition_fail(core::FailureCode::ResourceExceeded);
        remaining -= units;
    }
};

struct GramTerm {
    PgRankRow key;
    RankFraction<RankCell> raw;
    RankWide weighted;
};

[[nodiscard]] RankWide product(GramWork& work, const RankWide& x, const RankWide& y) {
    return rank_product(work, x, y);
}

[[nodiscard]] RankWide divide(GramWork& work, const RankWide& a, const RankWide& b) {
    work.charge(partition_mul(a.backend().size(), b.backend().size()));
    return RankWide{a / b};
}

[[nodiscard]] const ExactRational& full_flat_thickness(ReferenceArithmetic& math,
                                                       const PgReferenceSource& source,
                                                       const PgGraphView& graph) {
    for (const auto& term : graph.height.terms()) {
        math.visit(1);
        if (term.powers.first != 0 || term.powers.second != 0)
            partition_fail(core::FailureCode::Unsupported);
    }
    const auto full_area =
        math.div(math.constant("1"), math.square(index_value(math, 2 * graph.n)));
    const auto activity = partition_take(source.partition.activity_map());
    if (activity.size() != graph.cells.size())
        partition_fail(core::FailureCode::InvalidInput);
    const auto charge = [&](const ExactRational& x, const ExactRational& y) {
        math.visit(1);
        const auto words = [](std::uint64_t bits) { return bits / 32 + (bits % 32 != 0 ? 1 : 0); };
        math.visit(words(x.numerator_bits()) + words(x.denominator_bits()) +
                   words(y.numerator_bits()) + words(y.denominator_bits()));
    };
    for (const auto& cell : activity) {
        charge(cell.area, full_area);
        if (!cell.active || cell.area != full_area)
            partition_fail(core::FailureCode::Unsupported);
    }
    const auto thickness = partition_take(source.partition.thickness()).values;
    if (thickness.empty())
        partition_fail(core::FailureCode::InvalidInput);
    for (const auto& value : thickness) {
        charge(value, thickness.front());
        if (value != thickness.front())
            partition_fail(core::FailureCode::Unsupported);
    }
    return thickness.front();
}

[[nodiscard]] auto cell_key(const GramTerm& term) {
    return std::pair{term.key.cell.first, term.key.cell.second};
}

constexpr auto maximum_moment_degree = 2 * rank_degree[1];

struct GramScale {
    std::array<RankWide, maximum_moment_degree + 1> moment{};
    RankWide measure_numerator;
    RankWide denominator;
};

[[nodiscard]] GramScale setup_scale(GramWork& work, const ExactRational& h,
                                    const ExactRational& measure, const RankWide& common,
                                    std::span<GramTerm> terms) {
    if (h.numerator_bits() > 1536 || h.denominator_bits() > 1536 ||
        measure.numerator_bits() > 1536 || measure.denominator_bits() > 1536)
        partition_fail(core::FailureCode::ResourceExceeded);
    const auto half_thickness = rank_fraction<RankWide>(h);
    const auto volume = rank_fraction<RankWide>(measure);
    const auto p = half_thickness.numerator;
    const auto q = half_thickness.denominator;
    const std::array<RankWide, 3> field{product(work, RankWide{2}, product(work, q, q)),
                                        product(work, RankWide{2}, product(work, p, q)),
                                        product(work, p, p)};
    for (auto& term : terms) {
        work.charge();
        const auto ordinal = static_cast<std::size_t>(term.key.field);
        const auto integer = product(work, RankWide{term.raw.numerator},
                                     divide(work, common, RankWide{term.raw.denominator}));
        term.weighted = product(work, integer, field[ordinal]);
    }
    const auto value_denominator = product(work, common, field.front());
    auto moment_denominator = RankWide{1};
    for (std::size_t k = 0; k <= maximum_moment_degree; k += 2)
        moment_denominator = rank_lcm(work, moment_denominator, RankWide{k + 1});
    GramScale scale{};
    for (std::size_t k = 0; k < scale.moment.size(); k += 2)
        scale.moment[k] =
            product(work, RankWide{2}, divide(work, moment_denominator, RankWide{k + 1}));
    auto denominator = volume.denominator;
    for (const auto* factor : {&value_denominator, &value_denominator})
        denominator = product(work, denominator, *factor);
    for (std::size_t axis = 0; axis < 3; ++axis)
        denominator = product(work, denominator, moment_denominator);
    scale.measure_numerator = volume.numerator;
    scale.denominator = std::move(denominator);
    return scale;
}

[[nodiscard]] RankWide integrate(GramWork& work, const GramScale& scale,
                                 std::span<const GramTerm> a, std::span<const GramTerm> b) {
    RankWide total{0};
    while (!a.empty() && !b.empty()) {
        work.charge();
        const auto x = cell_key(a.front()), y = cell_key(b.front());
        if (x < y) {
            a = a.subspan(1);
            continue;
        }
        if (y < x) {
            b = b.subspan(1);
            continue;
        }
        std::size_t ac = 0, bc = 0;
        while (ac < a.size() && cell_key(a[ac]) == x) {
            work.charge();
            ++ac;
        }
        while (bc < b.size() && cell_key(b[bc]) == y) {
            work.charge();
            ++bc;
        }
        for (const auto& p : a.first(ac))
            for (const auto& q : b.first(bc)) {
                work.charge();
                const auto u = static_cast<std::size_t>(p.key.powers.first + q.key.powers.first);
                const auto v = static_cast<std::size_t>(p.key.powers.second + q.key.powers.second);
                const auto z =
                    static_cast<std::size_t>(p.key.field) + static_cast<std::size_t>(q.key.field);
                if (scale.moment[u] == 0 || scale.moment[v] == 0 || scale.moment[z] == 0)
                    continue;
                auto term = product(work, p.weighted, q.weighted);
                term = product(work, term, scale.moment[u]);
                term = product(work, term, scale.moment[v]);
                total += product(work, term, scale.moment[z]);
            }
        a = a.subspan(ac);
        b = b.subspan(bc);
    }
    return product(work, total, scale.measure_numerator);
}

[[nodiscard]] PgGramEntry enclose(GramWork& work, const RankWide& numerator,
                                  const RankWide& denominator, std::uint32_t bits) {
    work.charge();
    const RankWide unsigned_scaled = rank_abs(numerator) << bits;
    const auto quotient = divide(work, unsigned_scaled, denominator);
    const bool remainder = product(work, quotient, denominator) != unsigned_scaled;
    auto lower = RankWide{numerator < 0 ? -quotient : quotient};
    auto upper = lower;
    if (remainder) {
        if (numerator < 0)
            --lower;
        else
            ++upper;
    }
    if (rank_bits(lower) > 1536 || rank_bits(upper) > 1536)
        partition_fail(core::FailureCode::ResourceExceeded);
    auto& m = work.math;
    m.charge();
    auto lo = rank_materialize_dyadic(lower, bits, m.leaf_limits, m.plan);
    m.charge();
    auto hi = rank_materialize_dyadic(upper, bits, m.leaf_limits, m.plan);
    return {std::move(lo), std::move(hi)};
}

} // namespace

PgFlatGram::PgFlatGram(PgReferenceBasis&& basis, core::PlannedBuffer<PgGramEntry> entries,
                       std::uint32_t fractional_bits) noexcept
    : basis_(std::move(basis)), entries_(std::move(entries)), fractional_bits_(fractional_bits) {}

core::Result<PgFlatGram> PgFlatGram::prepare(PgReferenceBasis&& basis, std::uint32_t bits,
                                             const PgFlatGramLimits& limits, core::Workspace& ws,
                                             core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgFlatGram> {
        try {
            const auto lease = partition_take(ws.acquire());
            static_cast<void>(lease);
            auto math = ReferenceArithmetic{limits.arithmetic, plan};
            math.cancel();
            if (bits > 1535)
                partition_fail(core::FailureCode::ResourceExceeded);
            const auto reference = partition_take(basis.source());
            const auto graph = partition_take(reference.graph.view());
            const auto* rank = partition_take(basis.rank());
            const auto rank_count = partition_take(basis.size());
            const auto entry_count = partition_mul(rank_count, rank_count + 1) / 2;
            // The enclosure owns up to two rational values per packed entry.
            math.admit(partition_add(reference_live_rationals, partition_mul(entry_count, 2)));
            const auto& thickness = full_flat_thickness(math, reference, graph);
            auto work = GramWork{math, limits.maximum_integer_steps};
            const auto h = math.div(thickness, math.mul(math.constant("2"), basis.length_metres()));
            const auto four_n = index_value(math, 4 * graph.n);
            const auto measure = math.div(math.mul(graph.width, h), math.square(four_n));
            const auto rows = partition_take(rank->rows());
            std::uint64_t nnz = 0;
            for (std::uint64_t k = 0; k < rank_count; ++k) {
                work.charge();
                nnz = partition_add(
                    nnz, partition_take(rank->basis_indices({partition_index(k)})).size());
            }
            RankSizer forecast;
            forecast.add<GramTerm>(nnz);
            forecast.add<std::size_t>(rank_count + 1);
            forecast.add<PgGramEntry>(entry_count);
            if (forecast.bytes > plan.remaining_bytes())
                partition_fail(core::FailureCode::ResourceExceeded);
            auto terms = rank_allocate<GramTerm>(plan, nnz);
            auto offsets = rank_allocate<std::size_t>(plan, rank_count + 1);
            auto output = rank_allocate<PgGramEntry>(plan, entry_count);
            auto common = RankWide{1};
            std::size_t at = 0;
            for (std::size_t k = 0; k < rank_count; ++k) {
                const PgRankBasisColumn col{partition_index(k)};
                const auto support = partition_take(rank->basis_indices(col));
                for (const auto i : support) {
                    math.visit(static_cast<std::uint64_t>(std::bit_width(support.size())));
                    math.charge();
                    const auto value = partition_take(rank->basis(col, i, math.leaf_limits, plan));
                    auto raw = rank_fraction<RankCell>(value);
                    common = rank_lcm(work, common, RankWide{raw.denominator});
                    terms.values()[at++] = {rows[i], std::move(raw), {}};
                }
                const auto start = offsets.values()[k];
                std::ranges::sort(
                    terms.values().subspan(start, at - start), [&](const auto& a, const auto& b) {
                        work.charge();
                        const auto key = [](const GramTerm& term) {
                            return std::tuple{cell_key(term), term.key.field, term.key.powers.first,
                                              term.key.powers.second};
                        };
                        return key(a) < key(b);
                    });
                offsets.values()[k + 1] = at;
            }
            const auto scale = setup_scale(work, h, measure, common, terms.values());
            const auto column = [&](std::size_t k) {
                const auto begin = offsets.values()[k];
                return terms.values().subspan(begin, offsets.values()[k + 1] - begin);
            };
            std::size_t entry = 0;
            for (std::size_t j = 0; j < rank_count; ++j)
                for (std::size_t i = 0; i <= j; ++i) {
                    const auto value = integrate(work, scale, column(i), column(j));
                    output.values()[entry++] = enclose(work, value, scale.denominator, bits);
                }
            math.cancel();
            return PgFlatGram{std::move(basis), std::move(output), bits};
        } catch (const std::overflow_error&) {
            partition_fail(core::FailureCode::ResourceExceeded);
        } catch (const std::range_error&) {
            partition_fail(core::FailureCode::ResourceExceeded);
        }
    });
}

core::Result<PgDiscreteGramEnclosure> PgFlatGram::enclosure() const& noexcept {
    return partition_entry([&]() -> core::Result<PgDiscreteGramEnclosure> {
        static_cast<void>(partition_take(basis_.size()));
        return PgDiscreteGramEnclosure{entries_.values(), fractional_bits_};
    });
}

core::Result<PgConditioningCertificate> PgFlatGram::certify(PgConditioningCandidate candidate,
                                                            const PgConditioningLimits& limits,
                                                            core::Workspace& ws,
                                                            core::ResourcePlan& plan) && noexcept {
    return partition_entry([&]() -> core::Result<PgConditioningCertificate> {
        const auto input = partition_take(enclosure());
        auto accepted = PgConditioningCertificate::certify(std::move(basis_), input, candidate,
                                                           limits, ws, plan);
        if (accepted)
            entries_ = {};
        return accepted;
    });
}

} // namespace uvlas::structure
