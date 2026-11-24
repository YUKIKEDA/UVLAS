#include "pg_reference_hash.hpp"

#include <string_view>
#include <uvlas/structure/pg_prepared_patch.hpp>

namespace uvlas::structure {
using namespace detail;

struct PgPreparedPatch::Coordinates {
    core::PlannedBuffer<PgCartesianExactPosition> reference;
    core::PlannedBuffer<ExactRational> translation;
    PgPreparedPatchInfo info;
};

PgPreparedPatch::Coordinates PgPreparedPatch::derive(
    const PgConditioningCertificate& conditioning, std::span<const PgOwnedDiscreteVolumePoint> rule,
    const PgReferenceLimits& limits, core::Workspace& workspace, core::ResourcePlan& plan) {
    auto math = ReferenceArithmetic{limits, plan};
    math.cancel();
    const auto lease = partition_take(workspace.acquire());
    static_cast<void>(lease);
    const auto accepted = partition_take(conditioning.info());
    const auto* basis = partition_take(conditioning.reference_basis());
    const auto geometry = partition_take(basis->reference_identity());
    const auto frame = partition_take(basis->frame());
    const auto count = partition_index(accepted.rank);
    const auto result_bytes = partition_mul(accepted.rank, sizeof(PgCartesianExactPosition));
    const auto required_bytes =
        partition_mul(partition_mul(7, accepted.rank), sizeof(ExactRational));
    if (partition_add(result_bytes, required_bytes) > plan.remaining_bytes())
        partition_fail(core::FailureCode::ResourceExceeded);
    // q0, translation, one required coordinate column, and frame/accumulation work.
    math.admit(partition_add(reference_live_rationals + 16, partition_mul(5, accepted.rank)));
    math.visit(accepted.rank);
    auto reference =
        partition_take(plan.allocate<PgCartesianExactPosition>(count, partition_operation));
    auto translation = conditioning.materialize_required(math, {30});
    auto rows = reference.values();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        math.visit(1);
        for (std::size_t d = 0; d < frame.origin_metres.size(); ++d)
            rows[i][d] = math.mul(frame.origin_metres[d], translation.values()[i]);
    }

    constexpr core::Index normal_start = 27;
    for (const core::Index target : {core::Index{0}, core::Index{1}, core::Index{2}, normal_start,
                                     normal_start + 1, normal_start + 2}) {
        const auto axis = static_cast<std::size_t>(target % 3);
        std::array<ExactRational, 3> factor;
        for (std::size_t d = 0; d < factor.size(); ++d)
            factor[d] = math.mul(frame.length_metres, frame.rotation[3 * d + axis]);
        const auto coordinates = conditioning.materialize_required(math, {target});
        for (std::size_t i = 0; i < rows.size(); ++i) {
            math.visit(1);
            const auto& coefficient = coordinates.values()[i];
            if (coefficient.sign() == 0)
                continue;
            for (std::size_t d = 0; d < factor.size(); ++d)
                rows[i][d] = math.add(rows[i][d], math.mul(factor[d], coefficient));
        }
    }

    ReferenceHash stream{math};
    constexpr std::string_view domain = "uvlas.pg.prepared-patch.v1";
    stream.bytes(domain.data(), domain.size());
    conditioning.append_basis_identity(math, stream);
    stream.word(rule.size());
    for (const auto& p : rule) {
        math.visit(1);
        stream.word(static_cast<std::uint64_t>(p.volume));
        stream.rational(p.u);
        stream.rational(p.v);
        stream.rational(p.depth);
        stream.rational(p.weight);
    }
    math.cancel();
    return {std::move(reference),
            std::move(translation),
            {{stream.result().bytes}, geometry, accepted, math.info()}};
}

PgPreparedPatch::PgPreparedPatch(Conditioning conditioning, Coordinates coordinates) noexcept
    : conditioning_(std::move(conditioning)),
      reference_(std::move(coordinates.reference)),
      translation_(std::move(coordinates.translation)),
      info_(coordinates.info) {}

core::Result<PgPreparedPatch> PgPreparedPatch::prepare(PgConditioningCertificate&& source,
                                                       const PgReferenceLimits& limits,
                                                       core::Workspace& workspace,
                                                       core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgPreparedPatch> {
        auto coordinates = derive(source, {}, limits, workspace, plan);
        return PgPreparedPatch{std::move(source), std::move(coordinates)};
    });
}

core::Result<PgPreparedPatch> PgPreparedPatch::prepare(PgFixedDiscreteGramCertificate&& source,
                                                       const PgReferenceLimits& limits,
                                                       core::Workspace& workspace,
                                                       core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgPreparedPatch> {
        auto coordinates = derive(source.conditioning(), source.rule(), limits, workspace, plan);
        return PgPreparedPatch{std::move(source), std::move(coordinates)};
    });
}

core::Result<PgPreparedPatchInfo> PgPreparedPatch::info() const noexcept {
    if (const auto present = conditioning().info(); !present)
        return std::unexpected(present.error());
    return info_;
}

const PgConditioningCertificate& PgPreparedPatch::conditioning() const& noexcept {
    if (const auto* fixed = std::get_if<PgFixedDiscreteGramCertificate>(&conditioning_))
        return fixed->conditioning();
    // Both alternatives have nonthrowing move construction and assignment.
    return *std::get_if<PgConditioningCertificate>(&conditioning_);
}

std::optional<std::span<const PgOwnedDiscreteVolumePoint>>
PgPreparedPatch::fixed_rule() const& noexcept {
    if (const auto* fixed = std::get_if<PgFixedDiscreteGramCertificate>(&conditioning_))
        return fixed->rule();
    return std::nullopt;
}

core::Result<PgConditionedReferenceBatch>
PgPreparedPatch::evaluate_reference(std::span<const PgReferencePoint> queries,
                                    const PgReferenceLimits& limits,
                                    core::ResourcePlan& plan) const noexcept {
    return conditioning().evaluate_reference(queries, limits, plan);
}

} // namespace uvlas::structure
