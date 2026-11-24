#pragma once

#include <uvlas/structure/pg_prepared_patch.hpp>
#include <uvlas/structure/pg_volume_position.hpp>

namespace uvlas::structure {

struct PgPreparedVolumeSample {
    PgVolumeMassPoint reference;
    std::size_t rule_ordinal = 0;
    ExactRational weight;
};

/// Owns a contiguous window of a patch's certified finite rule and projected N/J.
class PgPreparedVolumeReference {
public:
    PgPreparedVolumeReference(PgPreparedVolumeReference&&) noexcept = default;
    PgPreparedVolumeReference& operator=(PgPreparedVolumeReference&&) noexcept = default;
    PgPreparedVolumeReference(const PgPreparedVolumeReference&) = delete;
    PgPreparedVolumeReference& operator=(const PgPreparedVolumeReference&) = delete;

    [[nodiscard]] static core::Result<PgPreparedVolumeReference>
    prepare(const PgPreparedPatch&, std::size_t first, std::size_t count,
            const PgVolumeMassProjection&, const PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] PgPreparedPatchIdentity source() const noexcept {
        return source_;
    }
    [[nodiscard]] std::size_t column_count() const noexcept {
        return points_.size() == 0 ? 0 : columns_;
    }
    [[nodiscard]] PgVolumeMassProjection projection() const noexcept {
        return projection_;
    }
    [[nodiscard]] PgReferenceCallInfo call_info() const noexcept {
        return call_info_;
    }
    [[nodiscard]] std::span<const PgPreparedVolumeSample> points() const& noexcept {
        return points_.values();
    }
    [[nodiscard]] std::span<const PgPreparedVolumeSample> points() const&& = delete;
    [[nodiscard]] PgVolumeMassWords words() const& noexcept {
        return {high_.values(), low_.values()};
    }
    [[nodiscard]] PgVolumeMassWords words() const&& = delete;

private:
    PgPreparedVolumeReference(PgPreparedPatchIdentity, std::size_t, PgVolumeMassProjection,
                              PgReferenceCallInfo, core::PlannedBuffer<PgPreparedVolumeSample>,
                              core::PlannedBuffer<double>, core::PlannedBuffer<double>) noexcept;

    PgPreparedPatchIdentity source_;
    std::size_t columns_;
    PgVolumeMassProjection projection_;
    PgReferenceCallInfo call_info_;
    core::PlannedBuffer<PgPreparedVolumeSample> points_;
    core::PlannedBuffer<double> high_, low_;
};

[[nodiscard]] core::Result<void>
apply_pg_prepared_volume_position(const PgPreparedVolumeReference&, PgPreparedPatchIdentity,
                                  PgVolumePositionDirection, PgCartesianWords,
                                  const PgVolumePositionLimits&, core::Workspace&,
                                  PgVolumePositionResult&) noexcept;

} // namespace uvlas::structure
