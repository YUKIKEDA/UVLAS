#pragma once

#include <uvlas/structure/pg_discrete_gram.hpp>

namespace uvlas::structure {

/// Reference, exact selected basis, coordinate transform, and retained finite rule.
struct PgPreparedPatchIdentity {
    std::array<std::byte, 32> bytes{};
    bool operator==(const PgPreparedPatchIdentity&) const = default;
};

struct PgPreparedPatchInfo {
    PgPreparedPatchIdentity identity;
    PgReferenceIdentity geometry;
    PgConditioningInfo basis;
    PgReferenceCallInfo preparation;
};

using PgCartesianExactPosition = std::array<ExactRational, 3>;

/// Immutable reference patch. q0 is Cartesian metres, by column and axis.
class PgPreparedPatch {
public:
    PgPreparedPatch(PgPreparedPatch&&) noexcept = default;
    PgPreparedPatch& operator=(PgPreparedPatch&&) noexcept = default;
    PgPreparedPatch(const PgPreparedPatch&) = delete;
    PgPreparedPatch& operator=(const PgPreparedPatch&) = delete;

    /// Preserve the caller's accepted conditioning source on failure.
    [[nodiscard]] static core::Result<PgPreparedPatch> prepare(PgConditioningCertificate&&,
                                                               const PgReferenceLimits&,
                                                               core::Workspace&,
                                                               core::ResourcePlan&) noexcept;
    /// Preserve the exact finite rule alongside its certified basis.
    [[nodiscard]] static core::Result<PgPreparedPatch> prepare(PgFixedDiscreteGramCertificate&&,
                                                               const PgReferenceLimits&,
                                                               core::Workspace&,
                                                               core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgPreparedPatchInfo> info() const noexcept;

    [[nodiscard]] std::span<const PgCartesianExactPosition>
    reference_coordinates() const& noexcept {
        return reference_.values();
    }
    [[nodiscard]] std::span<const PgCartesianExactPosition>
    reference_coordinates() const&& = delete;

    /// Exact scalar coordinates of one, usable for translation along any Cartesian axis.
    [[nodiscard]] std::span<const ExactRational> translation() const& noexcept {
        return translation_.values();
    }
    [[nodiscard]] std::span<const ExactRational> translation() const&& = delete;

    [[nodiscard]] const PgConditioningCertificate& conditioning() const& noexcept;
    [[nodiscard]] const PgConditioningCertificate& conditioning() const&& = delete;

    [[nodiscard]] std::optional<std::span<const PgOwnedDiscreteVolumePoint>>
    fixed_rule() const& noexcept;
    [[nodiscard]] std::optional<std::span<const PgOwnedDiscreteVolumePoint>>
    fixed_rule() const&& = delete;

    [[nodiscard]] core::Result<PgConditionedReferenceBatch>
    evaluate_reference(std::span<const PgReferencePoint>, const PgReferenceLimits&,
                       core::ResourcePlan&) const noexcept;

private:
    struct Coordinates;
    using Conditioning = std::variant<PgConditioningCertificate, PgFixedDiscreteGramCertificate>;
    PgPreparedPatch(Conditioning, Coordinates) noexcept;

    [[nodiscard]] static Coordinates derive(const PgConditioningCertificate&,
                                            std::span<const PgOwnedDiscreteVolumePoint>,
                                            const PgReferenceLimits&, core::Workspace&,
                                            core::ResourcePlan&);

    Conditioning conditioning_;
    core::PlannedBuffer<PgCartesianExactPosition> reference_;
    core::PlannedBuffer<ExactRational> translation_;
    PgPreparedPatchInfo info_;
};

} // namespace uvlas::structure
