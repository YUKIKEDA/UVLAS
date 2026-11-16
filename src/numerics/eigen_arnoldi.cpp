#include "backend_memory_detail.hpp"
#include "eigen_arnoldi_data.hpp"
#include "eigen_arnoldi_entry.hpp"
#include "eigen_full_phase_private.hpp"

#include <boost/hash2/sha2.hpp>

#include <array>
#include <bit>
#include <limits>
#include <utility>
#include <uvlas/numerics/eigen_arnoldi.hpp>

namespace uvlas::numerics::eigen::detail {
core::Result<void> arnoldi_normalize(std::span<ArnoldiComplex> vector,
                                     std::span<const double> scales, ArnoldiSide side,
                                     std::span<ArnoldiComplex> backup,
                                     std::span<ArnoldiComplex> candidate) {
    double norm = 0;
    for (std::size_t i = 0; i < vector.size(); ++i) {
        const auto coordinate =
            side == ArnoldiSide::Right ? vector[i] / scales[i] : vector[i] * scales[i];
        if (!arnoldi_finite(coordinate)) {
            return std::unexpected(
                core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
        }
        norm = std::hypot(norm, std::abs(coordinate));
    }
    if (!(norm > 0) || !std::isfinite(norm)) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    for (auto& value : vector) {
        value /= norm;
        if (!arnoldi_finite(value)) {
            return std::unexpected(
                core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
        }
    }
    std::ranges::copy(vector, backup.begin());
    // Rotation rounds too. Reuse the fixed, audited exact magnitude comparison
    // from Full, then verify the anchor against the actual rounded output.
    // The real/complex provider I/O planes are dead until the next action and
    // supply these scratch spans without new allocation.
    for (std::size_t pivot = 0; pivot < vector.size(); ++pivot) {
        const auto z = backup[pivot];
        const auto largest = std::max(std::abs(z.real()), std::abs(z.imag()));
        if (largest == 0) {
            continue;
        }
        const auto normalized = z / largest;
        const auto rotation = std::conj(normalized) / std::abs(normalized);
        const auto magnitude = std::abs(z);
        if (!std::isfinite(magnitude)) {
            continue;
        }
        bool finite = true;
        for (std::size_t i = 0; i < vector.size(); ++i) {
            candidate[i] = backup[i] * rotation;
            finite = finite && arnoldi_finite(candidate[i]);
        }
        if (!finite) {
            continue;
        }
        candidate[pivot] = {magnitude, 0};
        std::size_t maximum = 0;
        for (std::size_t i = 1; i < vector.size(); ++i) {
            if (full_phase_coordinate_larger(
                    candidate[i], scales[i], candidate[maximum], scales[maximum],
                    side == ArnoldiSide::Right ? PairSide::Right : PairSide::Left)) {
                maximum = i;
            }
        }
        if (maximum == pivot) {
            std::ranges::copy(candidate, vector.begin());
            return {};
        }
    }
    return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
}

} // namespace uvlas::numerics::eigen::detail

namespace uvlas::numerics::eigen {
ArnoldiWorkspace::ArnoldiWorkspace(std::unique_ptr<detail::ArnoldiWorkspaceData> data) noexcept
    : data_(std::move(data)) {}
ArnoldiWorkspace::ArnoldiWorkspace(ArnoldiWorkspace&&) noexcept = default;
ArnoldiWorkspace& ArnoldiWorkspace::operator=(ArnoldiWorkspace&&) noexcept = default;
ArnoldiWorkspace::~ArnoldiWorkspace() = default;
core::Index ArnoldiWorkspace::order() const noexcept {
    return data_ ? data_->order : 0;
}
core::Index ArnoldiWorkspace::capacity() const noexcept {
    return data_ ? data_->capacity : 0;
}
core::Result<ArnoldiWorkspace> ArnoldiWorkspace::create(core::Index order, core::Index capacity,
                                                        core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        core::Operation::ExtractEigen, core::FailureCode::InvalidInput, {},
        [&]() -> core::Result<ArnoldiWorkspace> {
            if (order < 0 || capacity < 0) {
                return std::unexpected(
                    core::Failure{core::FailureCode::InvalidInput, core::Operation::ExtractEigen});
            }
            const std::array shape{order, capacity}, square{capacity, capacity};
            const auto entries = core::checked_extent_product(shape, core::Operation::ExtractEigen);
            const auto small = core::checked_extent_product(square, core::Operation::ExtractEigen);
            if (!entries) {
                return std::unexpected(entries.error());
            }
            if (!small) {
                return std::unexpected(small.error());
            }
            auto data = std::make_unique<detail::ArnoldiWorkspaceData>();
            data->order = order;
            data->capacity = capacity;
            const auto allocate = [&](auto& destination, core::Index count) {
                return detail::arnoldi_allocate(destination, count, plan);
            };
            if (auto r = allocate(data->real_input, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->real_input_low, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->real_high, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->real_low, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->real_error, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->input, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->input_low, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->high, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->low, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->error, order); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->candidate_values, capacity); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->candidate_vectors, *entries); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->cluster_action, *entries); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->cluster_error, *entries); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->triangular_error, *small); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->assigned, capacity); !r) {
                return std::unexpected(r.error());
            }
            // The left candidate has only scratch values, no certificates.
            if (auto r = allocate(data->left_candidates.eigenvalues, capacity); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->left_candidates.left, *entries); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->left_candidates.left_t, *small); !r) {
                return std::unexpected(r.error());
            }
            if (auto r = allocate(data->left_candidates.left_mapping, capacity); !r) {
                return std::unexpected(r.error());
            }
            return ArnoldiWorkspace{std::move(data)};
        });
}

ArnoldiSpectrum::ArnoldiSpectrum(std::unique_ptr<detail::ArnoldiData> data) noexcept
    : data_(std::move(data)) {}
ArnoldiSpectrum::ArnoldiSpectrum(ArnoldiSpectrum&&) noexcept = default;
ArnoldiSpectrum& ArnoldiSpectrum::operator=(ArnoldiSpectrum&&) noexcept = default;
ArnoldiSpectrum::~ArnoldiSpectrum() = default;
bool ArnoldiSpectrum::active() const noexcept {
    return data_ != nullptr;
}
std::optional<ArnoldiCandidateView> ArnoldiSpectrum::retained_values() const& noexcept {
    if (!data_) {
        return std::nullopt;
    }
    return detail::arnoldi_view(*data_, data_->retained,
                                detail::retained_generation(data_->returned));
}
std::optional<ArnoldiCandidateView> ArnoldiSpectrum::values() const& noexcept {
    if (!data_) {
        return std::nullopt;
    }
    if (const auto* projection = std::get_if<linear::Binary64Projection>(&data_->returned)) {
        return detail::arnoldi_view(*data_, data_->projection, projection->target_generation);
    }
    return retained_values();
}
std::optional<ArnoldiMetadataView> ArnoldiSpectrum::metadata() const& noexcept {
    if (!data_) {
        return std::nullopt;
    }
    const auto& d = *data_;
    return ArnoldiMetadataView{d.context,
                               d.original_stamp,
                               d.state,
                               d.components.values(),
                               d.units.values(),
                               d.eigenvalue_scale,
                               d.eigenvalue_units,
                               detail::arnoldi_matrix(d.initial, d.order, 1),
                               d.association,
                               d.selection,
                               d.backend,
                               d.returned,
                               d.settings,
                               d.limits,
                               d.required_checks.size()};
}
std::optional<ArnoldiQualificationState> ArnoldiSpectrum::qualification() const noexcept {
    if (!data_) {
        return std::nullopt;
    }
    if (data_->subspace_requested) {
        return ArnoldiQualificationState{ArnoldiRequestedSubspaceQualified{}};
    }
    return ArnoldiQualificationState{ArnoldiRequestedSubsetQualified{}};
}
std::optional<PairMetrics> ArnoldiSpectrum::pair_metric(ArnoldiSide side,
                                                        core::Index index) const noexcept {
    if (!data_ || data_->subspace_requested || index < 0 ||
        (side != ArnoldiSide::Right && side != ArnoldiSide::Left) ||
        (side == ArnoldiSide::Left && !data_->left_requested)) {
        return std::nullopt;
    }
    const auto& packet = data_->projected ? data_->projection : data_->retained;
    const auto& proofs = side == ArnoldiSide::Right ? packet.right_pairs : packet.left_pairs;
    if (index >= packet.count) {
        return std::nullopt;
    }
    return proofs.values()[static_cast<std::size_t>(index)];
}
std::optional<BasisMetrics> ArnoldiSpectrum::subspace_metric(ArnoldiSide side) const noexcept {
    if (!data_ || !data_->subspace_requested ||
        (side != ArnoldiSide::Right && side != ArnoldiSide::Left) ||
        (side == ArnoldiSide::Left && !data_->left_requested)) {
        return std::nullopt;
    }
    const auto& packet = data_->projected ? data_->projection : data_->retained;
    return side == ArnoldiSide::Right ? packet.right_basis : packet.left_basis;
}
std::span<const ArnoldiSensitivity> ArnoldiSpectrum::sensitivities() const& noexcept {
    if (!data_ || data_->subspace_requested) {
        return {};
    }
    return data_->sensitivities.values().first(static_cast<std::size_t>(data_->retained.count));
}
std::optional<Context> ArnoldiSpectrum::context() const noexcept {
    if (!data_) {
        return std::nullopt;
    }
    return data_->context;
}
std::optional<StateDiagnosticContext> ArnoldiSpectrum::state() const noexcept {
    if (!data_) {
        return std::nullopt;
    }
    return data_->state;
}
std::optional<linear::ReturnedPrecision> ArnoldiSpectrum::returned_precision() const noexcept {
    if (!data_) {
        return std::nullopt;
    }
    return data_->returned;
}
std::optional<ArnoldiDiagnostics> ArnoldiSpectrum::diagnostics() const noexcept {
    if (!data_) {
        return std::nullopt;
    }
    return data_->diagnostics;
}
ArnoldiPhysicalCoverageUnresolved ArnoldiSpectrum::physical_coverage() const noexcept {
    return {};
}
std::optional<ArnoldiMeasurements> ArnoldiSpectrum::measurements() const noexcept {
    if (!data_) {
        return std::nullopt;
    }
    return data_->measurements;
}
std::optional<linear::BackendAdmission> ArnoldiSpectrum::backend_admission() const noexcept {
    if (!data_) {
        return std::nullopt;
    }
    return data_->admission;
}
std::size_t ArnoldiSpectrum::inverse_certificate_count() const noexcept {
    return data_ ? data_->inverse_count : 0;
}
std::optional<ArnoldiInverseAssessment>
ArnoldiSpectrum::inverse_assessment(core::Index index) const noexcept {
    if (!data_ || index < 0 || std::cmp_greater_equal(index, data_->inverse_count)) {
        return std::nullopt;
    }
    return data_->inverse_assessments.values()[static_cast<std::size_t>(index)];
}
std::optional<core::ErrorCertificate>
ArnoldiSpectrum::inverse_certificate(core::Index index) const noexcept {
    if (!data_ || index < 0 || std::cmp_greater_equal(index, data_->inverse_count)) {
        return std::nullopt;
    }
    const auto assessment = inverse_assessment(index);
    if (!assessment) {
        return std::nullopt;
    }
    return assessment->linear;
}
std::optional<linear::ConsumerCheck>
ArnoldiSpectrum::required_check(core::Index index) const noexcept {
    if (!data_ || index < 0 || std::cmp_greater_equal(index, data_->required_checks.size())) {
        return std::nullopt;
    }
    return data_->required_checks.values()[static_cast<std::size_t>(index)];
}
std::size_t ArnoldiSpectrum::consumer_certificate_count() const noexcept {
    return data_ ? data_->consumer_count : 0;
}
std::optional<core::ErrorCertificate>
ArnoldiSpectrum::consumer_certificate(core::Index index) const noexcept {
    if (!data_ || index < 0 || std::cmp_greater_equal(index, data_->consumer_count)) {
        return std::nullopt;
    }
    return data_->consumer_proofs
        .values()[static_cast<std::size_t>(index) + (data_->projected ? data_->consumer_count : 0)];
}

template <StorageScalar Scalar>
Result<ArnoldiSpectrum>
arnoldi_spectrum(const ActionProvider<Scalar>& provider, const ArnoldiScales& scales,
                 const ArnoldiRequest<Scalar>& request, const ArnoldiLimits& limits,
                 ArnoldiWorkspace& workspace, core::ResourcePlan& result_plan) noexcept {
    const Context entry_context = scales.context;
    const auto entry_state = request.state;
    return guard_entry<ArnoldiSpectrum>(
        entry_context, entry_state, {}, [&]() -> Result<ArnoldiSpectrum> {
            auto* scratch = detail::ArnoldiWorkspaceAccess::data(workspace);
            if (scratch == nullptr || scratch->busy.test_and_set(std::memory_order_acquire)) {
                return std::unexpected(complete_failure(
                    core::Failure{core::FailureCode::InvalidInput, core::Operation::ExtractEigen},
                    entry_context, entry_state));
            }
            struct Release {
                detail::ArnoldiWorkspaceData& data;
                ~Release() noexcept {
                    data.busy.clear(std::memory_order_release);
                }
            };
            const Release release{*scratch};
            auto result =
                detail::arnoldi_compute(provider, scales, request, limits, *scratch, result_plan);
            if (!result) {
                return std::unexpected(result.error());
            }
            return ArnoldiSpectrum{std::move(*result)};
        });
}
template Result<ArnoldiSpectrum> arnoldi_spectrum<double>(const ActionProvider<double>&,
                                                          const ArnoldiScales&,
                                                          const ArnoldiRequest<double>&,
                                                          const ArnoldiLimits&, ArnoldiWorkspace&,
                                                          core::ResourcePlan&) noexcept;
template Result<ArnoldiSpectrum> arnoldi_spectrum<std::complex<double>>(
    const ActionProvider<std::complex<double>>&, const ArnoldiScales&,
    const ArnoldiRequest<std::complex<double>>&, const ArnoldiLimits&, ArnoldiWorkspace&,
    core::ResourcePlan&) noexcept;

} // namespace uvlas::numerics::eigen
