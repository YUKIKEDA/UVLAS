#pragma once

#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/Util/Version.h>

#include <complex>

namespace uvlas::numerics::eigen::detail {

// Audited against the manifest-pinned Spectra 1.2.0. Only protected factorization
// access is adapted; no private vendor implementation is copied or modified.
static_assert(EIGEN_MAJOR_VERSION == 5 && EIGEN_MINOR_VERSION == 0 && EIGEN_PATCH_VERSION == 1);
static_assert(SPECTRA_MAJOR_VERSION == 1 && SPECTRA_MINOR_VERSION == 2 &&
              SPECTRA_PATCH_VERSION == 0);
template <class Operator>
class ArnoldiNative final : public Spectra::GenEigsSolver<Operator> {
    using Base = Spectra::GenEigsSolver<Operator>;

public:
    ArnoldiNative(Operator& operation, Eigen::Index requested, Eigen::Index width)
        : Base(operation, requested, width) {}
    [[nodiscard]] Eigen::Index actual_dimension() const noexcept {
        return this->m_fac.subspace_dim();
    }
    [[nodiscard]] const auto& active_basis_storage() const noexcept {
        return this->m_fac.matrix_V();
    }
    [[nodiscard]] const auto& hessenberg_storage() const noexcept {
        return this->m_fac.matrix_H();
    }
};

} // namespace uvlas::numerics::eigen::detail
