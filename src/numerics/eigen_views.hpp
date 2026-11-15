#pragma once

// Private adapter: dependency types never appear in the product's public headers.
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <limits>
#include <type_traits>
#include <uvlas/numerics/storage.hpp>

namespace uvlas::numerics::detail {

template <StorageScalar Scalar>
using EigenMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <StorageScalar Scalar>
using EigenReadMap =
    Eigen::Map<const EigenMatrix<Scalar>, Eigen::Unaligned, Eigen::OuterStride<Eigen::Dynamic>>;

template <StorageScalar Scalar>
[[nodiscard]] core::Result<EigenReadMap<Scalar>>
eigen_read_map(const DenseMatrix<Scalar>& input, core::Index expected_rows,
               core::Index expected_columns, core::Operation operation) noexcept {
    const DenseShape shape = input.shape();
    if (shape.rows != expected_rows || shape.columns != expected_columns) {
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    }
    if (const auto check = input.validate(operation); !check) {
        return std::unexpected(check.error());
    }
    if (!std::in_range<Eigen::Index>(shape.rows) || !std::in_range<Eigen::Index>(shape.columns) ||
        !std::in_range<Eigen::Index>(shape.leading_dimension)) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, operation});
    }
    return EigenReadMap<Scalar>{
        input.values().data(), static_cast<Eigen::Index>(shape.rows),
        static_cast<Eigen::Index>(shape.columns),
        Eigen::OuterStride<Eigen::Dynamic>{static_cast<Eigen::Index>(shape.leading_dimension)}};
}
template <StorageScalar Scalar>
core::Result<EigenReadMap<Scalar>> eigen_read_map(const DenseMatrix<Scalar>&&, core::Index,
                                                  core::Index, core::Operation) = delete;

// Check the original object before Eigen::Ref could materialize a temporary.
// Explicit Scalar selects the requested product representation; no cast or copy is made.
template <StorageScalar Scalar, class Derived>
[[nodiscard]] core::Result<std::span<const Scalar>>
borrow_eigen_storage(const Derived& input, DenseShape expected,
                     core::Operation operation) noexcept {
    if constexpr (!std::is_same_v<typename Derived::Scalar, Scalar> || Derived::IsRowMajor ||
                  !(Derived::Flags & Eigen::DirectAccessBit)) {
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    } else {
        if (input.rows() != expected.rows || input.cols() != expected.columns ||
            input.innerStride() != 1 || input.outerStride() != expected.leading_dimension ||
            expected.rows < 0 || expected.columns < 0 ||
            expected.leading_dimension < expected.rows) {
            return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
        }
        // The returned span covers logical values only for packed storage. Padded
        // foreign storage has no verifiable allocation extent, so reject it rather
        // than claiming padding bytes that might lie beyond a block/map's buffer.
        if (expected.leading_dimension != expected.rows) {
            return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
        }
        const std::array<core::Index, 2> extents{expected.rows, expected.columns};
        const auto count = core::checked_extent_product(extents, operation);
        if (!count) {
            return std::unexpected(count.error());
        }
        const auto bytes = core::checked_byte_count(*count, sizeof(Scalar), operation);
        if (!bytes) {
            return std::unexpected(bytes.error());
        }
        const std::span<const Scalar> values{input.data(), *bytes / sizeof(Scalar)};
        for (const auto& value : values) {
            if constexpr (std::is_same_v<Scalar, double>) {
                if (!std::isfinite(value)) {
                    return std::unexpected(core::Failure{core::FailureCode::Nonfinite, operation});
                }
            } else {
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
                    return std::unexpected(core::Failure{core::FailureCode::Nonfinite, operation});
                }
            }
        }
        return values;
    }
}
template <StorageScalar Scalar, class Derived>
    requires(!std::is_lvalue_reference_v<Derived>)
core::Result<std::span<const Scalar>> borrow_eigen_storage(Derived&&, DenseShape,
                                                           core::Operation) = delete;

using EigenSparse = Eigen::SparseMatrix<double, Eigen::ColMajor, core::Index>;
using EigenSparseReadMap = Eigen::Map<const EigenSparse>;
[[nodiscard]] inline core::Result<EigenSparseReadMap>
eigen_read_map(const SparseReal& input, core::Index expected_rows, core::Index expected_columns,
               core::Operation operation) noexcept {
    if (input.rows() != expected_rows || input.columns() != expected_columns) {
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    }
    if (!std::in_range<Eigen::Index>(input.rows()) ||
        !std::in_range<Eigen::Index>(input.columns()) ||
        !std::in_range<Eigen::Index>(input.values().size())) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, operation});
    }
    return EigenSparseReadMap{static_cast<Eigen::Index>(input.rows()),
                              static_cast<Eigen::Index>(input.columns()),
                              static_cast<Eigen::Index>(input.values().size()),
                              input.column_offsets().data(),
                              input.row_indices().data(),
                              input.values().data()};
}
core::Result<EigenSparseReadMap> eigen_read_map(const SparseReal&&, core::Index, core::Index,
                                                core::Operation) = delete;
} // namespace uvlas::numerics::detail
