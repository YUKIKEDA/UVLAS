#pragma once

#include <array>
#include <complex>
#include <concepts>
#include <cstdint>
#include <optional>
#include <span>
#include <uvlas/core/resource_plan.hpp>
#include <variant>

namespace uvlas::numerics {

/// Matrix dimensions and column stride, in elements rather than bytes.
struct DenseShape {
    core::Index rows;              ///< Logical row count.
    core::Index columns;           ///< Logical column count.
    core::Index leading_dimension; ///< Storage distance between successive columns.
};

/// Element types supported by owned binary64 matrix storage.
template <class Scalar>
concept StorageScalar = std::same_as<Scalar, double> || std::same_as<Scalar, std::complex<double>>;

/// Move-only column-major matrix, including padding. Writable access requires validate()
/// before external consumption. Moving invalidates views and leaves an empty source.
template <StorageScalar Scalar>
class DenseMatrix {
public:
    /// Construct an empty 0-by-0 matrix without allocation.
    DenseMatrix() noexcept = default;
    /// Transfer storage and leave the source empty.
    DenseMatrix(DenseMatrix&& other) noexcept;
    /// Transfer storage, invalidating previous destination views.
    DenseMatrix& operator=(DenseMatrix&& other) noexcept;
    /// Copying is forbidden because storage has one owner.
    DenseMatrix(const DenseMatrix&) = delete;
    /// Copy assignment is forbidden.
    DenseMatrix& operator=(const DenseMatrix&) = delete;
    /// Release owned storage.
    ~DenseMatrix() = default;
    /// Explicitly copy column-major storage; shape must have nonnegative dimensions
    /// and leading_dimension >= rows. Input has leading_dimension * columns elements,
    /// or zero for any empty logical shape. Padding must also be finite.
    /// @return Owned matrix, InvalidInput for shape/length mismatch, Nonfinite for
    /// nonfinite real/imaginary parts, or ResourceExceeded. Inputs live only during the call.
    [[nodiscard]] static core::Result<DenseMatrix> copy(DenseShape shape,
                                                        std::span<const Scalar> values,
                                                        core::ResourcePlan& plan,
                                                        core::Operation operation) noexcept;
    /// Allocate zero-initialized column-major storage, including padding.
    [[nodiscard]] static core::Result<DenseMatrix>
    zeroed(DenseShape shape, core::ResourcePlan& plan, core::Operation operation) noexcept;
    /// Return dimensions and column stride.
    [[nodiscard]] DenseShape shape() const noexcept {
        return shape_;
    }
    /// Borrow writable elements until assignment, moving from, or destruction.
    [[nodiscard]] std::span<Scalar> values() & noexcept {
        return values_.values();
    }
    /// Borrow read-only elements until assignment, moving from, or destruction.
    [[nodiscard]] std::span<const Scalar> values() const& noexcept {
        return values_.values();
    }
    /// A temporary cannot supply a borrowed view.
    [[nodiscard]] std::span<Scalar> values() && = delete;
    /// A const temporary cannot supply a borrowed view.
    [[nodiscard]] std::span<const Scalar> values() const&& = delete;
    /// Check finiteness after writable access, without copying or allocating.
    [[nodiscard]] core::Result<void> validate(core::Operation operation) const noexcept;

private:
    DenseMatrix(DenseShape shape, core::PlannedBuffer<Scalar> values) noexcept
        : shape_(shape), values_(std::move(values)) {}
    DenseShape shape_{0, 0, 0};
    core::PlannedBuffer<Scalar> values_;
};
/// Owned column-major binary64 matrix.
using DenseReal = DenseMatrix<double>;
/// Owned column-major complex-binary64 matrix.
using DenseComplex = DenseMatrix<std::complex<double>>;

/// One uncompressed contribution. Duplicate coordinates are legal input.
struct Triplet {
    core::Index row;    ///< Zero-based row.
    core::Index column; ///< Zero-based column.
    double value;       ///< Finite contribution.
};

/// Owns CSC storage with signed 64-bit indices and sorted unique rows per column.
/// Explicit zeros remain structural entries. Borrowed arrays expire on move/assignment/destruction.
class SparseReal {
public:
    /// Construct an empty sparse matrix without allocation.
    SparseReal() noexcept = default;
    /// Transfer arrays and leave the source empty.
    SparseReal(SparseReal&& other) noexcept;
    /// Replace arrays by moving.
    SparseReal& operator=(SparseReal&& other) noexcept;
    /// Copying is forbidden.
    SparseReal(const SparseReal&) = delete;
    /// Copy assignment is forbidden.
    SparseReal& operator=(const SparseReal&) = delete;
    /// Release owned arrays.
    ~SparseReal() = default;
    /// Aggregate duplicate contributions in input order with compensated summation.
    /// Sorting workspace is charged to plan. InvalidInput rejects shape/coordinates,
    /// Nonfinite rejects inputs or aggregation overflow, and ResourceExceeded rejects capacity.
    [[nodiscard]] static core::Result<SparseReal>
    from_triplets(core::Index rows, core::Index columns, std::span<const Triplet> input,
                  core::ResourcePlan& plan, core::Operation operation) noexcept;
    /// Explicitly copy canonical CSC. Offsets have columns+1 elements, start at zero,
    /// end at values.size(), and rows are sorted/unique. No implicit compression occurs.
    /// @return Owned matrix, InvalidInput for pattern mismatch, Nonfinite, or ResourceExceeded.
    [[nodiscard]] static core::Result<SparseReal>
    copy_csc(core::Index rows, core::Index columns, std::span<const core::Index> offsets,
             std::span<const core::Index> indices, std::span<const double> values,
             core::ResourcePlan& plan, core::Operation operation) noexcept;
    /// Logical row count.
    [[nodiscard]] core::Index rows() const noexcept {
        return rows_;
    }
    /// Logical column count.
    [[nodiscard]] core::Index columns() const noexcept {
        return columns_;
    }
    /// Borrow compressed column offsets.
    [[nodiscard]] std::span<const core::Index> column_offsets() const& noexcept;
    /// Borrow compressed row indices.
    [[nodiscard]] std::span<const core::Index> row_indices() const& noexcept {
        return indices_.values();
    }
    /// Borrow numerical values.
    [[nodiscard]] std::span<const double> values() const& noexcept {
        return values_.values();
    }
    /// A temporary cannot supply borrowed arrays.
    [[nodiscard]] std::span<const core::Index> column_offsets() const&& = delete;
    /// A temporary cannot supply borrowed arrays.
    [[nodiscard]] std::span<const core::Index> row_indices() const&& = delete;
    /// A temporary cannot supply borrowed arrays.
    [[nodiscard]] std::span<const double> values() const&& = delete;

private:
    core::Index rows_ = 0, columns_ = 0;
    core::PlannedBuffer<core::Index> offsets_, indices_;
    core::PlannedBuffer<double> values_;
};

namespace detail {
struct PointTag;
struct PatchTag;
struct MemberTag;
} // namespace detail
/// Physical point identity, distinct from an array offset.
using PointId = core::Identifier<detail::PointTag>;
/// Physical surface-patch identity.
using PatchId = core::Identifier<detail::PatchTag>;
/// Physical structural-member identity.
using MemberId = core::Identifier<detail::MemberTag>;

/// Owns equal-length finite x/y/z arrays in metres and unique point identities.
class PointSoA {
public:
    /// Construct an empty point set.
    PointSoA() noexcept = default;
    /// Transfer equal-length arrays; leave the source empty.
    PointSoA(PointSoA&&) noexcept = default;
    /// Replace arrays by moving; previous destination views expire.
    PointSoA& operator=(PointSoA&&) noexcept = default;
    /// Copying is forbidden.
    PointSoA(const PointSoA&) = delete;
    /// Copy assignment is forbidden.
    PointSoA& operator=(const PointSoA&) = delete;
    /// Release arrays.
    ~PointSoA() = default;
    /// Explicitly copy coordinates and identities in input order.
    /// @return InvalidInput for unequal lengths/duplicate IDs, Nonfinite for coordinates,
    /// ResourceExceeded for capacity, or owned points. Input is borrowed only during the call.
    [[nodiscard]] static core::Result<PointSoA>
    copy(std::span<const double> x, std::span<const double> y, std::span<const double> z,
         std::span<const PointId> ids, core::ResourcePlan& plan,
         core::Operation operation) noexcept;
    /// Point count.
    [[nodiscard]] std::size_t size() const noexcept {
        return x_.size();
    }
    /// Borrow x coordinates until assignment, moving from, or destruction.
    [[nodiscard]] std::span<const double> x() const& noexcept {
        return x_.values();
    }
    /// Borrow y coordinates until assignment, moving from, or destruction.
    [[nodiscard]] std::span<const double> y() const& noexcept {
        return y_.values();
    }
    /// Borrow z coordinates until assignment, moving from, or destruction.
    [[nodiscard]] std::span<const double> z() const& noexcept {
        return z_.values();
    }
    /// Return an identity by value, or InvalidInput for an invalid index.
    [[nodiscard]] core::Result<PointId> id(core::Index index,
                                           core::Operation operation) const noexcept;
    /// A temporary cannot supply a borrowed coordinate array.
    [[nodiscard]] std::span<const double> x() const&& = delete;
    /// A temporary cannot supply a borrowed coordinate array.
    [[nodiscard]] std::span<const double> y() const&& = delete;
    /// A temporary cannot supply a borrowed coordinate array.
    [[nodiscard]] std::span<const double> z() const&& = delete;

private:
    core::PlannedBuffer<double> x_, y_, z_;
    core::PlannedBuffer<std::uint64_t> ids_;
};

/// SI exponents in metre, kilogram, second, ampere, kelvin, mole, candela order.
struct SiDimension {
    std::array<std::int8_t, 7> powers; ///< Signed base-unit exponents.
};
/// Distinguish position, ANCF gradients and beam rotation-chart coordinates.
enum class CoordinateKind : std::uint8_t { Position, AncfGradient, BeamRotation, OtherScalar };
/// Metadata for a full coordinate; array order defines component order.
struct DofDescriptor {
    core::DofId id;              ///< Unique physical degree-of-freedom identity.
    core::ComponentId component; ///< Component identity defined by the model.
    CoordinateKind kind;         ///< Coordinate interpretation, distinct from its unit.
    SiDimension unit;            ///< SI dimension without an implicit unit conversion.
    double scale;                ///< Positive finite physical scale in the declared unit.
};
/// Patch or member associated with a coordinate; several owners are legal.
using DofOwner = std::variant<PatchId, MemberId>;
/// Associate a full coordinate index with a patch or member.
struct DofAssociation {
    core::Index full_index = 0; ///< Full-coordinate index.
    DofOwner owner;             ///< Physical patch/member identity.
};

/// Owns coordinate metadata, topology and independent-coordinate index selection.
/// Dependent coordinates are recovered by the model's nonlinear constraint map, not
/// set to zero. These index maps do not evaluate a physical coordinate transformation.
class DofLayout {
public:
    /// Construct an empty layout.
    DofLayout() noexcept = default;
    /// Transfer all arrays and leave the source empty.
    DofLayout(DofLayout&&) noexcept = default;
    /// Replace metadata by moving.
    DofLayout& operator=(DofLayout&&) noexcept = default;
    /// Copying is forbidden.
    DofLayout(const DofLayout&) = delete;
    /// Copy assignment is forbidden.
    DofLayout& operator=(const DofLayout&) = delete;
    /// Release arrays.
    ~DofLayout() = default;
    /// Explicitly copy metadata. independent_full has distinct full indices in
    /// independent-coordinate order. Every full coordinate needs an owner association.
    /// @return InvalidInput for invalid indices, duplicate IDs/associations, unknown kind,
    /// missing owners or nonpositive scales; Nonfinite for scales; ResourceExceeded otherwise.
    [[nodiscard]] static core::Result<DofLayout> copy(std::span<const DofDescriptor> full,
                                                      std::span<const core::Index> independent_full,
                                                      std::span<const DofAssociation> associations,
                                                      core::ResourcePlan& plan,
                                                      core::Operation operation) noexcept;
    /// Full-coordinate count.
    [[nodiscard]] std::size_t full_size() const noexcept {
        return full_.size();
    }
    /// Independent-coordinate count.
    [[nodiscard]] std::size_t independent_size() const noexcept {
        return independent_.size();
    }
    /// Unique association count, sorted by full index, owner kind and identity.
    [[nodiscard]] std::size_t association_size() const noexcept {
        return owners_.size();
    }
    /// Return metadata by value, or InvalidInput for an invalid index.
    [[nodiscard]] core::Result<DofDescriptor> descriptor(core::Index index,
                                                         core::Operation operation) const noexcept;
    /// Return an association by value, or InvalidInput for an invalid index.
    [[nodiscard]] core::Result<DofAssociation>
    association(core::Index index, core::Operation operation) const noexcept;
    /// Return the selected full index, or InvalidInput for an invalid independent index.
    [[nodiscard]] core::Result<core::Index>
    independent_to_full(core::Index index, core::Operation operation) const noexcept;
    /// Return an independent index, nullopt for a dependent coordinate, or InvalidInput.
    [[nodiscard]] core::Result<std::optional<core::Index>>
    full_to_independent(core::Index index, core::Operation operation) const noexcept;

private:
    struct StoredDof {
        std::uint64_t id = 0, component = 0;
        CoordinateKind kind = CoordinateKind::OtherScalar;
        SiDimension unit{};
        double scale = 1;
    };
    struct StoredOwner {
        core::Index full_index = 0;
        std::uint64_t id = 0;
        std::uint8_t kind = 0;
    };
    core::PlannedBuffer<StoredDof> full_;
    core::PlannedBuffer<core::Index> independent_, inverse_;
    core::PlannedBuffer<StoredOwner> owners_;
};

} // namespace uvlas::numerics
