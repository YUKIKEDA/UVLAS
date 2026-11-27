#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>
#include <uvlas/numerics/storage.hpp>

namespace uvlas::numerics {
namespace {
using core::Failure;
using core::FailureCode;
using core::Index;
using core::Operation;
using core::Result;

auto fail(FailureCode code, Operation operation) noexcept {
    return std::unexpected(Failure{code, operation});
}
Result<Index> length(std::size_t size, Operation operation) noexcept {
    if (!std::in_range<Index>(size)) {
        return fail(FailureCode::ResourceExceeded, operation);
    }
    return static_cast<Index>(size);
}
Result<Index> dense_extent(DenseShape shape, Operation operation) noexcept {
    if (shape.rows < 0 || shape.columns < 0 || shape.leading_dimension < shape.rows)
        return fail(FailureCode::InvalidInput, operation);
    const std::array<Index, 3> dimensions{shape.rows == 0 ? 0 : 1, shape.leading_dimension,
                                          shape.columns};
    return core::checked_extent_product(dimensions, operation);
}
bool valid_index(Index index, std::size_t size) noexcept {
    return index >= 0 && std::cmp_less(index, size);
}
bool finite(double value) noexcept {
    return std::isfinite(value);
}
bool finite(std::complex<double> value) noexcept {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}
template <class T>
Result<core::PlannedBuffer<T>> copy_buffer(std::span<const T> input, core::ResourcePlan& plan,
                                           Operation operation) noexcept {
    const auto count = length(input.size(), operation);
    if (!count) {
        return std::unexpected(count.error());
    }
    auto result = plan.allocate<T>(*count, operation);
    if (!result) {
        return std::unexpected(result.error());
    }
    std::copy(input.begin(), input.end(), result->values().begin());
    return result;
}
struct OrderedEntry {
    Index row = 0, column = 0, order = 0;
    double value = 0;
};
Result<Index> offset_count(Index rows, Index columns, Operation operation) noexcept {
    if (rows < 0 || columns < 0) {
        return fail(FailureCode::InvalidInput, operation);
    }
    if (columns == std::numeric_limits<Index>::max()) {
        return fail(FailureCode::ResourceExceeded, operation);
    }
    return columns + 1;
}
template <class T>
Result<void> check_finite(std::span<const T> values, Operation operation) noexcept {
    for (const auto& value : values) {
        if (!finite(value)) {
            return fail(FailureCode::Nonfinite, operation);
        }
    }
    return {};
}
} // namespace

template <StorageScalar Scalar>
DenseMatrix<Scalar>::DenseMatrix(DenseMatrix&& other) noexcept
    : shape_(std::exchange(other.shape_, DenseShape{0, 0, 0})), values_(std::move(other.values_)) {}
template <StorageScalar Scalar>
DenseMatrix<Scalar>& DenseMatrix<Scalar>::operator=(DenseMatrix&& other) noexcept {
    if (this != &other) {
        shape_ = std::exchange(other.shape_, DenseShape{0, 0, 0});
        values_ = std::move(other.values_);
    }
    return *this;
}
template <StorageScalar Scalar>
Result<DenseMatrix<Scalar>>
DenseMatrix<Scalar>::copy(DenseShape shape, std::span<const Scalar> values,
                          core::ResourcePlan& plan, Operation operation) noexcept {
    const auto elements = dense_extent(shape, operation);
    if (!elements) {
        return std::unexpected(elements.error());
    }
    if (!std::cmp_equal(*elements, values.size())) {
        return fail(FailureCode::InvalidInput, operation);
    }
    if (const auto checked = check_finite(values, operation); !checked) {
        return std::unexpected(checked.error());
    }
    auto owned = copy_buffer(values, plan, operation);
    if (!owned) {
        return std::unexpected(owned.error());
    }
    return DenseMatrix{shape, std::move(*owned)};
}

template <StorageScalar Scalar>
Result<DenseMatrix<Scalar>> DenseMatrix<Scalar>::zeroed(DenseShape shape, core::ResourcePlan& plan,
                                                        Operation operation) noexcept {
    const auto elements = dense_extent(shape, operation);
    if (!elements)
        return std::unexpected(elements.error());
    auto owned = plan.allocate<Scalar>(*elements, operation);
    if (!owned)
        return std::unexpected(owned.error());
    return DenseMatrix{shape, std::move(*owned)};
}
template <StorageScalar Scalar>
Result<void> DenseMatrix<Scalar>::validate(Operation operation) const noexcept {
    return check_finite<Scalar>(values_.values(), operation);
}
template class DenseMatrix<double>;
template class DenseMatrix<std::complex<double>>;

SparseReal::SparseReal(SparseReal&& other) noexcept
    : rows_(std::exchange(other.rows_, 0)),
      columns_(std::exchange(other.columns_, 0)),
      offsets_(std::move(other.offsets_)),
      indices_(std::move(other.indices_)),
      values_(std::move(other.values_)) {}
SparseReal& SparseReal::operator=(SparseReal&& other) noexcept {
    if (this != &other) {
        rows_ = std::exchange(other.rows_, 0);
        columns_ = std::exchange(other.columns_, 0);
        offsets_ = std::move(other.offsets_);
        indices_ = std::move(other.indices_);
        values_ = std::move(other.values_);
    }
    return *this;
}
std::span<const Index> SparseReal::column_offsets() const& noexcept {
    static constexpr Index zero = 0;
    if (columns_ == 0) {
        return {&zero, 1};
    }
    return offsets_.values();
}
Result<SparseReal> SparseReal::copy_csc(Index rows, Index columns, std::span<const Index> offsets,
                                        std::span<const Index> indices,
                                        std::span<const double> values, core::ResourcePlan& plan,
                                        Operation operation) noexcept {
    const auto count = offset_count(rows, columns, operation);
    if (!count) {
        return std::unexpected(count.error());
    }
    if (!std::cmp_equal(*count, offsets.size()) || indices.size() != values.size()) {
        return fail(FailureCode::InvalidInput, operation);
    }
    const auto nnz = length(values.size(), operation);
    if (!nnz) {
        return std::unexpected(nnz.error());
    }
    if (offsets.front() != 0 || offsets.back() != *nnz) {
        return fail(FailureCode::InvalidInput, operation);
    }
    for (Index column = 0; column < columns; ++column) {
        const auto c = static_cast<std::size_t>(column);
        if (offsets[c] < 0 || offsets[c] > offsets[c + 1] || offsets[c + 1] > *nnz) {
            return fail(FailureCode::InvalidInput, operation);
        }
        Index previous = -1;
        for (Index at = offsets[c]; at < offsets[c + 1]; ++at) {
            const Index row = indices[static_cast<std::size_t>(at)];
            if (row < 0 || row >= rows || row <= previous) {
                return fail(FailureCode::InvalidInput, operation);
            }
            previous = row;
        }
    }
    if (const auto checked = check_finite(values, operation); !checked) {
        return std::unexpected(checked.error());
    }
    SparseReal result;
    auto row_data = copy_buffer(indices, plan, operation);
    if (!row_data) {
        return std::unexpected(row_data.error());
    }
    auto value_data = copy_buffer(values, plan, operation);
    if (!value_data) {
        return std::unexpected(value_data.error());
    }
    if (columns != 0) {
        auto column_data = copy_buffer(offsets, plan, operation);
        if (!column_data) {
            return std::unexpected(column_data.error());
        }
        result.offsets_ = std::move(*column_data);
    }
    result.rows_ = rows;
    result.columns_ = columns;
    result.indices_ = std::move(*row_data);
    result.values_ = std::move(*value_data);
    return result;
}
Result<SparseReal> SparseReal::from_triplets(Index rows, Index columns,
                                             std::span<const Triplet> input,
                                             core::ResourcePlan& plan,
                                             Operation operation) noexcept {
    const auto count = offset_count(rows, columns, operation);
    if (!count) {
        return std::unexpected(count.error());
    }
    const auto n = length(input.size(), operation);
    if (!n) {
        return std::unexpected(n.error());
    }
    for (const Triplet& entry : input) {
        if (entry.row < 0 || entry.column < 0 || entry.row >= rows || entry.column >= columns) {
            return fail(FailureCode::InvalidInput, operation);
        }
        if (!finite(entry.value)) {
            return fail(FailureCode::Nonfinite, operation);
        }
    }
    auto sorted = plan.allocate<OrderedEntry>(*n, operation);
    if (!sorted) {
        return std::unexpected(sorted.error());
    }
    for (std::size_t i = 0; i < input.size(); ++i) {
        sorted->values()[i] = {input[i].row, input[i].column, static_cast<Index>(i),
                               input[i].value};
    }
    // The original order is part of the key: equal coordinates retain a fixed
    // reduction order even though std::sort is not stable (N4-02).
    std::sort(sorted->values().begin(), sorted->values().end(),
              [](const OrderedEntry& a, const OrderedEntry& b) {
                  return std::tie(a.column, a.row, a.order) < std::tie(b.column, b.row, b.order);
              });
    Index unique = 0;
    for (std::size_t i = 0; i < input.size();) {
        const Index row = sorted->values()[i].row, column = sorted->values()[i].column;
        double sum = 0, correction = 0;
        std::size_t end = i;
        while (end < input.size() && sorted->values()[end].row == row &&
               sorted->values()[end].column == column) {
            const double value = sorted->values()[end].value;
            const double next = sum + value;
            if (!finite(next)) {
                return fail(FailureCode::Nonfinite, operation);
            }
            correction +=
                std::abs(sum) >= std::abs(value) ? (sum - next) + value : (value - next) + sum;
            if (!finite(correction)) {
                return fail(FailureCode::Nonfinite, operation);
            }
            sum = next;
            ++end;
        }
        const double value = sum + correction;
        if (!finite(value)) {
            return fail(FailureCode::Nonfinite, operation);
        }
        sorted->values()[static_cast<std::size_t>(unique)] = {row, column, 0, value};
        ++unique;
        i = end;
    }
    SparseReal result;
    auto indices = plan.allocate<Index>(unique, operation);
    if (!indices) {
        return std::unexpected(indices.error());
    }
    auto values = plan.allocate<double>(unique, operation);
    if (!values) {
        return std::unexpected(values.error());
    }
    if (columns != 0) {
        auto offsets = plan.allocate<Index>(*count, operation);
        if (!offsets) {
            return std::unexpected(offsets.error());
        }
        Index cursor = 0;
        for (Index column = 0; column < columns; ++column) {
            offsets->values()[static_cast<std::size_t>(column)] = cursor;
            while (cursor < unique &&
                   sorted->values()[static_cast<std::size_t>(cursor)].column == column) {
                ++cursor;
            }
        }
        offsets->values()[static_cast<std::size_t>(columns)] = unique;
        result.offsets_ = std::move(*offsets);
    }
    for (Index i = 0; i < unique; ++i) {
        const auto at = static_cast<std::size_t>(i);
        indices->values()[at] = sorted->values()[at].row;
        values->values()[at] = sorted->values()[at].value;
    }
    result.rows_ = rows;
    result.columns_ = columns;
    result.indices_ = std::move(*indices);
    result.values_ = std::move(*values);
    return result;
}

Result<PointSoA> PointSoA::copy(std::span<const double> x, std::span<const double> y,
                                std::span<const double> z, std::span<const PointId> ids,
                                core::ResourcePlan& plan, Operation operation) noexcept {
    if (x.size() != y.size() || x.size() != z.size() || x.size() != ids.size()) {
        return fail(FailureCode::InvalidInput, operation);
    }
    for (const auto coordinates : {x, y, z}) {
        if (const auto checked = check_finite(coordinates, operation); !checked) {
            return std::unexpected(checked.error());
        }
    }
    const auto n = length(ids.size(), operation);
    if (!n) {
        return std::unexpected(n.error());
    }
    auto identity = plan.allocate<std::uint64_t>(*n, operation);
    if (!identity) {
        return std::unexpected(identity.error());
    }
    auto scratch = plan.allocate<std::uint64_t>(*n, operation);
    if (!scratch) {
        return std::unexpected(scratch.error());
    }
    for (std::size_t i = 0; i < ids.size(); ++i) {
        identity->values()[i] = ids[i].value();
        scratch->values()[i] = ids[i].value();
    }
    std::sort(scratch->values().begin(), scratch->values().end());
    if (std::adjacent_find(scratch->values().begin(), scratch->values().end()) !=
        scratch->values().end()) {
        return fail(FailureCode::InvalidInput, operation);
    }
    auto xs = copy_buffer(x, plan, operation);
    if (!xs) {
        return std::unexpected(xs.error());
    }
    auto ys = copy_buffer(y, plan, operation);
    if (!ys) {
        return std::unexpected(ys.error());
    }
    auto zs = copy_buffer(z, plan, operation);
    if (!zs) {
        return std::unexpected(zs.error());
    }
    PointSoA result;
    result.x_ = std::move(*xs);
    result.y_ = std::move(*ys);
    result.z_ = std::move(*zs);
    result.ids_ = std::move(*identity);
    return result;
}
Result<PointId> PointSoA::id(Index index, Operation operation) const noexcept {
    if (!valid_index(index, ids_.size())) {
        return fail(FailureCode::InvalidInput, operation);
    }
    return PointId{ids_.values()[static_cast<std::size_t>(index)]};
}

Result<DofLayout> DofLayout::copy(std::span<const DofDescriptor> full,
                                  std::span<const Index> independent_full,
                                  std::span<const DofAssociation> associations,
                                  core::ResourcePlan& plan, Operation operation) noexcept {
    const auto n = length(full.size(), operation), m = length(independent_full.size(), operation),
               a = length(associations.size(), operation);
    if (!n) {
        return std::unexpected(n.error());
    }
    if (!m) {
        return std::unexpected(m.error());
    }
    if (!a) {
        return std::unexpected(a.error());
    }
    for (const auto& dof : full) {
        if (!finite(dof.scale)) {
            return fail(FailureCode::Nonfinite, operation);
        }
        if (dof.scale <= 0 ||
            (dof.kind != CoordinateKind::Position && dof.kind != CoordinateKind::AncfGradient &&
             dof.kind != CoordinateKind::BeamRotation && dof.kind != CoordinateKind::OtherScalar)) {
            return fail(FailureCode::InvalidInput, operation);
        }
    }
    for (const Index index : independent_full) {
        if (!valid_index(index, full.size())) {
            return fail(FailureCode::InvalidInput, operation);
        }
    }
    for (const auto& owner : associations) {
        if (!valid_index(owner.full_index, full.size())) {
            return fail(FailureCode::InvalidInput, operation);
        }
    }
    auto identifiers = plan.allocate<std::uint64_t>(*n, operation);
    if (!identifiers) {
        return std::unexpected(identifiers.error());
    }
    for (std::size_t i = 0; i < full.size(); ++i) {
        identifiers->values()[i] = full[i].id.value();
    }
    std::sort(identifiers->values().begin(), identifiers->values().end());
    if (std::adjacent_find(identifiers->values().begin(), identifiers->values().end()) !=
        identifiers->values().end()) {
        return fail(FailureCode::InvalidInput, operation);
    }
    auto metadata = plan.allocate<StoredDof>(*n, operation);
    if (!metadata) {
        return std::unexpected(metadata.error());
    }
    auto independent = copy_buffer(independent_full, plan, operation);
    if (!independent) {
        return std::unexpected(independent.error());
    }
    auto inverse = plan.allocate<Index>(*n, operation);
    if (!inverse) {
        return std::unexpected(inverse.error());
    }
    std::fill(inverse->values().begin(), inverse->values().end(), Index{-1});
    for (std::size_t i = 0; i < independent_full.size(); ++i) {
        const auto at = static_cast<std::size_t>(independent_full[i]);
        if (inverse->values()[at] != -1) {
            return fail(FailureCode::InvalidInput, operation);
        }
        inverse->values()[at] = static_cast<Index>(i);
    }
    for (std::size_t i = 0; i < full.size(); ++i) {
        const auto& dof = full[i];
        metadata->values()[i] = {dof.id.value(), dof.component.value(), dof.kind, dof.unit,
                                 dof.scale};
    }
    auto owners = plan.allocate<StoredOwner>(*a, operation);
    if (!owners) {
        return std::unexpected(owners.error());
    }
    for (std::size_t i = 0; i < associations.size(); ++i) {
        const auto& owner = associations[i];
        if (const auto* patch = std::get_if<PatchId>(&owner.owner)) {
            owners->values()[i] = {owner.full_index, patch->value(), 0};
        } else if (const auto* member = std::get_if<MemberId>(&owner.owner)) {
            owners->values()[i] = {owner.full_index, member->value(), 1};
        } else {
            return fail(FailureCode::InvalidInput, operation);
        }
    }
    auto key = [](const StoredOwner& owner) {
        return std::tie(owner.full_index, owner.kind, owner.id);
    };
    std::sort(owners->values().begin(), owners->values().end(),
              [&](const auto& x, const auto& y) { return key(x) < key(y); });
    for (std::size_t i = 1; i < owners->size(); ++i) {
        if (key(owners->values()[i - 1]) == key(owners->values()[i])) {
            return fail(FailureCode::InvalidInput, operation);
        }
    }
    Index cursor = 0;
    for (Index i = 0; i < *n; ++i) {
        if (cursor == *a || owners->values()[static_cast<std::size_t>(cursor)].full_index != i) {
            return fail(FailureCode::InvalidInput, operation);
        }
        while (cursor < *a && owners->values()[static_cast<std::size_t>(cursor)].full_index == i) {
            ++cursor;
        }
    }
    DofLayout result;
    result.full_ = std::move(*metadata);
    result.independent_ = std::move(*independent);
    result.inverse_ = std::move(*inverse);
    result.owners_ = std::move(*owners);
    return result;
}
Result<DofDescriptor> DofLayout::descriptor(Index index, Operation operation) const noexcept {
    if (!valid_index(index, full_.size())) {
        return fail(FailureCode::InvalidInput, operation);
    }
    const auto& dof = full_.values()[static_cast<std::size_t>(index)];
    return DofDescriptor{core::DofId{dof.id}, core::ComponentId{dof.component}, dof.kind, dof.unit,
                         dof.scale};
}
Result<DofAssociation> DofLayout::association(Index index, Operation operation) const noexcept {
    if (!valid_index(index, owners_.size())) {
        return fail(FailureCode::InvalidInput, operation);
    }
    const auto& owner = owners_.values()[static_cast<std::size_t>(index)];
    const DofOwner id =
        owner.kind == 0 ? DofOwner{PatchId{owner.id}} : DofOwner{MemberId{owner.id}};
    return DofAssociation{owner.full_index, id};
}
Result<Index> DofLayout::independent_to_full(Index index, Operation operation) const noexcept {
    if (!valid_index(index, independent_.size())) {
        return fail(FailureCode::InvalidInput, operation);
    }
    return independent_.values()[static_cast<std::size_t>(index)];
}
Result<std::optional<Index>> DofLayout::full_to_independent(Index index,
                                                            Operation operation) const noexcept {
    if (!valid_index(index, inverse_.size())) {
        return fail(FailureCode::InvalidInput, operation);
    }
    const Index mapped = inverse_.values()[static_cast<std::size_t>(index)];
    return mapped == -1 ? std::optional<Index>{} : std::optional<Index>{mapped};
}
} // namespace uvlas::numerics
