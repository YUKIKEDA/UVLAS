#include "backend_memory_detail.hpp"
#include "compensated.hpp"
#include "eigen_views.hpp"
#include "linear_solve_detail.hpp"

#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <type_traits>
#include <uvlas/numerics/dense_solve.hpp>

namespace uvlas::numerics::linear {
namespace {
using core::FailureCode;
using core::Index;
using core::Result;
using Complex = std::complex<double>;
using RealMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using ComplexMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using RealLu = Eigen::PartialPivLU<RealMatrix, Index>;
using ComplexLu = Eigen::PartialPivLU<ComplexMatrix, Index>;
constexpr auto operation = core::Operation::SolveLinear;
auto fail(FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}
template <class Value>
Result<Value> with_context(Result<Value> result, const Context& context) noexcept {
    if (!result) {
        auto& details = result.error().details;
        // An inner failure may identify a more specific evaluation or resource request.
        // Complete only unknown identities at the public boundary; never relabel known ones.
        if (!details.model)
            details.model = context.cache.model;
        if (!details.trial)
            details.trial = context.trial;
        if (!details.numerical_generation)
            details.numerical_generation = context.cache.matrix;
    }
    return result;
}
bool finite(double value) noexcept {
    return std::isfinite(value);
}
bool finite(Complex value) noexcept {
    return finite(value.real()) && finite(value.imag());
}
template <StorageScalar Scalar>
const DenseMatrix<Scalar>& high(const PrecisionBuffer<Scalar>& values) noexcept {
    if (const auto* binary = std::get_if<Binary64<Scalar>>(&values.storage())) {
        return binary->values;
    }
    return std::get_if<TwoWord<Scalar>>(&values.storage())->hi;
}
template <StorageScalar Scalar>
const DenseMatrix<Scalar>* low(const PrecisionBuffer<Scalar>& values) noexcept {
    const auto* pair = std::get_if<TwoWord<Scalar>>(&values.storage());
    return pair == nullptr ? nullptr : &pair->lo;
}
template <StorageScalar Scalar>
Scalar entry(const DenseMatrix<Scalar>& matrix, Index row, Index column) noexcept {
    return matrix
        .values()[static_cast<std::size_t>(column * matrix.shape().leading_dimension + row)];
}
Result<Index> count(Index rows, Index columns) noexcept {
    return core::checked_extent_product(std::array<Index, 2>{rows, columns}, operation);
}
Result<std::uint64_t> sum_bytes(std::uint64_t a, std::uint64_t b) noexcept {
    if (b > std::numeric_limits<std::uint64_t>::max() - a) {
        return fail(FailureCode::ResourceExceeded);
    }
    return a + b;
}
Result<BackendMemoryPlan> backend_plan(Index n, Index rank) noexcept {
    if (n == 0) {
        return detail::BackendPlanFactory::make(0, 0);
    }
    const auto square = count(n, n), small = count(rank, rank);
    if (!square) {
        return std::unexpected(square.error());
    }
    if (!small) {
        return std::unexpected(small.error());
    }
    // Fixed Eigen 5 forecast: retained LU/permutations plus concurrent packing,
    // small-Schur storage and estimator work. This is deliberately not a proof.
    const auto base = core::checked_byte_count(*square, 24, operation);
    const auto rows = core::checked_byte_count(n, 64, operation);
    const auto schur = core::checked_byte_count(*small, 48, operation);
    const auto update = core::checked_byte_count(rank, 128, operation);
    if (!base) {
        return std::unexpected(base.error());
    }
    if (!rows) {
        return std::unexpected(rows.error());
    }
    if (!schur) {
        return std::unexpected(schur.error());
    }
    if (!update) {
        return std::unexpected(update.error());
    }
    auto total = sum_bytes(*base, *rows);
    if (!total) {
        return std::unexpected(total.error());
    }
    total = sum_bytes(*total, *schur);
    if (!total) {
        return std::unexpected(total.error());
    }
    total = sum_bytes(*total, *update);
    if (!total) {
        return std::unexpected(total.error());
    }
    return detail::BackendPlanFactory::make(*total);
}
template <StorageScalar Scalar>
bool same_matrix(const DenseMatrix<Scalar>& a, const DenseMatrix<Scalar>& b) noexcept {
    if (a.shape().rows != b.shape().rows || a.shape().columns != b.shape().columns) {
        return false;
    }
    for (Index column = 0; column < a.shape().columns; ++column) {
        for (Index row = 0; row < a.shape().rows; ++row) {
            if constexpr (std::same_as<Scalar, double>) {
                if (std::bit_cast<std::uint64_t>(entry(a, row, column)) !=
                    std::bit_cast<std::uint64_t>(entry(b, row, column))) {
                    return false;
                }
            } else {
                const auto first = entry(a, row, column), second = entry(b, row, column);
                if (std::bit_cast<std::uint64_t>(first.real()) !=
                        std::bit_cast<std::uint64_t>(second.real()) ||
                    std::bit_cast<std::uint64_t>(first.imag()) !=
                        std::bit_cast<std::uint64_t>(second.imag())) {
                    return false;
                }
            }
        }
    }
    return true;
}
bool same_scales(const Scales& scales, const detail::MetadataBuffers& snapshot) noexcept {
    if (scales.unknown.size() != snapshot.unknown_scales.size()) {
        return false;
    }
    for (std::size_t i = 0; i < scales.unknown.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(scales.unknown[i]) !=
                std::bit_cast<std::uint64_t>(snapshot.unknown_scales.values()[i]) ||
            std::bit_cast<std::uint64_t>(scales.residual[i]) !=
                std::bit_cast<std::uint64_t>(snapshot.residual_scales.values()[i]) ||
            scales.unknown_units[i] != snapshot.unknown_units.values()[i] ||
            scales.residual_units[i] != snapshot.residual_units.values()[i]) {
            return false;
        }
    }
    return true;
}
template <class Factor>
Result<void> check_factor(const Factor& factor) noexcept {
    if (!factor.matrixLU().allFinite()) {
        return fail(FailureCode::Nonfinite);
    }
    for (Index i = 0; i < factor.rows(); ++i) {
        if (std::abs(factor.matrixLU()(i, i)) == 0) {
            return fail(FailureCode::RankUnresolved);
        }
    }
    return {};
}
template <class Factor, StorageScalar Scalar>
Result<void> triangular_solve(const Factor& factor, std::span<const Scalar> rhs,
                              std::span<Scalar> destination) noexcept {
    const Index n = static_cast<Index>(rhs.size());
    if (n == 0) {
        return {};
    }
    const typename Factor::MatrixType& lu = factor.matrixLU();
    const typename Factor::PermutationType::IndicesType& permutation =
        factor.permutationP().indices();
    // Eigen stores P's destination row per source row. The input/output do not alias.
    for (Index i = 0; i < n; ++i) {
        destination[static_cast<std::size_t>(permutation(i))] = rhs[static_cast<std::size_t>(i)];
    }
    for (Index i = 0; i < n; ++i) {
        Scalar value = destination[static_cast<std::size_t>(i)];
        for (Index j = 0; j < i; ++j) {
            value -= lu(i, j) * destination[static_cast<std::size_t>(j)];
        }
        if (!finite(value)) {
            return fail(FailureCode::Nonfinite);
        }
        destination[static_cast<std::size_t>(i)] = value;
    }
    for (Index end = n; end > 0;) {
        const Index i = --end;
        Scalar value = destination[static_cast<std::size_t>(i)];
        for (Index j = i + 1; j < n; ++j) {
            value -= lu(i, j) * destination[static_cast<std::size_t>(j)];
        }
        value /= lu(i, i);
        if (!finite(value)) {
            return fail(FailureCode::Nonfinite);
        }
        destination[static_cast<std::size_t>(i)] = value;
    }
    return {};
}
} // namespace
namespace detail {
template <StorageScalar Scalar>
struct CommonData {
    CommonData(Context identity, Limits limit, BackendAdmission admitted) noexcept
        : context(identity),
          capacity(limit),
          admission(admitted),
          lock(core::PlannedBuffer<std::byte>{}) {}
    Context context;
    Limits capacity;
    BackendAdmission admission;
    MetadataBuffers scales;
    core::Workspace lock;
    core::PlannedBuffer<Scalar> rhs_hi, rhs_lo, residual, correction, proposal;
    core::PlannedBuffer<double> residual_error, denominator_lower;
    core::PlannedBuffer<MagnitudeBounds> x_magnitude;
};
struct DenseData : CommonData<double> {
    using CommonData<double>::CommonData;
    DenseReal original;
    RealLu factor;
    ConditionDiagnostic condition = NotEvaluated{};
};
struct UpdatedData : CommonData<Complex> {
    using CommonData<Complex>::CommonData;
    DenseReal original;
    DenseComplex u, v, y;
    RealLu factor;
    ComplexLu schur;
    ConditionDiagnostic schur_diagnostic = NotEvaluated{};
    core::PlannedBuffer<Complex> small_rhs, small_solution;
};
} // namespace detail
namespace {
template <StorageScalar Scalar>
Result<void> prepare_arrays(detail::CommonData<Scalar>& data, Index n, const Scales& scales,
                            core::ResourcePlan& plan) noexcept {
    auto metadata = detail::SolutionFactory::prepare_metadata(n, 0, scales, {}, plan);
    if (!metadata) {
        return std::unexpected(metadata.error());
    }
    data.scales = std::move(*metadata);
    const auto size = count(n, data.capacity.maximum_rhs_columns);
    if (!size) {
        return std::unexpected(size.error());
    }
    for (auto* buffer :
         {&data.rhs_hi, &data.rhs_lo, &data.residual, &data.correction, &data.proposal}) {
        auto storage = plan.allocate<Scalar>(*size, operation);
        if (!storage) {
            return std::unexpected(storage.error());
        }
        *buffer = std::move(*storage);
    }
    for (auto* buffer : {&data.residual_error, &data.denominator_lower}) {
        auto storage = plan.allocate<double>(*size, operation);
        if (!storage) {
            return std::unexpected(storage.error());
        }
        *buffer = std::move(*storage);
    }
    for (auto* buffer : {&data.x_magnitude}) {
        auto storage = plan.allocate<detail::MagnitudeBounds>(*size, operation);
        if (!storage) {
            return std::unexpected(storage.error());
        }
        *buffer = std::move(*storage);
    }
    return {};
}
Result<void> prepare_input(const DenseReal& matrix, Context context, const Scales& scales,
                           const Limits& limits, Index rank) noexcept {
    const auto valid = matrix.validate(operation);
    if (!valid) {
        return valid;
    }
    if (matrix.shape().rows != matrix.shape().columns) {
        return fail(FailureCode::InvalidInput);
    }
    const auto identity = detail::check_context_inputs(context);
    if (!identity) {
        return identity;
    }
    if (context.cache.representation != core::Precision::Binary64) {
        return fail(FailureCode::PrecisionUnsupported);
    }
    const auto units = detail::check_scale_inputs(matrix.shape().rows, scales, context);
    if (!units) {
        return units;
    }
    return detail::check_capacity_inputs(matrix.shape().rows, limits.maximum_rhs_columns, rank,
                                         limits);
}
} // namespace
DensePrepared::DensePrepared(std::unique_ptr<detail::DenseData> data) noexcept
    : data_(std::move(data)) {}
DensePrepared::DensePrepared(DensePrepared&&) noexcept = default;
DensePrepared& DensePrepared::operator=(DensePrepared&&) noexcept = default;
DensePrepared::~DensePrepared() = default;
std::optional<BackendMemoryPlan> DensePrepared::memory_plan() const noexcept {
    return data_ ? std::optional<BackendMemoryPlan>{data_->admission.plan()} : std::nullopt;
}
Result<DensePrepared> DensePrepared::prepare(const DenseOperator& op, const Scales& scales,
                                             const Limits& limits,
                                             core::ResourcePlan& plan) noexcept {
    auto result = core::guard_entry(
        operation, FailureCode::FactorizationFailed, {}, [&]() -> Result<DensePrepared> {
            const auto valid = prepare_input(op.values, op.context, scales, limits, 0);
            if (!valid) {
                return std::unexpected(valid.error());
            }
            const Index n = op.values.shape().rows;
            const auto memory = backend_plan(n, 0);
            if (!memory) {
                return std::unexpected(memory.error());
            }
            const auto admitted = admit_backend_memory(*memory, limits.backend_memory);
            if (!admitted) {
                return std::unexpected(admitted.error());
            }
            auto data = std::make_unique<detail::DenseData>(op.context, limits, *admitted);
            auto original = DenseReal::copy(op.values.shape(), op.values.values(), plan, operation);
            if (!original) {
                return std::unexpected(original.error());
            }
            data->original = std::move(*original);
            const auto arrays = prepare_arrays(*data, n, scales, plan);
            if (!arrays) {
                return std::unexpected(arrays.error());
            }
            if (n > 0) {
                const auto map =
                    uvlas::numerics::detail::eigen_read_map(data->original, n, n, operation);
                if (!map) {
                    return std::unexpected(map.error());
                }
                data->factor.compute(*map);
                const auto factor = check_factor(data->factor);
                if (!factor) {
                    return std::unexpected(factor.error());
                }
                // Diagnostic for raw original A only. A failed/zero estimate never decides rank.
                const double estimate = data->factor.rcond();
                if (finite(estimate) && estimate >= 0) {
                    data->condition = ReciprocalOneNormEstimate{estimate};
                }
            }
            return DensePrepared{std::move(data)};
        });
    return with_context(std::move(result), op.context);
}
UpdatedPrepared::UpdatedPrepared(std::unique_ptr<detail::UpdatedData> data) noexcept
    : data_(std::move(data)) {}
UpdatedPrepared::UpdatedPrepared(UpdatedPrepared&&) noexcept = default;
UpdatedPrepared& UpdatedPrepared::operator=(UpdatedPrepared&&) noexcept = default;
UpdatedPrepared::~UpdatedPrepared() = default;
std::optional<BackendMemoryPlan> UpdatedPrepared::memory_plan() const noexcept {
    return data_ ? std::optional<BackendMemoryPlan>{data_->admission.plan()} : std::nullopt;
}
ConditionDiagnostic UpdatedPrepared::schur_condition() const noexcept {
    return data_ ? data_->schur_diagnostic : ConditionDiagnostic{NotEvaluated{}};
}
Result<UpdatedPrepared> UpdatedPrepared::prepare(const ColumnUpdatedOperator& op,
                                                 const Scales& scales, const Limits& limits,
                                                 core::ResourcePlan& plan) noexcept {
    auto result = core::guard_entry(
        operation, FailureCode::FactorizationFailed, {}, [&]() -> Result<UpdatedPrepared> {
            const auto base =
                prepare_input(op.base, op.context, scales, limits, op.u.shape().columns);
            if (!base) {
                return std::unexpected(base.error());
            }
            const auto valid_u = op.u.validate(operation), valid_v = op.v.validate(operation);
            if (!valid_u) {
                return std::unexpected(valid_u.error());
            }
            if (!valid_v) {
                return std::unexpected(valid_v.error());
            }
            const Index n = op.base.shape().rows, rank = op.u.shape().columns;
            if (op.u.shape().rows != n || op.v.shape().rows != n || op.v.shape().columns != rank) {
                return fail(FailureCode::InvalidInput);
            }
            const auto memory = backend_plan(n, rank);
            if (!memory) {
                return std::unexpected(memory.error());
            }
            const auto admitted = admit_backend_memory(*memory, limits.backend_memory);
            if (!admitted) {
                return std::unexpected(admitted.error());
            }
            auto data = std::make_unique<detail::UpdatedData>(op.context, limits, *admitted);
            auto original = DenseReal::copy(op.base.shape(), op.base.values(), plan, operation);
            auto u = DenseComplex::copy(op.u.shape(), op.u.values(), plan, operation);
            auto v = DenseComplex::copy(op.v.shape(), op.v.values(), plan, operation);
            if (!original) {
                return std::unexpected(original.error());
            }
            if (!u) {
                return std::unexpected(u.error());
            }
            if (!v) {
                return std::unexpected(v.error());
            }
            data->original = std::move(*original);
            data->u = std::move(*u);
            data->v = std::move(*v);
            const auto arrays = prepare_arrays(*data, n, scales, plan);
            if (!arrays) {
                return std::unexpected(arrays.error());
            }
            if (n > 0) {
                const auto map =
                    uvlas::numerics::detail::eigen_read_map(data->original, n, n, operation);
                if (!map) {
                    return std::unexpected(map.error());
                }
                data->factor.compute(*map);
                const auto factor = check_factor(data->factor);
                if (!factor) {
                    return std::unexpected(factor.error());
                }
                auto y = DenseComplex::copy(data->u.shape(), data->u.values(), plan, operation);
                if (!y) {
                    return std::unexpected(y.error());
                }
                data->y = std::move(*y);
                for (Index column = 0; column < rank; ++column) {
                    const auto input = data->u.values().subspan(
                        static_cast<std::size_t>(column * data->u.shape().leading_dimension),
                        static_cast<std::size_t>(n));
                    auto output = data->y.values().subspan(
                        static_cast<std::size_t>(column * data->y.shape().leading_dimension),
                        static_cast<std::size_t>(n));
                    const auto solved =
                        triangular_solve<RealLu, Complex>(data->factor, input, output);
                    if (!solved) {
                        return std::unexpected(solved.error());
                    }
                }
                if (rank > 0) {
                    const auto size = count(rank, rank);
                    if (!size) {
                        return std::unexpected(size.error());
                    }
                    auto coefficients = plan.allocate<Complex>(*size, operation);
                    if (!coefficients) {
                        return std::unexpected(coefficients.error());
                    }
                    for (Index column = 0; column < rank; ++column) {
                        for (Index row = 0; row < rank; ++row) {
                            Complex value = row == column ? Complex{1, 0} : Complex{};
                            for (Index i = 0; i < n; ++i) {
                                value += entry(data->v, i, row) * entry(data->y, i, column);
                            }
                            if (!finite(value)) {
                                return fail(FailureCode::Nonfinite);
                            }
                            coefficients->values()[static_cast<std::size_t>(column * rank + row)] =
                                value;
                        }
                    }
                    const Eigen::Map<const ComplexMatrix> small{coefficients->values().data(), rank,
                                                                rank};
                    data->schur.compute(small);
                    const auto factor_small = check_factor(data->schur);
                    if (!factor_small) {
                        return std::unexpected(factor_small.error());
                    }
                    auto first = plan.allocate<Complex>(rank, operation),
                         second = plan.allocate<Complex>(rank, operation);
                    if (!first) {
                        return std::unexpected(first.error());
                    }
                    if (!second) {
                        return std::unexpected(second.error());
                    }
                    const double reciprocal = data->schur.rcond();
                    if (std::isfinite(reciprocal) && reciprocal >= 0)
                        data->schur_diagnostic = ReciprocalOneNormEstimate{reciprocal};
                    data->small_rhs = std::move(*first);
                    data->small_solution = std::move(*second);
                }
            }
            return UpdatedPrepared{std::move(data)};
        });
    return with_context(std::move(result), op.context);
}
} // namespace uvlas::numerics::linear

namespace uvlas::numerics::linear {
namespace {
double coefficient(const detail::DenseData& data, Index row, Index col) noexcept {
    return entry(data.original, row, col);
}
Complex coefficient(const detail::UpdatedData& data, Index row, Index col) noexcept {
    Complex value = entry(data.original, row, col);
    // This rounded coefficient, in increasing rank order, defines the original operator.
    for (Index k = 0; k < data.u.shape().columns; ++k) {
        const Complex u = entry(data.u, row, k), v = entry(data.v, col, k);
        const double rr = u.real() * v.real(), ii = u.imag() * v.imag();
        const double ri = u.real() * v.imag(), ir = u.imag() * v.real();
        const Complex product{rr - ii, ri + ir};
        value = {value.real() + product.real(), value.imag() + product.imag()};
    }
    return value;
}
Result<void> inverse(detail::DenseData& data, std::span<const double> rhs,
                     std::span<double> out) noexcept {
    return triangular_solve<RealLu, double>(data.factor, rhs, out);
}
Result<void> inverse(detail::UpdatedData& data, std::span<const Complex> rhs,
                     std::span<Complex> out) noexcept {
    auto base = triangular_solve<RealLu, Complex>(data.factor, rhs, out);
    if (!base)
        return base;
    if (rhs.empty())
        return {};
    const Index rank = data.u.shape().columns;
    for (Index k = 0; k < rank; ++k) {
        Complex value{};
        for (Index i = 0; i < static_cast<Index>(rhs.size()); ++i)
            value += entry(data.v, i, k) * out[static_cast<std::size_t>(i)];
        if (!finite(value))
            return fail(FailureCode::Nonfinite);
        data.small_rhs.values()[static_cast<std::size_t>(k)] = value;
    }
    auto small = triangular_solve<ComplexLu, Complex>(data.schur, data.small_rhs.values(),
                                                      data.small_solution.values());
    if (!small)
        return small;
    for (Index i = 0; i < static_cast<Index>(rhs.size()); ++i) {
        for (Index k = 0; k < rank; ++k)
            out[static_cast<std::size_t>(i)] -=
                entry(data.y, i, k) * data.small_solution.values()[static_cast<std::size_t>(k)];
        if (!finite(out[static_cast<std::size_t>(i)]))
            return fail(FailureCode::Nonfinite);
    }
    return {};
}
using Accumulator = uvlas::numerics::detail::CompensatedAccumulator;
Result<void> product(Accumulator& real, Accumulator& imag, double a, double x) noexcept {
    (void)imag;
    return real.product(a, x, -1);
}
Result<void> product(Accumulator& real, Accumulator& imag, Complex a, Complex x) noexcept {
    auto rr = real.product(a.real(), x.real(), -1);
    if (!rr)
        return rr;
    auto ii = real.product(a.imag(), x.imag(), 1);
    if (!ii)
        return ii;
    auto ri = imag.product(a.real(), x.imag(), -1);
    if (!ri)
        return ri;
    return imag.product(a.imag(), x.real(), -1);
}
template <StorageScalar Scalar>
Result<void> add_rhs(Accumulator& real, Accumulator& imag, Scalar value) noexcept {
    if constexpr (std::same_as<Scalar, double>) {
        (void)imag;
        return real.add(value);
    } else {
        auto r = real.add(value.real());
        if (!r)
            return r;
        return imag.add(value.imag());
    }
}
template <StorageScalar Scalar, class Data>
Result<bool> evaluate(Data& data, const DenseMatrix<Scalar>& hi, const DenseMatrix<Scalar>* lo,
                      const Scales& scales, const Accuracy& accuracy,
                      std::span<ColumnMetrics> metrics) noexcept {
    const Index n = hi.shape().rows, columns = hi.shape().columns;
    bool accepted = true;
    for (Index c = 0; c < columns; ++c) {
        const auto offset = static_cast<std::size_t>(c * n), size = static_cast<std::size_t>(n);
        for (Index j = 0; j < n; ++j) {
            auto x =
                detail::actual_magnitude_bounds(entry(hi, j, c), lo ? entry(*lo, j, c) : Scalar{});
            if (!x)
                return std::unexpected(x.error());
            data.x_magnitude.values()[offset + static_cast<std::size_t>(j)] = *x;
        }
        for (Index i = 0; i < n; ++i) {
            const auto at = offset + static_cast<std::size_t>(i);
            Accumulator real, imag;
            auto first = add_rhs(real, imag, data.rhs_hi.values()[at]);
            if (!first)
                return std::unexpected(first.error());
            auto second = add_rhs(real, imag, data.rhs_lo.values()[at]);
            if (!second)
                return std::unexpected(second.error());
            auto bmag =
                detail::actual_magnitude_bounds(data.rhs_hi.values()[at], data.rhs_lo.values()[at]);
            if (!bmag)
                return std::unexpected(bmag.error());
            auto denominator = *bmag;
            for (Index j = 0; j < n; ++j) {
                const Scalar a = coefficient(data, i, j);
                if (!finite(a))
                    return fail(FailureCode::Nonfinite);
                auto term = product(real, imag, a, entry(hi, j, c));
                if (!term)
                    return std::unexpected(term.error());
                if (lo) {
                    auto trailing = product(real, imag, a, entry(*lo, j, c));
                    if (!trailing)
                        return std::unexpected(trailing.error());
                }
                auto amag = detail::actual_magnitude_bounds(a, Scalar{});
                if (!amag)
                    return std::unexpected(amag.error());
                denominator = detail::add_magnitudes(
                    denominator,
                    detail::multiply_magnitudes(
                        *amag, data.x_magnitude.values()[offset + static_cast<std::size_t>(j)]));
            }
            auto r = real.result();
            if (!r)
                return std::unexpected(r.error());
            if constexpr (std::same_as<Scalar, double>) {
                data.residual.values()[at] = r->value;
                data.residual_error.values()[at] = r->error;
            } else {
                auto im = imag.result();
                if (!im)
                    return std::unexpected(im.error());
                data.residual.values()[at] = {r->value, im->value};
                // L1 of component error enclosures is a sufficient complex-modulus bound.
                data.residual_error.values()[at] =
                    detail::add_magnitudes({r->error, r->error}, {im->error, im->error}).upper;
            }
            data.denominator_lower.values()[at] = denominator.lower;
        }
        auto metric = detail::evaluate_column_metrics<Scalar>(
            n, data.residual.values().subspan(offset, size),
            data.residual_error.values().subspan(offset, size),
            data.denominator_lower.values().subspan(offset, size),
            data.rhs_hi.values().subspan(offset, size), data.rhs_lo.values().subspan(offset, size),
            scales.residual, accuracy);
        if (!metric)
            return std::unexpected(metric.error());
        const auto old_correction = metrics[static_cast<std::size_t>(c)].scaled_correction_upper;
        metric->scaled_correction_upper = std::isfinite(old_correction) ? old_correction : 0;
        metrics[static_cast<std::size_t>(c)] = *metric;
        accepted = accepted && detail::meets_accuracy(*metric, accuracy);
    }
    return accepted;
}
template <StorageScalar Scalar, class Data>
Result<Solution<Scalar>> solve(Data& data, const Rhs<Scalar>& rhs, const Scales& scales,
                               const Precision<Scalar>& precision, const Accuracy& accuracy,
                               const Limits& limits, core::ResourcePlan& plan, Index rank,
                               ConditionDiagnostic condition) noexcept {
    const Index n = data.original.shape().rows, columns = rhs.values().shape().columns;
    auto valid =
        detail::check_solve_inputs(n, rank, data.context, rhs, scales, precision, accuracy, limits);
    if (!valid)
        return std::unexpected(valid.error());
    auto capacity = detail::check_capacity_inputs(n, columns, rank, data.capacity);
    if (!capacity)
        return std::unexpected(capacity.error());
    if (!same_scales(scales, data.scales))
        return fail(FailureCode::StaleGeneration);
    auto admitted = admit_backend_memory(data.admission.plan(), limits.backend_memory);
    if (!admitted)
        return std::unexpected(admitted.error());
    auto lease = data.lock.acquire();
    if (!lease)
        return std::unexpected(lease.error());
    auto metadata = detail::SolutionFactory::prepare_metadata(n, columns, scales,
                                                              precision.required_checks, plan);
    if (!metadata)
        return std::unexpected(metadata.error());
    auto hi = DenseMatrix<Scalar>::copy(rhs.values().shape(), high(rhs.values()).values(), plan,
                                        operation);
    if (!hi)
        return std::unexpected(hi.error());
    auto lo = DenseMatrix<Scalar>::copy(rhs.values().shape(), high(rhs.values()).values(), plan,
                                        operation);
    if (!lo)
        return std::unexpected(lo.error());
    std::fill(lo->values().begin(), lo->values().end(), Scalar{});
    std::optional<DenseMatrix<Scalar>> projected;
    const auto* projection = std::get_if<Binary64Projection>(&precision.returned);
    const auto* retained_request = std::get_if<RetainedPrecision>(&precision.returned);
    if (!projection && !retained_request)
        return fail(FailureCode::PrecisionUnsupported);
    if (projection) {
        auto p = DenseMatrix<Scalar>::copy(rhs.values().shape(), high(rhs.values()).values(), plan,
                                           operation);
        if (!p)
            return std::unexpected(p.error());
        projected.emplace(std::move(*p));
    }
    bool two = rhs.values().precision() == core::Precision::TwoWord;
    RetentionPath path = two ? RetentionPath::InitiallyTwoWord : RetentionPath::Binary64;
    const auto* rhs_low = low(rhs.values());
    for (Index c = 0; c < columns; ++c) {
        const auto offset = static_cast<std::size_t>(c * n), size = static_cast<std::size_t>(n);
        for (Index i = 0; i < n; ++i) {
            const auto at = offset + static_cast<std::size_t>(i);
            data.rhs_hi.values()[at] = entry(high(rhs.values()), i, c);
            data.rhs_lo.values()[at] = rhs_low ? entry(*rhs_low, i, c) : Scalar{};
        }
        auto target =
            hi->values().subspan(static_cast<std::size_t>(c * hi->shape().leading_dimension), size);
        auto initial = inverse(
            data, std::span<const Scalar>{data.rhs_hi.values().subspan(offset, size)}, target);
        if (!initial)
            return std::unexpected(initial.error());
        if (two) {
            auto correction = data.correction.values().subspan(offset, size);
            auto second =
                inverse(data, std::span<const Scalar>{data.rhs_lo.values().subspan(offset, size)},
                        correction);
            if (!second)
                return std::unexpected(second.error());
            for (Index i = 0; i < n; ++i) {
                const auto at = static_cast<std::size_t>(c * hi->shape().leading_dimension + i);
                auto updated = uvlas::numerics::detail::apply_two_word_correction(
                    hi->values()[at], lo->values()[at], correction[static_cast<std::size_t>(i)]);
                if (!updated)
                    return std::unexpected(updated.error());
            }
        }
    }
    unsigned refinements = 0;
    while (true) {
        auto assessed = evaluate(data, *hi, two ? &*lo : nullptr, scales, accuracy,
                                 metadata->retained_columns.values());
        if (!assessed)
            return std::unexpected(assessed.error());
        if (*assessed)
            break;
        if (refinements == limits.maximum_refinements)
            return fail(FailureCode::LinearUnconverged);
        // Progress is judged in the original system, not only by unchanged rounded coordinates.
        if (!two && refinements > 0) {
            for (std::size_t c = 0; c < metadata->retained_columns.size(); ++c) {
                const auto now = metadata->retained_columns.values()[c],
                           before = metadata->returned_columns.values()[c];
                if (!detail::meets_accuracy(now, accuracy) &&
                    (now.normalized_original_residual_upper <= 1 ||
                     now.normalized_original_residual_upper >
                         0.5 * before.normalized_original_residual_upper) &&
                    (now.componentwise_backward_error_upper <=
                         accuracy.componentwise_backward_tolerance ||
                     now.componentwise_backward_error_upper >
                         0.5 * before.componentwise_backward_error_upper)) {
                    two = true;
                    path = RetentionPath::PromotedTwoWord;
                }
            }
        }
        std::copy(metadata->retained_columns.values().begin(),
                  metadata->retained_columns.values().end(),
                  metadata->returned_columns.values().begin());
        bool stagnated = false;
        for (Index c = 0; c < columns; ++c) {
            const auto offset = static_cast<std::size_t>(c * n), size = static_cast<std::size_t>(n);
            auto correction = data.correction.values().subspan(offset, size);
            auto correction_ok =
                inverse(data, std::span<const Scalar>{data.residual.values().subspan(offset, size)},
                        correction);
            if (!correction_ok)
                return std::unexpected(correction_ok.error());
            bool unchanged = true;
            for (Index i = 0; i < n; ++i) {
                const auto packed = offset + static_cast<std::size_t>(i),
                           at = static_cast<std::size_t>(c * hi->shape().leading_dimension + i);
                data.proposal.values()[packed] =
                    hi->values()[at] + correction[static_cast<std::size_t>(i)];
                if (!finite(data.proposal.values()[packed]))
                    return fail(FailureCode::Nonfinite);
                unchanged = unchanged && (data.proposal.values()[packed] == hi->values()[at]);
            }
            if (!detail::meets_accuracy(
                    metadata->retained_columns.values()[static_cast<std::size_t>(c)], accuracy) &&
                unchanged)
                stagnated = true;
        }
        if (!two && stagnated) {
            two = true;
            path = RetentionPath::PromotedTwoWord;
        }
        bool changed = false;
        for (Index c = 0; c < columns; ++c)
            for (Index i = 0; i < n; ++i) {
                const auto packed = static_cast<std::size_t>(c * n + i),
                           at = static_cast<std::size_t>(c * hi->shape().leading_dimension + i);
                const Scalar before_hi = hi->values()[at], before_lo = lo->values()[at];
                if (two) {
                    auto updated = uvlas::numerics::detail::apply_two_word_correction(
                        hi->values()[at], lo->values()[at], data.correction.values()[packed]);
                    if (!updated)
                        return std::unexpected(updated.error());
                } else
                    hi->values()[at] = data.proposal.values()[packed];
                changed = changed || before_hi != hi->values()[at] || before_lo != lo->values()[at];
            }
        // Record the actual solved correction in physical unknown scales, using the same
        // outward norm policy as residual qualification. These arrays are scratch here.
        for (Index c = 0; c < columns; ++c) {
            const auto offset = static_cast<std::size_t>(c * n), size = static_cast<std::size_t>(n);
            auto errors = data.residual_error.values().subspan(offset, size);
            auto denominators = data.denominator_lower.values().subspan(offset, size);
            std::fill(errors.begin(), errors.end(), 0.0);
            std::fill(denominators.begin(), denominators.end(), 1.0);
            auto correction_accuracy = accuracy;
            correction_accuracy.absolute_tolerance = 1;
            correction_accuracy.relative_tolerance = 0;
            auto metric = detail::evaluate_column_metrics<Scalar>(
                n, data.correction.values().subspan(offset, size), errors, denominators,
                data.rhs_hi.values().subspan(offset, size),
                data.rhs_lo.values().subspan(offset, size), scales.unknown, correction_accuracy);
            if (!metric)
                return std::unexpected(metric.error());
            metadata->retained_columns.values()[static_cast<std::size_t>(c)]
                .scaled_correction_upper = metric->normalized_original_residual_upper;
        }
        ++refinements;
        if (!changed)
            return fail(FailureCode::PrecisionUnresolved);
    }
    if (projected) {
        std::copy(metadata->retained_columns.values().begin(),
                  metadata->retained_columns.values().end(),
                  metadata->returned_columns.values().begin());
        for (Index c = 0; c < columns; ++c)
            for (Index i = 0; i < n; ++i) {
                const auto at = static_cast<std::size_t>(c * hi->shape().leading_dimension + i);
                projected->values()[at] = hi->values()[at] + (two ? lo->values()[at] : Scalar{});
                if (!finite(projected->values()[at]))
                    return fail(FailureCode::Nonfinite);
            }
        auto checked = evaluate(data, *projected, static_cast<const DenseMatrix<Scalar>*>(nullptr),
                                scales, accuracy, metadata->returned_columns.values());
        if (!checked)
            return std::unexpected(checked.error());
        if (!*checked)
            return fail(FailureCode::PrecisionUnresolved);
    }
    const auto generation =
        projection ? projection->retained_generation : retained_request->generation;
    auto retained =
        two ? PrecisionBuffer<Scalar>::two_word(std::move(*hi), std::move(*lo), generation)
            : PrecisionBuffer<Scalar>::binary64(std::move(*hi), generation);
    if (!retained)
        return std::unexpected(retained.error());
    std::optional<PrecisionBuffer<Scalar>> projected_buffer;
    if (projected) {
        auto p =
            PrecisionBuffer<Scalar>::binary64(std::move(*projected), projection->target_generation);
        if (!p)
            return std::unexpected(p.error());
        projected_buffer.emplace(std::move(*p));
    }
    return detail::SolutionFactory::finalize(
        std::move(*retained), std::move(projected_buffer), precision, std::move(*metadata),
        detail::Finalization{data.context, rhs.cache_request().rhs, rhs.values().generation(),
                             accuracy, data.admission, *admitted, refinements, path, condition});
}
} // namespace
Result<Solution<double>> solve_linear(const DenseOperator& op, const Rhs<double>& rhs,
                                      const Scales& scales, const Precision<double>& precision,
                                      const Accuracy& accuracy, const Limits& limits,
                                      DensePrepared& prepared, core::ResourcePlan& plan) noexcept {
    auto result = core::guard_entry(
        operation, FailureCode::FactorizationFailed, {}, [&]() -> Result<Solution<double>> {
            if (!prepared.data_)
                return fail(FailureCode::InvalidInput);
            auto finite_input = op.values.validate(operation);
            if (!finite_input)
                return std::unexpected(finite_input.error());
            if (op.context != prepared.data_->context ||
                !same_matrix(op.values, prepared.data_->original))
                return fail(FailureCode::StaleGeneration);
            return solve<double>(*prepared.data_, rhs, scales, precision, accuracy, limits, plan, 0,
                                 prepared.data_->condition);
        });
    return with_context(std::move(result), op.context);
}
Result<Solution<Complex>> solve_linear(const ColumnUpdatedOperator& op, const Rhs<Complex>& rhs,
                                       const Scales& scales, const Precision<Complex>& precision,
                                       const Accuracy& accuracy, const Limits& limits,
                                       UpdatedPrepared& prepared,
                                       core::ResourcePlan& plan) noexcept {
    auto result = core::guard_entry(
        operation, FailureCode::FactorizationFailed, {}, [&]() -> Result<Solution<Complex>> {
            if (!prepared.data_)
                return fail(FailureCode::InvalidInput);
            for (const auto* matrix : {&op.u, &op.v}) {
                auto valid = matrix->validate(operation);
                if (!valid)
                    return std::unexpected(valid.error());
            }
            auto valid = op.base.validate(operation);
            if (!valid)
                return std::unexpected(valid.error());
            if (op.context != prepared.data_->context ||
                !same_matrix(op.base, prepared.data_->original) ||
                !same_matrix(op.u, prepared.data_->u) || !same_matrix(op.v, prepared.data_->v))
                return fail(FailureCode::StaleGeneration);
            return solve<Complex>(*prepared.data_, rhs, scales, precision, accuracy, limits, plan,
                                  op.u.shape().columns, NotEvaluated{});
        });
    return with_context(std::move(result), op.context);
}
} // namespace uvlas::numerics::linear
