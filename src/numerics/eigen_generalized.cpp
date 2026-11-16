#include "backend_memory_detail.hpp"
#include "compensated.hpp"
#include "sparse_spd.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <uvlas/numerics/eigen_generalized.hpp>

namespace uvlas::numerics::eigen {
namespace {
constexpr auto operation = core::Operation::ExtractEigen;
using Buffer = core::PlannedBuffer<double>;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

core::Failure failure(core::FailureCode code) noexcept {
    return {code, operation};
}
core::Failure arithmetic_failure(core::Failure error) noexcept {
    error.operation = operation;
    if (error.code == core::FailureCode::Nonfinite)
        error.code = core::FailureCode::PrecisionUnresolved;
    return error;
}
// Construct only the rounded candidate. All its physical consequences are certified
// later from original matrices; exponent separation avoids spurious intermediate overflow.
core::Result<double> scaled(double value, double first, double second, double divisor,
                            double other_divisor = 1) noexcept {
    if (value == 0)
        return value;
    int exponent = 0, next = 0;
    double mantissa = std::frexp(value, &exponent);
    for (const double factor : {first, second}) {
        mantissa *= std::frexp(factor, &next);
        exponent += next;
    }
    for (const double factor : {divisor, other_divisor}) {
        mantissa /= std::frexp(factor, &next);
        exponent -= next;
    }
    const double result = std::scalbn(mantissa, exponent);
    if (!std::isfinite(result) || result == 0)
        return std::unexpected(failure(core::FailureCode::PrecisionUnresolved));
    return result;
}
// Compare |physical_i|/D_i without an overflowing product or rounded quotient.
// Normalized EFT cross-products preserve the smallest-index tie policy on actual values.
core::Result<bool> larger_coordinate(double first, double first_scale, double second,
                                     double second_scale) noexcept {
    if (first == 0)
        return false;
    if (second == 0)
        return true;
    const auto product = [](double value,
                            double scale) -> core::Result<std::pair<TwoWordValue, int>> {
        int a = 0, b = 0;
        const double x = std::frexp(std::abs(value), &a), y = std::frexp(scale, &b);
        auto words = two_product(x, y);
        if (!words)
            return std::unexpected(arithmetic_failure(words.error()));
        int exponent = a + b;
        if (words->hi < .5 || (words->hi == .5 && words->lo < 0)) {
            words->hi *= 2;
            words->lo *= 2;
            --exponent;
        }
        return std::pair{*words, exponent};
    };
    const auto left = product(first, second_scale), right = product(second, first_scale);
    if (!left)
        return std::unexpected(left.error());
    if (!right)
        return std::unexpected(right.error());
    if (left->second != right->second)
        return left->second > right->second;
    if (left->first.hi != right->first.hi)
        return left->first.hi > right->first.hi;
    return left->first.lo > right->first.lo;
}
class PhaseClock {
public:
    explicit PhaseClock(GeneralizedTimings& timings) noexcept : timings_(timings) {
        if (timings_.enabled)
            start_ = std::chrono::steady_clock::now();
    }
    ~PhaseClock() noexcept {
        finish();
    }
    PhaseClock(const PhaseClock&) = delete;
    PhaseClock& operator=(const PhaseClock&) = delete;
    void next(GeneralizedStage stage) noexcept {
        finish();
        stage_ = stage;
        if (timings_.enabled)
            start_ = std::chrono::steady_clock::now();
    }
    void finish() noexcept {
        if (!start_.has_value())
            return;
        const auto beginning = *start_;
        start_.reset();
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                 std::chrono::steady_clock::now() - beginning)
                                 .count();
        if (elapsed < 0) {
            healthy_ = false;
            return;
        }
        auto& value = timings_.stages[static_cast<std::size_t>(stage_)];
        const auto sample = static_cast<std::uint64_t>(elapsed);
        if (value.calls == std::numeric_limits<std::uint64_t>::max() ||
            sample > std::numeric_limits<std::uint64_t>::max() - value.elapsed_nanoseconds) {
            healthy_ = false;
            return;
        }
        ++value.calls;
        value.elapsed_nanoseconds += sample;
    }
    [[nodiscard]] bool healthy() const noexcept {
        return healthy_;
    }

private:
    GeneralizedTimings& timings_;
    GeneralizedStage stage_ = GeneralizedStage::Preparation;
    std::optional<std::chrono::steady_clock::time_point> start_;
    bool healthy_ = true;
};
struct Work {
    DenseReal stiffness, mass, basis;
    Buffer values, scales, stiffness_action, mass_action, stiffness_error, mass_error;
    core::PlannedBuffer<linear::UnitDimensions> units;
    core::PlannedBuffer<std::optional<core::ErrorCertificate>> consumers;
};
core::Result<Work> allocate_work(const GeneralizedOperator& stiffness,
                                 const GeneralizedOperator& mass, const GeneralizedRequest& request,
                                 core::Index n, core::Index square,
                                 core::ResourcePlan& plan) noexcept {
    auto k =
        DenseReal::copy(stiffness.values().shape(), stiffness.values().values(), plan, operation);
    if (!k)
        return std::unexpected(k.error());
    auto m = DenseReal::copy(mass.values().shape(), mass.values().values(), plan, operation);
    if (!m)
        return std::unexpected(m.error());
    auto raw = plan.allocate<double>(square, operation);
    if (!raw)
        return std::unexpected(raw.error());
    auto basis = DenseReal::copy({n, n, n}, raw->values(), plan, operation);
    if (!basis)
        return std::unexpected(basis.error());
    auto values = plan.allocate<double>(n, operation);
    if (!values)
        return std::unexpected(values.error());
    auto scales = plan.allocate<double>(n, operation);
    if (!scales)
        return std::unexpected(scales.error());
    auto units = plan.allocate<linear::UnitDimensions>(n, operation);
    if (!units)
        return std::unexpected(units.error());
    auto ka = plan.allocate<double>(square, operation);
    if (!ka)
        return std::unexpected(ka.error());
    auto ma = plan.allocate<double>(square, operation);
    if (!ma)
        return std::unexpected(ma.error());
    auto ke = plan.allocate<double>(square, operation);
    if (!ke)
        return std::unexpected(ke.error());
    auto me = plan.allocate<double>(square, operation);
    if (!me)
        return std::unexpected(me.error());
    if (std::cmp_greater(request.precision.required_checks.size(),
                         std::numeric_limits<core::Index>::max()))
        return std::unexpected(failure(core::FailureCode::ResourceExceeded));
    auto checks = plan.allocate<std::optional<core::ErrorCertificate>>(
        static_cast<core::Index>(request.precision.required_checks.size()), operation);
    if (!checks)
        return std::unexpected(checks.error());
    std::copy(request.scales.metric.components.begin(), request.scales.metric.components.end(),
              scales->values().begin());
    std::copy(request.scales.component_units.begin(), request.scales.component_units.end(),
              units->values().begin());
    return Work{std::move(*k),      std::move(*m),     std::move(*basis), std::move(*values),
                std::move(*scales), std::move(*ka),    std::move(*ma),    std::move(*ke),
                std::move(*me),     std::move(*units), std::move(*checks)};
}
core::Result<void> charge(std::optional<std::uint64_t>& counter, std::uint64_t limit,
                          core::ResourceKind kind) noexcept {
    const auto started = counter.value_or(0);
    if (started >= limit) {
        auto error = failure(core::FailureCode::ResourceExceeded);
        if (started < std::numeric_limits<std::uint64_t>::max())
            error.details.resources = core::ResourceUsage{kind, started + 1, started};
        return std::unexpected(error);
    }
    counter = started + 1;
    return {};
}
core::Result<void> actions(Work& work, const ActionLimits& limits,
                           AttemptDiagnostics& diagnostics) noexcept {
    const auto n = static_cast<std::size_t>(work.basis.shape().rows);
    const auto phi = work.basis.values();
    for (std::size_t column = 0; column < n; ++column) {
        const auto charged = charge(diagnostics.forward_applications, limits.forward,
                                    core::ResourceKind::OperatorApplications);
        if (!charged)
            return charged;
        for (std::size_t row = 0; row < n; ++row) {
            numerics::detail::CompensatedAccumulator stiffness, mass;
            for (std::size_t inner = 0; inner < n; ++inner) {
                const auto kindex =
                    inner * static_cast<std::size_t>(work.stiffness.shape().leading_dimension) +
                    row;
                const auto mindex =
                    inner * static_cast<std::size_t>(work.mass.shape().leading_dimension) + row;
                auto status =
                    stiffness.product(work.stiffness.values()[kindex], phi[column * n + inner], 1);
                if (!status)
                    return std::unexpected(arithmetic_failure(status.error()));
                status = mass.product(work.mass.values()[mindex], phi[column * n + inner], 1);
                if (!status)
                    return std::unexpected(arithmetic_failure(status.error()));
            }
            const auto kvalue = stiffness.result();
            const auto mvalue = mass.result();
            if (!kvalue)
                return std::unexpected(arithmetic_failure(kvalue.error()));
            if (!mvalue)
                return std::unexpected(arithmetic_failure(mvalue.error()));
            const auto index = column * n + row;
            work.stiffness_action.values()[index] = kvalue->value;
            work.mass_action.values()[index] = mvalue->value;
            work.stiffness_error.values()[index] = kvalue->error;
            work.mass_error.values()[index] = mvalue->error;
        }
    }
    return {};
}
core::Result<void> validate(const GeneralizedOperator& stiffness, const GeneralizedOperator& mass,
                            const GeneralizedRequest& request) noexcept {
    const auto ks = stiffness.values().shape(), ms = mass.values().shape();
    if (request.measurement != GeneralizedMeasurement::Disabled &&
        request.measurement != GeneralizedMeasurement::Enabled)
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    if (ks.rows != ks.columns || ms.rows != ms.columns || ks.rows != ms.rows ||
        request.limits.maximum_order < 0 || request.limits.maximum_results < 0)
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    if (ks.rows > request.limits.maximum_order || ks.rows > request.limits.maximum_results)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded));
    const auto finite_k = stiffness.values().validate(operation);
    if (!finite_k)
        return finite_k;
    const auto finite_m = mass.values().validate(operation);
    if (!finite_m)
        return finite_m;
    const auto kc = stiffness.context(), mc = mass.context();
    if (!std::isfinite(kc.time) || !std::isfinite(mc.time) ||
        !std::isfinite(request.scales.context.time))
        return std::unexpected(failure(core::FailureCode::Nonfinite));
    auto comparable = mc;
    comparable.cache.matrix = kc.cache.matrix;
    if (comparable != kc || request.scales.context != kc)
        return std::unexpected(failure(core::FailureCode::StaleGeneration));
    if (request.precision.input != core::Precision::Binary64 ||
        kc.cache.representation != core::Precision::Binary64 ||
        mc.cache.representation != core::Precision::Binary64)
        return std::unexpected(failure(core::FailureCode::PrecisionUnsupported));
    if (request.require_physical_coverage)
        return std::unexpected(failure(core::FailureCode::EigenCoverageUnresolved));
    if (ks.rows == 0 && !request.precision.required_checks.empty())
        return std::unexpected(failure(core::FailureCode::PrecisionUnresolved));
    const auto n = static_cast<std::size_t>(ks.rows);
    if (request.scales.metric.components.size() != n || request.scales.component_units.size() != n)
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    for (const auto x : request.scales.metric.components) {
        if (!std::isfinite(x))
            return std::unexpected(failure(core::FailureCode::Nonfinite));
        if (!(x > 0))
            return std::unexpected(failure(core::FailureCode::InvalidInput));
    }
    for (const auto x : {request.scales.metric.eigenvalue_scale, request.scales.mass_scale}) {
        if (!std::isfinite(x))
            return std::unexpected(failure(core::FailureCode::Nonfinite));
        if (!(x > 0))
            return std::unexpected(failure(core::FailureCode::InvalidInput));
    }
    const auto valid_metric = [](MetricRequirement item) -> core::Result<void> {
        if (!std::isfinite(item.budget))
            return std::unexpected(failure(core::FailureCode::Nonfinite));
        if (item.budget < 0)
            return std::unexpected(failure(core::FailureCode::InvalidInput));
        return {};
    };
    for (const auto item : {request.accuracy.residual, request.accuracy.normalization,
                            request.accuracy.orthogonality}) {
        const auto status = valid_metric(item);
        if (!status)
            return status;
    }
    if (request.accuracy.residual.budget > 1.e-10)
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    if (const auto* item = std::get_if<MetricRequirement>(&request.accuracy.backward)) {
        const auto status = valid_metric(*item);
        if (!status)
            return status;
    }
    if (!request.precision.required_checks.empty() && request.precision.evaluator == nullptr)
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    for (const auto& check : request.precision.required_checks) {
        const auto status = valid_metric({check.budget, check.observation, check.scope});
        if (!status)
            return status;
    }
    if (const auto* projection =
            std::get_if<linear::Binary64Projection>(&request.precision.returned)) {
        if (projection->retained_generation == projection->target_generation)
            return std::unexpected(failure(core::FailureCode::StaleGeneration));
    }
    for (core::Index column = 0; column < ks.rows; ++column) {
        for (core::Index row = 0; row < ks.rows; ++row) {
            const auto i = static_cast<std::size_t>(column * ks.leading_dimension + row);
            const auto j = static_cast<std::size_t>(row * ks.leading_dimension + column);
            const auto mi = static_cast<std::size_t>(column * ms.leading_dimension + row);
            const auto mj = static_cast<std::size_t>(row * ms.leading_dimension + column);
            if (stiffness.values().values()[i] != stiffness.values().values()[j] ||
                mass.values().values()[mi] != mass.values().values()[mj])
                return std::unexpected(failure(core::FailureCode::InvalidInput));
        }
    }
    return {};
}
} // namespace

GeneralizedSpectrum::GeneralizedSpectrum(
    DenseReal stiffness, DenseReal mass, DenseReal basis, core::PlannedBuffer<double> values,
    core::PlannedBuffer<double> scales, core::PlannedBuffer<linear::UnitDimensions> units,
    core::PlannedBuffer<std::optional<core::ErrorCertificate>> consumers, Context stiffness_context,
    Context mass_context, core::NumericalGeneration generation, GeneralizedMetrics metrics,
    double shift, double eigenvalue_scale, double mass_scale,
    linear::UnitDimensions eigenvalue_units, linear::UnitDimensions mass_units,
    AttemptDiagnostics diagnostics, linear::BackendAdmission admission,
    StateDiagnosticContext state, linear::ReturnedPrecision returned,
    GeneralizedLimits limits) noexcept
    : stiffness_(std::move(stiffness)),
      mass_(std::move(mass)),
      basis_(std::move(basis)),
      values_(std::move(values)),
      scales_(std::move(scales)),
      units_(std::move(units)),
      consumers_(std::move(consumers)),
      stiffness_context_(stiffness_context),
      mass_context_(mass_context),
      generation_(generation),
      metrics_(metrics),
      mass_shift_(shift),
      eigenvalue_scale_(eigenvalue_scale),
      mass_scale_(mass_scale),
      eigenvalue_units_(eigenvalue_units),
      mass_units_(mass_units),
      diagnostics_(diagnostics),
      admission_(admission),
      limits_(limits),
      state_(state),
      returned_(returned) {}
core::Result<core::ErrorCertificate>
GeneralizedSpectrum::consumer(std::size_t index) const noexcept {
    if (index >= consumers_.size())
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    const auto& proof = consumers_.values()[index];
    if (!proof.has_value())
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    return *proof;
}
namespace detail {
class GeneralizedFactory {
public:
    static void record_timings(GeneralizedSpectrum& result, GeneralizedTimings timings) noexcept {
        result.timings_ = timings;
    }
    static Result<GeneralizedSpectrum> packet(Work work, const GeneralizedRequest& request,
                                              Context stiffness, Context mass,
                                              core::NumericalGeneration generation, double shift,
                                              AttemptDiagnostics& diagnostics,
                                              linear::BackendAdmission admission) noexcept {
        const auto applied = actions(work, request.limits.actions, diagnostics);
        if (!applied)
            return std::unexpected(
                complete_failure(applied.error(), stiffness, request.state, diagnostics));
        const auto phi = WordView<double>::make(work.basis.values(), {}, core::Precision::Binary64);
        const auto values =
            WordView<double>::make(work.values.values(), {}, core::Precision::Binary64);
        const auto ka =
            WordView<double>::make(work.stiffness_action.values(), {}, core::Precision::Binary64);
        const auto ma =
            WordView<double>::make(work.mass_action.values(), {}, core::Precision::Binary64);
        if (!phi || !values || !ka || !ma)
            return std::unexpected(complete_failure(failure(core::FailureCode::PrecisionUnresolved),
                                                    stiffness, request.state, diagnostics));
        auto proof = assess_generalized_metrics(
            work.basis.shape().rows, *phi, *values, {*ka, work.stiffness_error.values()},
            {*ma, work.mass_error.values()}, request.scales.metric, request.scales.mass_scale,
            request.accuracy, {stiffness, request.state, generation});
        if (!proof)
            return std::unexpected(proof.error());
        return GeneralizedSpectrum{std::move(work.stiffness),
                                   std::move(work.mass),
                                   std::move(work.basis),
                                   std::move(work.values),
                                   std::move(work.scales),
                                   std::move(work.units),
                                   std::move(work.consumers),
                                   stiffness,
                                   mass,
                                   generation,
                                   *proof,
                                   shift,
                                   request.scales.metric.eigenvalue_scale,
                                   request.scales.mass_scale,
                                   request.scales.eigenvalue_units,
                                   request.scales.mass_units,
                                   diagnostics,
                                   admission,
                                   request.state,
                                   request.precision.returned,
                                   request.limits};
    }
    static Result<void> consumers(GeneralizedSpectrum& returned,
                                  const GeneralizedSpectrum& retained,
                                  const GeneralizedRequest& request) noexcept {
        for (std::size_t i = 0; i < request.precision.required_checks.size(); ++i) {
            const auto& check = request.precision.required_checks[i];
            const core::CertificateContext expected{
                check.observation, check.scope, returned.generation_, core::Precision::Binary64};
            const auto result =
                core::guard_entry(operation, core::FailureCode::PrecisionUnresolved, {}, [&] {
                    return request.precision.evaluator->assess(retained, returned, check, expected);
                });
            if (!result)
                return std::unexpected(complete_failure(result.error(), returned.stiffness_context_,
                                                        request.state, returned.diagnostics_));
            if (result->context() != expected ||
                std::bit_cast<std::uint64_t>(result->budget()) !=
                    std::bit_cast<std::uint64_t>(check.budget) ||
                !std::holds_alternative<core::Achieved>(result->assessment()))
                return std::unexpected(complete_failure(
                    failure(core::FailureCode::PrecisionUnresolved), returned.stiffness_context_,
                    request.state, returned.diagnostics_));
            returned.consumers_.values()[i] = *result;
        }
        return {};
    }
};
} // namespace detail

Result<GeneralizedSpectrum> symmetric_generalized_spectrum(const GeneralizedOperator& stiffness,
                                                           const GeneralizedOperator& mass,
                                                           const GeneralizedRequest& request,
                                                           core::ResourcePlan& plan) noexcept {
    AttemptDiagnostics diagnostics;
    return guard_entry<GeneralizedSpectrum>(
        stiffness.context(), request.state, diagnostics, [&]() -> Result<GeneralizedSpectrum> {
            const auto reject = [&](core::Failure error) -> Result<GeneralizedSpectrum> {
                return std::unexpected(
                    complete_failure(error, stiffness.context(), request.state, diagnostics));
            };
            GeneralizedTimings timings;
            timings.enabled = request.measurement == GeneralizedMeasurement::Enabled;
            PhaseClock timer(timings);
            const auto valid = validate(stiffness, mass, request);
            if (!valid)
                return reject(valid.error());
            const auto n = stiffness.values().shape().rows;
            const auto stopping_limit =
                core::checked_extent_product(std::array{n, core::Index{30}}, operation);
            if (!stopping_limit)
                return reject(stopping_limit.error());
            const auto square = core::checked_extent_product(std::array{n, n}, operation);
            if (!square)
                return reject(square.error());
            // Forecast only: dense candidates/factors/solver 16*n*n doubles;
            // original-M ordering, permuted/shifted CSC and sparse LLT fill 64*n*n;
            // backend temporaries 16*n*n. Owned arrays use plan separately.
            const auto backend_count =
                core::checked_extent_product(std::array{n, n, core::Index{96}}, operation);
            if (!backend_count)
                return reject(backend_count.error());
            const auto backend_bytes =
                core::checked_byte_count(*backend_count, sizeof(double), operation);
            if (!backend_bytes)
                return reject(backend_bytes.error());
            const auto admission = linear::admit_backend_memory(
                linear::detail::BackendPlanFactory::make(
                    *backend_bytes, n == 0 ? std::optional<std::uint64_t>{0} : std::nullopt),
                request.limits.backend_memory);
            if (!admission)
                return reject(admission.error());
            diagnostics.requested_count = static_cast<std::uint64_t>(n);
            diagnostics.extraction_attempts = 0;
            diagnostics.forward_applications = 0;
            diagnostics.adjoint_applications = 0;
            diagnostics.inverse_applications = 0;
            auto work = allocate_work(stiffness, mass, request, n, *square, plan);
            if (!work)
                return reject(work.error());
            std::optional<Work> projected;
            if (std::holds_alternative<linear::Binary64Projection>(request.precision.returned)) {
                auto next = allocate_work(stiffness, mass, request, n, *square, plan);
                if (!next)
                    return reject(next.error());
                projected = std::move(*next);
            }
            double mass_shift = 0;
            // Own the canonical original CSC before asking the rigorous sparse qualifier.
            auto offsets = plan.allocate<core::Index>(n == 0 ? 0 : n + 1, operation);
            if (!offsets)
                return reject(offsets.error());
            auto indices = plan.allocate<core::Index>(*square, operation);
            if (!indices)
                return reject(indices.error());
            auto coefficients = plan.allocate<double>(*square, operation);
            if (!coefficients)
                return reject(coefficients.error());
            for (core::Index column = 0; column < n; ++column) {
                offsets->values()[static_cast<std::size_t>(column)] = column * n;
                for (core::Index row = 0; row < n; ++row) {
                    const auto packed = static_cast<std::size_t>(column * n + row);
                    indices->values()[packed] = row;
                    coefficients->values()[packed] =
                        mass.values().values()[static_cast<std::size_t>(
                            column * mass.values().shape().leading_dimension + row)];
                }
            }
            if (n > 0) {
                offsets->values()[static_cast<std::size_t>(n)] = *square;
                timer.next(GeneralizedStage::MassQualification);
                auto canonical = SparseReal::copy_csc(n, n, offsets->values(), indices->values(),
                                                      coefficients->values(), plan, operation);
                if (!canonical)
                    return reject(canonical.error());
                const auto qualification =
                    linear::detail::qualify_sparse_positive_definite(*canonical, plan);
                if (!qualification) {
                    auto error = arithmetic_failure(qualification.error());
                    if (error.code == core::FailureCode::PrecisionUnresolved)
                        error.code = core::FailureCode::RankUnresolved;
                    return reject(error);
                }
                if (!qualification->certified || !qualification->shift.has_value())
                    return reject(failure(core::FailureCode::RankUnresolved));
                mass_shift = *qualification->shift;
                timer.next(GeneralizedStage::Backend);
                Matrix k(n, n), m(n, n);
                for (core::Index column = 0; column < n; ++column) {
                    for (core::Index row = column; row < n; ++row) {
                        const auto drow =
                            request.scales.metric.components[static_cast<std::size_t>(row)];
                        const auto dcolumn =
                            request.scales.metric.components[static_cast<std::size_t>(column)];
                        const auto kval = scaled(
                            stiffness.values().values()[static_cast<std::size_t>(
                                column * stiffness.values().shape().leading_dimension + row)],
                            drow, dcolumn, request.scales.mass_scale,
                            request.scales.metric.eigenvalue_scale);
                        const auto mval =
                            scaled(mass.values().values()[static_cast<std::size_t>(
                                       column * mass.values().shape().leading_dimension + row)],
                                   drow, dcolumn, request.scales.mass_scale);
                        if (!kval)
                            return reject(kval.error());
                        if (!mval)
                            return reject(mval.error());
                        k(row, column) = k(column, row) = *kval;
                        m(row, column) = m(column, row) = *mval;
                    }
                }
                // Candidate M is rounded and needs its own backend precondition check.
                const Eigen::LLT<Matrix> mass_factor(m);
                if (mass_factor.info() != Eigen::Success || !mass_factor.matrixLLT().allFinite())
                    return reject(failure(core::FailureCode::PrecisionUnresolved));
                for (core::Index i = 0; i < n; ++i)
                    if (!(mass_factor.matrixLLT()(i, i) > 0))
                        return reject(failure(core::FailureCode::PrecisionUnresolved));
                const auto charged =
                    charge(diagnostics.extraction_attempts, request.limits.actions.attempts,
                           core::ResourceKind::Iterations);
                if (!charged)
                    return reject(charged.error());
                const Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(k, m);
                if (solver.info() != Eigen::Success)
                    return reject(failure(core::FailureCode::EigenCoverageUnresolved));
                if (!solver.eigenvalues().allFinite() || !solver.eigenvectors().allFinite())
                    return reject(failure(core::FailureCode::FactorizationFailed));
                for (core::Index column = 0; column < n; ++column) {
                    const auto value = scaled(solver.eigenvalues()[column],
                                              request.scales.metric.eigenvalue_scale, 1, 1);
                    if (!value)
                        return reject(value.error());
                    work->values.values()[static_cast<std::size_t>(column)] = *value;
                    for (core::Index row = 0; row < n; ++row) {
                        const auto entry = scaled(
                            solver.eigenvectors()(row, column),
                            request.scales.metric.components[static_cast<std::size_t>(row)], 1, 1);
                        if (!entry)
                            return reject(entry.error());
                        work->basis.values()[static_cast<std::size_t>(column * n + row)] = *entry;
                    }
                    core::Index pivot = 0;
                    for (core::Index row = 1; row < n; ++row) {
                        const auto larger = larger_coordinate(
                            work->basis.values()[static_cast<std::size_t>(column * n + row)],
                            request.scales.metric.components[static_cast<std::size_t>(row)],
                            work->basis.values()[static_cast<std::size_t>(column * n + pivot)],
                            request.scales.metric.components[static_cast<std::size_t>(pivot)]);
                        if (!larger)
                            return reject(larger.error());
                        if (*larger)
                            pivot = row;
                    }
                    if (work->basis.values()[static_cast<std::size_t>(column * n + pivot)] < 0)
                        for (core::Index row = 0; row < n; ++row)
                            work->basis.values()[static_cast<std::size_t>(column * n + row)] *= -1;
                }
            }
            diagnostics.converged_count = static_cast<std::uint64_t>(n);
            const core::NumericalGeneration retained_generation = std::visit(
                [](const auto& selection) {
                    if constexpr (std::same_as<std::decay_t<decltype(selection)>,
                                               linear::RetainedPrecision>)
                        return selection.generation;
                    else
                        return selection.retained_generation;
                },
                request.precision.returned);
            if (projected.has_value()) {
                std::copy(work->values.values().begin(), work->values.values().end(),
                          projected->values.values().begin());
                std::copy(work->basis.values().begin(), work->basis.values().end(),
                          projected->basis.values().begin());
            }
            timer.next(GeneralizedStage::Qualification);
            auto retained = detail::GeneralizedFactory::packet(
                std::move(*work), request, stiffness.context(), mass.context(), retained_generation,
                mass_shift, diagnostics, *admission);
            if (!retained)
                return retained;
            if (projected.has_value()) {
                const auto target = std::get<linear::Binary64Projection>(request.precision.returned)
                                        .target_generation;
                auto result = detail::GeneralizedFactory::packet(
                    std::move(*projected), request, stiffness.context(), mass.context(), target,
                    mass_shift, diagnostics, *admission);
                if (!result)
                    return result;
                if (!request.precision.required_checks.empty())
                    timer.next(GeneralizedStage::Consumer);
                const auto consumed =
                    detail::GeneralizedFactory::consumers(*result, *retained, request);
                if (!consumed)
                    return std::unexpected(consumed.error());
                timer.finish();
                if (!timer.healthy())
                    return reject(failure(core::FailureCode::ResourceExceeded));
                detail::GeneralizedFactory::record_timings(*result, timings);
                return result;
            }
            if (!request.precision.required_checks.empty())
                timer.next(GeneralizedStage::Consumer);
            const auto consumed =
                detail::GeneralizedFactory::consumers(*retained, *retained, request);
            if (!consumed)
                return std::unexpected(consumed.error());
            timer.finish();
            if (!timer.healthy())
                return reject(failure(core::FailureCode::ResourceExceeded));
            detail::GeneralizedFactory::record_timings(*retained, timings);
            return retained;
        });
}
} // namespace uvlas::numerics::eigen
